#!/usr/bin/env perl
package AI::GMTK::EmConvergeParallel;
use warnings;
use strict;

use Cwd;
use File::Path;
use File::Copy;
use Carp;
use OS::Util qw(readLines writeLines tsSystem );
use Config::OptionsSet::OptionsSet qw(readOpts writeOpts prettyPrintOpts overwriteOpts 
					setDiff hashDiff hashSymDiff);
use AI::GMTK::Util  qw(	prepareGmtkCmdArgs checkOpts
					isRestartable notifyOfJobCompletion);
use AI::GMTK::EmTrainParallel qw(emtrainParallel isTrainingNecessary);

use Exporter;
our @ISA=qw(Exporter);
our @EXPORT_OK=qw(emConvergeParallel);

## @fn
#do em training to converge with gaussian split/merging.  If the execution is interrupted, it can 
#be safely rerun again and will start where it left of.
#See emConvergeParallel.pl use help for details.
#
#Return values:
#	0  success
#	-1 cofiguration error. no files written or modified.
# 	-2 sanity check failed.  No parallel jobs launched.
#	-3 some parallel tasks failed.
#
#	calling this function will also modify following keys of $userDistrOpt
#   $userDistrOpt->{'jobName'} will be set if it was previously null
sub emConvergeParallel{
    my $opt = shift;
	
	$| =1; #flush pipes right away
	
    my ($err, $gmtkOpt, $distrOpt) = checkOpts($opt, \&initDistrOptsEmConverge, \&fillInOptsEmConverge);

	if ($err){
		return $err;
	}
	
	
	#check to see if the job is restartable.  A job is restartable if the 
	#distrOpt and gmtkOpt exist in the working directory and are identical to
	#the distrOpt and gtmkOpt we've set up so far
	my @allowedDiffs=(
        'clearWdOnStart',
        'clearAccOnEnd',
        'clearWdOnEnd',
        'breakRangePath',
        'qsubPath',
        'numChunks',
        'uttRanges',
        'verbosity',
        'emailAddress',
        'sendMailOnCompletion',
        'jobName',
		'extraDistributeParameters');
	my $restartCode = isRestartable($gmtkOpt, $distrOpt, \@allowedDiffs);
	return 0 if (!isTrainingNecessary($gmtkOpt, $restartCode));

	print "Beginning split/vanish-converge steps at ". localtime(time).".\n";


	#from now on we are no longer read-only and may modify the fs.

	if($distrOpt->{'clearWdOnStart'} && -d $distrOpt->{'wd'}){
		print "clear working directory forced. deleting working directory $distrOpt->{'wd'}\n";
		rmtree($distrOpt->{'wd'});
        $restartCode = isRestartable($gmtkOpt, $distrOpt, \@allowedDiffs);
	}
	
	
	if ($restartCode==2){
		return -2;
	}

	
	#create dirs only if they don't exist already
	mkpath($distrOpt->{'wd'}, 1);

	if($restartCode == 0){	#first-time running
	
		#write out the config state so that the experiment is reproducible
		writeOpts($gmtkOpt, $distrOpt->{'_gmtkArgsBackup'},': ');
		writeOpts($distrOpt, $distrOpt->{'_distrArgsBackup'},' ');
	
	}

	my $ret = doAllSplitVanishSteps($gmtkOpt, $distrOpt);
	return $ret if ($ret);


	#clear workingDir at successful end
	if($distrOpt->{'clearWdOnEnd'}){	
		print "deleting working directory $distrOpt->{'wd'} per clearWdOnEnd flag\n";
		rmtree($distrOpt->{'wd'});
		
	}


	print "All split/vanish-converge steps finished at ". localtime(time).".\n";
	return 0;
}

## @fn
#guess some defaults for distributed environment options.
sub initDistrOptsEmConverge{
	my $opt = shift;

	#accept all options that emtrainParallel accepts
	AI::GMTK::EmTrainParallel::initDistrOptsEmTrain($opt);

	#the name of the job as it will appear in the SGE queue monitoring tools
    #the name is not used anywhere else.
	$opt->{'jobName'}='emConverge';

	#an array of split/vanishes for each step.  See gmtk doc for meanings of mcvr mcsr, etc.
	#by default do no splits or vanishes.  
	$opt->{'mcvrArray'} = '';
	$opt->{'mcsrArray'} = '';
	$opt->{'botForceVanishArray'} = '';
	$opt->{'topForceSplitArray'} = '';

	#before each ith split/vanish, iterate until log-likelihood 
	#difference is less than ith llRatioThreshold entry
	#the last llRatioThreshold entry is the convergence schedule after the last split
	#there should be one more llRatioThreshold entries than splitMergeSchedule entries 
	#an example of a legal entry $opt->{'llRatioThreshold'} = '.02 .02 .02 .02 .002'
	$opt->{'llRatioThresholdArray'} = '.002';

	#same rules as llRatioThreshold but the maximum iterations to perform for each converge step
	##an example of a legal entry $opt->{'maxEmIterations'} = '20 20 20 20 20'
	$opt->{'maxEmIterationsArray'} = '20';

	$opt->{'verbosity'} = 2;

}



sub fillInOptsEmConverge{

	my $gmtkOpt = shift;
	my $distrOpt = shift;
	my $ret=0;


	#check options
	my $badThingDefined =0;
	map {$badThingDefined |= defined($gmtkOpt->{$_})} ('mcvr',  'mcsr', 'botForceVanish', 'topForceSplit', 'botForceV', 'topForceS');
	if($badThingDefined)	
	{
		print "mcvr, mcsr, botForceVanish and topForceSplit (and their synonyms) are disallowed in gmtk args file for emConvergeParallel\n";
		print "use the split/vanish schedule instead.\n";
		$ret++;
	}

	#convert strings into lists for those that are needed 
	foreach('mcvr', 'mcsr','botForceVanish', 'topForceSplit'){
		my @vals = split(' ',$distrOpt->{$_.'Array'});
		unshift (@vals, undef);
		$distrOpt->{"_$_"}=\@vals;
	}
	#convert strings into lists for those that are needed 
	foreach('llRatioThreshold', 'maxEmIterations'){
		my @vals = split(' ',$distrOpt->{$_.'Array'});
		$distrOpt->{"_$_"}=\@vals;
	}

	if (@{$distrOpt->{"_mcvr"}} != @{$distrOpt->{"_mcsr"}} || 
		@{$distrOpt->{"_mcvr"}} != @{$distrOpt->{"_botForceVanish"}} ||
		@{$distrOpt->{"_mcvr"}} != @{$distrOpt->{"_topForceSplit"}}){
		print 'mcvrArray, mcsrArray, botForceVanishArray and topForceSplitArray must have the same number of elements but their values are' . 
                "(" .   scalar(@{$distrOpt->{"_mcvr"}}).", " . 
                        scalar(@{$distrOpt->{"_mcsr"}}).", " .
                        scalar(@{$distrOpt->{"_botForceVanish"}}).", " .
                        scalar(@{$distrOpt->{"_topForceSplit"}}).")\n" ;
		$ret++;
	}

	if (@{$distrOpt->{"_llRatioThreshold"}} != @{$distrOpt->{"_maxEmIterations"}} || 
		@{$distrOpt->{"_mcvr"}} != @{$distrOpt->{"_llRatioThreshold"}}) {
		print 'llRatioThresholdArray and maxEmIterationsArray must have the same number of elements and they must be one more than mcvrArray number of elements'."\n";
		$ret++;
	}

	if (!$distrOpt->{'wd'}){
		print "working directory must be specified\n";
		$ret++;
	}

	if (!$gmtkOpt->{'storeAccFile'} && !$gmtkOpt->{'outputTrainable'}){
		print "outputTrainable or storeAccFile is required in gmtk args file.\n";
		$ret++;
	}

	return $ret if ($ret>0);

	#filenames for saving state
	$distrOpt->{'_gmtkArgsBackup'} = "$distrOpt->{'wd'}/gmtk.args";
	$distrOpt->{'_distrArgsBackup'} = "$distrOpt->{'wd'}/distr.args";

	return $ret;
}


sub doAllSplitVanishSteps{
	my $gmtkOpt = shift;
	my $distrOpt = shift;

	print "performing ".scalar(@{$distrOpt->{"_mcvr"}})." split/vanish steps\n";


	my %stepGmtkOpt = %$gmtkOpt;
	my %stepDistrOpt = %$distrOpt;

	foreach my $svStep (0..$#{$distrOpt->{"_llRatioThreshold"}}){
		print "doing split/vanish step $svStep \n";
		%stepGmtkOpt = %$gmtkOpt;
		%stepDistrOpt = %$distrOpt;
		
		#settings specific to emConverge must not be passed down to emTrain
		my @badArgs = (grep (/^_/, keys %$distrOpt), 'mcvrArray', 'mcsrArray', 'botForceVanishArray', 'topForceSplitArray', 'llRatioThresholdArray', 'maxEmIterationsArray');
		delete @stepDistrOpt{@badArgs};

		$stepDistrOpt{'wd'} = $distrOpt->{'wd'}."/svStep$svStep";
		$stepDistrOpt{'jobName'}=$distrOpt->{'jobName'}."_svStep$svStep";

		#these need to be done once per execution of emConvergeParallel  
		$stepDistrOpt{'clearWdOnStart'}=0;
		$stepDistrOpt{'clearWdOnEnd'}=0;
		#$stepDistrOpt->{'clearAccOnEnd'}=0; might be good to do after every iteration
		$stepDistrOpt{'sendMailOnCompletion'}=0;
		if($gmtkOpt->{'outputTrainable'}){
			$stepGmtkOpt{'outputTrainable'}="$distrOpt->{'wd'}/learnedParams_svStep$svStep.gmtk";
		}
		$stepGmtkOpt{'llStoreFile'}="$distrOpt->{'wd'}/logLikelihood_svStep$svStep.txt";
		$stepDistrOpt{'verbosity'}=$distrOpt->{'verbosity'}-1;

		if ($svStep > 0){
			#$stepDistrOpt{'doSanityCheck'}=0;
			$stepGmtkOpt{'inputTrainable'}="$distrOpt->{'wd'}/learnedParams_svStep". ($svStep-1) .".gmtk";
			$stepDistrOpt{'verbosity'}=$distrOpt->{'verbosity'}-1;
 		}
		else{
			print "Converging the first time without doing any splits or vanishes\n";
		}

		my $svOpts = undef;
		if (defined( $distrOpt->{"_mcvr"}[$svStep])){
			$svOpts = 	{'mcvr' => $distrOpt->{"_mcvr"}[$svStep],
						'mcsr' => $distrOpt->{"_mcsr"}[$svStep],
						'botForceVanish' => $distrOpt->{"_botForceVanish"}[$svStep],
						'topForceSplit' => $distrOpt->{"_topForceSplit"}[$svStep]};
 		}
		my $convergenceThresholds = {'llRatioThreshold' => $distrOpt->{"_llRatioThreshold"}[$svStep],
							 'maxEmIterations' => $distrOpt->{"_maxEmIterations"}[$svStep]};




		my $ret = doSplitVanishStep(\%stepGmtkOpt, \%stepDistrOpt, $convergenceThresholds, $svOpts);

		if($ret == 0){
			print "split/vanish step $svStep succeeded.\n";

		}
		else{
			print "split/vanish step $svStep failed. Exiting\n";
			return $ret;
		}
	}


	if($gmtkOpt->{'outputTrainable'}){
		copy($stepGmtkOpt{'outputTrainable'}, $gmtkOpt->{'outputTrainable'}) || confess("copy failed");
	}
	if($gmtkOpt->{'llStoreFile'}){
		copy($stepGmtkOpt{'llStoreFile'}, $gmtkOpt->{'llStoreFile'}) || confess("copy failed");
	}

	print "Finished all split/vanish steps.\n";
	return 0;

}

## @fn
#split according to svOpts and converge to $convergenceThresholds
# 
#   @param distrOpt,gmtkOpt settings to use on each step
#	@param convergenceThresholds rules for when to stop iterating must have llRatioThreshold and maxEmIterations keys
#	@param svOpts split/vanish gmtk options that will be executed only on first iteration
# @return
#	same error codes as emtrainParallel
sub doSplitVanishStep{
	my ($gmtkOpt, $distrOpt,$convergenceThresholds,	$svOpts) = @_;

	if (!isTrainingNecessary($gmtkOpt, 1)){
		my @lines = readLines($gmtkOpt->{'llStoreFile'});
		my $curLL = $lines[0]; 
		chomp $curLL;
		print "log-likelihood for this step is: $curLL.\n";
		return 0;
	}
	
	print "iterating EM until one of two things happens: either llratio < $convergenceThresholds->{'llRatioThreshold'}";
	print " or number of iterations is > $convergenceThresholds->{'maxEmIterations'}.\n";

	my $it=0;
	print "Starting iteration $it.\n";

	my %stepGmtkOpt = %$gmtkOpt;
	my %stepDistrOpt = %$distrOpt;
	#first iter
	if($svOpts){
		print "Splitting/vanishing gaussians in iteration $it according to following options:\n";
		prettyPrintOpts($svOpts);
		@stepGmtkOpt{keys %$svOpts}= values %$svOpts;
	}
	$stepDistrOpt{'wd'} = "$distrOpt->{'wd'}/it$it";
	$stepDistrOpt{'jobName'}="$distrOpt->{'jobName'}_it$it";
	if($gmtkOpt->{'outputTrainable'}){
		$stepGmtkOpt{'outputTrainable'}="$distrOpt->{'wd'}/learnedParams_it$it.gmtk";
	}
	
	$stepGmtkOpt{'llStoreFile'}="$distrOpt->{'wd'}/logLikelihood_it$it.txt";

	my %emTrainArgs=(%stepGmtkOpt,%stepDistrOpt);
	my $ret = emtrainParallel(\%emTrainArgs);
	return $ret if($ret);
	my @lines = readLines("$distrOpt->{'wd'}/logLikelihood_it$it.txt");
	my $curLL = $lines[0]; 
	chomp $curLL;
	print "Total log-likelihood after iteration $it: $curLL.\n";
	if($svOpts){
		print "number of gaussians may have changed due to splits/vanishes:\n";
		print "ignoring log-likelihood of iteration $it for purposes of llRatioThreshold\n";
	}
	my $llratio = 1e999999999;

	while ($llratio > $convergenceThresholds->{'llRatioThreshold'} && $it+1<$convergenceThresholds->{'maxEmIterations'}){
		$it++;
		print "Starting iteration $it.\n";

		%stepGmtkOpt = %$gmtkOpt;
		%stepDistrOpt = %$distrOpt;
		
		$stepDistrOpt{'wd'} = "$distrOpt->{'wd'}/it$it";
		$stepDistrOpt{'jobName'}="$distrOpt->{'jobName'}_it$it";
		$stepDistrOpt{'doSanityCheck'}=0;
		$stepDistrOpt{'verbosity'}=$distrOpt->{'verbosity'}-1;

		$stepGmtkOpt{'inputTrainable'}="$distrOpt->{'wd'}/learnedParams_it".($it-1).".gmtk";
		$stepGmtkOpt{'outputTrainable'}="$distrOpt->{'wd'}/learnedParams_it$it.gmtk";
		$stepGmtkOpt{'llStoreFile'}="$distrOpt->{'wd'}/logLikelihood_it$it.txt";

		%emTrainArgs=(%stepGmtkOpt,%stepDistrOpt);
		$ret = emtrainParallel(\%emTrainArgs);
		return $ret if($ret);

		@lines = readLines("$distrOpt->{'wd'}/logLikelihood_it$it.txt");
		my $lastLL = $curLL;
		$curLL = $lines[0]; 
		chomp $curLL;
		print "total log-likelihood after iteration $it: $curLL\n";
		if ($it>2 || $it>1 && !$svOpts){
			$llratio = 1-$curLL/$lastLL;
			print "log-likelihood ratio is $llratio=1-(current log-likelihood/previous log-likelihood)\n";
		}
		else{
			$llratio = 1e999999999;
			print "Not enough iterations to compute log-likelihood ratio.  Iterating again.\n";
		}
	}

	if($gmtkOpt->{'outputTrainable'}){
		copy($stepGmtkOpt{'outputTrainable'}, $gmtkOpt->{'outputTrainable'}) || confess("copy failed");
	}
	if($gmtkOpt->{'llStoreFile'}){
		copy($stepGmtkOpt{'llStoreFile'}, $gmtkOpt->{'llStoreFile'}) || confess("copy failed");
	}

	print "finished iterating EM. Log-likelihood ratio is $llratio.  ".($it+1)." iterations performed. \n";

	return 0;
}

1;

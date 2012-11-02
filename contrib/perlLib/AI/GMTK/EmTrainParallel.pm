#!/usr/bin/env perl
package AI::GMTK::EmTrainParallel;
use warnings;
use strict;

use Cwd;
use File::Path;
use Carp;

use Parallel::Distribute qw(distribute);
use OS::Util qw(readLines writeLines tsSystem );
use Config::OptionsSet::OptionsSet qw(readOpts writeOpts prettyPrintOpts overwriteOpts 
					setDiff hashDiff hashSymDiff);
use AI::GMTK::Util  qw(	prepareGmtkCmdArgs checkOpts
					isRestartable notifyOfJobCompletion);

use Exporter;
our @ISA=qw(Exporter);
our @EXPORT_OK=qw(emtrainParallel isTrainingNecessary);

## @fn
#do em training in parallel.  If the execution is interrupted, it can 
#be safely rerun again and will start where it left of.
#See emtrainParallel.pl use help for details.
#
# @return
#	0  success
#	-1 cofiguration error. no files written or modified.
# 	-2 sanity check failed.  No parallel jobs launched.
#	-3 some parallel tasks failed.
sub emtrainParallel{
    my $opt = shift;
	
	$| =1; #flush pipes right away
	
    my ($err, $gmtkOpt, $distrOpt) = checkOpts($opt, \&initDistrOptsEmTrain, \&fillInOptsEmTrain);

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
        'verbosity',
        'emailAddress',
        'sendMailOnCompletion',
        'jobName');
	my $restartCode = isRestartable($gmtkOpt, $distrOpt, \@allowedDiffs);
	return 0 if (!isTrainingNecessary($gmtkOpt, $restartCode));

	print "single EM iteration started on ". localtime(time)."\n";


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
	mkpath([$distrOpt->{'_taskOutDir'},$distrOpt->{'_logDir'},$distrOpt->{'_scriptDir'},$distrOpt->{'_SGEDir'}], 1);
	
	if($restartCode == 0){	#first-time running
	
		#write out the config state so that the experiment is reproducible
		writeOpts($gmtkOpt, $distrOpt->{'_gmtkArgsBackup'},': ');
		writeOpts($distrOpt, $distrOpt->{'_distrArgsBackup'},' ');
	
	}
	
	
	#prepare task argfiles and command lines.
	prepareDistrEnv($gmtkOpt, $distrOpt);
	
	
	#do sanity check if requested and not previously performed
	if($distrOpt->{'doSanityCheck'}){	

		if(-e $distrOpt->{'_sanityCheckStatusFile'}){
			print "Sanity check succeeded previously. Not doing it again.\n";
		}
		else{
			print "Performing sanity check with command\n$distrOpt->{'_sanityCheckCommand'}\n";
			tsSystem ($distrOpt->{'_sanityCheckCommand'});
			if($?){
				print "Sanity check failed with exit code $?\n";
				return -2;
			}
			else{
				`touch $distrOpt->{'_sanityCheckStatusFile'}`;
				print "Sanity check succeeded\n";
			}
		}
	}
	
	#filter out successfully completed tasks, and write the SGE script that will be handled by distribute
	my @succeededTasks = succeededTasks($gmtkOpt, $distrOpt);
	my @remainingTasks = sort(setDiff([1..$distrOpt->{'numChunks'}],\@succeededTasks));
	print scalar(@succeededTasks)." tasks completed on previous run. ".scalar(@remainingTasks)." tasks remaining.\n";
	open (TASKS, ">$distrOpt->{'_scriptDir'}/taskCmds.txt") || confess("cannot open >$distrOpt->{'_scriptDir'}/taskCmds.txt");
	map {print TASKS "@{$distrOpt->{'_taskCommands'}}[$_-1]\n"} @remainingTasks;
	close TASKS;
	

	if(@remainingTasks){
		#run the jobs in parallel
		#($commands_file, $jobName, $sge) = @_;
		if(distribute("$distrOpt->{'_scriptDir'}/taskCmds.txt", , $distrOpt->{'jobName'}, $distrOpt->{'_SGEDir'},  $distrOpt->{'runLocally'}, $distrOpt->{'extraDistributeParameters'})){
			print "At least one of the parallel jobs did not finish successfully. \n";
			return -3;
		}
		else{
			print "All remaining tasks finished. \n";
		}
	}
	
	#now estimate the new model
	my $errCode = postParallelStep($gmtkOpt, $distrOpt);
	if ($errCode){
		return -2;
	}


	#clear accumulators at successful end
	#they are also cleared if clearWdOnEnd is true
	if($distrOpt->{'clearAccOnEnd'}){	
		print "deleting accumulator dir $distrOpt->{'_taskOutDir'} per clearAccOnEnd flag\n";
		rmtree($distrOpt->{'_taskOutDir'});
	}

	#clear workingDir at successful end
	if($distrOpt->{'clearWdOnEnd'}){	
		print "deleting working directory $distrOpt->{'wd'} per clearWdOnEnd flag\n";
		rmtree($distrOpt->{'wd'});
		
	}


	print "single EM iteration succeeded on ". localtime(time)."\n";
	return 0;
}

## @fn
#guess some defaults for distributed environment options.
sub initDistrOptsEmTrain{
	my $opt = shift;
	

	#these MUST be set - no defaults for them

	#the working Dir path
	$opt->{'wd'} = undef;


	#these CAN be set and have defaults

	#paths of needed programs
	my $buf;
	$buf = `which qsub 2> /dev/null`;
	chomp $buf;
	$opt->{'qsubPath'}= $?==0 ? $buf : "";
	
	$buf = `which breakRange.pl 2> /dev/null`;
	chomp $buf;
	$opt->{'breakRangePath'}= $?==0 ? $buf : "";
	
	#The actual executable for gmtkEMtrainNew
	$buf = `which gmtkEMtrainNew 2> /dev/null`;
	chomp $buf;
	$opt->{'gmtkEMtrainNewPath'}= $?==0 ? $buf : "";

	#the name of the job as it will appear in the SGE queue monitoring tools
    #the name is not used anywhere else
	$opt->{'jobName'}='emtrain';

	#extra paramaters to pass to the batch queueing command (qsub, sbatch, etc...). 
	#The parameters are just passed into the distribute function without being interpreted
	$opt->{'extraDistributeParameters'} = '';


	#the following probably will be common to many parallel scripts
 
	#perform a sanityCheck before firing off a ton of parallel jobs
	$opt->{'doSanityCheck'}=1;

	#consider job restartable even if config files differ
	$opt->{'overrideRestartable'}=0;

	#number of chunks.  Will be set in fillInOpts
	$opt->{'numChunks'}='';

	#actual utterance Ranges.  Will be set in fillInOpts
	$opt->{'uttRanges'}='';
	
	#clear workingDir at start.  All the compute work will be repeated.
	$opt->{'clearWdOnStart'}=0;

	#clear workingDir at successful end
	$opt->{'clearWdOnEnd'}=0;

	#clear accumulators at successful end
	#they are also cleared if clearWdOnEnd is true
	$opt->{'clearAccOnEnd'}=0;


	$opt->{'sendMailOnCompletion'}=0;
	$opt->{'emailAddress'}='';

	#greater than 0 means print detailed config info. 
	#The higher the verbosity the more info is printed.
	#some messages are printed even at verbosity=0 
	$opt->{'verbosity'}=1;

	#if true, don't submit chunks to cluster but run them localy sequentialy, via bash
	$opt->{'runLocally'}='0';

}

sub postParallelStep{
	my $gmtkOpt = shift;
	my $distrOpt = shift;

	$gmtkOpt->{'loadAccFile'}= $distrOpt->{'_taskOutDir'}.'/acc_file_@D.data';
	$gmtkOpt->{'loadAccRange'}= '1:'.$distrOpt->{'numChunks'};
	$gmtkOpt->{'trrng'}= 'nil';

	my $emEstArgs=prepareGmtkCmdArgs($gmtkOpt,"$distrOpt->{'_scriptDir'}/emEst.args");
	my $emEstCmd= "$distrOpt->{'gmtkEMtrainNewPath'} $emEstArgs > $distrOpt->{'_logDir'}/emEst.log 2>&1 ";
	print "performing the EM estimation step with command $emEstCmd.\n";
	tsSystem($emEstCmd);
 	if ($?){
		print "EM estimation step failed with code $?\n";
	}
	
	if($gmtkOpt->{'storeAccFile'} && $gmtkOpt->{'outputTrainable'}){
		print "estimation step completed for storeAccFile $gmtkOpt->{'storeAccFile'}.  Now doing EM step for outputTrainable $gmtkOpt->{'outputTrainable'}\n";
		$gmtkOpt->{'loadAccFile'}= $gmtkOpt->{'storeAccFile'};
		delete $gmtkOpt->{'storeAccFile'};
		delete $gmtkOpt->{'loadAccRange'};
		$emEstArgs=prepareGmtkCmdArgs($gmtkOpt,"$distrOpt->{'_scriptDir'}/emEst2.args");
		$emEstCmd= "$distrOpt->{'gmtkEMtrainNewPath'} $emEstArgs > $distrOpt->{'_logDir'}/emEst2.log 2>&1 ";
		print "performing the EM estimation step with command $emEstCmd.\n";
		tsSystem($emEstCmd);
		if ($?){
			print "EM estimation step failed with code $?\n";
		}
	}
	return $?;
}

## @fn
#returns a list of tasknumbers that have succeeded
sub succeededTasks{
	my $gmtkOpt = shift;
	my $distrOpt = shift;

	my @successlines=();
	foreach my $i (1..$distrOpt->{'numChunks'}){
		system("(tail -n 5 $distrOpt->{'_logDir'}/emacc.$i.log | grep -H 'PROGRAM\ ENDED\ SUCCESSFULLY WITH STATUS 0' )> /dev/null 2>&1");
		confess("@?") if @?;
		my $status = $? >> 8;
		#print "$i: $status $distrOpt->{'_logDir'}/emacc.$i.log\n";
		push @successlines, $i if $status ==0;
	}
	@successlines = grep {-e "$distrOpt->{'_taskOutDir'}/acc_file_$_.data"} @successlines;
	return (@successlines);
}



sub fillInOptsEmTrain{

	my $gmtkOpt = shift;
	my $distrOpt = shift;

	my $ret=0;

	#check options
	if (!$distrOpt->{'wd'}){
		print "working directory must be specified\n";
		$ret++;
	}

	if (!$distrOpt->{'runLocally'} && (!$distrOpt->{'qsubPath'} || ! -x $distrOpt->{'qsubPath'})){
		print "cannot find qsub at path $distrOpt->{'qsubPath'}\n";
		$ret++;
	}
	if (!$distrOpt->{'gmtkEMtrainNewPath'} || ! -x $distrOpt->{'gmtkEMtrainNewPath'}){
		print "cannot find gmtkEMtrainNew at path $distrOpt->{'gmtkEMtrainNewPath'}\n" ;
		$ret++;
	}
	
	
	if (! $distrOpt->{'breakRangePath'} || ! -x $distrOpt->{'breakRangePath'}){
		print "cannot find breakRange.pl at path $distrOpt->{'breakRangePath'}\n";
		$ret++;
	}

	if ($distrOpt->{'sendMailOnCompletion'} && ! $distrOpt->{'emailAddress'}){
		print "emailAddress param is required if sendMailOnCompletion is true\n";
		$ret++;
	}

	if (!$gmtkOpt->{'storeAccFile'} && !$gmtkOpt->{'outputTrainable'}){
		print "outputTrainable or storeAccFile is required in gmtk args file.\n";
		$ret++;
	}

	if (!$gmtkOpt->{'of1'}){
		print "of1 is required in gmtk args file\n";
		$ret++;
	}

	return $ret if ($ret);

	#fill in unspecified values with reasonable defaults

	#let the number of chunks be as many as there SGE queues
	if(!$distrOpt->{'numChunks'}){
		my @queueCounts= `qhost  -F num_proc | grep num_proc `;
		$distrOpt->{'numChunks'}=0;
		foreach (@queueCounts){
			chomp;
			s/^.*=//;
			$distrOpt->{'numChunks'} += $_;
		}

		print "numChunks defaults to the number of SGE queues ($distrOpt->{'numChunks'}).\n";

	}

	#training range is the entire training file
	if(!$gmtkOpt->{'trrng'}){
		#my $lastIndex=`wc -l $gmtkOpt->{'of1'}`;
		#(undef,$lastIndex, undef)=split(',', $lastIndex, 2);
		my $cmd="obs-info -ifmt1 $gmtkOpt->{'fmt1'} -i1 $gmtkOpt->{'of1'}";
		if(defined($gmtkOpt->{'nf1'})){
			$cmd .= " -nf1 $gmtkOpt->{'nf1'}";
		}
		if(defined($gmtkOpt->{'ni1'})){
			$cmd .= " -ni1 $gmtkOpt->{'ni1'}";
		}
		my $lastIndex=`$cmd`;
		$lastIndex =~ /([0-9]+) sentences/;
		$lastIndex=$1;
		
		$lastIndex--;
		$gmtkOpt->{'trrng'}="0:$lastIndex";
		print "trrng defaults to all the utterances in the first observation file ($gmtkOpt->{'trrng'}).\n";
	}

	#determine the number of utterances
	my ($start,$end) = split(/:/,$gmtkOpt->{'trrng'});
	$end = $start if (!$end);
	my $items = $end-$start+1;

	if ($items < $distrOpt->{'numChunks'}){
		$distrOpt->{'numChunks'}=$items;
		print "There are only $items utterances. Decreasing numChunks to $distrOpt->{'numChunks'}.\n";
	} 

	#calc the utterance ranges for each chunk
	if(!$distrOpt->{'uttRanges'}){
		$distrOpt->{'uttRanges'} = `breakRange.pl -r $gmtkOpt->{'trrng'} -n $distrOpt->{'numChunks'} 2>&1`;
		if ($?){
			confess("Error executing breakRange.pl\n");
		}
		chomp $distrOpt->{'uttRanges'};
		$distrOpt->{'uttRanges'} =~ s/\s+$//;

	}
	else{ #if uttRanges were specified manually, set the numChunks to match 
		my @rangesStr = split ' ', $distrOpt->{'uttRanges'};
		if (@rangesStr>$items){
			print "there are more ranges in uttRanges (".scalar(@rangesStr).") than there are items ($items).\n";
			$ret++;
		}
		else{ 
			$distrOpt->{'numChunks'} = @rangesStr;
			print "numchunks set to match uttRanges $distrOpt->{'numChunks'} chunks generated.\n";
		}

	}

	#filenames for saving state
	$distrOpt->{'_gmtkArgsBackup'} = "$distrOpt->{'wd'}/gmtk.args";
	$distrOpt->{'_distrArgsBackup'} = "$distrOpt->{'wd'}/distr.args";

	#other paths and filenames 
	$distrOpt->{'_taskOutDir'} = "$distrOpt->{'wd'}/taskOutputs";
	$distrOpt->{'_logDir'} = "$distrOpt->{'wd'}/logs";
	$distrOpt->{'_scriptDir'} = "$distrOpt->{'wd'}/scripts";
	$distrOpt->{'_SGEDir'} = "$distrOpt->{'wd'}/SGE";
	$distrOpt->{'_sanityCheckStatusFile'} = "$distrOpt->{'wd'}/sanityCheckSucceeded";

	return $ret; 
}

## @fn 
#  @return 0 if training is not necessary, (all requested files already exist), or 1 if it is necessary
sub isTrainingNecessary{
	my ($gmtkOpt, $restartCode) = @_;
	my $retrain = 0;
	if($gmtkOpt->{'outputTrainable'}){
		if(-e "$gmtkOpt->{'outputTrainable'}"){
			print "trainable parameters file $gmtkOpt->{'outputTrainable'} already exists.\n";
		}
		else{
			$retrain =1;
		}
	}
	if($gmtkOpt->{'storeAccFile'}){
		if(-e "$gmtkOpt->{'storeAccFile'}"){
			print "accumulator file $gmtkOpt->{'storeAccFile'} already exists.\n";
		}
		else{
			$retrain =1;
		}
	}
	if(!$retrain){
		print "Not performing training.\n";
        if ($restartCode==2){
            print "WARNING: The config has changed.  The parameters file $gmtkOpt->{'outputTrainable'} or acc file $gmtkOpt->{'storeAccFile'} may not be generated from current settings.\n";
        }
	}
	return $retrain;
}

## @fn
#create the working dir, sge command lines, individual task files etc,
#Unless clearWdOnStart is true, this does not overwrite partially completed files
sub prepareDistrEnv{
	my $gmtkOpt = shift;
	my $distrOpt = shift;

	#construct the sanity check arg file and command-line args
	my %taskParams = %$gmtkOpt;
    $taskParams{'trrng'} =~ s/:.*//; #use only the first utterance from the range
	#$taskParams{'ckbeam'}=10;
	$taskParams{'outputTrainable'}="$distrOpt->{'wd'}/sanityCheck.gmtk_";
	my $sanityCheckArgs = prepareGmtkCmdArgs(\%taskParams, "$distrOpt->{'_scriptDir'}/taskSanityCheck.args");
	$distrOpt->{'_sanityCheckCommand'} = "$distrOpt->{'gmtkEMtrainNewPath'} $sanityCheckArgs";


	#construct the arg files for individual tasks
	my @ranges = split ' ', $distrOpt->{'uttRanges'};
	for (1..@ranges){
		my %taskParams = %$gmtkOpt;
		$taskParams{'trrng'}=$ranges[$_-1];
		$taskParams{'storeAccFile'}="$distrOpt->{'_taskOutDir'}/acc_file_$_.data";
		delete $taskParams{'outputTrainable'};
		my $taskArgs = prepareGmtkCmdArgs(\%taskParams, "$distrOpt->{'_scriptDir'}/task$_.args");
		my $taskCmd = "$distrOpt->{'gmtkEMtrainNewPath'} $taskArgs  > $distrOpt->{'_logDir'}/emacc.$_.log 2>&1";
		push(@{$distrOpt->{'_taskCommands'}}, $taskCmd);
	}


}

1;

#!/usr/bin/env perl
package AI::GMTK::ViterbiParallel;
use warnings;
use strict;

use File::Path;
use Carp;
use Data::Dumper;
use Parallel::Distribute qw(distribute);
use OS::Util qw(readLines writeLines tsSystem );
use Config::OptionsSet::OptionsSet qw(readOpts writeOpts prettyPrintOpts overwriteOpts 
					setDiff hashDiff hashSymDiff getNamespaceOpts);
use AI::GMTK::Util  qw(	prepareGmtkCmdArgs checkOpts
					isRestartable notifyOfJobCompletion);
use Exporter;
our @ISA=qw(Exporter);
our @EXPORT_OK=qw(viterbiParallel checkAccuracy viterbiAndAccuracy);


## @fn
#do viterbi decoding in parallel.  If the execution is interrupted, it can 
#be safely rerun again and will start where it left of.
#See viterbiParallel.pl use help for details.
#
# @return
#	0  success
#	-1 cofiguration error. no files written or modified.
# 	-2 sanity check failed.  No parallel jobs launched.
#	-3 some parallel tasks failed.
sub viterbiParallel{
    my $opt = shift;
	
	$| =1; #flush pipes right away
	
    my ($err, $gmtkOpt, $distrOpt) = checkOpts($opt, \&initDistrOptsViterbi, \&fillInViterbiOpts);

	if ($err){
		return $err;
	}
	

	#check to see if the job is restartable.  A job is restartable if the 
	#distrOpt and gmtkOpt exist in the working directory and are identical to
	#the distrOpt and gtmkOpt we've set up so far
	my @allowedDiffs=(
        'clearWdOnStart',
        'clearWdOnEnd',
        'breakRangePath',
        'qsubPath',
        'decoderOutFile',
        'verbosity',
        'emailAddress',
        'sendMailOnCompletion',
        'jobName');
	my $restartCode = isRestartable($gmtkOpt, $distrOpt, \@allowedDiffs);
    
	if(-e "$distrOpt->{'decoderOutFile'}"){
		my $lastLine = `tail -n 1 $distrOpt->{'decoderOutFile'}`;
		if ($lastLine =~ /^Total data log prob/){
			print "decoder output file $distrOpt->{'decoderOutFile'} already exists.\n";
			print "Not performing decoding.\n";
            if ($restartCode==2){
		        print "WARNING: The config has changed.  The decoder output file $distrOpt->{'decoderOutFile'} may not be generated from current settings.\n";
            }
			return 0;
		}
		else{
			print "decoder output file $distrOpt->{'decoderOutFile'} already exists, but seems corrupt.\n";
			print "Performing decoding again.\n";
		}
	}

	print "viterbi decoding started on ". localtime(time)."\n";


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
		if(distribute("$distrOpt->{'_scriptDir'}/taskCmds.txt", $distrOpt->{'jobName'}, $distrOpt->{'_SGEDir'},  $distrOpt->{'runLocally'},  $distrOpt->{'extraDistributeParameters'})){
			print "At least one of the parallel jobs did not finish successfully. \n";
			return -3;
		}
		else{
			print "All remaining tasks finished. \n";
		}
	}
	
	my $errCode = postParallelStep($gmtkOpt, $distrOpt);
	if ($errCode){
		return -2;
	}


	#clear workingDir at successful end
	if($distrOpt->{'clearWdOnEnd'}){	
		print "deleting working directory $distrOpt->{'wd'} per clearWdOnEnd flag\n";
		rmtree($distrOpt->{'wd'});
		
	}

	print "viterbi decoding job $distrOpt->{'jobName'} succeeded on ". localtime(time)."\n";

	return 0;
}


sub postParallelStep{
	my $gmtkOpt = shift;
	my $distrOpt = shift;


	print "concatinating the task outputs together in sequence\n";
	open(OUT,">$distrOpt->{'decoderOutFile'}") || confess("cannot open >$distrOpt->{'decoderOutFile'}");
	foreach(1..$distrOpt->{'numChunks'}){
		my @chunkLines =readLines("$distrOpt->{'_taskOutDir'}/vit_$_.decode");
		#for some reason it looks like STDERR goes to STDOUT
		@chunkLines = grep(!/Junction Tree/, @chunkLines);
		@chunkLines = grep(!/### Final time/, @chunkLines);
		@chunkLines = grep(!/PROGRAM ENDED SUCCESSFULLY/, @chunkLines);
		print OUT @chunkLines;
	}
	close(OUT) || confess("cannot close >$distrOpt->{'decoderOutFile'}");


 	if ($?){
		print "EM estimation step failed with code $?\n";
	}
	return $?;
}

#returns a list of tasknumbers that have succeeded
sub succeededTasks{
	my $gmtkOpt = shift;
	my $distrOpt = shift;

	my @successlines=();
	foreach my $i (1..$distrOpt->{'numChunks'}){
		system("(tail -n 5 $distrOpt->{'_taskOutDir'}/vit_$i.decode | grep -H 'PROGRAM\ ENDED\ SUCCESSFULLY WITH STATUS 0' )> /dev/null 2>&1");
		confess("@?") if @?;
		my $status = $? >> 8;
		push @successlines, $i if $status ==0;
	}

	return (@successlines);
}



sub fillInViterbiOpts{

	my $gmtkOpt = shift;
	my $distrOpt = shift;

	my $ret=0;

	#check options
	if (!$distrOpt->{'runLocally'} && (!$distrOpt->{'qsubPath'} || ! -x $distrOpt->{'qsubPath'})){
		print "cannot find qsub at path $distrOpt->{'qsubPath'}\n";
		$ret++;
	}
	if (!$distrOpt->{'gmtkViterbiPath'} || ! -x $distrOpt->{'gmtkViterbiPath'}){
		print "cannot find gmtkViterbiNew at path $distrOpt->{'gmtkViterbiPath'}\n" ;
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

	if (!$distrOpt->{'decoderOutFile'}){
		print "decoderOutFile is required in distr args file\n";
		$ret++;
	}

	if (!$distrOpt->{'wd'}){
		print "wd is required in distr args file\n";
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
	if(!$gmtkOpt->{'dcdrng'}){

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
		$gmtkOpt->{'dcdrng'}="0:$lastIndex";
		print "dcdrng defaults to all the utterances in the first obeservation file ($gmtkOpt->{'dcdrng'}).\n";
	}

	#determine the number of utterances
	my ($start,$end) = split(/:/,$gmtkOpt->{'dcdrng'});
	$end = $start if (!$end);
	my $items = $end-$start+1;

	if ($items < $distrOpt->{'numChunks'}){
		$distrOpt->{'numChunks'}=$items;
		print "There are only $items utterances. Decreasing numChunks to $distrOpt->{'numChunks'}.\n";
	} 

	#calc the utterance ranges for each chunk
	if(!$distrOpt->{'uttRanges'}){
		$distrOpt->{'uttRanges'} = `breakRange.pl -r $gmtkOpt->{'dcdrng'} -n $distrOpt->{'numChunks'} 2>&1`;
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
#create the working dir, sge command lines, individual task files etc,
#Unless clearWdOnStart is true, this does not overwrite partially completed files
sub prepareDistrEnv{
	my $gmtkOpt = shift;
	my $distrOpt = shift;

	#construct the sanity check arg file and command-line args
	my %taskParams = %$gmtkOpt;
	$taskParams{'dcdrng'} =~ s/:.*//; #use only the first utterance from the range
	$taskParams{'ckbeam'}=10;
	my $sanityCheckArgs = prepareGmtkCmdArgs(\%taskParams, "$distrOpt->{'_scriptDir'}/taskSanityCheck.args");
	$distrOpt->{'_sanityCheckCommand'} = "$distrOpt->{'gmtkViterbiPath'} $sanityCheckArgs ";



	#construct the arg files for individual tasks
	my @ranges = split ' ', $distrOpt->{'uttRanges'};
	for (1..@ranges){
		my %taskParams = %$gmtkOpt;
		$taskParams{'dcdrng'}=$ranges[$_-1];
		
		my $taskArgs = prepareGmtkCmdArgs(\%taskParams, "$distrOpt->{'_scriptDir'}/task$_.args");
		my $taskCmd = "$distrOpt->{'gmtkViterbiPath'} $taskArgs  > $distrOpt->{'_taskOutDir'}/vit_$_.decode 2>$distrOpt->{'_logDir'}/vit_$_.log";
		push(@{$distrOpt->{'_taskCommands'}}, $taskCmd);
	}


}

## @fn
#guess some defaults for distributed environment options.
sub initDistrOptsViterbi{
	my $opt = shift;
	
	#these MUST be set - no defaults for them

	#the working Dir path
	$opt->{'wd'} = undef;

	#the decoder output file
	$opt->{'decoderOutFile'}=undef;


	#these CAN be set and have defaults


	#paths of needed programs
	my $buf;
	$buf = `which qsub 2> /dev/null`;
	chomp $buf;
	$opt->{'qsubPath'}= $?==0 ? $buf : '';
	
	$buf = `which breakRange.pl 2> /dev/null`;
	chomp $buf;
	$opt->{'breakRangePath'}= $?==0 ? $buf : '';
	
	#The actual executable for gmtkViterbiNew
	$buf = `which gmtkViterbiNew 2> /dev/null`;
	chomp $buf;
	$opt->{'gmtkViterbiPath'}= $?==0 ? $buf : '';

	#the name of the job as it will appear in the SGE queue monitoring tools
    #the name is not used anywhere else
	$opt->{'jobName'}='viterbi';

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


	$opt->{'sendMailOnCompletion'}=0;
	$opt->{'emailAddress'}='';

	#greater than 0 means print detailed config info. 
	#The higher the verbosity the more info is printed.
	#some messages are printed even at verbosity=0 
	$opt->{'verbosity'}=1;

	#if true, don't submit chunks to cluster but run them localy sequentialy, via bash
	$opt->{'runLocally'}='0';

}


## @fn
#score acuracy of the decoding for data in switchboard format
#@pre the $opt (listed below) are valid, and the necessary files exist 
#runs sclite $distrOpt->{'decoderOutFile'} and stores results in
#$opt->{'wd'}/accuracy
#
#the transcriptions to be tested
#$opt->{'decoderOutFile'}
#
#the working dir.  (An accuracy subdirectory will be created under it) 
#$opt->{'wd'}
#
#The reference transcription file in sclite format.
#A perfect transcription will look identical to this reference transcription in this file
#$opt->{'referenceTrnFile'}
#
#The sentence ids file.  It is a list of  utterance ids (in switchboard format)
#where the ith sentence id (on the ith line in that file) is describing the ith utterance 
#in the decoderOutFile.
#$opt->{'sentIdsFile'}
#
#@return
# a list with first element being the error code (0 success -4 error)
# second element being the Word Error Rate
sub checkAccuracy{
	my $inOpt = shift;
	
	
    my ($err, $dummy, $opt) = checkOpts($inOpt, \&initOptsAccuracy, \&fillOptsAccuracy);
	prettyPrintOpts($opt) if ($opt->{'verbosity'} && $opt->{'verbosity'}  >=1);

	if ($err){
		return ($err,-1);
	}

	mkpath($opt->{'_accDir'},1);
    	my $cmd = "gmtk_output_to_trn_file.pl $opt->{'decoderOutFile'} $opt->{'sentIdsFile'} > $opt->{'_accDir'}/preOut.nosil.trn";
	tsSystem($cmd);
	
	#load words that should be deleted.
	my %delVocab;
	if($opt->{'deleteVocabFile'}){
		my @lines = readLines($opt->{'deleteVocabFile'});
		chomp @lines;
		@lines = grep {!/^#/} @lines;
		@delVocab{@lines}=@lines;
		print "read ".scalar(keys %delVocab)." deletable words from deleteVocabFile $opt->{'deleteVocabFile'}\n";
	}
	my %multWords=();
	if($opt->{'multiwordFile'}){
		my @lines = readLines($opt->{'multiwordFile'});
		chomp @lines;
		@lines = grep {!/^#/} @lines;
		foreach (@lines){
			my ($k, $v) =split(' ', $_, 2);
			$multWords{$k}=$v;
		}
		print "read ".scalar(keys %multWords)." multiwords from multiwordFile $opt->{'multiwordFile'}\n";
		#print Dumper(\%multWords);
	}
	
	my @lines = readLines("$opt->{'_accDir'}/preOut.nosil.trn");
	foreach(@lines){
		chomp;
		s/<s>//ig;
		s/<\/s>//ig;
		my @words = split;
		my @delWords=map {$delVocab{$_} ? '%HESITATION' : $_} @words;
		my @mWords=map {$multWords{$_} ? $multWords{$_} : $_} @delWords;
		#my @mWords= @delWords;
		$_ = join (" ", @mWords);
		#$_ = join (" ", grep {!/<UNK>/} split);
		$_ = $_."\n";
	}
	writeLines("$opt->{'_accDir'}/unfilteredOut.nosil.trn", \@lines); 

	if($opt->{'glmFile'}){
		$cmd = "csrfilt.sh -dh -t hyp $opt->{'glmFile'} < $opt->{'_accDir'}/unfilteredOut.nosil.trn  > $opt->{'_accDir'}/out.nosil.trn";
		print "$cmd\n";
		tsSystem($cmd);
	}
	else{
		tsSystem("cp -v $opt->{'_accDir'}/unfilteredOut.nosil.trn $opt->{'_accDir'}/out.nosil.trn ");
	}
	
   	$cmd = "sclite $opt->{'extraScliteParams'} -h $opt->{'_accDir'}/out.nosil.trn -i swb -o dtl -r $opt->{'referenceTrnFile'} > $opt->{'_accDir'}/out.sclite\n";
	print "$cmd\n";
	tsSystem($cmd);
	
	$cmd= "sclite $opt->{'extraScliteParams'} -h $opt->{'_accDir'}/out.nosil.trn -i swb -o all -o stdout -r $opt->{'referenceTrnFile'} > $opt->{'_accDir'}/out.sclite.snt";
	print "$cmd\n";
	tsSystem($cmd);

	#this is for the sc_stats tool
	$cmd= "sclite $opt->{'extraScliteParams'} -h $opt->{'_accDir'}/out.nosil.trn -i swb -o sgml -o stdout -r $opt->{'referenceTrnFile'} > $opt->{'_accDir'}/out.sclite.sgml";
	print "$cmd\n";
	tsSystem($cmd);

	@lines = readLines("$opt->{'_accDir'}/out.nosil.trn.dtl");
	@lines = grep(/Percent Total Error/,@lines);
	my $wer = $lines[0];
	chomp $wer;
	$wer =~ s/.*=[ ]*//;
	$wer =~ s/\%.*//;
	if ($wer !~ /[.-0-9]+/){
		print "error computing accuracy: WER=$wer. Check $opt->{'_accDir'}/out.nosil.trn.dtl\n";
		return -4;
	}
	return (0, $wer);
}

sub initOptsAccuracy{
	my $opt = shift;
	
	#these MUST be set - no defaults for them
	#the working Dir path
	$opt->{'wd'} = undef;

	#the decoder output file
	$opt->{'decoderOutFile'}=undef;

    #if non empty, calculate WER against the reference transcription file
	$opt->{'referenceTrnFile'}=undef;
    
    #the sentence Ids file - must be set if referenceTrnFile is not empty
    $opt->{'sentIdsFile'}=undef;

	#greater than 0 means print detailed config info. 
	#The higher the verbosity the more info is printed.
	#some messages are printed even at verbosity=0 
	$opt->{'verbosity'}=1;

	#extra parameters to be passed to sclite, e.g. -D -F
	#by default, don't penalize optionally deletable words, if they are specified in the reference transcription
	$opt->{'extraScliteParams'}='-D';
	
	#if specified, go through the transcription and delete the words in this file.
	#good for deleting non-speech and filled pauses
	$opt->{'deleteVocabFile'}='';
	
	#if specified, replace the multiword in the transcription with the multiple words that define it
	#each line of the file has form: multiword word1 word2 word3 ... 
	$opt->{'multiwordFile'}='';

	#if specified, use it in tranfilt to filter the hypothesis trn file in hyp mode
	$opt->{'glmFile'}='';
}

sub fillOptsAccuracy{

	my $leftoverOpt = shift;
	my $opt = shift;

	my $ret=0;

	#check options
	if (!$opt->{'wd'}){
		print "wd is required in distr args file\n";
		$ret++;
	}
	
	if ($opt->{'referenceTrnFile'} && ! -r $opt->{'referenceTrnFile'}){
		print "referenceTrnFile at path $opt->{'referenceTrnFile'} is not readable\n";
		$ret++;
	}
    
	if ($opt->{'deleteVocabFile'} && ! -r $opt->{'deleteVocabFile'}){
		print "deleteVocabFile at path $opt->{'deleteVocabFile'} is not readable\n";
		$ret++;
	}
	
	if ($opt->{'multiwordFile'} && ! -r $opt->{'multiwordFile'}){
		print "multiwordFile at path $opt->{'multiwordFile'} is not readable\n";
		$ret++;
	}
	
	if ($opt->{'sentIdsFile'} && ! -r $opt->{'sentIdsFile'}){
		print "sentIdsFile at path $opt->{'sentIdsFile'} is not readable\n";
		$ret++;
	}
	
    if ($opt->{'sentIdsFile'} && ! $opt->{'referenceTrnFile'}){
		print "referenceTrnFile param is required if sentIdsFile is specified\n";
		$ret++;
	}

	if (!$opt->{'decoderOutFile'}){
		print "decoderOutFile is required in distr args file\n";
		$ret++;
	}
	
	if (keys %$leftoverOpt){
		my @badkeys = keys %$leftoverOpt;
		print "The following accuracy arguments are not recognized: @badkeys\n";
	}
	
	$opt->{'_accDir'} = "$opt->{'wd'}/accuracy" if (!$opt->{'_accDir'});
	
	return $ret;
}

## @fn
#does viterbi and check accuracy, 
# @return ($errCode, $wer)
# @deprecated just use viterbiParallel followed by checkAccuracy. It's not that hard.
sub viterbiAndAccuracy{
	my ($allOpt)= @_;

	#we want any configuration errors to come back right away, so check accuracy stuff
	my $accOpts = $allOpt->{'accuracy'};
	if($accOpts->{'referenceTrnFile'}){
    	(my $err, my $dummy, $accOpts) = checkOpts($accOpts, \&initOptsAccuracy, \&fillOptsAccuracy);
		if ($err){
			return ($err,-1);
		}
	}

	my $vitOpts=$allOpt;
    delete $vitOpts->{'accuracy'};
	my $ret = viterbiParallel($vitOpts);
	my $wer;
	if ($ret ==0){
		if($accOpts->{'referenceTrnFile'}){
			print "checking accuracy against a reference transcription\n";
    		($ret, $wer) = checkAccuracy($accOpts);

			if($ret==0){
				print "Word Error rate: $wer\n";
			}
			else{
				print "checking accuracy failed\n";
			}
		}
		else{
			print "referenceTrnFile not specified.  Not checking transcription accuracy.\n";
		}
	}
	return ($ret, $wer);
}
1;

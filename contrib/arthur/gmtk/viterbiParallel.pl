#!/usr/bin/env perl
use warnings;
use strict;

## @file 
# 
#  @author Arthur Kantor 

use Getopt::Long;
use AI::GMTK::ViterbiParallel qw(viterbiParallel);
use Config::OptionsSet::OptionsSet qw(readOpts writeOpts prettyPrintOpts overwriteOpts 
					setDiff hashDiff hashSymDiff);
use AI::GMTK::Util  qw(	prepareGmtkCmdArgs checkOpts
					isRestartable notifyOfJobCompletion);

# use Gmtk::gmtkUtilParallel  qw(	prepareGmtkCmdArgs readOpts writeOpts prettyPrintOpts overwriteOpts 
# 					setDiff hashDiff hashSymDiff
# 					readLines writeLines
# 					isRestartable notifyOfJobCompletion
# 					tsSystem );

my $usage = <<EOS;
Usage: viterbiParallel.pl [-h] [-d <distributedConfigFile>] <gmtkArgFile> <workingDir> <decoderOutFile>

Runs gmtkViterbiNew in parallel via SGE.

Required:
<gmtkArgFile>  
	a file with arguments to gmtkViterbiNew.  It must be in the 
	same format as a file used by gmtkViterbiNew -gmtkArgFile flag,
	with additional requirement of one flag per line.
	
<workingDir> 
	a dir for all temprary files.

<decoderOutFile>
	A file where to store the final decoded transcriptions of the utterances.  
	Equivalent to redirecting  the non-parallel gmtkViterbiNew's STDOUT to this file
	
Options:
-d --distributedConfig 
	<distributedConfigFile>  a file containing distributed 
	settings in format of one VAR VALUE pair per line.
	To see the possible options just run viterbiParallel.pl without the -d 
	option and all legal settings will be printed.
	Note that any setting starting with _ cannot be set in 
	the distributedConfigFile on the command line, and 
	will be overwritten.

-h --help 
	print this help message

Return values:
	0  success
	-1 cofiguration error. no files written or modified.
 	-2 sanity check failed.  No parallel jobs launched.
	-3 some parallel tasks failed.
	-4 transcription accuracy checking failed.

This script is safe to rerun if it is interrupted by a user, cluster failure, etc.
It will safely restart, repeating minimum work necessary to finish.
A dry run is performed with a narrow beam on one utterance before 
sending out the job in parallel to make sure emtraining will at least start.

EOS

#parse the arguments
{#brackets to keep my variables private from subroutines 
	#read command-line options
	my $distrConfigFile;
	my $gmtkArgFile;
	my $workingDir;
	my $decoderOutFile;
	my $help;
	
	GetOptions (
				"distributedConfig|d=s" => \$distrConfigFile,
				"help|h" => \$help);
	

	local $SIG{INT} = \&killChildrenHandler; 
	local $SIG{QUIT} = \&killChildrenHandler;
	local $SIG{TERM} = \&killChildrenHandler;
	
	my $commands_file;
	if (@ARGV < 3 || $help) {
		die("$usage");
	} else {
		$gmtkArgFile = $ARGV[0];
		$workingDir = $ARGV[1];
		$decoderOutFile =$ARGV[2];
	}
	
	
	#read gmtk and distributed options
	my %gmtkOpt;
	readOpts($gmtkArgFile, \%gmtkOpt);

	my %userDistrOpt=();
	$userDistrOpt{'wd'} = $workingDir;
	$userDistrOpt{'decoderOutFile'} = $decoderOutFile;
	if($distrConfigFile){
		readOpts($distrConfigFile, \%userDistrOpt);
	}

	
	# %userDistrOpt can contains a mixture of options for 
	#checking accuracy and for the actual viterbi.
	#seperate them out.
	my @accOnlyKeys = ('referenceTrnFile','sentIdsFile');
	my @accKeys=('wd','decoderOutFile', @accOnlyKeys);
	my %accOpts;
	@accOpts{@accKeys}=@userDistrOpt{@accKeys};
	delete @userDistrOpt{@accOnlyKeys};
	prettyPrintOpts(\%userDistrOpt);
	my $ret = viterbiParallel(\%gmtkOpt,\%userDistrOpt);

	if ($ret ==0){
		if($accOpts{'referenceTrnFile'}){
			print "checking accuracy against a reference transcription\n";
			($ret, my $wer) = checkAccuracy(\%accOpts);
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
	notifyOfJobCompletion(\%userDistrOpt, $ret);

	print "Exiting with code $ret.\n";
	exit ($ret);
}



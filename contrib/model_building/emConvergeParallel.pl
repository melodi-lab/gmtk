#!/usr/bin/env perl
## @file 
# 
#  @author Arthur Kantor 
use warnings;
use strict;

use Getopt::Long;
use Config::OptionsSet::OptionsSet qw(readOpts writeOpts prettyPrintOpts overwriteOpts 
					setDiff hashDiff hashSymDiff);
use AI::GMTK::Util  qw(	prepareGmtkCmdArgs checkOpts
					isRestartable notifyOfJobCompletion);
use AI::GMTK::EmConvergeParallel qw(emConvergeParallel);

my $usage = <<EOS;
Usage: emConvergeParallel [-h] [-d <distributedConfigFile>] <gmtkArgFile> <workingDir>

Splits/Merges gaussians and emtrains to convergence according to 
some split/merge and convergence schedules.

The spliting/merging is governed by the split/merge schedule.  
Each split/merge step is followed by emtraining to convergence
, which has its own convergence schedule which can be 
different from step to step.

Required:
<gmtkArgFile>  
	a file with arguments to gmtkEmtrainNew.  It must be in the 
	same format as a file used by gmtkEmtrainNew -gmtkArgFile flag,
	with additional requirement of one flag per line.
	
<workingDir> a dir for all temprary files.
Options:
-d --distributedConfig 
	<distributedConfigFile>  a file containing distributed 
	settings in format of one VAR VALUE pair per line.
	To see the possible options just run emConvergeParallel without the -d 
	option and all legal settings will be printed.
	Note that any setting starting with _ cannot be set in 
	the distributedConfigFile.
	
	This file also contains the split/vanish schedule and convergence schedules 
	as in the following example. This will split 7 times, doubling number of gaussians each time
	#and afterwards converge to 1-llratio .002 (So if you started with one component, you will end with 128 components)
	botForceVanishArray 0 0 0 0 0 0 0
	topForceSplitArray 0 0 0 0 0 0 0
	mcsrArray 1e-200 1e-200 1e-200 1e-200 1e-200 1e-200 1e-200
	mcvrArray 2e200 2e200 2e200 2e200 2e200 2e200 2e200
	llRatioThresholdArray .02 .02 .02 .02 .02 .02 .02 .002
	maxEmIterationsArray  20  20  20  20  20  20  20  20

	llRatio is calculated as 1-currentLogLikelihood/previousLogLikelihood

-h --help 
	print this help message

Return values:
	0  success
	-1 cofiguration error. no files written or modified.
 	-2 sanity check failed.  No parallel jobs launched.
	-3 some parallel tasks failed.

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
	my $help;
	
	GetOptions (
				"distributedConfig|d=s" => \$distrConfigFile,
				"help|h" => \$help);
	
	
	my $commands_file;
	if (@ARGV < 2 || $help) {
		die("$usage");
	} else {
		$gmtkArgFile = $ARGV[0];
		$workingDir = $ARGV[1];
	}
	
	
	#read gmtk and distributed options
	my %gmtkOpt;
	readOpts($gmtkArgFile, \%gmtkOpt);

	my %userDistrOpt=();
	$userDistrOpt{'wd'} = $workingDir;
	if($distrConfigFile){
		readOpts($distrConfigFile, \%userDistrOpt);
	}

	my $ret = emConvergeParallel(\%gmtkOpt,\%userDistrOpt);
	notifyOfJobCompletion(\%userDistrOpt, $ret);

	print "Exiting with code $ret.\n";
	exit ($ret);
}


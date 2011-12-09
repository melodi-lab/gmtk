#!/usr/bin/env perl
## @file 
# 
#  @author Arthur Kantor 
use warnings;
use strict;

use Getopt::Long;
use AI::GMTK::EmTrainParallel qw(emtrainParallel);
use Config::OptionsSet::OptionsSet qw(readOpts writeOpts prettyPrintOpts overwriteOpts 
					setDiff hashDiff hashSymDiff);
use AI::GMTK::Util  qw(	prepareGmtkCmdArgs checkOpts
					isRestartable notifyOfJobCompletion);


my $usage = <<EOS;
Usage: emtrainParallel [-h] [-d <distributedConfigFile>] <gmtkArgFile> <workingDir>

Runs gmtkEmtrainNew in parallel via SGE. A single EM iteration is performed.

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
	To see the possible options just run emtrainParallel without the -d 
	option and all legal settings will be printed.
	Note that any setting starting with _ cannot be set in 
	the distributedConfigFile on the command line.

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
	$userDistrOpt{'gmtkArgFile'} = $gmtkArgFile; 
	if($distrConfigFile){
		readOpts($distrConfigFile, \%userDistrOpt);
	}

	my $ret = emtrainParallel(\%gmtkOpt,\%userDistrOpt);
	notifyOfJobCompletion(\%userDistrOpt, $ret);
	exit ($ret);
}


#!/usr/bin/env perl
## @file 
# 
#  @author Arthur Kantor 

use warnings;
use strict;
#use File::Spec;
#use File::Basename;
#use File::Temp qw/ tempfile tempdir /;
#use Cwd;
use Parallel::Distribute qw(distribute);
use Getopt::Lazier;

#parses options and checks if the files/dir specified as parameters exist
sub parseOpts{
    my $usage = <<'EOS';

Usage: $prog [--help] [$prog options] [--commandFile file] [qsub/sbatch specific opts]

    Launches jobs from commandFile in parallel via SGE or SLURM.
    
Options:

$opts
    
Unrecognized options are passed on qsub or sbatch or whatever command is actually used to 
submit the jobs.  The  qsub or sbatch options must follow the options to $prog. 


Based on distribute.pl from JHU CLSP workshop 2006.

EOS
    my $optSpec = [
        ['help|h',
        undef,
        'print this help message and exit'],
        ['commandsFile|c=s',
        undef,
        'A file with one command to be executed per per line.',
        'required fileExists'],
        ['jobName|N=s',
        undef,
        'job name (see -N in qsub). If omited defaults to name of executable in the commandFile'],
        ['sgeWorkDir|s=s',
        'SGE',
        'path for the working directory of this script. The dir is created if it does not exist. Task stdout streams are put here, and temporary script files as well.'],
        ['runLocally|l',
        undef,
        'Do not submit to cluster and just run each command in a bash shell instead.'],
    ];
    
    
    my %opts;
    Getopt::Lazier::Configure('pass_through');
	exit(-1) if (!parseOptions(\%opts, $optSpec));
	if($opts{'help'}){
		standardUsage($optSpec, $usage);
		exit(0);
	}	
	
	my @err=validateOptions(\%opts, $optSpec);
	if (@err){
		print "@err";
		standardUsage($optSpec, $usage);
		exit(-1);
	}

	#standardUsage($optSpec, $usage);
    return %opts;
}


my %opts = parseOpts();
my $extraOptions = join (' ', @ARGV);
#print "$extraOptions \n";

if ($opts{'sgeWorkDir'} !~ /^\//){ #relative path
    my $cwd = `echo pwd | bash`;
    chomp $cwd;
    $opts{'sgeWorkDir'} = "$cwd/$opts{'sgeWorkDir'}";
}


my $ret = distribute($opts{'commandsFile'}, $opts{'jobName'}, $opts{'sgeWorkDir'}, $opts{'runLocally'}, $extraOptions);

exit $ret;


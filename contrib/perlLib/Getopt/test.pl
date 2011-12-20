#!/usr/bin/env perl
use warnings;
use strict;

use Carp;
use Getopt::Lazier;
## @file 
#test and example of Lazier.pm
#  @author Arthur Kantor 
#  @date 10/10/08



#parses options and checks if the files/dir specified as parameters exist
sub parseOpts{

	my $usage = <<'EOS';

Usage: 
$prog [--help] [options]  

An example documentation for a program.  Write and format things here as you like, in any language.
Two special variables are recognized here: 
\$prog is the base name of the script
\$opts is the return value of docOptions(...)

Options:

$opts

Return values:
	0  success
	-1 error

EOS

    my $optSpec = [
		['help|h', 
		 undef,
		 'print this help message and exit'],
		['minimalisticOption|m=s'],
		['optionalStringArg|o=s',
		 undef, 
		 'This string option is documented here.'],
		['optionalStringArg2|aliasForOptionalStringArg2|o=s',
		 'Hello', 
		 'This string pption is has a default.'],
		['requiredStringArg|r=s', 
		 'hideho', 
		 'This is a required string argument as indicated by the "r" validator. If the parameter is not specified, validateOptions() will return an error.',
		 'required'],
		['someRequiredFile=s', 
		 undef,
		 'This argument is required to be present, and must be an existing filename as indeicated by the "r" and "f" validators.  If either of these validators fail, validateOptions() will return an error.',
		 'required fileExists'],
		['someOptionalDir|d=s', 
		 undef,
		 'This argument is not required to be present, but if it is present, it must be the name of an existing dir.',
		 'dirExists'],
	];


    my %opts;
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


#parse and check the arguments
my %opts=parseOpts();	









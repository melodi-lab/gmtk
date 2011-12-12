#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Lazier;
use OS::Util  qw(tsSystem );
use Carp;

#  Arthur Kantor 7/2/09



#parses options and checks if the files/dir specified as parameters exist
sub parseOpts{
    my $usage = <<'EOS';
$prog [--help] [options]  
 
Converts some float columns from an in-pfile into discrete features in the out-pfile.
This is to get around bugs in obs-print, pfile_create, labcat and feacat, which seem to 
stem from inability to handle more than one label column (labcat, feacat) and inability 
of obs-print to read the ascii <uttId> <frameId> <ftr1> <ftr2> ... <label1> <label2> ... format.

Options:

$opts
EOS
    
    my $optSpec = [
		["help|h", 
		 undef,
		 "print this help message and exit"],
		["inPfile|i=s", 
		 undef, 
		 "File containing the discrete features",
		 'required fileExists'],
		["outPfile|o=s", 
		 undef, 
		 "Where the discrete features will be stored",
		 'required'],
		["startFeatureRange|s=i",
		 0,
		 "The start of the discrete feature (label) range.",
		 'required'],
		["endFeatureRange|e=i",
		 undef,
 		 "The end of the discrete feature (label) range. Default: same as startFeatureRange.",
		],
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

	$opts{'endFeatureRange'} = $opts{'startFeatureRange'} if (!defined($opts{'endFeatureRange'}));
	#print "Script started on " . localtime(time)." with args: @ARGV\n";

    return %opts;
}

my $cmdLine=join(' ',@ARGV);
my %opts = parseOpts();
print "Script started on " . localtime(time)." with args: $cmdLine\n";

my $cnt=1;
my $mergeCmd = "obs-print -q -ofmt pfile -o $opts{outPfile} ";
for (my $i=$opts{'startFeatureRange'}; $i <= $opts{'endFeatureRange'}; $i++){
	print "generating $opts{outPfile}.$i\n";
	tsSystem("feacat $opts{inPfile} -op ascii -fr $i | labcat -ip ascii -op pfile -o $opts{outPfile}.$i");
	$mergeCmd .= " -i$cnt $opts{outPfile}.$i";
	$cnt++;
}
print "merging...\n";
tsSystem($mergeCmd);
for (my $i=$opts{'startFeatureRange'}; $i <= $opts{'endFeatureRange'}; $i++){
	print "deleting $opts{outPfile}.$i\n";
	tsSystem("rm $opts{outPfile}.$i &");
}

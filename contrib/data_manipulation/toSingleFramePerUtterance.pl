#!/usr/bin/env perl
use warnings;
use strict;

use Fatal qw(:void open close); #die on failed open and close if the return values are not checked

use Getopt::Lazier;
use OS::Util  qw(readLines writeLines tsSystem );
use List::Util qw[min max];
use File::Temp qw/ :mktemp  /;
use Carp;

#  Arthur Kantor 11/03/09

#parses options and checks if the files/dir specified as parameters exist
sub parseOpts{
    my $usage = <<'EOS';
$prog [--help] [options]  --inPfile <FILE> --outPfile <FILE>
 
 Turns every frame in inPfile into its own utterance and stores it in outPfile.

Options:

$opts

EOS

    my $optSpec = [
		['help|h', 
		 undef,
		 'print this help message and exit'],
		['inPfile|i=s', 
		 undef, 
		 'input Pfile. - means STDIN.',
		 'required'],
		['outPfile|o=s', 
		 undef, 
		 'output Pfile.',
		 'required'],
		['sourceUttFrameIds|f=s', 
		 undef, 
		 'output text files where to store the uttId frameId pair for each frame in the source pfile.',
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

    return %opts;
}

my %opts = parseOpts();

open(my $IN, "obs-print -i1 $opts{'inPfile'} -q |");  
open(my $IDS_FH, "> $opts{'sourceUttFrameIds'}") if ($opts{'sourceUttFrameIds'});
my ($out, $outfn) = mkstemps( "SingleFrameAsciiXXXX", '.txt');

my ($uttNo,$frameNo,$c);
while(<$IN>){
	($uttNo,$frameNo,$c)=split(" ", $_, 3); 
	print $IDS_FH "$uttNo $frameNo\n" if ($opts{'sourceUttFrameIds'});
	print $out ( $.-1, ' 0 ', $c);
}
close($out);
close($IDS_FH) if ($opts{'sourceUttFrameIds'});
close($IN);

my @n =split(' ', $c);
tsSystem("feacat -o $opts{'outPfile'} -i $outfn -ip ascii -width ".scalar(@n));
unlink $outfn;
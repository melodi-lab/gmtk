#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;
use File::Spec;

## @file

my $exe = basename $0;
my $usage = "$exe *.cm *.dg1 *.dg2 *.glo *.nas *.phn *.pl1 *.pl2 *.rd *.vow *.wd\n";
die $usage if $#ARGV < 0;

foreach my $file (@ARGV) {
    my $string = "";
    open FILE, "<$file";
    while (<FILE>) {
        chomp;
        if (/^\s*\d+\s+\d+\s+\S+\s*$/) {
            my ($startSamp, $stopSamp, $label) = split;
            $string .= sprintf("%.07f %.07f %s\n",$startSamp/8000,$stopSamp/8000,$label);
        }
        else {
            $string .= "$_\n";
        }
    }
    close FILE;
    
    open FILE, ">$file";
    print FILE $string;
    close FILE;
}

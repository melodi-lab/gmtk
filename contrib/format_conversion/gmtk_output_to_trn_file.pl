#!/usr/bin/env perl

## @file

# parse arguments
if (@ARGV != 2) {
    die("Usage: gmtk_output_to_trn_file.pl gmtk-output basenames-file\n");
}

$gmtk_output = $ARGV[0];
$basenames = $ARGV[1];

open(BNAMES, $basenames) || die "Couldn't open $basenames for reading.  Exiting.\n";
@basenames_lines = <BNAMES>;
close(BNAMES);

$syscall = "cat $gmtk_output | perl -pe 's/" . '\n/ /g' . "' | perl -pe 's/ (Segment|Total data)/" . '\n \1/g' . "' | grep 'Segment' | perl -pe 's/ " . '\(\d+\-\d+\)//g' . "' | perl -pe 's/" . '\s*Segment .+ per numUFrams = \S+\s+//' . "' | perl -pe 's/(" . '\<SILENCE\>)|(\<\\\?S\>)//g' . "' | perl -pe 's/ +/ /g' | perl -pe 's/^ //' | perl -pe 'tr/[a-z]/[A-Z]/' | perl -pe 's/SKIPPING SEGMENT SINCE PROBABILITY IS ESSENTIALLY ZERO//' > $gmtk_output.tmp";

system($syscall);

open(GMTK, "$gmtk_output.tmp") || die "Couldn't open $gmtk_output.tmp for reading.  Exiting.\n";
@gmtk_lines = <GMTK>;
close(GMTK);

foreach $i (0..$#gmtk_lines) {
    $bline = $basenames_lines[$i];
    chomp($bline);

    $gline = $gmtk_lines[$i];
    chomp($gline);

    print STDOUT "$gline ($bline)\n";
}

system("rm $gmtk_output.tmp");

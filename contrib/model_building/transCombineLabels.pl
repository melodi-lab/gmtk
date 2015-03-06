#!/usr/bin/perl

## @file

$0 =~ s%^.*/%%;

if (@ARGV != 1) {
    print STDERR "Usage:  $0 trans\n";
    exit;
}

$last_lab = "";
$last_start = 0;
open FILE, "<$ARGV[0]";
@lines = <FILE>;
close FILE;
open FILE, ">$ARGV[0]";
foreach $line (@lines) {
    chomp($line);
    ($start, $end, $lab) = split(/\s+/, $line);
#    if (($last_lab ne "") && ($lab ne $last_lab)) {
#    if ($last_start lt $start) {
    if ($lab ne $last_lab) {
	if ($last_lab ne "") {
	    print FILE "$last_start $start $last_lab\n";
	    $printed_something = 1;
	}
	$last_start = $start;
	$last_lab = $lab;
    }
}

if ($printed_something) {
    print FILE "$last_start $end $lab\n";
} else {
    print FILE "0 $end $lab\n";
}
close FILE;

#!/usr/bin/env perl

## @file
#################################################################
##
## preprocess_file.pl
##
## Make a set of variable replacements in the standard input,
## and output result to standard output.
##
## -KL 5/26/02
##
#################################################################

while (@ARGV) {
    $curr_arg = pop(@ARGV);
    if ($curr_arg =~ /^(\S+)=(\"?.*\"?)$/) {
	$var = $1;
	$value = $2;
	$replace{$var} = $value;
    } else {
	die 'Arguments must take the form VAR=["]value["].  Exiting.\n';
    }
}

while ($line = <STDIN>) {
    foreach $var (keys %replace) {
	$line =~ s/\*$var\*/$replace{$var}/;
    }
    print $line;
}


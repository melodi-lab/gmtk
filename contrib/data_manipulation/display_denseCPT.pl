#!/usr/bin/perl

## @file
# Prettier view of a denseCPT, showing parent values (and labels, if requested)
# The varmap-file argument is a file consisting of one line per parent variable,
# consisting of a filename that maps that parent's values to labels; use "ID"
# instead of the filename to use an identity mapping (i.e. just display the values
# themselves).
#
# See example.denseCPT and example.parmap for example input files.

$0 =~ s%^.*/%%;

if (@ARGV != 2) {
    die "Usage:  $0 denseCPT-file varmap-list-file\n (where 1st line of denseCPT-file is CPT name)\n";
}

open(CPT, $ARGV[0]) || die "Couldn't open denseCPT file $ARGV[0].\n";
@cpt_lines = <CPT>;
close(CPT);

# remove comments from CPT
foreach $i (0..$#cpt_lines) {
    $cpt_lines[$i] =~ s/\s+\%.*\n$//;
    chomp($cpt_lines[$i]);
}
@cpt = split(' ', join(' ', @cpt_lines));  # this join might be very large

open(MAPLIST, $ARGV[1]) || die "Couldn't open var map list file $ARGV[1].\n";
@map_lines = <MAPLIST>;
close(MAPLIST);

foreach $i (0..$#map_lines) {
    $file = $map_lines[$i];
    chomp($file);
    open(MAP, $file) || die "Couldn't open var map file $file.\n";
    @map = <MAP>;
    foreach $line (@map) {
	chomp($line);
	if ($line =~ /^(\d+) (\S+)$/) {
	    $varmap{$i . "_" . $1} = $2;
	} else {
	    die "ERROR:  Incorrectly formatted line in map file $file.  Offending line: $line\n";
	}
    }
}

# header of denseCPT must consist of denseCPT name, num parents, cardinalities
$name = shift(@cpt);
$num_parents = shift(@cpt);
if ($num_parents !~ /^[1-9][0-9]*$/) { # check that it is an int
    die "ERROR: Found $num_parents while looking for an int.\n";
}

@cards = ();
foreach $i (1..$num_parents) {
    push(@cards, shift(@cpt));
}

$child_card = shift(@cpt);

print STDOUT "CPT name = $name\n";
print STDOUT "num parents = $num_parents\n";
print STDOUT "parent cardinalities = ", join(' ' , @cards), "\n";
print STDOUT "child cardinality = $child_card\n";
print STDOUT "num probs in CPT = ", $#cpt+1, "\n";

foreach $i (0..@cpt/$child_card-1) {
    @pmf = splice(@cpt, 0, $child_card);
    $parent_vals[$num_parents-1] = $i - $cards[$num_parents-1] * int($i/$cards[$num_parents-1]);  # $i mod $cards[$num_parents-1]
    $divider = 1;
    foreach $p (reverse 1..$num_parents-1) {
	$divider = $divider*$cards[$p];
#	print "p = $p, divider = $divider, card = $cards[$p]\n";
	$parent_vals[$p-1] = int($i/$divider);
    }
    foreach $p (0..$num_parents-1) {
	print STDOUT "p$p=", $varmap{$p . "_" . $parent_vals[$p]}, " ";
    }
    print STDOUT join(' ', @pmf), "\n";
}

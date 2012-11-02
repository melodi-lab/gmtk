#!/usr/bin/env perl
#
# $Header: /export/ws06afsr/src/parallel/breakRange.pl,v 1.1.1.1 2006/06/30 21:46:08 plal Exp $
# Jeff Bilmes <bilmes@ee.washington.edu>

# breakRange: breaks a range in to n equal sub-ranges, the union
# of the n ranges is the original range.


sub usage {
 printf STDERR (
  "Description:
   breakRange: breaks a range in to n equal sub-ranges, the union
   of the n ranges is the original range.\n" .
  "Usage: $0 [options]\n" .
  "   -r rng     Range\n" .
  "   -n #       Number of subranges\n" .
  "");
}

sub checkfloat {
    ($opt,$arg) = @_;
    if (( $arg !~ /^[0-9.]+$/ )) {
	print STDERR ("Error: Value \"", $arg, "\" invalid for option ",
		      $opt, " (floating point value expected)\n");
	&usage; exit(-1);
    }
}

sub checknnint {
    ($opt,$arg) = @_;
    if (( $arg !~ /^[0-9]+$/ ) || ($arg < 0)) {
	print STDERR ("Error: Value \"", $arg, "\" invalid for option ",
		      $opt, " (positive integer expected)\n");
	&usage; exit(-1);
    }
}

# parse arguments
if ( $#ARGV < 0 ) {
    &usage;
    exit(-1);
}
while ( $#ARGV >= 0 ) {
    $opt = shift(@ARGV);
    $opt =~ tr/A-Z/a-z/; # ignore case
    if ( $opt =~ /^-r/ ) {
	if ( $#ARGV < 0 ) {
	    print STDERR ("Error: Option ", $opt, " requires an argument\n");
	    &usage; exit(-1);
	}
	$arg = shift(@ARGV);
	$rng = $arg;
    } elsif ( $opt =~ /^-n/ ) {
	if ( $#ARGV < 0 ) {
	    print STDERR ("Error: Option ", $opt, " requires an argument\n");
	    &usage; exit(-1);
	}
	$arg = shift(@ARGV);
	&checknnint($opt,$arg);
	$n = $arg;
    } elsif ( $opt =~ /^-help/ ) {
	&usage;
	exit(0);
    } else {
	printf STDERR ("Error: Unknown option (%s)\n",$opt);
	&usage;
	exit(-1);
    }
}

if (!(defined($rng) && defined($n))) {
    printf STDERR ("Argument error\n");
    &usage;
    exit(-1);
}

($rng =~ /^-?\d+:-?\d+$/) || die("malformed range $rng, need form n1:n2");

($start,$end) = split(/:/,$rng);

$items = $end-$start+1;
($items >= $n) || die("Number of subranges ($n) must be <= $items");

$num_per_sub = int($items/$n);
$rem = $items % $n;

$srend = $start-1;
for ($i=0;$i<$n;$i++) {
    $srstart = $srend + 1;
    $srend = $srstart + $num_per_sub - 1;
    if ($rem > 0) {
	$srend ++;
	$rem --;
    }
    if ($i == ($n-1)) {
	$srend = $end;
    }
    print "$srstart:$srend "
}
print "\n";

#!/usr/bin/perl

## @file
#################################################################
##
## generate_phoneState2featureDTs.pl
##
## This script takes in a phone-state-to-feature mapping file 
## and generates a DT for a given feature
##
#################################################################


sub usage {
 printf STDERR (
  "Description: Generate a DT mapping from phone state to a given feature.\n" .
  "Usage: progname [options]\n" .
   "   -p file        phone-state-to-feature mapping filename\n" .
   "   -f string      feature-name\n"
);
}

sub checkfile {
    ($opt,$arg) = @_;
    if (( $arg !~ /^\S+$/ )) {
        print STDERR ("Error: Value \"", $arg, "\" invalid for option ",
                      $opt, " (filename expected)\n");
        &usage; exit(-1);
    }
}

sub checkstring {
    ($opt,$arg) = @_;
    if (( $arg !~ /^\S+$/ )) {
        print STDERR ("Error: Value \"", $arg, "\" invalid for option ",
                      $opt, " (string with no whitespace expected)\n");
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
    if ( $opt =~ /^-p/ ) {
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        &checkfile($opt,$arg);
        $phone2feature = $arg;
    } elsif ( $opt =~ /^-f/ ) {
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        &checkstring($opt,$arg);
        $feature = $arg;
    } elsif ( $opt =~ /^-help/ ) {
        &usage;
        exit(0);
    } else {
        printf STDERR ("Error: Unknown option (%s)\n",$opt);
        &usage;
        exit(-1);
    }
}


# read in phone2feature mapping
open(PH2FEAT, $phone2feature) || die "Couln't open phone-to-feature file $phone2feature for reading.  Exiting.\n";

while ($line = <PH2FEAT>) {
    if ($line =~ /^\#\s+\S+(\s+)Index.+$/) {  # this is the comment line defining the order of the features
	# figure out the index of the feature of interest
	($hash, $phone, $index, @feature_names) = split(/\s+/, $line);

	foreach $i (0..$#feature_names) {
	    if ($feature_names[$i] eq $feature) {
		$feature_index = $i;
	    }
	}
    } elsif (($line =~ /^\s*$/) || ($line =~ /^\#.*$/)) {  # another comment or empty line--skip it
        next;
    } elsif ($line =~ /^([a-z\d_]+)\s+(\d+)[\.\,\(\)\s\d\*]+$/) {
        $num_phonemes++;
        ($ph_name, $ph_index, @feat_values) = split(/\s+/, $line);
	$phoneTable{$ph_name} = $ph_index;
        foreach $i (0..$#feat_values) {
	    if ($feat_values[$i] =~ /^\d+$/) {
		$phone2feat{$ph_index."_".$i} = $feat_values[$i];
	    } else {
		$phone2feat{$ph_index."_".$i} = 0;  # assuming that 0 is always a valid feature value
	    }
        }
    } else {
        die "ERROR: Incorrectly formatted line in phone2artic file: $line\n";
    }
}

close(PH2FEAT);

print "1 % number of DTs in this file\n";
print "0 % DT number\n";

print "phoneState2", $feature, "DT % DT name\n";
print "1 % num parents\n";
print "0 $num_phonemes ";

for ($ph=0;$ph<$num_phonemes-1;$ph++) {
    print "$ph ";
}
print "default\n";

for ($ph=0;$ph<$num_phonemes;$ph++) {
#    print STDERR "phone = $ph, feat index = $feature_index, val = ", $phone2feat{$ph."_".$feature_index}, "\n";
    print "\t -1 ", $phone2feat{$ph."_".$feature_index}, "\n";
}

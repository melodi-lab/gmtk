#!/usr/bin/perl -w

## @file 
# Based on genPhonePhonePos2WholePhoneStateDTs.pl from the gmtk Aurora tutorial

sub usage {
 printf STDERR (
  "Description: writes out a DT which maps from two RVs, cardinalities A and B, to one with cardinality AxB.  Handy for going from stateCounter and phone to state, say.  \n".
  "Usage: progname [options]\n" .
   "   -n        name of the DT\n" .
   "   -a        cardinality of first parent\n" .
   "   -b        cardinality of second parent\n" 
);
}

sub checknnint {
    ($opt,$arg) = @_;
    if (( $arg !~ /^[0-9]+$/ ) || ($arg < 0)) {
        print STDERR ("Error: Value \"", $arg, "\" invalid for option ",
                      $opt, " (positive integer expected)\n");
        &usage; exit(-1);
    }
}

#default values:
$numA = 3;
$numB = 61;
$dtName = "stateCounterAndPhoneToState"; 

# parse arguments
while ( $#ARGV >= 0 ) {
    $opt = shift(@ARGV);
    $opt =~ tr/A-Z/a-z/; # ignore case
    if ( $opt =~ /^-b/ ) {
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        &checknnint($opt,$arg);
        $numB = $arg;
    }elsif ( $opt =~ /^-n/ ) {
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        $dtName = $arg;
    }elsif ( $opt =~ /^-a/ ) {
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        &checknnint($opt,$arg);
        $numA = $arg;
    } elsif ( $opt =~ /^-help/ ) {
        &usage;
        exit(0);
    } else {
        printf STDERR ("Error: Unknown option (%s)\n",$opt);
        &usage;
        exit(-1);
    }
}

#header information for this DT file
print "1 % number of DTs in this file\n";

#obviously this is the first and only DT
print "0 % DT number\n";
print "$dtName % DT name\n";
print "2 % num parents\n";

#the root node queries the second parent
print "1 $numB ";
if($numB > 3){
    printf "0 ... %d ", $numB-2;
}else{
    for ($b=0 ; $b<$numB-1 ; $b++) {
        print "$b ";
    }
}
print "default\n";

#now make children for each of the options described in the root
$c = 0;
for ($b=0 ; $b<$numB ; $b++) {
    #next, the tree should query the first parent - there'll be $numA options
    print "\t 0 $numA ";
    if($numA > 3){
        printf "0 ... %d ", $numA-2;
    }else{
        for ($a=0 ; $a<$numA-1 ; $a++) {
            print "$a ";
        }
    }
    print "default\n";

    #these are the leaves, one per first parent value
    for ($a=0;$a<$numA;$a++) {
        print "\t\t-1 $c\n";
        $c++;
    }
}

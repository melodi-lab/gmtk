#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Config::OptionsSet::OptionsSet2 qw(readOpts writeOpts prettyPrintOpts expandOpts cppArgsToCppBlock cppBlockToCppArgs );

my $opts={};
my %con=(MODEL_DIR => 'modelDir/dir', DATA_DIR=>'dataDir');
readOpts($ARGV[0],$opts,\%con);
#print Dumper($opts);
#prettyPrintOpts($opts);
my $waOpts=expandOpts($opts,'wordAligned.trainTriUnit.splitConverge');
print("**********\n");
prettyPrintOpts($waOpts);
print("**********\n");
my $waOpts2 = cppBlockToCppArgs($waOpts);
prettyPrintOpts($waOpts);
prettyPrintOpts($waOpts2);
prettyPrintOpts(cppArgsToCppBlock($waOpts2));
writeOpts($opts,$ARGV[1]);
 
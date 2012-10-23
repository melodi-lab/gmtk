#!/usr/bin/perl -w

## @file
#Takes an htk scp format file (i.e. one or two filenames per line) and changes the path to make it work on the local filesystem (i.e. prepends $ROOTDIR)
#Assumes absolute filenames
#Takes input on STDIN and writes to STDOUT

$rootdir = $ENV{'ROOTDIR'};
print STDERR "Using ROOTDIR \"$rootdir\".  \n";

while(<STDIN>){
   @files = split;
   for($i=0 ; $i-1 < $#files ; $i++) {
      print $rootdir.$files[$i];
      if($i < $#files){print " "};
   }
   print "\n";
}

#!/usr/bin/perl -w
## @file 
# Modified version of genCounterToWordMapDTs.pl included with the gmtk Aurora tutorial
# Used for getting the AVICAR DTs

sub usage {
 printf STDERR (
  "Description: Generate a file containing a list of DTs, mapping from counter to word\n" . 
  "Usage: genCounterToWordMapDTs [options]\n" .
  "   -v file     vocabulary file, one word per line, remainder of line after the first word is ignored.\n" .
  "   -m N        maximum number of utterances to use\n" .
  "   -e ext      file extension of items in the MLF file\n" .
  "   -f file     HTK file (MLF) containing the training data\n");
}

#checks that the argument $arg for option $opt is a positive integer
sub checknnint {
    ($opt,$arg) = @_;
    if (( $arg !~ /^[0-9]+$/ ) || ($arg < 0)) { #the second half of the || is redundant, I think.  Swapping them around might help with any lazy evaluation but this is a benign mistake as it is.  - Partha
        print STDERR ("Error: Value \"", $arg, "\" invalid for option ",
                      $opt, " (positive integer expected)\n");
        &usage; exit(-1);
    }
}


#bail out if there are no arguments
if ( $#ARGV < 0 ) {
    &usage;
    exit(-1);
}

# parse arguments
$file="/dev/stdin/";
$ext="lab";
$maxuts = 1e8;
$eou = ".";#end of utterance, could be . sometimes, make this an arg

while ( $#ARGV >= 0 ) {
    $opt = shift(@ARGV);
    $opt =~ tr/A-Z/a-z/; # ignore case
    if ( $opt =~ /^-v/ ) {#the vocab file
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        $vocab = $arg;
    } elsif ( $opt =~ /^-e/ ) {#file extension
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        $ext = $arg;
    } elsif ( $opt =~ /^-f/ ) {#the MLF file
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        $file = $arg;
    } elsif ( $opt =~ /^-m/ ) {#max num utts
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        &checknnint($opt,$arg);
        $maxuts = $arg;
    } elsif ( $opt =~ /^-help/ ) {
        &usage;
        exit(0);
    } else {
        printf STDERR ("Error: Unknown option (%s)\n",$opt);
        &usage;
        exit(-1);
    }
}

#open the MLF file to filehandle F
open(F,"$file") || die("cannot open file $file");

#move the filepointer past the header
while(<F>) {
    chomp;
    # print "read ",$_,"\n";
    if ($_ eq "#!MLF!#") {
	last;
    }
}

# The result of this while loop is that, upon termination, $numUts holds the lesser value out of $maxuts and the number of utterances in the MLF file (as determined by the number of file extensions)
$numUts = 0;
while (($record = <F>)) {
   if($record =~ /$ext/){
       $numUts++;
   }
   if ($numUts > $maxuts) {
       $numUts = $maxuts;#this is probably unnecessary
       last;
   }
}


# go to beginning
seek(F,0,0);

print "% generated on ".localtime()." from: \n%mlf file $file, \n%vocab file $vocab \n\n";
print "\n\n", $numUts," % Number of DTs\n\n";

#skip the header again
while(<F>) {
    chomp;
    if ($_ eq "#!MLF!#") {
	last;
    }
}


%wordTable = load_vocab($vocab);
my %inversewordTable=reverse(%wordTable);

#iterate through the MLF file 
for ( $utt=0; $utt<$numUts; $utt++) {
    $_ = <F>;

    #grab the speakerid and sentenceid and create the decision tree header
    ($sentenceid) = /\/([^\/]*)\.$ext/;
    if ($utt==0){
        print $utt," % DT number\n";
        print "$sentenceid % DT name\n";
        print "1 % number of parents\n";
    }
    else{
        #don't print comments to keep the DT files small (with a better chance of being under the 2gb limit)
        print "$utt\n$sentenceid\n1\n";
    }    
    
    #this while loop will populate a map from index to word id
    $wordPos = 0;
    @wordMap = ();
    while (1) {
		$field = <F>;
		chomp($field),
		$field =~ s/\s+//g;
		if ($field  eq $eou) {#i.e. end of utterance
			last;
		}
		if (!defined($wordTable{$field})) {
			printf STDERR "ERROR, undefined word $field in file $sentenceid\n";
		}
		$wordMap[$wordPos] = $wordTable{$field};
		$wordPos++;
    }

    #create the information for this node
    $out_of_bounds = $wordPos;
    if($out_of_bounds>1){
        printf "0 $out_of_bounds ";
        print " 0 ... ".($out_of_bounds-2) if ($out_of_bounds-2 >=0);
        printf " default\n";
    }
    #else This is a degenerate decision tree, where the parent can take only a single value.
    #     So don't depend on the parent and just return the value always

    
    #now print the leaves
    for ($i=0;$i<$wordPos;$i++) {
        print "\t-1 $wordMap[$i] \% $inversewordTable{$wordMap[$i]}\n";
    }

    print "\n\n";
}


sub load_vocab{
    my $vocab_file = shift;
    my %table = ();
    my $index = 0;

    open(VOCAB,$vocab_file) or die "Couldn't open vocabulary file $vocab_file\n";

    while(<VOCAB>){
        chomp;
        if(/^#/){
           next;
        }
		s/[\t ].*//;
        $word = $_;
        $table{$word} = $index++;
    }
    close(VOCAB);

    print STDERR "$index items in vocab\n";

    return %table;
}

#!/usr/bin/perl

## @file 

sub usage {
 printf STDERR (
  "Description: Generate a file containing a decision tree that determines what phone is present, given the word, pronunciation variant and phoneCounter.\n " .
  "Usage: genDictionary [options]\n" .
  "   -v N	  vocabulary size\n" .
  "   -n str 	  DT name\n" .
  "   -p file     phoneme map file\n" .
  "   -d file     HTK file (DCT) containing pronunciation dictionary.  Format as described in htkbook, also assumes the file is sorted by word.  \n");
}

$dict = "/dev/stdin";
$dtName = "dictionary";

#bail out if there are no arguments
if ( $#ARGV < 0 ) {
    &usage;
    exit(-1);
}

# parse arguments
while ( $#ARGV >= 0 ) {
    $opt = shift(@ARGV);
    $opt =~ tr/A-Z/a-z/; # ignore case
    if ( $opt =~ /^-v/ ) {#the vocab file
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        $vocab_size = $arg;
    } elsif ( $opt =~ /^-p/ ) {#the name of the DT
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        $phoneme_map = $arg;
    } elsif ( $opt =~ /^-n/ ) {#the name of the DT
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        $dtName = $arg;
    } elsif ( $opt =~ /^-d/ ) {#the dictionary file
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        $dict = $arg;
    } elsif ( $opt =~ /^-help/ ) {
        &usage;
        exit(0);
    } else {
        printf STDERR ("Error: Unknown option (%s)\n",$opt);
        &usage;
        exit(-1);
    }
}

if (! defined $dict){
   print STDERR "Error: you need to specify a dictionary with -d\n";
   &usage;
   exit(-1);
}

#loading the phoneset into a table from phone to index (starts at 0)
# read in phoneme map table (swiped verbatim from generate_phoneStateDT.pl)

open(PH, $phoneme_map) || die "Couldn't open phoneme map file $phoneme_map for reading.  Exiting.\n";

while ($line = <PH>) {
    if ($line =~ /^(\d+)\s+(\S+)\s*$/) {
        $phoneTable{$2} = $1;
    }
}

close(PH);

#write the top of the DT file
print "1 % number of DTs in this file\n";
print "0 % DT number\n";
print "$dtName % DT name\n";
print "3 % three parents - word(0), pronunciation(0) and phoneCounter(0)\n";

#the decision tree first queries on word
printf("%s\n",branching_node(0,$vocab_size));

#now is a good time to open the dictionary
open(DICT,$dict) or die "Couldn't open dictionary $dict\n";

#loop through the dictionary
$previous_word = "";
@pronunciations = ("sil");
while(<DICT>){
    if(/^(\S+)[\ \[\]0-9\.]*(.+)/){
        $word = $1;
        $pronunciation = $2;

        if($word ne $previous_word && $previous_word ne ""){
             printf("\t%s %% %s\n",branching_node(1,scalar @pronunciations),$previous_word);

             foreach $p (@pronunciations){
                 @current_pronunciation =  split(/\ +/,$p);
                 printf("\t\t%s\n",branching_node(2,scalar @current_pronunciation));
                 
                 foreach $phone (@current_pronunciation){
                     $index = $phoneTable{$phone};
                     if(! defined $index){ print STDERR "$phone appears in $dict but not in $phoneme_map\n"; exit(-1);}
                     printf("\t\t\t-1 %d %% %s\n",$index,$phone);
                 }
             }

             @pronunciations = ($pronunciation);
        }else{
             push @pronunciations,$pronunciation;
        }
        $previous_word = $word;
    }else{
        #will end up here if the dictionary line contained an OUTSYM or a pronunciation probability (valid according to the htk spec of dct files but unsupported here)
        print STDERR "Didn't fit:$_";
    }
}

printf("\t%s %% %s\n",branching_node(1,scalar @pronunciations),$previous_word);

foreach $p (@pronunciations){
    @current_pronunciation =  split(/\ +/,$p);
    printf("\t\t%s\n",branching_node(2,scalar @current_pronunciation));
 
    foreach $phone (@current_pronunciation){
        $index = $phoneTable{$phone};
        if(! defined $index){ print STDERR "$phone appears in $dict but not in $phoneme_map\n"; exit(-1);}
        printf("\t\t\t-1 %d %% %s\n",$index,$phone);
    }
}

close(DICT);

sub branching_node{
    my $query_parent = shift;
    my $num_children = shift;

    my $node = "$query_parent $num_children";
    if($num_children > 3){
        $node .= " 0 ... ";
        $node .= $num_children-2;
    }else{
        for($i = 0; $i+1 < $num_children; $i++){
            $node .= " $i";
        }
    }
    $node .= " default";

    return $node;
}

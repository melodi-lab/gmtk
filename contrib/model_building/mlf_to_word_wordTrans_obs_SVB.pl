#!/usr/bin/perl

## @file
# Modified version of genCounterToWordMapDTs.pl included with the gmtk Aurora tutorial

sub usage {
 printf STDERR ("Description: Generate files of word and wordTransition observations for SVitchboard utterances\n" .
  "Usage: mlf_to_word_wordTrans_obs.pl [options]\n" .
  "   -w file     word_map file, each line consisting of '<index> <word>'\n" .
  "   -f file     HTK file (MLF) containing the training data\n" .
  "   -m N        maximum number of utterances to traverse in MLF file\n" .
  "   -e ext      file extension of items in the MLF file\n" .
  "   -u file     file containining list of utts to do\n" .
  "   -d dir      output directory in which to put the observation files\n");
}

#bail out if there are no arguments
if ( $#ARGV < 0 ) {
    &usage;
    exit(-1);
}

# parse arguments
$PLP_dir="/export/ws06afsr/data/SVB/PLP/ICSI/sw1-norms";
$file="/dev/stdin/";
$ext="lab";
$maxuts = 1e8;
$eou = ".";#end of utterance, could make this an arg

while ( $#ARGV >= 0 ) {
    $opt = shift(@ARGV);
    $opt =~ tr/A-Z/a-z/; # ignore case
    if ( $opt =~ /^-w/ ) {#the word-map file
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
    } elsif ( $opt =~ /^-u/ ) {#the MLF file
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        $utt_file = $arg;
    } elsif ( $opt =~ /^-d/ ) {#the output dir
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        $outdir = $arg;
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


#read in the utt list 
open(U,"$utt_file") || die("cannot open file $utt_file");
foreach $utt (<U>) {
    chomp($utt);

#    if ($utt =~ /^(sw\d\d\d\d[AB]-)ms98-\w-(\d\d\d\d-\d+)$/) {
#	$utt = $1 . $2;
#    } else {
#	die "Incorrectly formatted utt tag:  $utt\n";
#    }

    $do_utt{$utt} = 1;
}
close(U);

#open the MLF file
open(F,"$file") || die("cannot open file $file");

#move the filepointer past the header
while(<F>) {
    chomp;
    # print "read ",$_,"\n";
    if ($_ eq "#!MLF!#") {
	last;
    }
}

# The result of this while loop is that, upon termination, $numUts holds the lesser value of $maxuts and the number of utterances in the MLF file (as determined by the number of file extensions)
$numUts = 0;
while (($record = <F>)) {
   if($record =~ /\.$ext\"$/){
       $numUts++;
   }
   if ($numUts > $maxuts) {
       $numUts = $maxuts;
       last;
   }
}


# go to beginning
seek(F,0,0);

print STDERR "Generating from file $file\n\n";

#skip the header again
while(<F>) {
    chomp;
    if ($_ eq "#!MLF!#") {
	last;
    }
}


%wordTable = load_vocab($vocab);


#iterate through the MLF file 
for ( $utt=0; $utt<$numUts; $utt++) {
    $_ = <F>;

    #grab the sentenceid and create the output .wd and .wdTr files
    ($sentenceid) = /\/([^\/]*)\.$ext/;

    if ($sentenceid =~ /^(sw\d\d\d\d[AB])-(\d\d\d\d-\d+)$/) {
	$sentenceid = $1 . "-ms98-a-" . $2;
    } else {
	die "Incorrectly formatted sentence id in MLF file:  $sentenceid\n";
    }
    
    #debug
    print STDERR "$utt $sentenceid\n";

    #this while loop will go through a sentence and populate the 
    #$words, $start_times, and $end_times arrays

    $wordnum = 0;
    @words = ();
    @start_times = ();
    @end_times = ();

    while (1) {
	$line = <F>;
	chomp($line);

	if ($line eq $eou) {  # done with utt--now write the obs files

	    if (! $do_utt{$sentenceid}) {
		print STDERR "Skipping utt $sentenceid\n";
		last;
	    }

	    if ($sentenceid =~ /^(sw\d\d\d\d[AB])+/) {
		$conv_id = $1;
	    } else {
		die "Invalid sentence id:  $sentenceid\n";
	    }

	    $return = system("mkdir -p $outdir/$conv_id");
	    if ($return) {
		die "Couldn't make directory $outdir/$conv_id\n";
	    }

	    open(WD, ">$outdir/$conv_id/$sentenceid.wd") || die "Couldn't open file $outdir/$sentenceid.wd for writing.\n";
	    open(WDTR, ">$outdir/$conv_id/$sentenceid.wdTr") || die "Couldn't open file $outdir/$sentenceid.wdTr for writing.\n";

	    # get the number of frames in the PLP file; this number is off from
	    # the MLF file by a non-constant amount, so we just match it
	    $num_frames_PLP = `HList -h -z $PLP_dir/$conv_id/$sentenceid.plp | grep 'Num Samples' | awk '{print \$3}'`;
	    chomp($num_frames_PLP);
	    
	    print STDERR "num_frames_PLP = $num_frames_PLP\n";

	    $total_num_frames_out = 0; 
	    foreach $i (0..$#words) {

		# round to the nearest start & end frames using sprintf; 
		# the reason for the "+1" is that sprintf apparently rounds 
		# .5 down to 0
		$start_frame = sprintf("%.0f", ($start_times[$i]+1)/100000);
		$end_frame = sprintf("%.0f", ($end_times[$i]+1)/100000);

		$num_frames = $end_frame - $start_frame;

		# if we're in the last word, match the number of frames in the PLP file
		if ($i == $#words) {
		    print STDERR "NOTE:  Changing num_frames from $num_frames to ", $num_frames_PLP-$total_num_frames_out, ".\n";
		    $num_frames = $num_frames_PLP - $total_num_frames_out;
		}

		print STDERR "num_frames for word $words[$i] = $num_frames\n";

		$total_num_frames_out += $num_frames;

		# sentence-initial and -final silences are mapped to <s>, </s>
		# all other silences are mapped to previous word
		if ($words[$i] eq "sil") {
		    if ($i == 0) {
			$words[$i] = "<s>";
		    } elsif ($i == $#words) {
			$words[$i] = "</s>";
		    } else {
			$words[$i] = $words[$i-1];
		    }
		}

		if (!defined($wordTable{$words[$i]})) {
		    die "ERROR: undefined word $words[$i] in utterance $sentenceid\n";
		}
		
		if ($num_frames > 0) { # this check because of sw2237B-0023-1
		    foreach $f (1..$num_frames-1) {
			print WD "$wordTable{$words[$i]}\n";
			print WDTR "0\n";
		    }
		    
		    # last frame
		    print WD "$wordTable{$words[$i]}\n";
		    print WDTR "1\n";
		}
		    
	    }

	    close(WD);
	    close(WDTR);
	    last;
        }

	($start, $end, $word) = split(/\s+/, $line);

	$word =~ tr/[A-Z]/[a-z]/;
	
	$words[$wordnum] = $word;
	$start_times[$wordnum] = $start;
	$end_times[$wordnum] = $end;
	$wordnum++;
    }
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
        ($index, $word) = split(/\s+/, $_);
	$word =~ tr/[A-Z]/[a-z]/; # lowercase everything

        $table{$word} = $index++;

    }
    close(VOCAB);

    return %table;
}

#checks that the argument $arg for option $opt is a positive integer
sub checknnint {
    ($opt,$arg) = @_;
    if ( $arg !~ /^[0-9]+$/ ) {
        print STDERR ("Error: Value \"", $arg, "\" invalid for option ",
                      $opt, " (positive integer expected)\n");
        &usage; exit(-1);
    }
}



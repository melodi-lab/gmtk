#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Lazier;
use List::Util qw(sum max);
use Data::Dumper;

## @file
# @author rewritten bigtime by Arthur Kantor
# @date 6/18/08

sub parseOpts{
    my $usage = <<'EOS';
$prog [--help] [options]

#################################################################
##
## generate_phoneStateDT.pl
##
## This script takes in a dictionary file and generates a DT 
## for the phoneState variable with conditioning parents word, 
## pronVariant, and subWordState.
##
#################################################################

the program checks that the phones used in dictFile are present in phoneStateFile, 
and that words in dictFile are present in word_map.

$opts
EOS


	my $optSpec = [
		['help|h', 
		undef,
		'print this help message and exit'
		],
		['dictFile|d=s',
		undef,
		'dictionary filename, one pronunciation variant per line, formatted as "word phone1 phone2 ...".  '.
		'Multiple pronunciation for the same word will have an identical word at the beginning of the line.',
		'required fileExists'
		],
		['phoneStateFile|p=s',
		undef,
		'file with one phone per line, formatted as "<int> <phone>", where the <int> specifies the number of sub-phone states',
		'required fileExists'
		],
		['wordmap|w=s', 
		undef, 
		'A file with one word per line, mapping a word to an index which will represent the word in the DT.  '.
		'It is assumed that wordmap is ordered, with first word in file represented by index 0.',
		'required fileExists'
		],
		['subphoneStateDt|t=s', 
		'word_pronunciation_wordStateCounter_2_state.dt', 
		'A file into which the word_pronunciation_wordStateCounter_2_state DT will be written .',
		'required'
		],
		['wordTransitionDt|c=s', 
		undef, 
		'A file into which the word_pronunciation_wordStateCounter_stateTransition_2_wordTransition DT will be written .'
		],
		['endOfWordPhone|e=s', 
		undef, 
		'If the option is specified, this phone will be appended to every word.  This phone must be present in the '.
		'phoneStateFile.  Use this to specify interword pauses. Also useful for replacement of the  wordTransition DT.'
		],
		['singlePronunciation|s', 
		undef, 
		'If the option is specified, use the first pronunciation for a word in dictFile. The generated decision trees do not'.
		'have pronunciation as a dependency.'
		],
	];
	
	my %opts;
	exit(-1) if (!parseOptions(\%opts, $optSpec));
	if($opts{'help'}){
		standardUsage($optSpec, $usage);
		exit(0);
	}	
	
	my @err=validateOptions(\%opts, $optSpec);
	if (@err){
		print "@err";
		standardUsage($optSpec, $usage);
		exit(-1);
	}

    return %opts;
}



# read in phoneme map table
my %opts = parseOpts();
my ($psr, $sp2ir)=readPhoneStateFile($opts{'phoneStateFile'});
my %phoneSets=%$psr;
my %subPhone2index=%$sp2ir;


if ($opts{'endOfWordPhone'} ){
	if(!$phoneSets{$opts{'endOfWordPhone'}}){
		die "End of word marker $opts{'endOfWordPhone'} is not in the phoneState file $opts{'phoneStateFile'}"; 
	}
	my @eowSet = @{$phoneSets{$opts{'endOfWordPhone'}}};
	print "endOfWordPhone/endOfWordPhoneIndex (END_OF_WORD_MARKER): $eowSet[$#eowSet] $subPhone2index{$eowSet[$#eowSet]}\n";
}


my $w2ir = readWordMap($opts{'wordmap'});
my %word2index =%$w2ir;
#print keys %word2index;
my $numWords=scalar(keys %word2index);
print "$numWords total words.\n";


# keep track of mapping from wd nums to word labels
my %index2word= reverse(%word2index);


my ($ptr,$maxSubPhoneStatesInWord) = readDict($opts{'dictFile'});
my %phonTranscriptions=%$ptr;
print "max subphone states in a word (WORDSTATE_COUNTER_CARD) $maxSubPhoneStatesInWord\n";


open(SUBPHONEDT,">$opts{'subphoneStateDt'}") || die "Couldn't open dict file >$opts{'subphoneStateDt'} for reading.  Exiting.\n";
if($opts{'wordTransitionDt'}){
	open(WORDTRANSITIONDT,">$opts{'wordTransitionDt'}") || die "Couldn't open dict file >$opts{'wordTransitionDt'}. Exiting.\n";
}
else{
	open(WORDTRANSITIONDT,'>/dev/null') || die "Couldn't open dict file >/dev/null for reading.  Exiting.\n";
}

my $optString = Dumper(\%opts);
my @optList = split("\n", $optString);
$optString = join ("\n",map {"% $_"} @optList);
my $header =<< "END_OF_HEADER";
% generated with: $0
% parameters:
$optString

1 % number of DTs in this file

0 % DT number

END_OF_HEADER


print SUBPHONEDT $header;
print WORDTRANSITIONDT $header;
if ($opts{'singlePronunciation'}){
	print SUBPHONEDT <<EOS;
word_wordStateCounter_2_state % DT name
2 % num parents
EOS

	print WORDTRANSITIONDT <<EOS;
word_wordStateCounter_stateTransition_2_wordTransition % DT name
3 % num parents
EOS
}
else{
	print SUBPHONEDT <<EOS;
word_pronunciation_wordStateCounter_2_state % DT name
3 % num parents
EOS

	print WORDTRANSITIONDT <<EOS;
word_pronunciation_wordStateCounter_stateTransition_2_wordTransition % DT name
4 % num parents
EOS
}



my @buf  = ('0' , $numWords , makeRange($numWords-1) , 'default');
print SUBPHONEDT "@buf\n";
print WORDTRANSITIONDT "@buf\n";
my $wordStateCounterParent = defined($opts{'singlePronunciation'}) ? 1 : 2;
my $stateTransitionParent = $wordStateCounterParent + 1;
foreach my $wdIndex (0..$numWords-1) {
	my $tabs=1;
	die "word $index2word{$wdIndex} at index $wdIndex has no phonetic transcription in dictionary" if (!defined($phonTranscriptions{$wdIndex}));
	my @pronSet = @{$phonTranscriptions{$wdIndex}};
#    print "word = $word, num_phones = $num_phones{$word}\n";
    my $label = $index2word{$wdIndex};
    $label =~ s/\'/\_APO_/g;

	@buf = ("\t"x$tabs); 
	if(!$opts{'singlePronunciation'}){
	    @buf =(@buf, " 1" ,scalar(@pronSet), makeRange(@pronSet-1),'default');
	    $tabs=2;
	}
    @buf =(@buf, "% wd index $wdIndex, word $label");
	print SUBPHONEDT "@buf\n";
	print WORDTRANSITIONDT "@buf\n";


    foreach  (@pronSet) {
		my @pron = @$_;
		my $spsCount=scalar(@pron);
		@buf = ("\t"x$tabs,"$wordStateCounterParent $spsCount ", makeRange($spsCount-1), "default");
    	print SUBPHONEDT "@buf\n";
		foreach (@pron) {
			print SUBPHONEDT ("\t"x($tabs+1))."-1  $subPhone2index{$_} % $_\n";
		}

		if($spsCount==1){#there is only subphone state: a word transition happens when a subphone state transition happens
			print WORDTRANSITIONDT ("\t"x$tabs) ." -1 {p$stateTransitionParent}\n";
		}
		elsif($spsCount>1){
    		print WORDTRANSITIONDT ("\t"x$tabs) ." $wordStateCounterParent 2 $#pron default\n";
    		print WORDTRANSITIONDT ("\t"x($tabs+1))." -1 {p$stateTransitionParent}\n";
    		print WORDTRANSITIONDT ("\t"x($tabs+1))." -1 0\n";
		}
		else{
			die "Should not be here. pron=@pron";
		}
    }
}

close(WORDTRANSITIONDT);
close(SUBPHONEDT);

#the range specification in GMTK seems to be buggy: e.g. 0 ... 0 does not work, 
#this returns an acceptable range of N consecutive elements starting at 0.
sub makeRange{
	my $N= shift;
	if ($N < 0){
		die "cannot make range into negative numbers";
	}
	elsif($N == 0){
		return '';
	}
	elsif($N == 1){
		return '0';
	}
	else{
		return '0 ... '.($N-1);
	}
}


sub readPhoneStateFile{
    my $phoneStateFile = shift;
    my %phoneSets=();
    my %subPhone2index=();
    open(PH, $phoneStateFile ) || die "Couldn't open phoneme map file $phoneStateFile  for reading.  Exiting.\n";
    my $subPhoneIdx=0;
    while (<PH>) {
        chomp;
        my($phone, $states)=split;
        my @subPhones= map "$phone$_", (0 .. $states-1);
        $phoneSets{$phone} = \@subPhones;
        foreach (@subPhones){
            $subPhone2index{$_}=$subPhoneIdx++;
        }
    }
    print "total phones: ".scalar(keys %phoneSets)."\n";
    print "total subphone states (STATE_CARD): $subPhoneIdx\n";
    close(PH);
    return (\%phoneSets, \%subPhone2index);
}


sub readWordMap{
    my $wordMapFile = shift;
    my %word2index;
    open (WORDMAP, "<$wordMapFile") || die "cannot open <$wordMapFile";
    while (<WORDMAP>) {
        chomp;
        s/^ *//;
        s/ .*//;
        die " wordmap has duplicate indexes for word $_" if($word2index{$_});
        $word2index{$_} = $numWords++;
    }
    close WORDMAP  || die "cannot close <$opts{'wordmap'}";
    #print Dumper(\%word2index);
    
    return \%word2index;
}


# get pron of each word from dict & num phones per word, and
# build (word, position) --> (phoneme) mapping table
#map from wdIndex => [<pron1>, <pron2>,....]
#where <pron1> = [subphone1, subphone2, ...]
sub readDict{
    my $dictFileName = shift;
    open(DICT, $dictFileName) || die "Couldn't open dict file $dictFileName for reading.  Exiting.\n";
    my %phonTranscriptions=();

    my $maxSubPhoneStatesInWord=0;
    foreach (<DICT>) {
        chomp;
    #    print $line, "\n";
        my ($wd, @phones) = split;
        if (!defined($word2index{$wd})){
            #die "The word $wd on line $. of $opts{'dictFile'} is not in the wordmap file" 
            next;
        }
        push @phones, $opts{'endOfWordPhone'} if ($opts{'endOfWordPhone'});
        foreach (@phones){
            die "The phone $_ on line $. of $opts{'dictFile'} is not in the phoneStateFile file" if (!$phoneSets{$_});
        }
        my $wdIndex = $word2index{$wd};
        #print Dumper(\%word2index);
        my @subPhoneStates = map {@$_} @phoneSets{@phones};
        $maxSubPhoneStatesInWord=max($maxSubPhoneStatesInWord, scalar(@subPhoneStates));
        $phonTranscriptions{$wdIndex} = [] if(!defined($phonTranscriptions{$wdIndex}));
        push @{$phonTranscriptions{$wdIndex}}, \@subPhoneStates;
    
    #	print "wd = $wd, wd num = $wd_num, pron = $pron_num, pos = $numWords, phoneme = ", $phoneTable{$phones[$numWords]}, "\n";
    }
    close DICT;
    
    return (\%phonTranscriptions, $maxSubPhoneStatesInWord);
}
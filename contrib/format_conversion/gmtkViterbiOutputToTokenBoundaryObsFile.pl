#!/usr/bin/env perl
use warnings;
use strict;

use Carp;
use Getopt::Long;
use File::Path;
use File::Basename;
use Getopt::Lazier;
use OS::Util  qw(	readLines writeLines tsSystem );
use Fatal qw(:void open close); #die on failed open and close if the return values are not checked

## @file
#  @author
#  @date Arthur Kantor 9/6/09

#parses options and checks if the files/dir specified as parameters exist
sub parseOpts{

	my $usage = <<'EOS';
Usage: 
$prog [--help] [options]  

Takes a GMTK viterbi output file and converts it 
to a token and boundary observation file in ascii format.  Tokens are usually words.

Options:

$opts

EOS

    my $optSpec = [
		['help|h', 
		 undef,
		 'print this help message and exit'],
		['alignmentFile|a=s', 
		 undef, 
		 'A file in the format of the output of gmtkViterbi.',
		 'required fileExists'],
		['obsFile|o=s', 
		 undef,
		 'Where to write observations.',
		 'required'],
		['varMap|v=s', 
		 undef,
		 'The varmap file used by gmtkViterbiNew during forced aligment.',
		 'required'],
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


{ #brackets to keep my variables private from subroutines 

	#parse and check the arguments
    my %opts=parseOpts();	
	my $cmdlineArgs = join(' ', @ARGV);
	print "Script started on " . localtime(time)." with args: $cmdlineArgs\n";
	
	my %varMap;
	if($opts{'varMap'}){
    	my @varMapLines = readLines($opts{'varMap'});
    	print "Read ".scalar(@varMapLines)." words from $opts{'varMap'}\n";
    	chomp @varMapLines;
    	%varMap=reverse(map {split} @varMapLines);
	}
	
	
	open(AL,"<$opts{'alignmentFile'}");
	open(OBS,">$opts{'obsFile'}");
	my @curUtt = readUTT(\*AL);
	my $uttCounter =0;
	while(@curUtt){
	    my $uttId = shift(@curUtt);
       	print '.' if ($uttCounter % 10000==0);
	    
	    #print "$fname, $uttStartFrame, $uttEndFrame $uttEndSec $uttStartSec $uttPrefix\n";
	    if (@curUtt){
			my ($word, $startFrame, $endFrame);
			while (@curUtt){
				($word, $startFrame, $endFrame)=splice(@curUtt,0,3);
				 my $wordId=$varMap{$word};
				#print CURCONV pack("NN",($varMap{$word}, 0))x($endFrame-$startFrame);
				#print CURCONV pack("NN",($varMap{$word}, 1));
				for my $i  ($startFrame..($endFrame-1)){
					print OBS "$uttCounter $i $wordId 0\n";
				}
				print OBS "$uttCounter $endFrame $wordId 1\n";
				
			}
           $uttCounter++;
           
        }
	    else{
			die "Empty Utterance $uttId\n";
	    }
        
	    @curUtt = readUTT(\*AL);
	}
	close(AL);
    close(OBS);
}

##@fn
# @return a list of form (uttId, word1, word1StartFrame, word1EndFrame, word2, ...)
sub readUTT{
	my ($ALfh) = @_;
	
	my @ret=();
	$_ = <AL>;
	if (!$_){
	    return @ret;
	}
    elsif(/^Segment ([0-9]+)/){
        push @ret, $1;
    }
    my $prevLine=tell $ALfh;
    
	while(<AL>){
		chomp;
		if(/(.*) \(([0-9]+)-([0-9]+)\)$/){
		    push @ret, ($1, $2, $3);
		}
		elsif(/^Segment ([0-9]+)/){
			seek $ALfh, $prevLine, 0;
			last;
		}
		#else{
		#    die "unrecognized line \"$_\" in decode file\n";
		#}
		$prevLine = tell $ALfh;
	}
	
	return 	@ret;
}








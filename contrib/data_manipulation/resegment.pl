#!/usr/bin/env perl
use warnings;
use strict;

use Carp;
use Getopt::Long;
use File::Path;
use File::Basename;
use Getopt::Lazier;
use OS::Util  qw(	readLines writeLines tsSystem );

## @file
#  @author
#  @date Arthur Kantor 10/10/08

#prints the usage information
sub usage{
	my ($optSpec) = @_;	
    my $usage = <<'EOS';

Usage: 
$prog [--help] [options]  

resegments the utterances so that there is a single utterance per 
each emitted token in the alignmentFile.  The original utterance boundaries
will also be present in the resulting boundaries.

if --obsDir is specified the output is written in big-endian (network) format

FIXME use the new style (NIST swb) utterance ids

$opts

Return values:
	0  success
	-1 error
EOS

	return $usage
}

#parses options and checks if the files/dir specified as parameters exist
sub parseOpts{
    my $optSpec = [
		['help|h', 
		 undef,
		 'print this help message and exit'],
		['outTranscriptionFile|t=s',
		 undef, 
		 'The file into which store the transcriptions of the out utterances, each line containing "uttId transcription".  REQUIRED.'],
		['alignmentFile|a=s', 
		 undef, 
		 'A file in the format of the output of gmtkViterbi. REQUIRED.'],
		['inUttIdsFile|u=s', 
		 undef,
		 'The file with utterance ids, one per line, uttId of form SOMENAME_START_END, with START and END being 2-digit-after-decimal-point floats of utterance start and end in seconds. REQUIRED.'],
		['inHtkFeatFile|i=s', 
		 undef,
		 'An observation file in htk scp format corresponding to original partitioning(a file listing the actual feature files, with each feature file per utterance). REQUIRED.'],
		['outHtkFeatFile|o=s', 
		 undef,
		 'An observation file in htk scp format to which the resulting partition will be written. REQUIRED.'],
		['obsDir|g=s', 
		 undef,
		 'The dir where the actual data (pointed to by outHtkFeatFile) will be written.  The dir will be created if missing.  In this case, varMap is required, and outHtkFeatFile and outTranscriptionFile contain the full utterances.  They may be a subset of the in-utterances because some utterances could not be aligned to their transcriptions '],
		['varMap|v=s', 
		 undef,
		 'The varmap file used by gmtkViterbiNew during forced aligment. Required with --obsDir.'],
	];


    my %opts;
	exit(-1) if (!parseOptions(\%opts, $optSpec));
	
	if($opts{'help'}){
		standardUsage($optSpec, usage());
		exit(0);
	}	
	
	my $parseError=0;
	foreach ('outTranscriptionFile', 'alignmentFile', 'inHtkFeatFile', 'outHtkFeatFile', 'inUttIdsFile'){
	    if (!$opts{$_}) {
    		print "--$_ option required\n";
	        $parseError=1;
	    }
	}
	
	if($parseError){
		standardUsage($optSpec, usage());
		exit(-1);
    }
    
	checkFileExists($opts{'alignmentFile'});
	checkFileExists($opts{'inHtkFeatFile'});
	checkFileExists($opts{'inUttIdsFile'});
    checkFileExists($opts{'varMap'}) if($opts{'varMap'});
    return %opts;
}






# die if dir does not exit, return otherwise
sub checkDirExists{
	if (! -d $_[0]){
		die "$_[0] does not exist.\n";
	}
}

# die if file does not exit, return otherwise
sub checkFileExists{
	if (! (-e $_[0] && -f $_[0]) ){
		die "$_[0] does not exist.\n";
	}
}

{ #brackets to keep my variables private from subroutines 

	#parse and check the arguments
    my %opts=parseOpts();	
	
	print "creating $opts{'obsDir'} directory heirarchy\n";
	for my $i (0 .. 116){
	    tsSystem(sprintf("mkdir -p $opts{'obsDir'}/%03d",$i));
	}
	
	my %varMap;
	if($opts{'varMap'}){
    	my @varMapLines = readLines($opts{'varMap'});
    	print "Read ".scalar(@varMapLines)." words from $opts{'varMap'}\n";
    	chomp @varMapLines;
    	%varMap=reverse(map {split} @varMapLines);
	}
	
	my @inScp = readLines($opts{'inHtkFeatFile'});
	print "Read ".scalar(@inScp)." utterances from $opts{'inHtkFeatFile'}\n";
	my @utts = readLines($opts{'inUttIdsFile'});
	print "Read ".scalar(@utts)." utterances from $opts{'inUttIdsFile'}\n";
	die "Number of utterance ids and scp rows must be equal" if (@inScp != @utts);
	
	open(AL,"<$opts{'alignmentFile'}") || die "cannot open <$opts{'alignmentFile'}";
	open(HTK,">$opts{'outHtkFeatFile'}") || die "cannot open >$opts{'outHtkFeatFile'}";
	open(TRANS,">$opts{'outTranscriptionFile'}") || die "cannot open >$opts{'outTranscriptionFile'}";
	#if($opts{'obsDir'}){
	#    open(OBS,">$opts{'obsDir'}") || die "cannot open >$opts{'obsDir'}";
	#}
	my @curUtt = readUTT(\*AL);
	my $uttCounter =0;
	my $prevOutFileName='';
	while(@curUtt){
	    my $uttId = shift(@curUtt);
       	print '.' if ($uttCounter % 10000==0);
	    
	    my $uttScpLine=$inScp[$uttId];
	    chomp($uttScpLine);
	    my ($fname, $uttStartFrame, $uttEndFrame) = split(/[\[:\]]/,$uttScpLine);
	    my ($uttIdPrefix) = fileparse($fname, '\..*');
	    
	    my $uttIdLine=$utts[$uttId];
	    chomp($uttIdLine);
	    my @recs = split('_',$uttIdLine);
	    #FIXME use the new style (NIST swb) utterance ids
	    my $uttEndSec= pop(@recs);
	    my $uttStartSec= pop(@recs);
	    my $uttPrefix=join('_',@recs);
	    my $d=int($recs[2]/100);
	    my $outFname=sprintf("$opts{'obsDir'}/%03d/$uttPrefix.dat", $d);
	    #print "$outFname\n";
	    die "utterance id $uttIdLine is inconsistant with scp line $uttScpLine at line $uttId\n" if ($uttPrefix ne $uttIdPrefix);
	    
	    #switch to a new conversation-side file if needed
        if ($outFname ne $prevOutFileName && $opts{'obsDir'}){
            if ($prevOutFileName){
                close(CURCONV) || die "cannot close >$prevOutFileName";
            }
            open(CURCONV, ">$outFname") || die "cannot open >$outFname";
            $prevOutFileName=$outFname;
        }
	    
	    #print "$fname, $uttStartFrame, $uttEndFrame $uttEndSec $uttStartSec $uttPrefix\n";
	    if (@curUtt){
           my ($word, $startFrame, $endFrame);
	       if($opts{'obsDir'}){
               print HTK "$uttScpLine\n";
               print TRANS "$uttIdLine ";
               while (@curUtt){
                    ($word, $startFrame, $endFrame)=splice(@curUtt,0,3);
                    print TRANS "$word ";
                    print CURCONV pack("NN",($varMap{$word}, 0))x($endFrame-$startFrame);
                    print CURCONV pack("NN",($varMap{$word}, 1));
                    #for my $i  ($startFrame..($endFrame-1)){
                    #    print OBS "$uttCounter $i $varMap{$word} 0\n";
                    #}
                    #print OBS "$uttCounter $endFrame $varMap{$word} 1\n";
                    
                }
                print TRANS "\n";
	       
	       }
	       else{
                while (@curUtt){
                    ($word, $startFrame, $endFrame)=splice(@curUtt,0,3);
                    my ($wStartFrame, $wEndFrame)=($uttStartFrame+$startFrame, $uttStartFrame+$endFrame);
                    print HTK "$fname\[$wStartFrame:$wEndFrame\]\n";
                    my ($wStartFrameSec, $wEndFrameSec)=($uttStartSec+$startFrame/100.0, $uttStartSec+$endFrame/100.0);
                    print TRANS sprintf("%s_%.2f_%.2f %s\n",$uttIdPrefix, $wStartFrameSec, $wEndFrameSec, $word);
                }
           }
           $uttCounter++;
           
           die "Utterance (in-id $uttId) is ". ($uttEndFrame-$uttStartFrame+1)." frames long  but its component words are ". ($endFrame+1)." frames long."  if ($uttEndFrame-$uttStartFrame != $endFrame);
        }
	    else{
	       my $uttFrames = $uttEndFrame-$uttStartFrame+1;
	       print "WARNING: Utterance Id $uttId $uttIdPrefix\[$uttStartFrame:$uttEndFrame\] with $uttFrames frames could not be aligned. (too few frames for the given transcription?)\n";
	       if($opts{'obsDir'}){
	           print "appending $uttFrames dummy frames to $outFname to keep the word alignment indexes in sync with the plp files.\n";
	           #999999 is an illegal value for wordTransition and should cause an error if this utterance is ever used
	           print CURCONV pack("NN",(0, 999999))x$uttFrames;
	       }
	    }
        
	    @curUtt = readUTT(\*AL);
	}
	close(AL) || die "cannot close <$opts{'alignmentFile'}";
	close(HTK) || die "cannot close >$opts{'outHtkFeatFile'}";
	close(TRANS) || die "cannot close >$opts{'outTranscriptionFile'}";
	if($opts{'obsDir'}){
        close(CURCONV) || die "cannot close >$prevOutFileName";
        #close(OBS) || die "cannot close >$opts{'outTranscriptionFile'}";
	}
}

sub ecSystem{
	print "executing '@_'\n";
	system(@_)==0 or confess("error with $?");
}

#create skip and deslenfile lists
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








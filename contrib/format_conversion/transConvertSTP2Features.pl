#!/usr/bin/perl

## @file
#  @author Ari Bezman
#
#This script converts STP transcripts to Wavesurfer-type feature tier files.
#It can handle either STP_agreement or STP CD transcripts (some differences in format between these two).
#To actually map STP phones (+ diacritics) to features, it relies on a lookup table/mapping file.  The one used to convert STP  to WS06AFSR is called agreement_mappings.
#Output is to n .feat files, where n is the number of features, and .feat is each feature name.  The output files are "chunked", where all the timestamps line up.  This means that each feature tier has the same boundaries as all others, and, consequently, that it's possible to have adjacent identical values for some feature tier.  transProc.pl can be used to fix this, if desired.


#use strict;
use warnings;


if (@ARGV < 1) {
    print STDERR "Usage:  $0 phone2feature-file input-file";
    print STDERR "      output will be to input-file.pl1, input-file.pl2, etc.\n";
    exit;
}

# read in phone2articulator mapping
open(PH2ARTIC, shift @ARGV) || die "Couln't open phone2articulator file for reading.  Exiting.\n";

my @feature_names;
my %map;

while (my $line = <PH2ARTIC>) {
    chomp($line);
    if ($line =~ "HEADER"){  #line is the HEADER line, containing the list of features
	@feature_names = split(/\t/,$line);
	shift @feature_names; #get rid of the word HEADER
    }
    elsif (($line =~ /^\#.*$/) || ($line =~ /^\s*$/)) {   # line is a comment or is empty--skip it
	next;
    }
    else { #line is a phone to feature mapping line, so read in the mapping
	my @mapping_values = split(/\t/, $line);
	my $phone_name = shift @mapping_values; #The first entry in the line is actually the STP phone to be mapped
	unless (@mapping_values == @feature_names){  #There MUST be as many features for each phone as there were in the header
	    die ("ERROR: Incorrect number of features in line\n $line\n");
	}
	for (my $i = 0; $i < @feature_names; $i++){
	    $map{$phone_name}{$feature_names[$i]} = $mapping_values[$i];
	}
    }
}
	
close(PH2ARTIC);

while(my $in_name = shift @ARGV){ #For each input file
    open(IN_FILE, $in_name);

    #Open the corresponding output files
    my %out_files;
    foreach (@feature_names){
	my $out_name = $in_name . "." . $_;
	open(OUT."_".$_, ">".$out_name);
	$out_files{$_} = "OUT_$_";
    }

    #Dump the input file header
    my $line = <IN_FILE>;
    until($line =~ /END OF HEADER|^\#$/){
	$line = <IN_FILE>;
    }

    #Figure out what sort of file it is based on the header
    #In the case of Intra-STP agreement files, set to AGR
    #For the other files, set to "STP_CD"
    my $filetype;
    if($line =~ /END OF HEADER/){
	$filetype = "AGR";
    }
    else {$filetype = "STP_CD"};

    my %phone_features;
    #The following variables are for tracking the previous line when its values
    #are unknown until the following line is read.
    #e.g [hh]
    my %prev_phone_features;
    my $prev_phone;
    my $prev_start = 0;
    my $prev_end = 0;

    #Read in each line (phone+diacritic), and parse it
    while ($line = <IN_FILE>) {
	chomp($line);
	my ($start, $end, $phone);
	if($filetype eq "AGR"){
	    ($start, $end, $phone) = split(' ', $line);
	    
	    #Convert from ms to sec
	    $start /= 1000;
	    $end /= 1000;
	}
	elsif ($filetype eq "STP_CD"){
	    my @line = split(' ', $line);
	    ($end, $phone) = ($line[0], $line[2]);
	    $start = $prev_end || 0;
	}
	else { die ("File type not STP_CD or AGR")}
	
	if ($start >= $end){
	    warn ("WTF! Start $start is bigger than end $end in \n \t $line \n");
	}
	
	#First, the diphthongs
	#If the phone1 and phone2 mappings exist, it's a diphthong. So use them, splitting in the middle.
	if ($map{$phone."1"} && $map{$phone."2"}){
	    my $midpoint = ($start + $end) / 2;
	    my $phone1 = $phone."1";
	    my $phone2 = $phone."2";

	    #Print out to the files here
	    foreach (@feature_names){
		print {$out_files{$_}} "$start $midpoint $map{$phone1}{$_}\n";
		print {$out_files{$_}} "$midpoint $end $map{$phone2}{$_}\n";
	    }
	}

	#If the filetype is STP_CD, stops similarly need to be split into closures and bursts, as those haven't been distinguished
	#Note: it is assumed (and confirmed observationally) that STP_CD files NEVER use stop_cl or stopcl phones.
	elsif (($filetype eq "STP_CD") && ($phone =~ /^[pbtdkg]/)){
	    my $burst = $phone;
	    #The first character indicates the type of stop
	    my $root = substr($phone, 0, 1);
	    #Diacritics, like _vl, etc.
	    my $diacritics = substr($phone, 1, length($phone)); 

	    my $closure = $root . "cl" . $diacritics;

	    #Try and form a closure with all the diacritics
	    #Then, failing that, with only an initial subset of them
	    until ($map{$closure}){
		#If $diacritics is empty AND we're in this loop,
		#there is actually no mapping for $root."cl",
		#i.e. no mapping for any closure for this stop.

		#So, give up and set the closure to match the burst
		if ((length $diacritics) == 0){
		    $closure = $burst;
		    warn ("Could not find a closure for $phone from $line\n");
		    last;
		}
		#Lose one diacritic character and try, try again.
		else {
		    chop $diacritics;
		    $closure = $root . "cl" . $diacritics;
		}
	    }
	    
	    #It's possible, but unlikely, that, if the full phone is something like p_epi, there exists a p_epicl
	    #So, we look for it, but then start giving up diacritics until we get down to pcl
      	    if ($map{substr($phone, 0, 1)."cl"}){
		$closure = substr($phone, 0, 1)."cl";
	    }
	    else{
		warn "Could not split $phone in $line into closure/burst \n";
	    }
	    
	    #FIXME: Should this really be halfway?  Are bursts longer or shorter than closures?
	    my $midpoint = ($start + $end) / 2;
	    #Print out to the files here
	    foreach (@feature_names){
		print {$out_files{$_}} "$start $midpoint $map{$closure}{$_}\n";
		print {$out_files{$_}} "$midpoint $end $map{$burst}{$_}\n";
	    }
	}

	#hh (and its variants) is a special case, since some of its features
	#should be gathered from the following vowel, if any
	elsif ($phone =~ /^h$|^hv$|^hh/){
	    %prev_phone_features = %phone_features;
	    ($prev_phone, $prev_start, $prev_end) = ($phone, $start, $end);
	}

	#Only try to map phones that exist in the mapping table
	elsif ($map{$phone}){
	    foreach (@feature_names){
		$phone_features{$_} = $map{$phone}{$_};
		#If the previous phone was an [hh], and the current phone has a vow feature,
		#give the current phone's vow, frt, ht features onto the previous phone.
		#FIXME: Hardcoded strings
		 if ($prev_phone && ($prev_phone =~ /^h$|^hv$|^hh/) && ($phone_features{"vow"} !~ "N\/A")){
 		    $prev_phone_features{"vow"} = $phone_features{"vow"};
 		    $prev_phone_features{"frt"} = $phone_features{"frt"};
 		    $prev_phone_features{"ht"} = $phone_features{"ht"};
 		}

		#If the PREVIOUS phone was context-dependent, it hasn't been printed yet, so do that now
		if ($prev_phone && $prev_phone =~ /^h$|^hv$|^hh/){
		    print {$out_files{$_}} "$prev_start $prev_end $prev_phone_features{$_}\n";
		}
		print {$out_files{$_}} "$start $end $phone_features{$_}\n";
	    }	       
	}
	
	#If the phone couldn't be mapped, print asterisk filler into the output files
	#And an error to STDERR
	else {
	    warn ("Phone $phone in $line could not be found in the phone to feature mapping.\n");
	    foreach (@feature_names){
		print {$out_files{$_}} "$start $end \*\n";
	    }
	    next;
	}
	
    
	%prev_phone_features = %phone_features;
	($prev_phone, $prev_start, $prev_end) = ($phone, $start, $end);
	
    }
    
    #Close everything
    foreach (values %out_files){
	close $_;
    }
    close(IN_FILE);
}

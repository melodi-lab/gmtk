#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Lazier;
#use List::Util qw[min max sum reduce];
use Fatal qw(open close);
#use Data::Dumper;
use Carp;

#parses options and checks if the files/dir specified as parameters exist
sub parseOpts{
    my $usage = <<'EOS';
$prog [--help] [options]  
 
	Replaces keys in a field of one file with values from a second key-value map file

Options:

$opts
EOS
    
    my $optSpec = [
		["help|h", 
		 undef,
		 "print this help message and exit"],
		["map|m=s", 
		 undef, 
		 "The key-value map file, one key-value pair per line.",
		 'required fileExists'],
		['keyColumn|k=i',
		 '1',
		 'The column containing the key in the map file.',
		 'required'],
		['valueColumn|v=i',
		 '2',
		 'The column containing the value in the map file.',
		 'required'],
		["source|s=s", 
		 undef, 
		 "The source file.  If unspecified, uses STDIN",
		 'fileExists'],
		["output|o=s", 
		 undef, 
		 "The output file.  If unspecified, uses STDOUT",
		 undef],
		['replaceColumn|c=i',
		 '1',
		 'The column in the source file to replace with values from the map.',
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

     
{

	my %opts = parseOpts();

	my %map=();
	my $kcol = $opts{'keyColumn'}-1;
	my $vcol =$opts{'valueColumn'}-1;
	open(my $mapf,$opts{'map'});
	while(<$mapf>){
		chomp;
		my @l = split;
		print STDERR "WARNING: overriding value $map{$l[$kcol]} for key $l[$kcol] with new value $l[$vcol]\n" if(defined($map{$l[$kcol]}));
		$map{$l[$kcol]}=$l[$vcol];
	}
	close($mapf);
	
	my ($inf, $outf);
	if($opts{'source'}){
		open($inf,$opts{'source'});
	}
	else{
		$inf=*STDIN;
	}
	if ($opts{'output'}){
		open($outf,">$opts{'output'}");
	}
	else{
		$outf=*STDOUT;
	}
	my $rcol =$opts{'replaceColumn'}-1;
	while(<$inf>){
		chomp;
		my @l = split;
		if(defined($l[$rcol]) && defined($map{$l[$rcol]})){
			$l[$rcol]=$map{$l[$rcol]};
		}
		else{
			if(!defined($l[$rcol])){
				print STDERR "WARNING: value missing in the replace column for input line. Not performing the replacement in line\n$.: $_\n";
			}
			else{
				print STDERR "WARNING: key $l[$rcol] has no value in map. Not performing the replacement in line\n$.: $_\n";
			}
		}
		print $outf join(' ', @l ,"\n");
	}
	close($inf);
	close($outf);

}
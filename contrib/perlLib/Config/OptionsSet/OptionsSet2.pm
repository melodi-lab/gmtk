#!/usr/bin/env perl
use warnings;
use strict;
package Config::OptionsSet::OptionsSet2;

## @class Config::OptionsSet::OptionsSet2
#utilities dealing with IO, display and manipulation of heirachies of sets of options
#based on Config::General 
#
#This is a replacement for Config::OptionsSet::OptionsSet 

use Carp;
use Cwd;
use Getopt::Long;
use Config::General qw(ParseConfig SaveConfig SaveConfigString);
use Data::Dumper;
use List::Util qw[min max];
use Storable qw(dclone);

use Exporter;
our @ISA=qw(Exporter);
our @EXPORT_OK=qw(  readOpts writeOpts prettyPrintOpts mergeOpts expandOpts
					cppArgsToCppBlock cppBlockToCppArgs parseCppCommArgs genCppCommArgs);


## @fn
#read parameter value pairs from optFile into optHash
#an optional third argument contextHash
#is a hash of variable names to variable values
#which will be used to interpolate the option values.
#options defined in containing blocks can be used as variables
#in inner blocks, same as in config::general
#
#eg.
# if $contexthash->{'COMPUTE_SITE'}=='UIUC' and
# and optFiles contains a line
#of1: some/path/$COMPUTE_SITE/somefile.txt
# optHash will contain
# $optHash->{'of1'} == some/path/UIUC/somefile.txt
sub readOpts{
	my $optFile = shift;
	my $optHash = shift;
	my $contextHash = shift;
	my $fieldSep = shift;
	

	if (!defined($optFile)){
		confess "readOpts() requires filename.";
	}
	if (! -f $optFile ){
		confess "file $optFile does not exist.";
	}

	local %ENV=();
	if ($contextHash){
		%ENV = %$contextHash;
	}
	#print Dumper(%ENV); 
	my @settings=(	'-ConfigFile' => $optFile,
					#'-ConfigHash' => $optHash,
					'-InterPolateVars' => 1,
					'-InterPolateEnv' => 1,
					'-StrictVars' => 1,
					);
	if($fieldSep){
		@settings=(	@settings,
					'-SplitPolicy' => 'custom',
					'-SplitDelimiter' => '\s*'.$fieldSep.'\s*',
		);
	}

 	my %config = ParseConfig(@settings);

 	%$optHash=(%$optHash,%config);
}

## @fn
#write out the parameters in $opt to file $optFile
sub writeOpts{
	my $opt = shift;
	my $optFile = shift;
	my $fieldSep = shift;
	
	my $conf;
	if($fieldSep){
		$conf = new Config::General('-ConfigHash' => $opt,
								'-SplitPolicy' => 'custom',
								'-SplitDelimiter' => '\s*'.$fieldSep.'\s*',
								'-StoreDelimiter' => $fieldSep);
	}
	else{
		$conf = new Config::General('-ConfigHash' => $opt);
	}
	$conf->save_file($optFile);
}

sub prettyPrintOpts{
	my $opt = shift;
	my $indent = shift;
	if (!$indent){
		$indent=0;
	}
	my $arg;
	my $val;
	#print "in prettyPrintOpts, indent=$indent\n";
	#print Dumper($opt);
	my $maxOptLength=max(map {length($_)} keys %$opt);
	if (!$maxOptLength){
		$maxOptLength=0;
	}
	my ($rows,$cols) = split(/ /,`stty size`);		
	my $format  = "format STDOUT = \n"
		. ' ' x $indent. '^' . '<' x $maxOptLength . '  ^'.'<'x($cols-$maxOptLength-$indent-5)." ~~\n"
		. '$arg, $val' . "\n"
		. ".\n";
	#print $format;
	my @args=sort keys %$opt;
	#print non-recursive things first
	if (keys %$opt ==0){
		print ' 'x$indent."EMPTY\n";
	}
	
	foreach $arg (grep {ref($opt->{$_}) ne 'HASH'} @args){
		if (ref($opt->{$arg}) eq 'ARRAY'){
			my @valList = map {defined($_) ? $_: '<UNDEF>'} @{$opt->{$arg}}; 
			$val = "[ @valList ]";
		}
		else{
			if (defined($opt->{$arg})){
				$val = "$opt->{$arg}";			
			}
			else{
				$val = "<UNDEF>";
			}
		}
		#truncate values that are too long
		if (length($val)>800){
			$val=substr($val,0,400)." <TRUNCATED>";
		}
		
		{
			local $SIG{__WARN__} = sub { };
			eval $format;
			die $@ if $@;
		}
		write ();
	}
	#now print recursive things
	foreach $arg (grep {ref($opt->{$_}) eq 'HASH'} @args){
		my $hashKey=$arg;
		$val = '==>';
		{
			local $SIG{__WARN__} = sub { };
			eval $format;
			die $@ if $@;
		}
		write ();
		#print "arg: $hashKey\n";
		#print Dumper($opt->{$hashKey});
		prettyPrintOpts($opt->{$hashKey},$indent+3);
	}
	
	#print "\n";
}


sub expandOpts{
	my $opt = shift;
	my $name = shift;
	my $soFar = {};
	return _expandOpts($opt,$soFar,$name);
}

sub _expandOpts{
	my ($opt, $soFar, $name) = @_;
	if (!defined($opt->{'copyContaining'})){
		$soFar={};
	}
	mergeOpts($soFar,$opt);
	#print"################ $name\n";
	#print Dumper($soFar);
	if ($name){
		my ($curName,$rest)=split(/\./,$name,2);
		if (!defined($opt->{$curName})){
			die("The suffix $name of the specified tree branch does not exist.");
		}
		#print "about to call ExpandOpts(curName,rest)=($curName,$rest)\n";
		return _expandOpts($opt->{$curName},$soFar,$rest);
	}
	else{
		return $soFar;
	}
	
}

##recursively copy options, recursing into sub-hashes.
#if the source hash has <copyContaining/> child it will not be copied 
#sub-hashes are merged: if a value exists already, it is overwritten
sub mergeOpts{
	my $targetOpt = shift;
	my $sourceOpt = shift;
	foreach my $k (keys %$sourceOpt){
		if(ref($sourceOpt->{$k}) eq 'HASH'){
			if ( ($k ne 'copyContaining') && !defined($sourceOpt->{$k}->{'copyContaining'})){
				if (!defined($targetOpt->{$k})){
					$targetOpt->{$k}={};
				}
				#print"calling MERGEOPTS with $k\n";
				mergeOpts($targetOpt->{$k},$sourceOpt->{$k});
				#prettyPrint($targetOpt);
			}
		}
		else{
			if ($sourceOpt->{$k} eq '__UNDEF'){
				if(defined($targetOpt->{$k})){
					delete $targetOpt->{$k};
				}
			}
			else{
				$targetOpt->{$k}=$sourceOpt->{$k};			
			}
		}
	}
}

## @fn
#parses the $cppCommString and returns references to a list and a hash
#the string may be either '-quoted or unquoted
#the list is a list of includes specified via the -I flag
#the hash is a map from defines to their value speciied via the -DVAR=value flag
sub parseCppCommArgs{
	my $cppCommString=shift;
	#print "$cppCommString\n";
	my (@includes, @defines);
	my %defines;

	if (defined($cppCommString)){
	
		#GetoptionsFromString from string isn't until Getopt-Long-2.36 and not 
		#yet in perl 5.8.5 distribution
		#$ret = GetoptionsFromString(
		#	$cppCommString, 
		#	"I=s" => \@includes,
		#	"D=s" => \@defines,
		#);
		my @tempARGV= @ARGV;
		$cppCommString =~ s/^\s*\'//;#get rid of the quotes
		$cppCommString =~ s/\'\s*$//;
		$cppCommString =~ s/\s-D/ -D /g;# put a space between -D and -I and param (cpp allows no space)
		$cppCommString =~ s/\s-I/ -I /g;
		@ARGV = split(' ',$cppCommString);
		GetOptions(
			"I=s" => \@includes,
			"D=s" => \@defines,
		);
		@ARGV = @tempARGV;
		foreach (@defines){
			(my $var , my $val) = split /=/;
			$defines{$var}=$val;
		}
	}
	return (\@includes, \%defines);
}

## @fn
#the inverse of parseCppCommArgs
#takes references to a list and a hash and returns a '-quoted string suitable for cppComm
#the list is a list of includes specified via the -I flag
#the hash is a map from defines to their value speciied via the -DVAR=value flag
sub genCppCommArgs{
	my ($includes, $defines)=@_;
	if (!ref($includes)){
		$includes=[$includes];
	}
	my $cppCommString=join(' ',
			       map("-I $_",@$includes), 
			       map("-D $_".(defined($defines->{$_})?"=$defines->{$_}":''),(keys %$defines)));
	#print "$cppCommString\nincludes, defines follow\n".Dumper($includes, $defines );

	return "'$cppCommString'";
}

## @fn
#The argument list is a source configuration optionHash and an optional destination key name.
#If omitted, destination key defaults to cppComm.
#If the option Hash contains the <cpp> ... </cpp> block,
#move the variables in the block into the cppComm line in format suitable for the cpp preprocessor.
#The <cpp> block has the form
#<cpp>
#	<defines>
#		arg1 value1
#		arg2 value2
#	</defines>
#	include include/Path1
#	include include/Path2
#</cpp>
#(using the -D and the -I flags for defines and includes)
#The order of the defines and includes may  not be the same as the order in the @cppComm@.
#However, if the cppComm already has an include or define that is being added, the one in the @cppComm@ instruction 
#overwrites the definition already present in cppComm.
#@return a copy of the optionHash with the <cpp> block converted.
#constant function. (no side effects, the arguments are not modified.)
sub cppBlockToCppArgs{
	my $source = shift;
	my $destKey = shift;
	if (!$destKey){
		$destKey = 'cppComm';
	}
	
	my $dest=dclone($source);
	#my %args = %$conf;
	if ($dest->{'cpp'} && ref($dest->{'cpp'}) eq 'HASH'){
		if (!defined($dest->{$destKey})){
			$dest->{$destKey}='';
		}
		my ($includes, $defines) = parseCppCommArgs($dest->{$destKey});
		#print "parsed\n".Dumper([$includes, $defines]);
		my %includesHash;	#to remove duplicates
		@includesHash{@$includes}=@$includes;
		my %cpp=%{$dest->{'cpp'}};
		delete $dest->{'cpp'};
		foreach (keys %cpp){
			if ($_ eq 'defines' && ref($cpp{$_}) eq 'HASH'){
				#print Dumper($defines);
				foreach (keys %{$cpp{'defines'}}){
					$defines->{$_}=$cpp{'defines'}->{$_};					
				}
				
			}
			elsif($_ eq 'include' && !ref($cpp{$_})){#an include directive
				#print "adding I $_=$args{$_}\n";
				@includesHash{$cpp{'include'}}=$cpp{'include'};
			}
			elsif($_ eq 'include' && ref($cpp{$_}) eq 'ARRAY'){#an include directive
				#print "adding I $_=$args{$_}\n";
				@includesHash{@{$cpp{'include'}}}=@{$cpp{'include'}};
			}
			else{#unrecognized directive
				die "illegal format of cpp block: key $_ not recognized. ref() is :".ref($cpp{$_});
			}			
		}
		#print "after add\n".Dumper([$includes, $defines]);
		$dest->{$destKey} = genCppCommArgs([keys %includesHash], $defines);
	}
	#print Dumper(\@resultConfigList)."\n";
	return $dest;
}

#@fn convert to cpp command args to a <cpp></cpp> block
#the inverse of cppBlockTocppArgs
sub cppArgsToCppBlock{
	my $source = shift;
	my $cppArgsKey = shift;
	if (!$cppArgsKey){
		$cppArgsKey = 'cppComm';
	}
	
	my $dest=dclone($source);
	
	my ($includes, $defines) = parseCppCommArgs($dest->{$cppArgsKey});
	delete $dest->{$cppArgsKey};
	my %cpp=('include'=>$includes,'defines'=>$defines);
	$dest->{'cpp'}=\%cpp;
	#print("##################\n");
	#prettyPrintOpts($dest);	
		
	return $dest;
}

1;
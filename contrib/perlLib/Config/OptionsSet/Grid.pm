#!/usr/bin/env perl
package Config::OptionsSet::Grid;
use warnings;
use strict;

use File::Path;
use Carp qw(confess);
use Getopt::Long;
use Config::OptionsSet::OptionsSet qw(readOpts writeOpts prettyPrintOpts overwriteOpts 
					setDiff hashDiff hashSymDiff);
use Data::Dumper;
use Exporter;
our @ISA=qw(Exporter);
our @EXPORT_OK=qw(	gridOpts convertCppCommArgs parseCppCommArgs genCppCommArgs);

## @class Config::OptionsSet::Grid
#this module helps with generating config files over a grid of parameters.
#i.e., it takes a template config file, a n-dimentional array (grid) of 
#parameters and generates a config file instance for each point in grid

## @fn
#The options template and the grid parameters are all stored in the hash opt.
#all keys starting with @ specify grid parameters and other instructions 
#all other keys are part of the template, and will be present unaltered in 
#all config instances.
#the format for grid parameters is:
#
#	@N@<optName> <value1>,<value2>,...
#
#where 
#
#N is the grid number.
#
#<optName> is the option name
#
#<value1> and <value2> are the values which will be assigned to the optName in various instances of the config files. They are seperated by commas
#
#For example, the following config template 
# @verbatim
#  arg1 1 2
#  cppComm '-I someDir -D someArg'
#  @0@arg2 1,3 4
#  @0@arg3 4,5
#  @1@arg3 someVal
#  @cppComm@ arg3=D
#
# Will generate the following 5 config instances
#  arg1 1 2
#  arg2 1 
#  cppComm '-I someDir -D someArg -D arg3=4'
#
#  arg1 1 2
#  arg2 3 4
#  cppComm '-I someDir -D someArg -D arg3=4'
#
#  arg1 1 2
#  arg2 1
#  cppComm '-I someDir -D someArg -D arg3=5'
#
#  arg1 1 2
#  arg2 3 4
#  cppComm '-I someDir -D someArg -D arg3=5'
#
#  arg1 1 2
#  cppComm '-I someDir -D someArg -D arg3=someVal'
#@endverbatim
# The first arg with multiple options will be iterated in the outer-most loop 
# (arg2 in the example above)
#
# The following additional instructions are allowed:
# \@shortNames@ instruction to generate a configuration name.
# It is your responsibility to include enough arguments to make the name unique, 
# otherwise you run the risk of having different configs having the same name
# e.g. \@shortNames@ arg2=a2 arg3=a3 
# generates a name a2_1.a3_4 for the first config, a2_3 4.a3_4 for second config and so on.
#
#The order of the args is the order in which the names will be constructed.
#If shortnames are unique, then identical instance names implies identical config instances.
#
# \@cppComm@ instruction to move the args named in the cppCommInstruction to the cppComm
#string.  It is not applied by gridOpts autmatically.  Instead the \@cppComm@ instruction will 
#be copied into every instance.  To actually move the parameters into the cppComm, call  
#convertCppCommArgs() on the instance.
#The params can be treated either as defines (=D) or includes (=I).
#In addition to =D and =I, there is also the erase (=E) flag, meaning that 
#the parameter will simply be erased when cppComm is constructed via convertCppCommArgs function
#This is useful if you want to use the optionHash for your own purposes before sending it to some 
#GMTK program which will choke if it sees an unfamiliar option.
#In the example above arg3 is treated as define.  The order in cppComm is the order
#in which the args will be appended to cppComm
#
#returns a pair of list references
#the first is a list of config instances
#the second is a list of config names generated from config assignments 
#
#If an arg has a special value __UNDEF, the option is removed 
#
sub gridOpts{
	my $opt = shift;
	
	#@cppComm@ is treated as a regular template arg
	my @templateArgs = grep(!/^@/ || /^\@cppComm@/,keys %$opt);
	my @gridArgs = grep(/^@/ && !/^\@cppComm@/,keys %$opt);
	

	my %template;
	@template{@templateArgs} = @$opt{@templateArgs};


	my %shortNames;
	my @grids=();

	foreach (@gridArgs){
		if(/^@([0-9]+)@(.*)/){	#a grid arg
			$shortNames{$2}=undef;	#initially no short names are defined
			$grids[$1]->{$2}=$opt->{$_};
		}
	}
	@grids = grep {$_} @grids;

	my @nameOrder=();
	if($opt->{'@shortNames@'}){ #shortnames
		my @names = split(' ',$opt->{'@shortNames@'});
		foreach(@names){
			my ($name, $shortName) = split(/=/,$_);
			$shortNames{$name} = $shortName;
			push (@nameOrder,$name); 
		}
	}

	#in case some names are not present in the shortNames
	my @allNames = keys %shortNames;
	push(@nameOrder, setDiff(\@allNames, \@nameOrder));
	#print "@nameOrder\n";
	#prettyPrintOpts(\%shortNames);

	my @instances;
	my @instanceNames;

	foreach my $grid (@grids){
		#print "creating a grid of configs over the following args:\n";
		#prettyPrintOpts($grid);

		my @curArgOrder = grep {defined($grid->{$_})} @nameOrder;
		my($curInst,$curInstName) = instantiate($grid,\%template,'',\%shortNames, \@curArgOrder);


		push (@instances, @$curInst);
		push (@instanceNames, @$curInstName);
	}

	
	#if($opt->{'@cppComm@'}){ #some grid args are really options to CPP
	#	convertCppCommArgs($opt->{'@cppComm@'},\@instances);
	#}

	if(@instanceNames==0){
		Carp::carp ("WARNING: no config instances generated");
	}

	return (\@instances,\@instanceNames);

}

## @fn
#parses the CppCommArgs and returns references to a list and a hash
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
	my $cppCommString=join(' ',
			       map("-I $_",@$includes), 
			       map("-D $_".(defined($defines->{$_})?"=$defines->{$_}":''),(keys %$defines)));
	#print "$cppCommString\nincludes, defines follow\n".Dumper($includes, $defines );

	return "'$cppCommString'";
}

## @fn
#The argument list is a list of configuration optionHashes.
#For each configuration in the argument list which contains the @cppComm@ instruction,
#move the variables in the instruction into the cppComm line in format suitable for the cpp preprocessor.
#(using the -D and the -I flags for defines and includes)
#The order of the defines and includes may  not be the same as the order in the @cppComm@.
#However, if the cppString already has an include or define that is being added, the one in the @cppComm@ instruction 
#overwrites the definition already present in cppComm.
#returns the list of optionHashes with the @cppComm instruction applied.
#constant function. (no side effects, the arguments are not modified.)
sub convertCppCommArgs{
	my @configList = @_;
	my @resultConfigList=();
	foreach (@configList){
		my %args = %$_;
		if ($args{'@cppComm@'}){
			my $cppCommInstruction='';
			$cppCommInstruction = $args{'@cppComm@'};
			delete $args{'@cppComm@'};	
			my @cppCommPairs = split(' ',$cppCommInstruction);
			my %cppCommArgs=();
			my @argOrder=(); #FIXME this argOrder seems unnecessary?
			foreach(@cppCommPairs){
				my ($cppCommArg, $type) = split(/=/,$_);
				$cppCommArgs{$cppCommArg} = $type;
				push (@argOrder,$cppCommArg); 
			}
	
			my @curCppOpts;
			my ($includes, $defines) = parseCppCommArgs($args{'cppComm'});
			#print "parsed\n".Dumper([$includes, $defines]);
			my %includesHash;	#to remove duplicates
			@includesHash{@$includes}=@$includes;
			foreach (@argOrder){
				if(defined($args{$_})){ #the cppComm Arg actually is present in current config instance
					if($cppCommArgs{$_} eq 'D'){#a define directive
						#print "adding D $_=$args{$_}\n";
						$defines->{$_}=$args{$_};
					}
					elsif($cppCommArgs{$_} eq 'I'){#an include directive
						#print "adding I $_=$args{$_}\n";
						$includesHash{$_}=$args{$_};
					}
					elsif($cppCommArgs{$_} eq 'E'){#an erase directive (not part of cpp)
						1;
					}
					else{#unrecognized directive
						die "don't know how to handle cpp option $cppCommArgs{$_} for arg $_\n";
					}
				}
			}
			#print "after add\n".Dumper([$includes, $defines]);
			$args{'cppComm'} = genCppCommArgs([keys %includesHash], $defines);
			delete @args{@argOrder};
		}
		push @resultConfigList, \%args;
	}
	#print Dumper(\@resultConfigList)."\n";
	return @resultConfigList;
}


## @fn
#  @return (\@instances,\@instanceNames) given a grid, a template and a namePrefix.
sub instantiate{
	my ($grid, $template,$namePrefix,$shortNames,$argOrder) = @_;
	#print "instantiate($grid, $template,$namePrefix,%{$shortNames},@$argOrder)\n";

	my @subArgOrder = @$argOrder;
	if (!@subArgOrder){	#the grid is empty, so return just the template
		return ([$template], [$namePrefix]);
	}

	my $arg = shift @subArgOrder;
	my %subGrid = %$grid;
	my @values = split(/\s*,\s*/,$subGrid{$arg});
	delete $subGrid{$arg};

	my @instances;
	my @instanceNames;

	if(scalar(@values)==0){
		print Dumper(@_);
    	confess("parameter $arg has no value assigned to it.  Use the special __UNDEF value if you really want to undefine $arg in this config instance.");
	}
	foreach(@values){
		my $subNamePrefix=$namePrefix;
		if(defined($shortNames->{$arg})){
			if ($subNamePrefix){
				$subNamePrefix = "$subNamePrefix.";
			}
			#print "$shortNames->{$arg}\n";
			$subNamePrefix="$subNamePrefix$shortNames->{$arg}_$_"; 
		}
		my %subTemplate=%$template;

		if ($_ eq '__UNDEF'){
			delete $subTemplate{$arg} 
		}
		else{
			$subTemplate{$arg}=$_;
		}

		my($curInst,$curInstName) = instantiate(\%subGrid, \%subTemplate,$subNamePrefix,$shortNames,\@subArgOrder);
		push (@instances, @$curInst);
		push (@instanceNames, @$curInstName);
	}
	return (\@instances,\@instanceNames);
}


1;

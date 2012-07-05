#!/usr/bin/env perl

## @class
#
#Getopt::Lazier A command-line parser that also generates --help documentation and
#and performs validation on options.
#
#
#@verbatim
#	use Getopt::Lazier;
#	my $usage = <<'EOS';
#
#	Usage: 
#	$prog [--help] [options]  
#
#	An example documentation for a program.  Write and format things here as you like, in any language.
#	Two special variables are recognized here: 
#	\$prog is the base name of the script
#	\$opts is the return value of docOptions(...)
#
#	Options:
#
#	$opts
#
#	Return values:
#		0  success
#		-1 error
#
#	EOS
#
#	my $optSpec = [
#		['help|h', 
#			undef,
#			'print this help message and exit'],
#		['minimalisticOption|m=s'],
#		['optionalStringArg|o=s',
#			undef, 
#			'This string option is documented here.'],
#		['optionalStringArg2|aliasForOptionalStringArg2|o=s',
#			'Hello', 
#			'This string pption is has a default.'],
#		['requiredStringArg|r=s', 
#			'hideho', 
#			'This is a required string argument as indicated by the "r" validator. If the parameter is not specified, validateOptions() will return an error.',
#			'required'],
#		['someRequiredFile=s', 
#			undef,
#			'This argument is required to be present, and must be an existing filename as indeicated by the "r" and "f" validators.  If either of these validators fail, validateOptions() will return an error.',
#			'required fileExists'],
#		['someOptionalDir|d=s', 
#			undef,
#			'This argument is not required to be present, but if it is present, it must be the name of an existing dir.',
#			'dirExists'],
#	];
#
#
#	my %opts;
#	exit(-1) if (!parseOptions(\%opts, $optSpec));
#
#	if($opts{'help'}){
#		standardUsage($optSpec, $usage);
#		exit(0);
#	}	
#
#	my @err=validateOptions(\%opts, $optSpec);
#	if (@err){
#		print "@err";
#		standardUsage($optSpec, $usage);
#		exit(-1);
#	}
#@endverbatim
#
#Generates the following output if called with no arguments:
#
#@verbatim
#    --someRequiredFile: required but not specified.
#
#    Usage: 
#    test.pl [--help] [options]  
#
#    An example documentation for a program.  Write and format things here as you like, in any language.
#    Two special variables are recognized here: 
#    $prog is the base name of the script
#    $opts is the return value of docOptions(...)
#
#    Options:
#
#    --help                        print this help message and exit
#    -h
#
#    --minimalisticOption=STR
#    -m
#
#    --optionalStringArg=STR       This string option is documented here.
#    -o
#
#    --optionalStringArg2=STR      This string pption is has a default.
#    --aliasForOptionalStringArg2  (Default: Hello)
#    -o
#
#    --requiredStringArg=STR       This is a required string argument as indicated by the "r" validator. If
#    -r                            the parameter is not specified, validateOptions() will return an error.
#                                  REQUIRED.
#                                  (Default: hideho)
#
#    --someRequiredFile=STR        This argument is required to be present, and must be an existing filename
#                                  as indeicated by the "r" and "f" validators.  If either of these
#                                  validators fail, validateOptions() will return an error.
#                                  REQUIRED.
#                                  File must exist if specified.
#
#    --someOptionalDir=STR         This argument is not required to be present, but if it is present, it
#    -d                            must be the name of an existing dir.
#                                  Directory must exist if specified.
#
#
#
#    Return values:
#        0  success
#        -1 error
#@endverbatim
package Getopt::Lazier;

use warnings;
use strict;
use List::Util qw[min max];
use Getopt::Long;
use File::Basename;

use Carp;

use Exporter;
our @ISA=qw(Exporter);
our @EXPORT=qw(	parseOptions docOptions validateOptions standardUsage);

## @fn
#optSpec is of form [ [optSpec1, defaultForOpt1, docForOpt1], [optSpec2, defaultForOpt2, docForOpt2], ...]
#
#
sub parseOptions{
	my ($opts, $optSpec) = @_;	

	my @specList;
	foreach my $o (@$optSpec){
		if(defined($o->[1])){ #default is present
			$o->[0] =~ /^([^|!+=:]+)/;
			my $optName=$1;
			$opts->{$optName}=$o->[1];
		}
		push @specList, $o->[0];
	}
	return GetOptions($opts,@specList);
}

			
sub docOptions{
	my ($optSpec) = @_;	
			

	my ($prettySpec,$maxOptLength)=prettySpec(map($_->[0],@$optSpec));
			
	my ($rows,$cols) = split(/ /,`stty size 2>/dev/null`);	#FIXME probably there is a more efficient way of getting the size other than stty
	if(!$cols){
	    $cols=80;
	}
	elsif($cols<40){
		$cols=40;
	}
	my $format  = '^' . '>' x $maxOptLength . '  ^'.'<'x($cols-$maxOptLength-5)." ~~\n\n";

	my $ret='';
	for (my $i=0; $i<@$optSpec; $i++){
		my ($oSpec, $oDefault, $oDesc, $validatorSpec) = @{$optSpec->[$i]};
		$oDesc='' if (!$oDesc);

		if(defined($validatorSpec)){ 
			my @funcValidatorList=toValidatorFunctionList($validatorSpec);
			map { $oDesc .= "\r".&$_()} @funcValidatorList;
		}
		
		if (defined($oDefault)){
			$oDesc .= "\r(Default: $oDefault)";
		}
		$oSpec=$prettySpec->[$i];
		$ret .= swrite($format, ($oSpec, $oDesc));
		
	}
	

	return $ret;
}

## @fn
#
#How validator functions should behave.
#
#A validator can be called with either a scalar or no arguments.  
#If a validator is called with a scalar, the scalar will be the argument to be validated.
#e.g. It can be an integer, or a string, or an undef, or a ref to a list, or anything else that getOpts::Long GetOptions() may return.
#If a validator is called with a scalar, it should return either an empty list if the argument is valid
#or it should return the error message describing why the the argument is not valid.
#
#If a validator is called with no arguments, the validator must return a short message describing the conditions for validation success.
#this message will be appended to the help documentation for that argument.
sub validateOptions{
	my ($opts, $optSpec) = @_;	

	my @errors=();
	foreach my $o (@$optSpec){
		if(defined($o->[3])){ #a validator list or string is present
			my @funcValidatorList=toValidatorFunctionList($o->[3]);
			#now @funcValidatorList is the complete list of validator functions.  Apply them.
			foreach my $func (@funcValidatorList){
				my @aliases =getAliases($o->[0]);
				my $argName = shift @aliases;
				#find the value in the opts hash.
				my $val=$opts->{$argName};
				my $err = &$func($val);	
				if ($err){
					push @errors, "--$argName: $err\n";
					last;
				}
			}
		}

	}
	return @errors;
}

## @fn
#prints the usage information
sub standardUsage{
	my ($optSpec, $usage) = @_;	
	my $prog =basename($0);
	my $opts =docOptions($optSpec);

	print eval ('"'.$usage.'"');	

}


## @fn
#configure parameters get passed on to     Getopt::Long::Configure
sub Configure{
    Getopt::Long::Configure(@_);
}

##########################################################################################################
# Pre-defined validators:
##########################################################################################################

my %commonValidators = (
 'required' => \&requiredValidator,
 'fileExists' => \&fileExistsValidator,
 'dirExists' => \&dirExistsValidator
);

## @fn
# return an error string if arg is defined no dir with such name exists
# If arg is undefined, or if the file exists, return an empty list
sub requiredValidator{
	my @value = @_;	
	return "REQUIRED." if(@value==0);

	if (!defined($value[0])){
		return "required but not specified.";
	}
	return ();
}

## @fn
# return an error string if arg is defined no dir with such name exists
# If arg is undefined, or if the file exists, return an empty list
sub dirExistsValidator{
	my @value = @_;	
	return "Directory must exist if specified." if(@value==0);
	if (defined ($value[0]) && ! -d $value[0]){
		return "dir $value[0] does not exist.";
	}
	return ();
}

## @fn
# return an error string if arg is defined no file with such name exists
# If arg is undefined, or if the file exists, return an empty list
sub fileExistsValidator{
	my @value = @_;	
	return "File must exist if specified." if(@value==0);
	
	if (defined ($value[0]) && ! -f $value[0] ){
		return "file $value[0] does not exist.";
	}
}


##########################################################################################################
# Utility functions, of no use outside this module.  These aren't the functions you are looking for.
##########################################################################################################


## @fn
#swrite() subroutine, which is to write() what sprintf() is to printf().
#taken from perlform doc
sub swrite {
	croak "usage: swrite PICTURE ARGS" unless @_;
	#print "*** @_\n";
	my $format = shift;
	local $^A = "";
	formline($format,@_);
	return $^A;
}

## @fn
#returns a list of argument aliases from an optionSpec.  The order of aliases is preserved
sub getAliases{
	my $oSpec = shift;
	my @aliases = split('\|',$oSpec);
	$aliases[$#aliases] =~ /([^!+=:]*)/;
	$aliases[$#aliases] = $1;	
	return @aliases;
}

## @fn
#returns the option kind (the stuff following option aliases from an optionSpec.
sub getOptionKind{
	my $oSpec = shift;
	$oSpec =~ /[^!+=:]*(.*)/;
	return $1;
}

## @fn
#Takes a string or a ref to a list of strings or function refs denoting validators, 
#and returns only function refs denoting validators, in the same
#order as the input
sub toValidatorFunctionList{
	my $validatorSpec = $_[0];

	my @validatorList=();
	if (!ref($validatorSpec)){ #the validator list is a string of commonValidators keys. put it in a list.
		@validatorList=($validatorSpec);
	}
	elsif(ref($validatorSpec) eq 'ARRAY'){ #it's a list 
		@validatorList=@{$validatorSpec};
	}
	else{
		croak "The validator must be either a string a or a list reference\n";
	}

	my @funcValidatorList=();

	foreach my $v (@validatorList){
		if (ref($v)){ #the validator is a function
			push @funcValidatorList, $v;
		}
		else{#the validator is a string of commonValidators keys. expand them.
			foreach (split(' ',$v)){
				my $func = $commonValidators{$_};
				if ($func){
					push @funcValidatorList, $func;
				}
				else{
					croak "Validation specification error: There is no validator for shortcut '$_'\n";
				}
			}
		}
	}

	return @funcValidatorList;

}


sub prettySpec{
	my @optSpec = @_;	
	my @prettyOptSpec;
	my $maxPrettySpecLen=0;
	foreach my $oSpec (@optSpec){
		my @aliases = getAliases($oSpec);
		#print "@aliases\n";
		for (my $i=0; $i<@aliases; $i++){
			if (length($aliases[$i])>1){
				$aliases[$i] = '--' . $aliases[$i];
			}
			else{
				$aliases[$i] = '-' . $aliases[$i];
			}
			
		}
		my $spec=getOptionKind($oSpec);
		$spec =~ s/ .*//;
		my %types =('s' => 'STR',
					'i' => 'INT',
					'f' => 'FLOAT',
					'o' => 'EXT_INT'
					);
		if ($spec =~ /^([:=])([sifo])$/){ #the other parameter specs seem too complicated to me. Not printing them nicely discourages their use.
			if ($1 eq ':'){
				$spec = "[=$types{$2}]";
			}
			else{
				$spec = "=$types{$2}";
			}
		}
		$aliases[0] .= $spec;
		$maxPrettySpecLen=max($maxPrettySpecLen, map {length($_)} @aliases);
		$oSpec =join("\r",@aliases);
		push @prettyOptSpec, $oSpec;
	}
	

	return (\@prettyOptSpec,$maxPrettySpecLen);
}



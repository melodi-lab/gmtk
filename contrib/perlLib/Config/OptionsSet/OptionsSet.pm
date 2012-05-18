#!/usr/bin/env perl
package Config::OptionsSet::OptionsSet;
use warnings;
use strict;

## @class Config::OptionsSet::OptionsSet
#utilities dealing with IO, display and manipulation of sets of options

use Carp;
use Cwd;
use String::Interpolate;
use Data::Dumper;
use List::Util qw[min max];
use Exporter;
our @ISA=qw(Exporter);
our @EXPORT_OK=qw(  readOpts writeOpts prettyPrintOpts overwriteOpts 
					getNamespaceOpts getRootNamespaceOpts
					setDiff hashDiff hashSymDiff);


## @fn
#read parameter value pairs from optFile into optHash
#an optional third argument contextHash
#is a hash of variable names to variable values
#which will be used to interpolate the option values
#
#eg.
# if $contexthash->{'COMPUTE_SITE'}=='UIUC' and
# and optFiles contains a line
#of1: some/path/$COMPUTE_SITE/somefile.txt
# optHash will contain
# $optHash->{'of1'} == some/path/UIUC/somefile.txt
#
#interpolation is done using perl's string interpolation, and follows perl's rules
sub readOpts{
	my $optFile = shift;
	my $optHash = shift;
	my $contextHash = shift;

	my $interp = safe String::Interpolate;
	#TODO: would be nice to throw exception if the opts file has an undefined variable. right now it maps undefined vars to ''
   	#$interp->pragma('FATAL','NO_TRAP','NO UNSAFE_SYMBOLS');
	if(defined($contextHash)){
    	$interp->($contextHash);

	}


	open(IN,$optFile) || die  "Couldn't open options file $optFile for reading\n";
	while(<IN>){
		chomp;
		s/^\s+//; #remove leading and trailing whitespace 
		s/\s+$//;
		
		next if (/^#.*/);
		next if (/^$/);
		(my $par, my $val) = split(' ',$_,2);
		$par =~ s/\:$//;
		$val = '' if(!defined($val));
		if(defined($contextHash)){
			my $interpolated_val = $interp->($val);
			$optHash->{$par} = $interpolated_val;
		}
		else{
			$optHash->{$par} = $val
		}
	}

	close(IN)  || confess("Couldn't close options  file $optFile\n");
}

## @fn
#write out the parameters in $opt to file $optFile, one arg-value pair per line, seperated by $fieldSep
sub writeOpts{
	my $opt = shift;
	my $optFile = shift;
	my $fieldSep = shift;
	open(OUT,">$optFile") || die "Couldn't open options file $optFile for writing\n";
	foreach (sort keys %$opt){
		print OUT "$_$fieldSep$opt->{$_}\n";
	}
	close(OUT)  || confess("Couldn't close options file $optFile\n");
	
}

sub prettyPrintOpts{
	my $opt = shift;
	my $arg;
	my $val;

	my $maxOptLength=max(map {length($_)} keys %$opt);
	my ($rows,$cols) = split(/ /,`stty size`);		
	my $format  = "format STDOUT = \n"
		. '^' . '>' x $maxOptLength . '  ^'.'<'x($cols-$maxOptLength-5)." ~~\n"
		. '$arg, $val' . "\n"
		. ".\n";
	#print $format;
	foreach $arg (sort keys %$opt){
		if (ref($opt->{$arg}) eq 'ARRAY'){
			my @valList = map {defined($_) ? $_: 'UNDEF'} @{$opt->{$arg}}; 
			$val = "[ @valList ]";
		}
		else{
			$val = "$opt->{$arg}";
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
	print "\n";
}


## @fn
#returns a subset of $opt, whose keys have have a prefix "$namespace."
#the prefix is stripped of the keys.
#The special option $namespace..import is interpreted to be an instruction
#to copy options from other namespaces into this namespace.
#e.g.
#if %opts={'opt1'=>'v1',
# 		   'ns1.opt1'=>'v2',
# 		   'ns1.nestedNamespace.opt1'=>'v3',
# 		   'ns2.opt1'=>'v4',
# 		   'ns1.import'=>'importedOpt=ns2.opt1 nestedNamespace.opt2=opt1'}
#
#then getNamespaceOpts(\%opts,'ns1') returns
#
# 	{'opt1'=>'v2',
# 	'nestedNamespace.opt1'=>'v3',
# 	'importedOpt'=>'v4',
# 	'nestedNamespace.opt2=v1'}
#
#It is an error to import an option into a namespace where the option is already defined.
#It is also an error to import a non-existing option that is undefined (importing an option which has an undef value is ok).
#
#Namespaces can be nested.  In that case getNamespaceOpts(\%opts,'ns1.nestedNamespace') is equivalent to 
#getNamespaceOpts(getNamespaceOpts(\%opts,'ns1'),'nestedNamespace')
sub getNamespaceOpts{
	my $opt = shift;
	my $namespace = shift;
	if (!$namespace){ #base case
		my %ret = %$opt;
		return \%ret;
	}
	else{
		my ($ns,$rest)=split(/\./,$namespace,2);
		#print "***($ns,$rest)\n";
		my %nsOpts= map {substr($_,0,length($ns)) eq $ns ? (substr($_,length($ns)+1), $opt->{$_}) : () } keys %$opt;
		if ($nsOpts{'.import'}){
			my @im = split(' ', $nsOpts{'.import'});
			for (@im){
				my ($targ, $source) = split('=');
				if (!exists($opt->{$source})){
					croak "error while importing into $ns.$targ: The source $source does not exist.\n";
				}
				if (exists($nsOpts{$targ})){
					croak "error while importing into $ns.$targ: The target exists.\n";
				}
				$nsOpts{$targ}=$opt->{$source};
			}
			delete $nsOpts{'.import'};
		}
		return getNamespaceOpts(\%nsOpts,$rest);
	}
	
}

# @fn 
# @return the options not having a . in the name, stripping away all inner namespaces.  
sub getRootNamespaceOpts{
	my $opt = shift;
	my $namespace = shift;
	my %ret;
	for (keys %$opt){
		if (!/\./){
			$ret{$_}=$opt->{$_}
		}
	}
	return \%ret;
}

## @fn
#overwrites the options in the target hash with the options in the source hash
#only if the target options are already present.
#
#	returns false if some source options were not present in the target hash
#				  In this case the target hash is left unmodified.
#	returns true if target options were succesfully overwritten.
sub overwriteOpts{
	my $targetOpt = shift;
	my $sourceOpt = shift;

	my @sOptArgs =keys %$sourceOpt;
	my @tOptArgs =keys %$targetOpt;

	my @badArgs = setDiff(\@sOptArgs, \@tOptArgs);
	if (@badArgs){
		print "the following distributed config options are not recognized:\n";
		map {print "$_\n"} @badArgs;
		return 0;
	}

	#everything's ok - override the defaults with provided user options
	@$targetOpt{@sOptArgs}= @$sourceOpt{@sOptArgs};
	
	return 1;
}

## @fn
#returns those elements in \@a that are not in \@b
#order is not preserved
sub setDiff{
	my $a = shift;
	my $b = shift;
	my %temp = ();
	@temp{@$a} = ();
	foreach (@$b) {
			delete $temp{$_};
	}
	return keys %temp;
}

## @fn
#returns those keys k of  %a with either of the following conditions:
# defined($a->{k}) and  not defined($b->{k}) 
# $a->{$k} ne $b->{$k}
sub hashDiff{
	my $a = shift;
	my $b = shift;
	my @diffKeys;
	foreach (keys %$a) {
		if ( defined($a->{$_}) && !defined($b->{$_}) || $b->{$_} ne $a->{$_}){
			push(@diffKeys, $_);
			#print "$_\t'$b->{$_}'\t'$a->{$_}'\n";
		}
	}
	return @diffKeys;

}

## @fn
#symmetric hash difference. 
#returns those keys of  %a that are not present in %b and vise versa
#or the values for those keys in %b is different from %a and vise versa 
sub hashSymDiff{
	my $a = shift;
	my $b = shift;

	my %temp;
	@temp{hashDiff($a, $b)}=0;
	@temp{hashDiff($b, $a)}=0;
	return keys %temp;
}

1;
#!/usr/bin/env perl
package AI::GMTK::Util;
use warnings;
use strict;

use Carp;
use Cwd;
use Data::Dumper;
use Exporter;
use Config::OptionsSet::OptionsSet qw(readOpts writeOpts prettyPrintOpts overwriteOpts 
					setDiff hashDiff hashSymDiff);

our @ISA=qw(Exporter);
our @EXPORT_OK=qw(	prepareGmtkCmdArgs checkOpts
					isRestartable notifyOfJobCompletion);

## @fn
#takes in a gmtk args hash, and a filename, writes the args to the filename 
#and returns the args that can be passed to a gmtk program.
#a workaround for broken parsing of args from file in gmtk.
#otherwise only the -argsfile arg would be necessary.
#constant function
#
#concatinated cppComm strings are treated correctly, ie 'some defines''other defines' is read as
#'some defines  other defines'

sub prepareGmtkCmdArgs{

	my %gmtkOpt=%{$_[0]};
	my $argsFname=$_[1];

	
	my $cmdArgs="-argsFile $argsFname";
	#this is a workaround for gmtk not being able to have array variables in the argsFile
	#also -cppComm doesn't work
	#hopefully it will be fixed soon by Chris or Jeff
	foreach(grep (/[0-9]$/, keys %gmtkOpt)){
		$cmdArgs = $cmdArgs . " -$_ $gmtkOpt{$_}";
		delete $gmtkOpt{$_};
	}
	if ($gmtkOpt{'cppComm'}){
		my $cpp = $gmtkOpt{'cppComm'};
		$cpp =~ s/'/ /g;
		$cmdArgs = $cmdArgs . " -cppComm '$cpp'";
		delete $gmtkOpt{'cppComm'};
	}

	writeOpts(\%gmtkOpt, $argsFname, ': ');
	#print "cmdArgs $cmdArgs";
	return $cmdArgs;

}

## @fn
#splits opts into gmtkOpts and distrOpts.
#gmtkOpts are passed untouched to gmtk executables, distrOpts manange behavior of 
#gmtkParallel scripts.
#arg $opt   a ref to a hash table with options
#arg $initDistrOptsFun a ref to a subroutine which returns the complete set of distrOpts
#   with their default values
#arg $fillInOptsFun a ref to a subroutine which verifies the options
sub checkOpts{
    my ($opt, $initDistrOptsFun, $fillInOptsFun) = @_;
    
	my $distrOpt={};
    my $gmtkOpt={};
	&$initDistrOptsFun($distrOpt); #initialize to default values
    my @gmtkArgs = setDiff([keys %$opt], [keys %$distrOpt]);
    #print Dumper($distrOpt, @gmtkArgs);
    @$gmtkOpt{@gmtkArgs}=@$opt{@gmtkArgs};
    my %userDistrOpt=%$opt;
    delete @userDistrOpt{@gmtkArgs};
    
	if(!overwriteOpts($distrOpt, \%userDistrOpt)){
		print "THIS IS A BUG AND SHOULD NOT HAPPEN. Check the distributed params\n";
		return -1;
	}

	#make sure that settings in distributed settings file from the command line are ok
	#and fill in missing settings from defaults if reasonable
    #print Dumper($gmtkOpt, $distrOpt);
    
	my $errCount=0;
	$errCount += &$fillInOptsFun($gmtkOpt, $distrOpt);
	if ($errCount){
		print "$errCount errors found in the settings.\n" ;
		return -1;
	}

    
	if($distrOpt->{'verbosity'} && $distrOpt->{'verbosity'}>0){
		print "using the following gtmk parameters:\n";
		prettyPrintOpts($gmtkOpt);
		
		print "using the following distributed parameters:\n";
		prettyPrintOpts($distrOpt);
	}

    return (0,$gmtkOpt, $distrOpt);
}

## @fn
#@pre  $distrOpt->{'_gmtkArgsBackup'} and $distrOpt->{'_distrArgsBackup'} must contain
#the options filenames from previous run of the job, if the job was run previously 
#
# @return
#returns 0 if the job has not been run before 
#returns 1 if the job has been run before and is restartable
#returns 2 if the job has been run before and is not restartable
sub isRestartable{
	my $gmtkOpt = shift;
	my $distrOpt = shift;
	my $diffAllowedOpts = shift; # a list of options that can differ
	my %gmtkOptOld;
	my %distrOptOld;

	my %allowedDiffs;
	map {$allowedDiffs{$_}=1} @$diffAllowedOpts;
	
	if(-e $distrOpt->{'_gmtkArgsBackup'} && -e $distrOpt->{'_distrArgsBackup'}){
		print "A previously started job detected.\n";
		readOpts($distrOpt->{'_gmtkArgsBackup'}, \%gmtkOptOld);
		readOpts($distrOpt->{'_distrArgsBackup'}, \%distrOptOld);
		
		my @gmtkDiff = hashSymDiff(\%gmtkOptOld, $gmtkOpt); 
		my @distrDiff = hashSymDiff(\%distrOptOld, $distrOpt);
		#some params can be different 
		@distrDiff = grep(!/^_.*/,@distrDiff);
		@distrDiff = grep {!$allowedDiffs{$_}} @distrDiff;

		if ( @gmtkDiff +@distrDiff){
			print "This previously started job is not restartable because the following parameters differ:\n";
			foreach (@gmtkDiff){
				my $old = defined($gmtkOptOld{$_}) ? $gmtkOptOld{$_} : 'UNDEF'; 
				my $new = defined($gmtkOpt->{$_}) ? $gmtkOpt->{$_} : 'UNDEF'; 
				print "gmtk param:\t$_\nold job:    \t'$old'\ncur job:    \t'$new'\n";
			}
			foreach (@distrDiff){
				my $old = defined($distrOptOld{$_}) ? $distrOptOld{$_} : 'UNDEF'; 
				my $new = defined($distrOpt->{$_}) ? $distrOpt->{$_} : 'UNDEF'; 
				print "distr param:\t$_\nold job:   \t'$old'\ncur job:    \t'$new'\n";
			}

			if($distrOpt->{'overrideRestartable'}){
				print "overrideRestartable flag is set. Allowing restart despite the config differences\n";
				return 1;
			}
			else{
				return 2;
			}
		}
		else{
			return 1;
		}
	}
	else{
		print "Starting job from beginning.\n";
		return 0;
	}
}

## @fn
#sends an email if necessary according to options in $distrOpt
#
sub notifyOfJobCompletion{
	my $distrOpt = shift;
	my $exitCode = shift;

	my %exitDesc=(0 => 'success',
				 -1 => 'configuration error. no files written or modified.',
 				 -2 => 'sanity check failed.  No parallel jobs launched.',
				 -3 => 'some parallel tasks failed.',
				 -4 => 'transcription accuracy checking failed.');


	if($distrOpt->{'sendMailOnCompletion'} && $distrOpt->{'emailAddress'}){	
		print "sending email notification of completed job\n";
		my $hn = `hostname`;
		chomp($hn);
		my $cwd = getcwd();
		my $execName = $0;
		my $time = localtime(time);
		my $descr = $exitDesc{$exitCode};
		my $jobName = defined($distrOpt->{'jobName'}) ? $distrOpt->{'jobName'} : "unknown";
my $emailbody = <<EOS;
$execName completed.
job name: $jobName
on machine: $hn
launched from dir: $cwd
wd: $distrOpt->{'wd'}
exit code: $exitCode ($descr).
message generated at: $time 
EOS
		#print "$emailbody\n";
		`echo -n'$emailbody' | mail $distrOpt->{'emailAddress'} -s 'job $jobName finished.'`;
	}

}

1;
#!/usr/bin/env perl
package OS::Util;
use warnings;
use strict;

use Carp;
use Cwd;
use IPC::Run qw( run timeout ) ;
use Exporter;
our @ISA=qw(Exporter);
our @EXPORT_OK=qw(	readLines writeLines tsSystem killChildrenHandler);

## @fn
#reads all the lines from the file and returns them in a list.
#confesses if the file cannot be opened or closed
sub readLines{
	open(IN, "<$_[0]") ||  confess("cannot open <$_[0]");
	my @lines = <IN>;
	close(IN) ||  confess("cannot close <$_[0]");
	return @lines; 
}

sub writeLines{
	open(OUT, ">$_[0]") ||  confess("cannot open >$_[0]");
	foreach (@{$_[1]}){
		print OUT $_;
	}
	close(OUT) ||  confess("cannot close >$_[0]");
	return 0; 
}


## @fn
#the following routine is a a "kill all children and threads and exit" signal handler
#example use (dynamic scope):
#local $SIG{INT} = \&killChildrenHandler; 
#local $SIG{QUIT} = \&killChildrenHandler;
#local $SIG{TERM} = \&killChildrenHandler;
sub killChildrenHandler{
    my $signame = shift;
	print "Caught a SIG$signame.\nKilling my process group with SIG$signame.\n";
	{
		local $SIG{$signame} = 'IGNORE';   # exempt myself
		kill($signame, -$$);               # signal my own process group
	}
	print "Killing self with SIG$signame.\n";
	$SIG{$signame} = 'DEFAULT';
	kill($signame, $$);

}

## @fn
#threadsafe system.  This one actually uses the STDOUT and STDERR that are set for the current thread,
#instead of the first thread in the process.  I think it's a bug in perl.  We do it via IPC::run.
#$? is set as usual
sub tsSystem{
	if(!run (['bash', '-c', $_[0]], \undef, \*STDOUT, \*STDERR)){
        print "Error running command: $_[0]\n";
        confess "error with $?";
    }
}

1;
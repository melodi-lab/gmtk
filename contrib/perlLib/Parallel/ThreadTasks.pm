#!/usr/bin/env perl
package Parallel::ThreadTasks;
use warnings;
use strict;

use threads;
use Thread::Queue; 
use Error ':try';

use Exporter;
our @ISA=qw(Exporter);
our @EXPORT_OK=qw(	spawnAndJoin );

#This module has routines for launching multiple jobs in their own threads 
#redirecting their STDIO and waiting for their completion.


#spawns a thread for every element in tasks, and waits for them all to finish before returning
#each element in tasks should be a list reference containing the subroutine and the arglist
#returns a list of list references, with the ith list reference pointing to the output of ith task.
my $doneQueue = Thread::Queue->new();
sub spawnAndJoin{


	my ($tasks, $taskNames, $logFileNames) = @_;
	my @threadHandles = ();
	foreach(0..$#{$tasks}){
		my @subAndArgs = @{$tasks->[$_]};
		push(@threadHandles, threads->new(\&taskWrapper, $taskNames->[$_], $logFileNames->[$_], @subAndArgs));
	}
	print "waiting for ".scalar(@threadHandles)." threads to finish\n";
	my @ret=();
	foreach(@threadHandles){
		my $tid;
		#spin-lock on the doneQueue, because blocking on join will make perl ignore signals
		#I think this is a perl bug.  Once it's fixed, we don't need the doneQueue at all.
 		while (!($tid = $doneQueue->dequeue_nb())) {
			 sleep(1);
		} 
	}
	foreach(@threadHandles){
		push(@ret, [$_->join()]);
	}

	print "All threads finished\n";
	#print "ignore Scalars leaked warnings here... some sort of a perl threads join bug\n";
	return @ret;
}

sub taskWrapper{
	my ($name, $logFileName, $sub, @args) = @_;
	try{
		print "starting thread ".threads->self->tid()." with STDOUT/STDERR redirected to $logFileName\n";
		my @ret;
		{
			#change stdout/stderr for this thread
			local *main::STDOUT;
			local *main::STDERR;
			open (STDOUT, ">$logFileName") || die "cannot open  >$logFileName for STDOUT";
			*STDERR=*STDOUT;
			#make the pipes flush immediately
  			my $oldfh = select(STDERR);
			$| = 1;
			select(STDOUT);
			$| = 1;
			select($oldfh);

			#run the subroutine
			@ret = &$sub(@args);
			close (STDOUT) || confess(" cannot close file >$logFileName");

		}
		print "thread (id ".threads->self->tid().") completed\n";
		return @ret;
	}
	otherwise{
		my $err = shift;
		print "thread ".threads->self->tid().": caught exception: $@ \n";
	}
	finally{
		$doneQueue->enqueue(threads->tid());
	};
}
1;

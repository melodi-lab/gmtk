#!/usr/bin/env perl
package Parallel::Distribute;

use warnings;
use strict;
use File::Spec;
use File::Basename;
use File::Temp qw( tempfile tempdir );
use Cwd;
use OS::Util qw(tsSystem readLines writeLines);
use Carp;
use List::Util;
use Data::Dumper;
use Exporter;
our @ISA=qw(Exporter);
our @EXPORT_OK=qw(distribute);


## @fn
#now distribute can be called directly without cluttering up the 'ps' output
sub distribute{
	my ($commands_file, $jobName, $sge, $runLocally, $extraArgs) = @_;

	$extraArgs = '' if (!defined($extraArgs));
		
	#where am i?
	my $location = `hostname -d`;
	chomp $location;
	
	#using this expression gives you the user-friendly name of the mounted directory we're in
	#getcwd gives /amd/nfs/... (at least at Edinburgh) which breaks things
	my $cwd = `echo pwd | bash`;
	chomp $cwd;
		
	my @commands = readLines($commands_file);
	chomp(@commands);
	
	if($runLocally){
		tsSystem("bash $commands_file");
	}
	else{
		#use the script name as job name if none is specified
		if (!$jobName){
			($jobName,undef) = split /\s+/, $commands[0];
			$jobName = basename $jobName;
		}
		
		mkdir $sge if (!-e $sge);
		(my $sfh, my $scriptName) = tempfile( "$jobName.XXXX", DIR => "$sge");
		close $sfh;
		            
		if($location eq 'ncsa.uiuc.edu'){
		    bsubLaunchAndWait(\@commands, $scriptName, $sge, $extraArgs, $cwd);
		}
		elsif($location eq 't-illiac.uiuc.edu'){#the trusted illiac cluster
			sbatchLaunchAndWait(\@commands, $scriptName, $sge, $extraArgs, $cwd);
		}
		elsif($location eq "clsp.jhu.edu"){
		
            my $queues = 'all.q@x14.clsp.jhu.edu,all.q@x15.clsp.jhu.edu,all.q@x16.clsp.jhu.edu,all.q@x17.clsp.jhu.edu,all.q@x18.clsp.jhu.edu,all.q@x19.clsp.jhu.edu,all.q@x20.clsp.jhu.edu,all.q@x21.clsp.jhu.edu,all.q@x22.clsp.jhu.edu,all.q@x23.clsp.jhu.edu,all.q@x24.clsp.jhu.edu';
		    $extraArgs .= " -q '$queues' ";
		    qsubLaunchAndWait(\@commands, $scriptName, $sge, $extraArgs, $cwd);
		}
		elsif($location eq "ifp.uiuc.edu"){
		    $extraArgs .= " -cwd ";
		    qsubLaunchAndWait(\@commands, $scriptName, $sge, $extraArgs, $cwd);
		}
		else{ #don't know where we are - assume SGE
		    qsubLaunchAndWait(\@commands, $scriptName, $sge, $extraArgs, $cwd);
		}
		
	}
	return $?;
}

sub bsubLaunchAndWait{
	my ($commandsRef, $scriptName, $sgeWorkDir, $extraArgs, $cwd) = @_;
    
    #the sysadmins at NCSA tell me that jobs with multiple subtasks are 
    #not supported on the tungsten cluster (4/5/07) -Arthur
    #So we don't submit the $scriptName as usual but queue up 
    #all the tasks directly
    my @commands = @$commandsRef;
    my $command_line;
    my @jobIds=();
    
	my $saltedJobName = basename $scriptName;

    foreach my $el (1..@commands) {
        my $c = $commands[$el-1];
        $command_line = "bsub $extraArgs -P kph -J $saltedJobName -e $scriptName.$el.err -o $scriptName.$el.out '$c '";
        #print "$command_line\n";
        my $subStatus =`$command_line`;
        $subStatus =~ /Job <([0-9]+)> is submitted/s;
        push (@jobIds, $1);
    }
    
    print "submitted commands to bsub. Last command follows:\n$command_line\n";
    #print "jobIds @jobIds\n";
    $command_line ="bsub $extraArgs -n 1 -P kph -q debug -J $saltedJobName.Blocker -K -w 'ended(\"$saltedJobName\")' echo 'all jobs finished' ";
    tsSystem($command_line);
    
}

{
	my $saltedJobName;
	sub sbatchHandleDeath {
		print "caught signal $_[0].  Terminating tasks...\n";
		tsSystem("scancel -n $saltedJobName"); return @_; 
	}
	
	sub sbatchLaunchAndWait{
		my ($commandsRef, $scriptName, $sgeWorkDir, $extraArgs, $cwd) = @_;
		my $numJobs = scalar(@$commandsRef);
		my $saltedJobName = basename $scriptName;
		my @jobIds=();
	
		my $bashPath = `which bash`;
		my @scriptLines=("#!$bashPath");
		my $command_line ;
	
		for my $s (qw(HUP  INT  PIPE  TERM)){
			$SIG{$s} = \&sbatchHandleDeath;
		}

		foreach my $el (1..@$commandsRef) {
			$scriptLines[1] = "$commandsRef->[$el-1]\n";
			writeLines("$scriptName.$el", \@scriptLines);
			$command_line = "sbatch $extraArgs -J $saltedJobName -e $scriptName.$el.err -o $scriptName.$el.out $scriptName.$el 2>&1";
			#print "$command_line\n";
			my $subStatus = `$command_line`;
			my $err = ($? >>8);
			croak "sbatch '$command_line' failed with error code $err" if ($err);
			$subStatus =~ /Submitted batch job ([0-9]+)/s;
			push (@jobIds, $1);
		}
		#print "jobIds @jobIds\n";
		print "Launched " .scalar(@jobIds). " tasks for job $saltedJobName.\n";# of the form: $command_line\n";
		#print "Check Status with: squeue -t all -o '%j %T' | grep $saltedJobName\n";
	
		#wait until no more useful work can be done: the state is non of PENDING, RUNNING, SUSPENDED or COMPLETING
		my $stillRunning =1;
		my %taskStatus;
		while ($stillRunning>0){
			%taskStatus = getStatusViaSqueue($saltedJobName);
			#print Dumper(\%taskStatus);
			$stillRunning = List::Util::sum( map {defined($_) ? $_ : 0} @taskStatus{('PENDING', 'RUNNING', 'SUSPENDED', 'COMPLETING')});
			sleep(2);
		}
	
		if (defined($taskStatus{'COMPLETED'}) && $taskStatus{'COMPLETED'}==$numJobs){
			print "SUCCESS: All $numJobs tasks for job $saltedJobName completed successfully.\n";
			return 0;
		}
		else{
			%taskStatus = getStatusViaSqueue($saltedJobName);
			print "ERROR: Some tasks for job $saltedJobName are not in COMPLETED state, and no further work will be done.\n";
			print "NumTasks\tState\n";
			print map("$taskStatus{$_}\t $_\n",keys %taskStatus);
			$?=-2;
			return -2;
		}
	}
}


sub qsubLaunchAndWait{
	my ($commandsRef, $scriptName, $sgeWorkDir, $extraArgs, $cwd) = @_;
    my @commands = @$commandsRef;
	my $numJobs = $#commands + 1;
	my $saltedJobName = basename $scriptName;

    writeScript($scriptName,"$scriptName.\$SGE_TASK_ID", $cwd);
    foreach my $el (0..$#commands) {
		#rewrite the commands to use absolute paths and not relative paths
        #$commands[$el] =~ s/\'/\\\'/g; 
        $commands[$el] =~ s/ \.\// $cwd\//g;
        $commands[$el] =~ s/ \.\.\// $cwd\/\.\.\//g;
        $commands[$el] =~ s/^\.\.\//$cwd\/\.\.\//g;
        #print $scriptName . "." . ($el+1)."\n";
        writeScript($scriptName . "." . ($el+1),$commands[$el],$cwd);
    }

    my $command_line;


    $command_line = "qsub  $extraArgs -t 1-$numJobs -sync y -j y -o $scriptName.\\\$TASK_ID.out -N $saltedJobName $scriptName";
    
    print "submitted: $command_line\n";
    my $oldflush = $|;
    $|=1;
    tsSystem($command_line);
    $|=$oldflush;
}


sub getStatusViaSqueue{
	my $saltedJobName = shift;
	my @lines = `squeue -h -t all -o '%j %T' | grep $saltedJobName | uniq -c`;
	my $err = ($? >>8);
	croak "unable to check job status with squeue -h -t all -o '%j %T' | grep $saltedJobName | uniq -c" if ($err);

	chomp(@lines);
	my %statuses;
	for (@lines){
		my ($count, $job, $status) = split;
		$statuses{$status}=$count;
	}

	return %statuses;
}

sub writeScript{
  my ($scriptName,$command,$cwd) = @_;

  chomp $command;
  open FILE, ">$scriptName" || die "cannot open >$scriptName";
  print FILE "#!/bin/bash\n\n";
  print FILE "cd $cwd\n";
  print FILE "$command\n";
  close FILE;
  chmod 0755, $scriptName;
}

#sub gen_unique_id{
#  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
#  my $datetime = sprintf("%04d-%02d-%02d.%02d-%02d-%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec);
#  my $count = 0;
#  my $filename = "$datetime.$$";
#  while (glob "*.$filename") {
#    $count += 1;
#    $filename = "$datetime.$$.$count";
#  }
#  return $filename;
#}

1;

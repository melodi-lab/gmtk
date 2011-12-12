#!/usr/bin/perl -w

## @file
# @author Partha
#Based on genGMParms.pl that came with the gmtk Aurora tutorial
#

sub usage {
 printf STDERR (
  "Description: Generate the gmtk parameters file for any bunch of mixture collections.  No longer makes a valid I/O trainable file, see the commented out bits to revive that.  .\n" .
  "The m, f, c, p, a, r and v options specify a collections of mixtures, specify them multiple times to get multiple collections. \n".
  "Usage: progname [options]\n" .
   "   -s \#   phone state file\n" .
   "   -m \#       number of mixtures\n" .
   "   -f \#       number of features\n" .
   "   -c \#       number of Gaussian components per mixture\n" .
   "   -p \#       prefix for mixture names e.g. gm.  Always followed by the mixture number e.g. gm5\n" .
   "   -a name     name of the current collection of mixtures\n" .
   "   -r true/false	random means will be used if true, zero otherwise\n" .
   "   -v file     set means and diag variances to those in file. Number of rows in file must\n" .
   "               match number of features.  the row format is <featureIndex> <mean> <var> <other stuff>\n" .
   "               which is what 'obs-print -q -stats ...' outputs\n" .
   "   -o file     file to output the I/O trainable file to\n" .
   "   -n file     file to output the name collections to\n" .
   "   -i          if present, will generate a valid I/O trainable file, just makes the mixture stuff otherwise\n"
);
}

sub checknnint {
    ($opt,$arg) = @_;
    if (( $arg !~ /^[0-9]+$/ ) || ($arg < 0)) {
        print STDERR ("Error: Value \"", $arg, "\" invalid for option ",
                      $opt, " (positive integer expected)\n");
        &usage; exit(-1);
    }
}

sub checktf {
    $arg = $_;
    if ( $arg ne "true" && $arg ne "false"){
	print STDERR "Option -r must be given either true or false, not $arg.  \n";
        &usage; exit(-1);
    }
}

$io = 0;
#default value for the number of features - overriden if the -f argument is used
@nFeats = ();
#default value for the number of mixtures - overriden if the -m argument is used. 
@numGaussianMixtures = ();
#default value for the number of components per Gaussian - overriden if the -c argument is used
@nComponentsPerGaussian = ();
#don't use random mean by default - overriden by -r
@random = ();
#don't read means and covars from file - overriden by -v
@initialValueFiles = ();
#default prefix for mixture names i.e. gm8
@prefix = ();

if(@ARGV == 0){
    &usage; exit(-1);
}

# parse arguments
while ( $#ARGV >= 0 ) {
    $opt = shift(@ARGV);
    $opt =~ tr/A-Z/a-z/; # ignore case
    if ( $opt =~ /^-m/ ) {
        if ( $#ARGV < 0 ) { #set the number of mixtures
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        &checknnint($opt,$arg);
        push(@numGaussianMixtures, $arg);
    } elsif ( $opt =~ /^-f/ ){ #set the number of features
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        &checknnint($opt,$arg);
        push(@nFeats,$arg);
    } elsif ( $opt =~ /^-i/ ){ #i/o trainable?
	$io = 1;
    } elsif ( $opt =~ /^-p/ ){ #set the mixture name prefix
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        push(@prefix,$arg);
    } elsif ( $opt =~ /^-n/ ){ #set the file containing the name collections
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
	$nameColsFile = $arg;
    } elsif ( $opt =~ /^-o/ ){ #set the name of the output file
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
	$paramsFile = $arg;
    } elsif ( $opt =~ /^-s/ ){ #set the name of the output file
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
	$stateFile = $arg;
    } elsif ( $opt =~ /^-a/ ){ #set the name of this collection of mixtures
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
	push(@names,$arg);
    } elsif ( $opt =~ /^-v/ ){ #name of the file from which to read the means and diag vars
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
	push(@initialValueFiles,$arg);
    } elsif ( $opt =~ /^-c/ ){ #set the number of components per Gaussian 
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        &checknnint($opt,$arg);
        push(@nComponentsPerGaussian,$arg);
    } elsif ( $opt =~ /^-r/ ) {#set random means
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
	#checktf($arg);
        push(@random,$arg);
    } elsif ( $opt =~ /^-help/ ) {
        &usage;
        exit(0);
    } else {
        printf STDERR ("Error: Unknown option (%s)\n",$opt);
        &usage;
        exit(-1);
    }
}

die " no phone state file specified" if(!$stateFile);

open (PHONE, "<$stateFile") || die "cannot open $stateFile";
$i          = 0;
@phone      = ();
%phone      = ();
%startPhone = ();
while(<PHONE>){
    chomp;
    @line            = split " ", $_;
    $phone[$i]       = $line[0];
    $phone{$line[0]} = $line[1];

    if($i > 0){
	$startPhone{$phone[$i]} = $startPhone{$phone[$i-1]} + $phone{$phone[$i-1]};
    } else{
	$startPhone{$phone[$i]} = 0;
    }
    $i++;
}
close PHONE;
#$number_of_states = $startPhone{$phone[$i-1]} + $phone{$phone[$i-1]};

#check args
$numCols = @names;
if($#numGaussianMixtures+1 != $numCols){
    &usage; exit(-1);
}elsif($#nComponentsPerGaussian+1 != $numCols){
    &usage; exit(-1);
}elsif($#prefix+1 != $numCols){
    &usage; exit(-1);
}elsif($#nFeats+1 != $numCols){
    &usage; exit(-1);
}elsif($#random+1 != $numCols){
    &usage; exit(-1);
}elsif(@initialValueFiles != $numCols && @initialValueFiles>0){
    &usage; exit(-1);
}

my @meansList=();
my @diagvarsList=();
foreach $vi (0 .. $#initialValueFiles){
 open(VALS,"<$initialValueFiles[$vi]") or die "Couldn't open $initialValueFiles[$vi] for reading\n";
 my @means=();
 my @diagvars=();
 while(<VALS>){
    my ($f, $m, $v, $rest) = split(' ',$_,4);
    push @means, $m;
    push @diagvars, $v;
 }
 close(VALS);
 die "number of features in file $initialValueFiles[$vi] (".scalar(@means).") does not match feature count on the command line ($nFeats[$vi])" if (@means != $nFeats[$vi]);
 push @meansList, join(" ",@means);
 push @diagvarsList, join(" ",@diagvars);
}

open(PARAMS,">$paramsFile") or die "Couldn't open $paramsFile for writing\n";

#tell the user
print("Tonight, Matthew, I'm going to be generating $numCols collections of mixtures and put the parameters in $paramsFile and the collections in $nameColsFile.  \n");
for($i = 0 ; $i < $numCols ; $i++){
    printf("$names[$i] will contain $numGaussianMixtures[$i] mixtures with $nComponentsPerGaussian[$i] components each and they'll be named $prefix[$i]0 etc. They model $nFeats[$i] features and are initialized %s\n",$random[$i] eq "true" ? "with random means" : $initialValueFiles[$i] ? "from file $initialValueFiles[$i]" : "to mean zero" );
    $totalComponents += $numGaussianMixtures[$i]*$nComponentsPerGaussian[$i];
    $totalMixtures += $numGaussianMixtures[$i];
}
# first do the dense 1d PMFs
$mixNo = 0;
printf PARAMS ("\n%d %% number of DPMFs\n",$totalMixtures);
for($i = 0 ; $i < $numCols ; $i++){
    foreach $phone(@phone){
	for($j = 0; $j < $phone{$phone}; $j++){
	    printf PARAMS ("%d ",$mixNo);
	    printf PARAMS ("mx:%s:%s:%d %d", $names[$i], $phone, $j, $nComponentsPerGaussian[$i]);
	    for ($k=0;$k<$nComponentsPerGaussian[$i];$k++) {
		printf PARAMS (" %f",1.0/$nComponentsPerGaussian[$i]);
	    }
	    print PARAMS "\n";
	    $mixNo++;
	}
    }
}
print PARAMS "\n\n";

#no sparse PMFs
if($io == 1){
   printf PARAMS ("%% no sparse PMFs\n0\n");
   print PARAMS "\n\n";
}

#now for the means
printf PARAMS ("\n%d %% number of Means\n",$totalComponents);
$meanNo = 0;
for ($i=0;$i<$numCols;$i++) {
    foreach $phone(@phone){
	for($j = 0; $j < $phone{$phone}; $j++){
	    for ($k=0;$k<$nComponentsPerGaussian[$i];$k++) {
		printf PARAMS ("%d\nmean:%s:%s:%d %d ", $meanNo, $names[$i], $phone, $j, $nFeats[$i]);
        if($initialValueFiles[$i]){
            print PARAMS $meansList[$i];
        }
        else{
            for ($l=0;$l<$nFeats[$i];$l++) {
                if($random[$i] eq "true"){
                        $val = rand(10000)/10000 - 0.5;
                }
                else{
                        $val = 0;
                }
                print PARAMS "$val ";
            }
        }
		print PARAMS "\n";
		$meanNo++;
	    }
	}
    }
}


print PARAMS "\n\n";

#then the diagonal covariances
printf PARAMS ("\n%d %% number of Covars\n",$totalComponents);
$covarNo = 0;
for ($i=0;$i<$numCols;$i++) {
    foreach $phone(@phone){
        for($j = 0; $j < $phone{$phone}; $j++){
            for ($k=0;$k<$nComponentsPerGaussian[$i];$k++) {
                printf PARAMS ("%d\ncovar:%s:%s:%d %d ", $covarNo, $names[$i], $phone, $j, $nFeats[$i]);
                if($initialValueFiles[$i]){
                    print PARAMS $diagvarsList[$i];
                }
                else{
                    for ($l=0;$l<$nFeats[$i];$l++) {
                    print PARAMS "10.0 ";
                    }
                }
                print PARAMS "\n";
                $covarNo++;
            }
        }
    }
}
if($io == 1){
   print PARAMS "\n%% dlink matrices\n0\n";
   print PARAMS "%% weight matrices\n0\n\n";

   #no Dense CPTs, leave that to gmtkParmConvert
   print PARAMS "%% dense CPTs\n0\n\n";
}

#the Gaussian components...
printf PARAMS ("\n%d %% number of gaussian components\n",$totalComponents);
$compNo = 0;
for ($i=0;$i<$numCols;$i++) {
    foreach $phone(@phone){
        for($j = 0; $j < $phone{$phone}; $j++){
            for ($k=0;$k<$nComponentsPerGaussian[$i];$k++) {
                printf PARAMS ("%d %d\n0 gc:%s:%s:%d \n",$compNo,$nFeats[$i],$names[$i], $phone, $j);
                printf PARAMS ("mean:%s:%s:%d covar:%s:%s:%d\n", $names[$i], $phone, $j, $names[$i], $phone, $j);
                $compNo++;
            }
        }
    }
}
print PARAMS "\n\n";

#finally the mixtures
printf PARAMS ("\n%d %% number of gaussian mixtures\n", $totalMixtures);
$mixNo = 0;
for ($i=0;$i<$numCols;$i++) {
    foreach $phone(@phone){
        for($j = 0; $j < $phone{$phone}; $j++){
            printf PARAMS ("%d %d gm:%s:%s:%d\n",$mixNo, $nFeats[$i], $prefix[$i], $phone, $j);
            printf PARAMS ("%d\n", $nComponentsPerGaussian[$i]);
            printf PARAMS ("mx:%s:%s:%d\n",$names[$i], $phone, $j);
            for ($k=0;$k<$nComponentsPerGaussian[$i];$k++) {
            printf PARAMS ("gc:%s:%s:%d\n",$names[$i], $phone, $j);
            }
            print PARAMS "\n";
            $mixNo++;
        }
	print PARAMS "\n";
    }
}
print PARAMS "\n";

if($io == 1){
   print PARAMS "%% Gaussian Switching Mixtures\n0\n\n";
   print PARAMS "%% Logistic-Regression-based Switching Mixtures\n0\n\n";
   print PARAMS "%% MLP-based Switching Mixtures\n0\n\n";
}

close(PARAMS);

open(COLLECTIONS,">$nameColsFile") or die "Couldn't open $nameColsFile for writing\n";
print COLLECTIONS "$numCols\n\n";
for ($i=0;$i<$numCols;$i++) {
    print COLLECTIONS "$i\n";
    print COLLECTIONS "$names[$i]\n";
    print COLLECTIONS "$numGaussianMixtures[$i]\n";
    foreach $phone(@phone){
        for($j = 0; $j < $phone{$phone}; $j++){
            printf COLLECTIONS  ("gm:%s:%s:%d\n",$prefix[$i], $phone, $j);
        }
    }
    print COLLECTIONS "\n";
}
close(COLLECTIONS);

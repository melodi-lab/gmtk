#!/usr/bin/perl
# written by Ozgur Cetin

sub skip{
    while( <HTK> ){ 
	chomp; 
	$_    =~ s/\"//g;
	$_    =~ tr/[A-Z]/[a-z]/;
	$line = join(" ", split(" ", $_)); 
	if( !( $line eq "" ) ){ last; }
    }
}

sub readTransP{
    @line = split(" ", $_[0]);
    $name = $line[1];

    skip(); @line = split(" ", $line); $size = $line[1]; 
    skip();
    for( $i = 0; $i < ($size-2); $i++ ){ 
	skip(); 
	@line              = split(" ", $line);
	$transP{$name}[$i] = join(" ", @line[($i+1)..($i+2)]); 
	$transPCard++;
    }
    skip();
}

sub readState{
    @line = split(" ", $_[0]);
    $name = "$line[1]\_$state"; $state++;;

    skip(); skip(); $state{$name}[0] = $line;
    skip(); skip(); $state{$name}[1] = $line;
}

sub readHmm{   
    @line = split(" ", $_[0]);
    $name = $line[1]; $hmm++;
    
    skip(); skip(); @line = split(" ", $line); $stateCard = $line[1] - 2;
    for( $i = 0; $i < $stateCard; $i++ ){
	skip();
	skip();
	if( ! ($line =~ /^~/ ) ){ print "Not implemented yet.\n"; exit(1); }
	@line = split(" ", $line); $hmm{$name}[$i] = $line[1];
    }
    skip();
    if( ! ($line =~ /^~/ ) ){ print "Not implemented yet.\n"; exit(1); }
    @line = split(" ", $line); $hmm{$name}[$i] = $line[1];
    skip();
}


$htk        = $ARGV[0];
$hmmMap     = $ARGV[1];
$gmtk       = $ARGV[2];

%hmm    = (); $hmm        = 0;
%state  = (); $state      = 0;
%transP = (); $transPCard = 0;

open HTK, $htk;
while( <HTK> ){
    chomp;
    $_ =~ s/\"//g;
    $_ =~ tr/[A-Z]/[a-z]/;
    $_ = join(" ", split(" ", $_));
    if(    $_ =~ /^~t/ ){ readTransP($_); }
    elsif( $_ =~ /^~s/ ){ readState($_);  }
    elsif( $_ =~ /^~h/ ){ readHmm($_);    }
    elsif( $_ =~ /^<vecsize>/ ){ 	
	$_ =~ s/[^0-9]//g; 
	$d = $_;
    }
}
close HTK;

@state  = sort keys %state;  %stateI = ();
@transP = sort keys %transP;

open GMTK, ">$gmtk";
print GMTK ($#state+1+$transPCard), " % DPMFs\n";
for( $i = 0; $i <= $#state; $i++ ){ 
    $state[$i]  =~ /(.*)\_(.*)/;
    $stateI{$1} = $i; 
    print GMTK "$i p\_$state[$i] 1 1.0\n";
}
foreach $transP( @transP ){
    @transPi =  @{$transP{$transP}};
    for( $j = 0; $j <= $#transPi; $j++ ){
	print GMTK "$i p\_$transP\_$j 2 $transPi[$j]\n";
        $i++;
    }
}
print GMTK "\n\n\n";

print GMTK "$transPCard \% SPMFs\n";
$i = 0;
foreach $transP( @transP ){
    @transPi =  @{$transP{$transP}};
    for( $j = 0; $j <= $#transPi; $j++ ){
	print GMTK "$i s\_$transP\_$j 2 2 0 1 p\_$transP\_$j\n";
        $i++;
    }
}
print GMTK "\n\n\n";

print GMTK ($#state+1), " % means\n";
for( $i = 0; $i <= $#state; $i++ ){ 
    print GMTK "$i m\_$state[$i] $d $state{$state[$i]}[0]\n\n";
}
print GMTK "\n\n";

print GMTK ($#state+1), " % vars\n";
for( $i = 0; $i <= $#state; $i++ ){ 
    print GMTK "$i v\_$state[$i] $d $state{$state[$i]}[1]\n\n";
}
print GMTK "\n\n";

print GMTK "0 \% dlink matrices\n\n\n\n"; 
print GMTK "0 \% weight matrices\n\n\n\n"; 
print GMTK "0 \% dense CPTs\n\n\n\n";

print GMTK ($#state+1), " % gc\n";
for( $i = 0; $i <= $#state; $i++ ){ 
    print GMTK "$i $d 0 gc\_$state[$i] m\_$state[$i] v\_$state[$i]\n";
}
print GMTK "\n\n\n";

print GMTK ($#state+1), " % gm\n";
for( $i = 0; $i <= $#state; $i++ ){ 
    print GMTK "$i $d gm\_$state[$i] 1 p\_$state[$i] gc\_$state[$i]\n";
}
print GMTK "\n\n\n";

print GMTK "0 \% smog\n\n\n\n";
print GMTK "0 \% logistic smog\n\n\n\n";
print GMTK "0 \% mlp smog\n";

close GMTK;

open HMM, ">$hmmMap";
foreach $hmm( keys %hmm ){
    @hmm = @{$hmm{$hmm}};
    print HMM "$hmm ";
    for( $i = 0; $i < $#hmm; $i++ ){ 	
	print HMM "$stateI{$hmm[$i]} "; 
    } 
    print HMM "$hmm[$i]\n";
}
close HMM;

#!/usr/bin/perl
# written by Ozgur Cetin

sub hmmMap{
    %hmm    = ();
    %transP = ();
    open M, $map;
    while( <M> ){
	chomp; 
	@line                  = split(" ", $_);
	$hmm{$line[0]}         = join(" ", @line[1..$#line]);
	$transP{$line[$#line]} = $#line - 1;

	$phn                   = $line[0]; 
	$phn                   =~ s/(.*)\-//;
	$phn                   =~ s/\+(.*)//;
	$trpeq{$line[$#line]}  = $phn;
    }
}

sub skip{
    while( <GMTK> ){
	chomp; 
	$_ = join(" ", split(" ", $_));
	if( !( /^$/ ) ){ if( !( /^\%/ ) ){ last; } }
    }
    s/(\%.*)//;
    $line = join(" ", split(" ", $_));
}

sub getEntity{
    %entity = ();

    skip(); $entityCard = $line;
    for( $i = 0; $i < $entityCard; $i++ ){
	skip(); 
	skip();	
	@line = split(" ", $line); 
	$entity{$line[0]} = join(" ", @line[2..$#line]); #print "$line[0]\n";
    }
}

sub readSpmf{
    %spmf = ();

    skip(); $noOfSpmf = $line;
    for( $i = 0; $i < $noOfSpmf; $i++ ){
	skip(); 
	skip(); @line = split(" ", $line); $name = $line[0];
	skip(); $spmf{$name} = $line;
    }
} 

sub readCPT{
    @cpt = ();
    skip(); $i     = $line;
    skip(); $name  = $line;
    skip(); 
    skip(); @cards = split(" ", $line);

    if( $#cards == 0 ){
	skip(); @cpt = split(" ", $line);
    } elsif( $#cards == 1 ){
	for( $i = 0; $i < $cards[0]; $i++ ){ 
	    skip(); $cpt[$i] = $line;
	}       
    } else{ 
	printf "Doesn't support reading of \> 2D CPT\n"; 
	exit(-1);
    }   
}

sub readGC{
    %gc = ();

    skip(); $noOfGC = $line;
    for( $i = 0; $i < $noOfGC; $i++ ){
	skip(); 
	skip(); 
	skip(); @line = split(" ", $line); $name = $line[1];
	skip(); $gc{$name} = $line;
    }
} 

sub readGM{
    %gm    = ();
    @gm    = ();

    skip(); $noOfGM = $line;
    for( $i = 0; $i < $noOfGM; $i++ ){
	skip(); 
	skip(); $d       = $line;
	skip(); $gm      = $line;
	skip(); @line    = split(" ", $line);
	skip(); $gm{$gm} = join(" ", $line[1], $line);
	$gm[$i]          = $gm;
    }
} 

sub writeState{
    $gm    = $_[0]; 

    @pmfGc = split(" ", $gm{$gm});
    @pmf   = split(" ", $dpmf{$pmfGc[0]});    
    @gc    = @pmfGc[1..$#pmfGc];
    
    print HTK "\~s \"$gm\"\n";
    print HTK "\<NUMMIXES\> ", ($#pmf+1), "\n";
    for( $i = 0; $i <= $#pmf; $i++ ){
	($mean, $var) = split(" ", $gc{$gc[$i]});
	print HTK "\<MIXTURE\> ", ($i+1), " ", $pmf[$i], "\n";
	print HTK "\<MEAN\> $d\n";
	print HTK "$mean{$mean}\n"; 
	print HTK "\<VARIANCE\> $d\n";
	print HTK "$var{$var}\n";
    }
}

sub writeHmm{
    $hmm = $_[0];
    @hmm = split(" ", $hmm{$hmm});

    print HTK "\~h \"$hmm\"\n";
    print HTK "\<BEGINHMM\>\n";
    print HTK "\<NUMSTATES\> ", ($#hmm+2), "\n";
    for( $i = 0; $i < $#hmm; $i++ ){
	print HTK "\<STATE\> ", ($i+2), "\n";
	print HTK "\~s \"$gm[$hmm[$i]]\"\n";
    }
    print HTK "\~t \"$hmm[$i]\"\n";
    print HTK "\<ENDHMM\>\n";
}


$gmtk = $ARGV[0];
$map  = $ARGV[1];
$info = $ARGV[2];
$htk  = $ARGV[3];

hmmMap();

open GMTK, $gmtk; 
getEntity(); %dpmf = %entity;                    # dmpf
readSpmf();                                      # sparse pmf
getEntity(); %mean = %entity;                    # mean
getEntity(); %var  = %entity;                    # var  
getEntity();                                     # dlink
getEntity();                                     # weight matrices

skip(); $noOfCPTs = $line;                       # cpt
for( $k = 0; $k < $noOfCPTs; $k++ ){ 
    readCPT();                                   
} 

readGC();                                        # gaussian components
readGM();                                        # gaussian mixtures
getEntity();                                     # switching mog (smog)
getEntity();                                     # logistic regression smog
getEntity();                                     # mlp-based smog
close GMTK;

open HTK, ">$htk";
print HTK "\~o\n";
print HTK "\<STREAMINFO\> 1 $d\n";
print HTK "\<VECSIZE\> $d \<NULLD\>\<$info\>\n\n";

$currentRow = 0;
foreach $transP( sort keys %transP ){    
    $Q = $transP{$transP};
    print HTK "\~t \"$transP\"\n";
    print HTK "\<TRANSP\> ", ($Q+2), "\n";
    print HTK "0.0 1.0 "; for( $i = 0; $i < $Q; $i++ ){ print HTK "0.0 "; } print HTK "\n";
    for( $i = 0; $i < $Q; $i++ ){
	for( $j = 0; $j < $i+1; $j++ ){	print HTK "0.0 "; }
#	print HTK $dpmf{"p\_$transP\_$i"}," ";
	print HTK $dpmf{"tr\:$trpeq{$transP}\_$i"}," ";
	for( $j = $i+1; $j < $Q; $j++ ){ print HTK " 0.0 "; } print HTK "\n";	
	$currentRow++;
    }
    for( $j = 0; $j < $Q+2; $j++ ){ print HTK "0.0 "; } print HTK "\n";
}
print HTK "\n";

foreach $state( sort keys %gm ){    
    writeState($state);
}
print HTK "\n";

foreach $hmm( sort keys %hmm ){
    writeHmm($hmm);
}
close HTK;

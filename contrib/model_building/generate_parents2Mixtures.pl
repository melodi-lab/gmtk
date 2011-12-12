#!/usr/bin/env perl

## @file
#################################################################
##
## generate_parents2Mixtures.pl
##
## This script takes in a list of parent variable names and
## cardinalities, and outputs initial Gaussians, name collection,
## and DT mapping parent values to mixture numbers.
##
#################################################################

# parse arguments
if (@ARGV < 7) {
    die("Usage: generate_parents2Mixtures.pl output-nameCol output-GMfile num-components num-dimensions u|r prefix parent1 card1 parent2 card2 ...\n");
}

$nameColFile = shift();
$GMFile = shift();
$num_comps = shift();
$num_dims = shift();
$unif_or_rand = shift();
$prefix = shift();

if (!($unif_or_rand eq "u") && !($unif_or_rand eq "r")) {
    die "6th argument must be either 'u' or 'r' (for uniform/random means).  Exiting.\n";
}

$num_parents = 0;
while (@ARGV) {
    $parent[$num_parents] = shift();
    $card[$num_parents++] = shift();
}

$num_mixtures = 1;
foreach $i (0..$num_parents-1) {
    $num_mixtures *= $card[$i];
}

#####################
#
# print the nameCol
#
#####################

open(COL, ">$nameColFile") || die "Couldn't open $nameColFile for writing.  Exiting.\n";

print COL "1\n0\n";
print COL $prefix . "Col\n";
print COL "$num_mixtures\n";

foreach $i (0..$num_mixtures-1) {
    print COL $prefix . "gm_" . &mixnum2parvalues($i, @card) . "\n";
}

close(COL);

#####################
#
# print the GMs
#
#####################

open(GM, ">$GMFile") || die "Couldn't open $GMFile for writing.  Exiting.\n";

print GM "% dense PMFs\n";
print GM "$num_mixtures  % number of DPMFs\n";

foreach $i (0..$num_mixtures-1) {
    print GM "$i  % DPMF num\n";
    print GM $prefix . "gmMx_" . &mixnum2parvalues($i, @card) . " $num_comps";
    foreach $j (0..$num_comps-1) {
        printf GM " %f", 1.0/$num_comps;
    }
    print GM "\n";
}
print GM "\n\n";

print GM "% means\n";
printf GM "%d  %% number of means\n", $num_mixtures*$num_comps;

$mean_num = 0;
foreach $i (0..$num_mixtures-1) {
    foreach $j (0..$num_comps-1) {
        print GM "$mean_num\n";
        printf GM "%smean_%s_%d $num_dims", $prefix, &mixnum2parvalues($i, @card), $j;
        foreach $k (0..$num_dims-1) {
            if ($unif_or_rand eq "u") {
                print GM " 0";
            } else {
                $rn = (rand(10000)/10000 - 0.5)/1000;
                print GM " $rn";
            }
        }
        print GM "\n";
        $mean_num++;
    }
}

print GM "\n\n";
print GM "% covars\n";
printf GM "%d  %% number of covars\n", $num_mixtures*$num_comps;

$covar_num = 0;
foreach $i (0..$num_mixtures-1) {
    foreach $j (0..$num_comps-1) {
        print GM "$covar_num\n";
        printf GM "%scovar_%s_%d $num_dims", $prefix, &mixnum2parvalues($i, @card), $j;
        foreach $k (0..$num_dims-1) {
            print GM " 10.0";
        }
        print GM "\n";
        $covar_num++;
    }
}
print GM "\n\n";
print GM "% Gaussian components\n";
printf GM "%d  %% number of components\n", $num_mixtures*$num_comps;

$comp_num = 0;
foreach $i (0..$num_mixtures-1) {
    foreach $j (0..$num_comps-1) {
        print GM "$comp_num\n";
        print GM "$num_dims\n";
        printf GM "0 %sgc_%s_%d  %% num, dim, type, name\n", $prefix, &mixnum2parvalues($i, @card), $j;
        printf GM "%smean_%s_%d %scovar_%s_%d\n", $prefix, &mixnum2parvalues($i, @card), $j, $prefix, &mixnum2parvalues($i, @card), $j;
        $comp_num++;
    }
}

print GM "\n\n";
print GM "% Mixtures of Gaussians\n";
print GM "$num_mixtures  % number of mixtures\n";

foreach $i (0..$num_mixtures-1) {
    print GM "$i $num_dims $prefix" . "gm_" . &mixnum2parvalues($i, @card) . "  % num, dim, name\n";
    print GM "$num_comps  % num components\n";
    print GM "$prefix" . "gmMx_" . &mixnum2parvalues($i, @card) . "\n";
    foreach $j (0..$num_comps-1) {
        printf GM "%sgc_%s_%d ", $prefix, &mixnum2parvalues($i, @card), $j;
    }
    print GM "\n";
}

close(GM);

sub mixnum2parvalues {
    my ($mixnum, @card) = @_;

    $num_parents = @card;

    $pvals[$num_parents-1] = mod($mixnum,$card[$num_parents-1]);
    
    $num = int(($mixnum - $pvals[$num_parents-1])/$card[$num_parents-1]);
    $divider = 1;
    foreach $p (reverse 0..$num_parents-2) {
	$pvals[$p] = mod($num,$card[$p]);
	$num = int(($num - $pvals[$p])/$card[$p]);
    }

    $pval_string = "$pvals[0]";
    foreach $p (1..$num_parents-1) {
        $pval_string .= "_$pvals[$p]";
    }

#    print $pval_string, "\n";

    return $pval_string;
}

sub mod {
    my ($n1, $n2) = @_;

    return $n1 - $n2 * int($n1/$n2);
}

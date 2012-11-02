#!/usr/bin/perl -w

## @file

sub usage {
 printf STDERR (
  "Description: Generate a file containing a DT for a wordTransition in the phone-free model.\n" .
  "Usage: progname [options]\n" .
   "   -d file            dictionary file\n" .
   "   -dtName dtname     name of the decision tree\n" .
   "   -n #               num features\n" .
   "   -m #               max. num positions in a word (for now, needs to be 1 more than the actual max length of any word in the dict)\n"
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

sub checkfile {
    ($opt,$arg) = @_;
    if (( $arg !~ /^\S+$/ )) {
        print STDERR ("Error: Value \"", $arg, "\" invalid for option ",
                      $opt, " (filename expected)\n");
        &usage; exit(-1);
    }
}

# parse arguments
if ( $#ARGV < 2 ) {
    &usage;
    exit(-1);
}
while ( $#ARGV >= 0 ) {
    $opt = shift(@ARGV);
    $opt =~ tr/A-Z/a-z/; # ignore case
    if ( $opt =~ /^-d$/ ) {
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        &checkfile($opt,$arg);
        $dict = $arg;
    } elsif ( $opt =~ /^-dtName$/i ) {
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        $dtName = $arg;
    } elsif ( $opt =~ /^-n/ ) {
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        &checknnint($opt,$arg);
        $num_features = $arg;
    } elsif ( $opt =~ /^-m/ ) {
        if ( $#ARGV < 0 ) {
            print STDERR ("Error: Option ", $opt, " requires an argument\n");
            &usage; exit(-1);
        }
        $arg = shift(@ARGV);
        &checknnint($opt,$arg);
        $max_num_pos_per_word = $arg;
    } elsif ( $opt =~ /^-help/ ) {
        &usage;
        exit(0);
    } else {
        printf STDERR ("Error: Unknown option (%s)\n",$opt);
        &usage;
        exit(-1);
    }
}


# get num phones in each word
open(DICT, $dict) || die "Couldn't open dict file $dict for reading.  Exiting.\n";

$num_words = 0;
$last_wd = '';
foreach $line (<DICT>) {
    chomp($line);
    ($wd, $pron) = split(/\s*\:\s*/, $line);
    @phones = split(/ /, $pron);

    if (! ($last_wd eq $wd)) {
	$num_words++;
	$last_wd = $wd;
    }

    $wd_num = $num_words-1;

    # keep track of how many prons each word has
    $num_prons[$wd_num] += 1;

    # word identities
    $wd_label[$wd_num] = $wd;
    $wd_label[$wd_num] =~ s/\'/\_/g;

    $pron_num = $num_prons[$wd_num] - 1;

    # keep track of how many phones each variant has
    $num_phones{$wd_num."_".$pron_num} = @phones;

}
close(DICT) || die "Couldn't close dict file $dict. Exiting.\n";

# calculate number of parents = 
# word + variant + position and transition vars for each feature
$num_parents = 2 + 2 * $num_features;  

# note: parents 0, ... $num_features-1 are feature transitions
#       parents $num_features, ... , 2*$num_features-1 are feature positions
#       parent 2*$num_features is word,
#       parent 2*$num_features+1 is wordVariant,

$wd_parent_num = 2*$num_features;
$variant_parent_num = $wd_parent_num+1;

print "%generated from dict $dict on ".`date`."\n";
print "1 % number of DTs in this file\n";
print "0 % DT number\n";

print "$dtName % DT name\n";
print "$num_parents % num parents\n";
$num_tabs = 1; # keep track of level of tree
foreach $i (0..$num_features-1) {               # now check that every feature has a transition
    foreach $j (1..$num_tabs) { # enter the right number of tabs
	print "  ";
    }
    print "$i 2 0 default\n";

    foreach $j (1..$num_tabs) { # enter the right number of tabs
	print "  ";
    }

    print "  -1 0\n";                           # if this feature doesn't have transition, then no wd transition
    $num_tabs++;
}

foreach $i (1..$num_tabs) { # enter the right number of tabs
    print "  ";
}

print "$wd_parent_num $num_words ";

foreach $i (0..$num_words-2) {
    print "$i ";
}

print "default\n";

$initial_num_tabs = $num_tabs;
foreach $i (0..$num_words-1) {
    
    $num_tabs = $initial_num_tabs;

    foreach $j (0..$num_prons[$i]-1) {
	$last_position_in_word[$j] = $num_phones{$i."_".$j} - 1;
	$one_before_last[$j] = $last_position_in_word[$j] - 1;
	$one_after_last[$j] = $last_position_in_word[$j] + 1;
    }

    # list branches for pron variants for this word
    foreach $j (1..$num_tabs) { # enter the right number of tabs
	print "  ";
    }

    print "$variant_parent_num $num_prons[$i] ";

    foreach $j (0..$num_prons[$i]-2) {
	print "$j ";
    }

    print "default  % wd num $i, wd $wd_label[$i]\n";

    foreach $j (0..$num_prons[$i]-1) {

	$num_tabs = $initial_num_tabs+1;

	foreach $k (0..$num_features-1) {
	    foreach $l (1..$num_tabs) { # enter the right number of tabs
		print "  ";
	    }
	    
	    $parent_num = $num_features + $k;
	    
	    # list branches for current feature position variable
	    if ($last_position_in_word[$j] > 0) {
		print "  $parent_num 3 0:$one_before_last[$j] $one_after_last[$j]:", $max_num_pos_per_word-1, " default\n";

# Failed attempt to fix bug that occurs when last_position_in_word >= max_num_pos_per_word-1
# Should try to implement this properly sometime
#
#		if ($last_position_in_word[$j] < $max_num_pos_per_word-1) {
#		    print "  $parent_num 3 0:$one_before_last[$j] $one_after_last[$j]:", $max_num_pos_per_word-1, " default\n";
#		} else {
#		    print "  $parent_num 3 0:$one_before_last[$j] default\n";
#		}

	    } else {
		print "  $parent_num 1 default\n";
	    }
	    
	    if ($last_position_in_word[$j] > 0) {
		foreach $l (1..$num_tabs) { # enter the right number of tabs
		    print "  ";
		}	
		print "    -1 0\n";
		foreach $l (1..$num_tabs) { # enter the right number of tabs
		    print "  ";
		}	
		print "    -1 0\n";
	    }
	    
	    $num_tabs++;
	}

	foreach $l (1..$num_tabs) { # enter the right number of tabs
	    print "  ";
	}
	print "  -1 1\n";
    }
}

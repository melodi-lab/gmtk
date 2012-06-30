#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;
use File::Path;

## @file

my $exe = basename $0;
my $usage = "$exe [options] <transcript directories>\n\n".
    "  OPTIONS:\n".
    "   -silence    - ensure silence consistency.  if any of pl1, pl2, dg1, or dg2 is marked as SIL,\n".
    "                 then make them all SIL.\n".
    "   -order      - ensure intra-transcriber constriction order.  swap if out of order.\n".
    "   -swap       - ensure inter-transcriber constriction alignment.  if swapping the constrictions\n".
    "                 of one transcriber will increase the constriction place alignments then do so.\n".
    "   -recombine  - post process the transcriptions using recombineLabels.pl to combine\n".
    "                 sequential labels if they are the same\n".
    "   -chunk      - ensure that all the feature tiers have the same time boundaries.  basically the\n".
    "                 opposite of recombine.\n".
    "   -agreement  - calculate percent agreement and Kappa statistic\n".
    "   -matrix     - display agreement/confusion matrix between transcribers\n".
    "   -unnormalize  don't normalize the agreement/confusion matrix\n".
    "   -save=<DIR> - save modified transcriptions to DIR (replace top level directory w/ what?)\n".
    "   -quiet      - turn off some debugging information\n";
die $usage if $#ARGV < 0;

my @constrictionFeatures = ('pl1','pl2','dg1','dg2');
my @otherFeatures = ('glo','nas','rd','vow');
my @allFeatures = (@otherFeatures,@constrictionFeatures);
my (%baseFiles, %validLabels, %transcriptions,%agreementMatrix,
    %labelTime, %observedAgreement, %expectedAgreement, %kappa);
my $recombineCommand = "/home/ws05/nborges/jhu/ws06/scripts/transCombineLabels.pl";
my $valuesDir = "/export/ws06afsr/data/manual_trans.newdir/values";
my (@directories,@transcribers);

# parse command line options
my ($optSilence,$optOrder,$optSwap,$optRecombine,$optChunk,
    $optAgreement,$optMatrix,$optUnnormalize,$optSave,$optDebug) = (0,0,0,0,0,0,0,0,0,1);
my $optSaveDir;
foreach (@ARGV) {
  if (/^-silence$/i) {
    $optSilence = 1;
  }
  elsif (/^-order$/i) {
    $optOrder = 1;
  }
  elsif (/^-swap$/i) {
    $optSwap = 1;
  }
  elsif (/^-recombine$/i) {
    $optRecombine = 1;
  }
  elsif (/^-chunk$/i) {
    $optChunk = 1;
  }
  elsif (/^-agreement$/i) {
    $optAgreement = 1;
  }
  elsif (/^-matrix$/i) {
    $optMatrix = 1;
  }
  elsif (/^-unnormalize$/i) {
    $optUnnormalize = 1;
  }
  elsif (/^-save=(.*)$/i) {
    $optSaveDir = $1;
    $optSave = 1;
  }
  elsif (/^-quiet$/i) {
    $optDebug = 0;
  }
  else {
    if (-e $_) {
      push @directories, $_;
    } else {
      die "$usage\nINVALID ARGUMENT: $_\n";
    }
  }
}

die "$usage\nINVALID OPTIONS: recombine and chunk are mutually exclusive\n" if ($optRecombine and $optChunk);
die "$usage\nINVALID OPTIONS: order and swap are mutually exclusive\n" if ($optOrder and $optSwap);

# load valid feature labels into hash table so we can
# easily test later using 'exists'
foreach my $feature (@allFeatures) {
  my $valueFile = "$valuesDir/values.$feature";
  die "$valueFile not found\n" unless (-e $valueFile);
  open FILE, "<$valueFile";
  while (<FILE>) {
    chomp;
    unless (/^\s*$/) {
      $validLabels{$feature}{$_} = 1;
    }
  }
  close FILE;
}

my %transcribersHash;
# store file information
foreach my $feature (@allFeatures) {
  foreach my $directory (@directories) {
    my @files = `find $directory -follow -name "*.$feature"`;
    foreach my $file (@files) {
      chomp $file;
      my $dir = dirname $file;
      $dir =~ /\/([^\/\.]{2})(?:\/|$)/;
      my $transcriber = $1;
      die "transcriber not found" if ($transcriber =~ /^\s*$/);
      $transcribersHash{$transcriber} = 1 if (not exists $transcribersHash{$transcriber});
      my $base = basename $file;
      $base =~ s/\.(\w+)$//;
      my $feature = $1;
      $baseFiles{$base}{$feature}{$transcriber} = $file if (not exists $baseFiles{$base}{$feature}{$transcriber});
    }
  }
}
@transcribers = sort keys %transcribersHash;

die "not enough transcribers found\n" if ($optAgreement and @transcribers < 2);
die "not two transcribers\n" if ($optMatrix and @transcribers != 2);

# read in transcription information
foreach my $baseFile (keys %baseFiles) {

  # check to make sure all files exist from both transcribers
  my $allFilesExist = 1;
  foreach my $transcriber (@transcribers) {
    foreach my $feature (@allFeatures) {
      if (not exists $baseFiles{$baseFile}{$feature}{$transcriber} or
          !-e $baseFiles{$baseFile}{$feature}{$transcriber}) {
        $allFilesExist = 0;
      }
    }
  }
  if (!$allFilesExist) {
    print STDERR "WARNING: skipping $baseFile due to incomplete transcriptions\n";
    next;
  }

  # store start,stop,label data in transcriptions
  foreach my $feature (keys %{$baseFiles{$baseFile}}) {
    foreach my $transcriber (keys %{$baseFiles{$baseFile}{$feature}}) {
      my $fullFeatureFile = $baseFiles{$baseFile}{$feature}{$transcriber};
      open FILE, "<$fullFeatureFile";
      while (<FILE>) {
        next if (m/^\s*$/);
        chomp;
        my ($start,$stop,$label) = split;
        if ((!$start && $start != 0) or (!$stop && $stop != 0)) {
          print STDERR "WARNING: failed to split '$_' in $fullFeatureFile\n";
          next;
        }
        # convert samples to seconds if there's only digits w/o a decimal
        if ( ($start =~ /^\d+$/) and ($stop =~ /^\d+$/) ) {
          $start /= 8000.0;
          $stop /= 8000.0
        }
        $start = sprintf("%.7f",$start);
        $stop = sprintf("%.7f",$stop);

        $label =~ s#N/A#N-A#;
        my @labels = split /\//, $label;
        #print STDERR "WARNING: multiple labels '$label' from $start-$stop in $fullFeatureFile\n" if ($label =~ /\//);
        foreach my $lab (@labels) {
          $lab =~ s#N-A#N/A#;
          if (!exists $validLabels{$feature}{$lab}) {
            print STDERR "WARNING: invalid label $lab from $start-$stop in $fullFeatureFile\n" if $optDebug;
          }
          $transcriptions{$baseFile}{$feature}{$transcriber}{$start}{LABELS}{$lab} = 1;
        }
        $transcriptions{$baseFile}{$feature}{$transcriber}{$start}{STOP} = $stop;
      }
      close FILE;
    }
  }
}

# split appropriate time boundaries
foreach my $baseFile (keys %baseFiles) {
  foreach my $feature (keys %{$baseFiles{$baseFile}}) {
    foreach my $transcriber (@transcribers) {
      my %newMarks;

      # go through all of the other features for this transcriber
      if ($optChunk) {
        foreach my $otherFeature (@allFeatures) {
          next if ($otherFeature eq $feature);
          die "missing $baseFile $otherFeature $transcriber\n" if not exists $transcriptions{$baseFile}{$otherFeature}{$transcriber};
          foreach my $start (sort {$a<=>$b} keys %{$transcriptions{$baseFile}{$otherFeature}{$transcriber}}) {
            $newMarks{$start} = 1 if not exists $newMarks{$start};
            my $stop = $transcriptions{$baseFile}{$otherFeature}{$transcriber}{$start}{STOP};
            $newMarks{$stop} = 1 if not exists $newMarks{$stop};
          }
        }
      }

      # go through all of the other constriction features for this transcriber
      if ($optOrder or $optSilence or $optSwap) {
        my $isConstriction = 0;
        foreach my $constriction (@constrictionFeatures) {
          $isConstriction = 1 if ($feature eq $constriction);
        }
        if ($isConstriction) {
          foreach my $constriction (@constrictionFeatures) {
            next if ($constriction eq $feature);
            die "missing $baseFile $constriction $transcriber\n" if not exists $transcriptions{$baseFile}{$constriction}{$transcriber};
            foreach my $start (sort {$a<=>$b} keys %{$transcriptions{$baseFile}{$constriction}{$transcriber}}) {
              $newMarks{$start} = 1 if not exists $newMarks{$start};
              my $stop = $transcriptions{$baseFile}{$constriction}{$transcriber}{$start}{STOP};
              $newMarks{$stop} = 1 if not exists $newMarks{$stop};
            }
          }
        }
      }

      # go through all of the other transcribers
      if ($optAgreement or $optMatrix or $optSwap) {
        foreach my $otherTranscriber (keys %{$baseFiles{$baseFile}{$feature}}) {
          next if ($otherTranscriber eq $transcriber);

          # go through all of the other transcribers' transcriptions of this feature
          die "missing $baseFile $feature $otherTranscriber\n" if not exists $transcriptions{$baseFile}{$feature}{$otherTranscriber};
          foreach my $otherStart (sort {$a<=>$b} keys %{$transcriptions{$baseFile}{$feature}{$otherTranscriber}}) {
            $newMarks{$otherStart} = 1 if not exists $newMarks{$otherStart};
            my $otherStop = $transcriptions{$baseFile}{$feature}{$otherTranscriber}{$otherStart}{STOP};
            $newMarks{$otherStop} = 1 if not exists $newMarks{$otherStop};
          }

          # go through all of the other transcribers' transcriptions of the other constriction features
          if ($optSwap) {
            my $isConstriction = 0;
            foreach my $constriction (@constrictionFeatures) {
              $isConstriction = 1 if ($feature eq $constriction);
            }
            if ($isConstriction) {
              foreach my $constriction (@constrictionFeatures) {
                die "missing $baseFile $constriction $otherTranscriber\n" if not exists $transcriptions{$baseFile}{$constriction}{$otherTranscriber};
                foreach my $otherStart (sort {$a<=>$b} keys %{$transcriptions{$baseFile}{$constriction}{$otherTranscriber}}) {
                  $newMarks{$otherStart} = 1 if not exists $newMarks{$otherStart};
                  my $otherStop = $transcriptions{$baseFile}{$constriction}{$otherTranscriber}{$otherStart}{STOP};
                  $newMarks{$otherStop} = 1 if not exists $newMarks{$otherStop};
                }
              }
            }
          }
        }
      }

      # actually do mark splitting
      foreach my $newMark (keys %newMarks) {
        # find our first and last mark (they might change as we go along)
        my @starts = sort {$a<=>$b} keys %{$transcriptions{$baseFile}{$feature}{$transcriber}};
        my $firstStart = $starts[0];
        my $lastStop = $transcriptions{$baseFile}{$feature}{$transcriber}{$starts[$#starts]}{STOP};
        if (not exists $transcriptions{$baseFile}{$feature}{$transcriber}{$newMark} and $newMark != $lastStop) {
          if ($newMark < $firstStart) {
            $transcriptions{$baseFile}{$feature}{$transcriber}{$newMark}{STOP} = $firstStart;
            $transcriptions{$baseFile}{$feature}{$transcriber}{$newMark}{LABELS}{"NULL"} = 1;
          } 
          elsif ($newMark > $lastStop) {
            $transcriptions{$baseFile}{$feature}{$transcriber}{$lastStop}{STOP} = $newMark;
            $transcriptions{$baseFile}{$feature}{$transcriber}{$lastStop}{LABELS}{"NULL"} = 1;
          }
          else {
            foreach my $start (@starts) {
              my $stop = $transcriptions{$baseFile}{$feature}{$transcriber}{$start}{STOP};
              if ($start < $newMark and $stop > $newMark) {
                $transcriptions{$baseFile}{$feature}{$transcriber}{$newMark}{STOP} = $stop;
                $transcriptions{$baseFile}{$feature}{$transcriber}{$start}{STOP} = $newMark;
                $transcriptions{$baseFile}{$feature}{$transcriber}{$newMark}{LABELS} =
                    \%{$transcriptions{$baseFile}{$feature}{$transcriber}{$start}{LABELS}};
              }
            }
          }
        }
      }
    }
  }
}

validateLabels();

if ($optSilence) {
  foreach my $baseFile (keys %transcriptions) {
    my $fPlaceHolder = $constrictionFeatures[0];
    foreach my $transcriber (keys %{$transcriptions{$baseFile}{$fPlaceHolder}}) {
      foreach my $start (sort {$a<=>$b} keys %{$transcriptions{$baseFile}{$fPlaceHolder}{$transcriber}}) {
        my $stop = $transcriptions{$baseFile}{$fPlaceHolder}{$transcriber}{$start}{STOP};
        die if (not exists $transcriptions{$baseFile}{'pl1'}{$transcriber}{$start}{STOP});
        die if (not exists $transcriptions{$baseFile}{'pl2'}{$transcriber}{$start}{STOP});
        die if (not exists $transcriptions{$baseFile}{'dg1'}{$transcriber}{$start}{STOP});
        die if (not exists $transcriptions{$baseFile}{'dg2'}{$transcriber}{$start});
        die if (not exists $transcriptions{$baseFile}{'dg2'}{$transcriber}{$start}{STOP});

        my @pl1Labels = keys %{$transcriptions{$baseFile}{'pl1'}{$transcriber}{$start}{LABELS}};
        my @pl2Labels = keys %{$transcriptions{$baseFile}{'pl2'}{$transcriber}{$start}{LABELS}};
        my @dg1Labels = keys %{$transcriptions{$baseFile}{'dg1'}{$transcriber}{$start}{LABELS}};
        my @dg2Labels = keys %{$transcriptions{$baseFile}{'dg2'}{$transcriber}{$start}{LABELS}};
        my @labels = (@pl1Labels,@pl2Labels,@dg1Labels,@dg2Labels);
        if ((($#pl1Labels == 0 and $pl1Labels[0] eq 'SIL') or
             ($#pl2Labels == 0 and $pl2Labels[0] eq 'SIL') or
             ($#dg1Labels == 0 and $dg1Labels[0] eq 'SIL') or
             ($#dg2Labels == 0 and $dg2Labels[0] eq 'SIL')) and
            !(($#pl1Labels == 0 and $pl1Labels[0] eq 'SIL') and
              ($#pl2Labels == 0 and $pl2Labels[0] eq 'SIL') and
              ($#dg1Labels == 0 and $dg1Labels[0] eq 'SIL') and
              ($#dg2Labels == 0 and $dg2Labels[0] eq 'SIL'))) {
          my %SILENCE = ('SIL' => 1);
          $transcriptions{$baseFile}{'pl1'}{$transcriber}{$start}{LABELS} = \%SILENCE;
          $transcriptions{$baseFile}{'pl2'}{$transcriber}{$start}{LABELS} = \%SILENCE;
          $transcriptions{$baseFile}{'dg1'}{$transcriber}{$start}{LABELS} = \%SILENCE;
          $transcriptions{$baseFile}{'dg2'}{$transcriber}{$start}{LABELS} = \%SILENCE;
          my @newlabels = (keys %{$transcriptions{$baseFile}{'pl1'}{$transcriber}{$start}{LABELS}},
                           keys %{$transcriptions{$baseFile}{'pl2'}{$transcriber}{$start}{LABELS}},
                           keys %{$transcriptions{$baseFile}{'dg1'}{$transcriber}{$start}{LABELS}},
                           keys %{$transcriptions{$baseFile}{'dg2'}{$transcriber}{$start}{LABELS}});
          print STDERR "WARNING: ${transcriber}'s inconsistent silence from $start-$stop in $baseFile has been corrected to @newlabels\n" if $optDebug;
        }
      }
    }
  }
}

if ($optOrder) {
  my %constrictionDepth = ('LAB' => 0,
                           'L-D' => 1,
                           'DEN' => 2,
                           'ALV' => 3,
                           'P-A' => 4,
                           'VEL' => 5,
                           'GLO' => 6,
                           'NONE'=> 7);
  foreach my $baseFile (keys %transcriptions) {
    my $fPlaceHolder = $constrictionFeatures[0];
    foreach my $transcriber (keys %{$transcriptions{$baseFile}{$fPlaceHolder}}) {
      foreach my $start (sort {$a<=>$b} keys %{$transcriptions{$baseFile}{$fPlaceHolder}{$transcriber}}) {
        my $stop = $transcriptions{$baseFile}{$fPlaceHolder}{$transcriber}{$start}{STOP};
        my @pl1Labels = keys %{$transcriptions{$baseFile}{'pl1'}{$transcriber}{$start}{LABELS}};
        my @pl2Labels = keys %{$transcriptions{$baseFile}{'pl2'}{$transcriber}{$start}{LABELS}};
        if ($#pl1Labels == 0 and
            $#pl2Labels == 0 and
            exists $constrictionDepth{$pl1Labels[0]} and
            exists $constrictionDepth{$pl2Labels[0]}) {
          if ($constrictionDepth{$pl2Labels[0]} - $constrictionDepth{$pl1Labels[0]} < 0 ) {
            print STDERR "WARNING: ${transcriber}'s pl1:$pl1Labels[0] and pl2:$pl2Labels[0] have been reordered from $start-$stop in $baseFile\n" if $optDebug;
            my $pl1Ref = $transcriptions{$baseFile}{'pl1'}{$transcriber}{$start}{LABELS};
            my $pl2Ref = $transcriptions{$baseFile}{'pl2'}{$transcriber}{$start}{LABELS};
            my $dg1Ref = $transcriptions{$baseFile}{'dg1'}{$transcriber}{$start}{LABELS};
            my $dg2Ref = $transcriptions{$baseFile}{'dg2'}{$transcriber}{$start}{LABELS};
            $transcriptions{$baseFile}{'dg1'}{$transcriber}{$start}{LABELS} = $dg2Ref;
            $transcriptions{$baseFile}{'dg2'}{$transcriber}{$start}{LABELS} = $dg1Ref;
            $transcriptions{$baseFile}{'pl1'}{$transcriber}{$start}{LABELS} = $pl2Ref;
            $transcriptions{$baseFile}{'pl2'}{$transcriber}{$start}{LABELS} = $pl1Ref;
          }
        }
      }
    }
  }
}

    # check if swapping the transcribers places of constriction would increase the alignment
if ($optSwap) {
  foreach my $baseFile (keys %transcriptions) {
    my $fPlaceHolder = $constrictionFeatures[0];
    my $tPlaceHolder = $transcribers[0];
    foreach my $start (sort {$a<=>$b} keys %{$transcriptions{$baseFile}{$fPlaceHolder}{$tPlaceHolder}}) {
      my $stop = $transcriptions{$baseFile}{$fPlaceHolder}{$tPlaceHolder}{$start}{STOP};

      # compute the normal agreement measures (sans swapping)
      my $normalTotalTimeInAgreement = 0;
      my (@normalAgreements,@normalNotAgreements);
      foreach my $el (0..1) {   # the first two constriction features are pl1 and pl2
        ($normalAgreements[$el],$normalNotAgreements[$el]) =
            computeAgreement($baseFile,$constrictionFeatures[$el],$constrictionFeatures[$el],$tPlaceHolder,$start);
        $normalTotalTimeInAgreement += $normalAgreements[$el];
      }

      # compute the agreement for transcribers whose features can be swapped
      my $swappedTotalTimeInAgreement = 0;
      my (@swappedAgreements,@swappedNotAgreements,@swappableTranscribers,@swappableTranscriberScores);
      foreach my $el (0..$#constrictionFeatures) {
        ($swappedAgreements[$el],$swappedNotAgreements[$el]) = (0,0);
      }
      foreach my $transcriber (@transcribers) {
        my @pl2Labels = keys %{$transcriptions{$baseFile}{'pl2'}{$transcriber}{$start}{LABELS}};
        my @dg2Labels = keys %{$transcriptions{$baseFile}{'dg2'}{$transcriber}{$start}{LABELS}};
        # if this transcriber's marks can be swapped, calculated what the agreement would be
        if ($#pl2Labels == 0 and $#dg2Labels == 0 and $pl2Labels[0] eq "NONE" and $dg2Labels[0] eq "VOW") {
          push @swappableTranscribers, $transcriber;
          ($swappedAgreements[0],$swappedNotAgreements[0]) =
              computeAgreement($baseFile,$constrictionFeatures[1],$constrictionFeatures[0],$transcriber,$start);
          ($swappedAgreements[1],$swappedNotAgreements[1]) =
              computeAgreement($baseFile,$constrictionFeatures[0],$constrictionFeatures[1],$transcriber,$start);
          $swappedTotalTimeInAgreement = $swappedAgreements[0] + $swappedAgreements[1];
          push @swappableTranscriberScores, $swappedTotalTimeInAgreement;
        }
      }

      # find the best transcriber to swap
      my $indexTranscriberToSwap = -1;
      foreach my $index (0..$#swappableTranscribers) {
        if ($swappableTranscriberScores[$index] > $normalTotalTimeInAgreement) {
          if ( ($indexTranscriberToSwap == -1) or
               ($indexTranscriberToSwap != -1 and $swappableTranscriberScores[$index] > $swappableTranscriberScores[$indexTranscriberToSwap]) ) {
            $indexTranscriberToSwap = $index;
          }
        }
      }

      # actually swap the labels
      if ($indexTranscriberToSwap != -1) {
        my $transcriber = $swappableTranscribers[$indexTranscriberToSwap];
        my $pl1Ref = $transcriptions{$baseFile}{'pl1'}{$transcriber}{$start}{LABELS};
        my $pl2Ref = $transcriptions{$baseFile}{'pl2'}{$transcriber}{$start}{LABELS};
        my $dg1Ref = $transcriptions{$baseFile}{'dg1'}{$transcriber}{$start}{LABELS};
        my $dg2Ref = $transcriptions{$baseFile}{'dg2'}{$transcriber}{$start}{LABELS};
        my $pl1 = (keys %{$pl1Ref})[0];
        print STDERR "WARNING: swapped ${transcriber}'s single constriction pl1:$pl1 to match the other transcribers pl2 from $start-$stop in $baseFile\n" if $optDebug;
        $transcriptions{$baseFile}{'dg1'}{$transcriber}{$start}{LABELS} = $dg2Ref;
        $transcriptions{$baseFile}{'dg2'}{$transcriber}{$start}{LABELS} = $dg1Ref;
        $transcriptions{$baseFile}{'pl1'}{$transcriber}{$start}{LABELS} = $pl2Ref;
        $transcriptions{$baseFile}{'pl2'}{$transcriber}{$start}{LABELS} = $pl1Ref;
      }
    }
  }
}

if ($optSave) {
  # print out hash tree of the marks w/ all the same time boundaries
  $" = "\/";
  foreach my $baseFile (keys %transcriptions) {
    foreach my $feature (keys %{$transcriptions{$baseFile}}) {
      foreach my $transcriber (keys %{$transcriptions{$baseFile}{$feature}}) {
        my $fullFile = $baseFiles{$baseFile}{$feature}{$transcriber};
        $fullFile =~ s/[^\/]*trans[^\/]*/$optSaveDir/;
        my $dir = dirname $fullFile;
        mkpath($dir,0,0775) if (!-e $dir);
        open FILE, ">$fullFile";
        foreach my $start (sort {$a<=>$b} keys %{$transcriptions{$baseFile}{$feature}{$transcriber}}) {
          my $stop = $transcriptions{$baseFile}{$feature}{$transcriber}{$start}{STOP};
          my @labels = keys %{$transcriptions{$baseFile}{$feature}{$transcriber}{$start}{LABELS}};
          next if ($#labels == 0 and $labels[0] eq "NULL");
          print FILE "$start ";
          print FILE "$stop ";
          print FILE "@labels\n";
        }
        close FILE;
        `$recombineCommand $fullFile` if $optRecombine;
      }
    }
  }
}

if ($optAgreement or $optMatrix) {
  # accumulate agreements
  foreach my $baseFile (keys %transcriptions) {
    foreach my $feature (keys %{$transcriptions{$baseFile}}) {
      foreach my $transcriber (keys %{$transcriptions{$baseFile}{$feature}}) {
        foreach my $start (sort {$a<=>$b} keys %{$transcriptions{$baseFile}{$feature}{$transcriber}}) {
          storeLabelTime($baseFile,$feature,$transcriber,$start);
          storeAgreementMatrix($baseFile,$feature,$start) if $optMatrix;
        }
      }
    }
    my $tPlaceHolder = $transcribers[0];
    foreach my $feature (@allFeatures) {
      foreach my $start (sort {$a<=>$b} keys %{$transcriptions{$baseFile}{$feature}{$tPlaceHolder}}) {
        my ($timeInAgreement,$timeNotInAgreement) = computeAgreement($baseFile,$feature,$feature,$tPlaceHolder,$start);
        storeObservedFeatureAgreement($feature,$timeInAgreement,$timeNotInAgreement);
      }
    }
    foreach my $feature (@allFeatures) {
      $observedAgreement{$feature}{RATIO} = $observedAgreement{$feature}{TIMEAGREE} / $observedAgreement{$feature}{TOTAL};
    }
  }


  printPercentAgreement() if $optAgreement;
  printAgreementMatrix() if $optMatrix;
}

sub printPercentAgreement {
  # print out final accumulations and percentages
  print "|  feature  |  timeInAgreement  |  timeNotInAgreement  |  obsAgreement  |  expAgreement  |  Kappa  |\n";
  my ($agreeTotal, $dontTotal, $timeTotal) = (0,0,0);
  $expectedAgreement{OVERALL}{TIME} = 0;
  foreach my $feature (sort keys %observedAgreement) {
    foreach my $label (keys %{$validLabels{$feature}}) {
      foreach my $transcriber (@transcribers) {
        $labelTime{$feature}{$label}{$transcriber} = 0 if (not exists $labelTime{$feature}{$label}{$transcriber});
        if (exists $expectedAgreement{$feature}{$label}) {
          $expectedAgreement{$feature}{$label} *= $labelTime{$feature}{$label}{$transcriber};
        } else {
          $expectedAgreement{$feature}{$label} = $labelTime{$feature}{$label}{$transcriber};
        }
      }
      $expectedAgreement{$feature}{$label} /= $observedAgreement{$feature}{TOTAL};
      if (exists $expectedAgreement{$feature}{TIME}) {
        $expectedAgreement{$feature}{TIME} += $expectedAgreement{$feature}{$label};
      } else {
        $expectedAgreement{$feature}{TIME} = $expectedAgreement{$feature}{$label};
      }
    }
    $expectedAgreement{$feature}{RATIO} = $expectedAgreement{$feature}{TIME} / $observedAgreement{$feature}{TOTAL};
    my $agree = (exists $observedAgreement{$feature}{TIMEAGREE}) ? $observedAgreement{$feature}{TIMEAGREE} : 0;
    my $dont = (exists $observedAgreement{$feature}{TIMEDONTAGREE}) ? $observedAgreement{$feature}{TIMEDONTAGREE} : 0;
    $agreeTotal += $agree;
    $dontTotal += $dont;
    if ($agree + $dont > 0) {
      $observedAgreement{$feature}{RATIO} = $agree / ($agree + $dont);
      $kappa{$feature} = ($observedAgreement{$feature}{TIMEAGREE} - $expectedAgreement{$feature}{TIME}) /
          ($observedAgreement{$feature}{TOTAL} - $expectedAgreement{$feature}{TIME});
      printf("|  %5s  |  %6.1f  |  %6.1f  |  %6.3f  |  %6.3f  |  %6.3f  |\n",$feature,$agree,$dont,
             $observedAgreement{$feature}{RATIO},$expectedAgreement{$feature}{RATIO},$kappa{$feature});
    }
    $expectedAgreement{OVERALL}{TIME} += $expectedAgreement{$feature}{TIME};
  }
  $timeTotal = $agreeTotal + $dontTotal;
  my $num = $agreeTotal - $expectedAgreement{OVERALL}{TIME};
  my $den = $timeTotal - $expectedAgreement{OVERALL}{TIME};
  my $totalKappa = ($den != 0) ? $num/$den : 0;
  if ($agreeTotal + $dontTotal > 0) {
    printf("|  TOTAL  |  %6.1f  |  %6.1f  |  %6.3f  |  %6.3f  |  %6.3f  |\n",$agreeTotal,$dontTotal,
           $agreeTotal/$timeTotal,$expectedAgreement{OVERALL}{TIME}/$timeTotal,$totalKappa);
  }
}

sub validateLabels {
  # print out hash tree of the marks w/ all the same time boundaries
  $" = "\/";
  foreach my $baseFile (keys %transcriptions) {
    foreach my $feature (keys %{$transcriptions{$baseFile}}) {
      foreach my $transcriber (keys %{$transcriptions{$baseFile}{$feature}}) {
        foreach my $start (sort {$a<=>$b} keys %{$transcriptions{$baseFile}{$feature}{$transcriber}}) {
          my $stop = $transcriptions{$baseFile}{$feature}{$transcriber}{$start}{STOP};
          die if $stop eq "";
        }
      }
    }
  }
}

sub computeAgreement {
  # count up the number of matching labels across the transcribers as well as the total number
  # of labels.  this allows for a 50% match if e.g. ll marked something 'APP/CLO' and xc marked only 'APP'
  my ($baseFile,$feature1,$feature2,$transcriber,$start) = @_;
  my ($comparisons,$matches) = (0,0);
  my $stop = $transcriptions{$baseFile}{$feature1}{$transcriber}{$start}{STOP};
  foreach my $label (keys %{$transcriptions{$baseFile}{$feature1}{$transcriber}{$start}{LABELS}}) {
    foreach my $otherTranscriber (keys %{$transcriptions{$baseFile}{$feature1}}) {
      next if ($otherTranscriber eq $transcriber);
      foreach my $otherLabel (keys %{$transcriptions{$baseFile}{$feature2}{$otherTranscriber}{$start}{LABELS}}) {
        $comparisons += 1;
        $matches += 1 if ($label eq $otherLabel and $label ne "NULL");
      }
    }
  }
  my $agreementRatio = ($comparisons != 0) ? $matches/$comparisons : 0;
  my $timeInAgreement = ($stop - $start) * $agreementRatio;
  my $timeNotInAgreement = ($stop - $start) * (1 - $agreementRatio);
  return ($timeInAgreement,$timeNotInAgreement);
}

sub storeLabelTime {
  my ($baseFile,$feature,$transcriber,$start) = @_;
  foreach my $label (keys %{$transcriptions{$baseFile}{$feature}{$transcriber}{$start}{LABELS}}) {
    my $stop = $transcriptions{$baseFile}{$feature}{$transcriber}{$start}{STOP};
    if (exists $labelTime{$feature}{$label}{$transcriber}) {
      $labelTime{$feature}{$label}{$transcriber} += ($stop - $start);
    } else {
      $labelTime{$feature}{$label}{$transcriber} = ($stop - $start);
    }
  }
}

sub storeAgreementMatrix {
  my ($baseFile,$feature,$start) = @_;
  my $numLabels = 0;
  my %labels;

  # calculate the total number of labels so we can adjust for multiple labels
  my $stop = $transcriptions{$baseFile}{$feature}{$transcribers[0]}{$start}{STOP};
  my $stop2 = $transcriptions{$baseFile}{$feature}{$transcribers[1]}{$start}{STOP};
  die "label lengths not the same for $baseFile $feature $start $transcribers[0]=$stop $transcribers[1]=$stop2\n" if (!$stop or !$stop2 or $stop != $stop2);
  foreach my $transcriber (@transcribers) {
    my @labels = keys %{$transcriptions{$baseFile}{$feature}{$transcriber}{$start}{LABELS}};
    $numLabels += $#labels + 1;
  }
  my $timePerCell = ($stop - $start)/$numLabels;

  # actually store information in the agreement/confusion matrix
  foreach my $label0 (keys %{$transcriptions{$baseFile}{$feature}{$transcribers[0]}{$start}{LABELS}}) {
    foreach my $label1 (keys %{$transcriptions{$baseFile}{$feature}{$transcribers[1]}{$start}{LABELS}}) {
      if (exists $agreementMatrix{$feature}{$label0}{$label1}) {
        $agreementMatrix{$feature}{$label0}{$label1} += $timePerCell;
      }
      else {
        $agreementMatrix{$feature}{$label0}{$label1} = $timePerCell;
      }
    }
  }
}

sub printAgreementMatrix {

  my @features;


  print "\n\n$transcribers[0] (rows) by $transcribers[1] (cols)\n\n";
  foreach my $feature (sort keys %agreementMatrix) {
    @features = ('VOW','APP','FLAP','FRIC','CLO','SIL') if ($feature =~ /^dg[12]$/);
    @features = ('A+VO','ASP','VL','VOI','IRR','STOP') if ($feature eq "glo");
    @features = ('+','-') if ($feature =~ /^(nas|rd)$/);
    @features = ('LAB','L-D','DEN','ALV','P-A','VEL','GLO','RHO','LAT','NONE','SIL') if ($feature eq "pl1");
    @features = ('L-D','DEN','ALV','P-A','VEL','GLO','RHO','LAT','NONE','SIL') if ($feature eq "pl2");
    @features = sort keys %{$validLabels{'vow'}} if ($feature eq "vow");

    # print row/column headings
    print "---+++++$feature\n";
    print "|       |";
    foreach my $label1 (@features) {
      print "   $label1   |";
    }
    print "\n";

    # calculate column totals
    foreach my $label0 (@features) {
      $agreementMatrix{$feature}{$label0}{TOTAL} = 0 if (not exists $agreementMatrix{$feature}{$label0}{TOTAL});
      foreach my $label1 (@features) {
        if (exists $agreementMatrix{$feature}{$label0}{$label1}) {
          $agreementMatrix{$feature}{$label0}{TOTAL} += $agreementMatrix{$feature}{$label0}{$label1};
        }
      }
    }

    # print confusion matrix normalized by columns
    foreach my $label0 (@features) {
      print "|  $label0  |";
      foreach my $label1 (@features) {
        if (exists $agreementMatrix{$feature}{$label0}{$label1}) {
          if ($optUnnormalize) {
            printf("  %.03f  |",$agreementMatrix{$feature}{$label0}{$label1});
          } else {
            printf("  %.03f  |",$agreementMatrix{$feature}{$label0}{$label1}/$agreementMatrix{$feature}{$label0}{TOTAL});
          }
        } else {
          print "    0    |";
        }
      }
      print "\n";
    }
    print "\n";
  }
}

sub storeObservedFeatureAgreement {
  my ($feature,$timeInAgreement,$timeNotInAgreement) = @_;
  if (not exists $observedAgreement{$feature}{TOTAL}) {
    $observedAgreement{$feature}{TIMEDONTAGREE} = 0;
    $observedAgreement{$feature}{TIMEAGREE} = 0;
    $observedAgreement{$feature}{TOTAL} = 0;
  }
  $observedAgreement{$feature}{TIMEAGREE} += $timeInAgreement;
  $observedAgreement{$feature}{TIMEDONTAGREE} += $timeNotInAgreement;
  $observedAgreement{$feature}{TOTAL} += $timeInAgreement + $timeNotInAgreement;
}

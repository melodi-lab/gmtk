#!/usr/bin/env perl

$strFile = $ARGV[0];
$defFile = $ARGV[1];

open FILE, "$defFile";
while (<FILE>) {
  s/(^\s+|\s+$)//g;
  ($f1,$f2,$f3) = split /\s+/;
  next if $f1 ne "#define";
  $defined{$f2} = $f3;
}
close FILE;

open FILE, "$strFile";
while ($line = <FILE>) {
  foreach $key (keys %defined) {
    $line =~ s/$key/$defined{$key}/g;
  }
  if ($line =~ /\s+([A-Z_]{2,})\s+/) {
    if (not exists $warnings{$1}) {
      warn "WARNING: $1 might need to be defined\n";
      $warnings{$1} = 1;
    }
  }
  print $line;
}
close FILE;

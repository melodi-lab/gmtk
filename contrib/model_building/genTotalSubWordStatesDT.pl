#!/usr/bin/env perl
use File::Basename;

## @file

$exe = basename $0;
$usage = "$exe <dict.sil> <word_map> <phone.state>\n";
die $usage if $#ARGV != 2;

$dict = $ARGV[0];
$wordMap = $ARGV[1];
$phoneStates = $ARGV[2];

%word2phones = readFile2HashOfArrays($dict,":");
%index2word = readFile2Hash($wordMap);
%phone2numStates = readFile2Hash($phoneStates);

%word2index = reverse %index2word;
delete $word2index{"THIS_IS_A_BUG"} if exists $word2index{"THIS_IS_A_BUG"};
%index2word = reverse %word2index;

@indices = sort {$a<=>$b} keys %index2word;
$numDTs = $#indices + 1;

print << "END_OF_HEADER";
% generated with: '$exe $dict $wordMap $phoneStates'

1                                       % number of decision trees

0                                       % index number
word_pronunciation_2_totalSubWordStates % name
2                                       % number of parents
END_OF_HEADER
print DTString(0,$numDTs);

foreach $index (@indices) {
   $word = $index2word{$index};
   @pronunciations = @{$word2phones{$word}};
   $numPronunciations = $#pronunciations+1;
   next if $numPronunciations < 1;
   print "\t% p0 = $index => '$word'\n";
   print "\t" . DTString(1,$numPronunciations);
   foreach $i (0..$#pronunciations) {
      $phones = $pronunciations[$i];
      $numStates = 0;
      foreach $phone (split /\s+/, $phones) {
         $numStates += $phone2numStates{$phone};
      }
      print "\t\t% p1 = $i => $phones\n";
      print "\t\t-1 $numStates\n";
   }
}

sub DTString {
   my ($parent,$numOptions) = @_;
   $string = "$parent $numOptions ";
   for ($i=0;$i<$numOptions-1;$i++) {
      $string .= "$i ";
   }
   $string .= "default\n";
   return $string;
}

sub readFile2HashOfArrays {
   $filename = shift;
   $char = ($#_ >= 0) ? shift : " ";
   my %hash;
   open FILE, "<$filename";
   while (<FILE>) {
      chomp;
      s/\s+/ /g;
      ($key,$value) = split $char;
      $key =~ s/(^\s+|\s+$)//g;
      $value =~ s/(^\s+|\s+$)//g;
      if (exists $hash{$key}) {
         push @{$hash{$key}}, $value;
      } else {
         $hash{$key} = [$value];
      }
   }
   close FILE;
   return %hash;
}

sub readFile2Hash {
   $filename = shift;
   $char = ($#_ >= 0) ? shift : " ";
   my %hash;
   open FILE, "<$filename";
   while (<FILE>) {
      chomp;
      s/\s+/ /g;
      ($key,$value) = split $char;
      $key =~ s/(^\s+|\s+$)//g;
      $value =~ s/(^\s+|\s+$)//g;
      if (exists $hash{$key}) {
         die "$key already exists\n";
      } else {
         $hash{$key} = $value;
      }
   }
   close FILE;
   return %hash;
}

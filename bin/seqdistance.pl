#!/usr/bin/perl -w

my $f1 = $ARGV[0];
my $f2 = $ARGV[1];

my $name1, $seq1;
open IN1, $f1 || die("Can't open file.\n");
while(<IN1>) {
  chomp;
  if( /^>/ ) {
    s/^>//;
    s/ /_/g;
    $name1 = $_;
    $seq1 = "";
  }
  else {
    $seq1 .= $_;
  }
}
close IN1;

my $name2, $seq2;
open IN2, $f2 || die("Can't open file.\n");
while(<IN2>) {
  chomp;
  if( /^>/ ) {
    s/^>//;
    s/ /_/g;
    $name2 = $_;
    $seq2 = "";
  }
  else {
    $seq2 .= $_;
  }
}
close IN2;

#get distance
my @seq1 = split($seq1, "");
my @seq2 = split($seq2, "");

my $dist = 0;
my $gap = 1;
my $match = 0;
my $mismatch = 1;
for( my $i=0; $i<$#seq1+1; $i++ ) {
   if( $seq1[ $i ] eq "-" || $seq2[ $i ] eq "-" ) {
     $dist += $gap;
   }
   elsif( $seq1[ $i ] ne $seq2[ $i ] ) {
     $dist += $mismatch;
   }
   else {
     $dist += $match;
   }
}

$dist /= $#seq1+1;
print $dist;


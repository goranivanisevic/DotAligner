#!/usr/bin/perl -w

use Cwd;

# get current directory
my $pwd = getcwd;

# process input file name
my @tmp = split /\//, $ARGV[1];
my $fname = $tmp[$#tmp];
my $len = length($ARGV[0]); 

# initialize probability matrix
my @p;
for( my $i=0; $i<$len; $i++ ) {
  for( my $j=0; $j<$len; $j++ ) {
    $p[ $i ][ $j ] = 0;
  }
}

# call RNAfold -p
#system "RNAfold -p -d2 --noLP < $fasta > /dev/null";
#system "RNAfold -p2 < $fasta > /dev/null";

# fill matrix with probabilities of folding energies predicted by RNAfold
open(file_local, $ARGV[1]);
# open(file_local, $pwd . "/" . $name . "_dp.ps");
my $start = 0;
while( <file_local> ) {
  if( /%start of base pair probability data/ ) {
    $start = 1;
    next;
  }
  $start = 0 if /showpage/;
  next unless $start;

  @l = split " ";
  # fetch small bug in RNAfold
  next if $l[0]>$len || $l[0]<1 || $l[1]>$len || $l[1]<1 || $l[2]==0;

  # RNAfold returns the square roots of the base pair probabilities sqrt{P(i,j)}
  $p[ $l[0]-1 ][ $l[1]-1 ] = $p[ $l[1]-1 ][ $l[0]-1 ] = $l[2]*$l[2] if $l[3]=~/ubox/;
}
close file_local;

# write output
print ">" . $fname . "\n" . $ARGV[0] . "\n";
for( my $i=0; $i<$len; $i++ ) {
  for( my $j=0; $j<$len; $j++ ) {
    print $p[ $i ][ $j ] . " ";
  }
  print "\n";
}
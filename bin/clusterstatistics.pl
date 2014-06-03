#!/usr/bin/perl -w
# Copyright (C) 2013 Stefan E Seemann <seemann@rth.dk>
#
# clusterstatistics.pl
# get sensitivity (SN) and specificity (SP) of a set of clusters 
# whereas the True Positives of one input family is the cluster with the largest amount of this family
#
# Example of 3 clusters (C1,C2,C3) and 3 input families (FA,FB,FC)
# Input:
# C1 FA
# C1 FA
# C1 FB
# C2 FA
# C2 FC
# C3 FB
# C3 FC
# Number of occurance of input families in each cluster:
# FC 0 1 1
# FB 1 0 1
# FA 2 1 0
# Cluster(s) where the input families are True Positives (largest occurance)
# FC 1 2
# FB 0 2
# FA 0
# Output:
# FC:	TP = 2	FP = 2	TN = 8	FN = 2	SP = 0.8000	SN = 0.5000
# FB:	TP = 2	FP = 3	TN = 7	FN = 2	SP = 0.7000	SN = 0.5000
# FA:	TP = 2	FP = 1	TN = 3	FN = 1	SP = 0.7500	SN = 0.6667
# SP = 0.7500
# SN = 0.5556

use strict;

sub help {
  print <<EOT;
  Usage: clusterstatistics.pl FILE1 FILE2
  Get sensitivity (SN) and specificity (SP) of a set of clusters 
  whereas the True Positives of one input family is the cluster with the largest amount of this family

  FILE1: list of clustered pairs - cluster family
  FILE2: list of non-clustered pairs - Number_items family
EOT
exit;
}

&help() if( $ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "-?" );


my (%h, %c, @l);
my $i = 0;

#read input
#1st column: cluster id
#2nd column: member id
my @input;
open IN, $ARGV[0] || die("Die!\n");
while( <IN> ) {
  chomp;
  push @input, $_;

  @l = split " ";
  if( ! defined $c{ $l[0] } ) {
    $c{ $l[0] } = $i;
    $i++;
  }
}
close IN;
 
my $nr = keys %c;

foreach( @input ) {
  @l = split " ";
  if( ! defined $h{ $l[1] } ) {
    for(my $k=0; $k<=$nr; $k++){
       ${$h{ $l[1] }}[ $k ] = 0;
    }
  }

  ${$h{ $l[1] }}[ $c{ $l[0] } ]++;
}

#read number of items without cluster for each family
my %noc;
open IN, $ARGV[1] || die("Die!\n");
while( <IN> ) {
  @l = split " ";
  if( ! defined $h{ $l[1] } ) {
    for(my $k=0; $k<=$nr; $k++){
       ${$h{ $l[1] }}[ $k ] = 0;
    }
  }
  ${$h{ $l[1] }}[ $nr ] = $l[0];
}
close IN;

#fill empty holes
#my $hl = keys %h;
#foreach my $key (keys %h) {
#  for(my $k=0; $k<$hl; $k++){
#    if( ! defined ${$h{ $key }}[ $k ] ){ ${$h{ $key }}[ $k ] = 0 };
#  }
#}
#map { print $_; map { print " " . $_} @{$h{$_}}; print "\n" } keys %h;

#calculate statistics
my (@tp, @fp, @tn, @fn, @sp, @sn);
my ($asp, $asn);
$i = 0;
foreach my $key (keys %h) {
  #find dominant cluster
  my $max = 0;
  my @maxi;
  for( my $k=0; $k<$#{$h{$key}}; $k++ ) {
    if( ${$h{$key}}[ $k ] == $max ) { push @maxi, $k; }
    elsif( ${$h{$key}}[ $k ] > $max ) { $max = ${$h{$key}}[ $k ]; @maxi = ( $k ); }
  }
  #print $key; map { print " " . $_ } @maxi; print "\n";

  #TP
  $tp[ $i ] = 0;
  foreach my $m ( @maxi ) {
    $tp[ $i ] += ${$h{$key}}[ $m ];
  }
  
  #FP
  $fp[ $i ] = 0;
  foreach my $m ( @maxi ) {
    foreach my $keys ( keys %h ) {
      if( $keys ne $key ) {	
        $fp[ $i ] += ${$h{$keys}}[ $m ];
      }
    }
  }
 
  #TN
  $tn[ $i ] = 0;
  foreach my $m ( @maxi ) {
    foreach my $keys ( keys %h ) {
      if( $keys ne $key ) {	
        for( my $k=0; $k<$#{$h{$key}}+1; $k++ ) {
          if( $k != $m ) {
            $tn[ $i ] += ${$h{$keys}}[ $k ];
          }
        }
      }
    }
  }
 
  #FN
  $fn[ $i ] = 0;
  foreach my $m ( @maxi ) {
    for( my $k=0; $k<$#{$h{$key}}+1; $k++ ) {
      if( $k != $m ) {	
        $fn[ $i ] += ${$h{$key}}[ $k ];
      }
    }
  }
 
  #Specificity SP
  my $denominator = $tn[ $i ] + $fp[ $i ];
  $sp[ $i ] = ( $denominator != 0 ) ? $tn[ $i ] / $denominator : 0;
  #avg SP
  $asp += $sp[ $i ]; 

  #Sensitivity SN
  $denominator = $tp[ $i ] + $fn[ $i ]; 
  $sn[ $i ] = ( $denominator != 0 ) ? $tp[ $i ] / $denominator : 0;
  #avg SN
  $asn += $sn[ $i ];

  $i++;
}

#average SP
$asp /= $i;

#average SN
$asn /= $i;

#output
my @names = keys %h;
for( my $k=0; $k<$i; $k++ ) {
  printf("%s:\tTP = %i\tFP = %i\tTN = %i\tFN = %i\tSP = %6.4f\tSN = %6.4f\n", $names[$k], $tp[$k], $fp[$k], $tn[$k], $fn[$k], $sp[$k], $sn[$k]);
}
printf("SP = %6.4f\nSN = %6.4f\n", $asp, $asn);


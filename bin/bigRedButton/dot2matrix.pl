#!/usr/bin/perl -w


use strict;
use warnings;


my %opts;


# This program takes the RNAfold (-m1) or RNAplfold (-m2) dotplot output in PS
# format to process and report the base pair probabilities in a matrix format,
# which can be used as input for RNAbound

if(scalar(@ARGV)==1){
    open(FIRSTDP,"$ARGV[0]") or die "Error in opening the input file \"$ARGV[0]\" \n";
    my @firstDP = <FIRSTDP>;
    my $seq;my $len;
    my @bpp;
  for(my $i=0; $i<scalar(@firstDP); $i++)
  {
    # find the sequence length
    if($firstDP[$i]=~/\/sequence\s/)
    {
      my $k=++$i;
      for(; $k<scalar(@firstDP); $k++)
        {
         if($firstDP[$k]=~/def/){last;}
         else{ $seq.=$firstDP[$k]; }
        }
    $seq=~s/[^\w]//g;
    $len=length $seq;

    # initialize a matrix to store the bp values to show in txt file
     for(my $i=1;$i<=$len;$i++)
     {
       for(my $j=1;$j<=$len;$j++)
       {
          $bpp[$i][$j]=0;
       }
     }
    }

    elsif($firstDP[$i]=~/^\d+.*ubox/){
                 # save the bpp values 
                 my @tt=split(/\s+/,$firstDP[$i]);
                 $bpp[$tt[0]][$tt[1]]=$tt[2]**2;
        }
  } # end of for loop 

    print $len,"\n";
    for(my $i=1;$i<=$len;$i++)
    {
      for(my $j=1;$j<=$len;$j++)
      {
        print $bpp[$i][$j]," ";
      }
      print "\n";
    }
}
else
{
 print "\n\tUsage: perl <postscript file>  >out.txt\n\n";
}

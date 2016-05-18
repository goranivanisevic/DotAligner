#!/usr/bin/perl -w
# -*-CPerl-*-
# Last changed Time-stamp: <2005-11-07 12:08:38 ivo>
# $Id: colorrna.pl,v 1.1 2005/11/07 12:42:27 ivo Exp $
# colorize a secondary structure plot with reliability annotation
# from positional entropy
use diagnostics;
use strict;
use Getopt::Std;
$main::VERSION = 1.0;
$Getopt::Std::STANDARD_HELP_VERSION=1;

our $opt_p;
getopts('p');

sub HELP_MESSAGE {
  print STDERR "\nusage: $0 FOO_ss.ps FOO_dp.ps > FOO_css.ps\n";
  print STDERR "For more details run\n\tperldoc -F $0\n";
}

HELP_MESSAGE() unless $#ARGV >0;
my $macro_seen= 0;
my %mfe = ();        # hash of mfe pairs
my @ss_ps = ('',''); # head and tail of the ss.ps file

my $aln=readClustal();
my $consStruc = swallow_ss_ps(length($aln->[0]->{seq}));    # read ss plot

print $ss_ps[0];     # print head
print "/cmark { % i cmark   draw circle around base i\n";
print "   newpath 1 sub coor exch get aload pop\n";
print "   fsize 2 div 0 360 arc stroke\n";
print "} bind def\n";

print $ss_ps[1];
if (!$macro_seen) {
  print <<_E_O_F_
/hsb {
    dup 0.3 mul 1 exch sub sethsbcolor
} bind def
/colorpair { % i j hue sat colorpair
   % draw basepair i,j in color
   % 1 index 0.00 ne {
   gsave
   newpath
   hsb
   fsize setlinewidth
   1 sub coor exch get aload pop moveto
   1 sub coor exch get aload pop lineto
   stroke
   grestore
   % } if
} bind def
_E_O_F_
}

my %hsb = getColors($aln,$consStruc);

foreach my $k (keys %hsb) {
  my ($i,$j) = split($;,$k);
  print "$i $j ", $hsb{$k}->[0], " ", $hsb{$k}->[1]," colorpair\n";
}

print $ss_ps[2];  # print main part

print "% Start Annotations\n";
foreach my $k (keys %hsb) {
  my ($i,$j) = split($;,$k);
  print "$i cmark\n" if $hsb{$k}->[2];
  print "$j cmark\n" if $hsb{$k}->[3];
}
print "\n% End Annotations\n";

print $ss_ps[3];  # print showpage etc

sub swallow_ss_ps {
  my ($l) = @_;

  # read the secondary structure plot
  my $length=0;
  my $tail=0;
  my @ss=();
  for my $i (1..$l) {
    push @ss,'.';
  }
  my $pairsFlag=0;
  my $sequenceFlag=0;

  open(ALIRNA,$ARGV[0]);
  while (<ALIRNA>) {
    $pairsFlag=1 if (/\/pairs/);
    if ($pairsFlag and /\[(\d+) (\d+)\]/) {
      $ss[$1-1]='(';
      $ss[$2-1]=')';
    }
    $pairsFlag=0 if ($pairsFlag and /def/);

    $macro_seen=1 if /colorpair/;
    $length ++ if /^\/coor/ .. /^\] def/;
    if (/^\/pairs/ .. /^\] def/) {
      $mfe{$1,$2}=1 if /(\d+)\s+(\d+)/;
    }
    next if /\d gmark$/;
    $tail++ if /^end/;
    $tail++ if /^drawpairs/;
    s/^drawpairs/% drawpairs/;
    $tail++ if /^showpage/;
    $ss_ps[$tail] .= $_;
  }
  close(ALIRNA);

  my $consStruc=join('',@ss);
  return $consStruc;
}


######################################################################
#
# readClustal(filehandle)
#
# Reads Clustal W formatted alignment file and returns it in list of
# hash references with keys "name" and "seq" for the name and the sequence,
# resp.
#
######################################################################

sub readClustal{
  #  my $fh=shift;
  my @out=();
  my (%order, $order, %alignments);

  open(ALIGN,$ARGV[1]);
  while (<ALIGN>) {
    next if ( /^\s+$/ );
    my ($seqname, $aln_line) = ('', '');
    if ( /^\s*(\S+)\s*\/\s*(\d+)-(\d+)\s+(\S+)\s*$/ ) {
      # clustal 1.4 format
      ($seqname,$aln_line) = ("$1/$2-$3",$4);
    } elsif ( /^(\S+)\s+([A-Za-z\-]+)\s*\d*$/ ) {
	  ($seqname,$aln_line) = ($1,$2);
	} else {
	  next;
  }
	  if ( !exists $order{$seqname} ) {
	    $order{$seqname} = $order++;
	  }
    $alignments{$seqname} .= $aln_line;
  }
  close(ALIGN);

  foreach my $name ( sort { $order{$a} <=> $order{$b} } keys %alignments ) {
    if ( $name =~ /(\S+):(\d+)-(\d+)/ ) {
      (my $sname,my $start, my $end) = ($1,$2,$3);
    } else {
      (my $sname, my $start) = ($name,1);
      my $str  = $alignments{$name};
      $str =~ s/[^A-Za-z]//g;
      my $end = length($str);
    }
    my $seq=$alignments{$name};
    push @out, {name=>$name,seq=>$seq};
  }
  return [@out];
}


sub getColors{

  my @aln=@{$_[0]};
  my $ss=$_[1];

  my $length=length($aln[0]->{seq});

  my $white="1 1 1";
  my $black="0 0 0" ;
  my $grey="1 1 1";

  my $red="0.0 1";
  my $ocre="0.16 1";
  my $green="0.32 1";
  my $turq="0.48 1";
  my $blue="0.65 1";
  my $violet="0.81 1";

  my $red1="0.0 0.6";
  my $ocre1="0.16 0.6";
  my $green1="0.32 0.6";
  my $turq1="0.48 0.6";
  my $blue1="0.65 0.6";
  my $violet1="0.81 0.6";

  my $red2="0.0 0.2";
  my $ocre2="0.16 0.2";
  my $green2="0.32 0.2";
  my $turq2="0.48 0.2";
  my $blue2="0.65 0.2";
  my $violet2="0.81 0.2";

  my @colorMatrix=([$red,$red1,$red2],
		   [$ocre,$ocre1,$ocre2],
		   [$green,$green1,$green2],
		   [$turq,$turq1,$turq2],
		   [$blue,$blue1,$blue2],
		   [$violet,$violet1,$violet2]);

  my @pairs=getPairs(\@aln,$ss);
  my %hsb;

  foreach my $pair (@pairs) {

    my $pairing = (@{$pair->{pairing}}>6) ? 6 : @{$pair->{pairing}};
    my $nonpairing = (@{$pair->{nonpairing}} > 2) ? 2 : @{$pair->{nonpairing}};
    

    my $color = $colorMatrix[$pairing-1][$nonpairing]; 
    my @F = split " ", $color;
    $hsb{$pair->{open}+1,$pair->{close}+1} = [$F[0], $F[1], $pair->{left}, $pair->{right}];
  }

  return %hsb;
}

######################################################################
#
# getPairs(\@aln alnref, $ss string)
#
# Evalutates the pairing of an alignment according to a given
# consensus secondary structure
#
# Returns list of all base pairs which is a hash with the following
# keys:
#
#  open ... column in the alignment of the first base in the pair "("
#  close ... column in the alignment of the second base in the pair ")"
#  all ... list of all basepairs in the different sequences in the alignment
#  pairing ... list of all different pairing basepairs
#  nonpairing ... list of all incompatible basepairs
#
######################################################################

sub getPairs{

  my @inputAln=@{$_[0]};
  my $ss=$_[1];

  # return nothing if there are no pairs
  if (!($ss=~tr/(/(/)) {
    return ();
  }

  my @aln=();
  foreach my $row (@inputAln) {
    my $seq=$row->{seq};
    $seq=uc($seq);
    $seq=~s/T/U/g;
    my @tmp=split(//,$seq);
    push @aln,\@tmp;
  }
  my @ss=split(//,$ss);

  my @pairs=();
  my @stack=();

  foreach my $column (0..$#ss) {

    my $currChar=$ss[$column];

    if ($currChar eq '(') {
      push @stack,$column;
    }

    if ($currChar eq ')') {
      my $openedCol=pop @stack;
      push @pairs,{open=>$openedCol,close=>$column};
    }
  }

  @pairs=sort {$a->{open} <=> $b->{open}} @pairs;

  foreach my $i (0..$#pairs) {
    #print "$i: $pairs[$i]->{open} - $pairs[$i]->{close}\n";

    my @all=();
    my @pairing=();
    my @nonpairing=();

    for my $j (0..$#aln) {
      my $currPair=$aln[$j][$pairs[$i]->{open}].$aln[$j][$pairs[$i]->{close}];
      push @all,$currPair;
    }

    for my $pair (@all) {
      if (($pair eq 'AU') or
	  ($pair eq 'UA') or
	  ($pair eq 'GC') or
	  ($pair eq 'CG') or
	  ($pair eq 'UG') or
	  ($pair eq 'GU')) {

	push @pairing, $pair;
      } elsif ($pair eq "--") {
	# do nothing
      } else {
	push @nonpairing,$pair;
      }
    }

    undef my %saw;
    my @uniquePairing = grep(!$saw{$_}++, @pairing);

    $pairs[$i]->{all}=[@all];
    $pairs[$i]->{pairing}=[@uniquePairing];
    $pairs[$i]->{nonpairing}=[@nonpairing];

    my $bl = "";
    my $br = "";
    my $compl = 0;
    my $compr = 0;
    foreach( @{$pairs[$i]->{pairing}} ) {
	my @b = split //;
	if( $b[0] ne $bl ) {
		$compl++;
	}
	if( $b[1] ne $br ) {
		$compr++;
	}
	$bl = $b[0];
	$br = $b[1];
    }
    $pairs[$i]->{left} = ($compl>1) ? 1 : 0;
    $pairs[$i]->{right} = ($compr>1) ? 1 : 0;
  }

  return @pairs;

}

=head1 NAME

colorrna - colorize an alirna.ps file

=head1 SYNOPSIS

   colorrna.pl alirna.ps alidot.ps > colorRNA.ps

=head1 DESCRIPTION

colorrna reads an RNA secondary structure plot and a dot plot
containing pair probabilities and covariance annotation as, produced
by C<RNAalifold -p>, and writes a new secondary structure plot with
color annotated seqence annotation to stdout.  The color annotation
is taken directly from the color dot plot file.

=head1 AUTHOR

Ivo L. Hofacker <ivo@tbi.univie.ac.at>
Stefan Washietl <wash@tbi.univie.ac.at>

=cut

#  End of file


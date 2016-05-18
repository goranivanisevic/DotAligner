#!/usr/bin/perl -w
# Copyright (C) 2014 Stefan E Seemann <seemann@rth.dk>

&usage() unless $#ARGV==0;

print "CLUSTAL W (1.83) multiple sequence alignment\n\n\n";

open IN, $ARGV[0] || die("Can not open the file!\n");

while(<IN>) {
	chomp $_;
	if(/^>/) {
		$aln{$k} = $v if defined $k;
		$v="";

		/>(.*)/;
		$k=$1;
		$k=~s/\s+/_/g;
	}
	else  {
		$_=~s/\./-/g;
		$v=$v.$_;
	}
}

close IN;

$aln{$k} = $v if defined $k;

map{ print "$_\t$aln{$_}\n" } keys %aln;
print "\n";

sub usage
{
   die("usage: " .
       "fa2aln.pl <file>\n");
}


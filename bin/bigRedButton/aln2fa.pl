#!/usr/bin/perl -w
# Copyright (C) 2014 Stefan E Seemann <seemann@rth.dk>

&usage() if defined $ARGV[0] && $ARGV[0] eq "--help";

while(<>) {
	chomp $_;
	if(/^CLUSTAL/ || /^$/ || /^\s+/) {
		next;
	}
	else {
		@l = split " ";
		$fa{$l[0]} = (defined $l[0] && defined $fa{$l[0]}) ? $fa{$l[0]}.$l[1] : $l[1];
	}
}

foreach $k (keys %fa) {
	print ">$k\n";
	for($i=0; $i<length($fa{$k}); $i+=60) {
		print substr($fa{$k},$i,60)."\n";
	}
}

sub usage
{
   die("usage: " .
       "aln2fa.pl <file>\n" .
       "       cat <file> | aln2fa.pl\n" . 
       "       aln2fa.pl --help\n");
}


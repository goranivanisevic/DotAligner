#!/usr/bin/perl -w
# Copyright (C) 2011 Stefan S Seemann <seemann@rth.dk>
#
# uniqCol1maxCol2.pl
# Return the unique list of the key column together with its
# maximal value in another column
use strict;

sub help {
    print <<EOT;
Usage: $0 FILE KEYCOLUMN[,KEYCOLUMN2,...] VALUECOLUMN
Returns unique list of the KEYCOLUMN together with its
maximal value in VALUECOLUMN

EOT
    exit;
}

# get options
help() if @ARGV!=3;

my @c = split /,/, $ARGV[1]; 
my (%h, %l);

open IN, $ARGV[0] || die("Cannot open $ARGV[0].\n");
while(<IN>)
{
	chomp;
	my @l = split " ";
	my $key;
	foreach my $k (@c) {
		$key .= $l[$k-1] . ":";
	}
	my $val = $l[$ARGV[2]-1];

	if(defined $h{$key})
	{
		if($val>$h{$key})
		{
			$h{$key}=$val;
			$l{$key}=$_;
		}
	}
	else
	{
		$h{$key}=$val;
		$l{$key}=$_;
	}
}

map{print "$l{$_}\n"}keys %h;

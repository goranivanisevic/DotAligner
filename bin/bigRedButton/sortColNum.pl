#!/usr/bin/perl -w
# Copyright (C) 2014 Stefan E Seemann <seemann@rth.dk>

$file=$ARGV[0];
$column=($ARGV[1]) ? $ARGV[1]-1 : 0;
$"="\t";

open IN, $file || die "GOODBYE\n";
while( <IN> ) {
	my @l = split " ";
	push @a, \@l;
}
close IN;

@sa = sort { $a->[$column] <=> $b->[$column] } @a;
foreach( @sa ) {
	print "@$_\n";
}

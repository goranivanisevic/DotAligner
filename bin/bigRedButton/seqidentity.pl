#!/usr/bin/perl -w
# Copyright (C) 2014 Stefan E Seemann <seemann@rth.dk>

&usage() unless $#ARGV==1;

open IN, $ARGV[0] || die("Can not open the file!\n");

$j=-1;
while(<IN>) {
	chomp $_;
	if(/^>/) {
		$j++;
		/>(.*)/;
		$name[$j] = $1;
	}
	else  {
		$_ =~ s/\./-/g;
		@tmp = split "";
		push @{$seq[$j]}, @tmp;
	}
}

close IN;

$sum=0;
$nr=0;
for( $tj=0;$tj<=$j;$tj++ ) {
    for( $ttj=$tj+1;$ttj<=$j;$ttj++ ) {
	#compare $seq[$tj] with $seq[$ttj]
	$e=0;
	$l=$#{$seq[$tj]}+1;
	$updatedl = $l;
	for( $i=0;$i<$l;$i++ ) {
	    if( ${$seq[$tj]}[$i] eq ${$seq[$ttj]}[$i] ) {
		if( ${$seq[$tj]}[$i] eq "-" ) {
			$updatedl--;
		}
		else {
			$e++ 	
		}
	    }
	}
	$id=$e/$updatedl;
	print "$name[$tj]\t$name[$ttj]\t$id\t$e\t$updatedl\n" if $id>=$ARGV[1]/100;
	$sum+=$id;
	$nr++;
    }
}

$id=$sum/$nr*100;
#print "average pairwise identity: $id\n";
print "$id\n";

sub usage
{
   die("usage: seqidentity.pl <fasta-file> <min-identity-in-percent>\n" .
       "The identities between all sequence pairs in the fasta-file are calculated.");
}


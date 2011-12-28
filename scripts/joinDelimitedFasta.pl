#!/usr/bin/perl

use strict;
use warnings;

# make sure delimiter is properly escaped, e.g. for "+" we need "\+".
my $delimiter = $ARGV[0];
my $infileA = $ARGV[1];
my $infileB = $ARGV[2];

open(INA, $infileA) or die "Can't read $infileA\n";
open(INB, $infileB) or die "Can't read $infileB\n";

my $fastaHeader = "";
	
while(<INA>)
	{
	chomp $a;
	my $b = <INB>;
	
	print $a . $delimiter . $b . "\n";

	}
close INA;
close INB;
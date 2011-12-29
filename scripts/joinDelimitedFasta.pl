#!/usr/bin/perl

use strict;
use warnings;

# make sure delimiter is properly escaped, e.g. for "+" we need "\+".
my $delimiter = $ARGV[0];
my $infileA = $ARGV[1];
my $infileB = $ARGV[2];
my $rc2 = $ARGV[3];

open(INA, $infileA) or die "Can't read $infileA\n";
open(INB, $infileB) or die "Can't read $infileB\n";

my $fastaHeader = "";

#sub revcompl {
#    return ((my $rev = reverse $_) =~ tr/ACGTacgt/TGCAtgca/); 
#}
	
while(<INA>)
	{
	my $a = $_;
	my $b = <INB>;
	chomp $a;
	chomp $b;
	if($rc2)
	    {
	    # $b = revcompl($b);
	    $b = reverse $b;
	    $b =~ tr/ACGTacgt/TGCAtgca/;
	    }
	
	print $a . $delimiter . $b . "\n";

	}
close INA;
close INB;
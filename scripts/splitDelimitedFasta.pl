#!/usr/bin/perl

use strict;
use warnings;

# make sure delimiter is properly escaped, e.g. for "+" we need "\+".
my $delimiter = $ARGV[0];
my $infile = $ARGV[1];

my $infilename = `basename $infile`;
chomp $infilename;

open(IN, $infile) or die "Can't read $infile\n";
open(A,">$infilename.a") or die "Can't write $infilename.a\n";
open(B,">$infilename.b") or die "Can't write $infilename.b\n";

my $fastaHeader = "";
my $outfile = *A;
	
while(<IN>)
	{
	chomp;
	if(/^>/)
		{
		$outfile = *A;
		$fastaHeader = $_;
		print A "$fastaHeader\n";
		print B "$fastaHeader\n";
		}
	else
		{
		if ( /(.*)$delimiter(.*)/)
		    {
		        print A "$1\n";
		        print B "$2\n";
		        $outfile = *B;
	        }
		else
		    {
		        print $outfile "$_\n";
	        }
		}
	}
close IN;
close A;
close B;
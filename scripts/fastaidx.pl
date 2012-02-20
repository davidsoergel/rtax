#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Temp;

use FindBin;
use lib "$FindBin::Bin";
use FastaIndex;

my ($fastaFileName, $idRegex) = @ARGV;

my $indexA = FastaIndex->new();    # '-filename' => "A.idx", '-write_flag' => 1 );
print STDERR "Indexing $fastaFileName...\n";
$indexA->make_index( $fastaFileName, $idRegex, $fastaFileName );

while(<STDIN>)
{
    chomp;
    print $indexA->fetch($_);
}
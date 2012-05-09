#!/usr/bin/env perl

# read a uclust result file; for each query sequence, collect the prokMSAids of the hits.

use strict;
use warnings;

my %clusters   = ();
my %worstPcids = ();

parseUclust( $ARGV[0] );

for my $key ( sort ( keys %clusters ) ) {
    my $worstPcid = $worstPcids{$key};

    #print STDERR "$key -> $clusters{$key}\n";
    print "$key\t$worstPcid\t" . ( join "\t", @{ $clusters{$key} } ) . "\n";
}

sub parseUclust {
    my ($ucFileName) = @_;
    open( UC, $ucFileName ) || die("Could not open $ucFileName");

    while (<UC>) {

        if (/^\s*#/) { next; }

        my ( $type, $cluster, $size, $percentid, $strand, $querystart, $targetstart, $alignment, $querylabel, $targetlabel ) = split /\t/;
        chomp $querylabel;
        chomp $targetlabel;

        if ( $type eq "S" ) {

            #print STDERR "S $querylabel\n";
            $clusters{$querylabel}    = [$querylabel];
            $worstPcids{$querylabel} = 100.0;
        }
        elsif ( $type eq "H" ) {

            #print STDERR "H $targetlabel $querylabel\n";
            push @{ $clusters{$targetlabel} }, $querylabel;

            if ( $percentid < $worstPcids{$targetlabel} ) {
                $worstPcids{$targetlabel} = $percentid;
            }

        }

        # ignore other types
    }
    close UC;
}

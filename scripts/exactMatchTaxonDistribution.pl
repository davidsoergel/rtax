#!/usr/bin/perl

# read the prokMSAids of the hits, as produced by exactMatchIds.
# Then grab the taxonomy info for those prokMSAids and make some classification info from the set.

use strict;
use warnings;

#use lib '.';
use FindBin;
use lib "$FindBin::Bin";

use LevelPrediction;

loadTaxonomy( $ARGV[0] );
processHits( $ARGV[1] );

my %taxonomies = ();

# my %taxonomyPrior = ();

sub loadTaxonomy {
    my ($taxFileName) = @_;
    print STDERR "loading taxonomy...\n";
    open( TAX, $taxFileName ) or die "Cannot open input taxonomy file $taxFileName: $!";
    my $lines = 0;
    while ( my $line = <TAX> ) {
		chomp $line;
        my ( $prokMSAid, @taxonomy ) = split( /\t/, $line );

	    if ( @taxonomy && ($taxonomy[0] eq "" || $taxonomy[0] =~ /^(\d+\.?\d*|\.\d+)$/ )) {
	        # value is numeric or empty, must be the pcid width of the target cluster
	        my $taxPcid = shift @taxonomy;
	        #print STDERR "Ignoring target cluster width: $taxPcid\n";
	    }

        $taxonomies{$prokMSAid} = \@taxonomy;

        my $taxString = "";
        for my $taxElem (@taxonomy) {
            $taxString .= "$taxElem; ";

            #			$taxonomyPrior{$taxString}++;
        }
        $lines++;
    }
    close TAX;
    print STDERR "...done loading $lines taxonomy lines\n";

    #	print STDERR "normalizing prior...\n";
    #	for my $taxKey (keys %taxonomyPrior)
    #	{
    #		$taxonomyPrior{$taxKey} /= $lines;
    #	}
    #	$taxonomyPrior{""} = 1;
    #    print STDERR "...done normalizing prior\n";
}

sub printLine {
    my ( $label, $bestPcid, @ids ) = @_;

    my $levelPredictions = makeTaxonomyPrediction(@ids);
    print "$label\t$bestPcid";
    for my $levelPrediction (@$levelPredictions) {
        print "\t" . $levelPrediction->toString();
    }

    #if(scalar(@$levelPredictions) == 0) {$unclassified++}
    #if(scalar(@$levelPredictions) == 0) { print "\tUNCLASSIFIED"; }

    print "\n";

    return ( scalar(@$levelPredictions) == 0 );
}

sub processHits {
    my ($hitsFileName) = @_;
    open( HITS, $hitsFileName ) || die("Could not open $hitsFileName");

    my $hit          = 0;
    my $noprimer     = 0;
    my $nohit        = 0;
    my $toomanyhits  = 0;
    my $unclassified = 0;

    while (<HITS>) {
        chomp;
        my ( $label, $bestPcid, @ids ) = split /\t/;

        if ( ( !@ids ) || ( $ids[0] eq "" ) )    #  an empty id list
		{
			print "$label\t\tNOHIT\n";
            $nohit++;
		}
		#elsif($ids[0] eq "NOHIT")     # this does not happen; use an empty id list instead up to this point
        #{
        #    print "$label\tNOHIT\n";
        #    $nohit++;
        #}
        elsif ( ( $ids[0] eq "NOPRIMER" ) ) {
            print "$label\t\tNOPRIMER\n";
            $noprimer++;
        }
        elsif ( ( $ids[0] eq "TOOMANYHITS" ) ) {
            print "$label\t\tTOOMANYHITS\n";
            $toomanyhits++;
        }
        else {
            $unclassified += printLine( $label, $bestPcid, @ids );
            $hit++;
        }
    }

    my $samples = $hit + $nohit + $toomanyhits + $noprimer;

    print STDERR "$samples items\n";
    if ( $samples != 0 ) {
        print STDERR "$noprimer (" . sprintf( "%.1f",     ( 100 * $noprimer / $samples ) ) . "%) had no primer match\n";
        print STDERR "$nohit (" . sprintf( "%.1f",        ( 100 * $nohit / $samples ) ) . "%) had no hits\n";
        print STDERR "$toomanyhits (" . sprintf( "%.1f",        ( 100 * $toomanyhits / $samples ) ) . "%) had too many hits\n";
        print STDERR "$unclassified (" . sprintf( "%.1f", ( 100 * $unclassified / $samples ) ) . "%) had hits but no classification\n";
        print STDERR "" . ( $hit - $unclassified ) . " (" . sprintf( "%.1f",          ( 100 * ( $hit - $unclassified ) / $samples ) ) . "%) were classified\n";
    }

    # classifications per level?  That depends on the stringency filter, which is downstream

    close HITS;
}

#sub normalizeByPrior {
#	my ($taxString, $taxonCounts) = @_;
#
#	print STDERR "Normalizing: \n";
#	while( my ($k, $v) = each %$taxonCounts) {
#        print STDERR "$k = $v ; ";
#    }
#	print STDERR "\n";
#
#	my $normalizer = $taxonomyPrior{$taxString};
#	if(!defined $normalizer)
#		{
#			print STDERR ("No prior: $taxString\n");
#		}
#
#	for my $taxon (keys %$taxonCounts)
#	{
#		$taxonCounts->{$taxon} = ($taxonCounts->{$taxon} / $taxonomyPrior{$taxString . $taxon . "; "}) * $normalizer;
#	}
#	print STDERR "       Done: \n";
#	while( my ($k, $v) = each %$taxonCounts) {
#        print STDERR "$k = $v ; ";
#    }
#	print STDERR "\n";
#}

sub makeTaxonomyPrediction {
    my (@ids) = @_;
	
    my @levelPredictions = ();

# TOOMANYHITS is now added in ucFilterBest.pl, so it should be caught at line 99 above

#	if(scalar(@ids) == 1000) {
#		# when we did the uclust at exactMatchOneReadSet:50, we used " --allhits --maxaccepts 1000 ".
#		# after that we filtered for the best pcid set (using ucFilterBest.pl).
#		# if 1000 hits remain, then the real set of best-pcid hits is larger, and we missed some.
#		# In that case we should bail because the set of 1000 hits we do have may not be representative.
#		# I think this is the reason why matching the E517F primer only (17nt reads) produced predictions, and at different levels to boot.
#		# That also depends on the classification thresholds.
#			
#		my $levelPrediction = LevelPrediction->new();	
#		$levelPrediction->label("TOOMANYHITS");
#		push @levelPredictions, $levelPrediction;
#		return \@levelPredictions;
#		}

    my @taxonomyVectors = map { $taxonomies{$_} } @ids;
#	my @taxonomyClusterSizes = map { $taxonomyWorstPcids{$_} } @ids;

    my $globalTotal = @taxonomyVectors;

    my $predictedTaxa    = "";

    my $globalUnknowns = 0;    # at all levels, incrementing as we go

    for my $level ( 0 .. 10 ) {

        my $levelPrediction = LevelPrediction->new();

        # assert the remaining taxonomyVectors are equal at higher levels

        my %taxonCounts   = ();
        my $localUnknowns = 0;
        my $localTotal    = @taxonomyVectors;

        # count up how often each label is seen descending from this node

        for my $vec (@taxonomyVectors) {
            my $taxon = $vec->[$level];

            # "Unknown" should occur only towards the leaves; an unspecified intermediate level followed by a populated level is a "skip".
            # Here, "skip" is counted as just another taxon.
            if   ( !defined $taxon || $taxon eq "" || $taxon =~ /unknown/i ) { $localUnknowns++; }
            else                                                             { $taxonCounts{$taxon}++; }
        }

        if ( $localUnknowns == $localTotal ) { last; }

        # this normalization makes no sense, because we don't want a uniform prior either
        # e.g., one Archaeal hit among dozens of Bacterial hits will win, because there are so few Archaea in GreenGenes to begin with
        # normalizeByPrior( $predictedTaxa, \%taxonCounts );

        # get the best label and check for ties

        $levelPrediction->numChildren( scalar( keys %taxonCounts ) );

        my @taxaInOrder = sort { $taxonCounts{$b} <=> $taxonCounts{$a} } keys %taxonCounts;
        my $bestTaxon = $taxaInOrder[0];

        # print STDERR "$bestLabel\n";
        $levelPrediction->label($bestTaxon);
        $predictedTaxa .= "$bestTaxon; ";

        my $bestCount = $taxonCounts{$bestTaxon};
        $levelPrediction->count($bestCount);

        my $secondBestTaxon = $taxaInOrder[1];
        if ( defined $secondBestTaxon ) {
            my $secondBestCount = $taxonCounts{$secondBestTaxon};

            if ( $levelPrediction->count() < 2 * $secondBestCount ) {
                # Declare a tie if the winning taxon doesn't have at least twice as many votes as the runner-up.  
				# just consider this level unknown and bail
                last;
            }
        }

        # compute various ratios of the prediction vs. alternatives

        $levelPrediction->localMinProportion( $bestCount / $localTotal );
        $levelPrediction->localMaxProportion( ( $bestCount + $localUnknowns ) / $localTotal );

        my $globalUnknowns += $localUnknowns;
        $levelPrediction->globalMinProportion( $bestCount / $globalTotal );
        $levelPrediction->globalMaxProportion( ( $bestCount + $globalUnknowns ) / $globalTotal );

        # what if all of the "unknown" matches should have been the same?  Then an "alternate" classification might have won
        $levelPrediction->alternateLocalProportion( $localUnknowns / $localTotal );
        $levelPrediction->alternateGlobalProportion( $globalUnknowns / $globalTotal )
            ;    # note we already added the local unknowns to the global unknowns

        # it's possible that a completely different path has a higher global probability than this one,
        # but we'd never know because we pick the max _per level_ and never explore the other paths.

        # decide whether to bother continuing
        # if ( $bestLocalMinProportion < 0.5 ) { last; }
        # for now, print all the best labels until everything is unknown or there is a tie; sort out the confidence later

        push @levelPredictions, $levelPrediction;

        # remove any non-matching paths from consideration at the next level

        my @newTaxonomyVectors = ();
        for my $vec (@taxonomyVectors) {
            my $taxon = $vec->[$level];
            if ( defined $taxon && $taxon eq $bestTaxon ) { push @newTaxonomyVectors, $vec; }
        }
        @taxonomyVectors = @newTaxonomyVectors;

    }
    return \@levelPredictions;

}

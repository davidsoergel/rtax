#!/usr/bin/perl

my $countThreshold               = $ARGV[0];
my $localMinProportionThreshold  = $ARGV[1];
my $globalMinProportionThreshold = $ARGV[2];

while (<STDIN>) {
    chomp;
    my ( $id, $bestPcid, @taxonomyEntries ) = split /\t/;

    if ( $taxonomyEntries[0] eq "NOPRIMER" ) {
        print "$id\t\tNOPRIMER\n";
    }
    elsif ( $taxonomyEntries[0] eq "NOHIT" ) {
        print "$id\t\tNOHIT\n";
    }
	elsif ( $taxonomyEntries[0] =~ /TOOMANYHITS/ ) {
	    print "$id\t\tTOOMANYHITS\n";
	}
    else {
        my @newTaxonomyEntries = ();
        for my $entry (@taxonomyEntries) {

            my (
                $label,               $numChildren,              $count,
                $localMinProportion,  $localMaxProportion,       $globalMinProportion,
                $globalMaxProportion, $alternateLocalProportion, $alternateGlobalProportion
            ) = split /;/, $entry;

            if (   $count >= $countThreshold
                && $localMinProportion >= $localMinProportionThreshold
                && $globalMinProportion >= $globalMinProportionThreshold )
            {
                push @newTaxonomyEntries, $entry;
            }
            else { last; }
        }
        #if ( @newTaxonomyEntries == 0 ) {
        #    push @newTaxonomyEntries, "Root";
        #}

		# if the trailing entries of the filtered taxonomy string are "SKIP", strip them
		while(substr(@newTaxonomyEntries[$#newTaxonomyEntries],0,4) eq "SKIP") { pop @newTaxonomyEntries; }

        print "$id\t$bestPcid\t" . ( join "\t", @newTaxonomyEntries ) . "\n";
    }

}

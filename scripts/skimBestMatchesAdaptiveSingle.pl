#!/usr/bin/perl

# read a uclust result file from STDIN; for each query sequence, collect the prokMSAids of those hits that are within $slop of the best.

use strict;
use warnings;
use Bio::Index::Fasta;

my $usearch                             = shift @ARGV;
my $slop                                = shift @ARGV;
my $minMaxAccepts                       = shift @ARGV;
my $maxMaxAccepts                       = shift @ARGV;
my $maxSinglePercentDifferenceThreshold = shift @ARGV;
my $databaseFile                        = shift @ARGV;
my $readAFileAll                        = shift @ARGV;

print STDERR "usearch = $usearch\n";
print STDERR "slop = $slop\n";
print STDERR "minMaxAccepts = $minMaxAccepts\n";
print STDERR "maxMaxAccepts = $maxMaxAccepts\n";
print STDERR "maxSinglePercentDifferenceThreshold = $maxSinglePercentDifferenceThreshold\n";
print STDERR "databaseFile = $databaseFile\n";
print STDERR "readAFileAll = $readAFileAll\n";

my $indexA;

# allow variable nember of reads?
# my @readFiles = @ARGV;

sub printLine {
    my ( $label, $idsByIdentity ) = @_;

    my @ids        = ();
    my @pcids      = sort { $b <=> $a } keys %$idsByIdentity;
    my $bestPcid   = $pcids[0];
    my $acceptPcid = $bestPcid - $slop;
    for my $pcid (@pcids) {
        if ( $pcid < $acceptPcid ) { last; }
        push @ids, @{ $idsByIdentity->{$pcid} };
    }

    my $idString = ( join "\t", @ids );
    print "$label\t$bestPcid\t" . $idString . "\n";
}

sub collectIds {

    # assume that the records for a given label are all grouped together, so we can collect the hits for each label and process them in turn
    # need to jump through some hoops to save the first line of the next block (since we can't rewind a pipe with seek()).

    my $fh        = shift;
    my $firstLine = shift;

    my %hits = ();

    # process the first line (cached from the previous call)
    #print STDERR "Processing firstLine: $firstLine\n";
    my ( $typeF, $clusterF, $sizeF, $percentIdF, $strandF, $queryStartF, $targetStartF, $alignmentF, $queryLabelF, $targetLabelF ) =
        split /\t/, $firstLine;
    chomp $targetLabelF;
    if ( $typeF eq "N" ) {
        #print STDERR "$queryLabelF -> N\n";

        # NOHIT represented as empty hash
        # note we still have to cache the next line, after filtering comments
        while (<$fh>) {
            if (/^\s*#/) { next; }
            if (/^\s*$/) { next; }
            my $line = $_;
            chomp $line;
            #print STDERR "$line\n";
            return ( $queryLabelF, \%hits, $line );
        }
    }
    elsif ( $typeF eq "H" ) {
        $hits{$targetLabelF} = $percentIdF;
    }

    # process the remaining lines
    while (<$fh>) {
        if (/^\s*#/) { next; }
        if (/^\s*$/) { next; }
        my $line = $_;
        chomp $line;
        #print STDERR "$line\n";

        my ( $type, $cluster, $size, $percentId, $strand, $queryStart, $targetStart, $alignment, $queryLabel, $targetLabel ) = split /\t/;
        chomp $targetLabel;
    
        if ( $queryLabel ne $queryLabelF ) {

            # advanced past the end of the current block
            print STDERR "End of block $queryLabelF -> block of " . scalar( keys %hits ) . " hits\n";
            #print STDERR "End of block, return\n";
            return ( $queryLabelF, \%hits, $line );
        }
        else {
            if ( $type eq "N" ) {
                die "impossible: $queryLabel reports no hit after it already had a hit";
            }
            elsif ( $type eq "H" ) {
                $hits{$targetLabel} = $percentId;
            }
        }
    }

    # end of the stream (the last block)
    print STDERR "End of stream $queryLabelF -> block of " . scalar( keys %hits ) . " hits\n";
    #print STDERR "End of stream, return\n";
    return ( $queryLabelF, \%hits, "" );
}

sub reconcileHitsAndPrint {
    my ( $queryLabel, $idsA, $singlePercentIdThreshold ) = @_;

    my $idsByIdentity = {};
    for my $targetLabel ( keys %$idsA ) {
        my $percentIdA = $idsA->{$targetLabel};
        if ( !defined $idsByIdentity->{$percentIdA} ) {
            $idsByIdentity->{$percentIdA} = [];
        }
        push @{ $idsByIdentity->{$percentIdA} }, $targetLabel;

    }

    if ( !%$idsByIdentity ) {
        print STDERR "$queryLabel -> no reconciled hits at $singlePercentIdThreshold%\n";

        # this is the NOHIT case, but we'll back off and try again
        return 0;
    }
    else {
        print STDERR "$queryLabel -> printing reconciled hits at $singlePercentIdThreshold%\n";
        printLine( $queryLabel, $idsByIdentity );
        return 1;
    }
}

sub doSearch {
    my $numQuerySequences                = shift;
    my $readAFile                        = shift;
    my $singlePercentDifferenceThreshold = shift;
    my $maxAccepts                       = shift;

    my $singlePercentIdThreshold = 1. - $singlePercentDifferenceThreshold;

    print STDERR
"$numQuerySequences query sequences remaining; searching $readAFile with %id $singlePercentDifferenceThreshold and maxAccepts $maxAccepts\n";

    my $nohitQueryIds      = [];
    my $tooManyHitQueryIds = [];

    print STDERR "$usearch --quiet --iddef 2 --query $readAFile --db $databaseFile --uc /dev/stdout --id $singlePercentIdThreshold --maxaccepts $maxAccepts --maxrejects 128 --nowordcountreject\n";

    # open the USEARCH streams
    open( UCA,
"$usearch --quiet --iddef 2 --query $readAFile --db $databaseFile --uc /dev/stdout --id $singlePercentIdThreshold --maxaccepts $maxAccepts --maxrejects 128 --nowordcountreject |"
    ) || die "can't fork usearch: $!";

    # Load the first non-comment line from each stream
    my $nextLineA;
    while (<UCA>) {
        if (/^\s*#/) { next; }
        if (/^\s*$/) { next; }
        $nextLineA = $_;
        last;
    }

    # read synchronized blocks from each stream
    while (1) {
        my ( $queryLabelA, $idsA );
        #print STDERR "reading next block...\n";

        # idsA is a reference to a hash from targetIds to %ids
        ( $queryLabelA, $idsA, $nextLineA ) = collectIds( *UCA, $nextLineA );

        if ( ( scalar keys %$idsA ) >= $maxAccepts ) {
            push @$tooManyHitQueryIds, $queryLabelA;
        }

        elsif ( !reconcileHitsAndPrint( $queryLabelA, $idsA, $singlePercentIdThreshold ) ) {
            push @$nohitQueryIds, $queryLabelA;
        }

        if ( !$nextLineA ) {
            #print STDERR "End of stream, close\n";
            last;
        }

    }

    close(UCA) || die "can't close usearch: $!";

    #print STDERR "Closed usearch stream.\n";

    # for the TOOMANYHITS cases, escalate maxAccepts and try again
    # Note this recurses, so no need to iterate here
    if ( scalar(@$tooManyHitQueryIds) && ( $maxAccepts * 2 <= $maxMaxAccepts ) ) {
        my $nohitQueryIdsB;

        print STDERR "Escalating maxAccepts to " . ( $maxAccepts * 2 ) . " for " . scalar(@$tooManyHitQueryIds) . " sequences.\n";

        # prepare input files with the remaining sequences
        $readAFile = extractFasta( $indexA, $tooManyHitQueryIds );

        ( $nohitQueryIdsB, $tooManyHitQueryIds ) =
            doSearch( scalar(@$tooManyHitQueryIds), $readAFile, $singlePercentDifferenceThreshold, $maxAccepts * 2 );
        if (@$nohitQueryIdsB) {
            die "A TOOMANYHITS case can't turn into a NOHITS case";
        }
    }

    # any tooManyHitQueryIds that remain had more than maxMaxAccepts hits
    return ( $nohitQueryIds, $tooManyHitQueryIds );
}

my $fastaNum;

sub extractFasta {
    my $index = shift;
    my $ids   = shift;
    my $name  = $fastaNum++;

    print STDERR "Extracting " . scalar(@$ids) . " to file $name\n";

    open( OUT, ">$name" );

    my $out = Bio::SeqIO->new( '-format' => 'Fasta', '-fh' => \*OUT );
    for my $id (@$ids) {
        my $seqobj = $index->fetch($id);
        if ( !defined $seqobj ) {
            print STDERR "Undefined: $id\n";
        }
        else {
            $out->write_seq($seqobj);
        }
    }
    close OUT;

    return $name;
}

sub main {

    my $nohitQueryIds = [];
    push @$nohitQueryIds, "ALL";
    my $tooManyHitQueryIds = [];

    my $singlePercentDifferenceThreshold = 0.005;    # this gets doubled to 1% before being used the first time

    # these are redundant between multiple runs, oh well
    $indexA = Bio::Index::Fasta->new( '-filename' => "A.idx", '-write_flag' => 1 );
    $indexA->make_index($readAFileAll);

    my $readAFile = $readAFileAll;

    while ( @$nohitQueryIds > 0 ) {

        # double the allowed %diff on every round
        $singlePercentDifferenceThreshold *= 2;
        if ( $singlePercentDifferenceThreshold > $maxSinglePercentDifferenceThreshold ) {
            last;
        }

        if ( $nohitQueryIds->[0] ne "ALL" ) {

            # prepare input files with the remaining sequences
            $readAFile = extractFasta( $indexA, $nohitQueryIds );
        }

        # within doSearch we escalate maxAccepts, so if a queryLabel is still marked tooManyHits at this point,
        # that means that there are more than maxMaxAccepts hits for this threshold--
        # so there's really no point in testing this query again at an even lower threshold
        my $tooManyHitQueryIdsThisRound;
        ( $nohitQueryIds, $tooManyHitQueryIdsThisRound ) =
            doSearch( scalar(@$nohitQueryIds), $readAFile, $singlePercentDifferenceThreshold, $minMaxAccepts );

        print STDERR "Finished round at threshold $singlePercentDifferenceThreshold; "
            . scalar(@$nohitQueryIds)
            . " NOHIT, "
            . scalar(@$tooManyHitQueryIdsThisRound)
            . " TOOMANYHIT.\n";
        push @$tooManyHitQueryIds, @$tooManyHitQueryIdsThisRound;
    }

    print STDERR scalar(@$nohitQueryIds) . " query sequences remaining with NOHIT\n";
    for my $queryLabel (@$nohitQueryIds) {
        print "$queryLabel\t\tNOHIT\n";
    }

    print STDERR scalar(@$tooManyHitQueryIds) . " query sequences remaining with TOOMANYHITS\n";
    for my $queryLabel (@$tooManyHitQueryIds) {
        print "$queryLabel\t\tTOOMANYHITS\n";
    }

}

main();

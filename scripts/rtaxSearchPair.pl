#!/usr/bin/perl

# RTAX: Rapid and accurate taxonomic classification of short paired-end
#       sequence reads from the 16S ribosomal RNA gene.
#
# David A. W. Soergel 1*, Rob Knight 2, and Steven E. Brenner 1
#
# 1 Department of Plant and Microbial Biology,
#   University of California, Berkeley
# 2 Howard Hughes Medical Institute and Department of Chemistry
#   and Biochemistry, University of Colorado at Boulder
# * Corresponding author: soergel@berkeley.edu
#
# http://www.davidsoergel.com/rtax
#
# Version 0.9303  (January 9, 2012)
#
# For usage instructions: just run the script with no arguments
#
#
# Copyright (c) 2011 Regents of the University of California
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the University of California, Berkeley nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

use strict;
use warnings;
use Getopt::Long;
use File::Temp;

use FindBin;
use lib "$FindBin::Bin";
use FastaIndex;

use vars qw (
    $usearch
    $slop
    $minMaxAccepts
    $maxMaxAccepts
    $maxPairPercentDifferenceThreshold
    $databaseFile
    $readAFileAll
    $readBFileAll
    $idRegex
);

# each file provides reads from a different primer extracted from full-length sequence, in the same order

sub init {
    Getopt::Long::Configure("bundling");
    GetOptions(
        "usearch=s"                           => \$usearch,
        "slop=s"                              => \$slop,
        "minMaxAccepts=s"                     => \$minMaxAccepts,
        "maxMaxAccepts=s"                     => \$maxMaxAccepts,
        "maxPairPercentDifferenceThreshold=s" => \$maxPairPercentDifferenceThreshold,
        "databaseFile=s"                      => \$databaseFile,
        "idRegex=s"                           => \$idRegex,
        "queryA=s"                            => \$readAFileAll,
        "queryB=s"                            => \$readBFileAll
    );

    if (   !defined $databaseFile
        || !defined $readAFileAll
        || !defined $readBFileAll )
    {
        print STDERR "Missing required argument!\n\n";
        usage();
    }

    if ( !defined $usearch ) {
        $usearch = `whereis usearch`;
        chomp $usearch;
        if ( !defined $usearch || $usearch eq "" ) {
            print STDERR "Could not find usearch.\n\n";
            usage();
        }
    }

    if ( !defined $slop )          { $slop          = 0.005; }
    if ( !defined $minMaxAccepts ) { $minMaxAccepts = 1000; }
    if ( !defined $maxMaxAccepts ) { $maxMaxAccepts = 16000; }
    if ( !defined $maxPairPercentDifferenceThreshold ) {
        $maxPairPercentDifferenceThreshold = 0.2;
    }
    if ( !defined $idRegex || $idRegex eq "") { $idRegex = "(\\S+)"; }
    print STDERR "Header Regex: $idRegex\n";
}

sub usage {
    print STDERR << "EOF";

RTAX: Rapid and accurate taxonomic classification of short paired-end
      sequence reads from the 16S ribosomal RNA gene.
      
David A. W. Soergel (1*), Rob Knight (2), and Steven E. Brenner (1)
1 Department of Plant and Microbial Biology, University of California, Berkeley 
2 Howard Hughes Medical Institute and Department of Chemistry and Biochemistry,
  University of Colorado at Boulder
* Corresponding author: soergel\@berkeley.edu

http://www.davidsoergel.com/rtax

Version 0.9  (September 7, 2011)

OPTIONS:

    --usearch <path>    location of usearch program
                        (defaults to result of "whereis usearch")
     
    --slop <number>     %id difference from maximum to accept
                        (default 0.005, i.e. accept hits to a 200nt query that
                        has 1 more mismatch than the best match, or 2 extra
                        mismatches out of 400nt.)
     
    --minMaxAccepts <number>
                        The initial maxaccepts value to pass to USEARCH,
                        controlling the maximum number of hits that will be
                        returned in the first iteration.  If this number of hits
                        is reached, the search is repeated with a doubled
                        maxaccepts value.  Thus, a smaller value causes the
                        first iteration to run faster, but increases the
                        probability of needing more iterations.
                        (Default 1000)
     
    --maxMaxAccepts <number>
                        The highest escalated maxaccepts value to allow.  The
                        maxaccepts value is doubled on every iteration, so this
                        value is effectively rounded down to a power of two 
                        factor of minMaxAccepts.  USEARCH runs with large 
                        maxaccepts values are very slow, so this controls when 
                        RTAX gives up on a query because there are too many
                        hits (i.e., a query that is too short and/or too highly
                        conserved to be informative).  (Default 16000) 

    --maxPairPercentDifferenceThreshold <number>
                        The largest percent difference to allow between a query
                        (considering the read pair jointly) and a database hit.
                        (Default 0.2, i.e. reference sequences of less than 80%
                        identity with the query sequence are never matched)

    --databaseFile <file>
                        A FASTA file containing a reference database.
                        
    --queryA <file>     A FASTA file containing one set of query reads.
    
    --queryB <file>     A FASTA file containing the other set of query reads.
    
    --idRegex <regex>   A regular expression for extracting IDs from fasta
                        headers.  The first capture group (aka $1) will be
                        used as the sequence ID; these should be unique.
                        This is useful when the ID needed to match mate
                        pairs is embedded in the header in some way.
                        "^>" is automatically prepended to the regex
                        provided here.  Take care that double escaping may
                        be required.  Default: "(\\S+)" (the first field)
    
    Note that the two query files must provide mate-paired reads with exactly
    matching identifiers (though they need not be in the same order).  Any ids
    present in one file but not the other are silently dropped.  Please contact
    us for help accomodating alternate naming schemes for the paired reads
    (e.g., "SOMEID.a" paired with "SOMEID.b", and so forth).
    
    Various indexes and derived FASTA files will be created in a temporary
    directory and are cleaned up afterwards.
    
    The output is tab-delimited text, one line per query pair, in the form
    
        <query ID>  <%id>   <list of reference IDs>
    
    where the second column gives the %id between the query and the best match,
    and the reference IDs provided are all those matches within "slop" of this
    best value.
    
    This output is suitable for piping directly into rtaxVote for taxonomy
    assignment.

EOF

    exit;
}

my $indexA;
my $indexB;

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

    # print STDERR "$label\t$bestPcid\t" . $idString . "\n";
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
    
    ($queryLabelF) = $queryLabelF =~ /$idRegex/;
    #$queryLabelF =~ s/\s.*//;    # primary ID is only the portion before whitespace

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

        my ( $type, $cluster, $size, $percentId, $strand, $queryStart, $targetStart, $alignment, $queryLabel, $targetLabel ) =
            split /\t/;
        chomp $targetLabel;
        
        ($queryLabelF) = $queryLabelF =~ /$idRegex/;
        #$queryLabel =~ s/\s.*//;    # primary ID is only the portion before whitespace

        if ( $queryLabel ne $queryLabelF ) {

            # advanced past the end of the current block
            #print STDERR "End of block $queryLabelF -> block of " . scalar( keys %hits ) . " hits\n";

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
    #print STDERR "End of stream $queryLabelF -> block of " . scalar( keys %hits ) . " hits\n";

    #print STDERR "End of stream, return\n";
    return ( $queryLabelF, \%hits, "" );
}

sub reconcileHitsAndPrint {
    my ( $queryLabel, $idsA, $idsB, $pairPercentIdThreshold ) = @_;

    my $idsByIdentity = {};
    for my $targetLabel ( keys %$idsA ) {
        my $percentIdA = $idsA->{$targetLabel};
        my $percentIdB = $idsB->{$targetLabel};
        if ($percentIdB) {

            # both reads hit the same target
            # compute the average percent id
            my $pairPercentId = ( $percentIdA + $percentIdB ) / 2.0;
            if ( $pairPercentId >= $pairPercentIdThreshold ) {
                if ( !defined $idsByIdentity->{$pairPercentId} ) {
                    $idsByIdentity->{$pairPercentId} = [];
                }
                push @{ $idsByIdentity->{$pairPercentId} }, $targetLabel;
            }
        }
    }

    if ( !%$idsByIdentity ) {

        #print STDERR "$queryLabel -> no reconciled hits at $pairPercentIdThreshold%\n";

        # this is the NOHIT case, but we'll back off and try again
        return 0;
    }
    else {

        #print STDERR "$queryLabel -> printing reconciled hits at $pairPercentIdThreshold%\n";
        printLine( $queryLabel, $idsByIdentity );
        return 1;
    }

}

sub doPairSearch {
    my $numQuerySequences              = shift;
    my $readAFile                      = shift;
    my $readBFile                      = shift;
    my $pairPercentDifferenceThreshold = shift;
    my $maxAccepts                     = shift;

    my $pairPercentIdThreshold = 1. - $pairPercentDifferenceThreshold;

    # because we're going to average the two %ids, we have to search each read with twice the %diff for the pair
    my $singlePercentDifferenceThreshold = $pairPercentDifferenceThreshold * 2;
    my $singlePercentIdThreshold         = 1. - $singlePercentDifferenceThreshold;

    print STDERR "doPairSearch $readAFile, $readBFile: $numQuerySequences query sequences remaining\n     searching with pair %id "
        . $pairPercentDifferenceThreshold
        . " and maxAccepts "
        . $maxAccepts . "\n";

    my $nohitQueryIds      = [];
    my $tooManyHitQueryIds = [];

# open the USEARCH streams
#    print STDERR
#"$usearch --quiet --iddef 2 --query $readAFile --db $databaseFile --uc /dev/stdout --id $singlePercentIdThreshold --maxaccepts $maxAccepts --maxrejects 128 --nowordcountreject\n";

#    open( UCA,
#"$usearch --quiet --iddef 2 --query $readAFile --db $databaseFile --uc /dev/stdout --id $singlePercentIdThreshold --maxaccepts $maxAccepts --maxrejects 128 --nowordcountreject |"
#    ) || die "can't fork usearch: $!";

#    print STDERR
#"$usearch --quiet --iddef 2 --query $readBFile --db $databaseFile --uc /dev/stdout --id $singlePercentIdThreshold --maxaccepts $maxAccepts --maxrejects 128 --nowordcountreject\n";

#    open( UCB,
#"$usearch --quiet --iddef 2 --query $readBFile --db $databaseFile --uc /dev/stdout --id $singlePercentIdThreshold --maxaccepts $maxAccepts --maxrejects 128 --nowordcountreject |"
#    ) || die "can't fork usearch: $!";

    my $dir = File::Temp->newdir();
    if ( system( 'mknod', "$dir/a", 'p' ) && system( 'mkfifo', "$dir/a" ) ) { die "mk{nod,fifo} $dir/a failed"; }
    if ( system( 'mknod', "$dir/b", 'p' ) && system( 'mkfifo', "$dir/b" ) ) { die "mk{nod,fifo} $dir/b failed"; }
    if ( !fork() ) {
        my $cmd =
"$usearch --quiet --iddef 2 --query $readAFile --db $databaseFile --uc $dir/a --id $singlePercentIdThreshold --maxaccepts $maxAccepts --maxrejects 128 --nowordcountreject";
        print STDERR $cmd . "\n";
        exec $cmd || die "can't fork usearch: $!";
    }

    if ( !fork() ) {
        my $cmd =
"$usearch --quiet --iddef 2 --query $readBFile --db $databaseFile --uc $dir/b --id $singlePercentIdThreshold --maxaccepts $maxAccepts --maxrejects 128 --nowordcountreject";
        print STDERR $cmd . "\n";
        exec $cmd || die "can't fork usearch: $!";
    }

    open( UCA, "$dir/a" );
    open( UCB, "$dir/b" );

    # Load the first non-comment line from each stream
    my $nextLineA;
    while (<UCA>) {
        if (/^\s*#/) { next; }
        if (/^\s*$/) { next; }
        $nextLineA = $_;
        last;
    }

    my $nextLineB;
    while (<UCB>) {
        if (/^\s*#/) { next; }
        if (/^\s*$/) { next; }
        $nextLineB = $_;
        last;
    }

    # read synchronized blocks from each stream
    while (1) {
        my ( $queryLabelA, $idsA, $queryLabelB, $idsB );

        # idsA is a reference to a hash from targetIds to %ids
        ( $queryLabelA, $idsA, $nextLineA ) = collectIds( *UCA, $nextLineA );
        ( $queryLabelB, $idsB, $nextLineB ) = collectIds( *UCB, $nextLineB );

        if ( !( $queryLabelA eq $queryLabelB ) ) {
            die "Usearch results desynchronized: $queryLabelA neq $queryLabelB";
        }

        my $numHitsA = ( scalar keys %$idsA );
        my $numHitsB = ( scalar keys %$idsB );

        # if either read got NOHITS, then it's definitely NOHITS for the pair.  This trumps TOOMANYHITS for the other read.
        # don't bother trying to reconcile hits in this case

        if ( ( $numHitsA == 0 ) || ( $numHitsB == 0 ) ) {
            push @$nohitQueryIds, $queryLabelA;
        }

        # if both reads got TOOMANYHITS, then it's definitely TOOMANYHITS for the pair.

        elsif ( ( $numHitsA >= $maxAccepts ) && ( $numHitsB >= $maxAccepts ) ) {
            push @$tooManyHitQueryIds, $queryLabelA;
        }

        # if neither read got TOOMANYHITS, then we're good to go

        elsif ( ( $numHitsA < $maxAccepts ) && ( $numHitsB < $maxAccepts ) ) {
            if ( !reconcileHitsAndPrint( $queryLabelA, $idsA, $idsB, $pairPercentIdThreshold ) ) {
                push @$nohitQueryIds, $queryLabelA;
            }

        }
        else {

            # if only one read got TOOMANYHITS, escalate if possible
            if ( $maxAccepts < $maxMaxAccepts ) {
                push @$tooManyHitQueryIds, $queryLabelA;
            }

            # but if only one read got TOOMANYHITS and we're already at maxMaxAccepts, fall back to the info provided by the other read.
            elsif ( $numHitsA < $maxAccepts ) {
                if ( !reconcileHitsAndPrint( $queryLabelA, $idsA, $idsA, $pairPercentIdThreshold ) ) {
                    push @$nohitQueryIds, $queryLabelA;
                }
            }
            elsif ( $numHitsB < $maxAccepts ) {
                if ( !reconcileHitsAndPrint( $queryLabelA, $idsB, $idsB, $pairPercentIdThreshold ) ) {
                    push @$nohitQueryIds, $queryLabelA;
                }
            }
            else { die "impossible"; }
        }
        if ( !$nextLineA || !$nextLineB ) {
            if ( !( !$nextLineA && !$nextLineB ) ) {
                die "Usearch results desynchronized at end:\nA: $nextLineA\nB: $nextLineB";
            }
            last;
        }

    }

    close(UCA) || die "can't close usearch: $!";
    close(UCB) || die "can't close usearch: $!";

    # for the TOOMANYHITS cases, escalate maxAccepts and try again
    # Note this recurses, so no need to iterate here
    if ( scalar(@$tooManyHitQueryIds) && ( $maxAccepts * 2 <= $maxMaxAccepts ) ) {
        my $nohitQueryIdsB;

        print STDERR "doPairSearch $readAFile, $readBFile: Escalating maxAccepts to "
            . ( $maxAccepts * 2 ) . " for "
            . scalar(@$tooManyHitQueryIds)
            . " sequences.\n";

        # prepare input files with the remaining sequences
        my $readAFileEsc = extractFasta( $indexA, $tooManyHitQueryIds );
        my $readBFileEsc = extractFasta( $indexB, $tooManyHitQueryIds );

        ( $nohitQueryIdsB, $tooManyHitQueryIds ) =
            doPairSearch( scalar(@$tooManyHitQueryIds), $readAFileEsc, $readBFileEsc, $pairPercentDifferenceThreshold, $maxAccepts * 2 );

        # A TOOMANYHITS case can certainly turn into a NOHITS case:
        # once one read is no longer TOOMANYHITS, it may turn out that nothing can be reconciled with the other read.

        #if (@$nohitQueryIdsB) {
        #    die(      "A TOOMANYHITS case can't turn into a NOHITS case: "
        #            . ( join ", ", @$nohitQueryIdsB )
        #            . " at pair threshold $pairPercentDifferenceThreshold and maxAccepts "
        #            . ( $maxAccepts * 2 ) );
        #  }

        push @$nohitQueryIds, @$nohitQueryIdsB;
    }

    print STDERR
        "doPairSearch $readAFile, $readBFile: Finished at pair threshold $pairPercentDifferenceThreshold and maxAccepts $maxAccepts\n";
    print STDERR "         NOHITS: " . scalar(@$nohitQueryIds) . "\n";
    print STDERR "    TOOMANYHITS: " . scalar(@$tooManyHitQueryIds) . "\n";

    # print STDERR "         NOHITS: " .      ( join ", ", @$nohitQueryIds ) . "\n";
    # print STDERR "    TOOMANYHITS: " . ( join ", ", @$tooManyHitQueryIds ) . "\n";

    # any tooManyHitQueryIds that remain had more than maxMaxAccepts hits
    return ( $nohitQueryIds, $tooManyHitQueryIds );
}

my $fastaNum;

sub extractFasta {
    my $index = shift;
    my $ids   = shift;
    my $name  = $fastaNum++;

    #print STDERR "Extracting " . scalar(@$ids) . " to file $name\n";

    open( OUT, ">$name" );

    # my $out = Bio::SeqIO->new( '-format' => 'Fasta', '-fh' => \*OUT );
    for my $id (@$ids) {
        my $seqobj = $index->fetch($id);
        if ( !defined $seqobj ) {
            print STDERR "Extracting from " . $index->fastaFileName() . ": Undefined: $id\n";
        }

        # elsif ($seqobj->primary_id() ne $id)
        #    {
        #    print STDERR "Extracting from " . $index->filename() . ": ID problem: " . ($seqobj->primary_id()) . " ne $id\n";
        #    }
        else {

            #$out->write_seq($seqobj);
            print OUT $seqobj;
        }
    }
    close OUT;

    return $name;
}

sub main {

    init();

    my $nohitQueryIds = [];
    push @$nohitQueryIds, "ALL";
    my $tooManyHitQueryIds = [];

    my $pairPercentDifferenceThreshold = 0.005;    # this gets doubled to 1% before being used the first time

    # these are redundant between multiple runs, oh well
    $indexA = FastaIndex->new();                   # '-filename' => "A.idx", '-write_flag' => 1 );
    $indexA->make_index($readAFileAll, $idRegex);

    $indexB = FastaIndex->new();                   # '-filename' => "B.idx", '-write_flag' => 1 );
    $indexB->make_index($readBFileAll, $idRegex);

    my @union        = ();
    my @intersection = ();
    my @difference   = ();
    my %count        = ();
    my @idsA         = keys %{ $indexA->db() };
    my @idsB         = keys %{ $indexB->db() };

    #print STDERR "idsA = " . (join " ", @idsA) . "\n";
    #print STDERR "idsB = " . (join " ", @idsB) . "\n";

    foreach my $element ( @idsA, @idsB ) {
        if ( ( $element =~ /[\t ]/ ) > 0 )    # $indexA->db() should return parsed IDs
        {
            print STDERR "Invalid FASTA id: $element\n";
            exit(1);
        }
        $count{$element}++;
    }
    foreach my $element ( keys %count ) {
        if ( !( $element =~ /^__/ ) ) {
            push @union, $element;
            push @{ $count{$element} > 1 ? \@intersection : \@difference }, $element;
        }
    }

    # this is where we could fall back to single-ended classification as needed

    # sequence ids mentioned in one of the two files, but not both
    foreach my $element (@difference) {
        print "$element\t\tNOPRIMER\n";
        print STDERR "$element\t\tNOPRIMER\n";
    }

    #print STDERR "intersection = " . ( join " ", @intersection ) . "\n";
    print STDERR "intersection = " . scalar(@intersection) . " sequences\n";

    my $readAFile = extractFasta( $indexA, \@intersection );
    my $readBFile = extractFasta( $indexB, \@intersection );

    while ( @$nohitQueryIds > 0 ) {

        # double the allowed %diff on every round
        $pairPercentDifferenceThreshold *= 2;
        if ( $pairPercentDifferenceThreshold > $maxPairPercentDifferenceThreshold ) {
            last;
        }

        if ( $nohitQueryIds->[0] ne "ALL" ) {

            # prepare input files with the remaining sequences
            $readAFile = extractFasta( $indexA, $nohitQueryIds );
            $readBFile = extractFasta( $indexB, $nohitQueryIds );
        }

        # within doPairSearch we escalate maxAccepts, so if a queryLabel is still marked tooManyHits at this point,
        # that means that there are more than maxMaxAccepts hits for this threshold--
        # so there's really no point in testing this query again at an even lower threshold
        my $tooManyHitQueryIdsThisRound;

        my $numRemaining = ( $nohitQueryIds->[0] eq "ALL" ) ? scalar(@intersection) : scalar(@$nohitQueryIds);

        ( $nohitQueryIds, $tooManyHitQueryIdsThisRound ) =
            doPairSearch( $numRemaining, $readAFile, $readBFile, $pairPercentDifferenceThreshold, $minMaxAccepts );

        print STDERR "MAIN: Finished round at threshold $pairPercentDifferenceThreshold; "
            . scalar(@$nohitQueryIds)
            . " NOHIT, "
            . scalar(@$tooManyHitQueryIdsThisRound)
            . " TOOMANYHIT.\n";

        push @$tooManyHitQueryIds, @$tooManyHitQueryIdsThisRound;
    }

    print STDERR "MAIN: " . scalar(@$nohitQueryIds) . " query sequences remaining with NOHIT\n";
    for my $queryLabel (@$nohitQueryIds) {
        print "$queryLabel\t\tNOHIT\n";
        print STDERR "$queryLabel\t\tNOHIT\n";
    }

    print STDERR "MAIN: " . scalar(@$tooManyHitQueryIds) . " query sequences remaining with TOOMANYHITS\n";
    for my $queryLabel (@$tooManyHitQueryIds) {
        print "$queryLabel\t\tTOOMANYHITS\n";
        print STDERR "$queryLabel\t\tTOOMANYHITS\n";
    }

}

main();

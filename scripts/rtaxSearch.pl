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
# Version 0.96  (February 19, 2012)
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

### command-line arguments get stored in globals

use vars qw (
    $usearch
    $slop
    $minMaxAccepts
    $maxMaxAccepts
    $maxPercentDifferenceThreshold
    $databaseFile
    $readAFileAll
    $readBFileAll
    $idRegex
    $idList
    $singleOK
);

# just store these as globals for now-- they're basically like command-line args
my $indexA;
my $indexB;

sub init {
    Getopt::Long::Configure("bundling");
    GetOptions(
        "usearch=s"                       => \$usearch,
        "slop=s"                          => \$slop,
        "minMaxAccepts=s"                 => \$minMaxAccepts,
        "maxMaxAccepts=s"                 => \$maxMaxAccepts,
        "maxPercentDifferenceThreshold=s" => \$maxPercentDifferenceThreshold,
        "databaseFile=s"                  => \$databaseFile,
        "idRegex=s"                       => \$idRegex,
        "queryA=s"                        => \$readAFileAll,
        "queryB=s"                        => \$readBFileAll,
        "idList=s"                        => \$idList,
        "singleOK"                        => \$singleOK
    );

    if (   !defined $databaseFile
        || !defined $readAFileAll )
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
    if ( !defined $maxPercentDifferenceThreshold ) {
        $maxPercentDifferenceThreshold = 0.2;
    }
    if ( !defined $idRegex || $idRegex eq "" ) { $idRegex = "(\\S+)"; }
    print STDERR "Header Regex: $idRegex\n";

    # these are redundant between multiple runs, oh well
    # but the DBM indexes should persist, in the same directory as the original files
    $indexA = FastaIndex->new();    # '-filename' => "A.idx", '-write_flag' => 1 );
    $indexA->make_index( $readAFileAll, $idRegex, $readAFileAll );

    $indexB = FastaIndex->new();
    if ( defined $readBFileAll ) {    # '-filename' => "B.idx", '-write_flag' => 1 );
        $indexB->make_index( $readBFileAll, $idRegex, $readAFileAll );
    }
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

Version 0.95  (February 10, 2012)

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

    --maxPercentDifferenceThreshold <number>
                        The largest percent difference to allow between a query
                        (considering the read pair jointly) and a database hit.
                        (Default 0.2, i.e. reference sequences of less than 80%
                        identity with the query sequence are never matched)

    --databaseFile <file>
                        A FASTA file containing a reference database.
                        
    --queryA <file>     A FASTA file containing one set of query reads.
    
    --queryB <file>     A FASTA file containing the other set of query reads
                        (if any).
                        
    --singleOK          Classify a sequence based on only one read.  Required
                        when queryB is absent.  Default: sequences present in
                        only one of the two input files are dropped.  When
                        enabled, paired-end sequences are classified using
                        both reads as usual, but any remaining single-ended
                        reads are also classified.
    
    --idRegex <regex>   A regular expression for extracting IDs from fasta
                        headers.  The first capture group (aka \$1) will be
                        used as the sequence ID; these should be unique.
                        This is useful when the ID needed to match mate
                        pairs is embedded in the header in some way.
                        "^>" is automatically prepended to the regex
                        provided here.  Take care that double escaping may
                        be required.  Default: "(\\S+)" (the first field)
    
    --idList <file>     A file containing a list of IDs to process (one per 
                        line).  By default all IDs found in the query files
                        are used.
    
    Note that the two query files must provide mate-paired reads with exactly
    matching identifiers (though they need not be in the same order).  Any ids
    present in one file but not the other are classified in single-ended mode.
    Alternate naming schemes for the paired reads (e.g., "SOMEID.a" paired with
    "SOMEID.b", and so forth) are handled via the --idRegex option.
    
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

sub main {

    init();

    my ( $pairedReadAFile, $pairedReadBFile, $pairedBothCount, $singleReadAFile,$singleReadACount, $singleReadBFile,$singleReadBCount ) = partitionReadFiles();

    processPairs( $pairedReadAFile, $pairedReadBFile, $pairedBothCount );
    
    if($singleOK)
        {
            processSingle($singleReadAFile, $singleReadACount);
            processSingle($singleReadBFile, $singleReadBCount);
        }
}

sub partitionReadFiles {

    my @idList = ();
    if ($idList) {
        open( IDLIST, "$idList" ) or die "Could not open $idList";
        while (<IDLIST>) { chomp; push @idList, $_; }
        close IDLIST;
    }

    my @idsA = keys %{ $indexA->startpos() };    
    my @idsB = defined $indexB->startpos() ? keys %{ $indexB->startpos() } : ();    # in the single-read case, there should still be an empty index

    my @bothAandB = ();
    my @aOnly     = ();
    my @bOnly     = ();

    # encode which of the input files contain which IDs in three bits

    my %count = ();

    foreach my $element (@idsA) {
        if ( ( $element =~ /[\t ]/ ) > 0 )    # $indexA->db() should return parsed IDs with no whitespace
        {
            die "Invalid FASTA id: $element\n";
        }
        $count{$element} += 1;
    }

    foreach my $element (@idsB) {
        if ( ( $element =~ /[\t ]/ ) > 0 )    # $indexA->db() should return parsed IDs with no whitespace
        {
            die "Invalid FASTA id: $element\n";
        }
        $count{$element} += 2;
    }

    foreach my $element (@idList) {
        if ( ( $element =~ /[\t ]/ ) > 0 ) { die "Invalid FASTA id: $element\n"; }
        $count{$element} += 4;
    }

    if ($idList) {
        foreach my $element ( keys %count ) {
            $count{$element} -= 4;
        }
    }

    foreach my $element ( keys %count ) {
        if ( !( $element =~ /^__/ ) ) {
            if    ( $count{$element} == 1 ) { push @aOnly,     $element }
            elsif ( $count{$element} == 2 ) { push @bOnly,     $element }
            elsif ( $count{$element} == 3 ) { push @bothAandB, $element }
            else                            { 
                # no problem; these are sequences not in the idList
                }
        }
    }

    print STDERR "file a = " . scalar(@aOnly) . " sequences\n";
    print STDERR "file b = " . scalar(@bOnly) . " sequences\n";
    print STDERR "both   = " . scalar(@bothAandB) . " sequences\n";

    my $pairedReadAFile = extractFasta( $indexA, \@bothAandB );
    my $pairedReadBFile = extractFasta( $indexB, \@bothAandB );

    my $singleReadAFile = extractFasta( $indexA, \@aOnly );
    my $singleReadBFile = extractFasta( $indexB, \@bOnly );

    return ( $pairedReadAFile, $pairedReadBFile, scalar(@bothAandB), $singleReadAFile, scalar(@aOnly), $singleReadBFile, scalar(@bOnly) );
}

sub processPairs {
    my ( $pairedReadAFile, $pairedReadBFile,$pairedBothCount ) = @_;
    if ($pairedBothCount  == 0) { return; }
    
    my $nohitQueryIds = [];
    push @$nohitQueryIds, "ALL";
    my $tooManyHitQueryIds = [];

    my $percentDifferenceThreshold = 0.005;    # this gets doubled to 1% before being used the first time

    while ( @$nohitQueryIds > 0 ) {

        # double the allowed %diff on every round
        $percentDifferenceThreshold *= 2;
        if ( $percentDifferenceThreshold > $maxPercentDifferenceThreshold ) {
            last;
        }

        if ( $nohitQueryIds->[0] ne "ALL" ) {

            # prepare input files with the remaining sequences
            $pairedReadAFile = extractFasta( $indexA, $nohitQueryIds );
            $pairedReadBFile = extractFasta( $indexB, $nohitQueryIds );
        }

        # within doPairSearch we escalate maxAccepts, so if a queryLabel is still marked tooManyHits at this point,
        # that means that there are more than maxMaxAccepts hits for this threshold--
        # so there's really no point in testing this query again at an even lower threshold
        my $tooManyHitQueryIdsThisRound;

        my $numRemaining = ( $nohitQueryIds->[0] eq "ALL" ) ? $pairedBothCount : scalar(@$nohitQueryIds);

        print STDERR "doPairSearch $pairedReadAFile, $pairedReadBFile: $numRemaining query sequences remaining\n     searching with pair %id "
            . $percentDifferenceThreshold
            . " and maxAccepts "
            . $minMaxAccepts . "\n";

        ( $nohitQueryIds, $tooManyHitQueryIdsThisRound ) =
            doPairSearch( $pairedReadAFile, $pairedReadBFile, $percentDifferenceThreshold, $minMaxAccepts );

        print STDERR "MAIN: Finished round at threshold $percentDifferenceThreshold; "
            . scalar(@$nohitQueryIds)
            . " NOHIT, "
            . scalar(@$tooManyHitQueryIdsThisRound)
            . " TOOMANYHIT.\n";

        push @$tooManyHitQueryIds, @$tooManyHitQueryIdsThisRound;
    }

    print STDERR "MAIN: " . scalar(@$nohitQueryIds) . " query sequences remaining with NOHIT\n";
    for my $queryLabel (@$nohitQueryIds) {
        print "$queryLabel\t\tNOHIT\n";
    }

    print STDERR "MAIN: " . scalar(@$tooManyHitQueryIds) . " query sequences remaining with TOOMANYHITS\n";
    for my $queryLabel (@$tooManyHitQueryIds) {
        print "$queryLabel\t\tTOOMANYHITS\n";
    }

}

sub processSingle {
    my ($singleReadFile, $singleReadCount) = @_;
    if ( !defined $singleReadFile || $singleReadFile eq "" || $singleReadCount == 0) { return; }

    my $nohitQueryIds = [];
    push @$nohitQueryIds, "ALL";
    my $tooManyHitQueryIds = [];

    my $singlePercentDifferenceThreshold = 0.005;    # this gets doubled to 1% before being used the first time

    while ( @$nohitQueryIds > 0 ) {

        # double the allowed %diff on every round
        $singlePercentDifferenceThreshold *= 2;
        if ( $singlePercentDifferenceThreshold > $maxPercentDifferenceThreshold ) {
            last;
        }

        if ( $nohitQueryIds->[0] ne "ALL" ) {

            # prepare input files with the remaining sequences
            $singleReadFile = extractFasta( $indexA, $nohitQueryIds );
        }

        my $numRemaining = ( $nohitQueryIds->[0] eq "ALL" ) ? $singleReadCount : scalar(@$nohitQueryIds);

        # within doSearch we escalate maxAccepts, so if a queryLabel is still marked tooManyHits at this point,
        # that means that there are more than maxMaxAccepts hits for this threshold--
        # so there's really no point in testing this query again at an even lower threshold
        my $tooManyHitQueryIdsThisRound;

        print STDERR "doSingleSearch $singleReadFile: $numRemaining query sequences remaining\n     searching with %id "
            . $singlePercentDifferenceThreshold
            . " and maxAccepts "
            . $minMaxAccepts . "\n";

        ( $nohitQueryIds, $tooManyHitQueryIdsThisRound ) =
            doSingleSearch( $singleReadFile, $singlePercentDifferenceThreshold, $minMaxAccepts );

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

sub doSingleSearch {
    my $readAFile                        = shift;
    my $singlePercentDifferenceThreshold = shift;
    my $maxAccepts                       = shift;

    my $singlePercentIdThreshold = 1. - $singlePercentDifferenceThreshold;

    my $nohitQueryIds      = [];
    my $tooManyHitQueryIds = [];

# open the USEARCH streams
#  print STDERR "$usearch --quiet --iddef 2 --query $readAFile --db $databaseFile --uc /dev/stdout --id $singlePercentIdThreshold --maxaccepts $maxAccepts --maxrejects 128 --nowordcountreject\n";
# open( UCA,
# "$usearch --quiet --iddef 2 --query $readAFile --db $databaseFile --uc /dev/stdout --id $singlePercentIdThreshold --maxaccepts $maxAccepts --maxrejects 128 --nowordcountreject |"
#    ) || die "can't fork usearch: $!";

    my $dir = File::Temp->newdir();
    if ( system( 'mknod', "$dir/a", 'p' ) && system( 'mkfifo', "$dir/a" ) ) { die "mk{nod,fifo} $dir/a failed"; }
    if ( !fork() ) {
        my $cmd =
"$usearch --quiet --iddef 2 --query $readAFile --db $databaseFile --uc $dir/a --id $singlePercentIdThreshold --maxaccepts $maxAccepts --maxrejects 128 --nowordcountreject";
        print STDERR $cmd . "\n";
        exec $cmd || die "can't fork usearch: $!";
    }

    open( UCA, "$dir/a" );

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

        elsif ( !reconcileSingleHitsAndPrint( $queryLabelA, $idsA, $singlePercentIdThreshold ) ) {
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
    
    print STDERR
        "doPairSearch $readAFile: Finished at pair threshold $singlePercentDifferenceThreshold and maxAccepts $maxAccepts\n";
    print STDERR "         NOHITS: " . scalar(@$nohitQueryIds) . "\n";
    print STDERR "    TOOMANYHITS: " . scalar(@$tooManyHitQueryIds) . "\n";

    # print STDERR "         NOHITS: " .      ( join ", ", @$nohitQueryIds ) . "\n";
    # print STDERR "    TOOMANYHITS: " . ( join ", ", @$tooManyHitQueryIds ) . "\n";

    # any tooManyHitQueryIds that remain had more than maxMaxAccepts hits
    return ( $nohitQueryIds, $tooManyHitQueryIds );
}

sub reconcileSingleHitsAndPrint {
    my ( $queryLabel, $idsA, $idsB, $singlePercentIdThreshold ) = @_;

    my $idsByIdentity = {};
    for my $targetLabel ( keys %$idsA ) {
        my $percentIdA = $idsA->{$targetLabel};
        if ( !defined $idsByIdentity->{$percentIdA} ) {
            $idsByIdentity->{$percentIdA} = [];
        }
        push @{ $idsByIdentity->{$percentIdA} }, $targetLabel;

    }

    if ( !%$idsByIdentity ) {

        #print STDERR "$queryLabel -> no reconciled hits at $singlePercentIdThreshold%\n";

        # this is the NOHIT case, but we'll back off and try again
        return 0;
    }
    else {

        #print STDERR "$queryLabel -> printing reconciled hits at $singlePercentIdThreshold%\n";
        printLine( $queryLabel, $idsByIdentity );
        return 1;
    }

}

sub doPairSearch {
    my $readAFile                      = shift;
    my $readBFile                      = shift;
    my $pairPercentDifferenceThreshold = shift;
    my $maxAccepts                     = shift;

    my $pairPercentIdThreshold = 1. - $pairPercentDifferenceThreshold;

    # because we're going to average the two %ids, we have to search each read with twice the %diff for the pair
    my $singlePercentDifferenceThreshold = $pairPercentDifferenceThreshold * 2;
    my $singlePercentIdThreshold         = 1. - $singlePercentDifferenceThreshold;

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
            if ( !reconcilePairedHitsAndPrint( $queryLabelA, $idsA, $idsB, $pairPercentIdThreshold ) ) {
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
            doPairSearch( $readAFileEsc, $readBFileEsc, $pairPercentDifferenceThreshold, $maxAccepts * 2 );

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

sub reconcilePairedHitsAndPrint {
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

##### Read a block of target IDs for one query from a USEARCH stream
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

    $queryLabelF =~ /$idRegex/;
    $queryLabelF = $1;

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

        $queryLabel =~ /$idRegex/;
        $queryLabel = $1;

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

my $fastaNum;

sub extractFasta {
    my $index = shift;
    my $ids   = shift;
    my $name  = $fastaNum++;

    # print STDERR "Extracting " . scalar(@$ids) . " to file $name\n";

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

main();

#!/usr/bin/perl

use strict;
use warnings;

my $fieldname = $ARGV[0];

while (<STDIN>) {

    # read in by blocks
    if ( $_ =~ /^BEGIN$/ ) {
        my $prokMSAid   = "NONE";
        my $prokMSAname = "NONE";
        my $organism    = "NONE";
        my $val         = "NONE";

        until ( $_ =~ /^END$/ ) {
            if (/^prokMSA_id=(.*)/) {
                $prokMSAid = $1;
            }
            elsif (/^prokMSAname=(.*)/) {
                $prokMSAname = $1;
            }
            elsif (/^organism=(.*)/) {
                $organism = $1;
            }
            elsif (/^$fieldname=(.*)/) {
                $val = $1;
            }
            $_ = <STDIN>;
        }

        my @taxArray = formatTaxonomy( $val, $prokMSAname, $organism );

        print "$prokMSAid\t" . ( join "\t", @taxArray ) . "\n";
    }

    # else ignore anything outside a BEGIN/END block

}

sub formatTaxonomy {
    my ( $semicolonSeparated, $name, $organism ) = @_;

    $semicolonSeparated =~ s/Unclassified//g;
    $semicolonSeparated =~ s/"//g;

    my @tax = split( /;\s*/, $semicolonSeparated );

    # from fillInTheUnclassifieds
    my $numberTaxDivs  = 8;
    my $numberElements = @tax;

    if ( $numberElements > $numberTaxDivs ) {
        print "ERROR: taxonomy string has more than $numberTaxDivs elements: $semicolonSeparated\n";
    }

    for ( my $i = $numberElements ; $i < $numberTaxDivs ; $i++ ) {
        $tax[$i] = "";    #"Unclassified";
    }

    my ( $sp, $str, $info ) = parseSpeciesStrainName( $name, $organism );
    $tax[$numberTaxDivs]       = $sp;
    $tax[ $numberTaxDivs + 1 ] = $str;
    $tax[ $numberTaxDivs + 2 ] = $info;

    correctSkipsAndUnknowns( $numberTaxDivs, \@tax );

    return @tax;

}

sub correctSkipsAndUnknowns {
    my ( $numberTaxDivs, $tax ) = @_;

    # "Unknown" should occur only towards the leaves; an unspecified intermediate level followed by a populated level is a "skip".
    my $maxIndex = $numberTaxDivs + 2;
    my $gotone   = 0;
    for my $i ( 0 .. $maxIndex ) {

        # reverse traverse
        my $idx = $maxIndex - $i;

        if ( !defined $tax->[$idx] ) {
            $tax->[$idx] = "";
        }

        $tax->[$idx] = trim( $tax->[$idx] );

        if ( $tax->[$idx] eq "SKIP" || $tax->[$idx] eq "UNKNOWN" ) {
            $tax->[$idx] = "";
        }

        if ( $tax->[$idx] ne "" ) {
            $gotone = 1;

            # print STDERR "$idx: $tax[$idx] unchanged\n";
        }
        elsif ($gotone) {

            #   print STDERR "$idx: $tax[$idx] = SKIP\n";
            $tax->[$idx] = "SKIP";    # a more specific classification was already seen
        }
        else {

            #  print STDERR "$idx: $tax[$idx] = UNKNOWN\n";
            $tax->[$idx] = "UNKNOWN";
        }
    }

}

sub trim {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

sub parseSpeciesStrainName {
    my ( $prokMSAname, $organism ) = @_;

    my $sp   = $organism;
    my $str  = "";
    my $info = "";

    $sp =~ s/uncultured\s*//g;
    $sp =~ s/bacterium\s*//g;

    while ( $sp =~ /^(.*?)\s*(clone|isolate.*)$/i ) {
        $sp   = $1;
        $info = "$info $2 ";
    }

    if ( $sp =~ /^(.*?)\s*(str\.|strain|subsp\.)(.*)$/i ) {
        $sp  = $1;
        $str = $3;
    }

    if ( $sp =~ /^(.*?)\s*(sp\.)\s*(.*)$/ ) {
        $sp  = "$1 sp. ";
        $str = "$str $3 ";
    }

    # Species is first two words
    if ( $sp =~ /^(.+?)\s(.+?)\s(.+)$/ ) {
        $sp  = "$1 $2 ";
        $str = "$3 $str ";
    }

    $sp   = trim($sp);
    $str  = trim($str);
    $info = trim($info);

    #   print STDERR "$organism->[ $sp | $str | $info ] ";
    #
    #   if ( $prokMSAname ne $organism ) {
    #       print STDERR " ($prokMSAname) ";
    #   }
    #   print STDERR " \n ";

    return ( $sp, $str, $info );
}

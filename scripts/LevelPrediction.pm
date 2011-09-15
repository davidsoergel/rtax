package LevelPrediction;

use strict;
use warnings;

sub new {
	my ($this) = @_;
    my $class = ref($this) || $this;
    my $self = {};

    bless $self, $class;

    return $self;
}

for my $field (
    qw(label numChildren count localMinProportion localMaxProportion globalMinProportion globalMaxProportion alternateLocalProportion alternateGlobalProportion)
    )
{
    my $slot = __PACKAGE__ . "::$field";
    no strict "refs";    # So symbolic ref to typeglob works.

    *$field = sub {
        my $self = shift;
        $self->{$slot} = shift if @_;
        return $self->{$slot};
    };
}

sub print_tab_delimited {
    my $self = shift;
    for my $field (
        qw(label numChildren count localMinProportion localMaxProportion globalMinProportion globalMaxProportion alternateLocalProportion alternateGlobalProportion)
        )
    {
        my $slot = __PACKAGE__ . "::$field";
        print "\t" . ( $self->{$slot} ? $self->{$slot} : "undef" );
    }
}

sub toString {
    my $self = shift;
	my $result = "";
	
	for my $field (
        qw(label)
        )
    {
        my $slot = __PACKAGE__ . "::$field";
        $result .=  ( $self->{$slot} ? $self->{$slot} : "" );
    }
	
	for my $field (
        qw(numChildren count)
        )
    {
        my $slot = __PACKAGE__ . "::$field";
        $result .= ";" . ( $self->{$slot} ? $self->{$slot} : "" );
    }
	
    for my $fieldb (
        qw(localMinProportion localMaxProportion globalMinProportion globalMaxProportion alternateLocalProportion alternateGlobalProportion)
        )
    {
        my $slot = __PACKAGE__ . "::$fieldb";
        $result .=  $self->{$slot} ? (sprintf ";%.2f", $self->{$slot}) : ";" ;
    }
	return $result;
}

1;

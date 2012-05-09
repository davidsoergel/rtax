#!/usr/bin/env perl

while (<STDIN>) {
	my $line = $_;
	$line =~ s/;[^\t\n]*//g;
	print $line;
}

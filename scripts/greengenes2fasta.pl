use strict;
use warnings;

my $prokmsa;
my $seq;

while(<STDIN>)
	{
	# rely on the order; ignore BEGIN and END tags entirely
	if(/prokMSA_id=(.*)/)
		{
		$prokmsa = $1;
		}
	elsif(/aligned_seq=(.*)/)
		{
		$seq = $1;
		if($seq ne "unaligned")
			{
			print ">$prokmsa\n$seq\n";
			}
		undef $prokmsa;
		}
	}
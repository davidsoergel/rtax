package FastaIndex;

use strict;
use warnings;

use IO::File;


BEGIN {
for my $field (qw( fastaFileName db fh ))
	{
    my $slot = __PACKAGE__ . "::$field";
    no strict "refs";          # So symbolic ref to typeglob works.

    *$field = sub
    	{
        my $self = shift;
        $self->{$slot} = shift if @_;
        return $self->{$slot};
    	};
    }
}

sub new {
	my ($this) = @_;
	my $class = ref($this) || $this;
	my $self = {};
	
	bless $self, $class;
	
	return $self;
}

sub make_index {
    
	my ($this, $fastaFileName) = @_;
    
	$this->fastaFileName($fastaFileName);
	
	# for now just store a hash in memory.  Could tie to a dbm instead?
	$this->db({});
	    
	my $lastpos = 0;
	my $pos = 0;
	#my $seq = "";
	my $id = "";
	my $numlines = 0;
	
	#open(IN, "$fastaFileName") or die "could not open $fastaFileName";
	my $in = IO::File->new("$fastaFileName") or die "could not open $fastaFileName";
	$this->fh($in);
	while(<$in>)
	    {
	    my $line = $_;
	    if($line =~ /^>(\S*)/)
	        {
	            # write the last record
	            my @rec =  ($lastpos, $numlines);        
	            $this->db()->{$id} = \@rec;
	            
	            # start a new record
	            $id = $1;
	            #$seq = $line;
	            $lastpos = $pos;
	            $numlines = 0;
	        }
	    # else { $seq .= $line; }
            
	    $pos += length($line);
	    $numlines++;
        }
}

sub count_records()
    {
    my ($this) = @_;
    return scalar(keys %{$this->db()});
    }

sub fetch {
    
	my ($this, $id) = @_;
    
    my $r = $this->db()->{$id};
    if(! defined $r) { die "no sequence with id: $id"; }
    my ($pos, $numlines) = @{$r};
    
    my $result = "";
    my $in = $this->fh();
    #print $this->fastaFileName() . ": " . $in . "\n";
    seek($in,$pos,0);
    for (1..$numlines) { $result .= <$in> }
    return $result;
}

sub close {
	my ($this) = @_;
    my $in = $this->fh();
    close($in);
    $this->db(undef);
}
1;
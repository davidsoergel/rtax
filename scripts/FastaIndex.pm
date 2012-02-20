package FastaIndex;

use strict;
use warnings;

use IO::File;


BEGIN {
for my $field (qw( fastaFileName startpos lines fh ))
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
    
	my ($this, $fastaFileName, $idregex, $dbmFileName) = @_;
	if(!defined $idregex) { $idregex = "(\\S+)"; }
    
	$this->fastaFileName($fastaFileName);
	
	#open(IN, "$fastaFileName") or die "could not open $fastaFileName";
	my $in = IO::File->new("$fastaFileName") or die "could not open $fastaFileName";
	$this->fh($in);
	
	if($dbmFileName)
	    {
	    my %openpos = ();
	    dbmopen( %openpos, $dbmFileName.".pos", 0666) || die ("Could not open DBM file $dbmFileName.pos");
	    $this->startpos(\%openpos);
	    
	     my %openlines = ();
	    dbmopen( %openlines, $dbmFileName.".lines", 0666) || die ("Could not open DBM file $dbmFileName.lines");
	    $this->lines(\%openlines);
	    
	    my $sizepos = scalar keys %openpos;
	    my $sizelines = scalar keys %openlines;
	    if($sizepos != $sizelines) { die "DBM files broken $dbmFileName\n"; }
	    if($sizepos > 0) {
	        print STDERR ("Index file $dbmFileName.pos already exists; assuming valid\n");
	        return;  
        };
        }
    else
        {
	    # just store a hash in memory.
	    $this->startpos({});
	    $this->lines({});
        }
	    
	my $lastpos = 0;
	my $pos = 0;
	#my $seq = "";
	my $id = "";
	my $numlines = 0;
	
	while(<$in>)
	    {
	    my $line = $_;
	    if($line =~ /^>$idregex/)
	        {
	            if(defined $id && $id ne "") {
	                # write the previous record        
	                print STDERR "$id -> $lastpos\n";
	                $this->startpos()->{$id} = $lastpos;
	                $this->lines()->{$id} = $numlines;
                   }
	            
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
    # write the final record       
	$this->startpos()->{$id} = $lastpos;
	$this->lines()->{$id} = $numlines;
    #delete($this->startpos()->{""});  # spurious entries from the start of the loop
    #delete($this->lines()->{""});  # spurious entries from the start of the loop
}

sub count_records()
    {
    my ($this) = @_;
    return scalar(keys %{$this->startpos()});
    }

sub fetch {
    
	my ($this, $id) = @_;
    
    my $pos = $this->startpos()->{$id};
    if(! defined $pos) { die "no sequence with id: $id"; }
    
    my $numlines = $this->lines()->{$id};
    
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
    $this->startpos(undef);
    $this->lines(undef);
}
1;
RTAX
====

_Rapid and accurate taxonomic classification of short paired-end sequence reads from the 16S ribosomal RNA gene_

Short-read technologies for microbial community profiling are increasingly popular, yet previous techniques for assigning taxonomy to paired-end reads perform poorly.  RTAX provides rapid taxonomic assignments of paired-end reads using a consensus algorithm.

David A. W. Soergel (1*), Rob Knight (2), and Steven E. Brenner (1)

1 Department of Plant and Microbial Biology, University of California, Berkeley 
2 Howard Hughes Medical Institute and Department of Chemistry and Biochemistry, University of Colorado at Boulder

\* Corresponding author (current address): [soergel@cs.umass.edu](mailto:soergel@cs.umass.edu)

Download
--------

[rtax-0.984.tgz](http://static.davidsoergel.com/rtax-0.984.tgz) (23 KB)  _August 7, 2013_
[change log](changelog)

[rtax.greengenes.20110311.tgz](http://static.davidsoergel.com/rtax.greengenes.20110311.tgz) (49 MB)  ''March 11, 2011''

Citation
--------

Soergel DAW, Dey N, Knight R, Brenner SE. (2012). [Selection of primers for optimal taxonomic classification of environmental 16S rRNA gene sequences](http://www.nature.com/ismej/journal/vaop/ncurrent/abs/ismej2011208a.html).  ISME J (6), 1440–1444


Use within QIIME
----------------

RTAX can be used within [QIIME](http://qiime.org) workflows; see this [tutorial](http://www.qiime.org/tutorials/rtax.html).

Requirements
------------

 * Hardware: memory requirements vary with the size of the reference database; for the version of GreenGenes we provide, ~4GB is needed.

 * USEARCH (http://www.drive5.com/usearch/).  RTAX has been tested with v4.1.93 and v5.0.144.  Make sure that the program is in your search path as "usearch", or provide its path explicitly in the "rtax" script.  Note that v5.x uses quite a bit more memory, so 4.x may be preferred.

* A reference database consisting of sequences for which taxonomy assignments are known.  We provide prepared a version of the GreenGenes database on the RTAX web site; see the GreenGenes section below for details.

    Two files are needed:
    1. a FASTA file containing the sequences, and 
    2. a file listing taxonomy strings for each sequence ID found in the FASTA file.  The format of this file is `sequenceid	pcidwidth	kingdom; phylum; class; order; ...` one entry per line, where sequenceid, pcidwidth, and the taxonomy string are separated by a single tab, and the taxonomy string itself is delimited either by semicolons or by tabs. The pcidwidth column is optional, and (if present) is ignored in the current version of RTAX anyway. (In our usage, we cluster the reference database into 99% id OTUs, representing each by a single "seed" sequence. This column then lists the largest pairwise %id difference between the cluster seed and any cluster member.)


Installation
------------

RTAX consists of a set of perl scripts and a shell script ("rtax") to wire them
together. No compilation is needed; just extract the package to a convenient
location.

The perl scripts must remain in the "scripts" directory below the "rtax" shell
script, but the latter can be symlinked anywhere for convenience. A common
pattern might be to place the whole package at /usr/local/rtax, and then symlink
the script into your path (e.g., "ln -s /usr/local/rtax/rtax /usr/local/bin").


Running RTAX
------------

Sequence classification is done as follows:

```
    rtax -r gg.nr.fasta -t gg.nr.taxonomy -a queryA.fasta [-b queryB.fasta] -o result.out 
```

Substitute a different reference database for the GreenGenes files if desired,
of course.

If two query files are provided, then they must contain mate-paired reads
with exactly matching identifiers (though they need not be in the same order).
Any ids present in one file but not the other are silently dropped. Alternate
naming schemes for the paired reads (e.g., "SOMEID.a" paired with "SOMEID.b",
and so forth) are handled via the -i option.

RTAX may be run for single-ended reads as well; simply exclude the queryB
parameter in this case.

Run "rtax" with no arguments for a help message describing additional options.

Various progress messages are provided on stderr, and the predicted
classification for each query is provided in the file given by the -o option.
The output format is essentially the same as the taxonomy file input format
above. The second column indicates the best %id observed between the query
sequence and any reference sequence, and the taxonomy string is tab-delimited.

Temporary files are created in /tmp/rtax.SOMENUMBER (or under the temp
directory given by -m).  These can be quite large; please use -m to select a
location with sufficient free space, as needed.  The temporary directory is
deleted upon successful termination, but may need to be manually cleaned up in
the event of an error of some sort.


Using the GreenGenes reference database
---------------------------------------

An RTAX-formatted reference database based on GreenGenes (version of March 11,
2011) is available from the Downloads section above. That
contains gg.nr.fasta and gg.nr.taxonomy, which are the result of clustering the
GreenGenes input file at 99% identity, finding consensus taxonomy strings for
each cluster, and formatting the result as needed for use with RTAX.

Please see the "greengenes" subdirectory in the RTAX distribution for
information on using newer versions of GreenGenes.


Preparing a non-GreenGenes reference database
----------------------------------------------

You can certainly use different databases, or GreenGenes clustered by different
means or with different parameters (or not at all, though this approach will
have poor performance). If you have any trouble producing the reference fasta
and taxonomy files in the required format, examine the prepare-greengenes script
for hints, and feel free to contact us for help.


Plumbing
--------

Taxonomy assignment proceeds in two phases, which we implement separately:

1.  For each query sequence (or mate pair), find an appropriate set of matching hits in a reference database.  This is implemented as rtaxSearchSingle and rtaxSearchPair for single and paired-end reads, respectively.

2.  Find consensus taxonomy among each set of reference hits.  This is implemented as rtaxVote.

3.  Finally the detailed rtaxVote results are filtered and cleaned up for output.


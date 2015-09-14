==============
 swarmwrapper
==============

Dereplicate, pool, and cluster reads using swarm
(https://github.com/torognes/swarm), producing output suitable for
subsequent analysis by pplacer.

Input sequence names should *not* contain abundance annotations (these
are added before the sequences are provided to swarm). Output sequence
names are unannotated as well. Input files may be compressed
(compression is detected according to a file sufffix of either .bz2 or
.gz).

The output file specified by -a/--abundances contains three columns:

 1. the name of a seed sequence representing cluster C;
 2. the name of a sequence included in cluster C representing some specimen S;
 3. the number of reads representing cluster C originating from specimen S.

Let's name these columns ('seed', 'read_from_S', and 'abundance'). In
principle, the sum of the values in the column containing abundances
should equal the total number of input reads. However, swarm does not
allow sequences containing ambiguities, so these are silently
discarded.

In the specific case of pplacer, this file can be provided as an
argument to ``guppy redup -d`` after placement of sequences in
seeds.fasta to generate a placefile reflecting the original read
masses before clustering.

But more generally, this output, along with a specimen map and
assignment of seed sequences to taxon names, can be used to construct
a table describing taxon abundance by specimen.

Consider these two other tables:

 * 'specimen_map', a table provided as input to this script with
   columns ('read_from_S', 'specimen')
 * 'assignments', a table assigning clusters to taxa with columns
   ('seed', 'taxon').

Given these three inputs, a taxon table can be constructed as follows
(using SQL to illustrate the relations)::

  select specimen, taxon, sum(abundance) as read_count
  from abundances join specimen_map using (read_from_S)
  join assignments using(seed)
  group by specimen, taxon;

::

  positional arguments:
    infile                Input file containing trimmed reads in fasta format

  optional arguments:
    -h, --help            show this help message and exit
    -s SPECIMEN_MAP, --specimen-map SPECIMEN_MAP
			  headless csv file with columns (read, specimen)
    -d DIFFERENCES, --differences DIFFERENCES
    -k, --keep-abundances
			  retain abundance annotation in sequence names
			  (-w/--seeds only)
    --min-mass N          drop OTUs with total mass less than N
    -t THREADS, --threads THREADS
    --version             show program's version number and exit

  output files:
    -w SEEDS, --seeds SEEDS
			  Output fasta file containing OTU representatives
    -a ABUNDANCES, --abundances ABUNDANCES
			  csv file providing abundances by specimen
    --dropped DROPPED     file containing sequences dropped due to ambiguities
			  or cluster size
    --tmpdir TMPDIR       optional directory name for creating temporary
			  intermediate files (created in system temp file and
			  discarded by default)

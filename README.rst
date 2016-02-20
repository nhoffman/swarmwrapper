==============
 swarmwrapper
==============

Dereplicate, pool, and cluster reads using swarm
(https://github.com/torognes/swarm), producing output suitable for
subsequent analysis by pplacer (http://matsen.github.io/pplacer/).

Why?

* mostly, for compatibility with pplacer-style abundance annotation
  for pooled specimens
* adds abundance annotation to raw reads
* drops reads with ambiguities
* optionally drops OTUs with a mass below some threshold
* reads and writes compressed sequence and data files

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

Installation
============

Recommended use is in a python virtualenv. Here's what I do to install
both swarmwrapper and the linux swam binaries to a virtualenv::

  virtualenv env
  source env/bin/activate
  pip install git+https://github.com/nhoffman/swarmwrapper.git
  swarmwrapper install --prefix $VIRTUAL_ENV


Command line options
====================

::

  positional arguments:
    {dereplicate,cluster,install}
      dereplicate         Perform strict dereplication.
      cluster             Cluster dereplicated reads.
      install             Install swarm binaries (linux only for now)

  optional arguments:
    -h, --help            show this help message and exit
    -v                    increase verbosity of screen output (eg, -v is
			  verbose, -vv more so)
    -q, --quiet           suppress screen output from pip commands
    -V, --version         Print the version number and exit
    --tmpdir TMPDIR       optional directory name for creating temporary
			  intermediate files (created in system temp file and
			  discarded by default)
    -t N, --threads N     number of threads [4]

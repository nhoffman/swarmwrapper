#!/usr/bin/env python

"""Dereplicate, pool, and cluster reads using swarm
(https://github.com/torognes/swarm), producing output suitable for
subsequent analysis by pplacer.

Why a wrapper for a perfectly nice program like swarm?

* mostly, for compatibility with pplacer-style abundance annotation
  with pooled specimens
* adds abundance annotation to raw reads
* drops reads with ambiguities
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

Command line options
====================

::
"""

from __future__ import print_function
from collections import namedtuple, defaultdict
from tempfile import NamedTemporaryFile
from itertools import ifilter, imap
from distutils.version import LooseVersion
from os import path
import textwrap
import argparse
import csv
import gzip
import os
import subprocess
import sys
import logging
try:
    import bz2
except ImportError:
    bz2 = None


try:
    with open(path.join(path.dirname(__file__), 'data', 'ver')) as f:
        __version__ = f.read().strip().replace('-', '+', 1).replace('-', '.')
except Exception, e:
    __version__ = ''

log = logging

SWARM_VERSION = '2.1.4'


def mkdir(pth):
    try:
        os.makedirs(pth)
    except OSError:
        pass

    return path.abspath(pth)


class Opener(object):

    """Factory for creating file objects

    Keyword Arguments:
        - mode -- A string indicating how the file is to be opened. Accepts the
            same values as the builtin open() function.
        - bufsize -- The file's desired buffer size. Accepts the same values as
            the builtin open() function.
    """

    def __init__(self, mode='r', bufsize=-1):
        self._mode = mode
        self._bufsize = bufsize

    def __call__(self, string):
        if string is sys.stdout or string is sys.stdin:
            return string
        elif string == '-':
            return sys.stdin if 'r' in self._mode else sys.stdout
        elif string.endswith('.bz2'):
            if bz2 is None:
                raise ImportError(
                    'could not import bz2 module - was python built with libbz2?')
            return bz2.BZ2File(
                string, self._mode, self._bufsize)
        elif string.endswith('.gz'):
            return gzip.open(
                string, self._mode, self._bufsize)
        else:
            return open(string, self._mode, self._bufsize)

    def __repr__(self):
        args = self._mode, self._bufsize
        args_str = ', '.join(repr(arg) for arg in args if arg != -1)
        return '{}({})'.format(type(self).__name__, args_str)


SeqLite = namedtuple('SeqLite', 'id, description, seq')


def fastalite(handle):
    """
    Lightweight fasta parser. Returns iterator of namedtuple instances
    with fields (id, description, seq) given file-like object ``handle``.
    """

    name, seq = '', ''
    for line in handle:
        if line.startswith('>'):
            if name:
                yield SeqLite(name.split()[0], name, seq)
            name, seq = line[1:].strip(), ''
        else:
            seq += line.strip()

    if name and seq:
        yield SeqLite(name.split()[0], name, seq)


def check_swarm_version(min_version):
    try:
        subprocess.check_output(['which', 'swarm'])
    except subprocess.CalledProcessError:
        sys.exit('Error: no swarm executable found')

    output = subprocess.check_output('swarm -v 2>&1; true', shell=True)
    version = output.split()[1].decode(encoding='UTF-8')
    if LooseVersion(version) < LooseVersion(min_version):
        sys.exit('Error: swarm version >= {} is required, '
                 'found swarm version {}'.format(min_version, version))


def add_abundance(seq, abundance=1):
    return SeqLite('{}_{}'.format(seq.id, abundance), seq.description, seq.seq)


def rm_abundance(seq):
    name, __ = get_abundance(seq.id)
    return SeqLite(name, seq.description, seq.seq)


def get_abundance(name):
    key, abundance = name.rsplit('_', 1)
    return key, int(abundance)


def fiter_ambiguities(seq):
    return 'N' not in seq.seq


def write_seqs(fobj, seqs):
    for seq in seqs:
        fobj.write(as_fasta(seq))
        fobj.flush()


def as_fasta(seq):
    return '>{}\n{}\n'.format(seq.id, seq.seq)


def swarm(infile, seeds, clusters=None, differences=0, threads=1, quiet=True):
    with open(os.devnull, 'w') as devnull:
        cmd = ['swarm',
               '--seeds', seeds.name,
               '-o', clusters.name if clusters else os.devnull,
               '--differences', str(differences),
               '-t', str(threads),
               infile.name]

        log.info(' '.join(cmd))
        subprocess.check_call(cmd, stderr=devnull if quiet else None)


def dereplicate(seqs, seeds, tmpdir=None, threads=1, quiet=True):
    """Given iterator of sequences ``seqs``, write dereplicated reads in
    fasta format to open file ``seeds``.

    """

    seqs = imap(add_abundance, seqs)
    with ntf('w', prefix='raw-', suffix='.fasta', dir=tmpdir) as raw, \
         ntf('rwb', prefix='d0-', suffix='.fasta', dir=tmpdir) as d0:
        write_seqs(raw, seqs)
        swarm(raw, d0, differences=0, threads=threads, quiet=quiet)
        seeds.write(d0.read())


def dereplicate_and_pool(seqs, specimen_map, seeds, tmpdir=None,
                         threads=1, quiet=True):
    """Given iterator of sequences ``seqs``, and dict ``specimen_map``
    (providing the mapping {read: specimen}), write dereplicated reads
    in fasta format to open file ``seeds``.

    """

    # Create an open file for each specimen. For each read, write a
    # sequences annotated with _1 (abundance = 1) if there are no
    # ambiguities.
    raw_reads = {specimen: ntf('w', prefix='{}-'.format(specimen),
                               suffix='.fasta', dir=tmpdir, delete=False)
                 for specimen in set(specimen_map.values())}

    for seq in seqs:
        raw_reads[specimen_map[seq.id]].write(as_fasta(add_abundance(seq)))

    # dereplicate each specimen and concatenate all remaining reads
    for specimen, infile in list(raw_reads.items()):
        with ntf(prefix='d0-', suffix='.fasta', dir=tmpdir) as d0:
            infile.close()  # flush any pending writes from above
            swarm(infile, d0, differences=0, threads=threads, quiet=quiet)

            # concatenate to pooled file
            seeds.write(d0.read())

            if not tmpdir:
                os.remove(infile.name)


def ntf(*args, **kwargs):
    tmpdir = kwargs.get('dir')
    kwargs['delete'] = kwargs['delete'] if 'delete' in kwargs else (tmpdir is None)

    if tmpdir is not None:
        try:
            os.makedirs(tmpdir)
        except OSError:
            pass

    return NamedTemporaryFile(*args, **kwargs)


class VersionAction(argparse._VersionAction):
    """Write the version string to stdout and exit"""
    def __call__(self, parser, namespace, values, option_string=None):
        formatter = parser._get_formatter()
        formatter.add_text(parser.version if self.version is None else self.version)
        sys.stdout.write(formatter.format_help())
        sys.exit(0)


class Subparser(object):
    def __init__(self, subparsers, name):
        self.subparser = subparsers.add_parser(
            name,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            help=self.__doc__.strip().split('\n')[0],
            description=textwrap.dedent(self.__doc__.rstrip()))
        self.subparser.set_defaults(func=self.action)
        self.add_arguments()


class Cluster(Subparser):
    """
    Cluster dereplicated reads.
    """

    def add_arguments(self):
        self.subparser.add_argument(
            'seqs', type=Opener('r'),
            help="input sequences in fasta format "
            "(dereplicated and with abundance annotations)")
        self.subparser.add_argument(
            '-m', '--specimen-map', type=Opener('r'), metavar='INFILE',
            help="csv file with columns (read_name, specimen)")
        self.subparser.add_argument(
            '-w', '--seeds', type=Opener('w'), metavar='OUTFILE',
            help="output seed sequences in fasta format")
        self.subparser.add_argument(
            '-a', '--abundances', type=Opener('w'), metavar='OUTFILE',
            help="csv file providing abundances by specimen")
        self.subparser.add_argument(
            '--dropped', type=Opener('w'), metavar='OUTFILE',
            help="sequences discarded due to --min-mass threshold")
        self.subparser.add_argument(
            '-d', '--differences', type=int, default=1, metavar='N',
            help='value for "swarm -d" [default %(default)s]')
        self.subparser.add_argument(
            '-M', '--min-mass', type=int, default=None, metavar='N',
            help="drop OTUs with total mass less than N")
        self.subparser.add_argument(
            '-k', '--keep-abundance', action='store_true', default=False,
            help="keep abundance annotation in seed names")

    def action(self, args):
        check_swarm_version(SWARM_VERSION)
        # identifies specimen of origin (values) for each read (keys)
        specimen_map = dict(csv.reader(args.specimen_map))

        if args.abundances:
            writer = csv.writer(args.abundances)

        with ntf(prefix='infile-', suffix='.fasta', dir=args.tmpdir) as infile, \
             ntf(prefix='clusters-', suffix='.txt', dir=args.tmpdir) as clusters, \
             ntf(prefix='seeds-', suffix='.fasta', dir=args.tmpdir) as seeds:

            # copy contents of args.seqs to infile before running
            # swarm to handle compressed files
            infile.write(args.seqs.read())
            infile.flush()

            swarm(infile, seeds, clusters, differences=args.differences,
                  threads=args.threads, quiet=args.verbosity <= 1)

            grand_total, keep_total = 0.0, 0.0
            for seq, line in zip(fastalite(seeds), clusters):
                otu_rep, total = get_abundance(seq.id)
                names_and_counts = list(map(get_abundance, line.split()))

                assert otu_rep == names_and_counts[0][0]
                assert total == sum(x[1] for x in names_and_counts)

                grand_total += total

                if args.min_mass is not None and total < args.min_mass:
                    if args.dropped:
                        args.dropped.write(as_fasta(seq))
                    continue

                keep_total += total

                # write the OTU seed
                if not args.keep_abundance:
                    seq = rm_abundance(seq)
                args.seeds.write(as_fasta(seq))

                if not args.abundances:
                    continue

                # associate each specimen with a list of read names and masses
                specimens = defaultdict(list)
                for name, count in names_and_counts:
                    specimens[specimen_map[name]].append((name, count))

                # start with the specimen corresponding to the OTU representative
                specimen = specimen_map[otu_rep]
                names, counts = list(zip(*specimens.pop(specimen)))
                writer.writerow([otu_rep, otu_rep, sum(counts)])

                # ... then iterate over other specimens with reads in this OTU
                for specimen, val in specimens.items():
                    names, counts = list(zip(*val))
                    writer.writerow([otu_rep, names[0], sum(counts)])

        log.warning('total yield: {}'.format(round(100.0 * keep_total/grand_total, 2)))


class Dereplicate(Subparser):
    """
    Perform strict dereplication.

    Remove sequences with ambiguities, add abundance annotation, and
    dereplicate with ``swarm -d 0``. If a specimen_map is provided,
    reads from each specimen will be dereplicated individually.

    """

    def add_arguments(self):
        self.subparser.add_argument(
            'seqs', type=Opener('r'), help="sequences in fasta format")
        self.subparser.add_argument(
            '-m', '--specimen-map', type=Opener('r'), metavar='INFILE',
            help="csv file with columns (read_name, specimen)")
        self.subparser.add_argument(
            '-w', '--seeds', type=Opener('w'), metavar='OUTFILE',
            help="output seed sequences in fasta format")

    def action(self, args):
        check_swarm_version(SWARM_VERSION)
        seqs = fastalite(args.seqs)
        seqs = ifilter(fiter_ambiguities, seqs)

        if args.specimen_map:
            # identifies specimen of origin (values) for each read (keys)
            specimen_map = dict(csv.reader(args.specimen_map))
            dereplicate_and_pool(seqs, specimen_map, args.seeds, tmpdir=args.tmpdir,
                                 threads=args.threads, quiet=args.verbosity <= 1)
        else:
            dereplicate(seqs, args.seeds, tmpdir=args.tmpdir,
                        threads=args.threads, quiet=args.verbosity <= 1)


def main(arguments=None):

    if arguments is None:
        arguments = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        '-v', action='count', dest='verbosity', default=1,
        help='increase verbosity of screen output (eg, -v is verbose, '
        '-vv more so)')
    parser.add_argument(
        '-q', '--quiet', action='store_const', dest='verbosity', const=0,
        help='suppress screen output from pip commands')
    parser.add_argument(
        '-V', '--version', action=VersionAction, version=__version__,
        help='Print the version number and exit')
    parser.add_argument(
        '--tmpdir', help="""optional directory name for creating
        temporary intermediate files (created in system temp file and
        discarded by default)""")
    parser.add_argument(
        '-t', '--threads', default=4, type=int, metavar='N',
        help='number of threads [%(default)s]')

    subparsers = parser.add_subparsers()
    Dereplicate(subparsers, name='dereplicate')
    Cluster(subparsers, name='cluster')

    args = parser.parse_args(arguments)

    # set up logging
    loglevel = {
        0: logging.ERROR,
        1: logging.WARNING,
        2: logging.INFO,
        3: logging.DEBUG,
    }.get(args.verbosity, logging.DEBUG)

    logformat = ('%(levelname)s %(funcName)s %(lineno)s %(message)s'
                 if args.verbosity > 1
                 else '%(message)s')
    logging.basicConfig(file=sys.stderr, format=logformat, level=loglevel)

    return args.func(args)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

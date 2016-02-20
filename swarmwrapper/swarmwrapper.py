#!/usr/bin/env python

"""Dereplicate, pool, and cluster reads using swarm
(https://github.com/torognes/swarm), producing output suitable for
subsequent analysis by pplacer.

Why a wrapper for a perfectly nice program like swarm?

* mostly, for compatibility with pplacer-style abundance annotation
  with pooled specimens
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
import urllib2

try:
    from bz2 import BZ2File
except ImportError, err:
    BZ2File = lambda x, *args, **kwargs: sys.exit(err)

try:
    with open(path.join(path.dirname(__file__), 'data', 'ver')) as f:
        __version__ = f.read().strip().replace('-', '+', 1).replace('-', '.')
except Exception, e:
    __version__ = ''

log = logging

SWARM_VERSION = '2.1.6'


def mkdir(pth):
    try:
        os.makedirs(pth)
    except OSError:
        pass

    return path.abspath(pth)


class Opener(object):
    """Factory for creating file objects. Transparenty opens compressed
    files for reading or writing based on suffix (.gz and .bz2 only).

    Example::

        with Opener()('in.txt') as infile, Opener('w')('out.gz') as outfile:
            outfile.write(infile.read())
    """

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.writable = 'w' in kwargs.get('mode', args[0] if args else 'r')

    def __call__(self, obj):
        if obj is sys.stdout or obj is sys.stdin:
            return obj
        elif obj == '-':
            return sys.stdout if self.writable else sys.stdin
        else:
            __, suffix = obj.rsplit('.', 1)
            opener = {'bz2': BZ2File,
                      'gz': gzip.open}.get(suffix, open)
            return opener(obj, *self.args, **self.kwargs)

Seq = namedtuple('Seq', 'id, description, seq')


def fastalite(handle):
    """Return a sequence of namedtuple objects with attributes (id,
    description, seq) given open file-like object ``handle``

    """

    header, seq = '', []
    for line in handle:
        if line.startswith('>'):
            if header:
                yield Seq(header.split()[0], header, ''.join(seq))
            header, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())

    if header and seq:
        yield Seq(header.split()[0], header, ''.join(seq))


def check_swarm_version(min_version, swarm=None):
    if swarm is None:
        try:
            swarm = subprocess.check_output(['which', 'swarm']).strip()
        except subprocess.CalledProcessError:
            sys.exit('Error: no swarm executable found')
    else:
        swarm = path.abspath(swarm)

    # swarm -v always returns non-zero exit status
    output = subprocess.check_output('"{}" -v 2>&1; true'.format(swarm), shell=True)
    if output.lower().startswith('swarm'):
        version = output.split()[1].decode(encoding='UTF-8')
    else:
        log.info('No swarm executable at path "{}"'.format(swarm))
        return False

    version_ok = LooseVersion(version) >= LooseVersion(min_version)
    if version_ok:
        log.info('{} version {}'.format(swarm, version))
    else:
        log.error('Error: swarm version >= {} is required, '
                  'found {} version {}'.format(min_version, swarm, version))

    return version_ok


def add_abundance(seq, abundance=1):
    return Seq('{}_{}'.format(seq.id, abundance), seq.description, seq.seq)


def rm_abundance(seq):
    name, __ = get_abundance(seq.id)
    return Seq(name, seq.description, seq.seq)


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
            '-D', '--dereplicate', action='store_true', default=False,
            help='dereplicate before clustering')
        self.subparser.add_argument(
            '-d', '--differences', type=int, default=1, metavar='N',
            help='value for "swarm -d" [default %(default)s]')
        self.subparser.add_argument(
            '-M', '--min-mass', type=int, default=None, metavar='N',
            help="drop OTUs with total mass less than N [default 1, ie, no dropping]")
        self.subparser.add_argument(
            '-k', '--keep-abundance', action='store_true', default=False,
            help="keep abundance annotation in seed names")

    def action(self, args):
        if not check_swarm_version(SWARM_VERSION):
            sys.exit(1)

        # identifies specimen of origin (values) for each read (keys)
        if args.specimen_map:
            specimen_map = dict(csv.reader(args.specimen_map))
        else:
            specimen_map = defaultdict(lambda: 'specimen')

        if args.abundances:
            writer = csv.writer(args.abundances)

        with ntf(prefix='infile-', suffix='.fasta', dir=args.tmpdir) as infile, \
             ntf(prefix='clusters-', suffix='.txt', dir=args.tmpdir) as clusters, \
             ntf(prefix='seeds-', suffix='.fasta', dir=args.tmpdir) as seeds:

            if args.dereplicate:
                seqs = ifilter(fiter_ambiguities, fastalite(args.seqs))
                if args.specimen_map:
                    dereplicate_and_pool(
                        seqs, specimen_map, infile, tmpdir=args.tmpdir,
                        threads=args.threads, quiet=args.verbosity <= 1)
                else:
                    dereplicate(
                        seqs, infile, tmpdir=args.tmpdir,
                        threads=args.threads, quiet=args.verbosity <= 1)
            else:
                # copy contents of args.seqs to infile before running
                # swarm to handle compressed files
                infile.write(args.seqs.read())

            # cluster
            infile.flush()
            swarm(infile, seeds, clusters, differences=args.differences,
                  threads=args.threads, quiet=args.verbosity <= 1)

            grand_total, keep_total = 0.0, 0.0
            for seq, line in zip(fastalite(seeds), clusters):
                otu_rep, total = get_abundance(seq.id)
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

                names_and_counts = list(map(get_abundance, line.split()))
                assert otu_rep == names_and_counts[0][0]
                assert total == sum(x[1] for x in names_and_counts)

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
        if not check_swarm_version(SWARM_VERSION):
            sys.exit(1)
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


class Install(Subparser):
    """
    Install swarm binaries (linux only for now)

    """

    def add_arguments(self):
        self.subparser.add_argument(
            '--prefix',
            help='install swarm to PREFIX/bin/swarm (destination directory must exist)')
        self.subparser.add_argument(
            '--path', default='./swarm',
            help='if no --prefix, install swarm to PATH [default "%(default)s]"')
        self.subparser.add_argument(
            '--version', default=SWARM_VERSION,
            help='swarm version [default %(default)s]')
        self.subparser.add_argument(
            '--force', action='store_true', default=False,
            help='redownload even if the binary exists')

    def action(self, args):
        url = ('http://github.com/torognes/swarm/releases/download/'
               'v{version}/swarm-{version}-linux-x86_64').format(version=args.version)
        dest = '{}/bin/swarm'.format(args.prefix) if args.prefix else args.path
        version_ok = check_swarm_version(args.version, swarm=dest)
        if not version_ok or args.force:
            log.info('downloading {} to {}'.format(url, dest))
            with open(dest, 'w') as binary:
                handle = urllib2.urlopen(url)
                binary.write(handle.read())
            os.chmod(dest, 0755)
            # confirm installation
            check_swarm_version(args.version, swarm=dest)


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
    Install(subparsers, name='install')

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

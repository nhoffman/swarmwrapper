#!/usr/bin/env python

"""Dereplicate, pool, and cluster reads using swarm
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
                raise ImportError('could not import bz2 module - was python built with libbz2?')
            return bz2.BZ2File(string, self._mode, self._bufsize)
        elif string.endswith('.gz'):
            return gzip.open(string, self._mode, self._bufsize)
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


def get_abundance(name):
    key, abundance = name.rsplit('_', 1)
    return key, int(abundance)


def fiter_ambiguities(seq):
    return 'N' not in seq.seq


def write_seqs(fobj, seqs):
    for seq in seqs:
        fobj.write(as_fasta)
        fobj.flush()


def as_fasta(seq):
    return '>{}\n{}\n'.format(seq.id, seq.seq)


def dereplicate(infile, outfile, threads=1):
    with open(os.devnull, 'w') as devnull:
        cmd = ['swarm',
               '--seeds', outfile.name,
               '-o', devnull.name,  # ordinarily to stdout, but we don't want it
               '--differences', '0',
               '-t', str(threads),
               infile.name]
        subprocess.check_call(cmd, stderr=devnull)


def ntf(*args, **kwargs):
    tmpdir = kwargs.get('dir')
    kwargs['delete'] = kwargs['delete'] if 'delete' in kwargs else tmpdir is None

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


class Nothing(Subparser):
    """


    """

    def add_arguments(self):
        self.subparser.add_argument(
            'pkg', help="name of a package")
        self.subparser.add_argument(
            '--venv', help="Path to a virtualenv")

    def action(self, args):
        print(args)


class Dereplicate(Subparser):
    """Pre-process (remove sequences with ambiguities and add abundance
    annotation), dereplicate, and cluster sequences in fasta
    format. If a specimen_map is provided, reads from each specimen
    will be dereplicated individually.

    """

    def add_arguments(self):
        self.subparser.add_argument(
            'seqs', type=Opener('r'), help="sequences in fasta format")
        self.subparser.add_argument(
            '-m', '--specimen-map', type=Opener('r'),
            help="csv file with columns (read_name, specimen)")
        self.subparser.add_argument(
            '-w', '--seeds', type=Opener('w'), help="output seed sequences in fasta format")

    def action(self, args):
        seqs = fastalite(args.seqs)
        seqs = ifilter(fiter_ambiguities, seqs)

        if not args.specimen_map:
            log.info('dereplicating {}'.format(args.seqs.name))
            seqs = imap(add_abundance, seqs)
            with ntf('w', prefix='raw-', suffix='.fasta', dir=args.tmpdir) as raw_reads:
                write_seqs(raw_reads, seqs)
                dereplicate(raw_reads, args.seeds, threads=args.threads)
            return

        # identifies specimen of origin (values) for each read (keys)
        specimen_map = dict(csv.reader(args.specimen_map))

        # Create an open file for each specimen. For each read, write a
        # sequences annotated with _1 (abundance = 1) if there are no
        # ambiguities.
        raw_reads = {specimen: ntf('w', prefix='{}-'.format(specimen),
                                   suffix='.fasta', dir=args.tmpdir, delete=False)
                     for specimen in set(specimen_map.values())}

        for seq in seqs:
            raw_reads[specimen_map[seq.id]].write(as_fasta(add_abundance(seq)))

        # dereplicate each specimen and concatenate all remaining reads
        for specimen, infile in list(raw_reads.items()):
            with ntf(prefix='d0-', suffix='.fasta', dir=args.tmpdir) as d0:
                infile.close()  # flush any pending writes from above
                dereplicate(infile, d0, threads=args.threads)

                # concatenate to pooled file
                cmd = 'cat "{}" >> "{}"'.format(d0.name, args.seeds.name)
                subprocess.check_call(cmd, shell=True)


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
    parser.add_argument('-t', '--threads', default=4, type=int)

    subparsers = parser.add_subparsers()
    Nothing(subparsers, name='do-nothing')
    Dereplicate(subparsers, name='dereplicate')

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

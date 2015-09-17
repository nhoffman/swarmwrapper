#!/bin/bash

set -e

for specimen in $(csvcut -H -c 2 seq_info.csv | tail -n+2 | sort | uniq); do
    grep $specimen seq_info.csv | cut -f 1 -d, > want.txt
    seqmagick convert --include-from-file want.txt seqs.fasta some_words_${specimen}_other_words.fasta
done

rm want.txt

#!/bin/bash

set -e

run_test(){
    # return 0 if first argument matches any of the subsequent ones,
    # or if only one argument is provided.
    if [[ -z "$2" ]]; then
	return 0
    fi
    local e
    for e in "${@:2}"; do
	if [[ "$e" == "$1" ]]; then
	    return 0
	fi
    done
    return 1
}

setup(){
    local d
    d=test_output/${1:?}
    >&2 echo "running $1"
    rm -rf $d
    mkdir -p $d
    echo $d
}

verbose="-v"

testname=derep1
if run_test $testname "$@"; then
    outdir=$(setup $testname)
    tmpdir="--tmpdir $outdir/tmp"
    ./swarmwrapper.py "$verbose" $tmpdir \
		      dereplicate \
		      testfiles/seqs.fasta \
		      --seeds $outdir/seeds.fasta.gz
fi

testname=derep2
if run_test $testname "$@"; then
    outdir=$(setup $testname)
    # tmpdir="--tmpdir $outdir/tmp"
    ./swarmwrapper.py "$verbose" $tmpdir \
		      dereplicate \
		      testfiles/seqs.fasta \
		      --specimen-map testfiles/seq_info.csv \
		      --seeds $outdir/seeds.fasta.gz
fi

testname=cluster1
if run_test $testname "$@"; then
    outdir=$(setup $testname)
    tmpdir="--tmpdir $outdir/tmp"
    ./swarmwrapper.py "$verbose" \
		      dereplicate \
		      testfiles/seqs.fasta \
		      --specimen-map testfiles/seq_info.csv \
		      --seeds $outdir/d0.fasta

    ./swarmwrapper.py "$verbose" $tmpdir \
		      cluster \
		      $outdir/d0.fasta \
		      --specimen-map testfiles/seq_info.csv \
		      --seeds $outdir/seeds.fasta
fi

testname=cluster2
if run_test $testname "$@"; then
    outdir=$(setup $testname)
    tmpdir="--tmpdir $outdir/tmp"
    ./swarmwrapper.py "$verbose" \
		      dereplicate \
		      testfiles/seqs.fasta \
		      --specimen-map testfiles/seq_info.csv \
		      --seeds $outdir/d0.fasta

    ./swarmwrapper.py "$verbose" $tmpdir \
		      cluster \
		      $outdir/d0.fasta \
		      --specimen-map testfiles/seq_info.csv \
		      --seeds $outdir/seeds.fasta.gz \
		      --abundances $outdir/weights.csv
fi

testname=cluster3
if run_test $testname "$@"; then
    outdir=$(setup $testname)
    tmpdir="--tmpdir $outdir/tmp"
    ./swarmwrapper.py "$verbose" \
		      dereplicate \
		      testfiles/seqs.fasta \
		      --specimen-map testfiles/seq_info.csv \
		      --seeds $outdir/d0.fasta

    ./swarmwrapper.py "$verbose" $tmpdir \
		      cluster \
		      $outdir/d0.fasta \
		      --specimen-map testfiles/seq_info.csv \
		      --seeds $outdir/seeds.fasta \
		      --abundances $outdir/weights.csv \
		      --min-mass 2
fi

testname=cluster4
if run_test $testname "$@"; then
    outdir=$(setup $testname)
    tmpdir="--tmpdir $outdir/tmp"
    ./swarmwrapper.py "$verbose" \
		      dereplicate \
		      testfiles/seqs.fasta \
		      --specimen-map testfiles/seq_info.csv \
		      --seeds $outdir/d0.fasta

    ./swarmwrapper.py "$verbose" $tmpdir \
		      cluster \
		      $outdir/d0.fasta \
		      --specimen-map testfiles/seq_info.csv \
		      --seeds $outdir/seeds.fasta \
		      --abundances $outdir/weights.csv \
		      --min-mass 2 \
		      --dropped $outdir/dropped.fasta \
		      --keep-abundance
fi

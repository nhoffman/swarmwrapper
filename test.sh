#!/bin/bash

outdir=test_output

run_test(){
    # return 0 if first argument matches any of the subsequent ones,
    # or if only one argument is provided.
    if [[ -z "$2" ]]; then
	return 0
    fi

    for e in "${@:2}"; do
	if [[ "$e" == "$1" ]]; then
	    return 0
	fi
    done
    return 1
}

setup(){
    outdir=${outdir:?}/${1:?}
    >&2 echo "running $1"
    rm -rf $outdir
    mkdir -p $outdir
    echo $outdir
}

verbose="-v"

testname=derep1
if run_test $testname "$@"; then
    outdir=$(setup $testname)
    ./swarmwrapper.py "$verbose" dereplicate testfiles/seqs.fasta --seeds $outdir/seeds.fasta
fi

testname=derep2
if run_test $testname "$@"; then
    outdir=$(setup $testname)
    # tmpdir="--tmpdir $outdir/tmp"
    ./swarmwrapper.py "$verbose" $tmpdir \
		      dereplicate \
		      testfiles/seqs.fasta \
		      --specimen-map testfiles/seq_info.csv \
		      --seeds $outdir/seeds.fasta
fi

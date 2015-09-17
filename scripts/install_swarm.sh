#!/bin/bash

# Install swarm to $prefix/bin

set -e

source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--version', default='2.1.4', help='Swarm version [%(default)s]')
parser.add_argument('--prefix', default='/usr/local', help='base dir for install [%(default)s]')
EOF

version="$VERSION"
prefix=$(readlink -f "$PREFIX")

if $prefix/bin/swarm --version 2>&1 | grep -q $version; then
    echo "Swarm version $version is already installed in $prefix/bin"
    exit 0
fi

cd "$prefix/bin"

# Download swarm binaries
rm -f swarm-${version}-linux-x86_64
wget --no-check-certificate \
     http://github.com/torognes/swarm/releases/download/v${version}/swarm-${version}-linux-x86_64
chmod +x swarm-${version}-linux-x86_64
ln -f swarm-${version}-linux-x86_64 swarm

# confirm success
$prefix/bin/swarm --version 2>&1 | grep $version


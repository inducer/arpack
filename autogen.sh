#! /bin/sh

set -e

libtoolize --automake
aclocal
autoconf
automake --add-missing

./configure "$@"

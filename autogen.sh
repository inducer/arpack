#! /bin/sh

LIBTOOLIZE=libtoolize
if ! which libtoolize &>/dev/null; then
  # mac osx
  LIBTOOLIZE=glibtoolize
fi

set -ex

$LIBTOOLIZE --automake
aclocal
autoconf
automake --add-missing

./configure "$@"

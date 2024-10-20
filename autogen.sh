#!/bin/sh
#
# Run 'autoreconf -i' to generate 'configure', 'Makefile.in', etc.
#
# The first time this is run on a new cloned git repo the configure
# script will not be present, only the configure.ac and
# Makefile.am. The rest must be generated by `autoreconf -i`.
#
# If you are working with a distribution (file ending with ".tar.gz"
# or similar) then this script should not be needed, and should not be
# present, as all the files should already exist. You should only run
# this script if you know what you are doing with autoreconf.
#
# This script will only work with an argument to confirm the help
# message has been read.

runautoreconf() {
    autoreconf -i;
}

if test -d .git && test "$(basename "${PWD}")" = "preseq"
then
    runautoreconf
    exit 0
else
    echo "  It seems you are either attempting to run this script       "
    echo "  from the wrong directory, or in a source tree that was      "
    echo "  not obtained by cloning the preseq git repo.                "
    echo "                                                              "
    echo "  ./autogen.sh generates the configure script. Only run       "
    echo "                                                              "
    echo "  Only run this if you know what you are doing with           "
    echo "  autoreconf and are simply avoiding doing that. If you       "
    echo "  just want to use the software, download a release and       "
    echo "  this script will not be needed.                             "
    exit 1
fi

#!/bin/sh


CORRECT_AC_VERSION="2.69"
CORRECT_AM_VERSION="1.14"


acversion=`autoconf --version | head -n 1 | awk '{print $4}'`
if [ $acversion != $CORRECT_AC_VERSION ]; then
  echo You seem to be using GNU Autoconf version $acversion at `which autoconf`
  echo You should be using version ${CORRECT_AC_VERSION}, which is the latest release as of September 2014.
  echo The correct version is available on ruby at /opt/local/bin/autoconf
  echo or /g/melodi/software/bin/as60-amd64 on the MELODI Linux machines.
  exit 1
else
  echo You are using the correct GNU Autoconf release $CORRECT_AC_VERSION
fi

amversion=`automake --version | head -n 1 | awk '{print $4}'`
if [ ${amversion%.[0-9]} != $CORRECT_AM_VERSION ]; then
  echo You seem to be using GNU Automake version $amversion at `which automake`
  echo You should be using version ${CORRECT_AM_VERSION}, which is the latest release as of September 2014.
  echo The correct version is available on ruby at /opt/local/bin/automake
  echo or /g/melodi/software/bin/as60-amd64 on the MELODI Linux machines.
  exit 1
else
  echo You are using the correct GNU Automake release $CORRECT_AM_VERSION
fi

autoreconf -i || exit 1
echo "Now run \"./configure && make\""
echo "To run the GMTK test suite, \"export TEST_AND_DEV=...\" then \"make check\""

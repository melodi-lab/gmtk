#! /bin/sh

echo PATH is $PATH
echo
if [ -d .hg ]; then
  echo autoconf is `which autoconf`
  autoconf --version
  echo 
  echo automake is `which automake`
  automake --version
  echo
fi
echo gcc is at `which gcc`
gcc --version
echo 
echo g++ is at `which g++`
g++ --version
echo
if [ x$CC != x ]; then
  echo CC is $CC at `which $CC`
  $CC --version
  echo
fi
if [ x$CXX != x ]; then
  echo CXX is $CXX at `which $CXX`
  $CXX --version
  echo
fi
echo wx-config is at `which wx-config`
wx-config --version

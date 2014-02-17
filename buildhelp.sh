#! /bin/sh

echo PATH is $PATH
echo autoconf is `which autoconf`
autoconf --version
echo 
echo automake is `which automake`
automake --version
echo 
echo gcc is `which gcc`
gcc --version
echo 
echo g++ is `which g++`
g++ --version
echo
echo CC is $CC
echo CXX is $CXX 

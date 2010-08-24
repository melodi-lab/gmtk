#!/bin/sh
aclocal
autoconf
autoheader
automake --foreign --add-missing --copy
automake --foreign

#! /bin/sh

set -x

for fd in ts te; do
  ./testFilter -of1 fdiffact1.lst -fmt1  ascii -nf1 11 -ni1 9 \
               -fdiffact1 $fd                                 \
               -of2 fdiffact2.lst -fmt2  ascii -nf2 11 -ni2 9 \
               -fdiffact2 $fd > new.$fd                                \

  ./obs-print  -i1  fdiffact1.lst -ifmt1 ascii -nf1 11 -ni1 9 \
               -fdiffact1 $fd                                 \
               -i2  fdiffact2.lst -ifmt2 ascii -nf2 11 -ni2 9 \
               -fdiffact2 $fd > old.$fd

  cmp {old,new}.$fd || diff {old,new}.$fd > diff.$fd
done

for fd in rl rf se; do
  ./testFilter -of1 fdiffact1.lst -fmt1  ascii -nf1 11 -ni1 9 \
               -of2 fdiffact2.lst -fmt2  ascii -nf2 11 -ni2 9 \
               -fdiffact2 $fd > new.$fd

  ./obs-print  -i1  fdiffact1.lst -ifmt1 ascii -nf1 11 -ni1 9 \
               -i2  fdiffact2.lst -ifmt2 ascii -nf2 11 -ni2 9 \
               -fdiffact2 $fd > old.$fd

  cmp {old,new}.$fd || diff {old,new}.$fd > diff.$fd
done


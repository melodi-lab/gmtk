#!/bin/bash
EXPECTED_ARGS=2
if [ $# -lt $EXPECTED_ARGS ]
then
  echo -e "Usage: `basename $0` <startLine> <endLine> [filename].   \n\
  Prints the lines in range, [<startLine> <endLine>], inclusive. \n\
  If no filename, reads from STDIN. \n\
  If <endLine> < <startLine>, prints all lines until the end of file."
  E_BADARGS=65
  exit $E_BADARGS
fi

COUNT=`echo $2-$1+1 | bc`
tail -n +$1 $3 | head -n $COUNT
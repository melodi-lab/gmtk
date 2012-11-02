#! /bin/sh

# This script rewrites the cardinalities of continuous variables
# to be zero in trifile IDs written by older versions of GMTK that
# had a bug that caused continuous variables to inherit the 
# cardinality of the previous discrete variable in the structure file.

if [ $# != 1 ]; then
  echo Fix non-zero cardinality of continuous variables in trifiles
  echo Usage: $1 trifile
  exit 1
fi
awk '
  BEGIN { in_id = 1 }
  $1=="@@@!!!TRIFILE_END_OF_ID_STRING!!!@@@" { in_id = 0 }
  in_id == 1 { if ($5=="C") $4=0; print }
  in_id == 0 { print }'  $1 > $1.fixed && mv $1.fixed $1 

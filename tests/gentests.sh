#! /bin/sh

# gentests.sh TEST_AND_DEV

if [ $# != 1 ]; then
  echo $0 TEST_AND_DEV
  exit 1
fi

TEST_AND_DEV=$1
TAD=test_and_dev_scripts

for c in `find $TEST_AND_DEV -name tricommand`; do
  ddd=${c%/tricommand}
  dd=${ddd#$TEST_AND_DEV/}
  d=`echo $dd | tr / _`

  if test ! -f $TAD/test_and_dev_${d}.at; then
    echo Generating tests for $dd
    cat <<EOF >$TAD/test_and_dev_${d}.at

# test and dev model $dd

EOF
    for com in tri info train jt vit; do
      if test -x $ddd/${com}command; then
	cat <<EOF >>$TAD/test_and_dev_${d}.at

AT_SETUP([${dd}: ${com}command])
AT_SKIP_IF([test ! -x \$TEST_AND_DEV/$dd/${com}command])
AT_CHECK([(cd \$TEST_AND_DEV/$dd && ulimit -v 8000000 && \
           ulimit -t 120 && ./${com}command)], [], [ignore], [ignore])
AT_CLEANUP
EOF
      fi
    done
  fi
done

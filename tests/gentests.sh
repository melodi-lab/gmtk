#! /bin/sh

# gentests.sh TEST_AND_DEV

if [ $# != 1 ]; then
  echo $0 TEST_AND_DEV
  exit 1
fi

for d in `ls $1`; do
  if test -x $1/$d/tricommand -a ! -f gmtk_test_and_dev_${d}.at; then
    cat <<EOF >gmtk_test_and_dev_${d}.at

# test and dev model $d

EOF
    for com in tri info train jt vit; do
      if test -x $1/$d/${com}command; then
	cat <<EOF >>gmtk_test_and_dev_${d}.at

AT_SETUP([${d}: ${com}command])
AT_SKIP_IF([test ! -x \$TEST_AND_DEV/${d}/${com}command])
AT_CHECK([(cd \$TEST_AND_DEV/${d} && ./${com}command)], [], [ignore])
AT_CLEANUP
EOF
      fi
    done
  fi
done

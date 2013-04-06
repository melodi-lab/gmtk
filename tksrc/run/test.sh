#! /bin/sh

set -x

# $0 output-dir

outdir=$1

#valopts="--tool=massif --stacks=yes"
valopts="--tool=massif --stacks=no"

datadir=/data/rprogers/size-graph/inputData
TKSRC=..
LST=16
BASE=32
SQRT=0.5
ISLOPTS="-base $BASE -lst $LST"


for t in 000001024; do
$TKSRC/gmtkViterbi -strF big.str -inputM big.mtr                   \
    -fileBufferSize 4 -fileWindowSize 1 -allocateDenseCpts 2       \
    -of1 ${datadir}/$t.pfile -constantSpace T               \
    -pVitVal /dev/null                                             \
    -memoryGrowth conservative >test_out.con 

$TKSRC/gmtkViterbi -strF big.str -inputM big.mtr                   \
    -fileBufferSize 4 -fileWindowSize 1 -allocateDenseCpts 2       \
    -of1 ${datadir}/$t.pfile -constantSpace T               \
    -pVitVal /dev/null                                             \
    -memoryGrowth default      >test_out.def 

$TKSRC/gmtkViterbi -strF big.str -inputM big.mtr                   \
    -fileBufferSize 4 -fileWindowSize 1 -allocateDenseCpts 2       \
    -of1 ${datadir}/$t.pfile -constantSpace T               \
    -pVitVal /dev/null                                             \
    -memoryGrowth aggressive   >test_out.agg 
done
exit

for t in  \
000000020 \
000000040 \
000000080 \
000000100 \
000000120 \
000000140 \
000000180 \
000000200 \
000000250 \
000000300 \
000000350 \
000000400 \
000000450 \
000000500 ; do
    (valgrind ${valopts} --massif-out-file=${outdir}/mas.con.$t       \
       $TKSRC/gmtkViterbi -strF big.str -inputM big.mtr               \
       -fileBufferSize 4 -fileWindowSize 1 -allocateDenseCpts 2       \
       -of1 ${datadir}/$t.pfile -constantSpace T -pVitVal /dev/null   \
       -memoryGrowth conservative -mmapViterbiValues F                \
     &> ${outdir}/out.con.$t)&
    (valgrind ${valopts} --massif-out-file=${outdir}/mas.def.$t       \
       $TKSRC/gmtkViterbi -strF big.str -inputM big.mtr               \
       -fileBufferSize 4 -fileWindowSize 1 -allocateDenseCpts 2       \
       -of1 ${datadir}/$t.pfile -constantSpace T -pVitVal /dev/null   \
       -memoryGrowth default -mmapViterbiValues F                     \
     &> ${outdir}/out.def.$t)&
    (valgrind ${valopts} --massif-out-file=${outdir}/mas.agg.$t       \
       $TKSRC/gmtkViterbi -strF big.str -inputM big.mtr               \
       -fileBufferSize 4 -fileWindowSize 1 -allocateDenseCpts 2       \
       -of1 ${datadir}/$t.pfile -constantSpace T -pVitVal /dev/null   \
       -memoryGrowth aggressive    -mmapViterbiValues F               \
     &> ${outdir}/out.agg.$t)&
wait
done

exit

\
000001024 \
000002048 \
000004096 \
000020480 

000001024 \
000002048 \
000004096 \
000006144 \
000008192 \
000010240 \
000012288 \
000014336 \
000016384 \
000018432 \
000020480 \
000022528 \
000024576 \
000026624 \
000028672 \
000030720 \
000032768 \
000034816 \
000036864 \
000038912 \
000040960 \
000051200 \
000065536 \
000080000 \
000100000 \
000120000 \
000140000 \
000160000 \
000180000 \
000200000 \


\
000040960 \
000051200 \
000065536 \
000080000 \
000100000 \
000120000 \
000140000 \
000160000 

000200000 \
002200000 \
002400000 \
003000000 \
003600000 \
010000000 \
016000000


002200000 \
002400000 \
003000000 \
003600000
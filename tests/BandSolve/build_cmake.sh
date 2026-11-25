BDIR=build
PROG=TestBandSolve
if [ ! -d ${BDIR} ]
then
  mkdir ${BDIR}
fi

cd ${BDIR}
NP=4
cmake ../
make -j ${NP}
cp ${PROG} ../
cd ../

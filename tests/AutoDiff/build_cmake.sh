BDIR=build
PROG=TestAutoDiff
if [ ! -d ${BDIR} ]
then
  mkdir ${BDIR}
fi

cd ${BDIR}
NP=1
cmake ../
make -j ${NP}
cp ${PROG} ../
cd ../

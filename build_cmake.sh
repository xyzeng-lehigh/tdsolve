BDIR=build
PROG=TumorDynRun
if [ ! -d ${BDIR} ]
then
  mkdir ${BDIR}
fi

cd ${BDIR}
cmake ../
make
cp ${PROG} ../
cd ../

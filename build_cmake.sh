BDIR=build
PROG=TumorDynRun
PROGPAR=ParAnalRun
if [ ! -d ${BDIR} ]
then
  mkdir ${BDIR}
fi

cd ${BDIR}
cmake ../
make
cp ${PROG} ../
cp ${PROGPAR} ../
cd ../

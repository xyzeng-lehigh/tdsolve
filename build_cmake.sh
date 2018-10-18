BDIR=build
PROG=TumorDynRun
PROGPAR=ParAnalRun
if [ ! -d ${BDIR} ]
then
  mkdir ${BDIR}
fi

cd ${BDIR}
cmake -DCMAKE_C_COMPILER=$HOME/bin/gcc -DCMAKE_CXX_COMPILER=$HOME/bin/g++ ../
make -j 3
cp ${PROG} ../
cp ${PROGPAR} ../
cd ../

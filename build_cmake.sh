BDIR=build
PROG=TumorDynRun
PROGPAR=ParAnalRun
if [ ! -d ${BDIR} ]
then
  mkdir ${BDIR}
fi

cd ${BDIR}
NP=1
case ${OSTYPE//[0-9.]/} in
  darwin)
    #cmake -DCMAKE_C_COMPILER=$HOME/bin/gcc -DCMAKE_CXX_COMPILER=$HOME/bin/g++ ../
    #$TDS/support/cmake/bin/cmake -DCMAKE_C_COMPILER=$HOME/bin/gcc -DCMAKE_CXX_COMPILER=$HOME/bin/g++ ../
    cmake ../
    NP=2
    ;;
  linux-gnu)
    cmake ../
    NP=8
    ;;
esac
make -j ${NP}
cp ${PROG} ../
#cp ${PROGPAR} ../
cd ../

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
    export 
    cmake \
      -DCMAKE_C_COMPILER=$CSUPPROOT/gcc-9.3.0/gcc/bin/gcc \
      -DCMAKE_CXX_COMPILER=$CSUPPROOT/gcc-9.3.0/gcc/bin/g++ \
      -DTrilinos_DIR:PATH=$CSUPPROOT/gcc-9.3.0/trilinos-12.7 \
      ../
      #-DTrilinos_DIR:PATH= \
    #$TDS/support/cmake/bin/cmake -DCMAKE_C_COMPILER=$HOME/bin/gcc -DCMAKE_CXX_COMPILER=$HOME/bin/g++ ../
    #cmake ../
    NP=2
    ;;
  linux-gnu)
    cmake3 -DCMAKE_C_COMPILER=$CSUPPROOT/gcc-10.1.0/gcc/bin/gcc -DCMAKE_CXX_COMPILER=$CSUPPROOT/gcc-10.1.0/gcc/bin/g++ ../
    NP=8
    ;;
esac
make -j ${NP}
cp ${PROG} ../
cp ${PROGPAR} ../
cd ../

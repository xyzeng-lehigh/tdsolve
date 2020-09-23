BDIR=build
PROG=TestSacadoDup
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
    ;;
  linux-gnu)
    cmake3 \
      -DCMAKE_C_COMPILER=$CSUPPROOT/gcc-9.3.0/gcc/bin/gcc \
      -DCMAKE_CXX_COMPILER=$CSUPPROOT/gcc-9.3.0/gcc/bin/g++ \
      -DTrilinos_DIR:PATH=$CSUPPROOT/gcc-9.3.0/trilinos-12.7 \
      ../
    ;;
esac
make -j ${NP}
cp ${PROG} ../
cd ../

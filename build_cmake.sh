BDIR=build
PROG=TumorDynRun
PROGPAR=ParAnalRun
if [ ! -d ${BDIR} ]
then
  mkdir ${BDIR}
fi

cd ${BDIR}
NP=2
case ${OSTYPE//[0-9.]/} in
  darwin)
    export 
    /opt/homebrew/bin/cmake \
      -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-11 \
      -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-11 \
      -DTrilinos_DIR:PATH=$CSUPPROOT/gcc-10.2.1/trilinos-13.0.1 \
      -DCMAKE_CXX_EXTENSIONS=Off \
      ../
    NP=4
    ;;
  linux-gnu)
    cmake3 \
      -DCMAKE_C_COMPILER=$CSUPPROOT/gcc-9.3.0/gcc/bin/gcc \
      -DCMAKE_CXX_COMPILER=$CSUPPROOT/gcc-9.3.0/gcc/bin/g++ \
      -DTrilinos_DIR:PATH=$CSUPPROOT/gcc-9.3.0/trilinos-12.7 \
      ../
    NP=8
    ;;
esac
make -j ${NP}
cp ${PROG} ../
cp ${PROGPAR} ../
cd ../

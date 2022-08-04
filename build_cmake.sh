BDIR=build
PROG=TumorDynRun
PROGPAR=ParAnalRun
if [ ! -d ${BDIR} ]
then
  mkdir ${BDIR}
fi

cd ${BDIR}
NP=4
case ${OSTYPE//[0-9.]/} in
  darwin)
    export 
    /opt/homebrew/bin/cmake \
      -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-12 \
      -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-12 \
      -DTrilinos_DIR:PATH=$CSUPPROOT/gcc-12/trilinos-13.4.0 \
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
cd ../
if [ ! -f ${PROG} ];
then
  ln -s ${BDIR}/${PROG} ./
fi
#if [ ! -f ${PROGPAR} ];
#then
#  ln -s ${BDIR}/${PROGPAR} ./
#fi
#cp ${PROG} ../
#cp ${PROGPAR} ../
cd ../

# Trilinos is no longer needed. Extra arguments are passed on to cmake, e.g.
# to use the Sacado package of Trilinos for automatic differentiation instead
# of the built-in one:
#   sh build_cmake.sh -DUSE_TRILINOS=ON -DTrilinos_DIR:PATH=$CSUPPROOT/$TDSUPPNAME/trilinos
BDIR=build
PROG=TumorDynRun
PROGPAR=ParAnalRun
if [ ! -d ${BDIR} ]
then
  mkdir ${BDIR}
fi

cd ${BDIR}
NP=8
case ${OSTYPE//[0-9.]/} in
  darwin)
    export 
    /opt/homebrew/bin/cmake \
      -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-12 \
      -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-12 \
      -DCMAKE_CXX_EXTENSIONS=Off \
      "$@" ../
    NP=4
    ;;
  linux-gnu)
    cmake "$@" ../
    #cmake3 \
    #  -DCMAKE_C_COMPILER=$CSUPPROOT/gcc-9.3.0/gcc/bin/gcc \
    #  -DCMAKE_CXX_COMPILER=$CSUPPROOT/gcc-9.3.0/gcc/bin/g++ \
    #  -DTrilinos_DIR:PATH=$CSUPPROOT/gcc-9.3.0/trilinos-12.7 \
    #  -DCMAKE_PREFIX_PATH=$CSUPPROOT/gcc-$HCLGCCVER/openblas \
    #  ../
    NP=8
    ;;
esac
make -j ${NP}
cd ../
if [ ! -f ${PROG} ];
then
  ln -s ${BDIR}/${PROG} ./
fi
if [ ! -f ${PROGPAR} ];
then
  ln -s ${BDIR}/${PROGPAR} ./
fi
#cp ${PROG} ../
#cp ${PROGPAR} ../
cd ../

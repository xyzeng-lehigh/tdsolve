#ifndef COMPUT_CHOICES_H
#define COMPUT_CHOICES_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
using namespace std;

#ifndef MAX_OPTION_STR_LENGTH
#define MAX_OPTION_STR_LENGTH 40
#endif

// Part fifteen: choice of Quadrature object type
enum ChoiceQuadRule {QUADRULE_OFF = -1, QUADRULE_GAUSSIAN_LEGENDRE = 0, QUADRULE_GAUSSIAN_LOBATTO = 1, QUADRULE_CLENSHAW_CURTIS = 2};
enum ChoiceQuadType {QUADTYPE_OFF = -1, QUADTYPE_EDGE = 0, QUADTYPE_TRIA = 1, QUADTYPE_RECT = 2, QUADTYPE_TETRA = 3, QUADTYPE_CUBOID = 4};
string ChoiceQuadRule2Str(ChoiceQuadRule quadRule);
string ChoiceQuadType2Str(ChoiceQuadType quadType);

// Part eighteen: linear solver
enum ChoiceLSysMat {LSYS_MAT_FULL=0, LSYS_MAT_SPARSE=1};
enum ChoiceLSysMatInfo {LSYS_MAT_INFO_NORMAL=0, LSYS_MAT_INFO_SYM=1, LSYS_MAT_INFO_SPD=2, LSYS_MAT_INFO_M=3};
// The class of projection method, like SD (steepest descent), MR (minimal residual) are not explored
enum ChoiceLSysType {LSYS_TYPE_DIR=0, LSYS_TYPE_ITER=1, LSYS_TYPE_KRYLOV=2};
// Direct solver of the Sacado framework: LU - banded LU factorization after renumbering (Alg_BandSolve.h);
//                                        INV - inverse of the full matrix with largest entry pivoting (HMAT::inv)
enum ChoiceLSysDirect {LSYS_DIRECT_LU=0, LSYS_DIRECT_INV=1};
enum ChoiceLSysIter {LSYS_ITER_JAC=0, LSYS_ITER_GS=1, LSYS_ITER_BGS=2, LSYS_ITER_SOR=3, LSYS_ITER_BSOR=4};
enum ChoiceLSysKrylov {LSYS_KRYLOV_CG=0, LSYS_KRYLOV_GMRES=1, LSYS_KRYLOV_LANCZOS=2, LSYS_KRYLOV_QMR=3};
enum ChoiceLSysPrec {LSYS_PREC_OFF=-1, LSYS_PREC_LINK=0, LSYS_PREC_EYE=1, LSYS_PREC_ILU=2, LSYS_PREC_BILU=3};
string ChoiceLSysMat2Str(ChoiceLSysMat lsysMat);
string ChoiceLSysMatInfo2Str(ChoiceLSysMatInfo lsysMatInfo);
string ChoiceLSysType2Str(ChoiceLSysType lsysType); 
string ChoiceLSysDirect2Str(ChoiceLSysDirect lsysDirect);
string ChoiceLSysIter2Str(ChoiceLSysIter lsysIter);
string ChoiceLSysKrylov2Str(ChoiceLSysKrylov lsysKrolov);
string ChoiceLSysPrec2Str(ChoiceLSysPrec lsysPrec);

#endif

#include "Comput_Choices.h"
#include <cstring>
using namespace std;

string ChoiceQuadRule2Str(ChoiceQuadRule quadRule) {
	char quadRule_str[MAX_OPTION_STR_LENGTH];
	switch (quadRule) {
		case QUADRULE_OFF:
			strcpy(quadRule_str,"Off");
			break;
		case QUADRULE_GAUSSIAN_LEGENDRE:
			strcpy(quadRule_str,"GaussianLegendre");
			break;
		case QUADRULE_GAUSSIAN_LOBATTO:
			strcpy(quadRule_str,"GaussianLobatto");
			break;
		case QUADRULE_CLENSHAW_CURTIS:
			strcpy(quadRule_str,"ClenshawCurtis");
			break;
		default:
			strcpy(quadRule_str,"");
			break;
	}
	return quadRule_str;
}

string ChoiceQuadType2Str(ChoiceQuadType quadType) {
	char quadType_str[MAX_OPTION_STR_LENGTH];
	switch (quadType) {
		case QUADTYPE_OFF:
			strcpy(quadType_str,"Off");
			break;
		case QUADTYPE_EDGE:
			strcpy(quadType_str,"Edge");
			break;
		case QUADTYPE_TRIA:
			strcpy(quadType_str,"Triangle");
			break;
		case QUADTYPE_RECT:
			strcpy(quadType_str,"Rectangle");
			break;
		case QUADTYPE_TETRA:
			strcpy(quadType_str,"Tetrahedron");
			break;
		case QUADTYPE_CUBOID:
			strcpy(quadType_str,"Cuboid");
			break;
		default:
			strcpy(quadType_str,"");
			break;
	}
	return quadType_str;
}

// Part eighteen: linear solver

string ChoiceLSysMat2Str(ChoiceLSysMat lsysMat) {
    char lsysMat_str[MAX_OPTION_STR_LENGTH];
    switch (lsysMat) {
        case LSYS_MAT_FULL:
            strcpy(lsysMat_str,"FullMatrix");
            break;
        case LSYS_MAT_SPARSE:
            strcpy(lsysMat_str,"SparseMatrix");
            break;
        default:
            strcpy(lsysMat_str,"");
            break;
    }
    return lsysMat_str;
}

string ChoiceLSysMatInfo2Str(ChoiceLSysMatInfo lsysMatInfo) {
    char lsysMatInfo_str[MAX_OPTION_STR_LENGTH];
    switch (lsysMatInfo) {
        case LSYS_MAT_INFO_NORMAL:
            strcpy(lsysMatInfo_str,"NormalMatrix");
            break;
        case LSYS_MAT_INFO_SYM:
            strcpy(lsysMatInfo_str,"SymmetricMatrix");
            break;
        case LSYS_MAT_INFO_SPD:
            strcpy(lsysMatInfo_str,"SymmetricPositiveDefiniteMatrix");
            break;
        case LSYS_MAT_INFO_M:
            strcpy(lsysMatInfo_str,"MMatrix");
            break;
        default:
            strcpy(lsysMatInfo_str,"");
            break;
    }
    return lsysMatInfo_str;
}


string ChoiceLSysType2Str(ChoiceLSysType lsysType) {
    char lsysType_str[MAX_OPTION_STR_LENGTH];
    switch (lsysType) {
        case LSYS_TYPE_DIR:
            strcpy(lsysType_str,"Direct");
            break;
        case LSYS_TYPE_ITER:
            strcpy(lsysType_str,"Iterative");
            break;
        case LSYS_TYPE_KRYLOV:
            strcpy(lsysType_str,"KrylovSpace");
            break;
        default:
            strcpy(lsysType_str,"");
            break;
    }
    return lsysType_str;
}

string ChoiceLSysIter2Str(ChoiceLSysIter lsysIter) {
    char lsysIter_str[MAX_OPTION_STR_LENGTH];
    switch (lsysIter) {
        case LSYS_ITER_JAC:
            strcpy(lsysIter_str,"Jacobi");
            break;
        case LSYS_ITER_GS:
            strcpy(lsysIter_str,"GaussSeidel");
            break;
        case LSYS_ITER_BGS:
            strcpy(lsysIter_str,"BackwardGaussSeidel");
            break;
        case LSYS_ITER_SOR:
            strcpy(lsysIter_str,"SOR");
            break;
        case LSYS_ITER_BSOR:
            strcpy(lsysIter_str,"BackwardSOR");
            break;
        default:
            strcpy(lsysIter_str,"");
            break;
    }
    return lsysIter_str;
}

string ChoiceLSysKrylov2Str(ChoiceLSysKrylov lsysKrylov) {
    char lsysKrylov_str[MAX_OPTION_STR_LENGTH];
    switch (lsysKrylov) {
        case LSYS_KRYLOV_CG:
            strcpy(lsysKrylov_str,"ConjugateGradient");
            break;
        case LSYS_KRYLOV_GMRES:
            strcpy(lsysKrylov_str,"GMRES");
            break;
        case LSYS_KRYLOV_LANCZOS:
            strcpy(lsysKrylov_str,"Lanczos");
            break;
        case LSYS_KRYLOV_QMR:
            strcpy(lsysKrylov_str,"QuasiMinimalResidual");
            break;
        default:
            strcpy(lsysKrylov_str,"");
            break;
    }
    return lsysKrylov_str;
}

string ChoiceLSysPrec2Str(ChoiceLSysPrec lsysPrec) {
    char lsysPrec_str[MAX_OPTION_STR_LENGTH];
    switch (lsysPrec) {
        case LSYS_PREC_OFF:
            strcpy(lsysPrec_str,"Off");
            break;
        case LSYS_PREC_LINK:
            strcpy(lsysPrec_str,"LinkFromOutside");
            break;
        case LSYS_PREC_EYE:
            strcpy(lsysPrec_str,"Identity");
            break;
        case LSYS_PREC_ILU:
            strcpy(lsysPrec_str,"IncompleteLU");
            break;
        case LSYS_PREC_BILU:
            strcpy(lsysPrec_str,"BlockIncompleteLU");
            break;
        default:
            strcpy(lsysPrec_str,"");
            break;
    }
    return lsysPrec_str;
}


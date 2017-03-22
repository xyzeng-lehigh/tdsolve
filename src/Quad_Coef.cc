#include <cassert>
#include "Quad_Coef.h"
#include "HighDim_Alg.h"

QuadCoef* newQuadCoef(ChoiceQuadRule quadRule, 
		ChoiceQuadType quadType, int order) {
	QuadCoef*		quadCoef = NULL;
	QuadCoefGLeg* 	quadCoefGLeg;
	QuadCoefGLob*	quadCoefGLob;
	QuadCoefCC*		quadCoefCC;
	switch (quadRule) {
		case QUADRULE_GAUSSIAN_LEGENDRE:
			quadCoefGLeg = new QuadCoefGLeg(quadType,order);
			quadCoef = &(*quadCoefGLeg);
			break;
		case QUADRULE_GAUSSIAN_LOBATTO:
			quadCoefGLob = new QuadCoefGLob(quadType,order);
			quadCoef = &(*quadCoefGLob);
			break;
		case QUADRULE_CLENSHAW_CURTIS:
			quadCoefCC = new QuadCoefCC(quadType,order);
			quadCoef = &(*quadCoefCC);
			break;
        default:
            break;
	}
	return quadCoef;
}

/*
 *	Base class
 */
QuadCoef::QuadCoef(ChoiceQuadType quadType) {
	switch (quadType) {
		case QUADTYPE_EDGE:
			type = EDGE;
			dim = 1;
			break;
		case QUADTYPE_TRIA:
			type = TRIA;
			dim = 3;
			break;
		case QUADTYPE_RECT:
			type = RECT;
			dim = 2;
			break;
		case QUADTYPE_TETRA:
			type = TETRA;
			dim = 4;
			break;
		case QUADTYPE_CUBOID:
			type = CUBOID;
			dim = 3;
			break;
        default:
            break;
	}
}

QuadCoef::~QuadCoef(void) {
	if (W!=NULL) delete[] W;
	if (Q!=NULL) {
		for (int n=0; n<N; n++)
			if (Q[n]!=NULL) delete[] Q[n];
		delete[] Q;
	}
}

void QuadCoef::setup(void) {
	assert(order>0);
	//	Compute N and Np according to order:
	//	Np+1 >= order
	switch (type) {
		case EDGE:
			this->setup_edge();
			break;
		case TRIA:
			this->setup_tria();
			break;
		case RECT:
			this->setup_rect();
			break;
		case TETRA:
			this->setup_tetra();
			break;
		case CUBOID:
			this->setup_cuboid();
			break;
	}
	W = new DType[N];
	Q = new DType*[N];
	for (int n=0; n<N; n++) Q[n] = new DType[dim];
	for (int n=0; n<N; n++) W[n] = 0.0;
	for (int n=0; n<N; n++) for (int d=0; d<dim; d++) Q[n][d] = 0.0;
	bool support = this->checkType();
	if (!support)
		fprintf(stderr,"Error: geometric type not supported for %s quadrature rules.\n",(this->toString()).c_str());
	else
		switch (type) {
			case EDGE:
				this->compute_edge();
				break;
			case TRIA:
				this->compute_tria();
				break;
			case RECT:
				this->compute_rect();
				break;
			case TETRA:
				this->compute_tetra();
				break;
			case CUBOID:
				this->compute_cuboid();
				break;
		}
}

/*
 *	Gaussian-Legendre
 */

QuadCoefGLeg::~QuadCoefGLeg(void) {
	// Do nothing
}

string QuadCoefGLeg::toString(void) {
	char classNameString[40];
	strcpy(classNameString,"Gaussian-Legendre");
	return classNameString;
}

void QuadCoefGLeg::setup_edge(void) {
	N = (order+1)/2;		// order <= Np+1 = 2N
	Np = 2*N-1;
	return;
}

void QuadCoefGLeg::setup_tria(void) {
	switch (order) {
		case 1:
		case 2:
			N = 1;
			Np = 1;
			break;
		default:
			N = 1;
			Np = 1;
			fprintf(stderr,"Warning: Order %d not supported for TRIA in quadrature rule %s.\n",order,toString().c_str());
			break;
	}
	return;
}

void QuadCoefGLeg::setup_rect(void) {
	int N1d = (order+1)/2;	// order <= Np+1 = 2*N1d
	Np = 2*N1d-1;
	N  = N1d*N1d;
	return;
}

void QuadCoefGLeg::setup_tetra(void) {
	switch (order) {
		case 1:
		case 2:
			N = 1;
			Np = 1;
			break;
		default:
			N = 1;
			Np = 1;
			fprintf(stderr,"Warning: Order %d not supported for TETRA in quadrature rule %s.\n",order,toString().c_str());
			break;
	}
	return;
}

void QuadCoefGLeg::setup_cuboid(void) {
	int N1d = (order+1)/2;	// order <= Np+1 = 2*N1d
	Np = 2*N1d-1;
	N  = N1d*N1d*N1d;
	return;
}

void QuadCoefGLeg::compute_edge(void) {
	switch (N) {
		case 1:
			Q[0][0] = 0.0;				W[0] = 1.0;
			break;
		case 2:
			Q[0][0] = -1.0/sqrt(3.0);	W[0] = 0.5;
			Q[1][0] = 1.0/sqrt(3.0);	W[1] = 0.5;
			break;
		case 3:
			Q[0][0] = -sqrt(3.0/5.0);	W[0] = 5.0/18.0;
			Q[1][0] = 0.0;				W[1] = 4.0/9.0;
			Q[2][0] = sqrt(3.0/5.0);	W[2] = 5.0/18.0;
			break;
		case 4:
			Q[0][0] = -sqrt((3.0+2.0*sqrt(1.2))/7.0);
			Q[1][0] = -sqrt((3.0-2.0*sqrt(1.2))/7.0);
			Q[2][0] = sqrt((3.0-2.0*sqrt(1.2))/7.0);
			Q[3][0] = sqrt((3.0+2.0*sqrt(1.2))/7.0);
			W[0] = 0.25-sqrt(30.0)/72.0;
			W[1] = 0.25+sqrt(30.0)/72.0;
			W[2] = 0.25+sqrt(30.0)/72.0;
			W[3] = 0.25-sqrt(30.0)/72.0;
			break;
		case 5:
			Q[0][0] = -sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0;
			Q[1][0] = -sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0;
			Q[2][0] = 0.0;
			Q[3][0] = sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0;
			Q[4][0] = sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0;
			W[0] = (322.0-13.0*sqrt(70.0))/1800.0;
			W[1] = (322.0+13.0*sqrt(70.0))/1800.0;
			W[2] = 64.0/225.0;
			W[3] = (322.0+13.0*sqrt(70.0))/1800.0;
			W[4] = (322.0-13.0*sqrt(70.0))/1800.0;
			break;
		default:
			fprintf(stderr,"Error: Quadrature rule of %s for %d nodes in the case of EDGE not implemented.\n",toString().c_str(),N);
			break;
	}
	return;
}

void QuadCoefGLeg::compute_tria(void) {
	switch (N) {
		case 1:
			Q[0][0] = 0.0;
			Q[0][1] = 0.0;
			Q[0][2] = 0.0;
			W[0] 	= 1.0;
			break;
		default:
			fprintf(stderr,"Error: Quadrature rule of %s for %d nodes in the case of TRIA not implemented.\n",toString().c_str(),N);
			break;
	}
	return;
}

void QuadCoefGLeg::compute_rect(void) {
	DType* Q1D = NULL;
	DType* W1D = NULL;
	int N1d = 0;
	switch (N) {
		case 1:
			N1d = 1;
			Q1D = new DType[N1d];	W1D = new DType[N1d];
			Q1D[0] = 0.0;				W1D[0] = 1.0;
			break;
		case 4:
			N1d = 2;
			Q1D = new DType[N1d];	W1D = new DType[N1d];
			Q1D[0] = -1.0/sqrt(3.0);	W1D[0] = 0.5;
			Q1D[1] = 1.0/sqrt(3.0);		W1D[1] = 0.5;
			break;
		case 9:
			N1d = 3;
			Q1D = new DType[N1d];	W1D = new DType[N1d];
			Q1D[0] = -sqrt(3.0/5.0);		W1D[0] = 5.0/18.0;
			Q1D[1] = 0.0;				W1D[1] = 4.0/9.0;
			Q1D[2] = sqrt(3.0/5.0);		W1D[2] = 5.0/18.0;
			break;
		case 16:
			N1d = 4;
			Q1D = new DType[N1d];	W1D = new DType[N1d];
			Q1D[0] = -sqrt((3.0+2.0*sqrt(1.2))/7.0);
			Q1D[1] = -sqrt((3.0-2.0*sqrt(1.2))/7.0);
			Q1D[2] = sqrt((3.0-2.0*sqrt(1.2))/7.0);
			Q1D[3] = sqrt((3.0+2.0*sqrt(1.2))/7.0);
			W1D[0] = 0.25-sqrt(30.0)/72.0;
			W1D[1] = 0.25+sqrt(30.0)/72.0;
			W1D[2] = 0.25+sqrt(30.0)/72.0;
			W1D[3] = 0.25-sqrt(30.0)/72.0;
			break;
		case 25:
			N1d = 5;
			Q1D = new DType[N1d];	W1D = new DType[N1d];
			Q1D[0] = -sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0;
			Q1D[1] = -sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0;
			Q1D[2] = 0.0;
			Q1D[3] = sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0;
			Q1D[4] = sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0;
			W1D[0] = (322.0-13.0*sqrt(70.0))/1800.0;
			W1D[1] = (322.0+13.0*sqrt(70.0))/1800.0;
			W1D[2] = 64.0/225.0;
			W1D[3] = (322.0+13.0*sqrt(70.0))/1800.0;
			W1D[4] = (322.0-13.0*sqrt(70.0))/1800.0;
			break;
		default:
			fprintf(stderr,"Error: Quadrature rule of %s for %d nodes in the case of RECT not implemented.\n",toString().c_str(),N);
			break;
	}
	for (int k1 = 0; k1<N1d; k1++)
		for (int k2 = 0; k2<N1d; k2++) {
			int n = k1*N1d+k2;
			Q[n][0] = Q1D[k1];
			Q[n][1] = Q1D[k2];
			W[n] = W1D[k1]*W1D[k2];
		}
	if (Q1D!=NULL) delete[] Q1D;
	if (W1D!=NULL) delete[] W1D;
	return;
}

void QuadCoefGLeg::compute_tetra(void) {
	switch (N) {
		case 1:
			Q[0][0] = 0.0;
			Q[0][1] = 0.0;
			Q[0][2] = 0.0;
			Q[0][3] = 0.0;
			W[0] 	= 1.0;
			break;
		default:
			fprintf(stderr,"Error: Quadrature rule of %s for %d nodes in the case of TETRA not implemented.\n",toString().c_str(),N);
			break;
	}
	return;
}

void QuadCoefGLeg::compute_cuboid(void) {
	DType* Q1D = NULL;
	DType* W1D = NULL;
	int N1d = 0;
	switch (N) {
		case 1:
			N1d = 1;
			Q1D = new DType[N1d];	W1D = new DType[N1d];
			Q1D[0] = 0.0;				W1D[0] = 1.0;
			break;
		case 8:
			N1d = 2;
			Q1D = new DType[N1d];	W1D = new DType[N1d];
			Q1D[0] = -1.0/sqrt(3.0);	W1D[0] = 0.5;
			Q1D[1] = 1.0/sqrt(3.0);		W1D[1] = 0.5;
			break;
		case 27:
			N1d = 3;
			Q1D = new DType[N1d];	W1D = new DType[N1d];
			Q1D[0] = -sqrt(3.0/5.0);		W1D[0] = 5.0/18.0;
			Q1D[1] = 0.0;				W1D[1] = 4.0/9.0;
			Q1D[2] = sqrt(3.0/5.0);		W1D[2] = 5.0/18.0;
			break;
		case 64:
			N1d = 4;
			Q1D = new DType[N1d];	W1D = new DType[N1d];
			Q1D[0] = -sqrt((3.0+2.0*sqrt(1.2))/7.0);
			Q1D[1] = -sqrt((3.0-2.0*sqrt(1.2))/7.0);
			Q1D[2] = sqrt((3.0-2.0*sqrt(1.2))/7.0);
			Q1D[3] = sqrt((3.0+2.0*sqrt(1.2))/7.0);
			W1D[0] = 0.25-sqrt(30.0)/72.0;
			W1D[1] = 0.25+sqrt(30.0)/72.0;
			W1D[2] = 0.25+sqrt(30.0)/72.0;
			W1D[3] = 0.25-sqrt(30.0)/72.0;
			break;
		case 125:
			N1d = 5;
			Q1D = new DType[N1d];	W1D = new DType[N1d];
			Q1D[0] = -sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0;
			Q1D[1] = -sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0;
			Q1D[2] = 0.0;
			Q1D[3] = sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0;
			Q1D[4] = sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0;
			W1D[0] = (322.0-13.0*sqrt(70.0))/1800.0;
			W1D[1] = (322.0+13.0*sqrt(70.0))/1800.0;
			W1D[2] = 64.0/225.0;
			W1D[3] = (322.0+13.0*sqrt(70.0))/1800.0;
			W1D[4] = (322.0-13.0*sqrt(70.0))/1800.0;
			break;
		default:
			fprintf(stderr,"Error: Quadrature rule of %s for %d nodes in the case of CUBOID not implemented.\n",toString().c_str(),N);
			break;
	}
	for (int k1 = 0; k1<N1d; k1++)
		for (int k2 = 0; k2<N1d; k2++)
			for (int k3 = 0; k3<N1d; k3++) {
				int n = (k1*N1d+k2)*N1d+k3;
				Q[n][0] = Q1D[k1];
				Q[n][1] = Q1D[k2];
				Q[n][2] = Q1D[k3];
				W[n] = W1D[k1]*W1D[k2]*W1D[k3];
			}
	if (Q1D!=NULL) delete[] Q1D;
	if (W1D!=NULL) delete[] W1D;
	return;
}

bool QuadCoefGLeg::checkType(void) {
	bool support;
	switch (type) {
		case EDGE:
		case RECT:
		case CUBOID:
		case TRIA:
		case TETRA:
			support = true;
			break;
		default:
			support = false;
			break;
	}
	return support;
}

/*
 *	Gaussian-Lobatto
 */

QuadCoefGLob::~QuadCoefGLob(void) {
	// Do nothing
}

string QuadCoefGLob::toString(void) {
	char classNameString[40];
	strcpy(classNameString,"Gaussian-Lobatto");
	return classNameString;
}

void QuadCoefGLob::setup_edge(void) {
	N = (order+3)/2;		// order <= Np+1 = 2N-2
	Np = 2*N-3;
	return;
}

void QuadCoefGLob::setup_tria(void) {
	switch (order) {
		case 1:
		case 2:
			N = 3;
			Np = 1;
			break;
        case 3:
            N = 4;
            Np = 2;
            break;
		default:
			N = 3;
			Np = 1;
			fprintf(stderr,"Warning: Order %d not supported for TRIA in quadrature rule %s.\n",order,toString().c_str());
			break;
	}
	return;
}

void QuadCoefGLob::setup_rect(void) {
	int N1d = (order+3)/2;	// order <= Np+1 = 2N-2
	Np = 2*N1d-3;
	N  = N1d*N1d;
	return;
}

void QuadCoefGLob::setup_tetra(void) {
	switch (order) {
		case 1:
		case 2:
			N = 4;
			Np = 1;
			break;
        case 3:
            N = 5;
            Np = 2;
            break;
		default:
			N = 4;
			Np = 1;
			fprintf(stderr,"Warning: Order %d not supported for TETRA in quadrature rule %s.\n",order,toString().c_str());
			break;
	}
	return;
}

void QuadCoefGLob::setup_cuboid(void) {
	int N1d = (order+3)/2;	// order <= Np+1 = 2N-2
	Np = 2*N1d-3;
	N  = N1d*N1d*N1d;
	return;
}

void QuadCoefGLob::compute_edge(void) {
	switch (N) {
		case 2:
			Q[0][0] = -1.0;				W[0] = 0.5;
			Q[1][0] = 1.0;				W[1] = 0.5;
			break;
		case 3:
			Q[0][0] = -1.0;				W[0] = 1.0/6.0;
			Q[1][0] = 0.0;				W[1] = 2.0/3.0;
			Q[2][0] = 1.0;				W[2] = 1.0/6.0;
			break;
		case 4:
			Q[0][0] = -1.0;				W[0] = 1.0/12.0;
			Q[1][0] = -sqrt(0.2);		W[1] = 5.0/12.0;
			Q[2][0] = sqrt(0.2);		W[2] = 5.0/12.0;
			Q[3][0] = 1.0;				W[3] = 1.0/12.0;
			break;
		case 5:
			Q[0][0] = -1.0;				W[0] = 0.05;
			Q[1][0] = -sqrt(3.0/7.0);	W[1] = 49.0/180.0;
			Q[2][0] = 0.0;				W[2] = 16.0/45.0;
			Q[3][0] = sqrt(3.0/7.0);	W[3] = 49.0/180.0;
			Q[4][0] = 1.0;				W[4] = 0.05;
			break;
		default:
			fprintf(stderr,"Error: Quadrature rule of %s for %d nodes in the case of EDGE not implemented.\n",toString().c_str(),N);
			break;
	}
	return;
}

void QuadCoefGLob::compute_tria(void) {
	switch (N) {
		case 3:
			Q[0][0] = 1.0; Q[0][1] = 0.0; Q[0][2] = 0.0; W[0] = 1.0/3.0;
			Q[1][0] = 0.0; Q[1][1] = 1.0; Q[1][2] = 0.0; W[1] = 1.0/3.0;
			Q[2][0] = 0.0; Q[2][1] = 0.0; Q[2][2] = 1.0; W[2] = 1.0/3.0;
			break;
        case 4:
			Q[0][0] = 1.0; Q[0][1] = 0.0; Q[0][2] = 0.0; W[0] = 1.0/12.0;
			Q[1][0] = 0.0; Q[1][1] = 1.0; Q[1][2] = 0.0; W[1] = 1.0/12.0;
			Q[2][0] = 0.0; Q[2][1] = 0.0; Q[2][2] = 1.0; W[2] = 1.0/12.0;
            Q[3][0] = 0.0; Q[3][1] = 0.0; Q[3][2] = 0.0; W[3] = 3.0/4.0;
            break;
		default:
			fprintf(stderr,"Error: Quadrature rule of %s for %d nodes in the case of TRIA not implemented.\n",toString().c_str(),N);
			break;
	}
	return;
}

void QuadCoefGLob::compute_rect(void) {
	DType* Q1D = NULL;
	DType* W1D = NULL;
	int N1d = 0;
	switch (N) {
		case 4:
			N1d = 2;
			Q1D = new DType[N1d];		W1D = new DType[N1d];
			Q1D[0] = -1.0;				W1D[0] = 0.5;
			Q1D[1] = 1.0;				W1D[1] = 0.5;
			break;
		case 9:
			N1d = 3;
			Q1D = new DType[N1d];		W1D = new DType[N1d];
			Q1D[0] = -1.0;				W1D[0] = 1.0/6.0;
			Q1D[1] = 0.0;				W1D[1] = 2.0/3.0;
			Q1D[2] = 1.0;				W1D[2] = 1.0/6.0;
			break;
		case 16:
			N1d = 4;
			Q1D = new DType[N1d];		W1D = new DType[N1d];
			Q1D[0] = -1.0;				W1D[0] = 1.0/12.0;
			Q1D[1] = -sqrt(0.2);		W1D[1] = 5.0/12.0;
			Q1D[2] = sqrt(0.2);			W1D[2] = 5.0/12.0;
			Q1D[3] = 1.0;				W1D[3] = 1.0/12.0;
			break;
		case 25:
			N1d = 5;
			Q1D = new DType[N1d];		W1D = new DType[N1d];
			Q1D[0] = -1.0;				W1D[0] = 0.05;
			Q1D[1] = -sqrt(3.0/7.0);	W1D[1] = 49.0/180.0;
			Q1D[2] = 0.0;				W1D[2] = 16.0/45.0;
			Q1D[3] = sqrt(3.0/7.0);		W1D[3] = 49.0/180.0;
			Q1D[4] = 1.0;				W1D[4] = 0.05;
			break;
		default:
			fprintf(stderr,"Error: Quadrature rule of %s for %d nodes in the case of RECT not implemented.\n",toString().c_str(),N);
			break;
	}
	for (int k1 = 0; k1<N1d; k1++)
		for (int k2 = 0; k2<N1d; k2++) {
			int n = k1*N1d+k2;
			Q[n][0] = Q1D[k1];
			Q[n][1] = Q1D[k2];
			W[n] = W1D[k1]*W1D[k2];
		}
	if (Q1D!=NULL) delete[] Q1D;
	if (W1D!=NULL) delete[] W1D;
	return;
}

void QuadCoefGLob::compute_tetra(void) {
	switch (N) {
		case 4:
			Q[0][0] = 1.0; Q[0][1] = 0.0; Q[0][2] = 0.0; Q[0][3] = 0.0;
			W[0] = 1.0/4.0;
			Q[1][0] = 0.0; Q[1][1] = 1.0; Q[1][2] = 0.0; Q[1][3] = 0.0;
			W[1] = 1.0/4.0;
			Q[2][0] = 0.0; Q[2][1] = 0.0; Q[2][2] = 1.0; Q[2][3] = 0.0;
			W[2] = 1.0/4.0;
			Q[3][0] = 0.0; Q[3][1] = 0.0; Q[3][2] = 0.0; Q[3][3] = 1.0;
			W[3] = 1.0/4.0;
			break;
        case 5:
			Q[0][0] = 1.0; Q[0][1] = 0.0; Q[0][2] = 0.0; Q[0][3] = 0.0;
			W[0] = 1.0/20.0;
			Q[1][0] = 0.0; Q[1][1] = 1.0; Q[1][2] = 0.0; Q[1][3] = 0.0;
			W[1] = 1.0/20.0;
			Q[2][0] = 0.0; Q[2][1] = 0.0; Q[2][2] = 1.0; Q[2][3] = 0.0;
			W[2] = 1.0/20.0;
			Q[3][0] = 0.0; Q[3][1] = 0.0; Q[3][2] = 0.0; Q[3][3] = 1.0;
			W[3] = 1.0/20.0;
			Q[4][0] = 0.0; Q[4][1] = 0.0; Q[4][2] = 0.0; Q[4][3] = 0.0;
            W[4] = 4.0/5.0;
            break;
		default:
			fprintf(stderr,"Error: Quadrature rule of %s for %d nodes in the case of TETRA not implemented.\n",toString().c_str(),N);
			break;
	}
	return;
}

void QuadCoefGLob::compute_cuboid(void) {
	DType* Q1D = NULL;
	DType* W1D = NULL;
	int N1d = 0;
	switch (N) {
		case 8:
			N1d = 2;
			Q1D = new DType[N1d];		W1D = new DType[N1d];
			Q1D[0] = -1.0;				W1D[0] = 0.5;
			Q1D[1] = 1.0;				W1D[1] = 0.5;
			break;
		case 27:
			N1d = 3;
			Q1D = new DType[N1d];		W1D = new DType[N1d];
			Q1D[0] = -1.0;				W1D[0] = 1.0/6.0;
			Q1D[1] = 0.0;				W1D[1] = 2.0/3.0;
			Q1D[2] = 1.0;				W1D[2] = 1.0/6.0;
			break;
		case 64:
			N1d = 4;
			Q1D = new DType[N1d];		W1D = new DType[N1d];
			Q1D[0] = -1.0;				W1D[0] = 1.0/12.0;
			Q1D[1] = -sqrt(0.2);		W1D[1] = 5.0/12.0;
			Q1D[2] = sqrt(0.2);			W1D[2] = 5.0/12.0;
			Q1D[3] = 1.0;				W1D[3] = 1.0/12.0;
			break;
		case 125:
			N1d = 5;
			Q1D = new DType[N1d];		W1D = new DType[N1d];
			Q1D[0] = -1.0;				W1D[0] = 0.05;
			Q1D[1] = -sqrt(3.0/7.0);	W1D[1] = 49.0/180.0;
			Q1D[2] = 0.0;				W1D[2] = 16.0/45.0;
			Q1D[3] = sqrt(3.0/7.0);		W1D[3] = 49.0/180.0;
			Q1D[4] = 1.0;				W1D[4] = 0.05;
			break;
		default:
			fprintf(stderr,"Error: Quadrature rule of %s for %d nodes in the case of CUBOID not implemented.\n",toString().c_str(),N);
			break;
	}
	for (int k1 = 0; k1<N1d; k1++)
		for (int k2 = 0; k2<N1d; k2++) 
			for (int k3 = 0; k3<N1d; k3++) {
				int n = (k1*N1d+k2)*N1d+k3;
				Q[n][0] = Q1D[k1];
				Q[n][1] = Q1D[k2];
				Q[n][2] = Q1D[k3];
				W[n] = W1D[k1]*W1D[k2]*W1D[k3];
			}
	if (Q1D!=NULL) delete[] Q1D;
	if (W1D!=NULL) delete[] W1D;
	return;
}

bool QuadCoefGLob::checkType(void) {
	bool support;
	switch (type) {
		case EDGE:
		case RECT:
		case CUBOID:
		case TRIA:
		case TETRA:
			support = true;
			break;
		default:
			support = false;
			break;
	}
	return support;
}

/*
 *	Clenshaw-Curtis
 */

QuadCoefCC::~QuadCoefCC(void) {
	// Do nothing
}

string QuadCoefCC::toString(void) {
	char classNameString[40];
	strcpy(classNameString,"Clenshaw-Curtis");
	return classNameString;
}

void QuadCoefCC::setup_edge(void) {
	Np = -1;
	N = 2*order-1;
	return;
}

void QuadCoefCC::setup_tria(void) {
	return;
}

void QuadCoefCC::setup_rect(void) {
	Np = -1;
	int N1d = 2*order-1;
	N = N1d*N1d;
	return;
}

void QuadCoefCC::setup_tetra(void) {
	return;
}

void QuadCoefCC::setup_cuboid(void) {
	Np = -1;
	int N1d = 2*order-1;
	N = N1d*N1d*N1d;
	return;
}

void QuadCoefCC::compute_edge(void) {
	if (N==1) {
		Q[0][0] = 0.0;	W[0] = 1.0;
	}
	else {
		DType theta = PI/(1.0*N-1.0);
		for (int n=0; n<N; n++)	Q[n][0]	= cos((N-1-n)*theta);
		HVEC vec_f((N+1)/2,0.0);
		HMAT mat_d((N+1)/2,0.0);
		for (int n=0; n<(N+1)/2; n++) 
			vec_f[n] = 2.0/(1.0-(2.0*n)*(2.0*n));
		vec_f[0] 		*= 0.5;
		vec_f[(N+1)/2]	*= 0.5;
		for (int n=0; n<(N+1)/2; n++)
			for (int k=0; k<(N+1)/2; k++)
				mat_d[n][k] = (2.0/(1.0*N-1.0))*cos(2.0*n*k*theta);
		for (int k=0; k<(N+1)/2; k++) {
			mat_d[0][k]			*= 0.5;
			mat_d[(N-1)/2][k]	*= 0.5;
		}
		HVEC vec_w((N+1)/2,0.0);
		vec_w = mat_d*vec_f;
		for (int n=0; n<(N-1)/2; n++) {
			W[n]		= 0.5*vec_w[n];
			W[N-1-n]	= 0.5*vec_w[n];
		}
		W[(N-1)/2] = vec_w[(N-1)/2];
	}
	return;
}

void QuadCoefCC::compute_tria(void) {
	return;
}

void QuadCoefCC::compute_rect(void) {
	int N1d = sqrt(N);
	assert(N==(N1d*N1d));
	DType* Q1D = new DType[N1d];
	DType* W1D = new DType[N1d];
	if (N1d==1) {
		Q1D[0] = 0.0;	W1D[0] = 1.0;
	}
	else {
		DType theta = PI/(1.0*N1d-1.0);
		for (int n=0; n<N1d; n++)	Q1D[n]	= cos((N1d-1-n)*theta);
		HVEC vec_f((N1d+1)/2,0.0);
		HMAT mat_d((N1d+1)/2,0.0);
		for (int n=0; n<(N1d+1)/2; n++) 
			vec_f[n] = 2.0/(1.0-(2.0*n)*(2.0*n));
		vec_f[0] 		*= 0.5;
		vec_f[(N1d+1)/2]	*= 0.5;
		for (int n=0; n<(N1d+1)/2; n++)
			for (int k=0; k<(N1d+1)/2; k++)
				mat_d[n][k] = (2.0/(1.0*N1d-1.0))*cos(2.0*n*k*theta);
		for (int k=0; k<(N1d+1)/2; k++) {
			mat_d[0][k]			*= 0.5;
			mat_d[(N1d-1)/2][k]	*= 0.5;
		}
		HVEC vec_w((N1d+1)/2,0.0);
		vec_w = mat_d*vec_f;
		for (int n=0; n<(N1d-1)/2; n++) {
			W1D[n]		= 0.5*vec_w[n];
			W1D[N1d-1-n]	= 0.5*vec_w[n];
		}
		W1D[(N1d-1)/2] = vec_w[(N1d-1)/2];
	}
	for (int k1 = 0; k1<N1d; k1++)
		for (int k2 = 0; k2<N1d; k2++) {
			int n = k1*N1d+k2;
			Q[n][0] = Q1D[k1];
			Q[n][1] = Q1D[k2];
			W[n] = W1D[k1]*W1D[k2];
		}
	if (Q1D!=NULL) delete[] Q1D;
	if (W1D!=NULL) delete[] W1D;
	return;
}

void QuadCoefCC::compute_tetra(void) {
	return;
}

void QuadCoefCC::compute_cuboid(void) {
	int N1d = sqrt(N);
	assert(N==(N1d*N1d*N1d));
	DType* Q1D = new DType[N1d];
	DType* W1D = new DType[N1d];
	if (N1d==1) {
		Q1D[0] = 0.0;	W1D[0] = 1.0;
	}
	else {
		DType theta = PI/(1.0*N1d-1.0);
		for (int n=0; n<N1d; n++)	Q1D[n]	= cos((N1d-1-n)*theta);
		HVEC vec_f((N1d+1)/2,0.0);
		HMAT mat_d((N1d+1)/2,0.0);
		for (int n=0; n<(N1d+1)/2; n++) 
			vec_f[n] = 2.0/(1.0-(2.0*n)*(2.0*n));
		vec_f[0] 		*= 0.5;
		vec_f[(N1d+1)/2]	*= 0.5;
		for (int n=0; n<(N1d+1)/2; n++)
			for (int k=0; k<(N1d+1)/2; k++)
				mat_d[n][k] = (2.0/(1.0*N1d-1.0))*cos(2.0*n*k*theta);
		for (int k=0; k<(N1d+1)/2; k++) {
			mat_d[0][k]			*= 0.5;
			mat_d[(N1d-1)/2][k]	*= 0.5;
		}
		HVEC vec_w((N1d+1)/2,0.0);
		vec_w = mat_d*vec_f;
		for (int n=0; n<(N1d-1)/2; n++) {
			W1D[n]		= 0.5*vec_w[n];
			W1D[N1d-1-n]	= 0.5*vec_w[n];
		}
		W1D[(N1d-1)/2] = vec_w[(N1d-1)/2];
	}
	for (int k1 = 0; k1<N1d; k1++)
		for (int k2 = 0; k2<N1d; k2++)
			for (int k3 = 0; k3<N1d; k3++) {
				int n = (k1*N1d+k2)*N1d+k3;
				Q[n][0] = Q1D[k1];
				Q[n][1] = Q1D[k2];
				Q[n][2] = Q1D[k3];
				W[n] = W1D[k1]*W1D[k2]*W1D[k3];
			}
	if (Q1D!=NULL) delete[] Q1D;
	if (W1D!=NULL) delete[] W1D;
	return;
}

bool QuadCoefCC::checkType(void) {
	bool support;
	switch (type) {
		case EDGE:
		case RECT:
		case CUBOID:
			support = true;
			break;
		case TRIA:
		case TETRA:
			support = false;
			break;
		default:
			support = false;
			break;
	}
	return support;
}

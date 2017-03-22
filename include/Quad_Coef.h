#ifndef QUAD_COEF_H
#define QUAD_COEF_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "Data_Type.h"
#include "Comput_Choices.h"
#include "LowDim_Alg.h"
using namespace std;

class QuadCoef;
QuadCoef* newQuadCoef(ChoiceQuadRule,ChoiceQuadType,int);

/*
 *	General remark
 *	- Though the name Quad(rature)Coef(ficients)
 *	  The following objects are also used to construct cubatures in high-
 *	  dimensions.
 *	- Strictly speaking, there is no general rule for either 
 *	  Gauss-Legendre, or Gauss-Lobatto, or Clenshaw-Curtis for simplex
 *	  in multiple dimensions, thus they are only included in the 
 *	  corresponding class conceptually. For example, Gauss-Lobatto for
 *	  triangles have points on corners as well as edges; whereas Gauss-
 *	  Legendre not
 *	- Due to the reason above, each implementation for TRIA or TETRA MUST
 *	  be supplemented by documentary, which is kepted in the file
 *	  docs/QuadCoef.pdf
 */

/*
 *	Base class for quadrature coefficients
 */
class QuadCoef {
	protected:
		enum Type {EDGE=0, TRIA=1, RECT=2, TETRA=3, CUBOID=4} type;
		int	dim;	// Number of coordinates:
					// - Relative to the center of the object
					// - Not necessary independent
					// 1 for EDGE
					// 2 for RECT
					// 3 for TRIA and CUBOID
					// 4 for TETRA
		int N;		// Number of quadrature points
		int Np;		// Allowed exact evaluation of polynomials of degree Np
		int order;	// Targeted order
		DType* 	W;	// A single array of weights of size N
		DType** Q;	// Quadrature coordinate: N X dim

		void setup(void);					// Setup N, Np given order
		virtual string toString(void) = 0;	// Object name to string
		virtual void setup_edge(void) = 0;
		virtual void setup_tria(void) = 0;
		virtual void setup_rect(void) = 0;
		virtual void setup_tetra(void) = 0;
		virtual void setup_cuboid(void) = 0;
		virtual void compute_edge(void) = 0;
		virtual void compute_tria(void) = 0;
		virtual void compute_rect(void) = 0;
		virtual void compute_tetra(void) = 0;
		virtual void compute_cuboid(void) = 0;
	public:
		QuadCoef(ChoiceQuadType);
		virtual ~QuadCoef(void);

		int getNquad(void) const {return N;}
		int getNpoly(void) const {return Np;}
		int getOrder(void) const {return order;}
		DType getWeight(int n) const {return W[n];}
		DType getQuad(int n, int d=0) const {return Q[n][d];}
		
		virtual bool checkType(void) = 0;
};

/*
 *	Gaussian-Legendre quadrature: 
 *	Minimum # of points for given desired order
 *	Accurate for polynomials with order 2n-1 given n points in 1D
 *	- EDGE -
 *	  region for Q : [-1, 1]
 *	  region for W : [-0.5, 0.5]
 *	- TRIA -
 *	- RECT -
 *	  region for Q : [-1, 1]^2
 *	  region for W : [-0.5, 0.5]^2
 *	- TETRA - 
 *	- CUBOID -
 *	  region for Q : [-1, 1]^3
 *	  region for W : [-0.5, 0.5]^3
 */
class QuadCoefGLeg : public QuadCoef {
	protected:
		string toString(void);
		void setup_edge(void);
		void setup_tria(void);
		void setup_rect(void);
		void setup_tetra(void);
		void setup_cuboid(void);
		void compute_edge(void);
		void compute_tria(void);
		void compute_rect(void);
		void compute_tetra(void);
		void compute_cuboid(void);
	public:
		QuadCoefGLeg(ChoiceQuadType quadType, int order_in) :
			QuadCoef(quadType) {
			order = order_in;
			setup();
		}
		~QuadCoefGLeg(void);

		bool checkType(void);
};

/*
 *	Gaussian-Lobatto
 *	Inclusion of end points of interval in 1D
 *	Accurate for polynomials with order 2n-3 given n points in 1D
 *	- EDGE -
 *	  region for Q : [-1, 1]
 *	  region for W : [-0.5, 0.5]
 *	- TRIA -
 *	- RECT -
 *	  region for Q : [-1, 1]^2
 *	  region for W : [-0.5, 0.5]^2
 *	- TETRA - 
 *	- CUBOID -
 *	  region for Q : [-1, 1]^3
 *	  region for W : [-0.5, 0.5]^3
 */
class QuadCoefGLob : public QuadCoef {
	protected:
		string toString(void);
		void setup_edge(void);
		void setup_tria(void);
		void setup_rect(void);
		void setup_tetra(void);
		void setup_cuboid(void);
		void compute_edge(void);
		void compute_tria(void);
		void compute_rect(void);
		void compute_tetra(void);
		void compute_cuboid(void);
	public:
		QuadCoefGLob(ChoiceQuadType quadType, int order_in) :
			QuadCoef(quadType) {
			order = order_in;
			setup();
		}
		~QuadCoefGLob(void);

		bool checkType(void);
};

/*
 *	Clenshaw-Curtis quadrature
 *	Use Chebyshev polynomial's roots:
 *	Allows extension to higher-order by expanding sample points: reuse existing nodes
 *	The error is not polynomial order; hence the input order represents:
 *	2*order-1 is the number of sample points in [-1,1]
 *	This method is available for all input order
 *	- EDGE -
 *	  region for Q : [-1, 1]
 *	  region for W : [-0.5, 0.5]
 *	- TRIA -
 *	- RECT -
 *	  region for Q : [-1, 1]^2
 *	  region for W : [-0.5, 0.5]^2
 *	- TETRA - 
 *	- CUBOID -
 *	  region for Q : [-1, 1]^3
 *	  region for W : [-0.5, 0.5]^3
 */
class QuadCoefCC : public QuadCoef {
	protected:
		string toString(void);
		void setup_edge(void);
		void setup_tria(void);
		void setup_rect(void);
		void setup_tetra(void);
		void setup_cuboid(void);
		void compute_edge(void);
		void compute_tria(void);
		void compute_rect(void);
		void compute_tetra(void);
		void compute_cuboid(void);
	public:
		QuadCoefCC(ChoiceQuadType quadType, int order_in) :
			QuadCoef(quadType) {
			order = order_in;
			fprintf(stdout,"Warning: Clenshaw-Curtis rule is employed -- the input order corresponds to the case of 2*order-1 sample points; and Np is always -1.\n");
			setup();
		}
		~QuadCoefCC(void);

		bool checkType(void);
};

#endif

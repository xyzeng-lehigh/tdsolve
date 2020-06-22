#include "Alg_Limiters.h"

// No limiter: linear reconstruction with no limiter
DType Alg_LimSlope_Linear(DType U_c, DType U_l, DType U_r, DType dx_c, DType dx_l, DType dx_r) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	return ALG_LIMSLOPE_LINEAR_L*2.0*(U_c-U_l)/(dx_c+dx_l)+
		   ALG_LIMSLOPE_LINEAR_R*2.0*(U_r-U_c)/(dx_r+dx_c)+
		   ALG_LIMSLOPE_LINEAR_C*2.0*(U_r-U_l)/(dx_r+2.0*dx_c+dx_l);
}

VEC3 Alg_LimSlope_Linear(VEC3 U_c, VEC3 U_l, VEC3 U_r, DType dx_c, DType dx_l, DType dx_r) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	int i;
	VEC3 slope;
	for (i=0; i<3; i++) {
		slope[i] = ALG_LIMSLOPE_LINEAR_L*2.0*(U_c[i]-U_l[i])/(dx_c+dx_l)+
				   ALG_LIMSLOPE_LINEAR_R*2.0*(U_r[i]-U_c[i])/(dx_r+dx_c)+
				   ALG_LIMSLOPE_LINEAR_C*2.0*(U_r[i]-U_l[i])/(dx_r+2.0*dx_c+dx_l);
	}
	return slope;
}

VEC4 Alg_LimSlope_Linear(VEC4 U_c, VEC4 U_l, VEC4 U_r, DType dx_c, DType dx_l, DType dx_r) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	int i;
	VEC4 slope;
	for (i=0; i<4; i++) {
		slope[i] = ALG_LIMSLOPE_LINEAR_L*2.0*(U_c[i]-U_l[i])/(dx_c+dx_l)+
				   ALG_LIMSLOPE_LINEAR_R*2.0*(U_r[i]-U_c[i])/(dx_r+dx_c)+
				   ALG_LIMSLOPE_LINEAR_C*2.0*(U_r[i]-U_l[i])/(dx_r+2.0*dx_c+dx_l);
	}
	return slope;
}


// Minmod limiter: non-uniform mesh implemented
DType Alg_LimSlope_Minmod(DType U_c, DType U_l, DType U_r, DType dx_c, DType dx_l, DType dx_r, DType coef_a, DType coef_b) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	// phi(a)=b must be satisfied for second order accuracy
//	DType coef_a, coef_b;
	if (coef_a==-1)
		coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	if (coef_b==-1)
		coef_b = 2.0*dx_c/(dx_c+dx_r);

	if (U_r>U_c)
		return coef_b*max(0.0,min( (U_r-U_c) , (U_c-U_l)/coef_a )/dx_c);
	else if (U_r<U_c)
		return coef_b*min(0.0,max( (U_r-U_c) , (U_c-U_l)/coef_a )/dx_c);
	else
		return 0.0;
}

VEC3 Alg_LimSlope_Minmod(VEC3 U_c, VEC3 U_l, VEC3 U_r, DType dx_c, DType dx_l, DType dx_r, DType coef_a, DType coef_b) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	// phi(a)=b must be satisfied for second order accuracy
//	DType coef_a, coef_b;
	if (coef_a==-1)
		coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	if (coef_b==-1)
		coef_b = 2.0*dx_c/(dx_c+dx_r);

	int i;
	VEC3 slope;
	for (i=0; i<3; i++) {
		if (U_r[i]>U_c[i])
			slope[i] = coef_b*max(0.0,min( (U_r[i]-U_c[i]) , (U_c[i]-U_l[i])/coef_a )/dx_c);
		else if (U_r[i]<U_c[i])
			slope[i] = coef_b*min(0.0,max( (U_r[i]-U_c[i]) , (U_c[i]-U_l[i])/coef_a )/dx_c);
		else
			slope[i] = 0.0;
	}
	return slope;
}

VEC4 Alg_LimSlope_Minmod(VEC4 U_c, VEC4 U_l, VEC4 U_r, DType dx_c, DType dx_l, DType dx_r, DType coef_a, DType coef_b) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	// phi(a)=b must be satisfied for second order accuracy
//	DType coef_a, coef_b;
	if (coef_a==-1)
		coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	if (coef_b==-1)
		coef_b = 2.0*dx_c/(dx_c+dx_r);

	int i;
	VEC4 slope;
	for (i=0; i<4; i++) {
		if (U_r[i]>U_c[i])
			slope[i] = coef_b*max(0.0,min( (U_r[i]-U_c[i]) , (U_c[i]-U_l[i])/coef_a )/dx_c);
		else if (U_r[i]<U_c[i])
			slope[i] = coef_b*min(0.0,max( (U_r[i]-U_c[i]) , (U_c[i]-U_l[i])/coef_a )/dx_c);
		else
			slope[i] = 0.0;
	}
	return slope;
}



// Superbee limiter: non-uniform mesh implemented
DType Alg_LimSlope_Superbee(DType U_c, DType U_l, DType U_r, DType dx_c, DType dx_l, DType dx_r, DType coef_a, DType coef_b) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	// phi(a)=b must be satisfied for second order accuracy
//	DType coef_a, coef_b;
	if (coef_a==-1)
		coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	if (coef_b==-1)
		coef_b = 2.0*dx_c/(dx_c+dx_r);

	if (U_r>U_c)
		return max(0.0,max(min(2.0*(U_c-U_l),coef_b*(U_r-U_c)),
						   min((coef_b/coef_a)*(U_c-U_l),2.0*(U_r-U_c))))/dx_c;
	else if (U_r<U_c)
		return min(0.0,min(max(2.0*(U_c-U_l),coef_b*(U_r-U_c)),
						   max((coef_b/coef_a)*(U_c-U_l),2.0*(U_r-U_c))))/dx_c;
	else
		return 0.0;
}

VEC3 Alg_LimSlope_Superbee(VEC3 U_c, VEC3 U_l, VEC3 U_r, DType dx_c, DType dx_l, DType dx_r, DType coef_a, DType coef_b) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	// phi(a)=b must be satisfied for second order accuracy
//	DType coef_a, coef_b;
	if (coef_a==-1)
		coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	if (coef_b==-1)
		coef_b = 2.0*dx_c/(dx_c+dx_r);

	int i;
	VEC3 slope;
	for (i=0; i<3; i++) {
		if (U_r[i]>U_c[i])
			slope[i] = max(0.0,max(min(2.0*(U_c[i]-U_l[i]),coef_b*(U_r[i]-U_c[i])),
							   min((coef_b/coef_a)*(U_c[i]-U_l[i]),2.0*(U_r[i]-U_c[i]))))/dx_c;
		else if (U_r[i]<U_c[i])
			slope[i] = min(0.0,min(max(2.0*(U_c[i]-U_l[i]),coef_b*(U_r[i]-U_c[i])),
							   max((coef_b/coef_a)*(U_c[i]-U_l[i]),2.0*(U_r[i]-U_c[i]))))/dx_c;
		else
			slope[i] = 0.0;
	}
	return slope;
}

VEC4 Alg_LimSlope_Superbee(VEC4 U_c, VEC4 U_l, VEC4 U_r, DType dx_c, DType dx_l, DType dx_r, DType coef_a, DType coef_b) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	// phi(a)=b must be satisfied for second order accuracy
//	DType coef_a, coef_b;
	if (coef_a==-1)
		coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	if (coef_b==-1)
		coef_b = 2.0*dx_c/(dx_c+dx_r);

	int i;
	VEC4 slope;
	for (i=0; i<4; i++) {
		if (U_r[i]>U_c[i])
			slope[i] = max(0.0,max(min(2.0*(U_c[i]-U_l[i]),coef_b*(U_r[i]-U_c[i])),
							   min((coef_b/coef_a)*(U_c[i]-U_l[i]),2.0*(U_r[i]-U_c[i]))))/dx_c;
		else if (U_r[i]<U_c[i])
			slope[i] = min(0.0,min(max(2.0*(U_c[i]-U_l[i]),coef_b*(U_r[i]-U_c[i])),
							   max((coef_b/coef_a)*(U_c[i]-U_l[i]),2.0*(U_r[i]-U_c[i]))))/dx_c;
		else
			slope[i] = 0.0;
	}
	return slope;
}


// MC limiter: non-uniform mesh implemented
DType Alg_LimSlope_MC(DType U_c, DType U_l, DType U_r, DType dx_c, DType dx_l, DType dx_r, DType coef_a, DType coef_b) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	// phi(a)=b must be satisfied for second order accuracy
//	DType coef_a, coef_b;
	if (coef_a==-1)
		coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	if (coef_b==-1)
		coef_b = 2.0*dx_c/(dx_c+dx_r);

	if (U_r>U_c)
		return max(0.0,min(2.0*(U_r-U_c),
						   min(2.0*(U_c-U_l),(coef_b/(coef_a+1.0))*(U_r-U_l))))/dx_c;
	else if (U_r<U_c)
		return min(0.0,max(2.0*(U_r-U_c),
						   max(2.0*(U_c-U_l),(coef_b/(coef_a+1.0))*(U_r-U_l))))/dx_c;
	else
		return 0.0;
}

VEC3 Alg_LimSlope_MC(VEC3 U_c, VEC3 U_l, VEC3 U_r, DType dx_c, DType dx_l, DType dx_r, DType coef_a, DType coef_b) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	// phi(a)=b must be satisfied for second order accuracy
//	DType coef_a, coef_b;
	if (coef_a==-1)
		coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	if (coef_b==-1)
		coef_b = 2.0*dx_c/(dx_c+dx_r);

	int i;
	VEC3 slope;
	for (i=0; i<3; i++) {
		if (U_r[i]>U_c[i])
			slope[i] = max(0.0,min(2.0*(U_r[i]-U_c[i]),
							   min(2.0*(U_c[i]-U_l[i]),(coef_b/(coef_a+1.0))*(U_r[i]-U_l[i]))))/dx_c;
		else if (U_r[i]<U_c[i])
			slope[i] = min(0.0,max(2.0*(U_r[i]-U_c[i]),
							   max(2.0*(U_c[i]-U_l[i]),(coef_b/(coef_a+1.0))*(U_r[i]-U_l[i]))))/dx_c;
		else
			slope[i] = 0.0;
	}
	return slope;
}

VEC4 Alg_LimSlope_MC(VEC4 U_c, VEC4 U_l, VEC4 U_r, DType dx_c, DType dx_l, DType dx_r, DType coef_a, DType coef_b) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	// phi(a)=b must be satisfied for second order accuracy
//	DType coef_a, coef_b;
	if (coef_a==-1)
		coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	if (coef_b==-1)
		coef_b = 2.0*dx_c/(dx_c+dx_r);

	int i;
	VEC4 slope;
	for (i=0; i<4; i++) {
		if (U_r[i]>U_c[i])
			slope[i] = max(0.0,min(2.0*(U_r[i]-U_c[i]),
							   min(2.0*(U_c[i]-U_l[i]),(coef_b/(coef_a+1.0))*(U_r[i]-U_l[i]))))/dx_c;
		else if (U_r[i]<U_c[i])
			slope[i] = min(0.0,max(2.0*(U_r[i]-U_c[i]),
							   max(2.0*(U_c[i]-U_l[i]),(coef_b/(coef_a+1.0))*(U_r[i]-U_l[i]))))/dx_c;
		else
			slope[i] = 0.0;
	}
	return slope;
}

// van Leer limiter
//
// calculate the parameter k_vl in generalized van Leer limiter
int Alg_LimSlope_vanLeer_calc_k(DType a_in, DType b_in, DType& denom_ot, DType& nomtr_ot) {
	int k_ot = 1;
	denom_ot = a_in+1.0;
	nomtr_ot = a_in;
	while (k_ot<K_VANLEER_MAX) {
		if (b_in<=(2.0*nomtr_ot/denom_ot))
			return k_ot;
		denom_ot = a_in*denom_ot+1.0;
		nomtr_ot = denom_ot-1.0;
		k_ot = k_ot+1;
	}
	return k_ot;
}

DType Alg_LimSlope_vanLeer(DType U_c, DType U_l, DType U_r, DType dx_c, DType dx_l, DType dx_r, int k_vl, DType denom, DType nomtr, DType coef_a, DType coef_b) {
	if ((dx_l==-1)&&(dx_r==-1))
		k_vl = 0;			// Indicate uniform meshes
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;
	if (k_vl<-1) {
		printf("Error in Alg_LimSlope_vanLeer, k_vl should be >= -1\n");
		exit(8);
	}
	else if ((k_vl>0)&&((denom<0)||(nomtr<0))) {
		printf("Error in Alg_LimSlope_vanLeer, for manual k_vl>0, denom and nomtr must also be provided\n");
		exit(8);
	}

	// phi(a)=b must be satisfied for second order accuracy
//	DType coef_a, coef_b;
	if (coef_a==-1)
		coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	if (coef_b==-1)
		coef_b = 2.0*dx_c/(dx_c+dx_r);

	if (k_vl==0) {			// Uniform meshes
		k_vl = 1;
		denom = 2.0;
		nomtr = 1.0;
	}
	else if (k_vl==-1)		// Automatical determination
		k_vl = Alg_LimSlope_vanLeer_calc_k(coef_a,coef_b,denom,nomtr);
							// Otherwise, use inputed value of k_vl

	DType par_theta;		// This can be either the theta described in header
							// or the inverse of the theta described there.
	int k;
	DType denom_theta, nomtr_theta;

	if (((U_c-U_l)*(U_r-U_c))<0.0)
		return 0.0;

	if (abs(U_c-U_l)<abs(U_r-U_c)) {
		par_theta = (U_c-U_l)/(U_r-U_c);
		if (k_vl==1) {
			if (U_r>U_c)
				return (0.5*coef_b*denom/nomtr)*((U_c-U_l)+abs(U_c-U_l))/((1.0+abs(par_theta))*dx_c);
			else if (U_r<U_c)
				return (0.5*coef_b*denom/nomtr)*((U_c-U_l)-abs(U_c-U_l))/((1.0+abs(par_theta))*dx_c);
			else
				return 0.0;
		}
		else {
			nomtr_theta = par_theta;
			denom_theta = par_theta+1.0;
			for (k=1; k<k_vl; k++) {
				denom_theta = par_theta*denom_theta+1.0;
				nomtr_theta = denom_theta-1.0;
			}
			if (par_theta>0.0)
				return ((U_r-U_c)/dx_c)*coef_b*(nomtr_theta/denom_theta)*(denom/nomtr);
			else
				return 0.0;
		}
	}
	else if (abs(U_c-U_l)>abs(U_r-U_c)) {
		par_theta = (U_r-U_c)/(U_c-U_l);
		if (k_vl==1) {
			if (U_c>U_l)
				return (0.5*coef_b*denom/nomtr)*((U_r-U_c)+abs(U_r-U_c))/((1.0+abs(par_theta))*dx_c);
			else if (U_c<U_l)
				return (0.5*coef_b*denom/nomtr)*((U_r-U_c)-abs(U_r-U_c))/((1.0+abs(par_theta))*dx_c);
			else
				return 0.0;
		}
		else {
			nomtr_theta = par_theta;
			denom_theta = par_theta+1.0;
			for (k=1; k<k_vl; k++) {
				denom_theta = par_theta*denom_theta+1.0;
				nomtr_theta = denom_theta-1.0;
			}
			if (par_theta>0.0)
				return ((U_c-U_l)/dx_c)*coef_b*(nomtr_theta/denom_theta)*(denom/nomtr);
			else
				return 0.0;
		}
	}
	else {
		if (U_r==U_l)
			return 0.0;
		else
			return (coef_b*(k_vl/(k_vl+1.0))*denom/nomtr)*(U_r-U_c)/dx_c;
	}
}

VEC3 Alg_LimSlope_vanLeer(VEC3 U_c, VEC3 U_l, VEC3 U_r, DType dx_c, DType dx_l, DType dx_r, int k_vl, DType denom, DType nomtr, DType coef_a, DType coef_b) {
	if ((dx_l==-1)&&(dx_r==-1))
		k_vl = 0;			// Indicate uniform meshes
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;
	if (k_vl<-1) {
		printf("Error in Alg_LimSlope_vanLeer, k_vl should be >= -1\n");
		exit(8);
	}
	else if ((k_vl>0)&&((denom<0)||(nomtr<0))) {
		printf("Error in Alg_LimSlope_vanLeer, for manual k_vl>0, denom and nomtr must also be provided\n");
		exit(8);
	}

	// phi(a)=b must be satisfied for second order accuracy
//	DType coef_a, coef_b;
	if (coef_a==-1)
		coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	if (coef_b==-1)
		coef_b = 2.0*dx_c/(dx_c+dx_r);

	if (k_vl==0) {			// Uniform meshes
		k_vl = 1;
		denom = 2.0;
		nomtr = 1.0;
	}
	else if (k_vl==-1)		// Automatical determination
		k_vl = Alg_LimSlope_vanLeer_calc_k(coef_a,coef_b,denom,nomtr);
							// Otherwise, use inputed value of k_vl

	DType par_theta;		// This can be either the theta described in header
							// or the inverse of the theta described there
	int i, k;
	DType denom_theta, nomtr_theta;
	VEC3 slope;
	for (i=0; i<3; i++) {
		if (((U_c[i]-U_l[i])*(U_r[i]-U_c[i]))<0.0)
			slope[i] = 0.0;

		if (abs(U_c[i]-U_l[i])<abs(U_r[i]-U_c[i])) {
			par_theta = (U_c[i]-U_l[i])/(U_r[i]-U_c[i]);
			if (k_vl==1) {
				if (U_r[i]>U_c[i])
					slope[i] = (0.5*coef_b*denom/nomtr)*((U_c[i]-U_l[i])+abs(U_c[i]-U_l[i]))/((1.0+abs(par_theta))*dx_c);
				else if (U_r[i]<U_c[i])
					slope[i] = (0.5*coef_b*denom/nomtr)*((U_c[i]-U_l[i])-abs(U_c[i]-U_l[i]))/((1.0+abs(par_theta))*dx_c);
				else
					slope[i] = 0.0;
			}
			else {
				nomtr_theta = par_theta;
				denom_theta = par_theta+1.0;
				for (k=1; k<k_vl; k++) {
					denom_theta = par_theta*denom_theta+1.0;
					nomtr_theta = denom_theta-1.0;
				}
				if (par_theta>0.0)
					slope[i] = ((U_r[i]-U_c[i])/dx_c)*coef_b*(nomtr_theta/denom_theta)*(denom/nomtr);
				else
					slope[i] = 0.0;
			}
		}
		else if (abs(U_c[i]-U_l[i])>abs(U_r[i]-U_c[i])) {
			par_theta = (U_r[i]-U_c[i])/(U_c[i]-U_l[i]);
			if (k_vl==1) {
				if (U_c[i]>U_l[i])
					slope[i] = (0.5*coef_b*denom/nomtr)*((U_r[i]-U_c[i])+abs(U_r[i]-U_c[i]))/((1.0+abs(par_theta))*dx_c);
				else if (U_c[i]<U_l[i])
					slope[i] = (0.5*coef_b*denom/nomtr)*((U_r[i]-U_c[i])-abs(U_r[i]-U_c[i]))/((1.0+abs(par_theta))*dx_c);
				else
					slope[i] = 0.0;
			}
			else {
				nomtr_theta = par_theta;
				denom_theta = par_theta+1.0;
				for (k=1; k<k_vl; k++) {
					denom_theta = par_theta*denom_theta+1.0;
					nomtr_theta = denom_theta-1.0;
				}
				if (par_theta>0.0)
					slope[i] = ((U_c[i]-U_l[i])/dx_c)*coef_b*(nomtr_theta/denom_theta)*(denom/nomtr);
				else
					slope[i] = 0.0;
			}
		}
		else {
			if (U_r[i]==U_l[i])
				slope[i] = 0.0;
			else
				slope[i] = (coef_b*(k_vl/(k_vl+1.0))*denom/nomtr)*(U_r[i]-U_c[i])/dx_c;
		}
	}
	return slope;
}

VEC4 Alg_LimSlope_vanLeer(VEC4 U_c, VEC4 U_l, VEC4 U_r, DType dx_c, DType dx_l, DType dx_r, int k_vl, DType denom, DType nomtr, DType coef_a, DType coef_b) {
	if ((dx_l==-1)&&(dx_r==-1))
		k_vl = 0;			// Indicate uniform meshes
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;
	if (k_vl<-1) {
		printf("Error in Alg_LimSlope_vanLeer, k_vl should be >= -1\n");
		exit(8);
	}
	else if ((k_vl>0)&&((denom<0)||(nomtr<0))) {
		printf("Error in Alg_LimSlope_vanLeer, for manual k_vl>0, denom and nomtr must also be provided\n");
		exit(8);
	}

	// phi(a)=b must be satisfied for second order accuracy
//	DType coef_a, coef_b;
	if (coef_a==-1)
		coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	if (coef_b==-1)
		coef_b = 2.0*dx_c/(dx_c+dx_r);

	if (k_vl==0) {			// Uniform meshes
		k_vl = 1;
		denom = 2.0;
		nomtr = 1.0;
	}
	else if (k_vl==-1)		// Automatical determination
		k_vl = Alg_LimSlope_vanLeer_calc_k(coef_a,coef_b,denom,nomtr);
							// Otherwise, use inputed value of k_vl

	DType par_theta;		// This can be either the theta described in header
							// or the inverse of the theta described there
	int i, k;
	DType denom_theta, nomtr_theta;
	VEC4 slope;
	for (i=0; i<4; i++) {
		if (((U_c[i]-U_l[i])*(U_r[i]-U_c[i]))<0.0)
			slope[i] = 0.0;

		if (abs(U_c[i]-U_l[i])<abs(U_r[i]-U_c[i])) {
			par_theta = (U_c[i]-U_l[i])/(U_r[i]-U_c[i]);
			if (k_vl==1) {
				if (U_r[i]>U_c[i])
					slope[i] = (0.5*coef_b*denom/nomtr)*((U_c[i]-U_l[i])+abs(U_c[i]-U_l[i]))/((1.0+abs(par_theta))*dx_c);
				else if (U_r[i]<U_c[i])
					slope[i] = (0.5*coef_b*denom/nomtr)*((U_c[i]-U_l[i])-abs(U_c[i]-U_l[i]))/((1.0+abs(par_theta))*dx_c);
				else
					slope[i] = 0.0;
			}
			else {
				nomtr_theta = par_theta;
				denom_theta = par_theta+1.0;
				for (k=1; k<k_vl; k++) {
					denom_theta = par_theta*denom_theta+1.0;
					nomtr_theta = denom_theta-1.0;
				}
				if (par_theta>0.0)
					slope[i] = ((U_r[i]-U_c[i])/dx_c)*coef_b*(nomtr_theta/denom_theta)*(denom/nomtr);
				else
					slope[i] = 0.0;
			}
		}
		else if (abs(U_c[i]-U_l[i])>abs(U_r[i]-U_c[i])) {
			par_theta = (U_r[i]-U_c[i])/(U_c[i]-U_l[i]);
			if (k_vl==1) {
				if (U_c[i]>U_l[i])
					slope[i] = (0.5*coef_b*denom/nomtr)*((U_r[i]-U_c[i])+abs(U_r[i]-U_c[i]))/((1.0+abs(par_theta))*dx_c);
				else if (U_c[i]<U_l[i])
					slope[i] = (0.5*coef_b*denom/nomtr)*((U_r[i]-U_c[i])-abs(U_r[i]-U_c[i]))/((1.0+abs(par_theta))*dx_c);
				else
					slope[i] = 0.0;
			}
			else {
				nomtr_theta = par_theta;
				denom_theta = par_theta+1.0;
				for (k=1; k<k_vl; k++) {
					denom_theta = par_theta*denom_theta+1.0;
					nomtr_theta = denom_theta-1.0;
				}
				if (par_theta>0.0)
					slope[i] = ((U_c[i]-U_l[i])/dx_c)*coef_b*(nomtr_theta/denom_theta)*(denom/nomtr);
				else
					slope[i] = 0.0;
			}
		}
		else {
			if (U_r[i]==U_l[i])
				slope[i] = 0.0;
			else
				slope[i] = (coef_b*(k_vl/(k_vl+1.0))*denom/nomtr)*(U_r[i]-U_c[i])/dx_c;
		}
	}
	return slope;
}


// van Albada limiter
// calculate the parameter k_va in generalized van Albada limiter
int Alg_LimSlope_vanAlbada_calc_k(DType a_in, DType b_in) {
	int k_ot = 2;
	DType var_dum;
	while (k_ot<K_VANALBADA_MAX) {
		var_dum = 1.0+pow((k_ot-1.0)/k_ot,k_ot-1)/k_ot;
		if ((b_in<=2.0*a_in/var_dum)&&(b_in<=2.0/var_dum))
			return k_ot;
		k_ot = k_ot+1;
	}
	return k_ot;
}

DType Alg_LimSlope_vanAlbada(DType U_c, DType U_l, DType U_r, DType dx_c, DType dx_l, DType dx_r, int k_va, DType coef_a, DType coef_b) {
	if ((dx_l==-1)&&(dx_r==-1))
		k_va = 0;			// Indicate uniform meshes
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;
	if ((k_va!=0)&&(k_va!=-1)&&(k_va<2)) {
		printf("Error in Alg_LimSlope_vanAlbada, k_va should be >2 or =-1\n");
		exit(8);
	}

	// phi(a)=b must be satisfied for second order accuracy
//	DType coef_a, coef_b;
	if (coef_a==-1)
		coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	if (coef_b==-1)
		coef_b = 2.0*dx_c/(dx_c+dx_r);

	if (k_va==0)
		k_va = 2;
	else if (k_va==-1)
		k_va = Alg_LimSlope_vanAlbada_calc_k(coef_a,coef_b);

	DType par_theta;		// This can be either the theta described in header
							// or the inverse of the theta described there
	if (((U_c-U_l)*(U_r-U_c))<0.0)
		return 0.0;

	if (abs(U_c-U_l)<abs(U_r-U_c)) {
		par_theta = (U_c-U_l)/(U_r-U_c);
		return ((coef_b*(1.0+pow(par_theta,k_va-1)))
				/(coef_a+pow(par_theta,k_va)))*(U_c-U_l)/dx_c;
	}
	else if (abs(U_c-U_l)>abs(U_r-U_c)) {
		par_theta = (U_r-U_c)/(U_c-U_l);
		return ((coef_b*(1.0+pow(par_theta,k_va-1)))
				/(1.0+coef_a*pow(par_theta,k_va)))*(U_r-U_c)/dx_c;
	}
	else {
		if (U_r==U_l)
			return 0.0;
		else
			return (2.0*coef_b/(1.0+coef_a))*(U_r-U_c)/dx_c;
	}
}

VEC3 Alg_LimSlope_vanAlbada(VEC3 U_c, VEC3 U_l, VEC3 U_r, DType dx_c, DType dx_l, DType dx_r, int k_va, DType coef_a, DType coef_b) {
	if ((dx_l==-1)&&(dx_r==-1))
		k_va = 0;			// Indicate uniform meshes
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;
	if ((k_va!=0)&&(k_va!=-1)&&(k_va<2)) {
		printf("Error in Alg_LimSlope_vanAlbada, k_va should be >2 or =-1\n");
		exit(8);
	}

	// phi(a)=b must be satisfied for second order accuracy
//	DType coef_a, coef_b;
	if (coef_a==-1)
		coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	if (coef_b==-1)
		coef_b = 2.0*dx_c/(dx_c+dx_r);

	if (k_va==0)
		k_va = 2;
	else if (k_va==-1)
		k_va = Alg_LimSlope_vanAlbada_calc_k(coef_a,coef_b);

	DType par_theta;		// This can be either the theta described in header
							// or the inverse of the theta described there
	int i;
	VEC3 slope;
	for (i=0; i<3; i++) {
		if (((U_c[i]-U_l[i])*(U_r[i]-U_c[i]))<0.0)
			slope[i] = 0.0;

		if (abs(U_c[i]-U_l[i])<abs(U_r[i]-U_c[i])) {
			par_theta = (U_c[i]-U_l[i])/(U_r[i]-U_c[i]);
			slope[i] = ((coef_b*(1.0+pow(par_theta,k_va-1)))
					/(coef_a+pow(par_theta,k_va)))*(U_c[i]-U_l[i])/dx_c;
		}
		else if (abs(U_c[i]-U_l[i])>abs(U_r[i]-U_c[i])) {
			par_theta = (U_r[i]-U_c[i])/(U_c[i]-U_l[i]);
			slope[i] = ((coef_b*(1.0+pow(par_theta,k_va-1)))
					/(1.0+coef_a*pow(par_theta,k_va)))*(U_r[i]-U_c[i])/dx_c;
		}
		else {
			if (U_r[i]==U_l[i])
				slope[i] = 0.0;
			else
				slope[i] = (2.0*coef_b/(1.0+coef_a))*(U_r[i]-U_c[i])/dx_c;
		}
	}
	return slope;
}

VEC4 Alg_LimSlope_vanAlbada(VEC4 U_c, VEC4 U_l, VEC4 U_r, DType dx_c, DType dx_l, DType dx_r, int k_va, DType coef_a, DType coef_b) {
	if ((dx_l==-1)&&(dx_r==-1))
		k_va = 0;			// Indicate uniform meshes
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;
	if ((k_va!=0)&&(k_va!=-1)&&(k_va<2)) {
		printf("Error in Alg_LimSlope_vanAlbada, k_va should be >2 or =-1\n");
		exit(8);
	}

	// phi(a)=b must be satisfied for second order accuracy
//	DType coef_a, coef_b;
	if (coef_a==-1)
		coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	if (coef_b==-1)
		coef_b = 2.0*dx_c/(dx_c+dx_r);

	if (k_va==0)
		k_va = 2;
	else if (k_va==-1)
		k_va = Alg_LimSlope_vanAlbada_calc_k(coef_a,coef_b);

	DType par_theta;		// This can be either the theta described in header
							// or the inverse of the theta described there
	int i;
	VEC4 slope;
	for (i=0; i<4; i++) {
		if (((U_c[i]-U_l[i])*(U_r[i]-U_c[i]))<0.0)
			slope[i] = 0.0;

		if (abs(U_c[i]-U_l[i])<abs(U_r[i]-U_c[i])) {
			par_theta = (U_c[i]-U_l[i])/(U_r[i]-U_c[i]);
			slope[i] = ((coef_b*(1.0+pow(par_theta,k_va-1)))
					/(coef_a+pow(par_theta,k_va)))*(U_c[i]-U_l[i])/dx_c;
		}
		else if (abs(U_c[i]-U_l[i])>abs(U_r[i]-U_c[i])) {
			par_theta = (U_r[i]-U_c[i])/(U_c[i]-U_l[i]);
			slope[i] = ((coef_b*(1.0+pow(par_theta,k_va-1)))
					/(1.0+coef_a*pow(par_theta,k_va)))*(U_r[i]-U_c[i])/dx_c;
		}
		else {
			if (U_r[i]==U_l[i])
				slope[i] = 0.0;
			else
				slope[i] = (2.0*coef_b/(1.0+coef_a))*(U_r[i]-U_c[i])/dx_c;
		}
	}
	return slope;
}

// Berger limiter
DType Alg_LimSlope_Berger1(DType U_c, DType U_l, DType U_r, DType dx_c, DType dx_l, DType dx_r) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	// phi(a)=b must be satisfied for second order accuracy
	DType coef_a, coef_b;
	coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	coef_b = 2.0*dx_c/(dx_c+dx_r);

	DType par_theta;
	if ((U_c-U_l)*(U_r-U_c)<=0.0)
		return 0.0;
	else {
		par_theta = (U_c-U_l)/(U_r-U_c);
		if (par_theta<=coef_a)
			return 2.0*((U_c-U_l)/dx_c)*(1.0-(1.0-0.5*coef_b/coef_a))*
				pow((par_theta/(1.0+par_theta))*((coef_a+1.0)/coef_a),
						coef_b/(2.0*coef_a-coef_b));
		else
			return 2.0*((U_r-U_c)/dx_c)*(1.0-(1.0-0.5*coef_b)*
					pow((1.0+coef_a)/(1.0+par_theta),coef_b/(2.0-coef_b)));
	}
}

VEC3 Alg_LimSlope_Berger1(VEC3 U_c, VEC3 U_l, VEC3 U_r, DType dx_c, DType dx_l, DType dx_r) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	// phi(a)=b must be satisfied for second order accuracy
	DType coef_a, coef_b;
	coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	coef_b = 2.0*dx_c/(dx_c+dx_r);

	DType par_theta;
	VEC3 slope;
	int i;
	for (i=0; i<3; i++) {
		if ((U_c[i]-U_l[i])*(U_r[i]-U_c[i])<=0.0)
			slope[i] = 0.0;
		else {
			par_theta = (U_c[i]-U_l[i])/(U_r[i]-U_c[i]);
			if (par_theta<=coef_a)
				slope[i] = 2.0*((U_c[i]-U_l[i])/dx_c)*(1.0-(1.0-0.5*coef_b/coef_a)*
					pow((par_theta/(1.0+par_theta))*((coef_a+1.0)/coef_a),
							coef_b/(2.0*coef_a-coef_b)));
			else
				slope[i] = 2.0*((U_r[i]-U_c[i])/dx_c)*(1.0-(1.0-0.5*coef_b)*
						pow((1.0+coef_a)/(1.0+par_theta),coef_b/(2.0-coef_b)));
		}
	}
	return slope;
}

VEC4 Alg_LimSlope_Berger1(VEC4 U_c, VEC4 U_l, VEC4 U_r, DType dx_c, DType dx_l, DType dx_r) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	// phi(a)=b must be satisfied for second order accuracy
	DType coef_a, coef_b;
	coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	coef_b = 2.0*dx_c/(dx_c+dx_r);

	DType par_theta;
	VEC4 slope;
	int i;
	for (i=0; i<4; i++) {
		if ((U_c[i]-U_l[i])*(U_r[i]-U_c[i])<=0.0)
			slope[i] = 0.0;
		else {
			par_theta = (U_c[i]-U_l[i])/(U_r[i]-U_c[i]);
			if (par_theta<=coef_a)
				slope[i] = 2.0*((U_c[i]-U_l[i])/dx_c)*(1.0-(1.0-0.5*coef_b/coef_a*
					pow((par_theta/(1.0+par_theta))*((coef_a+1.0)/coef_a),
							coef_b/(2.0*coef_a-coef_b))));
			else
				slope[i] = 2.0*((U_r[i]-U_c[i])/dx_c)*(1.0-(1.0-0.5*coef_b)*
						pow((1.0+coef_a)/(1.0+par_theta),coef_b/(2.0-coef_b)));
		}
	}
	return slope;
}



// Berger's second limiter
DType Alg_LimSlope_Berger2(DType U_c, DType U_l, DType U_r, DType dx_c, DType dx_l, DType dx_r) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	// phi(a)=b must be satisfied for second order accuracy
	DType coef_a, coef_b;
	coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	coef_b = 2.0*dx_c/(dx_c+dx_r);

	DType par_theta;
	if ((U_c-U_l)*(U_r-U_c)<=0.0)
		return 0.0;
	else {
		par_theta = (U_c-U_l)/(U_r-U_c);
		if (par_theta<=coef_a)
			return (coef_b/(1.0+coef_a))*((U_r-U_l)/dx_c)*(1.0-pow(
					1.0-(par_theta/(1.0+par_theta))*((1.0+coef_a)/coef_a),
					2.0*coef_a/coef_b));
		else
			return (coef_b/(1.0+coef_a))*((U_r-U_l)/dx_c)*(1.0-pow(
					1.0-(1.0+coef_a)/(1.0+par_theta),2.0/coef_b));
	}
}

VEC3 Alg_LimSlope_Berger2(VEC3 U_c, VEC3 U_l, VEC3 U_r, DType dx_c, DType dx_l, DType dx_r) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	// phi(a)=b must be satisfied for second order accuracy
	DType coef_a, coef_b;
	coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	coef_b = 2.0*dx_c/(dx_c+dx_r);

	DType par_theta;
	VEC3 slope;
	int i;
	for (i=0; i<3; i++) {
		if ((U_c[i]-U_l[i])*(U_r[i]-U_c[i])<=0.0)
			slope[i] = 0.0;
		else {
			par_theta = (U_c[i]-U_l[i])/(U_r[i]-U_c[i]);
			if (par_theta<=coef_a)
				slope[i] = (coef_b/(1.0+coef_a))*((U_r[i]-U_l[i])/dx_c)*(1.0-pow(
						1.0-(par_theta/(1.0+par_theta))*((1.0+coef_a)/coef_a),
						2.0*coef_a/coef_b));
			else
				slope[i] = (coef_b/(1.0+coef_a))*((U_r[i]-U_l[i])/dx_c)*(1.0-pow(
						1.0-(1.0+coef_a)/(1.0+par_theta),2.0/coef_b));
		}
	}
	return slope;
}

VEC4 Alg_LimSlope_Berger2(VEC4 U_c, VEC4 U_l, VEC4 U_r, DType dx_c, DType dx_l, DType dx_r) {
	if (dx_l==-1)
		dx_l = dx_c;
	if (dx_r==-1)
		dx_r = dx_c;

	// phi(a)=b must be satisfied for second order accuracy
	DType coef_a, coef_b;
	coef_a = (dx_c+dx_l)/(dx_c+dx_r);
	coef_b = 2.0*dx_c/(dx_c+dx_r);

	DType par_theta;
	VEC4 slope;
	int i;
	for (i=0; i<4; i++) {
		if ((U_c[i]-U_l[i])*(U_r[i]-U_c[i])<=0.0)
			slope[i] = 0.0;
		else {
			par_theta = (U_c[i]-U_l[i])/(U_r[i]-U_c[i]);
			if (par_theta<=coef_a)
				slope[i] = (coef_b/(1.0+coef_a))*((U_r[i]-U_l[i])/dx_c)*(1.0-pow(
						1.0-(par_theta/(1.0+par_theta))*((1.0+coef_a)/coef_a),
						2.0*coef_a/coef_b));
			else
				slope[i] = (coef_b/(1.0+coef_a))*((U_r[i]-U_l[i])/dx_c)*(1.0-pow(
						1.0-(1.0+coef_a)/(1.0+par_theta),2.0/coef_b));
		}
	}
	return slope;
}



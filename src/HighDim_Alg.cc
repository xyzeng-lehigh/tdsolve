#include "HighDim_Alg.h"

HVEC operator+ (const HVEC& vec_1, const HVEC& vec_2) {
	int dim = vec_1.size();
	if (dim!=vec_2.size()) {
		fprintf(stderr,"Error: HVEC size not the same in + oper.\n");
		return HVEC(dim,0.0);}
	HVEC vec(vec_1);
	for (int d=0; d<dim; d++) vec[d] += vec_2[d];
	return vec;
}

HVEC operator+ (const HVEC& vec_1) {
	return HVEC(vec_1);
}

HVEC operator+ (const HVEC& vec_in, const DType& sca_in) {
	int dim = vec_in.size();
	HVEC vec(vec_in);
	for (int d=0; d<dim; d++) vec[d] = sca_in+vec[d];
	return vec;
}

HVEC operator+ (const DType& sca_in, const HVEC& vec_in) {
	int dim = vec_in.size();
	HVEC vec(vec_in);
	for (int d=0; d<dim; d++) vec[d] = sca_in+vec[d];
	return vec;
}

HVEC operator- (const HVEC& vec_1, const HVEC& vec_2) {
	int dim = vec_1.size();
	if (dim!=vec_2.size()) {
		fprintf(stderr,"Error: HVEC size not the same in - oper.\n");
		return HVEC(dim,0.0);}
	HVEC vec(vec_1);
	for (int d=0; d<dim; d++) vec[d] -= vec_2[d];
	return vec;
}

HVEC operator- (const HVEC& vec_in) {
	int dim = vec_in.size();
	HVEC vec(vec_in);
	for (int d=0; d<dim; d++) vec[d] = -vec[d];
	return vec;
}

HVEC operator- (const HVEC& vec_in, const DType& sca_in) {
	int dim = vec_in.size();
	HVEC vec(vec_in);
	for (int d=0; d<dim; d++) vec[d] -= sca_in;
	return vec;
}

HVEC operator- (const DType& sca_in, const HVEC& vec_in) {
	int dim = vec_in.size();
	HVEC vec(vec_in);
	for (int d=0; d<dim; d++) vec[d] = sca_in-vec_in[d];
	return vec;
}

HVEC operator* (const DType& sca_in, const HVEC& vec_in) {
	int dim = vec_in.size();
	HVEC vec(vec_in);
	for (int d=0; d<dim; d++) vec[d] = sca_in*vec[d];
	return vec;
}

HVEC operator* (const HVEC& vec_in, const DType& sca_in) {
	int dim = vec_in.size();
	HVEC vec(vec_in);
	for (int d=0; d<dim; d++) vec[d] = sca_in*vec[d];
	return vec;
}

HVEC operator/ (const HVEC& vec_in, const DType& sca_in) {
	int dim = vec_in.size();
	HVEC vec(vec_in);
	for (int d=0; d<dim; d++) vec[d] = vec_in[d]/sca_in;
	return vec;
}


// This is the inner product
DType operator* (const HVEC& vec_1, const HVEC& vec_2) {
	int dim = vec_1.size();
	if (dim!=vec_2.size()) {
		fprintf(stderr,"Error: HVEC size not the same in * oper.\n");
		return 0.0;}
	DType sca = 0.0;
	for (int d=0; d<dim; d++)
		sca += vec_1[d]*vec_2[d];
	return sca;
}

// This is the exterior (cross) product (dim=3 required)
HVEC operator^ (const HVEC& vec_1, const HVEC& vec_2) {
	if ((vec_1.size()!=3)||(vec_2.size()!=3)) {
		fprintf(stderr,"Error: exterior (cross) product only available for HVEC with dim=3.\n");
		return HVEC(3,0.0);}
	HVEC vec(3);
	vec[0] = vec_1[1]*vec_2[2]-vec_1[2]*vec_2[1];
	vec[1] = vec_1[2]*vec_2[0]-vec_1[0]*vec_2[2];
	vec[2] = vec_1[0]*vec_2[1]-vec_1[1]*vec_2[0];
	return vec;
}

DType TriDMAT::det(void) {
    DType det_pre = 1.0;
    DType det_cur = mdiag[0];
    for (int i=1; i<dim; i++) {
        DType det_tmp = mdiag[i]*det_cur-upper[i-1]*lower[i-1]*det_pre;
        det_pre = det_cur;
        det_cur = det_tmp;
    }
    return det_cur;
}

HMAT TriDMAT::inv(void) {
    HMAT mat(dim,0.0);
    DType* dets = NULL;
    dets = new DType[dim+1];
    dets[0] = 1.0;
    dets[1] = mdiag[0];
    for (int i=2; i<dim+1; i++)
        dets[i] = mdiag[i-1]*dets[i-1]-upper[i-2]*lower[i-2]*dets[i-2];
    if (dets[dim]==0.0) {
        if (dets!=NULL) delete[] dets;
        return mat;
    }
    DType inv_det = 1.0/dets[dim];
    DType* phis = NULL;
    phis = new DType[dim+1];
    phis[dim]   = 1.0;
    phis[dim-1] = mdiag[dim-1];
    for (int i=dim-2; i>-1; i--)
        phis[i] = mdiag[i]*phis[i+1]-upper[i]*lower[i]*phis[i+2];
    for (int i=0; i<dim; i++)
        mat[i][i] = inv_det*dets[i]*phis[i+1];
    DType scal_i, scal_j;
    for (int i=0; i<dim-1; i++) {
        scal_i = inv_det*dets[i];
        scal_j = 1.0;
        for (int j=i+1; j<dim; j++) {
            scal_j *= -upper[j-1];
            mat[i][j] = scal_i*scal_j*phis[j+1];
        }
    }
    for (int i=1; i<dim; i++) {
        scal_i = inv_det*phis[i+1];
        scal_j = 1.0;
        for (int j=i-1; j>-1; j--) {
            scal_j *= -lower[j];
            mat[i][j] = scal_i*scal_j*dets[j];
        }
    }
    if (dets!=NULL) delete[] dets;
    if (phis!=NULL) delete[] phis;
    return mat;
}

TriDMAT operator+ (const TriDMAT& mat_1, const TriDMAT& mat_2) {
    int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: TriDMAT size not the same in + oper.\n");
		return TriDMAT(dim,0.0);}
    TriDMAT mat(mat_1);
    for (int i=0; i<dim-1; i++)
        mat.upper[i] += mat_2.upper[i];
    for (int i=0; i<dim-1; i++)
        mat.lower[i] += mat_2.lower[i];
    for (int i=0; i<dim; i++)
        mat.mdiag[i] += mat_2.mdiag[i];
    return mat;
}

TriDMAT operator+ (const TriDMAT& mat_in) {
	return TriDMAT(mat_in);
}

TriDMAT operator+ (const TriDMAT& mat_in, const DType sca_in) {
    int dim = mat_in.size();
    TriDMAT mat(mat_in);
    for (int i=0; i<dim-1; i++)
        mat.upper[i] += sca_in;
    for (int i=0; i<dim-1; i++)
        mat.lower[i] += sca_in;
    for (int i=0; i<dim; i++)
        mat.mdiag[i] += sca_in;
    return mat;
}

TriDMAT operator+ (const DType sca_in, const TriDMAT& mat_in) {
    int dim = mat_in.size();
    TriDMAT mat(mat_in);
    for (int i=0; i<dim-1; i++)
        mat.upper[i] += sca_in;
    for (int i=0; i<dim-1; i++)
        mat.lower[i] += sca_in;
    for (int i=0; i<dim; i++)
        mat.mdiag[i] += sca_in;
    return mat;
}

TriDMAT operator- (const TriDMAT& mat_1, const TriDMAT& mat_2) {
    int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: TriDMAT size not the same in - oper.\n");
		return TriDMAT(dim,0.0);}
    TriDMAT mat(mat_1);
    for (int i=0; i<dim-1; i++)
        mat.upper[i] -= mat_2.upper[i];
    for (int i=0; i<dim-1; i++)
        mat.lower[i] -= mat_2.lower[i];
    for (int i=0; i<dim; i++)
        mat.mdiag[i] -= mat_2.mdiag[i];
    return mat;
}

TriDMAT operator- (const TriDMAT& mat_in) {
    int dim = mat_in.size();
    TriDMAT mat(dim,0.0);
    for (int i=0; i<dim-1; i++)
        mat.upper[i] = -mat_in.upper[i];
    for (int i=0; i<dim-1; i++)
        mat.lower[i] = -mat_in.lower[i];
    for (int i=0; i<dim; i++)
        mat.mdiag[i] = -mat_in.mdiag[i];
    return mat;
}

TriDMAT operator- (const TriDMAT& mat_in, const DType sca_in) {
    int dim = mat_in.size();
    TriDMAT mat(mat_in);
    for (int i=0; i<dim-1; i++)
        mat.upper[i] -= sca_in;
    for (int i=0; i<dim-1; i++)
        mat.lower[i] -= sca_in;
    for (int i=0; i<dim; i++)
        mat.mdiag[i] -= sca_in;
    return mat;
}

TriDMAT operator- (const DType sca_in, const TriDMAT& mat_in) {
    int dim = mat_in.size();
    TriDMAT mat(dim,0.0);
    for (int i=0; i<dim-1; i++)
        mat.upper[i] = sca_in - mat_in.upper[i];
    for (int i=0; i<dim-1; i++)
        mat.lower[i] = sca_in - mat_in.lower[i];
    for (int i=0; i<dim; i++)
        mat.mdiag[i] = sca_in - mat_in.mdiag[i];
    return mat;
}

TriDMAT operator* (const TriDMAT& mat_in, const DType sca_in) {
    int dim = mat_in.size();
    TriDMAT mat(mat_in);
    for (int i=0; i<dim-1; i++)
        mat.upper[i] *= sca_in;
    for (int i=0; i<dim-1; i++)
        mat.lower[i] *= sca_in;
    for (int i=0; i<dim; i++)
        mat.mdiag[i] *= sca_in;
    return mat;
}

TriDMAT operator* (const DType sca_in, const TriDMAT& mat_in) {
    int dim = mat_in.size();
    TriDMAT mat(mat_in);
    for (int i=0; i<dim-1; i++)
        mat.upper[i] *= sca_in;
    for (int i=0; i<dim-1; i++)
        mat.lower[i] *= sca_in;
    for (int i=0; i<dim; i++)
        mat.mdiag[i] *= sca_in;
    return mat;
}

TriDMAT operator/ (const TriDMAT& mat_in, const DType sca_in) {
    int dim = mat_in.size();
    TriDMAT mat(mat_in);
    for (int i=0; i<dim-1; i++)
        mat.upper[i] /= sca_in;
    for (int i=0; i<dim-1; i++)
        mat.lower[i] /= sca_in;
    for (int i=0; i<dim; i++)
        mat.mdiag[i] /= sca_in;
    return mat;
}

HVEC operator* (const TriDMAT& mat_in, const HVEC& vec_in) {
	int dim = vec_in.size();
	if (dim!=mat_in.size()) {
		fprintf(stderr,"Error: TriDMAT and HVEC sizes not the same in * oper.\n");
		return HVEC(dim,0.0);}
	HVEC vec(dim,0.0);
    for (int i=0; i<dim; i++)
        vec[i] = mat_in.mdiag[i]*vec_in[i];
    for (int i=0; i<dim-1; i++) {
        vec[i]   += mat_in.upper[i]*vec_in[i+1];
        vec[i+1] += mat_in.lower[i]*vec_in[i];
    }
	return vec;
}

TriDMAT operator* (const TriDMAT& mat_1, const TriDMAT& mat_2) {
	int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: TriDMAT and TriDMAT sizes not the same in * oper.\n");
		return TriDMAT(dim,0.0);}
	TriDMAT mat(dim,0.0);
    for (int i=0; i<dim; i++)
        mat.mdiag[i] = mat_1.mdiag[i]*mat_2.mdiag[i];
    for (int i=0; i<dim-1; i++) {
        mat.mdiag[i]   += mat_1.upper[i]*mat_2.lower[i];
        mat.mdiag[i+1] += mat_1.lower[i]*mat_2.upper[i];
    }
    for (int i=0; i<dim-1; i++) {
        mat.upper[i] = mat_1.mdiag[i]*mat_2.upper[i] +
                       mat_1.upper[i]*mat_2.mdiag[i+1];
        mat.lower[i] = mat_1.lower[i]*mat_2.mdiag[i] +
                       mat_1.mdiag[i+1]*mat_2.lower[i];
    }
	return mat;
}

UTriMAT operator+ (const UTriMAT& mat_1, const UTriMAT& mat_2) {
    int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: UTriMAT size not the same in + oper.\n");
		return UTriMAT(dim,0.0);}
    int len = dim*(dim+1)/2;
    UTriMAT mat(mat_1);
    for (int i=0; i<len; i++)
        mat.data[i] += mat_2.data[i];
    return mat;
}

UTriMAT operator+ (const UTriMAT& mat_in) {
    return UTriMAT(mat_in);
}

UTriMAT operator+ (const UTriMAT& mat_in, const DType sca_in) {
    int dim = mat_in.size();
    UTriMAT mat(mat_in);
    int len = dim*(dim+1)/2;
    for (int i=0; i<len; i++)
        mat.data[i] += sca_in;
    return mat;
}

UTriMAT operator+ (const DType sca_in, const UTriMAT& mat_in) {
    int dim = mat_in.size();
    UTriMAT mat(mat_in);
    int len = dim*(dim+1)/2;
    for (int i=0; i<len; i++)
        mat.data[i] += sca_in;
    return mat;
}

UTriMAT operator- (const UTriMAT& mat_1, const UTriMAT& mat_2) {
    int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: UTriMAT size not the same in + oper.\n");
		return UTriMAT(dim,0.0);}
    int len = dim*(dim+1)/2;
    UTriMAT mat(mat_1);
    for (int i=0; i<len; i++)
        mat.data[i] -= mat_2.data[i];
    return mat;
}

UTriMAT operator- (const UTriMAT& mat_in) {
    int dim = mat_in.size();
    UTriMAT mat(dim,0.0);
    int len = dim*(dim+1)/2;
    for (int i=0; i<len; i++) mat.data[i] = -mat_in.data[i];
    return mat;
}

UTriMAT operator- (const UTriMAT& mat_in, const DType sca_in) {
    int dim = mat_in.size();
    UTriMAT mat(mat_in);
    int len = dim*(dim+1)/2;
    for (int i=0; i<len; i++)
        mat.data[i] -= sca_in;
    return mat;
}

UTriMAT operator- (const DType sca_in, const UTriMAT& mat_in) {
    int dim = mat_in.size();
    UTriMAT mat(dim,0.0);
    int len = dim*(dim+1)/2;
    for (int i=0; i<len; i++)
        mat.data[i] = sca_in - mat_in.data[i];
    return mat;
}

UTriMAT operator* (const UTriMAT& mat_in, const DType sca_in) {
    int dim = mat_in.size();
    UTriMAT mat(mat_in);
    int len = dim*(dim+1)/2;
    for (int i=0; i<len; i++)
        mat.data[i] *= sca_in;
    return mat;
}

UTriMAT operator* (const DType sca_in, const UTriMAT& mat_in) {
    int dim = mat_in.size();
    UTriMAT mat(mat_in);
    int len = dim*(dim+1)/2;
    for (int i=0; i<len; i++)
        mat.data[i] *= sca_in;
    return mat;
}

UTriMAT operator/ (const UTriMAT& mat_in, const DType sca_in) {
    int dim = mat_in.size();
    UTriMAT mat(mat_in);
    int len = dim*(dim+1)/2;
    for (int i=0; i<len; i++)
        mat.data[i] /= sca_in;
    return mat;
}

HVEC operator* (const UTriMAT& mat_in, const HVEC& vec_in) {
    int dim = vec_in.size();
	if (dim!=mat_in.size()) {
		fprintf(stderr,"Error: UTriMAT and HVEC sizes not the same in * oper.\n");
		return HVEC(dim,0.0);}
    HVEC vec(dim,0.0);
    for (int i=0; i<dim; i++) {
        for (int j=i; j<dim; j++)
            vec[i] += mat_in.data[mat_in.index(i,j)] * vec_in[j];
    }
    return vec;
}

UTriMAT operator* (const UTriMAT& mat_1, const UTriMAT& mat_2) {
	int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: UTriMAT and UTriMAT sizes not the same in * oper.\n");
		return UTriMAT(dim,0.0);}
    UTriMAT mat(dim,0.0);
    for (int i=0; i<dim; i++)
        for (int j=i; j<dim; j++) {
            DType& val = mat.data[mat.index(i,j)];
            for (int k=i; k<=j; k++)
                val += mat_1.data[mat_1.index(i,k)] *
                       mat_2.data[mat_2.index(k,j)];
        }
    return mat;
}

UTriMAT operator^ (const UTriMAT& mat_1, const UTriMAT& mat_2) {
    int dim1 = mat_1.size();
    int dim2 = mat_2.size();
    UTriMAT mat(dim1*dim2);
    for (int i1=0; i1<dim1; i1++) {
        int I = i1*dim2;
        for (int j1=i1; j1<dim1; j1++) {
            int J = j1*dim2;
            DType sca = mat_1.data[mat_1.index(i1,j1)];
            for (int i2=0; i2<dim2; i2++)
                for (int j2=i2; j2<dim2; j2++)
                    mat.data[mat.index(I+i2,J+j2)] = sca*mat_2.data[mat_2.index(i2,j2)];
        }
    }
    return mat;
}

void UTriMAT::print(void) {
    HMAT mat(*this);
    mat.print();
    return;
}

LTriMAT UTriMAT::trans(void) {
    LTriMAT mat(dim,0.0);
    for (int i=0; i<dim; i++)
        for (int j=i; j<dim; j++)
            mat.data[mat.index(j,i)] = data[index(i,j)];
    return mat;
}

UTriMAT UTriMAT::inv(void) {
    UTriMAT mat(dim,0.0);
    for (int i=0; i<dim; i++) {
        DType val = data[index(i,i)];
        if (val==0.0) {
            mat = 0.0;
            return mat;
        }
        else
            mat.data[index(i,i)] = 1.0/val;
    }
    for (int i=dim-2; i>=0; i--) {
        for (int j=i+1; j<dim; j++) {
            DType rhs = 0.0;
            for (int k=i+1; k<=j; k++)
                rhs += data[index(i,k)]*mat.data[index(k,j)];
            mat.data[index(i,j)] = -rhs/data[index(i,i)];
        }
    }
    return mat;
}

void UTriMAT::expand(DType val) {
    if (data==NULL) {
        fprintf(stderr,"Error: UTriMAT cannot expand a null matrix ... Exiting ...\n");
        exit(-1);
    }
    DType* data_tmp = new DType[len];
    for (int i=0; i<len; i++)
        data_tmp[i] = data[i];
    dim++;
    len = dim*(dim+1)/2;
    delete[] data;
    data = new DType[len];
    for (int i=0; i<dim-1; i++)         // old index
        for (int j=i; j<dim-1; j++) {   // old index
            data[index(i,j)] = data_tmp[index(i,j)-i];
        }
    for (int i=0; i<dim; i++)           // new column
        data[index(i,dim-1)] = val;
    delete[] data_tmp;
    return;
}

UTriMAT UTriMAT::sub(int sub_dim) {
    if (sub_dim<=0) sub_dim = dim;
    sub_dim = min(sub_dim,dim);
    UTriMAT mat(sub_dim,0.0);
    for (int i=0; i<sub_dim; i++)
        for (int j=i; j<sub_dim; j++)
            mat.data[mat.index(i,j)] = data[index(i,j)];
    return mat;
}

LTriMAT operator+ (const LTriMAT& mat_1, const LTriMAT& mat_2) {
    int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: LTriMAT size not the same in + oper.\n");
		return LTriMAT(dim,0.0);}
    int len = dim*(dim+1)/2;
    LTriMAT mat(mat_1);
    for (int i=0; i<len; i++)
        mat.data[i] += mat_2.data[i];
    return mat;
}

LTriMAT operator+ (const LTriMAT& mat_in) {
    return LTriMAT(mat_in);
}

LTriMAT operator+ (const LTriMAT& mat_in, const DType sca_in) {
    int dim = mat_in.size();
    LTriMAT mat(mat_in);
    int len = dim*(dim+1)/2;
    for (int i=0; i<len; i++)
        mat.data[i] += sca_in;
    return mat;
}

LTriMAT operator+ (const DType sca_in, const LTriMAT& mat_in) {
    int dim = mat_in.size();
    LTriMAT mat(mat_in);
    int len = dim*(dim+1)/2;
    for (int i=0; i<len; i++)
        mat.data[i] += sca_in;
    return mat;
}

LTriMAT operator- (const LTriMAT& mat_1, const LTriMAT& mat_2) {
    int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: LTriMAT size not the same in + oper.\n");
		return LTriMAT(dim,0.0);}
    int len = dim*(dim+1)/2;
    LTriMAT mat(mat_1);
    for (int i=0; i<len; i++)
        mat.data[i] -= mat_2.data[i];
    return mat;
}

LTriMAT operator- (const LTriMAT& mat_in) {
    int dim = mat_in.size();
    LTriMAT mat(dim,0.0);
    int len = dim*(dim+1)/2;
    for (int i=0; i<len; i++) mat.data[i] = -mat_in.data[i];
    return mat;
}

LTriMAT operator- (const LTriMAT& mat_in, const DType sca_in) {
    int dim = mat_in.size();
    LTriMAT mat(mat_in);
    int len = dim*(dim+1)/2;
    for (int i=0; i<len; i++)
        mat.data[i] -= sca_in;
    return mat;
}

LTriMAT operator- (const DType sca_in, const LTriMAT& mat_in) {
    int dim = mat_in.size();
    LTriMAT mat(dim,0.0);
    int len = dim*(dim+1)/2;
    for (int i=0; i<len; i++)
        mat.data[i] = sca_in - mat_in.data[i];
    return mat;
}

LTriMAT operator* (const LTriMAT& mat_in, const DType sca_in) {
    int dim = mat_in.size();
    LTriMAT mat(mat_in);
    int len = dim*(dim+1)/2;
    for (int i=0; i<len; i++)
        mat.data[i] *= sca_in;
    return mat;
}

LTriMAT operator* (const DType sca_in, const LTriMAT& mat_in) {
    int dim = mat_in.size();
    LTriMAT mat(mat_in);
    int len = dim*(dim+1)/2;
    for (int i=0; i<len; i++)
        mat.data[i] *= sca_in;
    return mat;
}

LTriMAT operator/ (const LTriMAT& mat_in, const DType sca_in) {
    int dim = mat_in.size();
    LTriMAT mat(mat_in);
    int len = dim*(dim+1)/2;
    for (int i=0; i<len; i++)
        mat.data[i] /= sca_in;
    return mat;
}

HVEC operator* (const LTriMAT& mat_in, const HVEC& vec_in) {
    int dim = vec_in.size();
	if (dim!=mat_in.size()) {
		fprintf(stderr,"Error: LTriMAT and HVEC sizes not the same in * oper.\n");
		return HVEC(dim,0.0);}
    HVEC vec(dim,0.0);
    for (int i=0; i<dim; i++) {
        for (int j=0; j<=i; j++)
            vec[i] += mat_in.data[mat_in.index(i,j)] * vec_in[j];
    }
    return vec;
}

LTriMAT operator* (const LTriMAT& mat_1, const LTriMAT& mat_2) {
	int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: LTriMAT and LTriMAT sizes not the same in * oper.\n");
		return LTriMAT(dim,0.0);}
    LTriMAT mat(dim,0.0);
    for (int i=0; i<dim; i++)
        for (int j=0; j<=i; j++) {
            DType& val = mat.data[mat.index(i,j)];
            for (int k=j; k<=i; k++)
                val += mat_1.data[mat_1.index(i,k)] *
                       mat_2.data[mat_2.index(k,j)];
        }
    return mat;
}

LTriMAT operator^ (const LTriMAT& mat_1, const LTriMAT& mat_2) {
    int dim1 = mat_1.size();
    int dim2 = mat_2.size();
    LTriMAT mat(dim1*dim2);
    for (int i1=0; i1<dim1; i1++) {
        int I = i1*dim2;
        for (int j1=0; j1<=i1; j1++) {
            int J = j1*dim2;
            DType sca = mat_1.data[mat_1.index(i1,j1)];
            for (int i2=0; i2<dim2; i2++)
                for (int j2=0; j2<=i2; j2++)
                    mat.data[mat.index(I+i2,J+j2)] = sca*mat_2.data[mat_2.index(i2,j2)];
        }
    }
    return mat;
}

void LTriMAT::print(void) {
    HMAT mat(*this);
    mat.print();
    return;
}

UTriMAT LTriMAT::trans(void) {
    UTriMAT mat(dim,0.0);
    for (int i=0; i<dim; i++)
        for (int j=0; j<=i; j++)
            mat.data[mat.index(j,i)] = data[index(i,j)];
    return mat;
}

LTriMAT LTriMAT::inv(void) {
    LTriMAT mat(dim,0.0);
    for (int i=0; i<dim; i++) {
        DType val = data[index(i,i)];
        if (val==0.0) {
            mat = 0.0;
            return mat;
        }
        else
            mat.data[index(i,i)] = 1.0/val;
    }
    for (int i=1; i<dim; i++) {
        for (int j=i-1; j>=0; j--) {
            DType rhs = 0.0;
            for (int k=j+1; k<=i; k++)
                rhs += mat.data[index(i,k)]*data[index(k,j)];
            mat.data[index(i,j)] = -rhs/data[index(j,j)];
        }
    }
    return mat;
}

void LTriMAT::expand(DType val) {
    if (data==NULL) {
        fprintf(stderr,"Error: LTriMAT cannot expand a null matrix ... Exiting ...\n");
        exit(-1);
    }
    DType* data_tmp = new DType[len];
    for (int i=0; i<len; i++)
        data_tmp[i] = data[i];
    dim++;
    int len_pre = len;
    len = dim*(dim+1)/2;
    delete[] data;
    data = new DType[len];
    for (int i=0; i<len_pre; i++)   data[i] = data_tmp[i];
    for (int i=len_pre; i<len; i++) data[i] = val;
    delete[] data_tmp;
    return;
}

LTriMAT LTriMAT::sub(int sub_dim) {
    if (sub_dim<=0) sub_dim = dim;
    sub_dim = min(sub_dim,dim);
    LTriMAT mat(sub_dim,0.0);
    int sub_len = sub_dim*(sub_dim+1)/2;
    for (int i=0; i<sub_len; i++)
        mat.data[i] = data[i];
    return mat;
}

/*
 * HMAT
 */

HMAT::HMAT(const HTENS2& tens_in) {
    dim = tens_in.size(1);
    assert( dim == tens_in.size(2) );
    data = new HVEC*[dim];
    for (int d=0; d<dim; d++)
        data[d] = new HVEC(dim);
    for (int i=0; i<dim; i++)
        for (int j=0; j<dim; j++)
            (*data[i])[j] = tens_in[i][j];
}

HMAT& HMAT::operator=(const HTENS2& tens_in) {
    assert( dim == tens_in.size(1) );
    assert( dim == tens_in.size(2) );
    for (int i=0; i<dim; i++)
        for (int j=0; j<dim; j++)
            (*data[i])[j] = tens_in[i][j];
    return *this;
}

HMAT operator+ (const HMAT& mat_1, const HMAT& mat_2) {
	int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: HMAT size not the same in + oper.\n");
		return HMAT(dim,0.0);}
	HMAT mat(mat_1);
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] += mat_2[i][j];
	return mat;
}

HMAT operator+ (const HMAT& mat_in) {
	return HMAT(mat_in);
}

HMAT operator+ (const HMAT& mat_in, const DType& sca_in) {
	int dim = mat_in.size();
	HMAT mat(mat_in);
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] += sca_in;
	return mat;
}

HMAT operator+ (const DType& sca_in, const HMAT& mat_in) {
	int dim = mat_in.size();
	HMAT mat(mat_in);
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] += sca_in;
	return mat;
}

HMAT operator- (const HMAT& mat_1, const HMAT& mat_2) {
	int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: HMAT size not the same in - oper.\n");
		return HMAT(dim,0.0);}
	HMAT mat(mat_1);
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] -= mat_2[i][j];
	return mat;
}

HMAT operator- (const HMAT& mat_in) {
	int dim = mat_in.size();
	HMAT mat(dim);
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] = -mat_in[i][j];
	return mat;
}

HMAT operator- (const HMAT& mat_in, const DType& sca_in) {
	int dim = mat_in.size();
	HMAT mat(mat_in);
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] -= sca_in;
	return mat;
}

HMAT operator- (const DType& sca_in, const HMAT& mat_in) {
	int dim = mat_in.size();
	HMAT mat(dim);
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] = sca_in-mat_in[i][j];
	return mat;
}

HMAT operator* (const DType& sca_in, const HMAT& mat_in) {
	int dim = mat_in.size();
	HMAT mat(mat_in);
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] *= sca_in;
	return mat;
}

HMAT operator* (const HMAT& mat_in, const DType& sca_in) {
	int dim = mat_in.size();
	HMAT mat(mat_in);
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] *= sca_in;
	return mat;
}

HMAT operator/ (const HMAT& mat_in, const DType& sca_in) {
	int dim = mat_in.size();
	HMAT mat(mat_in);
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] /= sca_in;
	return mat;
}

HVEC operator* (const HMAT& mat_in, const HVEC& vec_in) {
	int dim = vec_in.size();
	if (dim!=mat_in.size()) {
		fprintf(stderr,"Error: HMAT and HVEC sizes not the same in * oper.\n");
		return HVEC(dim,0.0);}
	HVEC vec(dim,0.0);
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			vec[i] += mat_in[i][j]*vec_in[j];
	return vec;
}

HMAT operator* (const HMAT& mat_1, const HMAT& mat_2) {
	int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: HMAT and HMAT sizes not the same in * oper.\n");
		return HMAT(dim,0.0);}
	HMAT mat(dim,0.0);
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			for (int d=0; d<dim; d++)
				mat[i][j] += mat_1[i][d]*mat_2[d][j];
	return mat;
}

// // This is the outer product
// HMAT operator| (const HVEC& vec_1, const HVEC& vec_2) {
// 	int dim = vec_1.size();
// 	if (dim!=vec_2.size()) {
// 		fprintf(stderr,"Error: HVEC and HVEC sizes not the same in | oper.\n");
// 		return HMAT(dim,0.0);}
// 	HMAT mat(dim);
// 	for (int i=0; i<dim; i++) {
// 		DType sca = vec_1[i];
// 		for (int j=0; j<dim; j++)
// 			mat[i][j] = sca*vec_2[j];
// 	}
// 	return mat;
// }

// This is the Kronecker product
HMAT operator^ (const HMAT& mat_1, const HMAT& mat_2) {
	int dim1 = mat_1.size();
	int dim2 = mat_2.size();
	HMAT mat(dim1*dim2);
	for (int i1=0; i1<dim1; i1++) {
		int I = i1*dim2;
		for (int j1=0; j1<dim1; j1++) {
			int J = j1*dim2;
			DType sca = mat_1[i1][j1];
			for (int i2=0; i2<dim2; i2++)
				for (int j2=0; j2<dim2; j2++)
					mat[I+i2][J+j2] = sca*mat_2[i2][j2];
		}
	}
	return mat;
}

HMAT operator+ (const HMAT& mat_1, const TriDMAT& mat_2) {
	int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: HMAT and TriDMAT size not the same in + oper.\n");
		return HMAT(dim,0.0);}
    HMAT mat(mat_1);
    for (int i=0; i<dim; i++)
        mat[i][i] += mat_2.mdiag[i];
    for (int i=0; i<dim-1; i++) {
        mat[i][i+1] += mat_2.upper[i];
        mat[i+1][i] += mat_2.lower[i];
    }
    return mat;
}

HMAT operator+ (const TriDMAT& mat_1, const HMAT& mat_2) {
	int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: TriDMAT and HMAT size not the same in + oper.\n");
		return HMAT(dim,0.0);}
    HMAT mat(mat_2);
    for (int i=0; i<dim; i++)
        mat[i][i] += mat_1.mdiag[i];
    for (int i=0; i<dim-1; i++) {
        mat[i][i+1] += mat_1.upper[i];
        mat[i+1][i] += mat_1.lower[i];
    }
    return mat;
}

HMAT operator- (const HMAT& mat_1, const TriDMAT& mat_2) {
	int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: HMAT and TriDMAT size not the same in - oper.\n");
		return HMAT(dim,0.0);}
    HMAT mat(mat_1);
    for (int i=0; i<dim; i++)
        mat[i][i] -= mat_2.mdiag[i];
    for (int i=0; i<dim-1; i++) {
        mat[i][i+1] -= mat_2.upper[i];
        mat[i+1][i] -= mat_2.lower[i];
    }
    return mat;
}

HMAT operator- (const TriDMAT& mat_1, const HMAT& mat_2) {
	int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: TriDMAT and HMAT size not the same in + oper.\n");
		return HMAT(dim,0.0);}
    HMAT mat(-mat_2);
    for (int i=0; i<dim; i++)
        mat[i][i] += mat_1.mdiag[i];
    for (int i=0; i<dim-1; i++) {
        mat[i][i+1] += mat_1.upper[i];
        mat[i+1][i] += mat_1.lower[i];
    }
    return mat;
}

HMAT operator* (const HMAT& mat_1, const TriDMAT& mat_2) {
	int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: HMAT and TriDMAT size not the same in * oper.\n");
		return HMAT(dim,0.0);}
    HMAT mat(dim,0.0);
    for (int i=0; i<dim; i++)
        for (int j=0; j<dim; j++) {
            if (j>0)
                mat[i][j] += mat_1[i][j-1]*mat_2.upper[j-1];
            mat[i][j] += mat_1[i][j]*mat_2.mdiag[j];
            if (j<dim-1)
                mat[i][j] += mat_1[i][j+1]*mat_2.lower[j];
        }
    return mat;
}

HMAT operator* (const TriDMAT& mat_1, const HMAT& mat_2) {
	int dim = mat_1.size();
	if (dim!=mat_2.size()) {
		fprintf(stderr,"Error: TriDMAT and HMAT size not the same in * oper.\n");
		return HMAT(dim,0.0);}
    HMAT mat(dim,0.0);
    for (int i=0; i<dim; i++)
        for (int j=0; j<dim; j++) {
            if (i>0)
                mat[i][j] += mat_1.lower[i-1]*mat_2[i-1][j];
            mat[i][j] += mat_1.mdiag[i]*mat_2[i][j];
            if (i<dim-1)
                mat[i][j] += mat_1.upper[i]*mat_2[i+1][j];
        }
    return mat;
}

HMAT operator^ (const HMAT& mat_1, const TriDMAT& mat_2) {
	int dim1 = mat_1.size();
	int dim2 = mat_2.size();
	HMAT mat(dim1*dim2,0.0);
    for (int i1=0; i1<dim1; i1++) {
        int I = i1*dim2;
        for (int j1=0; j1<dim1; j1++) {
            int J = j1*dim2;
            DType sca = mat_1[i1][j1];
            for (int k=1; k<dim2-1; k++) {
                mat[I+k][J+k+1] = sca*mat_2.lower[k];
                mat[I+k+1][J+k] = sca*mat_2.upper[k];
            }
            for (int k=1; k<dim2; k++)
                mat[I+k][I+k] = sca*mat_2.mdiag[k];
        }
    }
	return mat;
}

HMAT operator^ (const TriDMAT& mat_1, const HMAT& mat_2) {
	int dim1 = mat_1.size();
	int dim2 = mat_2.size();
	HMAT mat(dim1*dim2,0.0);
    for (int k=0; k<dim1-1; k++) {
        int I, J;
        I = k*dim2;
        J = (k+1)*dim2;
        DType sca_u = mat_1.upper[k];
        DType sca_l = mat_1.lower[k];
        for (int i2=0; i2<dim2; i2++)
            for (int j2=0; j2<dim2; j2++) {
                DType sca = mat_2[i2][j2];
                mat[I+i2][J+j2] = sca_u*sca;
                mat[J+i2][I+i2] = sca_l*sca;
            }
    }
    for (int k=0; k<dim1; k++) {
        int K = k*dim2;
        DType sca = mat_1.mdiag[k];
        for (int i2=0; i2<dim2; i2++)
            for (int j2=0; j2<dim2; j2++)
                mat[K+i2][K+j2] = sca*mat_2[i2][j2];
    }
	return mat;
}

HMAT operator^ (const TriDMAT& mat_1, const TriDMAT& mat_2) {
    int dim1 = mat_1.size();
    int dim2 = mat_2.size();
    HMAT mat(dim1*dim2,0.0);
    for (int i1=0; i1<dim1-1; i1++) {
        int I = i1*dim2;
        int J = (i1+1)*dim2;
        DType sca_u = mat_1.upper[i1];
        DType sca_l = mat_1.lower[i1];
        for (int i2=0; i2<dim2-1; i2++) {
            DType sca_u2 = mat_2.upper[i2];
            DType sca_l2 = mat_2.lower[i2];
            mat[I+i2][J+i2+1] = sca_u * sca_u2;
            mat[I+i2+1][J+i2] = sca_u * sca_l2;
            mat[J+i2][I+i2+1] = sca_l * sca_u2;
            mat[J+i2+1][I+i2] = sca_l * sca_l2;
        }
        for (int i2=0; i2<dim2; i2++) {
            DType sca = mat_2.mdiag[i2];
            mat[I+i2][J+i2] = sca_u*sca;
            mat[J+i2][I+i2] = sca_l*sca;
        }
    }
    for (int i1=0; i1<dim1; i1++) {
        int I = i1*dim2;
        DType sca_1 = mat_1.mdiag[i1];
        for (int i2=0; i2<dim2-1; i2++) {
            mat[I+i2][I+i2+1] = sca_1*mat_2.upper[i2];
            mat[I+i2+1][I+i2] = sca_1*mat_2.lower[i2];
        }
        for (int i2=0; i2<dim2; i2++)
            mat[I+i2][I+i2] = sca_1*mat_2.mdiag[i2];
    }
    return mat;
}

// Find the index of the entry with largest absolute value
void HMAT::findMatMax(int* r, int* c, HMAT& M, int I) {
	int dim_loc = M.size();
	int r_idx = I;
	int c_idx = I;
	DType max_val = abs(M[I][I]);
	for (int i=I; i<dim_loc; i++)
		for (int j=I; j<dim_loc; j++) {
			DType tmp_val = abs(M[i][j]);
			if (tmp_val>max_val) {
				r_idx = i; c_idx = j;
				max_val = tmp_val;
			}
		}
	(*r) = r_idx;
	(*c) = c_idx;
	return;
}

void HMAT::findRowMax(int* c, HVEC& V, int I) {
	int dim_loc = V.size();
	int c_idx = I;
	DType max_val = abs(V[I]);
	for (int j=I; j<dim_loc; j++) {
		DType tmp_val = abs(V[j]);
		if (tmp_val>max_val) {
			c_idx = j;
			max_val = tmp_val;
		}
	}
	(*c) = c_idx;
	return;
}

// Permute two rows or two columns
void HMAT::permRow(HMAT& M, int r1, int r2) {
	if (r1==r2) return;
	int dim_loc = M.size();
	DType* r_tmp;
	r_tmp = new DType[dim_loc];
	for (int j=0; j<dim_loc; j++) r_tmp[j] = M[r1][j];
	for (int j=0; j<dim_loc; j++) {
		M[r1][j] = M[r2][j];
		M[r2][j] = r_tmp[j];
	}
	delete[] r_tmp;
	return;
}

void HMAT::permCol(HMAT& M, int c1, int c2) {
	if (c1==c2) return;
	int dim_loc = M.size();
	DType* c_tmp;
	c_tmp = new DType[dim_loc];
	for (int i=0; i<dim_loc; i++) c_tmp[i] = M[i][c1];
	for (int i=0; i<dim_loc; i++) {
		M[i][c1] = M[i][c2];
		M[i][c2] = c_tmp[i];
	}
	delete[] c_tmp;
	return;
}

// Computed by PLU decomposition with largest entry pivoting
DType HMAT::det_stable(void) {
	HMAT M(*this);
	int r_max, c_max;
	DType det_val = 1.0;
	for (int i=0; i<dim; i++) {
		findMatMax(&r_max,&c_max,M,i);
		permRow(M,i,r_max);
		permCol(M,i,c_max);

		if (i!=r_max) det_val = -det_val;
		if (i!=c_max) det_val = -det_val;

		DType m11 = abs(M[i][i]);
		if (m11==0.0) {
			det_val = 0.0;
			break;
		}
		DType l = sqrt(m11);
		DType u = (M[i][i]>0.0) ? l : -l;
		det_val *= l*u;

		M[i][i] = u;
		for (int j=i+1; j<dim; j++) {
			M[i][j] = M[i][j]/l;
			M[j][i] = M[j][i]/u;
		}
		for (int k=i+1; k<dim; k++)
			for (int l=i+1; l<dim; l++)
				M[k][l] -= M[k][i]*M[i][l];
	}
	return det_val;
}

// Computed by PLU decomposition with non-zero pivoting
DType HMAT::det_fast(void) {
	HMAT M(*this);
	int c_max;
	DType det_val = 1.0;
	for (int i=0; i<dim; i++) {
		findRowMax(&c_max,M[i],i);
		permCol(M,i,c_max);

		if (i!=c_max) det_val = -det_val;

		if (M[i][i]==0.0) {
			det_val = 0.0;
			break;
		}
		DType l = 1.0;
		DType u = M[i][i]; 
		det_val *= u;
		for (int j=i+1; j<dim; j++) {
			M[i][j] = M[i][j]/l;
			M[j][i] = M[j][i]/u;
		}
		for (int k=i+1; k<dim; k++)
			for (int l=i+1; l<dim; l++)
				M[k][l] -= M[k][i]*M[i][l];
	}
	return det_val;
}

// Computed by PLU decomposition with largest entry pivoting
HMAT HMAT::inv_stable(void) {
	HMAT M(*this);
	int* Pr;	// Row permutation
	int* Pc;	// Column permutation
	Pr = new int[dim];
	Pc = new int[dim];
	for (int i=0; i<dim; i++) Pr[i] = i;
	for (int i=0; i<dim; i++) Pc[i] = i;
	int r_max, c_max;
	bool isSingular=false;
	// Swept one: compute PM=LU, store LU in M
	for (int i=0; i<dim; i++) {
		findMatMax(&r_max,&c_max,M,i);
		permRow(M,i,r_max); Pr[i] = r_max;
		permCol(M,i,c_max); Pc[i] = c_max;

		DType m11 = abs(M[i][i]);
		if (m11==0.0) {
			isSingular = true;
			M = 0.0;
			break;
		}
		DType l = sqrt(m11);
		DType u = (M[i][i]>0.0) ? l : -l;

		M[i][i] = u;
		for (int j=i+1; j<dim; j++) {
			M[i][j] = M[i][j]/l;
			M[j][i] = M[j][i]/u;
		}
		for (int k=i+1; k<dim; k++)
			for (int l=i+1; l<dim; l++)
				M[k][l] -= M[k][i]*M[i][l];
	}
	// Swept two: compute inverses of
	// 	U: start from last colume and moving leftward 
	// 	L: starf from last row and moving upward 
	if (!isSingular) 
		for (int i=(dim-1); i>-1; i--) {
			M[i][i] = 1.0/M[i][i];
			for (int j=i-1; j>-1; j--) {
				DType tmp = 0.0;
				for (int k=j+1; k<i+1; k++) tmp -= M[j][k]*M[k][i];
				M[j][i] = tmp/M[j][j];
			}
			DType scal = (M[i][i]>0) ? 1.0 : -1.0;
			M[i][i] = abs(M[i][i]);
			for (int j=i-1; j>-1; j--) {
				DType tmp = 0.0;
				for (int k=j+1; k<i+1; k++) tmp -= M[k][j]*M[i][k];
				M[i][j] = tmp/abs(M[j][j]);
			}
			M[i][i] = scal*M[i][i];
		}
	// Swept three: construct inverse of PM from inverses of L and U
	if (!isSingular)
		for (int i=0; i<dim; i++) {
			M[i][i] = abs(M[i][i])*M[i][i];
			for (int k=i+1; k<dim; k++)
				M[i][i] += M[i][k]*M[k][i];
			for (int j=i+1; j<dim; j++) {
				DType tmp_ij, tmp_ji;
				tmp_ij = 0.0;
				tmp_ji = 0.0;
				for (int k=j; k<dim; k++) {
					if (k==j)
						tmp_ij += M[i][k]*abs(M[k][j]);
					else
						tmp_ij += M[i][k]*M[k][j];
					tmp_ji += M[j][k]*M[k][i];
				}
				M[i][j] = tmp_ij;
				M[j][i] = tmp_ji;
			}
		}
	//	Finally, permute the resulting matrix M = MP
	if (!isSingular)
		for (int i=(dim-1); i>-1; i--) {
			int cj = Pr[i];
			int rj = Pc[i];
			permRow(M,i,rj);
			permCol(M,i,cj);
		}
	delete[] Pr;
	delete[] Pc;
	return M;
}

// Computed by PLU decomposition with non-zero pivoting
HMAT HMAT::inv_fast(void) {
	HMAT M(*this);
	int* Pr;	// Row permutation
	int* Pc;	// Column permutation
	Pr = new int[dim];
	Pc = new int[dim];
	for (int i=0; i<dim; i++) Pr[i] = i;
	for (int i=0; i<dim; i++) Pc[i] = i;
	int c_max;
	bool isSingular=false;
	// Swept one: compute PM=LU, store LU in M
	for (int i=0; i<dim; i++) {
		findRowMax(&c_max,M[i],i);
		permCol(M,i,c_max); Pc[i] = c_max;

		DType m11 = abs(M[i][i]);
		if (m11==0.0) {
			isSingular = true;
			M = 0.0;
			break;
		}
		DType l = sqrt(m11);
		DType u = (M[i][i]>0.0) ? l : -l;

		M[i][i] = u;
		for (int j=i+1; j<dim; j++) {
			M[i][j] = M[i][j]/l;
			M[j][i] = M[j][i]/u;
		}
		for (int k=i+1; k<dim; k++)
			for (int l=i+1; l<dim; l++)
				M[k][l] -= M[k][i]*M[i][l];
	}
	// Swept two: compute inverses of
	// 	U: start from last colume and moving leftward 
	// 	L: starf from last row and moving upward 
	if (!isSingular) 
		for (int i=(dim-1); i>-1; i--) {
			M[i][i] = 1.0/M[i][i];
			for (int j=i-1; j>-1; j--) {
				DType tmp = 0.0;
				for (int k=j+1; k<i+1; k++) tmp -= M[j][k]*M[k][i];
				M[j][i] = tmp/M[j][j];
			}
			DType scal = (M[i][i]>0) ? 1.0 : -1.0;
			M[i][i] = abs(M[i][i]);
			for (int j=i-1; j>-1; j--) {
				DType tmp = 0.0;
				for (int k=j+1; k<i+1; k++) tmp -= M[k][j]*M[i][k];
				M[i][j] = tmp/abs(M[j][j]);
			}
			M[i][i] = scal*M[i][i];
		}
	// Swept three: construct inverse of PM from inverses of L and U
	if (!isSingular)
		for (int i=0; i<dim; i++) {
			M[i][i] = abs(M[i][i])*M[i][i];
			for (int k=i+1; k<dim; k++)
				M[i][i] += M[i][k]*M[k][i];
			for (int j=i+1; j<dim; j++) {
				DType tmp_ij, tmp_ji;
				tmp_ij = 0.0;
				tmp_ji = 0.0;
				for (int k=j; k<dim; k++) {
					if (k==j)
						tmp_ij += M[i][k]*abs(M[k][j]);
					else
						tmp_ij += M[i][k]*M[k][j];
					tmp_ji += M[j][k]*M[k][i];
				}
				M[i][j] = tmp_ij;
				M[j][i] = tmp_ji;
			}
		}
	//	Finally, permute the resulting matrix M = MP
	if (!isSingular)
		for (int i=(dim-1); i>-1; i--) {
			int rj = Pc[i];
			permRow(M,i,rj);
		}
	delete[] Pr;
	delete[] Pc;
	return M;
}

/*
 * HTENS2
 */

HTENS2 operator+ (const HTENS2& tens_1, const HTENS2& tens_2) {
    int dim1 = tens_1.size(1);
    int dim2 = tens_1.size(2);
    assert( (dim1==tens_2.size(1)) && (dim2==tens_2.size(2)) );
    HTENS2 tens(tens_1);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
			tens[i][j] += tens_2[i][j];
	return tens;
}

HTENS2 operator+ (const HTENS2& tens_1, const DType& data_2) {
    int dim1 = tens_1.size(1);
    int dim2 = tens_1.size(2);
    HTENS2 tens(tens_1);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
			tens[i][j] += data_2;
	return tens;
}

HTENS2 operator+ (const DType& data_1, const HTENS2& tens_2) {
    int dim1 = tens_2.size(1);
    int dim2 = tens_2.size(2);
    HTENS2 tens(tens_2);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
			tens[i][j] += data_1;
	return tens;
}

HTENS2 operator+ (const HTENS2& tens_in) {
    return HTENS2(tens_in);
}

HTENS2 operator- (const HTENS2& tens_1, const HTENS2& tens_2) {
    int dim1 = tens_1.size(1);
    int dim2 = tens_1.size(2);
    assert( (dim1==tens_2.size(1)) && (dim2==tens_2.size(2)) );
    HTENS2 tens(tens_1);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
			tens[i][j] -= tens_2[i][j];
	return tens;
}

HTENS2 operator- (const HTENS2& tens_1, const DType& data_2) {
    int dim1 = tens_1.size(1);
    int dim2 = tens_1.size(2);
    HTENS2 tens(tens_1);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
			tens[i][j] -= data_2;
	return tens;
}

HTENS2 operator- (const DType& data_1, const HTENS2& tens_2) {
    int dim1 = tens_2.size(1);
    int dim2 = tens_2.size(2);
    HTENS2 tens(dim1,dim2);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
			tens[i][j] = data_1 - tens_2[i][j];
	return tens;
}

HTENS2 operator- (const HTENS2& tens_in) {
    int dim1 = tens_in.size(1);
    int dim2 = tens_in.size(2);
    HTENS2 tens(dim1,dim2);
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            tens[i][j] = -tens_in[i][j];
    return tens;
}

HTENS2 operator* (const HTENS2& tens_1, const DType& data_2) {
    int dim1 = tens_1.size(1);
    int dim2 = tens_1.size(2);
    HTENS2 tens(tens_1);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
			tens[i][j] *= data_2;
	return tens;
}

HTENS2 operator* (const DType& data_1, const HTENS2& tens_2) {
    int dim1 = tens_2.size(1);
    int dim2 = tens_2.size(2);
    HTENS2 tens(tens_2);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
			tens[i][j] *= data_1;
	return tens;
}

HTENS2 operator/ (const HTENS2& tens_1, const DType& data_2) {
    int dim1 = tens_1.size(1);
    int dim2 = tens_1.size(2);
    HTENS2 tens(tens_1);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
			tens[i][j] /= data_2;
	return tens;
}

HTENS2 operator* (const HTENS2& tens_1, const HTENS2& tens_2) {
    int dim1 = tens_1.size(1);
    int dim  = tens_1.size(2);
    assert( dim == tens_2.size(1) );
    int dim2 = tens_2.size(2);
    HTENS2 tens(dim1,dim2);
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            for (int k=0; k<dim; k++)
                tens[i][j] += tens_1[i][k]*tens_2[k][j];
    return tens;
}

HTENS2 operator* (const HTENS2& tens_1, const HMAT& mat_2) {
    int dim1 = tens_1.size(1);
    int dim2 = tens_1.size(2);
    assert( dim2 == mat_2.size() );
    HTENS2 tens(dim1,dim2);
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            for (int k=0; k<dim2; k++)
                tens[i][j] += tens_1[i][k]*mat_2[k][j];
    return tens;
}

HTENS2 operator* (const HMAT& mat_1, const HTENS2& tens_2) {
    int dim1 = mat_1.size();
    assert( dim1 == tens_2.size(1) );
    int dim2 = tens_2.size(2);
    HTENS2 tens(dim1,dim2);
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            for (int k=0; k<dim1; k++)
                tens[i][j] += mat_1[i][k]*tens_2[k][j];
    return tens;
}

HVEC operator* (const HTENS2& tens_1, const HVEC& vec_2) {
    int dim = tens_1.size(1);
    int dim2 = tens_1.size(2);
    assert( dim2 == vec_2.size() );
    HVEC vec(dim,0.0);
    for (int i=0; i<dim; i++)
        for (int j=0; j<dim2; j++)
            vec[i] += tens_1[i][j] * vec_2[j];
    return vec;
}

// This is the outer product
HTENS2 operator| (const HVEC& vec_1, const HVEC& vec_2) {
    int dim1 = vec_1.size();
    int dim2 = vec_2.size();
    HTENS2 tens(dim1,dim2);
    for (int i=0; i<dim1; i++) {
        DType val1 = vec_1[i];
        for (int j=0; j<dim2; j++)
            tens[i][j] = val1*vec_2[j];
    }
    return tens;
}

// This is the Kronecker product
HTENS2 operator^ (const HTENS2& tens_1, const HTENS2& tens_2) {
    int dim1_1 = tens_1.size(1);
    int dim1_2 = tens_1.size(2);
    int dim2_1 = tens_2.size(1);
    int dim2_2 = tens_2.size(2);
    HTENS2 tens(dim1_1*dim2_1,dim1_2*dim2_2);
    for (int i1=0; i1<dim1_1; i1++) {
        int I = i1*dim2_1;
        for (int j1=0; j1<dim1_2; j1++) {
            int J = j1*dim2_2;
            DType sca = tens_1[i1][j1];
            for (int i2=0; i2<dim2_1; i2++)
                for (int j2=0; j2<dim2_2; j2++)
                    tens[I+i2][J+j2] = sca*tens_2[i2][j2];
        }
    }
    return tens;
}

/*
 * HTENS3
 */

HTENS3 operator+ (const HTENS3& tens_1, const HTENS3& tens_2) {
    int dim1 = tens_1.size(1);
    int dim2 = tens_1.size(2);
    int dim3 = tens_1.size(3);
    assert( (dim1==tens_2.size(1)) && (dim2==tens_2.size(2)) && (dim3==tens_2.size(3)) );
    HTENS3 tens(tens_1);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
			    tens[i][j][k] += tens_2[i][j][k];
	return tens;
}

HTENS3 operator+ (const HTENS3& tens_1, const DType& data_2) {
    int dim1 = tens_1.size(1);
    int dim2 = tens_1.size(2);
    int dim3 = tens_1.size(3);
    HTENS3 tens(tens_1);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
			    tens[i][j][k] += data_2;
	return tens;
}

HTENS3 operator+ (const DType& data_1, const HTENS3& tens_2) {
    int dim1 = tens_2.size(1);
    int dim2 = tens_2.size(2);
    int dim3 = tens_2.size(3);
    HTENS3 tens(tens_2);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
			    tens[i][j][k] += data_1;
	return tens;
}

HTENS3 operator+ (const HTENS3& tens_in) {
    return HTENS3(tens_in);
}

HTENS3 operator- (const HTENS3& tens_1, const HTENS3& tens_2) {
    int dim1 = tens_1.size(1);
    int dim2 = tens_1.size(2);
    int dim3 = tens_1.size(3);
    assert( (dim1==tens_2.size(1)) && (dim2==tens_2.size(2)) && (dim3==tens_2.size(3)) );
    HTENS3 tens(tens_1);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
			    tens[i][j][k] -= tens_2[i][j][k];
	return tens;
}

HTENS3 operator- (const HTENS3& tens_1, const DType& data_2) {
    int dim1 = tens_1.size(1);
    int dim2 = tens_1.size(2);
    int dim3 = tens_1.size(3);
    HTENS3 tens(tens_1);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
			    tens[i][j][k] -= data_2;
	return tens;
}

HTENS3 operator- (const DType& data_1, const HTENS3& tens_2) {
    int dim1 = tens_2.size(1);
    int dim2 = tens_2.size(2);
    int dim3 = tens_2.size(3);
    HTENS3 tens(dim1,dim2,dim3);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
			    tens[i][j][k] = data_1 - tens_2[i][j][k];
	return tens;
}

HTENS3 operator- (const HTENS3& tens_in) {
    int dim1 = tens_in.size(1);
    int dim2 = tens_in.size(2);
    int dim3 = tens_in.size(3);
    HTENS3 tens(dim1,dim2,dim3);
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
                tens[i][j][k] = -tens_in[i][j][k];
    return tens;
}

HTENS3 operator* (const HTENS3& tens_1, const DType& data_2) {
    int dim1 = tens_1.size(1);
    int dim2 = tens_1.size(2);
    int dim3 = tens_1.size(3);
    HTENS3 tens(tens_1);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
			    tens[i][j][k] *= data_2;
	return tens;
}

HTENS3 operator* (const DType& data_1, const HTENS3& tens_2) {
    int dim1 = tens_2.size(1);
    int dim2 = tens_2.size(2);
    int dim3 = tens_2.size(3);
    HTENS3 tens(tens_2);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
			    tens[i][j][k] *= data_1;
	return tens;
}

HTENS3 operator/ (const HTENS3& tens_1, const DType& data_2) {
    int dim1 = tens_1.size(1);
    int dim2 = tens_1.size(2);
    int dim3 = tens_1.size(3);
    HTENS3 tens(tens_1);
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
			    tens[i][j][k] /= data_2;
	return tens;
}

HVEC operator* (const HTENS3& tens3, const HTENS2& tens2) {
    int dim  = tens3.size(1);
    int dim2 = tens3.size(2);
    int dim3 = tens3.size(3);
    assert( dim2 == tens2.size(1) );
    assert( dim3 == tens2.size(2) );
    HVEC vec(dim,0.0);
    for (int i=0; i<dim; i++)
        for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
                vec[i] += tens3[i][j][k] * tens2[j][k];
    return vec;
}

HTENS2 operator* (const HTENS3& tens3, const HVEC& vec) {
    int dim1 = tens3.size(1);
    int dim2 = tens3.size(2);
    int dim  = tens3.size(3);
    assert( dim == vec.size() );
    HTENS2 tens(dim1,dim2);
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            for (int k=0; k<dim; k++)
                tens[i][j] += tens3[i][j][k] * vec[k];
    return tens;
}


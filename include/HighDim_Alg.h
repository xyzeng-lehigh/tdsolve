#ifndef HIGHDIM_ALG_H
#define HIGHDIM_ALG_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include "Data_Type.h"
using namespace std;

class HVEC {
	private:
		int		dim;
		DType	*data;
	public:
		int size(void) const {return dim;}
		void	print(void) {
			for (int d=0; d<dim; d++) fprintf(stdout,"%10.20f\t",data[d]);
			fprintf(stdout,"\n"); 
      return;
		}

		HVEC(int in_dim) {dim = in_dim; data = new DType[dim];
			for (int d=0; d<dim; d++) data[d] = 0.0;}

		HVEC(int in_dim, DType in_val) {dim = in_dim; data = new DType[dim];
			for (int d=0; d<dim; d++) data[d] = in_val;}

		HVEC(const HVEC& in_vec) {
			dim = in_vec.size();	data = new DType[dim];
			for (int d=0; d<dim; d++) data[d] = in_vec[d];}

		~HVEC(void) {if (data!=NULL) delete[] data;}

		DType& operator[] (int i) {
      //assert( (i>=0) && (i<dim) ); // Somehow generate strict-overflow warning...
      assert( i>=0 );
      assert( i<dim );
			return data[i];
    }

		const DType& operator[] (int i) const {
      //assert( (i>=0) && (i<dim) );
      assert( i>=0 );
      assert( i<dim );
			return data[i];
    }

		HVEC& operator= (DType in_data) {
			for (int d=0; d<dim; d++) 
        data[d] = in_data; 
      return *this;
    }

		HVEC& operator= (const HVEC& in_vec) {
      assert( dim == in_vec.size() );
			for (int d=0; d<dim; d++) 
        data[d] = in_vec[d]; 
      return *this;
    }

		HVEC& operator+= (const HVEC& in_vec) {
      assert( dim == in_vec.size() );
			for (int d=0; d<dim; d++) 
        data[d] += in_vec[d]; 
      return *this;
    }

		HVEC& operator+= (DType in_data) {
			for (int d=0; d<dim; d++) 
        data[d] += in_data; 
      return *this;
    }

		HVEC& operator-= (const HVEC& in_vec) {
      assert( dim == in_vec.size() );
			for (int d=0; d<dim; d++) 
        data[d] -= in_vec[d]; 
      return *this;
    }

		HVEC& operator-= (DType in_data) {
			for (int d=0; d<dim; d++) 
        data[d] -= in_data; 
      return *this;
    }

		HVEC& operator*= (DType in_data) {
			for (int d=0; d<dim; d++) 
        data[d] *= in_data; 
      return *this;
    }

		HVEC& operator/= (DType in_data) {
			for (int d=0; d<dim; d++) 
        data[d] /= in_data; 
      return *this;
    }

		DType norm(void) const {	// Euclidean nrom
			DType entry, len;	len = 0.0;
			for (int d=0; d<dim; d++) {
				entry = data[d]; len += entry*entry;}
			return sqrt(len);}

		DType norm(int p) const {	// Lp norm, Linf if p<=0
			DType entry, len;	len = 0.0;
			for (int d=0; d<dim; d++) {
				entry = fabs(data[d]);
				if (p<=0) len = max(len,entry);
				else len += pow(entry,p);
			}
			if (p>0) len = pow(len,1.0/(1.0*p));
			return len;}
};

HVEC operator+ (const HVEC&, const HVEC&);
HVEC operator+ (const HVEC&);
HVEC operator+ (const HVEC&, const DType&);
HVEC operator+ (const DType&, const HVEC&);
HVEC operator- (const HVEC&, const HVEC&);
HVEC operator- (const HVEC&);
HVEC operator- (const HVEC&, const DType&);
HVEC operator- (const DType&, const HVEC&);
HVEC operator* (const DType&, const HVEC&);
HVEC operator* (const HVEC&, const DType&);
HVEC operator/ (const HVEC&, const DType&);
// This is the inner product
DType operator* (const HVEC&, const HVEC&);
// This is the exterior (cross) product (one must have dim=3)
HVEC operator^ (const HVEC&, const HVEC&);

class HMAT;

class TriDMAT {
    private:
        int     dim;
    public:
        DType*  upper;
        DType*  mdiag;
        DType*  lower;
        DType   dummy_zero;

        int size(void) const {return dim;}

        void print(void) {
            fprintf(stdout,"Upper diagonal:\t");
            for (int i=0; i<dim-1; i++)
                fprintf(stdout,"%10.20f\t",upper[i]);
            fprintf(stdout,"\n");
            fprintf(stdout,"Main diagonal :\t");
            for (int i=0; i<dim; i++)
                fprintf(stdout,"%10.20f\t",mdiag[i]);
            fprintf(stdout,"\n");
            fprintf(stdout,"Lower diagonal:\t");
            for (int i=0; i<dim-1; i++)
                fprintf(stdout,"%10.20f\t",lower[i]);
            fprintf(stdout,"\n");}

        TriDMAT(int in_dim) {
            dim = in_dim;
            dummy_zero = 0.0;
            upper = NULL;
            lower = NULL;
            mdiag = NULL;
            if (dim>1) upper = new DType[dim-1];
            if (dim>1) lower = new DType[dim-1];
            mdiag = new DType[dim];
            for (int i=0; i<dim-1; i++) upper[i] = 0.0;
            for (int i=0; i<dim-1; i++) lower[i] = 0.0;
            for (int i=0; i<dim; i++)   mdiag[i] = 0.0;}

        TriDMAT(int in_dim, DType in_val) {
            dim = in_dim;
            dummy_zero = 0.0;
            upper = NULL;
            lower = NULL;
            mdiag = NULL;
            if (dim>1) upper = new DType[dim-1];
            if (dim>1) lower = new DType[dim-1];
            mdiag = new DType[dim];
            for (int i=0; i<dim-1; i++) upper[i] = in_val;
            for (int i=0; i<dim-1; i++) lower[i] = in_val;
            for (int i=0; i<dim; i++)   mdiag[i] = in_val;}

        TriDMAT(const HVEC& in_vec) {
            dim = in_vec.size();
            dummy_zero = 0.0;
            upper = NULL;
            lower = NULL;
            mdiag = NULL;
            if (dim>1) upper = new DType[dim-1];
            if (dim>1) lower = new DType[dim-1];
            mdiag = new DType[dim];
            for (int i=0; i<dim-1; i++) upper[i] = 0.0;
            for (int i=0; i<dim-1; i++) lower[i] = 0.0;
            for (int i=0; i<dim; i++)   mdiag[i] = in_vec[i];}

        TriDMAT(const TriDMAT& in_mat) {
            dim = in_mat.size();
            dummy_zero = 0.0;
            upper = NULL;
            lower = NULL;
            mdiag = NULL;
            if (dim>1) upper = new DType[dim-1];
            if (dim>1) lower = new DType[dim-1];
            mdiag = new DType[dim];
            for (int i=0; i<dim-1; i++) 
                upper[i] = in_mat.upper[i];
            for (int i=0; i<dim-1; i++) 
                lower[i] = in_mat.lower[i];
            for (int i=0; i<dim; i++)   
                mdiag[i] = in_mat.mdiag[i];}

        ~TriDMAT(void) {
            if (upper!=NULL) delete[] upper;
            if (lower!=NULL) delete[] lower;
            if (mdiag!=NULL) delete[] mdiag;}

        // WARNING: this cannot be used to ASSIGN value
        // but no warning will be casted if doing so
        HVEC getCol(int i) {
			if ((i<0)||(i>=dim)) {
				fprintf(stderr,"Error in TriDMAT: index outside range %d!\n",dim);
				exit(1);
			}
            HVEC vec(dim,0.0);
            if (i>0) vec[i-1] = upper[i-1];
            vec[i] = mdiag[i];
            if (i<dim-1) vec[i+1] = lower[i];
            return vec;
        }

        // WARNING: this cannot be used to ASSIGN value
        // but no warning will be casted if doing so
        HVEC getRow(int i) {
			if ((i<0)||(i>=dim)) {
				fprintf(stderr,"Error in TriDMAT: index outside range %d!\n",dim);
				exit(1);
			}
            HVEC vec(dim,0.0);
            if (i>0) vec[i-1] = lower[i-1];
            vec[i] = mdiag[i];
            if (i<dim-1) vec[i+1] = upper[i];
            return vec;
        }

        DType& get(int i, int j) {
			if ((i<0)||(i>=dim)||(j<0)||(j>=dim)) {
				fprintf(stderr,"Error in TriDMAT: index outside range %d!\n",dim);
				exit(1);
			}
            if (i==j) return mdiag[i];
            else if (i==(j+1)) return upper[j];
            else if (j==(i+1)) return lower[i];
            else return dummy_zero;
        }
        const DType& get(int i, int j) const {
			if ((i<0)||(i>=dim)||(j<0)||(j>=dim)) {
				fprintf(stderr,"Error in TriDMAT: index outside range %d!\n",dim);
				exit(1);
			}
            if (i==j) return mdiag[i];
            else if (i==(j+1)) return upper[j];
            else if (j==(i+1)) return lower[i];
            else return dummy_zero;
        }

        // WARNING: this cannot be used to ASSIGN value
        // but no warning will be casted if doing so
        HVEC operator[] (int i) {
			if ((i<0)||(i>=dim)) {
				fprintf(stderr,"Error in TriDMAT: index outside range %d!\n",dim);
				exit(1);
			}
            HVEC vec(dim,0.0);
            if (i>0) vec[i-1] = lower[i-1];
            vec[i] = mdiag[i];
            if (i<dim-1) vec[i+1] = upper[i];
            return vec;
        }

        TriDMAT& operator= (DType in_data) {
            for (int i=0; i<dim-1; i++) upper[i] = in_data;
            for (int i=0; i<dim-1; i++) lower[i] = in_data;
            for (int i=0; i<dim; i++)   mdiag[i] = in_data;
            return *this;
        }

        TriDMAT& operator= (const TriDMAT& in_mat) {
            assert( dim == in_mat.size() );
            for (int i=0; i<dim-1; i++)
                upper[i] = in_mat.upper[i];
            for (int i=0; i<dim-1; i++)
                lower[i] = in_mat.lower[i];
            for (int i=0; i<dim; i++)
                mdiag[i] = in_mat.mdiag[i];
            return *this;
        }

        TriDMAT& operator= (const HVEC& in_vec) {
            assert( dim == in_vec.size() );
            for (int i=0; i<dim-1; i++)
                upper[i] = 0.0;
            for (int i=0; i<dim-1; i++)
                lower[i] = 0.0;
            for (int i=0; i<dim; i++)
                mdiag[i] = in_vec[i];
            return *this;
        }

        TriDMAT& operator+= (const TriDMAT& in_mat) {
            assert( dim == in_mat.size() );
            for (int i=0; i<dim-1; i++)
                upper[i] += in_mat.upper[i];
            for (int i=0; i<dim-1; i++)
                lower[i] += in_mat.lower[i];
            for (int i=0; i<dim; i++)
                mdiag[i] += in_mat.mdiag[i];
            return *this;
        }

        TriDMAT& operator+= (DType in_data) {
            for (int i=0; i<dim-1; i++)
                upper[i] += in_data;
            for (int i=0; i<dim-1; i++)
                lower[i] += in_data;
            for (int i=0; i<dim; i++)
                mdiag[i] += in_data;
            return *this;
        }

        TriDMAT& operator-= (const TriDMAT& in_mat) {
            assert( dim == in_mat.size() );
            for (int i=0; i<dim-1; i++)
                upper[i] -= in_mat.upper[i];
            for (int i=0; i<dim-1; i++)
                lower[i] -= in_mat.lower[i];
            for (int i=0; i<dim; i++)
                mdiag[i] -= in_mat.mdiag[i];
            return *this;
        }

        TriDMAT& operator-= (DType in_data) {
            for (int i=0; i<dim-1; i++)
                upper[i] -= in_data;
            for (int i=0; i<dim-1; i++)
                lower[i] -= in_data;
            for (int i=0; i<dim; i++)
                mdiag[i] -= in_data;
            return *this;
        }

        TriDMAT& operator*= (DType in_data) {
            for (int i=0; i<dim-1; i++)
                upper[i] *= in_data;
            for (int i=0; i<dim-1; i++)
                lower[i] *= in_data;
            for (int i=0; i<dim; i++)
                mdiag[i] *= in_data;
            return *this;
        }

        TriDMAT& operator/= (DType in_data) {
            for (int i=0; i<dim-1; i++)
                upper[i] /= in_data;
            for (int i=0; i<dim-1; i++)
                lower[i] /= in_data;
            for (int i=0; i<dim; i++)
                mdiag[i] /= in_data;
            return *this;
        }

        TriDMAT trans(void) {
            TriDMAT mat(dim);
            for (int i=0; i<dim; i++)
                mat.mdiag[i] = mdiag[i];
            for (int i=0; i<dim-1; i++) {
                mat.upper[i] = lower[i];
                mat.lower[i] = upper[i];
            }
            return mat;
        }

        DType det(void); 
        HMAT  inv(void); 

        DType norm(void) const { // Spectral norm
            return norm('m',2);
        }

        DType norm(char type, int p) const {
            if ((type=='f')||(type=='F'))
                return norm_frb(p);
            else if ((type=='m')||(type=='M'))
                return norm_mat(p);
            else {
                fprintf(stderr,"Warning: norm type not defined for TriDMAT, return Frob norm");
                return norm_frb(p);
            }
        }

        DType norm_frb(int p) const {
            DType entry, len;
            len = 0.0;
            for (int i=0; i<dim-1; i++) {
                entry = fabs(upper[i]);
                if (p<=0) len = max(len,entry);
                else len += pow(entry,p);
            }
            for (int i=0; i<dim-1; i++) {
                entry = fabs(lower[i]);
                if (p<=0) len = max(len,entry);
                else len += pow(entry,p);
            }
            for (int i=0; i<dim; i++) {
                entry = fabs(mdiag[i]);
                if (p<=0) len = max(len,entry);
                else len += pow(entry,p);
            }
            if (p>0) len = pow(len,1.0/(1.0*p));
            return fabs(len);
        }

        DType norm_mat(int p) const {
            fprintf(stderr,"Error: matrix norm not implemented for TriDMAT ... Exiting ...\n");
            exit(-1);
            return 0.0;
        }
};

TriDMAT operator+ (const TriDMAT&, const TriDMAT&);
TriDMAT operator+ (const TriDMAT&);
TriDMAT operator+ (const TriDMAT&, const DType&);
TriDMAT operator+ (const DType&, const TriDMAT&);
TriDMAT operator- (const TriDMAT&, const TriDMAT&);
TriDMAT operator- (const TriDMAT&);
TriDMAT operator- (const TriDMAT&, const DType&);
TriDMAT operator- (const DType&, const TriDMAT&);
TriDMAT operator* (const TriDMAT&, const DType&);
TriDMAT operator* (const DType&, const TriDMAT&);
TriDMAT operator/ (const TriDMAT&, const DType&);
HVEC    operator* (const TriDMAT&, const HVEC&);
TriDMAT operator* (const TriDMAT&, const TriDMAT&);

class LTriMAT;

class UTriMAT {
    private:
        int     dim;
        int     len;
    public:
        DType*  data; // Stored row by row
        DType   dummy_zero;

        int size(void) const {return dim;}

        void print(void);

        int index(int i, int j) const {
            assert(i<=j);
            return j-i+(i*(2*dim-i+1))/2;
        }

        UTriMAT(int in_dim) {
            dim = in_dim;
            len = dim*(1+dim)/2;
            dummy_zero = 0.0;
            data = NULL;
            data = new DType[len];
            for (int i=0; i<len; i++) data[i] = 0.0;
        }

        UTriMAT(int in_dim, DType in_val) {
            dim = in_dim;
            len = dim*(1+dim)/2;
            dummy_zero = 0.0;
            data = NULL;
            data = new DType[len];
            for (int i=0; i<len; i++) data[i] = in_val;
        }

        UTriMAT(const HVEC& in_vec) {
            dim = in_vec.size();
            len = dim*(1+dim)/2;
            dummy_zero = 0.0;
            data = NULL;
            data = new DType[len];
            for (int i=0; i<dim; i++) 
                data[index(i,i)] = in_vec[i];
        }

        UTriMAT(const UTriMAT& in_mat) {
            dim = in_mat.size();
            len = dim*(1+dim)/2;
            dummy_zero = 0.0;
            data = NULL;
            data = new DType[len];
            for (int i=0; i<len; i++) 
                data[i] = in_mat.data[i];
        }

        ~UTriMAT(void) {
            if (data!=NULL) delete[] data;
        }

        // WARNING: this cannot be used to ASSIGN value
        // but no warning will be casted if doing so
        HVEC getCol(int i) {
            assert( (i>=0) && (i<dim) );
            HVEC vec(dim,0.0);
            for (int j=0; j<=i; j++)
                vec[j] = data[index(j,i)];
            return vec;
        }

        // WARNING: this cannot be used to ASSIGN value
        // but no warning will be casted if doing so
        HVEC getRow(int i) {
            assert( (i>=0) && (i<dim) );
            HVEC vec(dim,0.0);
            for (int j=i; j<dim; j++)
                vec[j] = data[index(i,j)];
            return vec;
        }

        DType& get(int i, int j) {
            assert( (i>=0) && (i<dim) );
            assert( (j>=0) && (j<dim) );
            if (i<=j) return data[index(i,j)];
            else      return dummy_zero;
        }

        const DType& get(int i, int j) const {
            assert( (i>=0) && (i<dim) );
            assert( (j>=0) && (j<dim) );
            if (i<=j) return data[index(i,j)];
            else      return dummy_zero;
        }

        // WARNING: this cannot be used to ASSIGN value
        // but no warning will be casted if doing so
        HVEC operator[] (int i) {
            assert( (i>=0) && (i<dim) );
            HVEC vec(dim,0.0);
            for (int j=i; j<dim; j++)
                vec[j] = data[index(i,j)];
            return vec;
        }

        UTriMAT& operator= (DType in_data) {
            for (int i=0; i<len; i++) data[i] = in_data;
            return *this;
        }

        UTriMAT& operator= (const UTriMAT& in_mat) {
            assert( dim==in_mat.size() );
            for (int i=0; i<len; i++) data[i] = in_mat.data[i];
            return *this;
        }

        UTriMAT& operator= (const HVEC& in_vec) {
            assert( dim==in_vec.size() );
            for (int i=0; i<len; i++) data[i] = 0.0;
            for (int i=0; i<dim; i++)
                data[index(i,i)] = in_vec[i];
            return *this;
        }

        UTriMAT& operator+= (const UTriMAT& in_mat) {
            assert( dim==in_mat.size() );
            for (int i=0; i<len; i++)
                data[i] += in_mat.data[i];
            return *this;
        }

        UTriMAT& operator+= (DType in_data) {
            for (int i=0; i<len; i++) data[i] += in_data;
            return *this;
        }

        UTriMAT& operator-= (const UTriMAT& in_mat) {
            assert( dim==in_mat.size() );
            for (int i=0; i<len; i++)
                data[i] -= in_mat.data[i];
            return *this;
        }

        UTriMAT& operator-= (DType in_data) {
            for (int i=0; i<len; i++) data[i] -= in_data;
            return *this;
        }

        UTriMAT& operator*= (DType in_data) {
            for (int i=0; i<len; i++) data[i] *= in_data;
            return *this;
        }

        UTriMAT& operator/= (DType in_data) {
            for (int i=0; i<len; i++) data[i] /= in_data;
            return *this;
        }

        LTriMAT trans(void); 

        DType det(void) {
            DType val = 1.0;
            for (int i=0; i<dim; i++)
                val *= data[index(i,i)];
            return val;
        }

        UTriMAT inv(void);

        // Expand by one column at the end
        void expand(DType val=0.0);

        // Upper-left sub matrix
        UTriMAT sub(int);

        DType norm(void) const { // Spectral norm
            return norm('m',2);
        }

        DType norm(char type, int p) const {
            if ((type=='f')||(type=='F'))
                return norm_frb(p);
            else if ((type=='m')||(type=='M'))
                return norm_mat(p);
            else {
                fprintf(stderr,"Warning: norm type not defined for UTriMAT, return Frob norm");
                return norm_frb(p);
            }
        }

        DType norm_frb(int p) const {
            DType entry, l_nrm;
            l_nrm = 0.0;
            for (int i=0; i<len; i++) {
                entry = fabs(data[i]);
                if (p<=0) l_nrm = max(l_nrm,entry);
                else l_nrm += pow(entry,p);
            }
            if (p>0) l_nrm = pow(l_nrm,1.0/(1.0*p));
            return fabs(l_nrm);
        }

        DType norm_mat(int p) const {
            fprintf(stderr,"Error: matrix norm not implemented for UTriMAT ... Exiting ...\n");
            exit(-1);
            return 0.0;
        }
};

UTriMAT operator+ (const UTriMAT&, const UTriMAT&);
UTriMAT operator+ (const UTriMAT&);
UTriMAT operator+ (const UTriMAT&, const DType&);
UTriMAT operator+ (const DType&, const UTriMAT&);
UTriMAT operator- (const UTriMAT&, const UTriMAT&);
UTriMAT operator- (const UTriMAT&);
UTriMAT operator- (const UTriMAT&, const DType&);
UTriMAT operator- (const DType&, const UTriMAT&);
UTriMAT operator* (const UTriMAT&, const DType&);
UTriMAT operator* (const DType&, const UTriMAT&);
UTriMAT operator/ (const UTriMAT&, const DType&);
HVEC    operator* (const UTriMAT&, const HVEC&);
UTriMAT operator* (const UTriMAT&, const UTriMAT&);
UTriMAT operator^ (const UTriMAT&, const UTriMAT&);

class LTriMAT {
    private:
        int     dim;
        int     len;
    public:
        DType*  data; // Stored row by row
        DType   dummy_zero;

        int size(void) const {return dim;}

        void print(void);

        int index(int i, int j) const {
            assert(i>=j);
            return j+(i*(i+1))/2;
        }

        LTriMAT(int in_dim) {
            dim = in_dim;
            len = dim*(1+dim)/2;
            dummy_zero = 0.0;
            data = NULL;
            data = new DType[len];
            for (int i=0; i<len; i++) data[i] = 0.0;
        }

        LTriMAT(int in_dim, DType in_val) {
            dim = in_dim;
            len = dim*(1+dim)/2;
            dummy_zero = 0.0;
            data = NULL;
            data = new DType[len];
            for (int i=0; i<len; i++) data[i] = in_val;
        }

        LTriMAT(const HVEC& in_vec) {
            dim = in_vec.size();
            len = dim*(1+dim)/2;
            dummy_zero = 0.0;
            data = NULL;
            data = new DType[len];
            for (int i=0; i<dim; i++) 
                data[index(i,i)] = in_vec[i];
        }

        LTriMAT(const LTriMAT& in_mat) {
            dim = in_mat.size();
            len = dim*(1+dim)/2;
            dummy_zero = 0.0;
            data = NULL;
            data = new DType[len];
            for (int i=0; i<len; i++) 
                data[i] = in_mat.data[i];
        }

        ~LTriMAT(void) {
            if (data!=NULL) delete[] data;
        }

        // WARNING: this cannot be used to ASSIGN value
        // but no warning will be casted if doing so
        HVEC getCol(int i) {
            assert ( (i>=0) && (i<dim) );
            HVEC vec(dim,0.0);
            for (int j=0; j<=i; j++)
                vec[j] = data[index(j,i)];
            return vec;
        }

        // WARNING: this cannot be used to ASSIGN value
        // but no warning will be casted if doing so
        HVEC getRow(int i) {
            assert ( (i>=0) && (i<dim) );
            HVEC vec(dim,0.0);
            for (int j=i; j<dim; j++)
                vec[j] = data[index(i,j)];
            return vec;
        }

        DType& get(int i, int j) {
            assert ( (i>=0) && (i<dim) );
            assert ( (j>=0) && (j<dim) );
            if (i>=j) return data[index(i,j)];
            else      return dummy_zero;
        }

        const DType& get(int i, int j) const {
            assert ( (i>=0) && (i<dim) );
            assert ( (j>=0) && (j<dim) );
            if (i>=j) return data[index(i,j)];
            else      return dummy_zero;
        }

        // WARNING: this cannot be used to ASSIGN value
        // but no warning will be casted if doing so
        HVEC operator[] (int i) {
            assert ( (i>=0) && (i<dim) );
            HVEC vec(dim,0.0);
            for (int j=i; j<dim; j++)
                vec[j] = data[index(i,j)];
            return vec;
        }

        LTriMAT& operator= (DType in_data) {
            for (int i=0; i<len; i++) data[i] = in_data;
            return *this;
        }

        LTriMAT& operator= (const LTriMAT& in_mat) {
            assert( dim == in_mat.size() );
            for (int i=0; i<len; i++) data[i] = in_mat.data[i];
            return *this;
        }

        LTriMAT& operator= (const HVEC& in_vec) {
            assert( dim == in_vec.size() );
            for (int i=0; i<len; i++) data[i] = 0.0;
            for (int i=0; i<dim; i++)
                data[index(i,i)] = in_vec[i];
            return *this;
        }

        LTriMAT& operator+= (const LTriMAT& in_mat) {
            assert( dim == in_mat.size() );
            for (int i=0; i<len; i++)
                data[i] += in_mat.data[i];
            return *this;
        }

        LTriMAT& operator+= (DType in_data) {
            for (int i=0; i<len; i++) data[i] += in_data;
            return *this;
        }

        LTriMAT& operator-= (const LTriMAT& in_mat) {
            assert( dim == in_mat.size() );
            for (int i=0; i<len; i++)
                data[i] -= in_mat.data[i];
            return *this;
        }

        LTriMAT& operator-= (DType in_data) {
            for (int i=0; i<len; i++) data[i] -= in_data;
            return *this;
        }

        LTriMAT& operator*= (DType in_data) {
            for (int i=0; i<len; i++) data[i] *= in_data;
            return *this;
        }

        LTriMAT& operator/= (DType in_data) {
            for (int i=0; i<len; i++) data[i] /= in_data;
            return *this;
        }

        UTriMAT trans(void); 

        DType det(void) {
            DType val = 1.0;
            for (int i=0; i<dim; i++)
                val *= data[index(i,i)];
            return val;
        }

        LTriMAT inv(void);

        // Expand by one column at the end
        void expand(DType val=0.0);

        // Upper-left sub matrix
        LTriMAT sub(int);

        DType norm(void) const { // Spectral norm
            return norm('m',2);
        }

        DType norm(char type, int p) const {
            if ((type=='f')||(type=='F'))
                return norm_frb(p);
            else if ((type=='m')||(type=='M'))
                return norm_mat(p);
            else {
                fprintf(stderr,"Warning: norm type not defined for LTriMAT, return Frob norm");
                return norm_frb(p);
            }
        }

        DType norm_frb(int p) const {
            DType entry, l_nrm;
            l_nrm = 0.0;
            for (int i=0; i<len; i++) {
                entry = fabs(data[i]);
                if (p<=0) l_nrm = max(l_nrm,entry);
                else l_nrm += pow(entry,p);
            }
            if (p>0) l_nrm = pow(l_nrm,1.0/(1.0*p));
            return fabs(l_nrm);
        }

        DType norm_mat(int p) const {
            fprintf(stderr,"Error: matrix norm not implemented for LTriMAT ... Exiting ...\n");
            exit(-1);
            return 0.0;
        }
};

LTriMAT operator+ (const LTriMAT&, const LTriMAT&);
LTriMAT operator+ (const LTriMAT&);
LTriMAT operator+ (const LTriMAT&, const DType&);
LTriMAT operator+ (const DType&, const LTriMAT&);
LTriMAT operator- (const LTriMAT&, const LTriMAT&);
LTriMAT operator- (const LTriMAT&);
LTriMAT operator- (const LTriMAT&, const DType&);
LTriMAT operator- (const DType&, const LTriMAT&);
LTriMAT operator* (const LTriMAT&, const DType&);
LTriMAT operator* (const DType&, const LTriMAT&);
LTriMAT operator/ (const LTriMAT&, const DType&);
HVEC    operator* (const LTriMAT&, const HVEC&);
LTriMAT operator* (const LTriMAT&, const LTriMAT&);
LTriMAT operator^ (const LTriMAT&, const LTriMAT&);

class HTENS2;

class HMAT {
	private:
		int		dim;
		HVEC	**data;

		// Find the index of the entry with largest absolute value
		// in the sub-matrix [I,dim]X[I,dim]
		void findMatMax(int*,int*,HMAT&,int);
		void findRowMax(int*,HVEC&,int);
		// Permute two rows or two columns
		void permRow(HMAT&,int,int);
		void permCol(HMAT&,int,int);
	public:
		int size(void) const {return dim;}

		void print(void) {
			for (int i=0; i<dim; i++) {
				for (int j=0; j<dim; j++)
					fprintf(stdout,"%10.20f\t",(*data[i])[j]);
				fprintf(stdout,"\n");
			}
			fprintf(stdout,"\n");
			return;
		}
		
		HMAT(int in_dim, DType in_val = 0.0) {
			dim = in_dim;
			data = new HVEC*[dim];
			for (int d=0; d<dim; d++)
				data[d] = new HVEC(dim);
			for (int i=0; i<dim; i++) 
				for (int j=0; j<dim; j++) 
					(*data[i])[j] = in_val;}

		HMAT(const HVEC& in_vec) {
			dim = in_vec.size();
			data = new HVEC*[dim];
			for (int d=0; d<dim; d++)
				data[d] = new HVEC(dim,0.0);
			for (int i=0; i<dim; i++)
				(*data[i])[i] = in_vec[i];}

		HMAT(const HMAT& in_mat) {
			dim = in_mat.size();
			data = new HVEC*[dim];
			for (int d=0; d<dim; d++)
				data[d] = new HVEC(dim);
			for (int i=0; i<dim; i++)
				for (int j=0; j<dim; j++)
					(*data[i])[j] = in_mat[i][j];}

        HMAT(const TriDMAT& in_mat) {
            dim = in_mat.size();
			data = new HVEC*[dim];
			for (int d=0; d<dim; d++)
				data[d] = new HVEC(dim,0.0);
            for (int d=0; d<dim-1; d++)
                (*data[d])[d+1] = in_mat.upper[d];
            for (int d=0; d<dim; d++)
                (*data[d])[d] = in_mat.mdiag[d];
            for (int d=0; d<dim-1; d++)
                (*data[d+1])[d] = in_mat.lower[d];
        }

        HMAT(const UTriMAT& in_mat) {
            dim = in_mat.size();
            data = new HVEC*[dim];
            for (int d=0; d<dim; d++) {
                data[d] = new HVEC(dim,0.0);
                for (int i=d; i<dim; i++)
                    (*data[d])[i] = in_mat.data[in_mat.index(d,i)];
            }
        }

        HMAT(const LTriMAT& in_mat) {
            dim = in_mat.size();
            data = new HVEC*[dim];
            for (int d=0; d<dim; d++) {
                data[d] = new HVEC(dim,0.0);
                for (int i=0; i<=d; i++)
                    (*data[d])[i] = in_mat.data[in_mat.index(d,i)];
            }
        }

        HMAT(const HTENS2&);

		~HMAT(void) {
			if (data!=NULL) {
			    for (int d=0; d<dim; d++)
			    	if (data[d]!=NULL) delete data[d];
			    delete[] data;
			}
		}

        // WARNING: this cannot be used to ASSIGN value
        // but no warning will be casted if doing so
        HVEC getCol(int i) {
            assert( (i>=0) && (i<dim) );
            HVEC vec(dim,0.0);
            for (int d=0; d<dim; d++)
                vec[d] = (*data[d])[i];
            return vec;
        }

        // WARNING: this cannot be used to ASSIGN value
        // but no warning will be casted if doing so
        HVEC getRow(int i) {
            assert( (i>=0) && (i<dim) );
            HVEC vec(dim,0.0);
            vec = *data[i];
            return vec;
        }

        DType& get(int i, int j) {
            assert( (i>=0) && (i<dim) );
            assert( (j>=0) && (j<dim) );
            return (*data[i])[j];
        }

        const DType& get(int i, int j) const {
            assert( (i>=0) && (i<dim) );
            assert( (j>=0) && (j<dim) );
            return (*data[i])[j];
        }

		HVEC& operator[] (int i) {
            assert( (i>=0) && (i<dim) );
			return (*data[i]);}

		const HVEC& operator[] (int i) const {
            assert( (i>=0) && (i<dim) );
			return (*data[i]);}

		HMAT& operator= (DType in_data) {
			for (int i=0; i<dim; i++)
				for (int j=0; j<dim; j++)
					(*data[i])[j] = in_data;
			return *this;}

		HMAT& operator= (const HMAT& in_mat) {
            assert( dim == in_mat.size() );
			for (int i=0; i<dim; i++)
				for (int j=0; j<dim; j++)
					(*data[i])[j] = in_mat[i][j];
			return *this;}

		HMAT& operator= (const HVEC& in_vec) {
            assert( dim == in_vec.size() );
			for (int i=0; i<dim; i++) {
				(*data[i]) = 0.0;	(*data[i])[i] = in_vec[i];}
			return *this;
		}

        HMAT& operator= (const TriDMAT& in_mat) {
            assert( dim == in_mat.size() );
            for (int i=0; i<dim; i++) {
                *data[i] = 0.0;
                (*data[i])[i] = in_mat.mdiag[i];
            }
            for (int i=0; i<dim-1; i++) {
                (*data[i])[i+1] = in_mat.upper[i];
                (*data[i+1])[i] = in_mat.lower[i];
            }
            return *this;
        }

        HMAT& operator= (const UTriMAT& in_mat) {
            assert( dim == in_mat.size() );
            for (int i=0; i<dim; i++) *data[i] = 0.0;
            for (int i=0; i<dim; i++)
                for (int j=i; j<dim; j++)
                    (*data[i])[j] = in_mat.data[in_mat.index(i,j)];
            return *this;
        }

        HMAT& operator= (const LTriMAT& in_mat) {
            assert( dim == in_mat.size() );
            for (int i=0; i<dim; i++) *data[i] = 0.0;
            for (int i=0; i<dim; i++)
                for (int j=0; j<=i; j++)
                    (*data[i])[j] = in_mat.data[in_mat.index(i,j)];
            return *this;
        }

        HMAT& operator= (const HTENS2&);

		HMAT& operator+=(const HMAT& in_mat) {
            assert( dim == in_mat.size() );
			for (int i=0; i<dim; i++)	(*data[i]) += in_mat[i];
			return *this;
		}

		HMAT& operator+=(DType in_data) {
			for (int i=0; i<dim; i++)	(*data[i]) += in_data;
			return *this;
		}

        HMAT& operator+=(const TriDMAT& in_mat) {
            assert( dim == in_mat.size() );
            for (int i=0; i<dim; i++) {
                (*data[i])[i] += in_mat.mdiag[i];
            }
            for (int i=0; i<dim-1; i++) {
                (*data[i])[i+1] += in_mat.upper[i];
                (*data[i+1])[i] += in_mat.lower[i];
            }
            return *this;
        }

        HMAT& operator+=(const UTriMAT& in_mat) {
            assert( dim == in_mat.size() );
            for (int i=0; i<dim; i++)
                for (int j=i; j<dim; j++)
                    (*data[i])[j] += in_mat.data[in_mat.index(i,j)];
            return *this;
        }

        HMAT& operator+=(const LTriMAT& in_mat) {
            assert( dim == in_mat.size() );
            for (int i=0; i<dim; i++)
                for (int j=0; j<=i; j++)
                    (*data[i])[j] += in_mat.data[in_mat.index(i,j)];
            return *this;
        }

		HMAT& operator-=(const HMAT& in_mat) {
            assert( dim == in_mat.size() );
			for (int i=0; i<dim; i++)	(*data[i]) -= in_mat[i];
			return *this;
		}

		HMAT& operator-=(DType in_data) {
			for (int i=0; i<dim; i++)	(*data[i]) -= in_data;
			return *this;
		}

        HMAT& operator-=(const TriDMAT& in_mat) {
            assert( dim == in_mat.size() );
            for (int i=0; i<dim; i++) {
                (*data[i])[i] -= in_mat.mdiag[i];
            }
            for (int i=0; i<dim-1; i++) {
                (*data[i])[i+1] -= in_mat.upper[i];
                (*data[i+1])[i] -= in_mat.lower[i];
            }
            return *this;
        }

        HMAT& operator-=(const UTriMAT& in_mat) {
            assert( dim == in_mat.size() );
            for (int i=0; i<dim; i++)
                for (int j=i; j<dim; j++)
                    (*data[i])[j] -= in_mat.data[in_mat.index(i,j)];
            return *this;
        }

        HMAT& operator-=(const LTriMAT& in_mat) {
            assert( dim == in_mat.size() );
            for (int i=0; i<dim; i++)
                for (int j=0; j<=i; j++)
                    (*data[i])[j] -= in_mat.data[in_mat.index(i,j)];
            return *this;
        }

		HMAT& operator*=(DType in_data) {
			for (int i=0; i<dim; i++)	(*data[i]) *= in_data;
			return *this;
		}

		HMAT& operator/=(DType in_data) {
			for (int i=0; i<dim; i++)	(*data[i]) /= in_data;
			return *this;
		}

		HMAT trans(void) {
			HMAT mat(dim);
			for (int i=0; i<dim; i++)
				for (int j=0; j<dim; j++)
					mat[i][j] = (*data[j])[i];;
			return mat;}

		DType det_stable(void);	// PLU with largest entry pivoting
		DType det_fast(void);	// PLU with non-zero pivoting
		DType det(void) {return det_stable();} 

		HMAT inv_stable(void);	// PLU with largest entry pivoting
		HMAT inv_fast(void);	// PLU with non-zero pivoting
		HMAT inv(void) {return inv_stable();} 

        DType norm(void) const { // Spectral norm
            return norm('m',2);
        }

        DType norm(char type, int p) const {
            if ((type=='f')||(type=='F'))
                return norm_frb(p);
            else if ((type=='m')||(type=='M'))
                return norm_mat(p);
            else {
                fprintf(stderr,"Warning: norm type not defined for HMAT, return Frob norm");
                return norm_frb(p);
            }
        }

        DType norm_frb(int p) const {
            DType entry, len;
            len = 0.0;
            for (int i=0; i<dim; i++)
                for (int j=0; j<dim; j++) {
                    entry = fabs((*data[i])[j]);
                    if (p<=0) len = max(len,entry);
                    else len += pow(entry,p);
                }
            if (p>0) len = pow(len,1.0/(1.0*p));
            return fabs(len);
        }

        DType norm_mat(int p) const {
            fprintf(stderr,"Error: matrix norm not implemented for HMAT ... Exiting ...\n");
            exit(-1);
            return 0.0;
        }
};

HMAT operator+ (const HMAT&, const HMAT&);
HMAT operator+ (const HMAT&);
HMAT operator+ (const HMAT&, const DType&);
HMAT operator+ (const DType&, const HMAT&);
HMAT operator- (const HMAT&, const HMAT&);
HMAT operator- (const HMAT&);
HMAT operator- (const HMAT&, const DType&);
HMAT operator- (const DType&, const HMAT&);
HMAT operator* (const DType&, const HMAT&);
HMAT operator* (const HMAT&, const DType&);
HMAT operator/ (const HMAT&, const DType&);
HVEC operator* (const HMAT&, const HVEC&);
HMAT operator* (const HMAT&, const HMAT&);
// // This is the outer product
// // Should first compute HTENS2 and then convert
// HMAT operator| (const HVEC&, const HVEC&);
// This is the Kronecker product	
HMAT operator^ (const HMAT&, const HMAT&);

HMAT operator+ (const HMAT&, const TriDMAT&);
HMAT operator+ (const TriDMAT&, const HMAT&);
HMAT operator- (const HMAT&, const TriDMAT&);
HMAT operator- (const TriDMAT&, const HMAT&);
HMAT operator* (const HMAT&, const TriDMAT&);
HMAT operator* (const TriDMAT&, const HMAT&);
HMAT operator^ (const HMAT&, const TriDMAT&);
HMAT operator^ (const TriDMAT&, const HMAT&);
HMAT operator^ (const TriDMAT&, const TriDMAT&);

// For the operations involving UTriMAT and LTriMAT and returning HMAT
// The automatic conversion should be activated, which may
// take twice the optimal computation cost
// If not, these needs to be implemented

// Second-order tensor
class HTENS2 {
    private:
        int     dim1, dim2;
        HVEC    **data;
    public:
        int size(int i) const {
            if (i==1) 
                return dim1;
            if (i==2)
                return dim2;
            return -1;
        }

        void print(void) {
            for (int i=0; i<dim1; i++) {
                for (int j=0; j<dim2; j++)
                    fprintf(stdout,"%10.20f\t",(*data[i])[j]);
                fprintf(stdout,"\n");
            }
            fprintf(stdout,"\n");
            return;
        }

        HTENS2(int in_dim1, int in_dim2, DType in_val=0.0) {
            dim1 = in_dim1;
            dim2 = in_dim2;
            data = new HVEC*[dim1];
            for (int d=0; d<dim1; d++)
                data[d] = new HVEC(dim2);
            for (int i=0; i<dim1; i++)
                for (int j=0; j<dim2; j++)
                    (*data[i])[j] = in_val;
        }

        HTENS2(const HMAT& in_mat) {
            dim1 = in_mat.size();
            dim2 = in_mat.size();
            for (int d=0; d<dim1; d++)
                data[d] = new HVEC(dim2);
            for (int i=0; i<dim1; i++)
                for (int j=0; j<dim2; j++)
                    (*data[i])[j] = in_mat[i][j];
        }

        HTENS2(const HTENS2& in_tens) {
            dim1 = in_tens.size(1);
            dim2 = in_tens.size(2);
            for (int d=0; d<dim1; d++)
                data[d] = new HVEC(dim2);
            for (int i=0; i<dim1; i++)
                for (int j=0; j<dim2; j++)
                    (*data[i])[j] = in_tens[i][j];
        }

        ~HTENS2(void) {
            if (data) {
                for (int d=0; d<dim1; d++)
                    if (data[d]) delete data[d];
                delete[] data;
            }
        }

        DType& get(int i, int j) {
            assert( (i>=0) && (i<dim1) );
            assert( (j>=0) && (j<dim2) );
            return (*data[i])[j];
        }

        const DType& get(int i, int j) const {
            assert( (i>=0) && (i<dim1) );
            assert( (j>=0) && (j<dim2) );
            return (*data[i])[j];
        }

        HVEC& operator[] (int i) {
            assert( (i>=0) && (i<dim1) );
            return (*data[i]);
        }

        const HVEC& operator[] (int i) const {
            assert( (i>=0) && (i<dim1) );
            return (*data[i]);
        }

        HTENS2& operator= (DType in_data) {
            for (int i=0; i<dim1; i++)
                for (int j=0; j<dim2; j++)
                    (*data[i])[j] = in_data;
            return *this;
        }

        HTENS2& operator= (const HMAT& in_mat) {
            assert( dim1 == in_mat.size() );
            assert( dim2 == in_mat.size() );
            for (int i=0; i<dim1; i++)
                for (int j=0; j<dim2; j++)
                    (*data[i])[j] = in_mat[i][j];
            return *this;
        }

        HTENS2& operator= (const HTENS2& in_tens) {
            assert( dim1 == in_tens.size(1) );
            assert( dim2 == in_tens.size(2) );
            for (int i=0; i<dim1; i++)
                for (int j=0; j<dim2; j++)
                    (*data[i])[j] = in_tens[i][j];
            return *this;
        }

        HTENS2& operator+= (DType in_data) {
            for (int i=0; i<dim1; i++)
                (*data[i]) += in_data;
            return *this;
        }

        HTENS2& operator+= (const HMAT& in_mat) {
            assert( dim1 == in_mat.size() );
            assert( dim2 == in_mat.size() );
            for (int i=0; i<dim1; i++)
                (*data[i]) += in_mat[i];
            return *this;
        }

        HTENS2& operator+= (const HTENS2& in_tens) {
            assert( dim1 == in_tens.size(1) );
            assert( dim2 == in_tens.size(2) );
            for (int i=0; i<dim1; i++)
                (*data[i]) += in_tens[i];
            return *this;
        }

        HTENS2& operator-= (DType in_data) {
            for (int i=0; i<dim1; i++)
                (*data[i]) -= in_data;
            return *this;
        }

        HTENS2& operator-= (const HMAT& in_mat) {
            assert( dim1 == in_mat.size() );
            assert( dim2 == in_mat.size() );
            for (int i=0; i<dim1; i++)
                (*data[i]) -= in_mat[i];
            return *this;
        }

        HTENS2& operator-= (const HTENS2& in_tens) {
            assert( dim1 == in_tens.size(1) );
            assert( dim2 == in_tens.size(2) );
            for (int i=0; i<dim1; i++)
                (*data[i]) -= in_tens[i];
            return *this;
        }

        HTENS2& operator*= (DType in_data) {
            for (int i=0; i<dim1; i++)
                (*data[i]) *= in_data;
            return *this;
        }

        HTENS2& operator/= (DType in_data) {
            for (int i=0; i<dim1; i++)
                (*data[i]) /= in_data;
            return *this;
        }

};

HTENS2 operator+ (const HTENS2&, const HTENS2&);
HTENS2 operator+ (const HTENS2&, const DType&);
HTENS2 operator+ (const DType&, const HTENS2&);
HTENS2 operator+ (const HTENS2&);
HTENS2 operator- (const HTENS2&, const HTENS2&);
HTENS2 operator- (const HTENS2&, const DType&);
HTENS2 operator- (const DType&, const HTENS2&);
HTENS2 operator- (const HTENS2&);
HTENS2 operator* (const HTENS2&, const DType&);
HTENS2 operator* (const DType&, const HTENS2&);
HTENS2 operator/ (const HTENS2&, const DType&);
// This is usual product, not double contraction
HTENS2 operator* (const HTENS2&, const HTENS2&);
HTENS2 operator* (const HTENS2&, const HMAT&);
HTENS2 operator* (const HMAT&, const HTENS2&);
HVEC   operator* (const HTENS2&, const HVEC&);
// This is the outer product
HTENS2 operator| (const HVEC&, const HVEC&);
// This is the Kronecker product
HTENS2 operator^ (const HTENS2&, const HTENS2&);

// Third-order tensor
class HTENS3 {
    private:
        int      dim1, dim2, dim3;
        HTENS2** data;
    public:
        int size(int i) const {
            if (i==1)
                return dim1;
            if (i==2)
                return dim2;
            if (i==3)
                return dim3;
            return -1;
        }

        void print(void) {
            for (int i=0; i<dim1; i++) {
                (*data[i]).print();
                fprintf(stdout,"\n");
            }
            return;
        }

        HTENS3(int in_dim1, int in_dim2, int in_dim3, DType in_val=0.0) {
            dim1 = in_dim1;
            dim2 = in_dim2;
            dim3 = in_dim3;
            data = new HTENS2*[dim1];
            for (int d=0; d<dim1; d++)
                data[d] = new HTENS2(dim2,dim3,in_val);
        }

        HTENS3(const HTENS3& in_tens) {
            dim1 = in_tens.size(1);
            dim2 = in_tens.size(2);
            dim3 = in_tens.size(3);
            for (int d=0; d<dim1; d++)
                data[d] = new HTENS2(in_tens[d]);
        }

        ~HTENS3(void) {
            if (data) {
                for (int d=0; d<dim1; d++)
                    if (data[d]) delete data[d];
                delete[] data;
            }
        }

        DType& get(int i, int j, int k) {
            assert( (i>=0) && (i<dim1) );
            assert( (j>=0) && (j<dim2) );
            assert( (k>=0) && (k<dim3) );
            return (*data[i])[j][k];
        }

        const DType& get(int i, int j, int k) const {
            assert( (i>=0) && (i<dim1) );
            assert( (j>=0) && (j<dim2) );
            assert( (k>=0) && (k<dim3) );
            return (*data[i])[j][k];
        }

        HTENS2& operator[] (int i) {
            assert( (i>=0) && (i<dim1) );
            return (*data[i]);
        }

        const HTENS2& operator[] (int i) const {
            assert( (i>=0) && (i<dim1) );
            return (*data[i]);
        }

        HTENS3& operator= (DType in_data) {
            for (int i=0; i<dim1; i++)
                (*data[i]) = in_data;
            return *this;
        }

        HTENS3& operator= (const HTENS3& in_tens) {
            assert( dim1 = in_tens.size(1) );
            assert( dim2 = in_tens.size(2) );
            assert( dim3 = in_tens.size(3) );
            for (int i=0; i<dim1; i++)
                (*data[i]) = in_tens[i];
            return *this;
        }

        HTENS3& operator+= (DType in_data) {
            for (int i=0; i<dim1; i++)
                (*data[i]) += in_data;
            return *this;
        }

        HTENS3& operator+= (const HTENS3& in_tens) {
            assert( dim1 = in_tens.size(1) );
            assert( dim2 = in_tens.size(2) );
            assert( dim3 = in_tens.size(3) );
            for (int i=0; i<dim1; i++)
                (*data[i]) += in_tens[i];
            return *this;
        }

        HTENS3& operator-= (DType in_data) {
            for (int i=0; i<dim1; i++)
                (*data[i]) -= in_data;
            return *this;
        }

        HTENS3& operator-= (const HTENS3& in_tens) {
            assert( dim1 = in_tens.size(1) );
            assert( dim2 = in_tens.size(2) );
            assert( dim3 = in_tens.size(3) );
            for (int i=0; i<dim1; i++)
                (*data[i]) = -in_tens[i];
            return *this;
        }

        HTENS3& operator*= (DType in_data) {
            for (int i=0; i<dim1; i++)
                (*data[i]) *= in_data;
            return *this;
        }

        HTENS3& operator/= (DType in_data) {
            for (int i=0; i<dim1; i++)
                (*data[i]) /= in_data;
            return *this;
        }

};

HTENS3 operator+ (const HTENS3&, const HTENS3&);
HTENS3 operator+ (const HTENS3&, const DType&);
HTENS3 operator+ (const DType&, const HTENS3&);
HTENS3 operator+ (const HTENS3&);
HTENS3 operator- (const HTENS3&, const HTENS3&);
HTENS3 operator- (const HTENS3&, const DType&);
HTENS3 operator- (const DType&, const HTENS3&);
HTENS3 operator- (const HTENS3&);
HTENS3 operator* (const HTENS3&, const DType&);
HTENS3 operator* (const DType&, const HTENS3&);
HTENS3 operator/ (const HTENS3&, const DType&);
// (dim1 x dim2 x dim3) X (dim2 x dim3)
HVEC   operator* (const HTENS3&, const HTENS2&);
// (dim1 x dim2 x dim3) X (dim3)
HTENS2 operator* (const HTENS3&, const HVEC&);

// General tensors -- to think carefully how to implement

template <int ord> class HTENS {
    private:
        int dim;
        HTENS<ord-1>* sub_tens;
        DType* data;
    public:
        int getDim(int ord_in) {
            assert( (ord_in>=1) && (ord_in<=ord) );
            if (ord_in==1) 
                return dim;
            else
                return sub_tens[0].getDim(ord_in-1);
        }

        HTENS<ord-1>& operator[] (int i) {
            assert( (i>=0) && (i<dim) );
            return sub_tens[i];
        }
};

#endif


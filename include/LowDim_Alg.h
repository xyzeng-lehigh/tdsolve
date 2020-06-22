#ifndef LOWDIM_ALG_H
#define LOWDIM_ALG_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <cassert>
#include "Data_Type.h"
using namespace std;

/* 3D vectors and matrices
   VEC3
   MAT3
   */
class VEC3 {
	private:
		DType data[3];
	public:
		void print(void) {
			printf("%10.20f\t%10.20f\t%10.20f\n",data[0],data[1],data[2]); 
			return;
		}

		VEC3() {
			data[0]=0.0; 
			data[1]=0.0; 
			data[2]=0.0;
		}
		VEC3(DType in_val) {
			data[0]=in_val;
			data[1]=in_val;
			data[2]=in_val;
		}
		VEC3(DType in_1, DType in_2, DType in_3) {
			data[0]=in_1; 
			data[1]=in_2; 
			data[2]=in_3;
		}
		VEC3(const VEC3& in_vec) {
			data[0] = in_vec[0];
			data[1] = in_vec[1];
			data[2] = in_vec[2];
		}

		DType& operator[] (int i) {
            assert( (i>=0) && (i<=2) );
			return data[i];
		}

		const DType& operator[] (int i) const {
            assert( (i>=0) && (i<=2) );
			return data[i];
		}

		void operator= (DType in_data) {
			data[0]=in_data; 
			data[1]=in_data; 
			data[2]=in_data; 
			return;
		}
		void operator= (const VEC3& in_data) {
			data[0]=in_data[0];
			data[1]=in_data[1];
			data[2]=in_data[2];
			return;
		}
		void operator+= (const VEC3& in_data) {
			data[0]+=in_data[0];
			data[1]+=in_data[1];
			data[2]+=in_data[2];
			return;
		}
		void operator+= (DType in_data) {
			data[0]+=in_data;
			data[1]+=in_data;
			data[2]+=in_data;
			return;
		}
		void operator-= (const VEC3& in_data) {
			data[0]-=in_data[0];
			data[1]-=in_data[1];
			data[2]-=in_data[2];
			return;
		}
		void operator-= (DType in_data) {
			data[0]-=in_data;
			data[1]-=in_data;
			data[2]-=in_data;
			return;
		}
		void operator*= (DType in_data) {
			data[0]*=in_data;
			data[1]*=in_data;
			data[2]*=in_data;
			return;
		}
		void operator/= (DType in_data) {
			data[0]/=in_data;
			data[1]/=in_data;
			data[2]/=in_data;
			return;
		}

		DType length(void) const {
			return sqrt(data[0]*data[0]+data[1]*data[1]+data[2]*data[2]);
		}
};

VEC3 operator+ (const VEC3& input_1, const VEC3& input_2);
VEC3 operator+ (const VEC3& input_1);
VEC3 operator+ (const VEC3& input_v, const DType& input_s);
VEC3 operator+ (const DType& input_s, const VEC3& input_v);
VEC3 operator- (const VEC3& input_1, const VEC3& input_2);
VEC3 operator- (const VEC3& input_1);
VEC3 operator- (const VEC3& input_v, const DType& input_s);
VEC3 operator- (const DType& input_s, const VEC3& input_v);
VEC3 operator* (const DType& input_s, const VEC3& input_v);
VEC3 operator* (const VEC3& input_v, const DType& input_s);
VEC3 operator/ (const VEC3& input_v, const DType& input_s);

class MAT3 {
	private:
		// data[i] is the i-th row of the matrix
		VEC3 data[3];
	public:
		void print(void) {
			int i;
			for (i=0; i<3; i++)
				printf("%10.20f\t%10.20f\t%10.20f\n",data[i][0],data[i][1],data[i][2]); 
			printf("\n");
			return;
		}

		MAT3() {
			int i,j;
			for (i=0; i<3; i++)
				for (j=0; j<3; j++)
					data[i][j]=0.0;
		}
		MAT3(DType in_val) {
			int i,j;
			for (i=0; i<3; i++)
				for (j=0; j<3; j++)
					data[i][j]=in_val;
		}
		MAT3(DType in_1, DType in_2, DType in_3) {
			MAT3(0.0);
			data[0][0]=in_1; 
			data[1][1]=in_2; 
			data[2][2]=in_3;
		}
		MAT3(const VEC3& in_data) {
			MAT3(0.0);
			int i;
			for (i=0; i<3; i++)
				data[i][i]=in_data[i];
		}
		MAT3 (const MAT3& in_mat) {
			int i,j;
			for (i=0; i<3; i++)
				for (j=0; j<3; j++)
					data[i][j] = in_mat[i][j];
		}

		VEC3& operator[] (int i) {
            assert( (i>=0) && (i<=2) );
			return data[i];
		}

		const VEC3& operator[] (int i) const {
            assert( (i>=0) && (i<=2) );
			return data[i];
		}

		void operator= (DType in_data) {
			int i,j;
			for (i=0;i<3;i++)
				for (j=0;j<3;j++)
					data[i][j] = in_data;
			return;
		}
		void operator= (const MAT3& in_data) {
			int i,j;
			for (i=0;i<3;i++)
				for (j=0;j<3;j++)
					data[i][j]=in_data[i][j];
			return;
		}
		void operator+= (const MAT3& in_data) {
			int i,j;
			for (i=0;i<3;i++)
				for (j=0;j<3;j++)
					data[i][j]+=in_data[i][j];
			return;
		}
		void operator+= (DType in_data) {
			int i,j;
			for (i=0;i<3;i++)
				for (j=0;j<3;j++)
					data[i][j]+=in_data;
			return;
		}
		void operator-= (const MAT3& in_data) {
			int i,j;
			for (i=0;i<3;i++)
				for (j=0;j<3;j++)
					data[i][j]-=in_data[i][j];
			return;
		}
		void operator-= (DType in_data) {
			int i,j;
			for (i=0;i<3;i++)
				for (j=0;j<3;j++)
					data[i][j]-=in_data;
			return;
		}
		void operator*= (DType in_data) {
			int i,j;
			for (i=0;i<3;i++)
				for (j=0;j<3;j++)
					data[i][j]*=in_data;
			return;
		}
		void operator/= (DType in_data) {
			int i,j;
			for (i=0;i<3;i++)
				for (j=0;j<3;j++)
					data[i][j]/=in_data;
			return;
		}

		DType det(void) {
			return data[0][0]*(data[1][1]*data[2][2]-data[1][2]*data[2][1])
				  +data[0][1]*(data[1][2]*data[2][0]-data[1][0]*data[2][2])
				  +data[0][2]*(data[1][0]*data[2][1]-data[1][1]*data[2][0]);
		}

		MAT3 trans(void) {
			MAT3 temp;
			int i,j;
			for (i=0;i<3;i++)
				for (j=0;j<3;j++)
					temp[i][j] = data[j][i];
			return temp;
		}

		// Cramer's rule is used
		// Zero matrix is returned if non-reversible
		MAT3 inv(void);
};

MAT3 operator+ (const MAT3& input_1, const MAT3& input_2);
MAT3 operator+ (const MAT3& input_1);
MAT3 operator+ (const MAT3& input_v, const DType& input_s);
MAT3 operator+ (const DType& input_s, const MAT3& input_v);
MAT3 operator- (const MAT3& input_1, const MAT3& input_2);
MAT3 operator- (const MAT3& input_1);
MAT3 operator- (const MAT3& input_v, const DType& input_s);
MAT3 operator- (const DType& input_s, const MAT3& input_v);
MAT3 operator* (const DType& input_s, const MAT3& input_v);
MAT3 operator* (const MAT3& input_v, const DType& input_s);
MAT3 operator/ (const MAT3& input_v, const DType& input_s);
VEC3 operator* (const MAT3& mat_l, const VEC3& vec_r);
MAT3 operator* (const MAT3& mat_l, const MAT3& mat_r);

/* 4D vectors and matrices
   VEC4
   MAT4
   */

class VEC4 {
	private:
		DType data[4];
	public:
		void print(void) {
			printf("%10.20f\t%10.20f\t%10.20f\t%10.20f\n",data[0],data[1],data[2],data[3]); 
			return;
		}

		VEC4() {
			data[0]=0.0; 
			data[1]=0.0; 
			data[2]=0.0;
			data[3]=0.0;
		}
		VEC4(DType in_val) {
			data[0]=in_val;
			data[1]=in_val;
			data[2]=in_val;
			data[3]=in_val;
		}
		VEC4(DType in_1, DType in_2, DType in_3, DType in_4) {
			data[0]=in_1; 
			data[1]=in_2; 
			data[2]=in_3;
			data[3]=in_4;
		}
		VEC4(const VEC4& in_vec) {
			data[0] = in_vec[0];
			data[1] = in_vec[1];
			data[2] = in_vec[2];
			data[3] = in_vec[3];
		}

		DType& operator[] (int i) {
            assert( (i>=0) && (i<=3) );
			return data[i];
		}

		const DType& operator[] (int i) const {
            assert( (i>=0) && (i<=3) );
			return data[i];
		}

		void operator= (DType in_data) {
			data[0]=in_data; 
			data[1]=in_data; 
			data[2]=in_data;
		   	data[3]=in_data;	
			return;
		}
		void operator= (const VEC4& in_data) {
			data[0]=in_data[0];
			data[1]=in_data[1];
			data[2]=in_data[2];
			data[3]=in_data[3];
			return;
		}
		void operator+= (const VEC4& in_data) {
			data[0]+=in_data[0];
			data[1]+=in_data[1];
			data[2]+=in_data[2];
			data[3]+=in_data[3];
			return;
		}
		void operator+= (DType in_data) {
			data[0]+=in_data;
			data[1]+=in_data;
			data[2]+=in_data;
			data[3]+=in_data;
			return;
		}
		void operator-= (const VEC4& in_data) {
			data[0]-=in_data[0];
			data[1]-=in_data[1];
			data[2]-=in_data[2];
			data[3]-=in_data[3];
			return;
		}
		void operator-= (DType in_data) {
			data[0]-=in_data;
			data[1]-=in_data;
			data[2]-=in_data;
			data[3]-=in_data;
			return;
		}
		void operator*= (DType in_data) {
			data[0]*=in_data;
			data[1]*=in_data;
			data[2]*=in_data;
			data[3]*=in_data;
			return;
		}
		void operator/= (DType in_data) {
			data[0]/=in_data;
			data[1]/=in_data;
			data[2]/=in_data;
			data[3]/=in_data;
			return;
		}

		DType length(void) const {
			return sqrt(data[0]*data[0]+data[1]*data[1]
					   +data[2]*data[2]+data[3]*data[3]);
		}
};

VEC4 operator+ (const VEC4& input_1, const VEC4& input_2);
VEC4 operator+ (const VEC4& input_1);
VEC4 operator+ (const VEC4& input_v, const DType& input_s);
VEC4 operator+ (const DType& input_s, const VEC4& input_v);
VEC4 operator- (const VEC4& input_1, const VEC4& input_2);
VEC4 operator- (const VEC4& input_1);
VEC4 operator- (const VEC4& input_v, const DType& input_s);
VEC4 operator- (const DType& input_s, const VEC4& input_v);
VEC4 operator* (const DType& input_s, const VEC4& input_v);
VEC4 operator* (const VEC4& input_v, const DType& input_s);
VEC4 operator/ (const VEC4& input_v, const DType& input_s);

class MAT4 {
	private:
		// data[i] is the i-th row of the matrix
		VEC4 data[4];
	public:
		void print(void) {
			int i;
			for (i=0; i<4; i++)
				printf("%10.20f\t%10.20f\t%10.20f\t%10.20f\n",data[i][0],data[i][1],data[i][2],data[i][3]); 
			printf("\n");
			return;
		}

		MAT4() {
			int i,j;
			for (i=0; i<4; i++)
				for (j=0; j<4; j++)
					data[i][j]=0.0;
		}
		MAT4(DType in_val) {
			int i,j;
			for (i=0; i<4; i++)
				for (j=0; j<4; j++)
					data[i][j]=in_val;
		}
		MAT4(DType in_1, DType in_2, DType in_3, DType in_4) {
			MAT4(0.0);
			data[0][0]=in_1; 
			data[1][1]=in_2; 
			data[2][2]=in_3;
			data[3][3]=in_4;
		}
		MAT4(const VEC4& in_data) {
			MAT4(0.0);
			int i;
			for (i=0; i<4; i++)
				data[i][i]=in_data[i];
		}
		MAT4(const MAT4& in_mat) {
			int i,j;
			for (i=0; i<4; i++)
				for (j=0; j<4; j++)
					data[i][j] = in_mat[i][j];
		}

		VEC4& operator[] (int i) {
            assert( (i>=0) && (i<=3) );
			return data[i];
		}

		const VEC4& operator[] (int i) const {
            assert( (i>=0) && (i<=3) );
			return data[i];
		}

		void operator= (DType in_data) {
			int i,j;
			for (i=0;i<4;i++)
				for (j=0;j<4;j++)
					data[i][j] = in_data;
			return;
		}
		void operator= (const MAT4& in_data) {
			int i,j;
			for (i=0;i<4;i++)
				for (j=0;j<4;j++)
					data[i][j]=in_data[i][j];
			return;
		}
		void operator+= (const MAT4& in_data) {
			int i,j;
			for (i=0;i<4;i++)
				for (j=0;j<4;j++)
					data[i][j]+=in_data[i][j];
			return;
		}
		void operator+= (DType in_data) {
			int i,j;
			for (i=0;i<4;i++)
				for (j=0;j<4;j++)
					data[i][j]+=in_data;
			return;
		}
		void operator-= (const MAT4& in_data) {
			int i,j;
			for (i=0;i<4;i++)
				for (j=0;j<4;j++)
					data[i][j]-=in_data[i][j];
			return;
		}
		void operator-= (DType in_data) {
			int i,j;
			for (i=0;i<4;i++)
				for (j=0;j<4;j++)
					data[i][j]-=in_data;
			return;
		}
		void operator*= (DType in_data) {
			int i,j;
			for (i=0;i<4;i++)
				for (j=0;j<4;j++)
					data[i][j]*=in_data;
			return;
		}
		void operator/= (DType in_data) {
			int i,j;
			for (i=0;i<4;i++)
				for (j=0;j<4;j++)
					data[i][j]/=in_data;
			return;
		}

		DType det(void) {
			MAT3 temp;
			int indx[3];
			DType result = 0.0;
			DType mult = 0.0;
			int i,j,k;
			for (k=0; k<4; k++) {
				switch (k) {
					case 0:
						indx[0] = 1;	indx[1] = 2;	indx[2] = 3;
						mult = 1.0;
						break;
					case 1:	
						indx[0] = 0;	indx[1] = 2;	indx[2] = 3;
						mult = -1.0;
						break;
					case 2:
						indx[0] = 0;	indx[1] = 1;	indx[2] = 3;
						mult = 1.0;
						break;
					case 3:
						indx[0] = 0;	indx[1] = 1;	indx[2] = 2;
						mult = -1.0;
						break;
					default:
						break;
				}
				for (i=0; i<3; i++)
					for (j=0; j<3; j++)
						temp[i][j] = data[i][indx[j]];
				result += -mult*temp.det()*data[3][k];
			}		
		    return result;
		}

		MAT4 trans(void) {
			MAT4 temp;
			int i,j;
			for (i=0;i<4;i++)
				for (j=0;j<4;j++)
					temp[i][j] = data[j][i];
			return temp;
		}

		// Cramer's rule is used
		// Zero matrix is returned if non-reversible
		MAT4 inv(void);
};

MAT4 operator+ (const MAT4& input_1, const MAT4& input_2);
MAT4 operator+ (const MAT4& input_1);
MAT4 operator+ (const MAT4& input_v, const DType& input_s);
MAT4 operator+ (const DType& input_s, const MAT4& input_v);
MAT4 operator- (const MAT4& input_1, const MAT4& input_2);
MAT4 operator- (const MAT4& input_1);
MAT4 operator- (const MAT4& input_v, const DType& input_s);
MAT4 operator- (const DType& input_s, const MAT4& input_v);
MAT4 operator* (const DType& input_s, const MAT4& input_v);
MAT4 operator* (const MAT4& input_v, const DType& input_s);
MAT4 operator/ (const MAT4& input_v, const DType& input_s);
VEC4 operator* (const MAT4& mat_l, const VEC4& vec_r);
MAT4 operator* (const MAT4& mat_l, const MAT4& mat_r);

/* General dim-dimensional vectors and matrices
   VEC<d>
   MAT<d>
   */

template <int dim> class VEC {
	private:
		DType data[dim];
	public:
		void print(void) {
			for (int d=0; d<dim; d++) fprintf(stdout,"%10.20f\t",data[d]);
			fprintf(stdout,"\n"); return;
		}

		VEC() {for (int d=0; d<dim; d++) data[d] = 0.0;}

		VEC(DType in_val) {for (int d=0; d<dim; d++) data[d] = in_val;}

		VEC(const VEC<dim>& in_vec) {
			for (int d=0; d<dim; d++) data[d] = in_vec[d];}

		DType& operator[] (int i) {
            assert( (i>=0) && (i<dim) );
			return data[i];
		}

		const DType& operator[] (int i) const {
            assert( (i>=0) && (i<dim) );
			return data[i];
		}

		void operator= (DType in_data) {
			for (int d=0; d<dim; d++)
        data[d] = in_data; 
      return;
    }

		void operator= (const VEC<dim>& in_data) {
			for (int d=0; d<dim; d++) 
        data[d] = in_data[d]; 
      return;
    }

		void operator+= (const VEC<dim>& in_data) {
			for (int d=0; d<dim; d++) 
        data[d] += in_data[d]; 
      return;
    }

		void operator+= (DType in_data) {
			for (int d=0; d<dim; d++) 
        data[d] += in_data; 
      return;
    }

		void operator-= (const VEC<dim>& in_data) {
			for (int d=0; d<dim; d++) 
        data[d] -= in_data[d]; 
      return;
    }

		void operator-= (DType in_data) {
			for (int d=0; d<dim; d++) 
        data[d] -= in_data; 
      return;
    }

		void operator*= (DType in_data) {
			for (int d=0; d<dim; d++) 
        data[d] *= in_data; 
      return;
    }

		void operator/= (DType in_data) {
			for (int d=0; d<dim; d++) 
        data[d] /= in_data; 
      return;
    }

		int size(void) const {return dim;}

		DType length(void) const {
			DType entry, len;	len = 0.0;
			for (int d=0; d<dim; d++) {
				entry = data[d]; len += entry*entry;}
			return sqrt(len);
		}
};

template <int dim> VEC<dim> operator+ (const VEC<dim>&, const VEC<dim>&);
template <int dim> VEC<dim> operator+ (const VEC<dim>&);
template <int dim> VEC<dim> operator+ (const VEC<dim>&, const DType&);
template <int dim> VEC<dim> operator+ (const DType&, const VEC<dim>&);
template <int dim> VEC<dim> operator- (const VEC<dim>&, const VEC<dim>&);
template <int dim> VEC<dim> operator- (const VEC<dim>&);
template <int dim> VEC<dim> operator- (const VEC<dim>&, const DType&);
template <int dim> VEC<dim> operator- (const DType&, const VEC<dim>&);
template <int dim> VEC<dim> operator* (const DType&, const VEC<dim>&);
template <int dim> VEC<dim> operator* (const VEC<dim>&, const DType&);
template <int dim> VEC<dim> operator/ (const VEC<dim>&, const DType&);
// This is the inner product
template <int dim> DType operator* (const VEC<dim>&, const VEC<dim>&);
// This is the exterior (cross) product
VEC<3> operator^ (const VEC<3>&, const VEC<3>&);

template <int dim> class MAT {
	private:
		// data[i] is the i-th row of the matrix
		VEC<dim> data[dim];

		// Find the index of the entry with largest absolute value
		// in the sub-matrix [I,dim]X[I,dim]
		void findMatMax(int*,int*,MAT<dim>&,int);
		void findRowMax(int*,VEC<dim>&,int);
		// Permute two rows or two columns
		void permRow(MAT<dim>&,int,int);
		void permCol(MAT<dim>&,int,int);
	public:
		void print(void) {
			for (int i=0; i<dim; i++) {
				for (int j=0; j<dim; j++)
					fprintf(stdout,"%10.20f\t",data[i][j]);
				fprintf(stdout,"\n");
			}
			fprintf(stdout,"\n");
			return;
		}
		
		MAT() {
			for (int i=0; i<dim; i++) 
				for (int j=0; j<dim; j++) 
					data[i][j] = 0.0;
		}

		MAT(DType in_val) {
			for (int i=0; i<dim; i++)
				for (int j=0; j<dim; j++)
					data[i][j] = in_val;
		}

		MAT(const VEC<dim>& in_vec) {
			MAT(0.0);
			for (int i=0; i<dim; i++)
				data[i][i] = in_vec[i];
		}

		MAT(const MAT<dim>& in_mat) {
			for (int i=0; i<dim; i++)
				for (int j=0; j<dim; j++)
					data[i][j] = in_mat[i][j];
		}

		VEC<dim>& operator[] (int i) {
            assert( (i>=0) && (i<dim) );
			return data[i];
		}

		const VEC<dim>& operator[] (int i) const {
            assert( (i>=0) && (i<dim) );
			return data[i];
		}

		void operator= (DType in_data) {
			for (int i=0; i<dim; i++)
				for (int j=0; j<dim; j++)
					data[i][j] = in_data;
			return;
		}

		void operator= (const MAT<dim>& in_data) {
			for (int i=0; i<dim; i++)
				for (int j=0; j<dim; j++)
					data[i][j] = in_data[i][j];
			return;
		}

		void operator= (const VEC<dim>& in_vec) {
			MAT(0.0);
			for (int i=0; i<dim; i++)
				data[i][i] = in_vec[i];
			return;
		}
		
		void operator+= (const MAT<dim>& in_data) {
			for (int d=0; d<dim; d++) 
        data[d] += in_data[d]; 
      return;
    }

		void operator+= (DType in_data) {
			for (int d=0; d<dim; d++) 
        data[d] += in_data; 
      return;
    }

		void operator-= (const MAT<dim>& in_data) {
			for (int d=0; d<dim; d++) 
        data[d] -= in_data[d]; 
      return;
    }

		void operator-= (DType in_data) {
			for (int d=0; d<dim; d++) 
        data[d] -= in_data; 
      return;
    }

		void operator*= (DType in_data) {
			for (int d=0; d<dim; d++) 
        data[d] *= in_data; 
      return;
    }

		void operator/= (DType in_data) {
			for (int d=0; d<dim; d++) 
        data[d] /= in_data; 
      return;
    }


		MAT<dim> trans(void) {
			MAT<dim> mat;
			for (int i=0; i<dim; i++)
				for (int j=0; j<dim; j++)
					mat[i][j] = data[j][i];;
			return mat;
		}

		DType det_stable(void);	// PLU with largest entry pivoting
		DType det_fast(void);	// PLU with non-zero pivoting
		DType det(void) {return det_stable();} 

		MAT<dim> inv_stable(void);	// PLU with largest entry pivoting
		MAT<dim> inv_fast(void);	// PLU with non-zero pivoting
		MAT<dim> inv(void) {return inv_stable();} 
};

template <int dim> MAT<dim> operator+ (const MAT<dim>&, const MAT<dim>&);
template <int dim> MAT<dim> operator+ (const MAT<dim>&);
template <int dim> MAT<dim> operator+ (const MAT<dim>&, const DType&);
template <int dim> MAT<dim> operator+ (const DType&, const MAT<dim>&);
template <int dim> MAT<dim> operator- (const MAT<dim>&, const MAT<dim>&);
template <int dim> MAT<dim> operator- (const MAT<dim>&);
template <int dim> MAT<dim> operator- (const MAT<dim>&, const DType&);
template <int dim> MAT<dim> operator- (const DType&, const MAT<dim>&);
template <int dim> MAT<dim> operator* (const DType&, const MAT<dim>&);
template <int dim> MAT<dim> operator* (const MAT<dim>&, const DType&);
template <int dim> MAT<dim> operator/ (const MAT<dim>&, const DType&);
template <int dim> VEC<dim> operator* (const MAT<dim>&, const VEC<dim>&);
template <int dim> MAT<dim> operator* (const MAT<dim>&, const MAT<dim>&);
// This is the outer product
template <int dim> MAT<dim> operator| (const VEC<dim>&, const VEC<dim>&);
// This is the Kronecker product	
template <int dim1, int dim2>
MAT<dim1*dim2> operator^ (const MAT<dim1>&, const MAT<dim2>&);

//	Template function definitions

// Find the index of the entry with largest absolute value
template <int dim>
void MAT<dim>::findMatMax(int* r, int* c, MAT<dim>& M, int I) {
	int r_idx = I;
	int c_idx = I;
	DType max_val = abs(M[I][I]);
	for (int i=I; i<dim; i++)
		for (int j=I; j<dim; j++) {
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

template <int dim>
void MAT<dim>::findRowMax(int* c, VEC<dim>& V, int I) {
	int c_idx = I;
	DType max_val = abs(V[I]);
	for (int j=I; j<dim; j++) {
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
template <int dim>
void MAT<dim>::permRow(MAT<dim>& M, int r1, int r2) {
	if (r1==r2) return;
	DType r_tmp[dim];
	for (int j=0; j<dim; j++) r_tmp[j] = M[r1][j];
	for (int j=0; j<dim; j++) {
		M[r1][j] = M[r2][j];
		M[r2][j] = r_tmp[j];
	}
	return;
}

template <int dim>
void MAT<dim>::permCol(MAT<dim>& M, int c1, int c2) {
	if (c1==c2) return;
	DType c_tmp[dim];
	for (int i=0; i<dim; i++) c_tmp[i] = M[i][c1];
	for (int i=0; i<dim; i++) {
		M[i][c1] = M[i][c2];
		M[i][c2] = c_tmp[i];
	}
	return;
}

// Computed by PLU decomposition with largest entry pivoting
template <int dim>
DType MAT<dim>::det_stable(void) {
	MAT<dim> M(*this);
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
template <int dim>
DType MAT<dim>::det_fast(void) {
	MAT<dim> M(*this);
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
template <int dim>
MAT<dim> MAT<dim>::inv_stable(void) {
	MAT<dim> M(*this);
	int Pr[dim];	// Row permutation
	int Pc[dim];	// Column permutation
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
			//for (int j=i-1; j>-1; j--) { // Cause a strict-overflow warning
      for (int idx_j=0; idx_j<i; idx_j++) {
        int j = i-1-idx_j;
				DType tmp = 0.0;
				for (int k=j+1; k<i+1; k++) tmp -= M[j][k]*M[k][i];
				M[j][i] = tmp/M[j][j];
			}
			DType scal = (M[i][i]>0) ? 1.0 : -1.0;
			M[i][i] = abs(M[i][i]);
			//for (int j=i-1; j>-1; j--) { // Cause a strict-overflow warning
      for (int idx_j=0; idx_j<i; idx_j++) {
        int j = i-1-idx_j;
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
	return M;
}

// Computed by PLU decomposition with non-zero pivoting
template <int dim>
MAT<dim> MAT<dim>::inv_fast(void) {
	MAT<dim> M(*this);
	int Pr[dim];	// Row permutation
	int Pc[dim];	// Column permutation
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
	return M;
}

template <int dim> 
VEC<dim> operator+ (const VEC<dim>& vec_1, const VEC<dim>& vec_2) {
	VEC<dim> vec;
	for (int d=0; d<dim; d++) vec[d] = vec_1[d]+vec_2[d];
	return vec;
}

template <int dim> 
VEC<dim> operator+ (const VEC<dim>& vec_in) {
	VEC<dim> vec;
	for (int d=0; d<dim; d++) vec[d] = vec_in[d];
	return vec;
}

template <int dim> 
VEC<dim> operator+ (const VEC<dim>& vec_in, const DType& scal_in) {
	VEC<dim> vec;
	for (int d=0; d<dim; d++) vec[d] = vec_in[d]+scal_in;
	return vec;
}

template <int dim> 
VEC<dim> operator+ (const DType& scal_in, const VEC<dim>& vec_in) {
	VEC<dim> vec;
	for (int d=0; d<dim; d++) vec[d] = vec_in[d]+scal_in;
	return vec;
}

template <int dim> 
VEC<dim> operator- (const VEC<dim>& vec_1, const VEC<dim>& vec_2) {
	VEC<dim> vec;
	for (int d=0; d<dim; d++) vec[d] = vec_1[d]-vec_2[d];
	return vec;
}

template <int dim> 
VEC<dim> operator- (const VEC<dim>& vec_in) {
	VEC<dim> vec;
	for (int d=0; d<dim; d++) vec[d] = -vec_in[d];
	return vec;
}

template <int dim> 
VEC<dim> operator- (const VEC<dim>& vec_in, const DType& scal_in) {
	VEC<dim> vec;
	for (int d=0; d<dim; d++) vec[d] = vec_in[d]-scal_in;
	return vec;
}

template <int dim> 
VEC<dim> operator- (const DType& scal_in, const VEC<dim>& vec_in) {
	VEC<dim> vec;
	for (int d=0; d<dim; d++) vec[d] = scal_in-vec_in[d];
	return vec;
}

template <int dim> 
VEC<dim> operator* (const DType& sca_in, const VEC<dim>& vec_in) {
	VEC<dim> vec;
	for (int d=0; d<dim; d++) vec[d] = sca_in*vec_in[d];
	return vec;
}

template <int dim> 
VEC<dim> operator* (const VEC<dim>& vec_in, const DType& sca_in) {
	VEC<dim> vec;
	for (int d=0; d<dim; d++) vec[d] = sca_in*vec_in[d];
	return vec;
}

template <int dim> 
VEC<dim> operator/ (const VEC<dim>& vec_in, const DType& sca_in) {
	VEC<dim> vec;
	for (int d=0; d<dim; d++) vec[d] = vec_in[d]/sca_in;
	return vec;
}

// This is the inner product
template <int dim> 
DType operator* (const VEC<dim>& vec_1, const VEC<dim>& vec_2) {
	DType sca = 0.0;
	for (int d=0; d<dim; d++)
		sca += vec_1[d]*vec_2[d];
	return sca;
}

template <int dim> 
MAT<dim> operator+ (const MAT<dim>& mat_1, const MAT<dim>& mat_2) {
	MAT<dim> mat;
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] = mat_1[i][j]+mat_2[i][j];
	return mat;
}

template <int dim> 
MAT<dim> operator+ (const MAT<dim>& mat_in) {
	MAT<dim> mat;
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] = mat_in[i][j];
	return mat;
}

template <int dim> 
MAT<dim> operator+ (const MAT<dim>& mat_in, const DType& scal_in) {
	MAT<dim> mat;
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] = mat_in[i][j]+scal_in;
	return mat;
}

template <int dim> 
MAT<dim> operator+ (const DType& scal_in, const MAT<dim>& mat_in) {
	MAT<dim> mat;
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] = mat_in[i][j]+scal_in;
	return mat;
}

template <int dim> 
MAT<dim> operator- (const MAT<dim>& mat_1, const MAT<dim>& mat_2) {
	MAT<dim> mat;
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] = mat_1[i][j]-mat_2[i][j];
	return mat;
}

template <int dim> 
MAT<dim> operator- (const MAT<dim>& mat_in) {
	MAT<dim> mat;
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] = -mat_in[i][j];
	return mat;
}

template <int dim> 
MAT<dim> operator- (const MAT<dim>& mat_in, const DType& scal_in) {
	MAT<dim> mat;
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] = mat_in[i][j]-scal_in;
	return mat;
}

template <int dim> 
MAT<dim> operator- (const DType& scal_in, const MAT<dim>& mat_in) {
	MAT<dim> mat;
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] = scal_in-mat_in[i][j];
	return mat;
}

template <int dim> 
MAT<dim> operator* (const DType& sca_in, const MAT<dim>& mat_in) {
	MAT<dim> mat;
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] = sca_in*mat_in[i][j];
	return mat;
}

template <int dim> 
MAT<dim> operator* (const MAT<dim>& mat_in, const DType& sca_in) {
	MAT<dim> mat;
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] = sca_in*mat_in[i][j];
	return mat;
}

template <int dim> 
MAT<dim> operator/ (const MAT<dim>& mat_in, const DType& sca_in) {
	MAT<dim> mat;
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat[i][j] = mat_in[i][j]/sca_in;
	return mat;
}

template <int dim> 
VEC<dim> operator* (const MAT<dim>& mat_in, const VEC<dim>& vec_in) {
	VEC<dim> vec(0.0);
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			vec[i] += mat_in[i][j]*vec_in[j];
	return vec;
}

template <int dim> 
MAT<dim> operator* (const MAT<dim>& mat_1, const MAT<dim>& mat_2) {
	MAT<dim> mat(0.0);
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			for (int d=0; d<dim; d++)
				mat[i][j] += mat_1[i][d]*mat_2[d][j];
	return mat;
}

// This is the outer product
template <int dim> 
MAT<dim> operator| (const VEC<dim>& vec_1, const VEC<dim>& vec_2) {
	MAT<dim> mat;
	for (int i=0; i<dim; i++) {
		DType sca = vec_1[i];
		for (int j=0; j<dim; j++)
			mat[i][j] = sca*vec_2[j];
	}
	return mat;
}

// This is the Kronecker product
template <int dim1, int dim2>
MAT<dim1*dim2> operator^ (const MAT<dim1>& mat_1, const MAT<dim2>& mat_2) 
{
	MAT<dim1*dim2> mat;
	for (int i1=0; i1<dim1; i1++) {
		int I = i1*dim1;
		for (int j1=0; j1<dim1; j1++) {
			int J = j1*dim1;
			DType sca = mat_1[i1][j1];
			for (int i2=0; i2<dim2; i2++)
				for (int j2=0; j2<dim2; j2++)
					mat[I+i2][J+j2] = sca*mat_2[i2][j2];
		}
	}
	return mat;
}

#endif

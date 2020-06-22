#include "LowDim_Alg.h"

/* 3D vectors and matrices
   VEC3
   MAT3
   */
MAT3 MAT3::inv(void) {
	DType div = det();
	MAT3 temp;
	if (abs(div)<TOL_ZERO)
		return temp;
	temp[0][0] = (data[1][1]*data[2][2]-data[1][2]*data[2][1])/div;
	temp[1][0] = (data[1][2]*data[2][0]-data[1][0]*data[2][2])/div;
	temp[2][0] = (data[1][0]*data[2][1]-data[1][1]*data[2][0])/div;

    temp[0][1] = (data[2][1]*data[0][2]-data[2][2]*data[0][1])/div;
	temp[1][1] = (data[2][2]*data[0][0]-data[2][0]*data[0][2])/div;
	temp[2][1] = (data[2][0]*data[0][1]-data[2][1]*data[0][0])/div;
	
	temp[0][2] = (data[0][1]*data[1][2]-data[0][2]*data[1][1])/div;
	temp[1][2] = (data[0][2]*data[1][0]-data[0][0]*data[1][2])/div;
	temp[2][2] = (data[0][0]*data[1][1]-data[0][1]*data[1][0])/div;
	return temp;
}

VEC3 operator+ (const VEC3& input_1, const VEC3& input_2) {
	return VEC3(input_1[0]+input_2[0], 
				input_1[1]+input_2[1],
				input_1[2]+input_2[2]);
}

VEC3 operator+ (const VEC3& input_1) {
	return VEC3(input_1[0],input_1[1],input_1[2]);
}

VEC3 operator+ (const VEC3& input_v, const DType& input_s) {
	return VEC3(input_v[0]+input_s,
				input_v[1]+input_s,
				input_v[2]+input_s);
}

VEC3 operator+ (const DType& input_s, const VEC3& input_v) {
	return VEC3(input_v[0]+input_s,
				input_v[1]+input_s,
				input_v[2]+input_s);
}

VEC3 operator- (const VEC3& input_1, const VEC3& input_2) {
	return VEC3(input_1[0]-input_2[0], 
				input_1[1]-input_2[1],
				input_1[2]-input_2[2]);
}

VEC3 operator- (const VEC3& input_1) {
	return VEC3(-input_1[0],-input_1[1],-input_1[2]);
}

VEC3 operator- (const VEC3& input_v, const DType& input_s) {
	return VEC3(input_v[0]-input_s,
				input_v[1]-input_s,
				input_v[2]-input_s);
}

VEC3 operator- (const DType& input_s, const VEC3& input_v) {
	return VEC3(input_s-input_v[0],
				input_s-input_v[1],
				input_s-input_v[2]);
}

VEC3 operator* (const DType& input_s, const VEC3& input_v) {
	return VEC3(input_s*input_v[0],
				input_s*input_v[1],
				input_s*input_v[2]);
}

VEC3 operator* (const VEC3& input_v, const DType& input_s) {
	return VEC3(input_s*input_v[0],
				input_s*input_v[1],
				input_s*input_v[2]);
}

VEC3 operator/ (const VEC3& input_v, const DType& input_s) {
	return VEC3(input_v[0]/input_s,
				input_v[1]/input_s,
				input_v[2]/input_s);
}

MAT3 operator+ (const MAT3& input_1, const MAT3& input_2) {
	MAT3 temp;
	int i,j;
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			temp[i][j]=input_1[i][j]+input_2[i][j];
	return temp;
}

MAT3 operator+ (const MAT3& input_1) {
	MAT3 temp;
	int i,j;
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			temp[i][j]=input_1[i][j];
	return temp;
}

MAT3 operator+ (const MAT3& input_v, const DType& input_s) {
	MAT3 temp;
	int i,j;
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			temp[i][j]=input_v[i][j]+input_s;
	return temp;
}

MAT3 operator+ (const DType& input_s, const MAT3& input_v) {
	MAT3 temp;
	int i,j;
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			temp[i][j]=input_v[i][j]+input_s;
	return temp;
}

MAT3 operator- (const MAT3& input_1, const MAT3& input_2) {
	MAT3 temp;
	int i,j;
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			temp[i][j]=input_1[i][j]-input_2[i][j];
	return temp;
}

MAT3 operator- (const MAT3& input_1) {
	MAT3 temp;
	int i,j;
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			temp[i][j]=-input_1[i][j];
	return temp;
}

MAT3 operator- (const MAT3& input_v, const DType& input_s) {
	MAT3 temp;
	int i,j;
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			temp[i][j]=input_v[i][j]-input_s;
	return temp;
}

MAT3 operator- (const DType& input_s, const MAT3& input_v) {
	MAT3 temp;
	int i,j;
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			temp[i][j]=input_s-input_v[i][j];
	return temp;
}

MAT3 operator* (const DType& input_s, const MAT3& input_v) {
	MAT3 temp;
	int i,j;
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			temp[i][j] = input_s*input_v[i][j];
	return temp;
}

MAT3 operator* (const MAT3& input_v, const DType& input_s) {
	MAT3 temp;
	int i,j;
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			temp[i][j] = input_s*input_v[i][j];
	return temp;
}

MAT3 operator/ (const MAT3& input_v, const DType& input_s) {
	MAT3 temp;
	int i,j;
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			temp[i][j] = input_v[i][j]/input_s;
	return temp;
}

VEC3 operator* (const MAT3& mat_l, const VEC3& vec_r) {
	VEC3 temp;
	int i,j;
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			temp[i] += mat_l[i][j]*vec_r[j];
	return temp;
}

MAT3 operator* (const MAT3& mat_l, const MAT3& mat_r) {
	MAT3 temp;
	int i,j,k;
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			for (k=0; k<3; k++)
				temp[i][j] += mat_l[i][k]*mat_r[k][j];
	return temp;
}

/* 4D vectors and matrices
   VEC4
   MAT4
   */
MAT4 MAT4::inv(void) {
	DType div = det();
	MAT4 trgt;
	if (abs(div)<TOL_ZERO)
		return trgt;
	DType mult  = 0.0;
	DType mult2 = 0.0;
	int idxi[3];
	int idxj[3];
	MAT3 temp;
	int i,j,k,l;
	for (i=0; i<4; i++) {
		switch (i) {
			case 0:
				idxi[0] = 1;	idxi[1] = 2;	idxi[2] = 3;
				mult = 1.0;
				break;
			case 1:
				idxi[0] = 0;	idxi[1] = 2;	idxi[2] = 3;
				mult = -1.0;
				break;
			case 2:
				idxi[0] = 0;	idxi[1] = 1;	idxi[2] = 3;
				mult = 1.0;
				break;
			case 3:
				idxi[0] = 0;	idxi[1] = 1;	idxi[2] = 2;
				mult = -1.0;
				break;
			default:
				break;
		}
    	for (j=0; j<4; j++) {
			switch (j) {
				case 0:
					idxj[0] = 1;	idxj[1] = 2;	idxj[2] = 3;
					mult2 = mult;
					break;
				case 1:
					idxj[0] = 0;	idxj[1] = 2;	idxj[2] = 3;
					mult2 = -mult;
					break;
				case 2:
					idxj[0] = 0;	idxj[1] = 1;	idxj[2] = 3;
					mult2 = mult;
					break;
				case 3:
					idxj[0] = 0;	idxj[1] = 1;	idxj[2] = 2;
					mult2 = -mult;
					break;
				default:
					break;
			}
			for (k=0; k<3; k++)
				for (l=0; l<3; l++)
					temp[k][l] = data[idxi[k]][idxj[l]];
			trgt[j][i] = mult2*temp.det()/div;
		}
	}
	return trgt;
}

VEC4 operator+ (const VEC4& input_1, const VEC4& input_2) {
	return VEC4(input_1[0]+input_2[0], 
				input_1[1]+input_2[1],
				input_1[2]+input_2[2],
				input_1[3]+input_2[3]);
}

VEC4 operator+ (const VEC4& input_1) {
	return VEC4(input_1[0],input_1[1],input_1[2],input_1[3]);
}

VEC4 operator+ (const VEC4& input_v, const DType& input_s) {
	return VEC4(input_v[0]+input_s,
				input_v[1]+input_s,
				input_v[2]+input_s,
				input_v[3]+input_s);
}

VEC4 operator+ (const DType& input_s, const VEC4& input_v) {
	return VEC4(input_v[0]+input_s,
				input_v[1]+input_s,
				input_v[2]+input_s,
				input_v[3]+input_s);
}

VEC4 operator- (const VEC4& input_1, const VEC4& input_2) {
	return VEC4(input_1[0]-input_2[0], 
				input_1[1]-input_2[1],
				input_1[2]-input_2[2],
				input_1[3]-input_2[3]);
}

VEC4 operator- (const VEC4& input_1) {
	return VEC4(-input_1[0],-input_1[1],-input_1[2],-input_1[3]);
}

VEC4 operator- (const VEC4& input_v, const DType& input_s) {
	return VEC4(input_v[0]-input_s,
				input_v[1]-input_s,
				input_v[2]-input_s,
				input_v[3]-input_s);
}

VEC4 operator- (const DType& input_s, const VEC4& input_v) {
	return VEC4(input_s-input_v[0],
				input_s-input_v[1],
				input_s-input_v[2],
				input_s-input_v[3]);
}

VEC4 operator* (const DType& input_s, const VEC4& input_v) {
	return VEC4(input_s*input_v[0],
				input_s*input_v[1],
				input_s*input_v[2],
				input_s*input_v[3]);
}

VEC4 operator* (const VEC4& input_v, const DType& input_s) {
	return VEC4(input_s*input_v[0],
				input_s*input_v[1],
				input_s*input_v[2],
				input_s*input_v[3]);
}

VEC4 operator/ (const VEC4& input_v, const DType& input_s) {
	return VEC4(input_v[0]/input_s,
				input_v[1]/input_s,
				input_v[2]/input_s,
				input_v[3]/input_s);
}

MAT4 operator+ (const MAT4& input_1, const MAT4& input_2) {
	MAT4 temp;
	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			temp[i][j]=input_1[i][j]+input_2[i][j];
	return temp;
}

MAT4 operator+ (const MAT4& input_1) {
	MAT4 temp;
	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			temp[i][j]=input_1[i][j];
	return temp;
}

MAT4 operator+ (const MAT4& input_v, const DType& input_s) {
	MAT4 temp;
	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			temp[i][j]=input_v[i][j]+input_s;
	return temp;
}

MAT4 operator+ (const DType& input_s, const MAT4& input_v) {
	MAT4 temp;
	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			temp[i][j]=input_v[i][j]+input_s;
	return temp;
}

MAT4 operator- (const MAT4& input_1, const MAT4& input_2) {
	MAT4 temp;
	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			temp[i][j]=input_1[i][j]-input_2[i][j];
	return temp;
}

MAT4 operator- (const MAT4& input_1) {
	MAT4 temp;
	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			temp[i][j]=-input_1[i][j];
	return temp;
}

MAT4 operator- (const MAT4& input_v, const DType& input_s) {
	MAT4 temp;
	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			temp[i][j]=input_v[i][j]-input_s;
	return temp;
}

MAT4 operator- (const DType& input_s, const MAT4& input_v) {
	MAT4 temp;
	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			temp[i][j]=input_s-input_v[i][j];
	return temp;
}

MAT4 operator* (const DType& input_s, const MAT4& input_v) {
	MAT4 temp;
	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			temp[i][j] = input_s*input_v[i][j];
	return temp;
}

MAT4 operator* (const MAT4& input_v, const DType& input_s) {
	MAT4 temp;
	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			temp[i][j] = input_s*input_v[i][j];
	return temp;
}

MAT4 operator/ (const MAT4& input_v, const DType& input_s) {
	MAT4 temp;
	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			temp[i][j] = input_v[i][j]/input_s;
	return temp;
}

VEC4 operator* (const MAT4& mat_l, const VEC4& vec_r) {
	VEC4 temp;
	int i,j;
	for (i=0; i<4; i++)
		for (j=0; j<4; j++)
			temp[i] += mat_l[i][j]*vec_r[j];
	return temp;
}

MAT4 operator* (const MAT4& mat_l, const MAT4& mat_r) {
	MAT4 temp;
	int i,j,k;
	for (i=0; i<4; i++)
		for (j=0; j<4; j++)
			for (k=0; k<4; k++)
				temp[i][j] += mat_l[i][k]*mat_r[k][j];
	return temp;
}

/* General dim-dimensional vectors and matrices
   VEC<d>
   MAT<d>
   */

// This is the exterior (cross) product
VEC<3> operator^ (const VEC<3>& vec_1, const VEC<3>& vec_2) {
	VEC<3> vec;
	vec[0] = vec_1[1]*vec_2[2]-vec_1[2]*vec_2[1];
	vec[1] = vec_1[2]*vec_2[0]-vec_1[0]*vec_2[2];
	vec[2] = vec_1[0]*vec_2[1]-vec_1[1]*vec_2[0];
	return vec;
}

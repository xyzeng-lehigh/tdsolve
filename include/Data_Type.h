#ifndef DATA_TYPE_H
#define DATA_TYPE_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cfloat>
using namespace std;

#ifndef FLOAT_ALG
#ifndef DType
#define DType double
#define DTYPE_MAX DBL_MAX
#endif
#else
#ifndef DType
#define DType float
#define DTYPE_MAX FLT_MAX
#endif
#endif

#ifndef TOL_ZERO
#define TOL_ZERO 1e-12
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#endif

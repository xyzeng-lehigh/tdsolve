#ifndef ALLOC_UTILS_H
#define ALLOC_UTILS_H

#include "Data_Type.h"
#include "HighDim_Alg.h"

/*
 * DType
 */

bool erase_array(DType*& target); 
DType* alloc_array_full(DType*& target, int dim1); 
DType* alloc_array_semi(DType*& target); 

bool erase_array(DType**& target, int dim1); 
DType** alloc_array_full(DType**& target, int dim1, int dim2); 
DType** alloc_array_semi(DType**& target, int dim1); 

bool erase_array(DType***& target, int dim1, int dim2); 
DType*** alloc_array_full(DType***& target, int dim1, int dim2, int dim3); 
DType*** alloc_array_semi(DType***& target, int dim1, int dim2); 

bool erase_array(DType****& target, int dim1, int dim2, int dim3); 
DType**** alloc_array_full(DType****& target, int dim1, int dim2, int dim3, int dim4); 
DType**** alloc_array_semi(DType****& target, int dim1, int dim2, int dim3); 


/*
 * int
 */

bool erase_array(int*& target); 
int* alloc_array_full(int*& target, int dim1); 
int* alloc_array_semi(int*& target); 

bool erase_array(int**& target, int dim1); 
int** alloc_array_full(int**& target, int dim1, int dim2); 
int** alloc_array_semi(int**& target, int dim1); 

bool erase_array(int***& target, int dim1, int dim2); 
int*** alloc_array_full(int***& target, int dim1, int dim2, int dim3); 
int*** alloc_array_semi(int***& target, int dim1, int dim2); 

bool erase_array(int****& target, int dim1, int dim2, int dim3); 
int**** alloc_array_full(int****& target, int dim1, int dim2, int dim3, int dim4); 
int**** alloc_array_semi(int****& target, int dim1, int dim2, int dim3); 

/*
 * HVEC
 */

bool erase_array(HVEC**& target, int dim1);
HVEC** alloc_array_full(HVEC**& target, int dim, int dim1);
bool reset_array(HVEC**& target, int dim1, DType val=0.0);

bool erase_array(HVEC***& target, int dim1, int dim2);
HVEC*** alloc_array_full(HVEC***& target, int dim, int dim1, int dim2);
bool reset_array(HVEC***& target, int dim1, int dim2, DType val=0.0);

bool erase_array(HVEC****& target, int dim1, int dim2, int dim3);
HVEC**** alloc_array_full(HVEC****& target, int dim, int dim1, int dim2, int dim3);
bool reset_array(HVEC****& target, int dim1, int dim2, int dim3, DType val=0.0);

#endif

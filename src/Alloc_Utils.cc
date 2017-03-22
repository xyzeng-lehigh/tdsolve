#include "Alloc_Utils.h"
using namespace std;

/*
 * DType
 */

bool erase_array(DType*& target) {
    if (target) {
        delete[] target;
        target = NULL;
        return true;
    }
    return false;
}

DType* alloc_array_full(DType*& target, int dim1) {
    target = new DType[dim1];
    return target;
}

DType* alloc_array_semi(DType*& target) {
    return target;
}


bool erase_array(DType**& target, int dim1) {
    if (target) {
        for (int d1=0; d1<dim1; d1++)
            if (target[d1]) delete[] target[d1];
        delete[] target;
        target = NULL;
        return true;
    }
    return false;
}

DType** alloc_array_full(DType**& target, int dim1, int dim2) {
    target = new DType*[dim1];
    for (int d1=0; d1<dim1; d1++)
        target[d1] = new DType[dim2];
    return target;
}

DType** alloc_array_semi(DType**& target, int dim1) {
    target = new DType*[dim1];
    return target;
}


bool erase_array(DType***& target, int dim1, int dim2) {
    if (target) {
        for (int d1=0; d1<dim1; d1++)
            if (target[d1]) {
                for (int d2=0; d2<dim2; d2++)
                    if (target[d1][d2]) delete[] target[d1][d2];
                delete[] target[d1];
            }
        delete[] target;
        target = NULL;
        return true;
    }
    return false;
}

DType*** alloc_array_full(DType***& target, int dim1, int dim2, int dim3) {
    target = new DType**[dim1];
    for (int d1=0; d1<dim1; d1++) {
        target[d1] = new DType*[dim2];
        for (int d2=0; d2<dim2; d2++)
            target[d1][d2] = new DType[dim3];
    }
    return target;
}

DType*** alloc_array_semi(DType***& target, int dim1, int dim2) {
    target = new DType**[dim1];
    for (int d1=0; d1<dim1; d1++)
        target[d1] = new DType*[dim2];
    return target;
}


bool erase_array(DType****& target, int dim1, int dim2, int dim3) {
    if (target) {
        for (int d1=0; d1<dim1; d1++)
            if (target[d1]) {
                for (int d2=0; d2<dim2; d2++) 
                    if (target[d1][d2]) {
                        for (int d3=0; d3<dim3; d3++)
                            if (target[d1][d2][d3]) delete[] target[d1][d2][d3];
                        delete[] target[d1][d2];
                    }
                delete[] target[d1];
            }
        delete[] target;
        target = NULL;
        return true;
    }
    return false;
}

DType**** alloc_array_full(DType****& target, int dim1, int dim2, int dim3, int dim4) {
    target = new DType***[dim1];
    for (int d1=0; d1<dim1; d1++) {
        target[d1] = new DType**[dim2];
        for (int d2=0; d2<dim2; d2++) {
            target[d1][d2] = new DType*[dim3];
            for (int d3=0; d3<dim3; d3++)
                target[d1][d2][d3] = new DType[dim4];
        }
    }
    return target;
}

DType**** alloc_array_semi(DType****& target, int dim1, int dim2, int dim3) {
    target = new DType***[dim1];
    for (int d1=0; d1<dim1; d1++) {
        target[d1] = new DType**[dim2];
        for (int d2=0; d2<dim2; d2++)
            target[d1][d2] = new DType*[dim3];
    }
    return target;
}



/*
 * int
 */

bool erase_array(int*& target) {
    if (target) {
        delete[] target;
        target = NULL;
        return true;
    }
    return false;
}

int* alloc_array_full(int*& target, int dim1) {
    target = new int[dim1];
    return target;
}

int* alloc_array_semi(int*& target) {
    return target;
}


bool erase_array(int**& target, int dim1) {
    if (target) {
        for (int d1=0; d1<dim1; d1++)
            if (target[d1]) delete[] target[d1];
        delete[] target;
        target = NULL;
        return true;
    }
    return false;
}

int** alloc_array_full(int**& target, int dim1, int dim2) {
    target = new int*[dim1];
    for (int d1=0; d1<dim1; d1++)
        target[d1] = new int[dim2];
    return target;
}

int** alloc_array_semi(int**& target, int dim1) {
    target = new int*[dim1];
    return target;
}


bool erase_array(int***& target, int dim1, int dim2) {
    if (target) {
        for (int d1=0; d1<dim1; d1++) 
            if (target[d1]) {
                for (int d2=0; d2<dim2; d2++)
                    if (target[d1][d2]) delete[] target[d1][d2];
                delete[] target[d1];
            }
        delete[] target;
        target = NULL;
        return true;
    }
    return false;
}

int*** alloc_array_full(int***& target, int dim1, int dim2, int dim3) {
    target = new int**[dim1];
    for (int d1=0; d1<dim1; d1++) {
        target[d1] = new int*[dim2];
        for (int d2=0; d2<dim2; d2++)
            target[d1][d2] = new int[dim3];
    }
    return target;
}

int*** alloc_array_semi(int***& target, int dim1, int dim2) {
    target = new int**[dim1];
    for (int d1=0; d1<dim1; d1++)
        target[d1] = new int*[dim2];
    return target;
}


bool erase_array(int****& target, int dim1, int dim2, int dim3) {
    if (target) {
        for (int d1=0; d1<dim1; d1++)
            if (target[d1]) {
                for (int d2=0; d2<dim2; d2++)
                    if (target[d1][d2]) {
                        for (int d3=0; d3<dim3; d3++)
                            if (target[d1][d2][d3]) delete[] target[d1][d2][d3];
                        delete[] target[d1][d2];
                    }
                delete[] target[d1];
            }
        delete[] target;
        target = NULL;
        return true;
    }
    return false;
}

int**** alloc_array_full(int****& target, int dim1, int dim2, int dim3, int dim4) {
    target = new int***[dim1];
    for (int d1=0; d1<dim1; d1++) {
        target[d1] = new int**[dim2];
        for (int d2=0; d2<dim2; d2++) {
            target[d1][d2] = new int*[dim3];
            for (int d3=0; d3<dim3; d3++)
                target[d1][d2][d3] = new int[dim4];
        }
    }
    return target;
}

int**** alloc_array_semi(int****& target, int dim1, int dim2, int dim3) {
    target = new int***[dim1];
    for (int d1=0; d1<dim1; d1++) {
        target[d1] = new int**[dim2];
        for (int d2=0; d2<dim2; d2++)
            target[d1][d2] = new int*[dim3];
    }
    return target;
}


/*
 * HVEC
 */

bool erase_array(HVEC**& target, int dim1) {
    if (target) {
        for (int i=0; i<dim1; i++)
            if (target[i]) delete target[i];
        delete[] target;
        return true;
    }
    return false;
}

HVEC** alloc_array_full(HVEC**& target, int dim, int dim1) {
    target = new HVEC*[dim1];
    for (int i=0; i<dim1; i++)
        target[i] = new HVEC(dim);
    return target;
}

bool reset_array(HVEC**& target, int dim1, DType val) {
    if (target) {
        for (int i=0; i<dim1; i++)
            *target[i] = val;
        return true;
    }

    return false;
}

bool erase_array(HVEC***& target, int dim1, int dim2) {
    if (target) {
        for (int i=0; i<dim1; i++)
            if (target[i]) {
                for (int j=0; j<dim2; j++)
                    if (target[i][j]) delete target[i][j];
                delete[] target[i];
            }
        delete[] target;
        return true;
    }
    return false;
}

HVEC*** alloc_array_full(HVEC***& target, int dim, int dim1, int dim2) {
    target = new HVEC**[dim1];
    for (int i=0; i<dim1; i++) {
        target[i] = new HVEC*[dim2];
        for (int j=0; j<dim2; j++)
            target[i][j] = new HVEC(dim);
    }
    return target;
}

bool reset_array(HVEC***& target, int dim1, int dim2, DType val) {
    if (target) {
        for (int i=0; i<dim1; i++)
            for (int j=0; j<dim2; j++)
                *target[i][j] = val;
        return true;
    }

    return false;
}

bool erase_array(HVEC****& target, int dim1, int dim2, int dim3) {
    if (target) {
        for (int i=0; i<dim1; i++)
            if (target[i]) {
                for (int j=0; j<dim2; j++)
                    if (target[i][j]) {
                        for (int k=0; k<dim3; k++)
                            if (target[i][j][k])
                                delete target[i][j][k];
                        delete[] target[i][j];
                    }
                delete[] target[i];
            }
        delete[] target;
        return true;
    }
    return false;
}

HVEC**** alloc_array_full(HVEC****& target, int dim, int dim1, int dim2, int dim3) {
    for (int i=0; i<dim1; i++) {
        target[i] = new HVEC**[dim2];
        for (int j=0; j<dim2; j++) {
            target[i][j] = new HVEC*[dim3];
            for (int k=0; k<dim3; k++)
                target[i][j][k] = new HVEC(dim);
        }
    }
    return target;
}

bool reset_array(HVEC****& target, int dim1, int dim2, int dim3, DType val) {
    if (target) {
        for (int i=0; i<dim1; i++)
            for (int j=0; j<dim2; j++)
                for (int k=0; k<dim3; k++)
                    *target[i][j][k] = val;
        return true;
    }

    return false;
}


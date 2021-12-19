#ifndef MAT_ENTRY_H
#define MAT_ENTRY_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include "Data_Type.h"
#include "LowDim_Alg.h"
#include "HighDim_Alg.h"
using namespace std;

// Block dim X dim entry 
class Entry {
    protected:
        int     dim;
        HMAT*   data;
    public:
        // Can have several use, for example to check
        // if an entry is zero for sparse matrices
        bool    sig;

        // Zero matrix
        Entry(int dim_in) {
            dim  = dim_in;
            data = new HMAT(dim,0.0);
            sig  = false;
        }
        // Multiple of identity matrix
        Entry(int dim_in, DType in_val) {
            dim  = dim_in;
            data = new HMAT(HVEC(dim,in_val));
            sig  = false;
        }
        // A general matrix
        Entry(const HMAT& mat) {
            dim  = mat.size();
            data = new HMAT(mat);
            sig  = false;
        }
        // An entry
        Entry(const Entry& ent) {
            dim  = ent.getSize();
            data = new HMAT(ent.getVal());
            sig  = ent.sig;
        }
        ~Entry(void) {
            if (data!=NULL)
                delete data;
        }

        int getSize(void) const {
            return dim;}
        void setSize(int dim_in) {
            dim = dim_in;
            if (data!=NULL) delete data;
            data = new HMAT(dim,0.0);
            return;
        }
        void print(void) {
            data->print(); return;}

        DType& get(int i, int j) {
            return (*data)[i][j];
        }

        const DType& get(int i, int j) const {
            return (*data)[i][j];
        }

        Entry trans(void) {
            Entry entry(data->trans());
            return entry;
        }
        const Entry trans(void) const {
            Entry entry(data->trans());
            return entry;
        }

        DType getNorm(void) const {
            return data->norm('f',2);}

        HMAT& getVal(void) {
            return *data;}
        const HMAT& getVal(void) const {
            return *data;}
        HMAT getInv(void) {
            return data->inv();}
        Entry& operator= (const HMAT& in_data) {
            assert(dim==in_data.size());
            *data = in_data; 
            return *this;
        }
        Entry& operator= (const Entry& in_ent) {
            assert(dim==in_ent.getSize());
            *data = in_ent.getVal();
            sig   = in_ent.sig;
            return *this;
        }
        Entry& operator= (const DType& val) {
            *data = 0.0;
            for (int i=0; i<dim; i++)
                (*data)[i][i] = val;
            return *this;
        }
        Entry& operator+= (const Entry& in_ent) {
            assert(dim==in_ent.getSize());
            *data += in_ent.getVal();
            return *this;
        }
        Entry& operator+= (const DType& val) {
            for (int i=0; i<dim; i++)
                (*data)[i][i] += val;
            return *this;
        }
        Entry& operator-= (const Entry& in_ent) {
            assert(dim==in_ent.getSize());
            *data -= in_ent.getVal();
            return *this;
        }
        Entry& operator-= (const DType& val) {
            for (int i=0; i<dim; i++)
                (*data)[i][i] -= val;
            return *this;
        }
        Entry& operator*= (const Entry& in_ent) {
            assert(dim==in_ent.getSize());
            *data = *data*in_ent.getVal(); 
            return *this;
        }
        Entry& operator*= (const DType in_val) {
            *data *= in_val; 
            return *this;
        }
        Entry& operator/= (const DType in_val) {
            *data /= in_val; 
            return *this;
        }
};

Entry operator+ (const Entry& ent_1, const Entry& ent_2); 
Entry operator+ (const Entry& ent_1, const DType& val); 
Entry operator+ (const DType& val, const Entry& ent_2); 
Entry operator+ (const Entry& ent_1); 
Entry operator- (const Entry& ent_1, const Entry& ent_2);
Entry operator- (const Entry& ent_1, DType val); 
Entry operator- (DType val, const Entry& ent_2); 
Entry operator- (const Entry& ent_1); 
Entry operator* (const Entry& ent_1, const Entry& ent_2); 
Entry operator* (const DType& sca_1, const Entry& ent_2); 
Entry operator* (const Entry& ent_1, const DType& sca_2);
Entry operator/ (const Entry& ent_1, const DType& sca_2); 

// Block dim X 1 component
class VComp {
    protected:
        int     dim;
        HVEC*   data;
    public:
        // Can have several use
        bool    sig;

        // Multiple of ones
        VComp(int dim_in, DType in_val=0.0) {
            dim  = dim_in;
            data = new HVEC(dim,0.0);
            sig  = false;
        }
        // A general vector
        VComp(const HVEC& vec) {
            dim  = vec.size();
            data = new HVEC(vec);
            sig  = false;
        }
        // A vector component
        VComp(const VComp& comp) {
            dim  = comp.getSize();
            data = new HVEC(comp.getVal());
            sig  = comp.sig;
        }
        ~VComp(void) {
            if (data!=NULL)
                delete data;
        }

        int getSize(void) const {
            return dim;}
        void setSize(int dim_in) {
            dim = dim_in;
            if (data!=NULL) delete data;
            data = new HVEC(dim,0.0);
            return;
        }
        void print(void) {
            data->print(); return;}

        DType& get(int i) {
            return (*data)[i];
        }
        const DType& get(int i) const {
            return (*data)[i];
        }
        HVEC& getVal(void) {
            return *data;}
        const HVEC& getVal(void) const {
            return *data;}
        DType getNorm(void) {
            return data->norm();}
        DType getInfNorm(void) {
            return data->norm(0);}
        VComp& operator= (const HVEC& in_data) {
            assert(dim==in_data.size());
            *data = in_data;
            return *this;
        }
        VComp& operator= (const VComp& in_comp) {
            assert(dim==in_comp.getSize());
            *data = in_comp.getVal();
            sig   = in_comp.sig;
            return *this;
        }
        VComp& operator= (const DType& val) {
            *data = val;
            return *this;
        }
        VComp& operator+= (const VComp& in_comp) {
            assert(dim==in_comp.getSize());
            *data += in_comp.getVal();
            return *this;
        }
        VComp& operator+= (const DType& val) {
            *data += val;
            return *this;
        }
        VComp& operator-= (const VComp& in_comp) {
            assert(dim==in_comp.getSize());
            *data -= in_comp.getVal();
            return *this;
        }
        VComp& operator-= (const DType& val) {
            *data -= val;
            return *this;
        }
        VComp& operator*= (const DType& val) {
            *data *= val;
            return *this;
        }
        VComp& operator/= (const DType& val) {
            *data /= val;
            return *this;
        }
};

VComp operator+ (const VComp& comp_1, const VComp& comp_2);
VComp operator+ (const VComp& comp_1, const DType& val); 
VComp operator+ (const DType& val, const VComp& comp_2); 
VComp operator+ (const VComp& comp_1); 
VComp operator- (const VComp& comp_1, const VComp& comp_2);
VComp operator- (const VComp& comp_1, const DType& val);
VComp operator- (const DType& val, const VComp& comp_2);
VComp operator- (const VComp& comp_1);

// Inner product
DType operator* (const VComp& comp_1, const VComp& comp_2);
VComp operator* (const VComp& comp_1, const DType& val);
VComp operator* (const DType& val, const VComp& comp_2);
VComp operator* (const Entry& ent, const VComp& comp);
VComp operator/ (const VComp& comp_1, const DType& val);

#endif

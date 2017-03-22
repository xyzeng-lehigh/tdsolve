#include "Alg_Solve.h"
#include <cassert>
#include <algorithm>

/*
 * Column vector storage
 */

LinVec::LinVec(void) {
    size = -1;
    dim  = -1;
    vec  = NULL;
}

LinVec::~LinVec(void) {
    if (vec!=NULL) {
        for (int i=0; i<size; i++)
            if (vec[i]!=NULL) delete vec[i];
        delete[] vec;
    }
}

void LinVec::setSize(int size_in) {
    assert(dim>0);
    size = size_in;
    if (vec!=NULL) {
        for (int i=0; i<size; i++)
            if (vec[i]!=NULL) delete vec[i];
        delete[] vec;
    }
    vec = new VComp*[size];
    for (int i=0; i<size; i++) {
        vec[i] = NULL;
        vec[i] = new VComp(dim,0.0);
    }
    return;
}

void LinVec::set(int i, VComp& comp_in) {
    *vec[i] = comp_in;
    return;
}

void LinVec::reset(int i, DType val_in) {
    *vec[i] = val_in;
    return;
}

DType LinVec::getNorm(void) {
    DType val = 0.0;
    for (int i=0; i<size; i++) {
        DType loc_norm = (*vec[i]).getNorm();
        val += loc_norm * loc_norm;
    }
    return sqrt(val);
}

void LinVec::copy(LinVec& vec_in) {
    assert( &vec_in != this );
    assert( size == vec_in.getSize() );
    assert( dim  == vec_in.getBlckDim() );
    for (int i=0; i<size; i++)
        *vec[i] = vec_in.get(i);
    return;
}

void LinVec::scale(const DType& a_in) {
    for (int i=0; i<size; i++)
        *vec[i] *= a_in;
    return;
}

void LinVec::scale_inc(const DType& a_in, const DType& b_in, LinVec& vec_in) {
    assert( &vec_in != this );
    assert( size == vec_in.getSize() );
    assert( dim  == vec_in.getBlckDim() );
    for (int i=0; i<size; i++) {
        *vec[i] *= a_in;
        *vec[i] += b_in * (vec_in.get(i));
    }
    return;
}

void LinVec::proj(LinVec& vec_proj) {
    DType inner_prod = calc_inner_prod(vec_proj,*this);
    for (int i=0; i<size; i++)
        *vec[i] -= 2.0*inner_prod*(vec_proj.get(i));
    //scale_inc(1.0,-2.0*calc_inner_prod(vec_proj,*this),vec_proj);
    return;
}

void LinVec::reset(DType val_in) {
    for (int i=0; i<size; i++)
        *vec[i] = val_in;
    return;
}

// Inner product
DType calc_inner_prod(const LinVec& vec_1, const LinVec& vec_2) {
    int size = vec_1.getSize();
    int dim  = vec_1.getBlckDim();
    assert( size == vec_2.getSize()    );
    assert( dim  == vec_2.getBlckDim() );
    DType result = 0.0;
    for (int i=0; i<size; i++)
        result += vec_1.get(i) * vec_2.get(i);
    return result;
}

/*
 * Matrix storage -- LinMat
 */

void LinMat::setSize(int size_r_in, int size_c_in) {
    assert(size_r_in>0);
    assert(size_c_in>0);
    size_r = size_r_in;
    size_c = size_c_in;
    return;
}

void LinMat::xMat_new(LinMat*& mat_in, LinMat*& mat_ot) {
    mat_ot->reset(0.0);
    int size_i = size_r;
    int size_j = mat_in->getSizeCol();
    int size_k = size_c; // == mat_in.getSizeRow();
    for (int i=0; i<size_i; i++)
        for (int j=0; j<size_j; j++) {
            Entry entry(dim,0.0);
            for (int k=0; k<size_k; k++)
                entry += this->get(i,k) * mat_in->get(k,j);
            mat_ot->set(i,j,entry);
        }
    return;
}

void LinMat::xMat_pst(LinMat*& mat_io) {
    assert( size_c == size_r );
    int size_i = mat_io->getSizeCol();
    // Create a temporary array of entries to store the current column
    Entry** cur_col = NULL;
    cur_col = new Entry*[size_r];
    for (int r=0; r<size_r; r++) {
        cur_col[r] = NULL;
        cur_col[r] = new Entry(dim,0.0);
    }

    for (int i=0; i<size_i; i++) {
        for (int r=0; r<size_r; r++)
            *cur_col[r] = mat_io->get(r,i);
        for (int r=0; r<size_r; r++) {
            Entry entry(dim,0.0);
            for (int r1=0; r1<size_r; r1++)
                entry += this->get(r,r1) * (*cur_col[r1]);
            mat_io->set(r,i,entry);
        }
    }

    if (cur_col!=NULL) {
        for (int r=0; r<size_r; r++)
            if (cur_col[r]!=NULL) delete cur_col[r];
        delete[] cur_col;
    }
    return;
}

void LinMat::xMat_pre(LinMat*& mat_in) {
    assert ( mat_in->getSizeRow() == mat_in->getSizeCol() );
    int size_i = mat_in->getSizeRow();
    // Create a temporary array of entries to store the current row
    Entry** cur_row = NULL;
    cur_row = new Entry*[size_c];
    for (int c=0; c<size_c; c++) {
        cur_row[c] = NULL;
        cur_row[c] = new Entry(dim,0.0);
    }

    for (int i=0; i<size_i; i++) {
        for (int c=0; c<size_c; c++)
            *cur_row[c] = this->get(i,c);
        for (int c=0; c<size_c; c++) {
            Entry entry(dim,0.0);
            for (int c1=0; c1<size_c; c1++)
                entry += (*cur_row[c1]) * mat_in->get(c1,c);
            this->set(i,c,entry);
        }
    }

    if (cur_row!=NULL) {
        for (int c=0; c<size_c; c++)
            if (cur_row[c]!=NULL) delete cur_row[c];
        delete[] cur_row;
    }
    return;
}

void LinMat::xVec_new(LinVec*& vec_in, LinVec*& vec_ot) {
    vec_ot->reset(0.0);
    int size_i = size_r;
    int size_j = size_c; // == vec_in.getSize();
    for (int i=0; i<size_i; i++) {
        VComp comp(dim,0.0);
        for (int j=0; j<size_j; j++)
            comp += this->get(i,j) * vec_in->get(j);
        vec_ot->set(i,comp);
    }
    return;
}

void LinMat::xVec_t_new(LinVec*& vec_in, LinVec*& vec_ot) {
    vec_ot->reset(0.0);
    int size_i = size_c;
    int size_j = size_r; // == vec_in.getSize();
    for (int i=0; i<size_i; i++) {
        VComp comp(dim,0.0);
        for (int j=0; j<size_j; j++)
            comp += (this->get(j,i)).trans() * vec_in->get(j);
        vec_ot->set(i,comp);
    }
    return;
}

void LinMat::xVec_rpl(LinVec*& vec_io) {
    assert( size_c == size_r );
    LinVec* vec_tmp = new LinVec(dim,vec_io->getSize());
    this->xVec_new(vec_io,vec_tmp);
    vec_io->copy(*vec_tmp);
    delete vec_tmp;
    return;
}

void LinMat::xVec_t_rpl(LinVec*& vec_io) {
    assert( size_c == size_r );
    LinVec* vec_tmp = new LinVec(dim,vec_io->getSize());
    this->xVec_t_new(vec_io,vec_tmp);
    vec_io->copy(*vec_tmp);
    delete vec_tmp;
    return;
}

void LinMat::xMat(LinMat*& mat_in, LinMat*& mat_ot) {
    assert( size_c == mat_in->getSizeRow() );
    assert( size_r == mat_ot->getSizeRow() );
    assert( mat_in->getSizeCol() == 
            mat_ot->getSizeCol() );
    assert( (dim==mat_in->getBlckDim()) && 
            (dim==mat_ot->getBlckDim()) );
    assert( mtype == mat_in->getType() );
    assert( mtype == mat_ot->getType() );
    if ((mat_ot!=mat_in) && (mat_ot!=this))
        xMat_new(mat_in,mat_ot);
    else if ((mat_ot!=this) && (mat_ot==mat_in))
        // In this case, this must be a square matrix
        xMat_pst(mat_in);
    else if ((mat_ot==this) && (mat_ot!=mat_in))
        // In this case, mat_in must be square
        xMat_pre(mat_in);
    else
        // Compute the square of current matrix
        xMat_sqr();
    return;
}

void LinMat::xVec(LinVec*& vec_in, LinVec*& vec_ot) {
    assert( size_c == vec_in->getSize() );
    assert( size_r == vec_ot->getSize() );
    assert( (dim==vec_in->getBlckDim()) &&
            (dim==vec_ot->getBlckDim()) );
    if (vec_in!=vec_ot)
        xVec_new(vec_in,vec_ot);
    else
        xVec_rpl(vec_in);
    return;
}

void LinMat::xVec_t(LinVec*& vec_in, LinVec*& vec_ot) {
    assert( size_r == vec_in->getSize() );
    assert( size_c == vec_ot->getSize() );
    assert( (dim==vec_in->getBlckDim()) &&
            (dim==vec_ot->getBlckDim()) );
    if (vec_in!=vec_ot)
        xVec_t_new(vec_in,vec_ot);
    else
        xVec_t_rpl(vec_in);
    return;
}

/*
 * Matrix storage -- LinSqMat
 */

void LinSqMat::setEqualSize(int size_in) {
    assert(size_in>0.0);
    size = size_in;
    LinMat::setSize(size,size);
    return;
}

/*
 * Matrix storage -- LinFullMat
 */

LinFullMat::LinFullMat(void) {
    mat     = NULL;
    mtype   = LINMAT_FULL;
    dtype   = LINMAT_GEN;
}

LinFullMat::~LinFullMat(void) {
    if (mat!=NULL) {
        for (int i=0; i<size_r*size_c; i++)
            if (mat[i]!=NULL) delete mat[i];
        delete[] mat;
    }
}

void LinFullMat::xMat_sqr(void) {
    LinMat* mat_tmp = new LinFullMat;
    LinMat* mat_ths = new LinFullMat(*this);
    mat_tmp->setBlckDim(dim);
    mat_tmp->setSize(size_r,size_c);
    xMat_new(mat_ths,mat_tmp);
    for (int i=0; i<size_r; i++)
        for (int j=0; j<size_c; j++)
            *mat[i*size_c+j] = mat_tmp->get(i,j);
    delete mat_tmp;
    delete mat_ths;
    return;
}

void LinFullMat::setSize(int size_r_in, int size_c_in) {
    LinMat::setSize(size_r_in,size_c_in);
    assert(dim>0);
    if (mat!=NULL) {
        for (int i=0; i<size_r*size_c; i++)
            if (mat[i]!=NULL) delete mat[i];
        delete[] mat;
    }
    mat = new Entry*[size_r*size_c];
    for (int i=0; i<size_r*size_c; i++) {
        mat[i] = NULL;
        mat[i] = new Entry(dim,0.0);
    }
    return;
}

Entry LinFullMat::get(int i, int j) {
    // assert((0<=i)&&(i<size_r));
    // assert((0<=j)&&(j<size_c));
    return *mat[i*size_c+j];
}

void LinFullMat::set(int i, int j, Entry& entry) {
    // assert((0<=i)&&(i<size_r));
    // assert((0<=j)&&(j<size_c));
    *mat[i*size_c+j] = entry;
    return;
}

void LinFullMat::reset(int i, int j, DType val) {
    // assert((0<=i)&&(i<size_r));
    // assert((0<=j)&&(j<size_c));
    *mat[i*size_c+j] = val;
    return;
}

void LinFullMat::reset(int i, DType val) {
    int offset = i*size_c;
    for (int j=0; j<size_c; j++)
        *mat[offset+j] = val;
    return;
}

void LinFullMat::reset(DType val) {
    for (int i=0; i<size_r*size_c; i++)
        *mat[i] = 0.0;
    if (val!=0.0) {
        int size_i = min(size_r,size_c);
        for (int i=0; i<size_i; i++)
            *mat[i*size_c+i] = val;
    }
    return;
}

/*
 * Matrix storage -- LinTriMat
 */
LinTriMat::LinTriMat(void) {
    mtype   = LINMAT_TRID;
    dtype   = LINMAT_SQ;
    lower = NULL;
    upper = NULL;
    mdiag = NULL;
}

LinTriMat::~LinTriMat(void) {
    if (lower!=NULL) {
        for (int i=0; i<size-1; i++)
            if (lower[i]!=NULL) delete lower[i];
        delete[] lower;
    }
    if (upper!=NULL) {
        for (int i=0; i<size-1; i++)
            if (upper[i]!=NULL) delete upper[i];
        delete[] upper;
    }
    if (mdiag!=NULL) {
        for (int i=0; i<size; i++)
            if (mdiag[i]!=NULL) delete mdiag[i];
        delete[] mdiag;
    }
}

void LinTriMat::xMat_new(LinMat*& mat_in, LinMat*& mat_ot) {
    mat_ot->reset(0.0);
    Entry** etemp = NULL;
    etemp = new Entry*[size];
    for (int i=0; i<size; i++) etemp[i] = NULL;
    for (int i=0; i<size; i++) etemp[i] = new Entry(dim);
    // Main diagonal
    for (int i=0; i<size; i++)
        *etemp[i] = (*mdiag[i])*mat_in->get(i,i);
    for (int i=0; i<size-1; i++)
        *etemp[i] += (*upper[i])*mat_in->get(i+1,i);
    for (int i=0; i<size-1; i++)
        *etemp[i+1] += (*lower[i])*mat_in->get(i,i+1);
    for (int i=0; i<size; i++)
        mat_ot->set(i,i,*etemp[i]);
    for (int i=0; i<size-1; i++) {
        Entry elem(dim);
        // Upper diagonal
        elem = (*mdiag[i])*mat_in->get(i,i+1)+
               (*upper[i])*mat_in->get(i+1,i+1);
        mat_ot->set(i,i+1,elem);
        // Lower diagonal
        elem = (*lower[i])*mat_in->get(i,i)+
               (*mdiag[i+1])*mat_in->get(i+1,i);
        mat_ot->set(i+1,i,elem);
    }
    if (etemp!=NULL) {
        for (int i=0; i<size; i++)
            if (etemp[i]!=NULL) delete etemp[i];
        delete[] etemp;
    }
    return;
}

void LinTriMat::xMat_pst(LinMat*& mat_io) {
    Entry elower(dim,0.0);
    Entry emdiag(dim,0.0);
    Entry eupper(dim,0.0);
    Entry m_left(dim,0.0);
    Entry m_diag(dim,0.0);
    Entry m_rght(dim,0.0);
    for (int i=0; i<size; i++) {
        Entry elem(dim);
        elower = (i<size-1) ? mat_io->get(i+1,i) : 0.0;
        emdiag = mat_io->get(i,i);
        eupper = (i>0) ? mat_io->get(i-1,i) : 0.0;
        // Upper diagonal
        if (i>0) {
            elem = (*mdiag[i-1])*eupper+(*upper[i-1])*emdiag;
            mat_io->set(i-1,i,elem);
        }
        // Main diagonal
        m_left = (i>0) ? *lower[i-1] : 0.0;
        m_diag = *mdiag[i];
        m_rght = (i<size-1) ? *upper[i] : 0.0;
        elem = m_left*eupper+m_diag*emdiag+m_rght*elower;
        mat_io->set(i,i,elem);
        // Lower diagonal
        if (i<size-1) {
            elem = (*lower[i])*emdiag+(*mdiag[i+1])*elower;
            mat_io->set(i+1,i,elem);
        }
    }
    return;
}

void LinTriMat::xMat_pre(LinMat*& mat_in) {
    Entry m_left(dim,0.0);
    Entry m_diag(dim,0.0);
    Entry m_rght(dim,0.0);
    Entry elower(dim,0.0);
    Entry emdiag(dim,0.0);
    Entry eupper(dim,0.0);
    for (int i=0; i<size; i++) {
        Entry elem(dim);
        m_left = (i>0) ? *lower[i-1] : 0.0;
        m_diag = *mdiag[i];
        m_rght = (i<size-1) ? *upper[i] : 0.0;
        // Lower diagonal
        if (i>0) {
            elem = m_left*mat_in->get(i-1,i-1)+m_diag*mat_in->get(i,i-1);
            *lower[i-1] = elem;
        }
        // Main diagonal
        elower = (i<size-1) ? mat_in->get(i+1,i) : 0.0;
        emdiag = mat_in->get(i,i);
        eupper = (i>0) ? mat_in->get(i-1,i) : 0.0;
        elem = m_left*eupper+m_diag*emdiag+m_rght*elower;
        *mdiag[i] = elem;
        // Upper diagonal
        if (i<size-1) {
            elem = m_diag*mat_in->get(i,i+1)+m_rght*mat_in->get(i+1,i+1);
            *upper[i] = elem;
        }
    }
    return;
}

void LinTriMat::xMat_sqr(void) {
    Entry** mdiag_copy = NULL;
    mdiag_copy = new Entry*[size];
    for (int i=0; i<size; i++) mdiag_copy[i] = NULL;
    for (int i=0; i<size; i++)
        mdiag_copy[i] = new Entry(*mdiag[i]);
    for (int i=0; i<size; i++)
        *mdiag[i] = (*mdiag_copy[i])*(*mdiag_copy[i]);
    for (int i=0; i<size-1; i++)
        *mdiag[i] += (*upper[i])*(*lower[i]);
    for (int i=0; i<size-1; i++)
        *mdiag[i+1] += (*lower[i])*(*upper[i]);
    for (int i=0; i<size-1; i++) {
        Entry diag_m(*mdiag_copy[i]);
        Entry diag_p(*mdiag_copy[i+1]);
        Entry elower(*lower[i]);
        Entry eupper(*upper[i]);
        *lower[i] = elower*diag_m+diag_p*elower;
        *upper[i] = diag_m*eupper+eupper*diag_p;
    }
    if (mdiag_copy!=NULL) {
        for (int i=0; i<size; i++)
            if (mdiag_copy[i]!=NULL) delete mdiag_copy[i];
        delete[] mdiag_copy;
    }
    return;
}

void LinTriMat::xVec_new(LinVec*& vec_in, LinVec*& vec_ot) {
    vec_ot->reset(0.0);
    for (int i=0; i<size; i++) {
        VComp comp(dim,0.0);
        comp = (*mdiag[i])*vec_in->get(i);
        vec_ot->set(i,comp);
    }
    for (int i=0; i<size-1; i++) {
        VComp comp(dim,0.0);
        comp  = vec_ot->get(i);
        comp += (*upper[i])*vec_in->get(i+1);
        vec_ot->set(i,comp);
        comp  = vec_ot->get(i+1);
        comp += (*lower[i])*vec_in->get(i);
        vec_ot->set(i+1,comp);
    }
    return;
}

void LinTriMat::xVec_t_new(LinVec*& vec_in, LinVec*& vec_ot) {
    vec_ot->reset(0.0);
    for (int i=0; i<size; i++) {
        VComp comp(dim,0.0);
        comp = (*mdiag[i]).trans()*vec_in->get(i);
        vec_ot->set(i,comp);
    }
    for (int i=0; i<size-1; i++) {
        VComp comp(dim,0.0);
        comp  = vec_ot->get(i);
        comp += (*lower[i]).trans()*vec_in->get(i+1);
        vec_ot->set(i,comp);
        comp  = vec_ot->get(i+1);
        comp += (*upper[i]).trans()*vec_in->get(i);
        vec_ot->set(i+1,comp);
    }
    return;
}

void LinTriMat::setEqualSize(int size_in) {
    LinSqMat::setEqualSize(size_in);
    assert(dim>0);
    if (lower!=NULL) {
        for (int i=0; i<size-1; i++)
            if (lower[i]!=NULL) delete lower[i];
        delete[] lower;
    }
    if (upper!=NULL) {
        for (int i=0; i<size-1; i++)
            if (upper[i]!=NULL) delete upper[i];
        delete[] upper;
    }
    if (mdiag!=NULL) {
        for (int i=0; i<size; i++)
            if (mdiag[i]!=NULL) delete mdiag[i];
        delete[] mdiag;
    }
    if (size>1) lower = new Entry*[size-1];
    if (size>1) upper = new Entry*[size-1];
    mdiag = new Entry*[size];
    for (int i=0; i<size-1; i++) {
        lower[i] = NULL;
        lower[i] = new Entry(dim,0.0);
        upper[i] = NULL;
        upper[i] = new Entry(dim,0.0);
    }
    for (int i=0; i<size; i++) {
        mdiag[i] = NULL;
        mdiag[i] = new Entry(dim,0.0);
    }
    return;
}

Entry LinTriMat::get(int i, int j) {
    // assert((0<=i)&&(i<size));
    // assert((0<=j)&&(j<size));
    if (i==j)
        return *mdiag[i];
    else if (i==(j+1))
        return *lower[j];
    else if ((i+1)==j)
        return *upper[i];
    else
        return Entry(dim,0.0);
}

void LinTriMat::set(int i, int j, Entry& entry) {
    // assert((0<=i)&&(i<size));
    // assert((0<=j)&&(j<size));
    if (i==j)
        *mdiag[i] = entry;
    else if (i==(j+1))
        *lower[j] = entry;
    else if ((i+1)==j)
        *upper[i] = entry;
    else
        fprintf(stderr,"Warning: trying to assign entry to a zero element of a tri-diagonal matrix, will do nothing");
    return;
}

void LinTriMat::reset(int i, int j, DType val) {
    // assert((0<=i)&&(i<size));
    // assert((0<=j)&&(j<size));
    if (i==j)
        *mdiag[i] = val;
    else if (i==(j+1))
        *lower[j] = val;
    else if ((i+1)==j)
        *upper[i] = val;
    else
        fprintf(stderr,"Warning: trying to reset a zero element of a tri-diagonal matrix, will do nothing");
    return;
}

void LinTriMat::reset(int i, DType val) {
    *mdiag[i] = val;
    if (i>0)
        *lower[i-1] = val;
    if (i<size-1)
        *upper[i+1] = val;
    return;
}

void LinTriMat::reset(DType val) {
    for (int i=0; i<size; i++)
        *mdiag[i] = val;
    for (int i=0; i<size-1; i++) {
        *lower[i] = 0.0;
        *upper[i] = 0.0;
    }
    return;
}

/*
 * Matrix storage -- LinSparseMat
 */

LinSparseMat::LinSparseMat(void) {
    mtype   = LINMAT_SPARSE;
    dtype   = LINMAT_GEN;
    mat.clear();
}

void LinSparseMat::xMat_new(LinMat*& mat_in, LinMat*& mat_ot) {
    mat_ot->reset(0.0);
    int size_k = mat_ot->getSizeCol();
    for (std::map<sm_idx,Entry>::iterator it=mat.begin();
         it!=mat.end(); it++) {
        sm_idx idx  = (*it).first;
        Entry entry = (*it).second;
        for (int k=0; k<size_k; k++) {
            Entry ent_in = mat_in->get(idx.j,k);
            if (!ent_in.sig) {
                Entry ent_ot = mat_ot->get(idx.i,k);
                ent_ot += entry * ent_in;
                mat_ot->set(idx.i,k,ent_ot);
            }
        }
    }
    return;
}

void LinSparseMat::xMat_pst(LinMat*& mat_io) {
    assert( size_c == size_r );
    int size_i = mat_io->getSizeCol();
    // A temporary map of a column
    std::map<sm_idx,Entry> cur_col;
    for (int i=0; i<size_i; i++)
        if ((mat_io->getCol(cur_col,i))>0) {
            for (int r=0; r<size_i; r++) {
                sm_idx idx(r,i);
                Entry  entry(dim,0.0);
                bool   is_zero = true;
                for (std::map<sm_idx,Entry>::iterator it=cur_col.begin();
                     it!=cur_col.end(); it++) {
                    Entry ent_this = this->get(r,(*it).first.i);
                    if (!ent_this.sig) {
                        is_zero = false;
                        entry += ent_this * ((*it).second);
                    }
                }
                if (!is_zero)
                    mat_io->set(r,i,entry);
            }
        }
    return;
}

void LinSparseMat::xMat_pre(LinMat*& mat_in) {
    assert ( mat_in->getSizeRow() == mat_in->getSizeCol() );
    int size_i = mat_in->getSizeRow();
    // A temporary map of a row
    std::map<sm_idx,Entry> cur_row;
    for (int i=0; i<size_i; i++)
        if ((this->getRow(cur_row,i))>0) {
            for (int c=0; c<size_i; c++) {
                sm_idx idx(i,c);
                Entry entry(dim,0.0);
                bool is_zero = true;
                for (std::map<sm_idx,Entry>::iterator it=cur_row.begin();
                     it!=cur_row.end(); it++) {
                    Entry ent_in = mat_in->get(c,(*it).first.j);
                    if (!ent_in.sig) {
                        is_zero = false;
                        entry += ((*it).second) * ent_in;
                    }
                }
                if (!is_zero)
                    this->set(i,c,entry);
            }
        }
    return;
}

void LinSparseMat::xMat_sqr(void) {
    LinMat* mat_tmp = new LinSparseMat;
    LinMat* mat_ths = new LinSparseMat(*this);
    mat_tmp->setBlckDim(dim);
    mat_tmp->setSize(size_r,size_c);
    xMat_new(mat_ths,mat_tmp);
    mat.clear();
    for (std::map<sm_idx,Entry>::iterator it=mat_tmp->it_beg();
         it!=mat_tmp->it_end(); it++)
        mat.insert(pair<sm_idx,Entry>((*it).first,(*it).second));
    delete mat_tmp;
    delete mat_ths;
    return;
}

void LinSparseMat::xVec_new(LinVec*& vec_in, LinVec*& vec_ot) {
    vec_ot->reset(0.0);
    for (std::map<sm_idx,Entry>::iterator it=mat.begin();
         it!=mat.end(); it++) {
        sm_idx idx  = (*it).first;
        Entry entry = (*it).second;
        VComp comp = vec_ot->get(idx.i);
        comp += entry * vec_in->get(idx.j);
        vec_ot->set(idx.i,comp);
    }
    return;
}

void LinSparseMat::xVec_t_new(LinVec*& vec_in, LinVec*& vec_ot) {
    vec_ot->reset(0.0);
    for (std::map<sm_idx,Entry>::iterator it=mat.begin();
         it!=mat.end(); it++) {
        sm_idx idx  = (*it).first;
        Entry entry = (*it).second; // mat[idx.i][idx.j]
        VComp comp = vec_ot->get(idx.j);
        comp += entry.trans() * vec_in->get(idx.i);
        vec_ot->set(idx.j,comp);
    }
    return;
}

void LinSparseMat::setSize(int size_r_in, int size_c_in) {
    LinMat::setSize(size_r_in,size_c_in);
    mat.clear();
    return;
}

Entry LinSparseMat::get(int i, int j) {
    sm_idx idx(i,j);
    std::map<sm_idx,Entry>::iterator iter;
    iter = mat.find(idx);
    if (iter==mat.end()) {
        Entry entry(dim,0.0);
        entry.sig = true;
        return entry;
    }
    else
        return (*iter).second;
}

void LinSparseMat::set(int i, int j, Entry& entry) {
    sm_idx idx(i,j);
    std::map<sm_idx,Entry>::iterator iter;
    iter = mat.find(idx);
    Entry loc_entry = entry;
    loc_entry.sig = false;
    if (iter!=mat.end())
        mat.erase(iter);
    mat.insert(pair<sm_idx,Entry>(idx,loc_entry));
    return;
}

void LinSparseMat::reset(int i, int j, DType val) {
    Entry entry(dim,val);
    if (val==0.0)
        mat.erase(sm_idx(i,j));
    else {
        entry.sig = false;
        this->set(i,j,entry);
    }
    return;
}

void LinSparseMat::reset(int i, DType val) {
    // Reset all registered entries in the i-th row to val
    sm_idx row_beg(i,  0);
    sm_idx row_end(i+1,0);
    for (std::map<sm_idx,Entry>::iterator it=mat.lower_bound(row_beg);
         (*it).first < row_end; it++)
        (*it).second = val;
    return;
}

void LinSparseMat::reset(DType val) {
    mat.clear();
    int size_i = min(size_c,size_r);
    if (val!=0.0)
        for (int i=0; i<size_i; i++) {
            sm_idx idx(i,i);
            Entry entry(dim,val);
            entry.sig = false;
            mat.insert(pair<sm_idx,Entry>(idx,entry));
        }
    return;
}

int LinSparseMat::getCol(std::map<sm_idx,Entry>& col, int col_idx) {
    col.clear();
    // Since sm_idx is sorted first according to row
    // number, a slower approach is adopted here
    for (std::map<sm_idx,Entry>::iterator it=mat.begin();
         it!=mat.end(); it++) {
        sm_idx idx = (*it).first;
        if (idx.j==col_idx)
            col.insert(pair<sm_idx,Entry>(idx,(*it).second));
    }
    return col.size();
}

int LinSparseMat::getRow(std::map<sm_idx,Entry>& row, int row_idx) {
    row.clear();
    
    // // The following is a slower approach
    // for (std::map<sm_idx,Entry>::iterator it=mat.begin();
    //      it!=mat.end(); it++) {
    //     sm_idx idx = (*it).first;
    //     if (idx.i==row_idx)
    //         row.insert(pair<sm_idx,Entry>(idx,(*it).second));
    // }

    // Using the specific ordering of sm_idx, the
    // following is a faster approach
    sm_idx row_beg(row_idx,  0);
    sm_idx row_end(row_idx+1,0);
    for (std::map<sm_idx,Entry>::iterator it=mat.lower_bound(row_beg);
         (*it).first < row_end; it++)
        row.insert(pair<sm_idx,Entry>((*it).first,(*it).second));
    return row.size();
}

/*
 * Linear solve -- LinSolve
 */

LinSolve::LinSolve(void) {
    size = -1;
    dim  = -1;
    rhs  = NULL;
    sol  = NULL;
    mat  = NULL;
    prc  = NULL;
    res  = NULL;
    preCond  = false;
    max_iter = -1;
    prec_ilu_fill = 0;
    coef_abs = 1.0;
    coef_ini = 0.0;
    coef_rhs = 0.0;
    print_stat  = false;
    target_res  = TOL_ZERO;
    rhs_nrm_ini = 1.0;
    abs_res.clear();
    rel_res_ini.clear();
    rel_res_rhs.clear();
    res_hist.clear();

    lsysType   = LSYS_TYPE_KRYLOV;
    lsysIter   = LSYS_ITER_GS;
    lsysKrylov = LSYS_KRYLOV_GMRES;
    lsysPrec   = LSYS_PREC_OFF;

    apply_prec          = NULL;
    solve_iter_single   = NULL;
    solve_krylov_pre    = NULL;
    solve_krylov_single = NULL;

    ilu_nzpat     = NULL;
    ilu_idx       = NULL;

    sor_par       = 1.0; // So SOR is effectively GS
    kspace_dim    = 0;   // Full GMRES
    kspace_res    = true;
    kspace_signal = false;

    vec_v.clear();
    vec_w.clear();
    vec_p.clear();

    gmres_g = NULL;
    gmres_z = NULL;
    gmres_v = NULL;
    gmres_h = NULL;
    gmres_cos.clear();
    gmres_sin.clear();

    lanczos_alpha.clear();
    lanczos_delta.clear();
    lanczos_beta.clear();
    lanczos_sol0 = NULL;

    vec_cg_r = NULL;
    vec_cg_p = NULL;

}

LinSolve::~LinSolve(void) {
    clear_krylov_spaces(vec_v);
    clear_krylov_spaces(vec_w);
    clear_krylov_spaces(vec_p);
    if (gmres_g!=NULL)      delete gmres_g;
    if (gmres_z!=NULL)      delete gmres_z;
    if (gmres_v!=NULL)      delete gmres_v;
    if (gmres_h!=NULL)      delete gmres_h;
    if (lanczos_sol0!=NULL) delete lanczos_sol0;
    if (vec_cg_r!=NULL)     delete vec_cg_r;
    if (vec_cg_p!=NULL)     delete vec_cg_p;
    if (ilu_nzpat!=NULL)    delete[] ilu_nzpat;
    if (ilu_idx!=NULL)      delete[] ilu_idx;
}

DType LinSolve::compute_res(void) {
    // res = 0.0
    res->reset(0.0); 
    // res = mat * sol
    mat->xVec(sol,res);         
    // res = - mat * sol + rhs
    res->scale_inc(-1.0,1.0,*rhs);
    return res->getNorm();
}

DType LinSolve::compute_tot_res(DType res_in) {
    DType cur_abs_res = res_in;
    DType cur_rel_ini = (fabs(rel_res_ini[0])>0.0) ? 
        res_in/fabs(rel_res_ini[0]) : 0.0;
    DType cur_rel_rhs = (fabs(rhs_nrm_ini)>0.0) ?
        res_in/fabs(rhs_nrm_ini) : 0.0;

    DType cur_tot_res = coef_abs*cur_abs_res +
                        coef_ini*cur_rel_ini +
                        coef_rhs*cur_rel_rhs;
    return cur_tot_res;
}

bool LinSolve::update_res(DType res_in) {
    DType cur_res = (res_in < 0.0) ? res->getNorm()
                                   : res_in;

    DType cur_abs_res = cur_res;
    DType cur_rel_ini = (fabs(rel_res_ini[0])>0.0) ? 
        cur_res/fabs(rel_res_ini[0]) : 0.0;
    DType cur_rel_rhs = (fabs(rhs_nrm_ini)>0.0) ?
        cur_res/fabs(rhs_nrm_ini) : 0.0;

    abs_res.push_back(cur_abs_res);
    rel_res_ini.push_back(cur_rel_ini);
    rel_res_rhs.push_back(cur_rel_rhs);

    DType cur_tot_res = coef_abs*cur_abs_res +
                        coef_ini*cur_rel_ini +
                        coef_rhs*cur_rel_rhs;

    res_hist.push_back(cur_tot_res);

    return (cur_tot_res>=target_res);

}

bool LinSolve::check_max_iter(int cur_iter) {
    return (max_iter<0) || (abs(cur_iter)<=abs(max_iter));
}

void LinSolve::apply_prec_link(int dummy) {
    prc->xMat(mat,mat);
    prc->xVec(rhs,rhs);
    return;
}

void LinSolve::apply_prec_eye(int semi_dummy) {
    prc->reset(1.0*((double)semi_dummy+1.0));
    prc->xMat(mat,mat);
    prc->xVec(rhs,rhs);
    return;
}

void LinSolve::apply_prec_ilu_0(int dummy_level) {
    DType l_tol = 1.0e-6;
    // Compute the incomplete LU preconditioner first
    calc_nz_pattern(mat,ilu_nzpat,ilu_idx,l_tol);
    
    // Construct P=prc=(LU)^{-1}, L is unit lower, U is invertible upper
    // LU is the incomplete LU factorization of A=mat

    UTriMAT matU(dim*size,0.0);
    for (int i=0; i<dim*size; i++) {
        DType dval = mat->getGlobVal(i,i);
        matU.get(i,i) = (fabs(dval)<=l_tol) ? 1.0 : dval;
        for (int j=i+1; j<dim*size; j++)
            matU.get(i,j) = mat->getGlobVal(i,j);
    }
    LTriMAT matL(dim*size,0.0);
    for (int i=0; i<dim*size; i++)
        matL.get(i,i) = 1.0;
    for (int i=1; i<dim*size; i++)
        for (int j=0; j<i; j++)
            matL.get(i,j) = mat->getGlobVal(i,j);

    for (int i=1; i<dim*size; i++) {
        std::vector<int>& pat = ilu_nzpat[i];
        for (int idx_k=0; idx_k<ilu_idx[i]; idx_k++) {
            int k = pat[idx_k];    // k < i
            matL.get(i,k) /= matU.get(k,k);
            for (unsigned idx_j=idx_k+1; idx_j<pat.size(); idx_j++) {
                int j = pat[idx_j];
                DType& val = (i>j) ? matL.get(i,j) : matU.get(i,j);
                val -= matL.get(i,k) * matU.get(k,j);
            }
        }
    }

    HMAT p_mat(dim*size,0.0);
    p_mat = matU.inv() * matL.inv();
    convert(prc,p_mat,dim,size,l_tol);

    prc->xMat(mat,mat);
    prc->xVec(rhs,rhs);
    return;
}

void LinSolve::calc_nz_pattern(LinMat*& mat_in, std::vector<int>*& pat_in, int*& idx_in, DType tol_in) {
    if (pat_in!=NULL) delete[] pat_in;
    pat_in = new std::vector<int>[dim*size];
    for (int i=0; i<dim*size; i++) pat_in[i].clear();
    if (idx_in!=NULL) delete[] idx_in;
    idx_in = new int[dim*size];
    for (int i=0; i<dim*size; i++) idx_in[i] = 0;

    tol_in = max(fabs(tol_in),TOL_ZERO);

    if (mat_in->isSparse()) {
        for (std::map<sm_idx,Entry>::iterator it=mat_in->it_beg(); it!=mat_in->it_end(); it++) {
            sm_idx l_idx = (*it).first;
            Entry  entry = (*it).second;
            int I_offset = l_idx.i * dim;
            int J_offset = l_idx.j * dim;
            for (int i=0; i<dim; i++) 
                for (int j=0; j<dim; j++) {
                    int I = I_offset+i;
                    int J = J_offset+j;
                    if (fabs(entry.get(i,j))>tol_in) {
                        pat_in[I].push_back(J);
                        if (I>J) idx_in[J]++;
                    }
                }
        }
        for (int i=0; i<dim*size; i++) {
            std::vector<int>& pat = pat_in[i];
            if (std::find(pat.begin(),pat.end(),i)==pat.end()) pat.push_back(i);
            std::sort(pat.begin(),pat.end());
        }
    }
    if (mat_in->isFull()) {
        for (int i=0; i<dim*size; i++)
            for (int j=0; j<dim*size; j++) {
                if (j==i)
                    pat_in[i].push_back(j);
                else {
                    if (fabs(mat_in->getGlobVal(i,j))>tol_in) {
                        pat_in[i].push_back(j);
                        if (i>j) idx_in[j]++;
                    }
                }
            }
    }
    if (mat_in->isTriD()) {
        fprintf(stderr,"LinSolve::calc_nz_pattern -- tridiagonal solver not implemented ... Exiting ...\n");
        exit(-1);
    }
    return;
}

void LinSolve::apply_prec_ilu_p(int fill_level) {
    // Compute the incomplete LU preconditioner first
    fprintf(stderr,"Error: LinSolve ILU preconditioner with nonzero fill-in levels not implememnted ... Exiting ...\n");
    exit(-1);
    return;
}

void LinSolve::apply_prec_bilu_0(int dummy_level) {
    DType l_tol = 1.0e-6;
    // Compute the incomplete LU preconditioner first
    calc_nz_bpattern(mat,ilu_nzpat,ilu_idx,l_tol);

    // Construct P=prc=(LU)^{-1}, L is unit lower, U is invertible upper
    // LU is the incomplete LU factorization of A=mat

    LinSqFullMat matLU;
    matLU.setBlckDim(dim);
    matLU.setEqualSize(size);
    
    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++) {
            Entry elem = mat->get(i,j);
            matLU.set(i,j,elem);
        }
    for (int i=0; i<size; i++) {
        Entry elem = mat->get(i,i);
        if (fabs((elem.getVal()).det())<=l_tol) {
            elem = 1.0;
            matLU.set(i,i,elem);
        }
    }

    // Construct LU and store it to matLU
    for (int i=1; i<size; i++) {
        std::vector<int>& pat = ilu_nzpat[i];
        for (int idx_k=0; idx_k<ilu_idx[i]; idx_k++) {
            int k = pat[idx_k]; // k < i
            Entry elem = matLU.get(i,k);
            HMAT  mate = elem.getVal();
            HMAT  matd = (matLU.get(k,k)).getVal();
            elem = mate*matd.inv();
            matLU.set(i,k,elem);
            for (unsigned idx_j=idx_k+1; idx_j<pat.size(); idx_j++) {
                int j = pat[idx_j];
                Entry l_elem = matLU.get(i,j);
                HMAT  l_emat = l_elem.getVal();
                HMAT  l_mat1 = (matLU.get(i,k)).getVal();
                HMAT  l_mat2 = (matLU.get(k,j)).getVal();
                l_elem = l_emat - l_mat1*l_mat2;
                matLU.set(i,j,l_elem);
            }
        }
    }
    
    // Construct P = U^{-1}L^{-1}
    // Step 1: replace the upper triangular part (including the diagonal) of matLU with U^{-1}
    for (int i=0; i<size; i++) {
        Entry elem = matLU.get(i,i);
        elem = (elem.getVal()).inv();
        matLU.set(i,i,elem);
    }
    for (int i=0; i<size-1; i++) {
        for (int j=i+1; j<size; j++) {
            HMAT  emat(dim,0.0);
            for (int k=i; k<j; k++)
                emat -= (matLU.get(i,k)).getVal() *
                        (matLU.get(k,j)).getVal();
            emat = emat*((matLU.get(j,j)).getVal());
            Entry elem = matLU.get(i,j); // Value is redundant, size matters
            elem = emat;
            matLU.set(i,j,elem);
        }
    }
    // Step 2: replace the lower triangular part (excluding the diagonal) of matLU with L^{-1}
    for (int i=1; i<size; i++) {
        for (int j=0; j<i; j++) {
            HMAT emat = - (matLU.get(i,j)).getVal();
            for (int k=j+1; k<i; k++)
                emat -= (matLU.get(i,k)).getVal() *
                        (matLU.get(k,j)).getVal();
            Entry elem = matLU.get(i,j); // Value is redundant, size matters
            elem = emat;
            matLU.set(i,j,elem);
        }
    }
    // Step 4: put U^{-1}L^{-1} into prc
    bool spse = prc->isSparse();
    prc->reset();
    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++) {
            HMAT emat(dim,0.0);
            for (int k=max(i,j); k<size; k++) {
                if (k==j) 
                    emat += (matLU.get(i,k)).getVal();
                else
                    emat += (matLU.get(i,k)).getVal() *
                            (matLU.get(k,j)).getVal();
            }
            Entry elem(dim);
            elem = emat;
            if (spse) {
                if (elem.getNorm()>l_tol)
                    prc->set(i,j,elem);
            }
            else
                prc->set(i,j,elem);
        }

    prc->xMat(mat,mat);
    prc->xVec(rhs,rhs);
    return;
}

void LinSolve::calc_nz_bpattern(LinMat*& mat_in, std::vector<int>*& pat_in, int*& idx_in, DType tol_in) {
    if (pat_in!=NULL) delete[] pat_in;
    pat_in = new std::vector<int>[size];
    for (int i=0; i<size; i++) pat_in[i].clear();
    if (idx_in!=NULL) delete[] idx_in;
    idx_in = new int[size];
    for (int i=0; i<size; i++) idx_in[i] = 0;

    tol_in = max(fabs(tol_in),TOL_ZERO);

    if (mat_in->isSparse()) {
        // Search for diagonal entries of mat_in
        // if it is not tracked by the map 
        // put a zero entry here
        // !! This may change the input matrix mat
        for (int i=0; i<size; i++) {
            Entry entry = mat_in->get(i,i);
            if (entry.sig) {
                Entry elem(dim,0.0);
                mat_in->set(i,i,elem);
            }
        }
        for (std::map<sm_idx,Entry>::iterator it=mat_in->it_beg(); it!=mat_in->it_end(); it++) {
            sm_idx l_idx = (*it).first;
            Entry  entry = (*it).second;
            if (!entry.sig) {
                pat_in[l_idx.i].push_back(l_idx.j);
                if (l_idx.i>l_idx.j) idx_in[l_idx.j]++;
            }
        }
        for (int i=0; i<size; i++) {
            std::vector<int>& pat = pat_in[i];
            std::sort(pat.begin(),pat.end());
        }
    }
    if (mat_in->isFull()) {
        for (int i=0; i<size; i++)
            for (int j=0; j<size; j++) {
                if (j==i)
                    pat_in[i].push_back(j);
                else {
                    Entry elem = mat_in->get(i,j);
                    if (elem.getNorm()>tol_in) {
                        pat_in[i].push_back(j);
                        if (i>j) idx_in[j]++;
                    }
                }
            }
    }
    if (mat_in->isTriD()) {
        fprintf(stderr,"LinSolve::calc_nz_bpattern -- tridiagonal solver not implemented ... Exiting ...\n");
        exit(-1);
    }
    return;
}

void LinSolve::apply_prec_bilu_p(int dummy_level) {
    // Compute the incomplete LU preconditioner first
    fprintf(stderr,"Error: LinSolve Block-ILU preconditioner with nonzero fill-in levels not implememnted ... Exiting ...\n");
    exit(-1);
}

bool LinSolve::solve_start(void) {
    abs_res.clear();
    rel_res_ini.clear();
    rel_res_rhs.clear();
    res_hist.clear();
    
    if (lsysType==LSYS_TYPE_KRYLOV) {
        if (max_iter<0) max_iter = dim*size-2;
        max_iter = min(max_iter,dim*size-2);
    }

    DType cur_abs_res = compute_res();
    abs_res.push_back(cur_abs_res);
    DType cur_rel_ini = (fabs(cur_abs_res)>0.0) ?
        1.0 : 0.0;
    rel_res_ini.push_back(cur_rel_ini);
    rhs_nrm_ini = rhs->getNorm();
    DType cur_rel_rhs = (fabs(rhs_nrm_ini)>0.0) ?
        cur_abs_res/rhs_nrm_ini : 0.0;
    rel_res_rhs.push_back(cur_rel_rhs);

    DType cur_tot_res = coef_abs*cur_abs_res +
                        coef_ini*cur_rel_ini +
                        coef_rhs*cur_rel_rhs;
    res_hist.push_back(cur_tot_res);

    int init_iter_num = 0;

    return check_max_iter(init_iter_num) && 
           (cur_tot_res >= target_res);
}

void LinSolve::solve_dir(bool toContinue, int cur_iter) {
    if (toContinue) {
        HVEC vec_rhs = convert(rhs);
        HMAT mat_lhs = convert(mat);
        HVEC vec_sol = mat_lhs.inv()*vec_rhs;
        convert(sol,vec_sol,dim,size);
        DType cur_res = compute_res();
        update_res(cur_res);
        cur_iter++;
        if (print_stat)
            fprintf(stdout,"Iter %d: Current Residual %10.20f\n",cur_iter,res_hist.back());
        toContinue = false;
    }
    return;
}

void LinSolve::solve_iter(bool toContinue, int cur_iter) {
    while (toContinue) {
        (this->*solve_iter_single)(sor_par);
        DType cur_res = compute_res();
        toContinue = update_res(cur_res);
        cur_iter ++;
        toContinue = toContinue && check_max_iter(cur_iter);
        if (print_stat)
            fprintf(stdout,"Iter %d: Current Residual %10.20f\n",cur_iter,res_hist.back());
    }
    return;
}

void LinSolve::solve_iter_jac(DType dummy) {
    // A = mat
    // x = sol
    // b = rhs
    // r = res = b - A x
    // A = D + L + U
    // x <- D.inv() [ - (L + U) x + b ]
    //    = D.inv() [ D x - A x + b ]
    //    = x + D.inv() r
    Entry entry(dim);
    for (int i=0; i<size; i++) {
        VComp comp(dim);
        comp = sol->get(i);
        entry = mat->get(i,i);
        comp += entry.getInv() * res->get(i);
        sol->set(i,comp);
    }
    return;
}

void LinSolve::solve_iter_gs(DType dummy) {
    // A = mat
    // x = sol
    // b = rhs
    // r = res = b - A x
    // A = D + L + U
    // x <- (D+L).inv() [ - U x + b ]
    //    = (D+L).inv() [ (D+L) x - A x + b]
    //    = x + (D+L).inv() r
    // D+L is lower triangular
    Entry entry(dim);
    LinVec vec_tmp(dim,size);
    for (int i=0; i<size; i++) {
        VComp sol_comp(dim);
        VComp res_comp(dim);
        res_comp = res->get(i);
        entry = mat->get(i,i);
        // assert(entry.sig==false);
        for (int j=0; j<i; j++) {
            Entry l_entry(dim);
            l_entry = mat->get(i,j);
            if (!l_entry.sig)
                res_comp -= l_entry * vec_tmp.get(j);
        }
        sol_comp = entry.getInv()*res_comp;
        vec_tmp.set(i,sol_comp);
        sol_comp += sol->get(i);
        sol->set(i,sol_comp);
    }
    return;
}

void LinSolve::solve_iter_bgs(DType dummy) {
    // A = mat
    // x = sol
    // b = rhs
    // r = res = b - A x
    // A = D + L + U
    // x <- (D+U).inv() [ - L x + b ]
    //    = (D+U).inv() [ (D+U) x - A x + b]
    //    = x + (D+U).inv() r
    // D+U is upper triangular
    Entry entry(dim);
    LinVec vec_tmp(dim,size);
    for (int i=size-1; i>-1; i--) {
        VComp sol_comp(dim);
        VComp res_comp(dim);
        res_comp = res->get(i);
        entry = mat->get(i,i);
        // assert(entry.sig==false);
        for (int j=i+1; j<size; j++) {
            Entry l_entry(dim);
            l_entry = mat->get(i,j);
            if (!l_entry.sig)
                res_comp -= l_entry * vec_tmp.get(j);
        }
        sol_comp = entry.getInv()*res_comp;
        vec_tmp.set(i,sol_comp);
        sol_comp += sol->get(i);
        sol->set(i,sol_comp);
    }
    return;
}

void LinSolve::solve_iter_sor(DType omega) {
    // A = mat
    // x = sol
    // b = rhs
    // r = res = b - A x
    // A = D + L + U
    // w = omega
    // x <- (D+wL).inv() [ (-wU + (1-w)D) x + wb ]
    //    = (D+wL).inv() [ (D+wL) x + w(-Ax + b) ]
    //    = x + w (D+wL).inv() r
    // D+wL is lower triangular
    Entry entry(dim);
    LinVec vec_tmp(dim,size);
    for (int i=0; i<size; i++) {
        VComp sol_comp(dim);
        VComp res_comp(dim);
        res_comp = res->get(i);
        entry = mat->get(i,i);
        // assert(entry.sig==false);
        for (int j=0; j<i; j++) {
            Entry l_entry(dim);
            l_entry = mat->get(i,j);
            if (!l_entry.sig)
                res_comp -= omega * l_entry * vec_tmp.get(j);
        }
        sol_comp = entry.getInv()*res_comp;
        vec_tmp.set(i,sol_comp);
        sol_comp *= omega;
        sol_comp += sol->get(i);
        sol->set(i,sol_comp);
    }
    return;
}

void LinSolve::solve_iter_bsor(DType omega) {
    // A = mat
    // x = sol
    // b = rhs
    // r = res = b - A x
    // A = D + L + U
    // w = omega
    // x <- (D+wU).inv() [ (-wL + (1-w)D)  x + wb ]
    //    = (D+wU).inv() [ (D+wU) x + w(-A x + b) ]
    //    = x + w (D+wU).inv() r
    // D+wU is upper triangular
    Entry entry(dim);
    LinVec vec_tmp(dim,size);
    for (int i=size-1; i>-1; i--) {
        VComp sol_comp(dim);
        VComp res_comp(dim);
        res_comp = res->get(i);
        entry = mat->get(i,i);
        // assert(entry.sig==false);
        for (int j=i+1; j<size; j++) {
            Entry l_entry(dim);
            l_entry = mat->get(i,j);
            if (!l_entry.sig)
                res_comp -= omega * l_entry * vec_tmp.get(j);
        }
        sol_comp = entry.getInv()*res_comp;
        vec_tmp.set(i,sol_comp);
        sol_comp *= omega;
        sol_comp += sol->get(i);
        sol->set(i,sol_comp);
    }
    return;
}

void LinSolve::solve_krylov(bool toContinue, int cur_iter) {
    // Temporary storage for vectors that span the
    // Krylov space
    if (toContinue) (this->*solve_krylov_pre)();
    while (toContinue) {
        (this->*solve_krylov_single)();
        DType cur_res = (kspace_res) ? res->getNorm() : compute_res();
        toContinue = update_res(cur_res);
        cur_iter ++;
        toContinue = toContinue && check_max_iter(cur_iter);
        if (kspace_signal)  toContinue = false;
        if (print_stat)
            fprintf(stdout,"Iter %d: Current Residual %10.20f\n",cur_iter,res_hist.back());
    }
    return;
}

void LinSolve::clear_krylov_spaces(std::vector<LinVec*>& vecs) {
    for (unsigned i=0; i<vecs.size(); i++) delete vecs[i];
    vecs.clear();
    return;
}

void LinSolve::solve_krylov_gmres_pre(void) {
    clear_krylov_spaces(vec_w);
    gmres_cos.clear();
    gmres_sin.clear();
    if (gmres_g!=NULL) delete gmres_g;
    gmres_g = new LinVec(dim,size);
    if (gmres_z!=NULL) delete gmres_z;
    gmres_z = new LinVec(dim,size);
    if (gmres_v!=NULL) delete gmres_v;
    gmres_v = new LinVec(dim,size);
    if (gmres_sol0!=NULL) delete gmres_sol0;
    gmres_sol0 = new LinVec(dim,size);
    if (gmres_h!=NULL) delete gmres_h;
    gmres_h = NULL; // Be cautious

    gmres_z->copy(*res);

    // w_1 or P_1
    LinVec* loc_w = new LinVec(dim,size);
    loc_w->copy(*gmres_z);
    DType res0_nrm = loc_w->getNorm();
    if (loc_w->getGlob(0)>=0.0)
        loc_w->getGlob(0) += res0_nrm;
    else
        loc_w->getGlob(0) -= res0_nrm;
    res0_nrm = loc_w->getNorm();
    loc_w->scale(1.0/res0_nrm);
    vec_w.push_back(loc_w);

    // h_0
    LinVec loc_h(dim,size);
    loc_h.copy(*gmres_z);
    loc_h.proj(*loc_w);

    // g
    gmres_g->getGlob(0) = loc_h.getGlob(0);

    // v
    gmres_v->getGlob(0) = 1.0;
    gmres_v->proj(*loc_w);

    // z
    mat->xVec(gmres_v,gmres_z);
    gmres_z->proj(*loc_w);

    // x0
    gmres_sol0->copy(*sol);
    return;
}

void LinSolve::solve_krylov_gmres(void) {
    int dim_v;
    if (gmres_h==NULL) {
        dim_v = 1;
        gmres_h = new UTriMAT(dim_v,0.0);
    }
    else {
        gmres_h->expand();
        dim_v = gmres_h->size();
    }

    // w_j
    LinVec* loc_w = new LinVec(dim,size);
    loc_w->copy(*gmres_z);
    for (int i=0; i<dim_v; i++)
        loc_w->getGlob(i) = 0.0;
    DType w_nrm = loc_w->getNorm();
    if (loc_w->getGlob(dim_v)>=0.0)
        loc_w->getGlob(dim_v) += w_nrm;
    else
        loc_w->getGlob(dim_v) -= w_nrm;
    w_nrm = loc_w->getNorm();
    loc_w->scale(1.0/w_nrm);
    vec_w.push_back(loc_w);

    // h_{j-1}
    LinVec loc_h(dim,size);
    loc_h.copy(*gmres_z);
    loc_h.proj(*loc_w);
    for (unsigned i=0; i<gmres_cos.size(); i++) {
        DType loc_cos = gmres_cos[i];
        DType loc_sin = gmres_sin[i];
        DType val_m = loc_h.getGlob(i);
        DType val_p = loc_h.getGlob(i+1);
        loc_h.getGlob(i)   = loc_cos*val_m+loc_sin*val_p;
        loc_h.getGlob(i+1) =-loc_sin*val_m+loc_cos*val_p;
    }
    DType h_m = loc_h.getGlob(dim_v-1);
    DType h_p = loc_h.getGlob(dim_v);
    DType l_cos = h_m/sqrt(h_m*h_m+h_p*h_p);
    DType l_sin = h_p/sqrt(h_m*h_m+h_p*h_p);
    gmres_cos.push_back(l_cos);
    gmres_sin.push_back(l_sin);
    loc_h.getGlob(dim_v-1) = l_cos*h_m+l_sin*h_p;
    loc_h.getGlob(dim_v)   =-l_sin*h_m+l_cos*h_p; // = 0.0
    DType val_m = gmres_g->getGlob(dim_v-1);
    DType val_p = gmres_g->getGlob(dim_v);
    gmres_g->getGlob(dim_v-1) = l_cos*val_m+l_sin*val_p;
    gmres_g->getGlob(dim_v)   =-l_sin*val_m+l_cos*val_p;
    
    // Fill last column of matrix h
    for (int i=0; i<dim_v; i++)
        gmres_h->get(i,dim_v-1) = loc_h.getGlob(i);
    HVEC g_trunc(dim_v,0.0);
    for (int i=0; i<dim_v; i++)
        g_trunc[i] = gmres_g->getGlob(i);
    HVEC loc_y(dim_v,0.0);
    loc_y = (gmres_h->inv())*g_trunc;
    LinVec loc_u(dim,size);
    for (int i=dim_v-1; i>=0; i--) {
        loc_u.getGlob(i) += loc_y[i];
        loc_u.proj(*vec_w[i]);
    }
    sol->copy(*gmres_sol0);
    sol->scale_inc(1.0,1.0,loc_u);
    DType cur_res = compute_res();
    DType cur_tot_res = compute_tot_res(cur_res);
    if (cur_tot_res<=target_res) return;

    // v
    gmres_v->reset(0.0);
    gmres_v->getGlob(dim_v) = 1.0;      // e_j
    for (int i=dim_v; i>=0; i--)
        gmres_v->proj(*vec_w[i]);

    // z
    mat->xVec(gmres_v,gmres_z);
    for (int i=0; i<=dim_v; i++)
        gmres_z->proj(*vec_w[i]);
    return;
}

void LinSolve::solve_krylov_dqgmres_pre(void) {
    clear_krylov_spaces(vec_v);
    clear_krylov_spaces(vec_p);
    gmres_cos.clear();
    gmres_sin.clear();

    gmres_gamma = res->getNorm();
    LinVec* loc_v = new LinVec(dim,size);
    loc_v->copy(*res);
    loc_v->scale(1.0/gmres_gamma);
    vec_v.push_back(loc_v);
    gmres_cos.push_back(1.0);
    gmres_sin.push_back(0.0);
    LinVec* loc_p = new LinVec(dim,size);
    vec_p.push_back(loc_p);
    return;
}

void LinSolve::solve_krylov_dqgmres(void) {
    LinVec* loc_v = vec_v.back();
    LinVec* loc_v_new = new LinVec(dim,size);
    LinVec* loc_w = new LinVec(dim,size);

    mat->xVec(loc_v,loc_w);
    int size_v = vec_v.size();
    HVEC loc_h(size_v+2,0.0);
    for (int i=0; i<size_v; i++) {
        loc_h[i+1] = calc_inner_prod(*loc_w,*vec_v[i]);
        loc_w->scale_inc(1.0,-loc_h[i+1],*vec_v[i]);
    }
    loc_h[size_v+1] = loc_w->getNorm();
    if (loc_h[size_v+1]==0.0) {
        kspace_signal = true;
        delete loc_v_new;
        delete loc_w;
        return;
    }
    loc_v_new->copy(*loc_w);
    loc_v_new->scale(1.0/loc_h[size_v+1]);
    vec_v.push_back(loc_v_new);
    if (vec_v.size()>(unsigned)kspace_dim) {
        delete vec_v[0];
        vec_v.erase(vec_v.begin());
    }   // The size of vec_v is no larger than kspace_dim

    // Rotate loc_h
    int size_rot = gmres_cos.size();
    for (int i=0; i<size_rot; i++) {
        DType val_m = loc_h[i];
        DType val_p = loc_h[i+1];
        DType l_cos = gmres_cos[i];
        DType l_sin = gmres_sin[i];
        loc_h[i]   = l_cos * val_m + l_sin * val_p;
        loc_h[i+1] =-l_sin * val_m + l_cos * val_p;
    }
    // New rotation
    DType val_m = loc_h[size_v];
    DType val_p = loc_h[size_v+1];
    DType l_cos = val_m/sqrt(val_m*val_m+val_p*val_p);
    DType l_sin = val_p/sqrt(val_m*val_m+val_p*val_p);
    loc_h[size_v]   = l_cos*val_m+l_sin*val_p;
    loc_h[size_v+1] = 0.0;
    gmres_cos.push_back(l_cos);
    gmres_sin.push_back(l_sin);
    if (gmres_cos.size()>(unsigned)kspace_dim)
        gmres_cos.erase(gmres_cos.begin());
    if (gmres_sin.size()>(unsigned)kspace_dim)
        gmres_sin.erase(gmres_sin.begin());

    LinVec* loc_p = new LinVec(dim,size);
    loc_p->copy(*loc_v);
    for (unsigned i=0; i<vec_p.size(); i++)
        loc_p->scale_inc(1.0,-loc_h[i],*vec_p[i]);
    loc_p->scale(1.0/loc_h[size_v]);
    sol->scale_inc(1.0,l_cos*gmres_gamma,*loc_p);
    gmres_gamma *= -l_sin;  // Update gamma
    vec_p.push_back(loc_p);
    if (vec_p.size()>(unsigned)kspace_dim) {
        delete vec_p[0];
        vec_p.erase(vec_p.begin());
    }
    delete loc_w;
    return;
}

void LinSolve::solve_krylov_cg_pre(void) {
    if (vec_cg_r==NULL) vec_cg_r = new LinVec(dim,size);
    vec_cg_r->copy(*res);          // r0_star

    if (vec_cg_p==NULL) vec_cg_p = new LinVec(dim,size);
    vec_cg_p->copy(*res);          // p0
    return;
}

void LinSolve::solve_krylov_cg(void) {
    LinVec* vec_s = new LinVec(dim,size);
    LinVec* vec_Ap = new LinVec(dim,size);
    LinVec* vec_As = new LinVec(dim,size);

    DType a0 = calc_inner_prod(*res,*vec_cg_r);
    // vec_Ap = Ap
    mat->xVec(vec_cg_p,vec_Ap);
    DType a = a0/calc_inner_prod(*vec_Ap,*vec_cg_r);

    vec_s->copy(*res);
    vec_s->scale_inc(1.0,-a,*vec_Ap);
    mat->xVec(vec_s,vec_As);
    DType norm_As = vec_As->getNorm();
    DType w = calc_inner_prod(*vec_As,*vec_s)/(norm_As*norm_As);

    sol->scale_inc(1.0,a,*vec_cg_p);
    sol->scale_inc(1.0,w,*vec_s);    // Update sol

    res->copy(*vec_s);
    res->scale_inc(1.0,-w,*vec_As);  // Update res

    DType b = calc_inner_prod(*res,*vec_cg_r)/a0;
    b *= (a/w);

    vec_cg_p->scale_inc(1.0,-w,*vec_Ap);
    vec_cg_p->scale(b);
    vec_cg_p->scale_inc(1.0,1.0,*res); // Update vec_cg_p
    delete vec_s;
    delete vec_As;
    delete vec_Ap;
    return;
}

void LinSolve::solve_krylov_lanczos_pre(void) {
    clear_krylov_spaces(vec_v);
    // Size of vec_w is 2
    clear_krylov_spaces(vec_w);
    lanczos_alpha.clear();
    lanczos_delta.clear();
    lanczos_beta.clear();

    lanczos_beta0 = res->getNorm();

    LinVec* loc_v0 = new LinVec(dim,size);
    LinVec* loc_w0 = new LinVec(dim,size);

    vec_v.push_back(loc_v0);
    vec_w.push_back(loc_w0);

    LinVec* loc_v = new LinVec(dim,size);
    LinVec* loc_w = new LinVec(dim,size);

    loc_v->copy(*res);
    loc_v->scale(1.0/lanczos_beta0);
    loc_w->copy(*loc_v);

    vec_v.push_back(loc_v);
    vec_w.push_back(loc_w);

    lanczos_beta.push_back(0.0);
    lanczos_delta.push_back(0.0);

    if (lanczos_sol0==NULL)
        lanczos_sol0 = new LinVec(dim,size);
    lanczos_sol0->copy(*sol);
    return;
}

void LinSolve::solve_krylov_lanczos(void) {
    LinVec* loc_v = vec_v.back();
    LinVec* loc_w = vec_w[1];
    LinVec* loc_v_pre = vec_v[vec_v.size()-2];
    LinVec* loc_w_pre = vec_w[0];
    LinVec* loc_v_new = new LinVec(dim,size);
    LinVec* loc_w_new = new LinVec(dim,size);
    DType   beta_pre  = lanczos_beta.back();
    DType   delta_pre = lanczos_delta.back();

    LinVec*  vec_Av = new LinVec(dim,size);
    mat->xVec(loc_v,vec_Av);
    DType alpha = calc_inner_prod(*vec_Av,*loc_w);
    lanczos_alpha.push_back(alpha);

    // v hat
    loc_v_new->copy(*vec_Av);
    loc_v_new->scale_inc(1.0,-alpha,*loc_v);
    loc_v_new->scale_inc(1.0,-beta_pre,*loc_v_pre);
    delete vec_Av;

    // w hat
    mat->xVec_t(loc_w,loc_w_new);
    loc_w_new->scale_inc(1.0,-alpha,*loc_w);
    loc_w_new->scale_inc(1.0,-delta_pre,*loc_w_pre);

    // delta
    DType inner_vw = calc_inner_prod(*loc_v_new,*loc_w_new);
    DType delta = sqrt(fabs(inner_vw));
    lanczos_delta.push_back(delta);

    int dim_t = lanczos_alpha.size();
    TriDMAT mat_t(dim_t,0.0);
    for (int i=0; i<dim_t; i++)
        mat_t.mdiag[i] = lanczos_alpha[i];
    for (int i=0; i<dim_t-1; i++) {
        mat_t.upper[i] = lanczos_beta[i+1];
        mat_t.lower[i] = lanczos_delta[i+1];
    }
    HMAT inv_mat_t = mat_t.inv();
    // y
    HVEC vec_y = lanczos_beta0 * inv_mat_t.getCol(0);
    sol->copy(*lanczos_sol0);
    for (int i=0; i<dim_t; i++)
        sol->scale_inc(1.0,vec_y[i],*vec_v[i+1]);
    DType cur_res = compute_res();
    DType cur_tot_res = compute_tot_res(cur_res);
    if (cur_tot_res<=target_res) {
        delete loc_v_new;
        delete loc_w_new;
        return;
    }
    if ((fabs(delta)<=TOL_ZERO)||
        ((kspace_dim>0)&&(vec_v.size()>=(unsigned)kspace_dim))) {    // Restart
        delete loc_v_new;
        delete loc_w_new;
        solve_krylov_lanczos_pre();
        return;
    }
    // beta
    DType beta = inner_vw/delta;
    lanczos_beta.push_back(beta);
    // v_new and w_new
    loc_v_new->scale(1.0/delta);
    loc_w_new->scale(1.0/beta);
    vec_v.push_back(loc_v_new);

    // Maintain the size of vec_w to be 2
    delete vec_w[0];
    vec_w.erase(vec_w.begin());
    vec_w.push_back(loc_w_new);
    return;
}

void LinSolve::solve_krylov_qmr_pre(void) {
    // The size of vec_w is 2
    clear_krylov_spaces(vec_v);
    // The size of vec_w is 2
    clear_krylov_spaces(vec_w);
    // In the case of QMR, vec_p has size 2
    clear_krylov_spaces(vec_p);
    
    qmr_gamma = res->getNorm(); // gamma_1
    LinVec* loc_v0 = new LinVec(dim,size);
    LinVec* loc_w0 = new LinVec(dim,size);
    vec_v.push_back(loc_v0);
    vec_w.push_back(loc_w0);

    LinVec* loc_v = new LinVec(dim,size);
    LinVec* loc_w = new LinVec(dim,size);
    loc_v->copy(*res);
    loc_v->scale(1.0/qmr_gamma);
    loc_w->copy(*loc_v);
    vec_v.push_back(loc_v);
    vec_w.push_back(loc_w);

    qmr_beta  = 0.0;
    qmr_delta = 0.0;

    LinVec* loc_p0 = new LinVec(dim,size);
    LinVec* loc_p  = new LinVec(dim,size);
    vec_p.push_back(loc_p0);
    vec_p.push_back(loc_p);

    qmr_cos[0] = 1.0;
    qmr_cos[1] = 1.0;
    qmr_sin[0] = 0.0;
    qmr_sin[1] = 0.0;

    qmr_iter = 0;
    return;
}

void LinSolve::solve_krylov_qmr(void) {
    LinVec* loc_v = vec_v[1];
    LinVec* loc_w = vec_w[1];
    LinVec* loc_p = vec_p[1];
    LinVec* loc_v_pre = vec_v[0];
    LinVec* loc_w_pre = vec_w[0];
    LinVec* loc_p_pre = vec_p[0];
    LinVec* loc_v_new = new LinVec(dim,size);
    LinVec* loc_w_new = new LinVec(dim,size);
    LinVec* loc_p_new = new LinVec(dim,size);
    DType   beta_pre  = qmr_beta;
    DType   delta_pre = qmr_delta;
    
    qmr_iter++;

    LinVec* vec_Av = new LinVec(dim,size);
    mat->xVec(loc_v,vec_Av);
    DType alpha = calc_inner_prod(*vec_Av,*loc_w);

    // v hat
    loc_v_new->copy(*vec_Av);
    loc_v_new->scale_inc(1.0,-alpha,*loc_v);
    loc_v_new->scale_inc(1.0,-beta_pre,*loc_v_pre);
    
    // w hat
    mat->xVec_t(loc_w,loc_w_new);
    loc_w_new->scale_inc(1.0,-alpha,*loc_w);
    loc_w_new->scale_inc(1.0,-delta_pre,*loc_w_pre);

    // delta
    DType inner_vw = calc_inner_prod(*loc_v_new,*loc_w_new);
    qmr_delta = sqrt(fabs(inner_vw));

    DType last_t[4] = {0.0, 0.0, 0.0, 0.0};
    int dim_t = qmr_iter;       // m = dim_t
    last_t[2] = alpha;          // unrotated t_{;m  , m}
    last_t[3] = qmr_delta;      // unrotated t_{;m+1, m}
    if (dim_t==1) {
        // do nothing about the rotation
    }
    if (dim_t==2) {
        // rotate using qmr_cos[1] and qmr_sin[1]
        //last_t[1] = beta_pre; // unrotated t_{;m-1, m}
        last_t[1] =  qmr_cos[1]*beta_pre+qmr_sin[1]*alpha;
        last_t[2] = -qmr_sin[1]*beta_pre+qmr_cos[1]*alpha;
    }
    else {
        // first rotate using qmr_cos[0] and qmr_sin[0]
        // then  rotate using qmr_cos[1] and qmr_sin[1]
        //last_t[0] = 0.0;      // unrotated t_{;m-2, m}
        //last_t[1] = beta_pre; // unrotated t_{;m-1, m}
        last_t[0] = qmr_sin[0]*beta_pre;
        last_t[1] = qmr_cos[1]*qmr_cos[0]*beta_pre+
                    qmr_sin[1]*alpha;
        last_t[2] =-qmr_sin[1]*qmr_cos[0]*beta_pre+
                    qmr_cos[1]*alpha;
    }
    qmr_cos[0] = qmr_cos[1];
    qmr_sin[0] = qmr_sin[1];
    DType temp_last_t = sqrt(last_t[2]*last_t[2]+last_t[3]*last_t[3]);
    qmr_cos[1] = last_t[2]/temp_last_t;
    qmr_sin[1] = last_t[3]/temp_last_t;

    DType val_m = last_t[2];
    DType val_p = last_t[3];
    last_t[2] = qmr_cos[1]*val_m+qmr_sin[1]*val_p;
    last_t[3] =-qmr_sin[1]*val_m+qmr_cos[1]*val_p;

    // Update vec_p
    loc_p_new->copy(*loc_v);
    loc_p_new->scale_inc(1.0,-last_t[0],*loc_p_pre);
    loc_p_new->scale_inc(1.0,-last_t[1],*loc_p);
    loc_p_new->scale(1.0/last_t[2]);

    // Maintain the size of vec_p
    delete vec_p[0];
    vec_p.erase(vec_p.begin());
    vec_p.push_back(loc_p_new);

    // Update sol
    sol->scale_inc(1.0,qmr_cos[1]*qmr_gamma,*loc_p_new);
    qmr_gamma = -qmr_sin[1]*qmr_gamma;
    DType cur_res = compute_res();
    DType cur_tot_res = compute_tot_res(cur_res);
    if (cur_tot_res<=target_res) {
        delete loc_v_new;
        delete loc_w_new;
        delete vec_Av;
        return;
    }
    if (fabs(qmr_delta)<=TOL_ZERO) {    // Restart
        delete loc_v_new;
        delete loc_w_new;
        delete vec_Av;
        solve_krylov_qmr_pre();
        return;
    }
    // beta
    qmr_beta = inner_vw/qmr_delta;
    // v_new and w_new
    loc_v_new->scale(1.0/qmr_delta);
    loc_w_new->scale(1.0/qmr_beta);
    delete vec_v[0];
    vec_v.erase(vec_v.begin());
    vec_v.push_back(loc_v_new);
    delete vec_w[0];
    vec_w.erase(vec_w.begin());
    vec_w.push_back(loc_w_new);
    delete vec_Av;
    return;
}

void LinSolve::linkMat(LinMat*& mat_in) {
    assert(mat_in->isSquare());
    if (dim<1)
        dim = mat_in->getBlckDim();
    else
        assert(dim==mat_in->getBlckDim());
    if (size<1)
        size = mat_in->getSize();
    else
        assert(size==mat_in->getSize());
    mat = mat_in;
    return;
}

void LinSolve::linkPrc(LinMat*& prc_in) {
    assert(prc_in->isSquare());
    if (dim<1)
        dim = prc_in->getBlckDim();
    else
        assert(dim==prc_in->getBlckDim());
    if (size<1)
        size = prc_in->getSize();
    else
        assert(size==prc_in->getSize());
    prc = prc_in;
    return;
}

void LinSolve::linkSol(LinVec*& sol_in) {
    if (dim<1)
        dim = sol_in->getBlckDim();
    else
        assert(dim==sol_in->getBlckDim());
    if (size<1)
        size = sol_in->getSize();
    else
        assert(size==sol_in->getSize());
    sol = sol_in;
    return;
}

void LinSolve::linkRhs(LinVec*& rhs_in) {
    if (dim<1)
        dim = rhs_in->getBlckDim();
    else
        assert(dim==rhs_in->getBlckDim());
    if (size<1)
        size = rhs_in->getSize();
    else
        assert(size==rhs_in->getSize());
    rhs = rhs_in;
    return;
}

void LinSolve::linkRes(LinVec*& res_in) {
    if (dim<1)
        dim = res_in->getBlckDim();
    else
        assert(dim==res_in->getBlckDim());
    if (size<1)
        size = res_in->getSize();
    else
        assert(size==res_in->getSize());
    res = res_in;
    return;
}

bool LinSolve::solve(void) {
    if (preCond) (this->*apply_prec)(prec_ilu_fill);
    bool toContinue = solve_start();
    int  cur_iter = 0;
    if (print_stat)
        fprintf(stdout,"Iter 0: Current Residual %10.20f\n",res_hist.back());
    switch (lsysType) {
        case LSYS_TYPE_DIR:
            solve_dir(toContinue,cur_iter);
            break;
        case LSYS_TYPE_ITER:
            solve_iter(toContinue,cur_iter);
            break;
        case LSYS_TYPE_KRYLOV:
            solve_krylov(toContinue,cur_iter);
            break;
        default:
            break;
    }
    return (res_hist.back()<=target_res);
}

HVEC LinSolve::convert(LinVec*& vec_in) {
    int dim  = vec_in->getBlckDim();
    int size = vec_in->getSize();
    HVEC vec(dim*size,0.0);
    for (int i=0; i<dim*size; i++)
        vec[i] = vec_in->getGlobVal(i);
    return vec;
}

HMAT LinSolve::convert(LinMat*& mat_in) {
    assert(mat_in->isSquare());
    int dim  = mat_in->getBlckDim();
    int size = mat_in->getSize();
    HMAT mat(dim*size,0.0);
    for (int i=0; i<dim*size; i++)
        for (int j=0; j<dim*size; j++)
            mat[i][j] = mat_in->getGlobVal(i,j);
    return mat;
}

void LinSolve::convert(LinVec*& vec_ot, HVEC vec_in, int dim, int size) {
    assert(vec_in.size()==dim*size);
    assert(vec_ot->getBlckDim()==dim);
    assert(vec_ot->getSize()==size);
    for (int n=0; n<size; n++) {
        VComp v_comp(dim,0.0);
        for (int i=0; i<dim; i++)
            v_comp.get(i) = vec_in[n*dim+i];
        vec_ot->set(n,v_comp);
    }
    return;
}

void LinSolve::convert(LinMat*& mat_ot, HMAT mat_in, int dim, int size, DType conv_tol) {
    assert(mat_ot->isSquare());
    assert(mat_in.size()==dim*size);
    assert(mat_ot->getBlckDim()==dim);
    assert(mat_ot->getSize()==size);
    bool spse = mat_ot->isSparse();
    if (spse) mat_ot->reset();
    for (int n1=0; n1<size; n1++)
        for (int n2=0; n2<size; n2++) {
            Entry m_entry(dim,0.0);
            for (int i1=0; i1<dim; i1++)
                for (int i2=0; i2<dim; i2++) {
                    m_entry.get(i1,i2) = mat_in[n1*dim+i1][n2*dim+i2];
                }
            if (spse) {
                if (m_entry.getNorm()>fabs(conv_tol))
                    mat_ot->set(n1,n2,m_entry);
            }
            else
                mat_ot->set(n1,n2,m_entry);
        }
}


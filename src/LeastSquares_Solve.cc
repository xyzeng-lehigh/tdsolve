#include <cassert>
#include "LeastSquares_Solve.h"
#include "Quad_Coef.h"

/*
 * General least squares solve
 */

LstSqr::LstSqr(void) {
    sol = NULL;
    rhs = NULL;
    mat = NULL;
    wgt = NULL;
    coef_mat = NULL;
    coef_inv = NULL;
    coef_rhs = NULL;
    c_mat = NULL;
    c_rhs = NULL;
    c_lag = NULL;
    cstr_mat = NULL;
    cstr_inv = NULL;
    cstr_rhs = NULL;
    reset();
}

LstSqr::~LstSqr(void) {
    clear();
}

void LstSqr::clear(void) {
    if (sol!=NULL) {delete sol; sol = NULL;}
    if (rhs!=NULL) {delete rhs; rhs = NULL;}
    if (mat!=NULL) {
        for (int i=0; i<Neqns; i++)
            if (mat[i]!=NULL) delete mat[i];
        delete[] mat;
        mat = NULL;
    }
    if (wgt!=NULL) {delete wgt; wgt = NULL;}

    if (coef_mat!=NULL) {delete coef_mat; coef_mat = NULL;}
    if (coef_inv!=NULL) {delete coef_inv; coef_inv = NULL;}
    if (coef_rhs!=NULL) {delete coef_rhs; coef_rhs = NULL;}

    if (c_mat!=NULL) {
        for (int i=0; i<Ncstr; i++)
            if (c_mat[i]!=NULL) delete c_mat[i];
        delete[] c_mat;
        c_mat = NULL;
    }
    if (c_rhs!=NULL) {delete c_rhs; c_rhs = NULL;}
    if (c_lag!=NULL) {delete c_lag; c_lag = NULL;}

    if (cstr_mat!=NULL) {delete cstr_mat; cstr_mat = NULL;}
    if (cstr_inv!=NULL) {delete cstr_inv; cstr_inv = NULL;}
    if (cstr_rhs!=NULL) {delete cstr_rhs; cstr_rhs = NULL;}

    return;
}

void LstSqr::cstr_solve(void) {
    // Construct cstr_mat = c_mat * coef_inv * c_mat^T
    *cstr_mat = 0.0;
    for (int i=0; i<Ncstr; i++)
        for (int j=0; j<Ncstr; j++) {
            DType loc_var = 0.0; 
            for (int k1=0; k1<Nvars; k1++) {
                DType loc_multip = (*c_mat[i])[k1];
                for (int k2=0; k2<Nvars; k2++)
                    loc_var += loc_multip * (*coef_inv)[k1][k2] * (*c_mat[j])[k2];
            }
            (*cstr_mat)[i][j] = loc_var;
        }
    *cstr_inv = cstr_mat->inv();
    // Construct cstr_rhs = c_mat * coef_inv * coef_rhs - c_rhs
    *cstr_rhs = 0.0;
    for (int i=0; i<Ncstr; i++) {
        DType loc_var = -(*c_rhs)[i];
        for (int k1=0; k1<Nvars; k1++) {
            DType loc_multip = (*c_mat[i])[k1];
            for (int k2=0; k2<Nvars; k2++)
                loc_var += loc_multip * (*coef_inv)[k1][k2] * (*coef_rhs)[k2];
        }
        (*cstr_rhs)[i] = loc_var;
    }
    // Compute Lagrangian multipliers
    *c_lag = (*cstr_inv) * (*cstr_rhs);
    return;
}

void LstSqr::reset(void) {
    Nvars    = 0;
    Neqns    = 0;
    Ncstr    = 0;
    use_cstr = false;
    clear();
}

void LstSqr::solve(void) {
    // Construct coef_mat = mat^T * wgt * wgt * mat
    *coef_mat = 0.0;
    for (int k=0; k<Neqns; k++) {
        DType loc_wgt = (*wgt)[k] * (*wgt)[k];
        *coef_mat += loc_wgt * ((*mat[k])|(*mat[k]));
    }
    *coef_inv = coef_mat->inv();
    // Construct coef_rhs = mat^T * wgt * wgt * rhs
    *coef_rhs = 0.0;
    for (int k=0; k<Neqns; k++) {
        DType loc_wgt = (*wgt)[k] * (*wgt)[k];
        *coef_rhs += loc_wgt * (*rhs)[k] * (*mat[k]);
    }
    // Construct constraints and solve for 
    // Lagrangian multipliers if needed
    if (use_cstr)   cstr_solve();
    // Copmute the solution to the least squares
    // problem
    if (use_cstr) {
        // Update coef_rhs for multipliers
        for (int i=0; i<Nvars; i++) {
            DType loc_var = (*coef_rhs)[i];
            for (int k=0; k<Ncstr; k++)
                loc_var -= (*c_mat[k])[i] * (*c_lag)[k];
            (*coef_rhs)[i] = loc_var;
        }
    }
    *sol = (*coef_inv) * (*coef_rhs);
    return;
}

void LstSqr::setup(void) {
    clear();

    assert(Nvars>0);
    sol      = new HVEC(Nvars,0.0);
    coef_mat = new HMAT(Nvars,0.0);
    coef_inv = new HMAT(Nvars,0.0);
    coef_rhs = new HVEC(Nvars,0.0);

    assert(Neqns>0);
    rhs = new HVEC(Neqns,0.0);
    mat = new HVEC*[Neqns];
    for (int i=0; i<Neqns; i++)
        mat[i] = new HVEC(Nvars,0.0);
    wgt = new HVEC(Neqns,1.0);

    if (use_cstr) {
        assert(Ncstr>0);
        c_mat    = new HVEC*[Ncstr];
        for (int i=0; i<Ncstr; i++)
            c_mat[i] = new HVEC(Nvars,0.0);
        c_rhs    = new HVEC(Ncstr,0.0);
        c_lag    = new HVEC(Ncstr,0.0);
        cstr_mat = new HMAT(Ncstr,0.0);
        cstr_inv = new HMAT(Ncstr,0.0);
        cstr_rhs = new HVEC(Ncstr,0.0);
    }

    return;
}

void LstSqr::setEqnsWeight(int i, DType wgt_in) {
    (*wgt)[i] = wgt_in;
    return;
}

void LstSqr::setEqnsCoef(int i, int j, DType coef_in) {
    (*mat[i])[j] = coef_in;
    return;
}

void LstSqr::setEqnsCoef(int i, HVEC* coef_vec) {
    *mat[i] = *coef_vec;
    return;
}

void LstSqr::setEqnsRhs(int i, DType val_in) {
    (*rhs)[i] = val_in;
    return;
}

void LstSqr::setEqns(int i, HVEC* coef_vec, DType val_in) {
    *mat[i]   = *coef_vec;
    (*rhs)[i] = val_in;
    return;
}

void LstSqr::setCstrCoef(int i, int j, DType coef_in) {
    (*c_mat[i])[j] = coef_in;
    return;
}

void LstSqr::setCstrCoef(int i, HVEC* coef_vec) {
    *c_mat[i] = *coef_vec;
    return;
}

void LstSqr::setCstrVal(int i, DType val_in) {
    (*c_rhs)[i] = val_in;
    return;
}

void LstSqr::setCstr(int i, HVEC* coef_vec, DType val_in) {
    *c_mat[i]   = *coef_vec;
    (*c_rhs)[i] = val_in;
    return;
}

/*
 * Local polynomial coefficients
 */

// Base case

PolyCoef::PolyCoef(void) {
    p_order = -1;
    p_coef  = NULL;
    n_coef  = NULL;
}

PolyCoef::~PolyCoef(void) {
    clear();
}

void PolyCoef::clear(void) {
    if (p_coef!=NULL) {
        assert(p_order>=0);
        for (int i=0; i<p_order+1; i++)
            if (p_coef[i]!=NULL)
                delete[] p_coef[i];
        delete[] p_coef;
    }
    p_coef = NULL;
    if (n_coef!=NULL) delete[] n_coef;
    n_coef = NULL;
    return;
}

void PolyCoef::setup(int p_order_in) {
    assert(p_order_in>=0);
    clear();
    p_order = p_order_in;
    p_coef = new DType*[p_order+1];
    n_coef = new int[p_order+1];
    for (int i=0; i<p_order+1; i++) n_coef[i] = 0;
    return;
}

// set the loc_idx'th coefficient in the order'th array
void PolyCoef::setLocCoef(
        DType coef, int order, int loc_idx) {
    assert((order>=0)&&(order<=p_order));
    assert((loc_idx>=0)&&(loc_idx<n_coef[order]));
    p_coef[order][loc_idx] = coef;
    return;
}

// set the global glob_idx'th coefficient
void PolyCoef::setGlobCoef(
        DType coef, int glob_idx) {
    assert((glob_idx>=0)&&(glob_idx<Ncoef));
    int order = 0;
    while ((1+glob_idx)>n_coef[order]) {
        glob_idx -= n_coef[order];
        order++;
    }
    p_coef[order][glob_idx] = coef;
    return;
}

// 1D case

void PolyCoef1D::setup(int p_order_in) {
    PolyCoef::setup(p_order_in);
    for (int i=0; i<p_order+1; i++) {
        p_coef[i]    = new DType[1];
        n_coef[i]    = 1;
        p_coef[i][0] = 0.0;
    }
    Ncoef = p_order+1;
    return;
}

DType PolyCoef1D::getCoef(int order) {
    assert((order>=0)&&(order<=p_order));
    return p_coef[order][0];
}

int PolyCoef1D::getNumCoef(int order) {
    assert((order>=0)&&(order<=p_order));
    return n_coef[order];
}

int PolyCoef1D::getLocCoefIdx(int order) {
    assert((order>=0)&&(order<=p_order));
    return 0;
}

int PolyCoef1D::getGlobCoefIdx(int order) {
    assert((order>=0)&&(order<=p_order));
    return order;
}

void PolyCoef1D::globIdx2Der(int idx, int& o_x) {
    assert((idx>=0)&&(idx<Ncoef));
    o_x = idx;
    return;
}

void PolyCoef1D::setCoef(DType coef, int order) {
    assert((order>=0)&&(order<=p_order));
    p_coef[order][0] = coef;
    return;
}

// 2D case

void PolyCoef2D::setup(int p_order_in) {
    PolyCoef::setup(p_order_in);
    Ncoef = 0;
    for (int i=0; i<p_order+1; i++) {
        p_coef[i] = new DType[i+1];
        n_coef[i] = i+1;
        Ncoef    += i+1;
        for (int j=0; j<i+1; j++) p_coef[i][j] = 0.0;
    }
    return;
}

DType PolyCoef2D::getCoef(int o_x, int o_y) {
    int order = o_x+o_y;
    assert((o_x>=0)&&(o_y>=0)&&(order<=p_order));
    return p_coef[order][getLocCoefIdx(o_x,o_y)];
}

int PolyCoef2D::getLocCoefIdx(int o_x, int o_y) {
    assert((o_x>=0)&&(o_y>=0)&&(o_x+o_y<=p_order));
    // local index in p_coef[o_x+o_y]:
    // p_coef[o_x+o_y][o_y] ~ (d/dx)^o_x
    //                        (d/dy)^o_y
    return o_y;
}

int PolyCoef2D::getGlobCoefIdx(int o_x, int o_y) {
    int order = o_x+o_y;
    assert((o_x>=0)&&(o_y>=0)&&(order<=p_order));
    return ((order*(order+1))/2) + 
           getLocCoefIdx(o_x,o_y);
}

void PolyCoef2D::globIdx2Der(
        int idx, int& o_x, int& o_y) {
    assert((idx>=0)&&(idx<Ncoef));
    o_x = 0;
    while (getGlobCoefIdx(0,o_x)<idx) o_x++;
    o_y = idx-getGlobCoefIdx(o_x,0);
    o_x = o_x - o_y;
    return;
}

void PolyCoef2D::setCoef(
        DType coef, int o_x, int o_y) {
    int order = o_x+o_y;
    assert((o_x>=0)&&(o_y>=0)&&(order<=p_order));
    p_coef[order][getLocCoefIdx(o_x,o_y)] = coef;
    return;
}

// 3D case

void PolyCoef3D::setup(int p_order_in) {
    PolyCoef::setup(p_order_in);
    Ncoef = 0;
    for (int i=0; i<p_order+1; i++) {
        n_coef[i] = ((i+1)*(i+2))/2;
        p_coef[i] = new DType[n_coef[i]];
        Ncoef    += n_coef[i];
        for (int j=0; j<n_coef[i]; j++) p_coef[i][j] = 0.0;
    }
    return;
}

DType PolyCoef3D::getCoef(
        int o_x, int o_y, int o_z) {
    int order = o_x+o_y+o_z;
    assert((o_x>=0)&&(o_y>=0)&&(o_z>=0)&&
           (order<=p_order));
    return p_coef[order][getLocCoefIdx(o_x,o_y,o_z)];
}

int PolyCoef3D::getLocCoefIdx(
        int o_x, int o_y, int o_z) {
    int o_yz = o_y+o_z;
    assert((o_x>=0)&&(o_y>=0)&&(o_z>=0)&&
           (o_x+o_yz<=p_order));
    // local index in p_ceof[o_x+o_y+o_z]:
    // p_coef[o_x+o_y+o_z]
    //       [(o_y+o_z)(o_y+o_z+1)/2+o_z] ~ 
    //           (d/dx)^o_x
    //           (d/dy)^o_y
    //           (d/dz)^o_z
    return ((o_yz*(o_yz+1))/2) + o_z;
}

int PolyCoef3D::getGlobCoefIdx(
        int o_x, int o_y, int o_z) {
    int o_yz  = o_y+o_z;
    int order = o_x+o_yz;
    assert((o_x>=0)&&(o_y>=0)&&(o_z>=0)&&
           (order<=p_order));
    return ((order*(order+1)*(order+2))/6) +
           getLocCoefIdx(o_x,o_y,o_z);
}

void PolyCoef3D::globIdx2Der(int idx, 
        int& o_x, int& o_y, int& o_z) {
    assert((idx>=0)&&(idx<Ncoef));
    o_x = 0;
    while (getGlobCoefIdx(0,0,o_x)<idx) o_x++;
    o_y = 0;
    while (getGlobCoefIdx(o_x,0,o_y)<idx) o_y++;
    o_z = idx-getGlobCoefIdx(o_x,o_y,0);
    o_y = o_y - o_z;
    o_x = o_x - o_y - o_z;
    return;
}

void PolyCoef3D::setCoef(
        DType coef, int o_x, int o_y, int o_z) {
    p_coef[o_x+o_y+o_z][getLocCoefIdx(o_x,o_y,o_z)] = coef;
    return;
}

/*
 * Local least squares solver for polynomial coefficients
 */

// Base case

LocLstSqr::LocLstSqr(void) {
    p_order = -1;
    Npars   = 0;
    use_val = false;
    coef_rhs = NULL;
    coef_par = NULL;
    coef_mat = NULL;
    coef_inv = NULL;

    Ncstr    = 0;
    use_cstr = false;
    cstr_vars.clear();
    cstr_vals.clear();
    cstr_mat = NULL;
    cstr_inv = NULL;
    cstr_lag = NULL;
    cstr_rhs = NULL;
}

void LocLstSqr::clear(void) {
    p_order = -1;
    Npars   = 0;
    use_val = false;
    if (coef_rhs!=NULL) delete coef_rhs;
    if (coef_par!=NULL) delete coef_par;
    if (coef_mat!=NULL) delete coef_mat;
    if (coef_inv!=NULL) delete coef_inv;
    coef_rhs = NULL;
    coef_par = NULL;
    coef_mat = NULL;
    coef_inv = NULL;

    cstr_reset();
    return;
}

void LocLstSqr::add2Arrays(DType loc_val, 
        DType* coef_int, DType weight) {
    for (int i=0; i<Npars; i++) {
        (*coef_rhs)[i] += loc_val * coef_int[i] * weight;
        for (int j=0; j<Npars; j++)
            (*coef_mat)[i][j] += coef_int[i] * coef_int[j] * weight;
    }
    return;
}

void LocLstSqr::add2Mat(
        DType* coef_int, DType weight) {
    for (int i=0; i<Npars; i++)
        for (int j=0; j<Npars; j++)
            (*coef_mat)[i][j] += coef_int[i] * coef_int[j] * weight;
    return;
}

void LocLstSqr::add2Rhs(DType loc_val, 
        DType* coef_int, DType weight) {
    for (int i=0; i<Npars; i++)
        (*coef_rhs)[i] += loc_val * coef_int[i] * weight;
    return;
}

void LocLstSqr::add2CstrLists(DType loc_val, DType* coef_int) {
    cstr_vals.push_back(loc_val);
    for (int i=0; i<Npars; i++)
        cstr_vars.push_back(coef_int[i]);
    return;
}

void LocLstSqr::resetMats(void) {
    if (coef_mat!=NULL) (*coef_mat) = 0.0;
    if (coef_par!=NULL) (*coef_par) = 0.0;
    if (coef_inv!=NULL) delete coef_inv;
}

void LocLstSqr::resetRhs(void) {
    if (coef_rhs!=NULL) (*coef_rhs) = 0.0;
    if (coef_par!=NULL) (*coef_par) = 0.0;
    cstr_reset();
}

void LocLstSqr::cstr_reset(void) {
    Ncstr    = 0;
    use_cstr = false;
    cstr_vars.clear();
    cstr_vals.clear();
    if (cstr_mat!=NULL) delete cstr_mat;
    if (cstr_inv!=NULL) delete cstr_inv;
    if (cstr_lag!=NULL) delete cstr_lag;
    if (cstr_rhs!=NULL) delete cstr_rhs;
    cstr_mat = NULL;
    cstr_inv = NULL;
    cstr_lag = NULL;
    cstr_rhs = NULL;
    return;
}

void LocLstSqr::cstr_start(void) {
    cstr_reset();
    use_cstr = true;
    return;
}

void LocLstSqr::cstr_end(void) {
    if (!use_cstr || (Ncstr==0)) {
        use_cstr = false;
        return;
    }

    if (Npars<Ncstr)
        fprintf(stderr,"Warning: number of constraints larger than number of unknowns, the least squares solver will return zero vector.\n");

    return;
}

void LocLstSqr::cstr_solve(void) {
    // Construct cstr_mat
    cstr_mat  = new HMAT(Ncstr);
    *cstr_mat = 0.0;
    for (int n1=0; n1<Ncstr; n1++)
        for (int n2=0; n2<Ncstr; n2++)
            for (int k1=0; k1<Npars; k1++)
                for (int k2=0; k2<Npars; k2++)
                    (*cstr_mat)[n1][n2] += cstr_vars[n1*Npars+k1] * (*coef_inv)[k1][k2] * cstr_vars[n2*Npars+k2];

    // Construct cstr_rhs
    cstr_rhs  = new HVEC(Ncstr);
    *cstr_rhs = 0.0;
    for (int n=0; n<Ncstr; n++) {
        (*cstr_rhs)[n] = -cstr_vals[n];
        for (int k1=0; k1<Npars; k1++)
            for (int k2=0; k2<Npars; k2++)
                (*cstr_rhs)[n] += cstr_vars[n*Npars+k1] * (*coef_inv)[k1][k2] * (*coef_rhs)[k2];
    }

    // Compute cstr_lag
    cstr_inv    = new HMAT(Ncstr);
    (*cstr_inv) = cstr_mat->inv();
    cstr_lag    = new HVEC(Ncstr);
    (*cstr_lag) = (*cstr_inv)*(*cstr_rhs);

    return;
}

// 1D case

void LocLstSqr1D::setup(int order_in, 
        PolyCoef1D* polyCoef_in, bool use_val_in) {
    LocLstSqr::clear();
    LocLstSqr::setOrder(order_in, use_val_in);
    polyCoef = polyCoef_in;
    polyCoef->setup(order_in);
    Npars = polyCoef->getNcoef();
    if (use_val) Npars = Npars-1;

    coef_rhs = new HVEC(Npars); (*coef_rhs) = 0.0;
    coef_par = new HVEC(Npars); (*coef_par) = 0.0;
    coef_mat = new HMAT(Npars); (*coef_mat) = 0.0;

    x_c   = 0.0;
    val_c = 0.0;

    return;
}

void LocLstSqr1D::setCenterV(DType val_c_in) {
    if (!use_val) {
        fprintf(stderr,"Warning: LocLstSqr1D::setCenterV is used but will not be active.\n");
        return;
    }
    val_c = val_c_in;
    polyCoef->setCoef(val_c,0);
    return;
}

void LocLstSqr1D::addNodalMat(
        DType x_n, DType weight) {
    DType dx = x_n - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    int   idx_inc = (use_val) ? 1 : 0;
    for (int i=0; i<Npars; i++)
        coef_int[i] = pow(dx,i+idx_inc);
    add2Mat(coef_int,weight);
    
    if (coef_int!=NULL) delete[] coef_int;
    return;
}

void LocLstSqr1D::addCellMat(
        DType x_l, DType x_r, DType weight) {
    assert(x_l<x_r);
    DType dxl = x_l - x_c;
    DType dxr = x_r - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    int   idx_inc = (use_val) ? 1 : 0;
    for (int i=0; i<Npars; i++)
        coef_int[i] = ((pow(dxr,i+1+idx_inc)-pow(dxl,i+1+idx_inc))/(dxr-dxl))*(1.0/(i+1.0+1.0*idx_inc));
    add2Mat(coef_int,weight);

    if (coef_int!=NULL) delete[] coef_int;
    return;
}

void LocLstSqr1D::addNodalVal(
        DType val, DType x_n, DType weight) {
    DType dx = x_n - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    DType loc_val = (use_val) ? (val - val_c) : val;
    int   idx_inc = (use_val) ? 1 : 0;
    for (int i=0; i<Npars; i++)
        coef_int[i] = pow(dx,i+idx_inc);
    add2Rhs(loc_val,coef_int,weight);
    
    if (coef_int!=NULL) delete[] coef_int;
    return;
}

void LocLstSqr1D::addCellVal(DType val, 
        DType x_l, DType x_r, DType weight) {
    assert(x_l<x_r);
    DType dxl = x_l - x_c;
    DType dxr = x_r - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    DType loc_val = (use_val) ? (val - val_c) : val;
    int   idx_inc = (use_val) ? 1 : 0;
    for (int i=0; i<Npars; i++)
        coef_int[i] = ((pow(dxr,i+1+idx_inc)-pow(dxl,i+1+idx_inc))/(dxr-dxl))*(1.0/(i+1.0+1.0*idx_inc));
    add2Rhs(loc_val,coef_int,weight);

    if (coef_int!=NULL) delete[] coef_int;
    return;
}

void LocLstSqr1D::addNodalArrays(
        DType val, DType x_n, DType weight) {
    DType dx = x_n - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    DType loc_val = (use_val) ? (val - val_c) : val;
    int   idx_inc = (use_val) ? 1 : 0;
    for (int i=0; i<Npars; i++)
        coef_int[i] = pow(dx,i+idx_inc);
    add2Arrays(loc_val,coef_int,weight);
    
    if (coef_int!=NULL) delete[] coef_int;
    return;
}

void LocLstSqr1D::addCellArrays(DType val, 
        DType x_l, DType x_r, DType weight) {
    assert(x_l<x_r);
    DType dxl = x_l - x_c;
    DType dxr = x_r - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    DType loc_val = (use_val) ? (val - val_c) : val;
    int   idx_inc = (use_val) ? 1 : 0;
    for (int i=0; i<Npars; i++)
        coef_int[i] = ((pow(dxr,i+1+idx_inc)-pow(dxl,i+1+idx_inc))/(dxr-dxl))*(1.0/(i+1.0+1.0*idx_inc));
    add2Arrays(loc_val,coef_int,weight);

    if (coef_int!=NULL) delete[] coef_int;
    return;
}

void LocLstSqr1D::addNodalCstr(DType val, DType x_n) {
    DType dx = x_n - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    DType loc_val = (use_val) ? (val - val_c) : val;
    int   idx_inc = (use_val) ? 1 : 0;
    for (int i=0; i<Npars; i++)
        coef_int[i] = pow(dx,i+idx_inc);
    add2CstrLists(loc_val,coef_int);

    Ncstr++;

    if (coef_int!=NULL) delete[] coef_int;
    return;
}

void LocLstSqr1D::addCellCstr(DType val, DType x_l, DType x_r) {
    assert(x_l<x_r);
    DType dxl = x_l - x_c;
    DType dxr = x_r - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    DType loc_val = (use_val) ? (val - val_c) : val;
    int   idx_inc = (use_val) ? 1 : 0;
    for (int i=0; i<Npars; i++)
        coef_int[i] = ((pow(dxr,i+1+idx_inc)-pow(dxl,i+1+idx_inc))/(dxr-dxl))*(1.0/(i+1.0+1.0*idx_inc));
    add2CstrLists(loc_val,coef_int);

    Ncstr++;

    if (coef_int!=NULL) delete[] coef_int;
    return;
}

void LocLstSqr1D::solve(void) {
    // The coef_inv has been reset, compute again 
    // Otherwise, the coef_inv is reused
    if (coef_inv==NULL) {
        coef_inv    = new HMAT(Npars);
        (*coef_inv) = coef_mat->inv();
    }
    if (use_cstr) {
        this->cstr_solve();
        for (int n=0; n<Npars; n++)
            for (int k=0; k<Ncstr; k++)
                (*coef_rhs)[n] -= cstr_vars[k*Npars+n]*(*cstr_lag)[k];
    }
    (*coef_par) = (*coef_inv)*(*coef_rhs);
    if (use_val)
        for (int i=0; i<Npars; i++)
            polyCoef->setGlobCoef((*coef_par)[i],i+1);
    else
        for (int i=0; i<Npars; i++)
            polyCoef->setGlobCoef((*coef_par)[i],i);
    return;
}

// 2D case

void LocLstSqr2D::calcOrd(void) {
    if (use_val)
        for (int i=0; i<Npars; i++)
            polyCoef->globIdx2Der(i+1,o_x[i],o_y[i]);
    else
        for (int i=0; i<Npars; i++)
            polyCoef->globIdx2Der(i,o_x[i],o_y[i]);
    return;
}

void LocLstSqr2D::setup(int order_in,
        PolyCoef2D* polyCoef_in, bool use_val_in) {
    LocLstSqr::clear();
    LocLstSqr::setOrder(order_in, use_val_in);
    if (o_x!=NULL) delete[] o_x;
    if (o_y!=NULL) delete[] o_y;
    polyCoef = polyCoef_in;
    polyCoef->setup(order_in);
    Npars = polyCoef->getNcoef();
    if (use_val) Npars = Npars-1;

    o_x = new int[Npars];
    o_y = new int[Npars];
    calcOrd();

    coef_rhs = new HVEC(Npars); (*coef_rhs) = 0.0;
    coef_par = new HVEC(Npars); (*coef_par) = 0.0;
    coef_mat = new HMAT(Npars); (*coef_mat) = 0.0;

    x_c   = 0.0;
    val_c = 0.0;

    return;
}

void LocLstSqr2D::setCenterV(DType val_c_in) {
    if (!use_val) {
        fprintf(stderr,"Warning: LocLstSqr2D::setCenterV is used but will not be active.\n");
        return;
    }
    val_c = val_c_in;
    polyCoef->setCoef(val_c,0,0);
    return;
}

void LocLstSqr2D::addNodalMat(
        VEC<2> x_n, DType weight) {
    VEC<2> dx = x_n - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];
  
    for (int i=0; i<Npars; i++)
        coef_int[i] = pow(dx[0],o_x[i]) * pow(dx[1],o_y[i]);
    add2Mat(coef_int,weight);
    
    if (coef_int!=NULL) delete[] coef_int;
    return;
}

void LocLstSqr2D::addCellMat(
        VEC<2> x_ll, VEC<2> x_ur, DType weight) {
    VEC<2> dx_ll = x_ll - x_c;
    VEC<2> dx_ur = x_ur - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    for (int i=0; i<Npars; i++)
        coef_int[i] = 
            (((pow(dx_ur[0],o_x[i]+1)-pow(dx_ll[0],o_x[i]+1))/(dx_ur[0]-dx_ll[0]))*(1.0/(o_x[i]+1.0))) *
            (((pow(dx_ur[1],o_y[i]+1)-pow(dx_ll[1],o_y[i]+1))/(dx_ur[1]-dx_ll[1]))*(1.0/(o_y[i]+1.0)));
    add2Mat(coef_int,weight);

    if (coef_int!=NULL) delete[] coef_int;
    return;
}

void LocLstSqr2D::addCellMat(VEC<2> x_1,
        VEC<2> x_2, VEC<2> x_3, DType weight) {
    VEC<2> dx_1 = x_1 - x_c;
    VEC<2> dx_2 = x_2 - x_c;
    VEC<2> dx_3 = x_3 - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    VEC<2> dx12 = dx_2 - dx_1;
    VEC<2> dx13 = dx_3 - dx_1;
    DType area = 0.5*(dx12[0]*dx13[1]-dx12[1]*dx13[0]);
    if (area<0.0) {
        VEC<2> dx_tmp = dx_2;
        dx_2 = dx_3;
        dx_3 = dx_tmp;
        dx12 = dx_2 - dx_1;
        dx13 = dx_3 - dx_1;
        area = -area;
    }

    // right now at most third-order accuracy
    int int_order = 3;
    QuadCoef* quadCoef = newQuadCoef(
            QUADRULE_GAUSSIAN_LOBATTO,
            QUADTYPE_TRIA, int_order);
    if (quadCoef->getNpoly()<polyCoef->getOrder())
        fprintf(stderr,"Warning: LocLstSqr2D::addCellMat: Cell-average evaluation order of accuracy less than requested polynomial order for triangular cell.\n");
    int Nquad = quadCoef->getNquad();

    VEC<2> Xc  = (dx_1+dx_2+dx_3)/3.0;
    VEC<2> dX1 = dx_1 - Xc;
    VEC<2> dX2 = dx_2 - Xc;
    VEC<2> dX3 = dx_3 - Xc;

    for (int i=0; i<Npars; i++) {
        coef_int[i] = 0.0; 
        for (int k=0; k<Nquad; k++) {
            VEC<2> Xloc = Xc 
                + quadCoef->getQuad(k,0) * dX1
                + quadCoef->getQuad(k,1) * dX2
                + quadCoef->getQuad(k,2) * dX3;
            coef_int[i] += quadCoef->getWeight(k) *
                pow(Xloc[0],o_x[i]) * 
                pow(Xloc[1],o_y[i]);
        }
    }
    add2Mat(coef_int,weight);

    if (coef_int!=NULL) delete[] coef_int;
    if (quadCoef!=NULL) delete quadCoef;
    return;
}

void LocLstSqr2D::addNodalVal(
        DType val, VEC<2> x_n, DType weight) {
    VEC<2> dx = x_n - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];
  
    DType loc_val = (use_val) ? (val - val_c) : val;
    for (int i=0; i<Npars; i++)
        coef_int[i] = pow(dx[0],o_x[i]) * pow(dx[1],o_y[i]);
    add2Rhs(loc_val,coef_int,weight);
    
    if (coef_int!=NULL) delete[] coef_int;
    return;
}

void LocLstSqr2D::addCellVal(DType val, 
        VEC<2> x_ll, VEC<2> x_ur, DType weight) {
    VEC<2> dx_ll = x_ll - x_c;
    VEC<2> dx_ur = x_ur - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    DType loc_val = (use_val) ? (val - val_c) : val;
    for (int i=0; i<Npars; i++)
        coef_int[i] = 
            (((pow(dx_ur[0],o_x[i]+1)-pow(dx_ll[0],o_x[i]+1))/(dx_ur[0]-dx_ll[0]))*(1.0/(o_x[i]+1.0))) *
            (((pow(dx_ur[1],o_y[i]+1)-pow(dx_ll[1],o_y[i]+1))/(dx_ur[1]-dx_ll[1]))*(1.0/(o_y[i]+1.0)));
    add2Rhs(loc_val,coef_int,weight);

    if (coef_int!=NULL) delete[] coef_int;
    return;
}

void LocLstSqr2D::addCellVal(DType val, VEC<2> x_1,
        VEC<2> x_2, VEC<2> x_3, DType weight) {
    VEC<2> dx_1 = x_1 - x_c;
    VEC<2> dx_2 = x_2 - x_c;
    VEC<2> dx_3 = x_3 - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    VEC<2> dx12 = dx_2 - dx_1;
    VEC<2> dx13 = dx_3 - dx_1;
    DType area = 0.5*(dx12[0]*dx13[1]-dx12[1]*dx13[0]);
    if (area<0.0) {
        VEC<2> dx_tmp = dx_2;
        dx_2 = dx_3;
        dx_3 = dx_tmp;
        dx12 = dx_2 - dx_1;
        dx13 = dx_3 - dx_1;
        area = -area;
    }

    // right now at most third-order accuracy
    int int_order = 3;
    QuadCoef* quadCoef = newQuadCoef(
            QUADRULE_GAUSSIAN_LOBATTO,
            QUADTYPE_TRIA, int_order);
    if (quadCoef->getNpoly()<polyCoef->getOrder())
        fprintf(stderr,"Warning: LocLstSqr2D::addCellVal: Cell-average evaluation order of accuracy less than requested polynomial order for triangular cell.\n");
    int Nquad = quadCoef->getNquad();

    VEC<2> Xc  = (dx_1+dx_2+dx_3)/3.0;
    VEC<2> dX1 = dx_1 - Xc;
    VEC<2> dX2 = dx_2 - Xc;
    VEC<2> dX3 = dx_3 - Xc;

    DType loc_val = (use_val) ? (val-val_c) : val;
    for (int i=0; i<Npars; i++) {
        coef_int[i] = 0.0; 
        for (int k=0; k<Nquad; k++) {
            VEC<2> Xloc = Xc 
                + quadCoef->getQuad(k,0) * dX1
                + quadCoef->getQuad(k,1) * dX2
                + quadCoef->getQuad(k,2) * dX3;
            coef_int[i] += quadCoef->getWeight(k) *
                pow(Xloc[0],o_x[i]) * 
                pow(Xloc[1],o_y[i]);
        }
    }
    add2Rhs(loc_val,coef_int,weight);

    if (coef_int!=NULL) delete[] coef_int;
    if (quadCoef!=NULL) delete quadCoef;
    return;
}

void LocLstSqr2D::addNodalArrays(
        DType val, VEC<2> x_n, DType weight) {
    VEC<2> dx = x_n - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];
  
    DType loc_val = (use_val) ? (val - val_c) : val;
    for (int i=0; i<Npars; i++)
        coef_int[i] = pow(dx[0],o_x[i]) * pow(dx[1],o_y[i]);
    add2Arrays(loc_val,coef_int,weight);
    
    if (coef_int!=NULL) delete[] coef_int;
    return;
}

void LocLstSqr2D::addCellArrays(DType val, 
        VEC<2> x_ll, VEC<2> x_ur, DType weight) {
    VEC<2> dx_ll = x_ll - x_c;
    VEC<2> dx_ur = x_ur - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    DType loc_val = (use_val) ? (val - val_c) : val;
    for (int i=0; i<Npars; i++)
        coef_int[i] = 
            (((pow(dx_ur[0],o_x[i]+1)-pow(dx_ll[0],o_x[i]+1))/(dx_ur[0]-dx_ll[0]))*(1.0/(o_x[i]+1.0))) *
            (((pow(dx_ur[1],o_y[i]+1)-pow(dx_ll[1],o_y[i]+1))/(dx_ur[1]-dx_ll[1]))*(1.0/(o_y[i]+1.0)));
    add2Arrays(loc_val,coef_int,weight);

    if (coef_int!=NULL) delete[] coef_int;
    return;
}

void LocLstSqr2D::addCellArrays(DType val, VEC<2> x_1,
        VEC<2> x_2, VEC<2> x_3, DType weight) {
    VEC<2> dx_1 = x_1 - x_c;
    VEC<2> dx_2 = x_2 - x_c;
    VEC<2> dx_3 = x_3 - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    VEC<2> dx12 = dx_2 - dx_1;
    VEC<2> dx13 = dx_3 - dx_1;
    DType area = 0.5*(dx12[0]*dx13[1]-dx12[1]*dx13[0]);
    if (area<0.0) {
        VEC<2> dx_tmp = dx_2;
        dx_2 = dx_3;
        dx_3 = dx_tmp;
        dx12 = dx_2 - dx_1;
        dx13 = dx_3 - dx_1;
        area = -area;
    }

    // right now at most third-order accuracy
    int int_order = 3;
    QuadCoef* quadCoef = newQuadCoef(
            QUADRULE_GAUSSIAN_LOBATTO,
            QUADTYPE_TRIA, int_order);
    if (quadCoef->getNpoly()<polyCoef->getOrder())
        fprintf(stderr,"Warning: LocLstSqr2D::addCellArrays: Cell-average evaluation order of accuracy less than requested polynomial order for triangular cell.\n");
    int Nquad = quadCoef->getNquad();

    VEC<2> Xc  = (dx_1+dx_2+dx_3)/3.0;
    VEC<2> dX1 = dx_1 - Xc;
    VEC<2> dX2 = dx_2 - Xc;
    VEC<2> dX3 = dx_3 - Xc;

    DType loc_val = (use_val) ? (val-val_c) : val;
    for (int i=0; i<Npars; i++) {
        coef_int[i] = 0.0; 
        for (int k=0; k<Nquad; k++) {
            VEC<2> Xloc = Xc 
                + quadCoef->getQuad(k,0) * dX1
                + quadCoef->getQuad(k,1) * dX2
                + quadCoef->getQuad(k,2) * dX3;
            coef_int[i] += quadCoef->getWeight(k) *
                pow(Xloc[0],o_x[i]) * 
                pow(Xloc[1],o_y[i]);
        }
    }
    add2Arrays(loc_val,coef_int,weight);

    if (coef_int!=NULL) delete[] coef_int;
    if (quadCoef!=NULL) delete quadCoef;
    return;
}

void LocLstSqr2D::addNodalCstr(DType val, VEC<2> x_n) {
    VEC<2> dx = x_n - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    DType loc_val = (use_val) ? (val - val_c) : val;
    for (int i=0; i<Npars; i++)
        coef_int[i] = pow(dx[0],o_x[i]) * pow(dx[1],o_y[i]);
    add2CstrLists(loc_val,coef_int);
    
    Ncstr++;

    if (coef_int!=NULL) delete[] coef_int;
    return;
}

void LocLstSqr2D::addCellCstr(DType val, 
        VEC<2> x_ll, VEC<2> x_ur) {
    VEC<2> dx_ll = x_ll - x_c;
    VEC<2> dx_ur = x_ur - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    DType loc_val = (use_val) ? (val - val_c) : val;
    for (int i=0; i<Npars; i++)
        coef_int[i] = 
            (((pow(dx_ur[0],o_x[i]+1)-pow(dx_ll[0],o_x[i]+1))/(dx_ur[0]-dx_ll[0]))*(1.0/(o_x[i]+1.0))) *
            (((pow(dx_ur[1],o_y[i]+1)-pow(dx_ll[1],o_y[i]+1))/(dx_ur[1]-dx_ll[1]))*(1.0/(o_y[i]+1.0)));
    add2CstrLists(loc_val,coef_int);

    Ncstr++;

    if (coef_int!=NULL) delete[] coef_int;
    return;
}

void LocLstSqr2D::addCellCstr(DType val, VEC<2> x_1,
        VEC<2> x_2, VEC<2> x_3) {
    VEC<2> dx_1 = x_1 - x_c;
    VEC<2> dx_2 = x_2 - x_c;
    VEC<2> dx_3 = x_3 - x_c;
    DType *coef_int;
    coef_int = new DType[Npars];

    VEC<2> dx12 = dx_2 - dx_1;
    VEC<2> dx13 = dx_3 - dx_1;
    DType area = 0.5*(dx12[0]*dx13[1]-dx12[1]*dx13[0]);
    if (area<0.0) {
        VEC<2> dx_tmp = dx_2;
        dx_2 = dx_3;
        dx_3 = dx_tmp;
        dx12 = dx_2 - dx_1;
        dx13 = dx_3 - dx_1;
        area = -area;
    }

    // right now at most third-order accuracy
    int int_order = 3;
    QuadCoef* quadCoef = newQuadCoef(
            QUADRULE_GAUSSIAN_LOBATTO,
            QUADTYPE_TRIA, int_order);
    if (quadCoef->getNpoly()<polyCoef->getOrder())
        fprintf(stderr,"Warning: LocLstSqr2D::addCellArrays: Cell-average evaluation order of accuracy less than requested polynomial order for triangular cell.\n");
    int Nquad = quadCoef->getNquad();

    VEC<2> Xc  = (dx_1+dx_2+dx_3)/3.0;
    VEC<2> dX1 = dx_1 - Xc;
    VEC<2> dX2 = dx_2 - Xc;
    VEC<2> dX3 = dx_3 - Xc;

    DType loc_val = (use_val) ? (val-val_c) : val;
    for (int i=0; i<Npars; i++) {
        coef_int[i] = 0.0; 
        for (int k=0; k<Nquad; k++) {
            VEC<2> Xloc = Xc 
                + quadCoef->getQuad(k,0) * dX1
                + quadCoef->getQuad(k,1) * dX2
                + quadCoef->getQuad(k,2) * dX3;
            coef_int[i] += quadCoef->getWeight(k) *
                pow(Xloc[0],o_x[i]) * 
                pow(Xloc[1],o_y[i]);
        }
    }
    add2CstrLists(loc_val,coef_int);

    Ncstr++;

    if (coef_int!=NULL) delete[] coef_int;
    if (quadCoef!=NULL) delete quadCoef;
    return;
}

void LocLstSqr2D::solve(void) {
    // The coef_inv has been reset, compute again 
    // Otherwise, the coef_inv is reused
    if (coef_inv==NULL) {
        coef_inv    = new HMAT(Npars);
        (*coef_inv) = coef_mat->inv();
    }
    if (use_cstr) {
        this->cstr_solve();
        for (int n=0; n<Npars; n++)
            for (int k=0; k<Ncstr; k++)
                (*coef_rhs)[n] -= cstr_vars[k*Npars+n]*(*cstr_lag)[k];
    }
    (*coef_par) = (*coef_inv)*(*coef_rhs);
    if (use_val)
        for (int i=0; i<Npars; i++)
            polyCoef->setGlobCoef((*coef_par)[i],i+1);
    else
        for (int i=0; i<Npars; i++)
            polyCoef->setGlobCoef((*coef_par)[i],i);
    return;
}

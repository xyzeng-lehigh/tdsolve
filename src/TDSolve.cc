#include "TDSolve.h"
#include "Alloc_Utils.h"
#include "TDInput.h"
#include "Alg_Solve.h"
using namespace std;

TDSolve::TDSolve(TDInput* tdInput_ptr) {
  tdInput = tdInput_ptr;
  tdCoef  = &(tdInput->tdCoef);
  tdMesh  = &(tdInput->tdMesh);
  tdSetup = &(tdInput->tdSetup);
  tdLin   = &(tdInput->tdLin);
  Nx      = tdMesh->Nx;
  dx      = tdMesh->dx;
  x_f     = tdMesh->x_f;
  x_c     = tdMesh->x_c;
  Na_m    = tdMesh->Na_m;
  Na      = tdMesh->Na;
  x_a_f   = tdMesh->x_a_f;
  x_a_c   = tdMesh->x_a_c;

  dtn     = tdSetup->dTinit;
  dtn_old = tdSetup->dTinit;
  tsol    = tdSetup->Tinit;

  reach_term_now  = tdSetup->Tterm <= tsol;
  reach_term_next = false;

  mat = NULL;
  prc = NULL;
  sol = NULL;
  rhs = NULL;
  res = NULL;
  linSolve = NULL;

  mat_a = NULL;
  prc_a = NULL;
  sol_a = NULL;
  rhs_a = NULL;
  res_a = NULL;
  linSolve_a = NULL;
  iter = 0;

  // OPN Model
  mat_p = NULL;
  prc_p = NULL;
  sol_p = NULL;
  rhs_p = NULL;
  res_p = NULL;
  linSolve_p = NULL;

  use_tfreq = tdSetup->Tfreq > 0.0;
  reach_output_now  = false;
  reach_output_next = false;
  if ( use_tfreq ) 
    t_output = tdSetup->Tfreq;
}

TDSolve::~TDSolve(void) {
  erase_array( data_g, 3 );
  erase_array( data_h, 3 );
  erase_array( data_n, 3 );
  erase_array( data_a, 3 );
  erase_array( data_v, 3 );
  // OPN Model
  erase_array(data_p1, 3);
  erase_array(data_p2, 3);

  if (mat!=NULL) delete mat;
  if (prc!=NULL) delete prc;
  if (rhs!=NULL) delete rhs;
  if (res!=NULL) delete res;
  if (sol!=NULL) delete sol;
  if (linSolve!=NULL) delete linSolve;

  if (mat_a!=NULL) delete mat_a;
  if (prc_a!=NULL) delete prc_a;
  if (rhs_a!=NULL) delete rhs_a;
  if (res_a!=NULL) delete res_a;
  if (sol_a!=NULL) delete sol_a;
  if (linSolve_a!=NULL) delete linSolve_a;

  // OPN Model
  if (mat_p!=NULL) delete mat_p;
  if (prc_p!=NULL) delete prc_p;
  if (rhs_p!=NULL) delete rhs_p;
  if (res_p!=NULL) delete res_p;
  if (sol_p!=NULL) delete sol_p;
  if (linSolve_p!=NULL) delete linSolve_p;
}


void TDSolve::setupLinSolve( 
    LinMat*& mat_in, LinMat*& prc_in, LinVec*& sol_in, LinVec*& rhs_in, 
    LinVec*& res_in, LinSolve*& linSolve_in, int Nsize ) {
  mat_in = new LinSqSparseMat;
  mat_in->setBlckDim(1);
  mat_in->setEqualSize(Nsize);

  prc_in = new LinSqSparseMat;
  prc_in->setBlckDim(1);
  prc_in->setEqualSize(Nsize);

  sol_in = new LinVec;
  sol_in->setBlckDim(1);
  sol_in->setSize(Nsize);

  rhs_in = new LinVec;
  rhs_in->setBlckDim(1);
  rhs_in->setSize(Nsize);

  res_in = new LinVec;
  res_in->setBlckDim(1);
  res_in->setSize(Nsize);

  linSolve_in = new LinSolve;
  linSolve_in->setBlckDim(1);
  linSolve_in->setSize(Nsize);

  linSolve_in->linkMat(mat_in);
  linSolve_in->linkPrc(prc_in);
  linSolve_in->linkSol(sol_in);
  linSolve_in->linkRhs(rhs_in);
  linSolve_in->linkRes(res_in);

  linSolve_in->setPrintStatus(tdLin->Print);
  linSolve_in->setSolvePrec(tdLin->Prec,tdLin->PrecPar);
  switch (tdLin->Type) {
    case LSYS_TYPE_DIR:
      linSolve_in->setSolveDir();
      break;
    case LSYS_TYPE_ITER:
      linSolve_in->setMaxIter(tdLin->MaxIter);
      linSolve_in->setTargetRes(tdLin->TargetRes);
      linSolve_in->setSolveIter(tdLin->Iter);
      break;
    case LSYS_TYPE_KRYLOV:
      linSolve_in->setMaxIter(tdLin->MaxIter);
      linSolve_in->setTargetRes(tdLin->TargetRes);
      linSolve_in->setSolveKrylov(tdLin->Krylov,tdLin->KrylovCap);
      break;
  }

  return;
}

void TDSolve::average_sol(const DType coef_1, const int fr_1,
    const DType coef_2, const int fr_2, const int to, 
    const DType shift_dt_to) {
  
  if ( tdSetup->frame == 2 && tdSetup->scell >= 0 )
  {
    downscale_var( *data_g[fr_1], *data_g[fr_1], x_c, Nx, data_r[fr_1] );
    downscale_var( *data_h[fr_1], *data_h[fr_1], x_c, Nx, data_r[fr_1] );
    downscale_var( *data_n[fr_1], *data_n[fr_1], x_c, Nx, data_r[fr_1] );

    downscale_var( *data_g[fr_2], *data_g[fr_2], x_c, Nx, data_r[fr_2] );
    downscale_var( *data_h[fr_2], *data_h[fr_2], x_c, Nx, data_r[fr_2] );
    downscale_var( *data_n[fr_2], *data_n[fr_2], x_c, Nx, data_r[fr_2] );
  }

  *data_g[to] = coef_1 * (*data_g[fr_1]) + coef_2 * (*data_g[fr_2]);
  *data_h[to] = coef_1 * (*data_h[fr_1]) + coef_2 * (*data_h[fr_2]);
  *data_n[to] = coef_1 * (*data_n[fr_1]) + coef_2 * (*data_n[fr_2]);
  *data_a[to] = coef_1 * (*data_a[fr_1]) + coef_2 * (*data_a[fr_2]);
  
  // OPN model -- only implemented for FVM or FVM_Mod
  if ( tdSetup->prob_id == 1 && 
       ( tdSetup->frame == 1 || tdSetup->frame == 2 ) ) {
    *data_p1[to] = coef_1 * (*data_p1[fr_1]) + coef_2 * (*data_p1[fr_2]);
    *data_p2[to] = coef_1 * (*data_p2[fr_1]) + coef_2 * (*data_p2[fr_2]);
  }

  data_r[to]  = coef_1 * data_r[fr_1]    + coef_2 * data_r[fr_2];

  if ( tdSetup->frame == 2 && tdSetup->scell >= 0 )
  {
    upscale_var( *data_g[fr_1], *data_g[fr_1], x_c, Nx, data_r[fr_1] );
    upscale_var( *data_h[fr_1], *data_h[fr_1], x_c, Nx, data_r[fr_1] );
    upscale_var( *data_n[fr_1], *data_n[fr_1], x_c, Nx, data_r[fr_1] );

    upscale_var( *data_g[fr_2], *data_g[fr_2], x_c, Nx, data_r[fr_2] );
    upscale_var( *data_h[fr_2], *data_h[fr_2], x_c, Nx, data_r[fr_2] );
    upscale_var( *data_n[fr_2], *data_n[fr_2], x_c, Nx, data_r[fr_2] );

    upscale_var( *data_g[to], *data_g[to], x_c, Nx, data_r[to] );
    upscale_var( *data_h[to], *data_h[to], x_c, Nx, data_r[to] );
    upscale_var( *data_n[to], *data_n[to], x_c, Nx, data_r[to] );
  }

  if ( tdSetup->re_vel ) // Do not use this option for FVM_Mod
    this->update_vel(to, 0, shift_dt_to, true);
  else
    *data_v[to] = coef_1 * (*data_v[fr_1]) + coef_2 * (*data_v[fr_2]);

  return;
}

void TDSolve::diagnostic(const int nt) {
  fprintf( stdout, "Iter %3d: Time %6.8f; Time step %6.8f; Radius %6.8f; Front-velocity %6.8f.\n",
      iter, tsol, dtn, data_r[nt], rprime );
  return;
}

void TDSolve::update_rad(const int fr, const int to, const int v_nt) {
  int vnt = (v_nt>=0) ? v_nt : fr;
  rprime     = (*data_v[vnt])[Nx];           // dr_opt == 0
  data_r[to] = data_r[fr] + dtn * rprime;   //  r_opt == 0
  if ( tdSetup->vCoup ) {
    if ( tdSetup->dr_opt == 1 )
      rprime = rprime * ( 1.0 - 0.5 * dtn * rprime / data_r[to] );
    if ( tdSetup->dr_opt == 2 )
      rprime = rprime / ( 1.0 - 0.25 * dx * dx );

    if ( tdSetup->r_opt == 1 )
      data_r[to] = std::sqrt( data_r[fr]*data_r[fr] + dtn*dtn*rprime*rprime ) + dtn*rprime;
  }
  else {
    if ( tdSetup->dr_opt == 1 )
      rprime = rprime * ( 1.0 + 0.5 * dtn * rprime / data_r[fr] );
    if ( tdSetup->dr_opt == 2 )
      rprime = rprime / ( 1.0 - 0.25 * dx * dx );

    if ( tdSetup->r_opt == 1 )
      data_r[to] = std::sqrt( data_r[fr]*data_r[fr] + 2.0*dtn*data_r[fr]*rprime );
  }
  return;
}

void TDSolve::upscale_var(
    HVEC& data_fr, HVEC& data_to, DType*& x_ptr, 
    int Nsize, DType radius, int scal_pow) {
  for ( int i = 0; i < Nsize; i++ ) {
    DType loc_r = x_ptr[i] * radius;
    data_to[i] = data_fr[i] * pow(loc_r, scal_pow);
  }
  return;
}

void TDSolve::downscale_var(
    HVEC& data_fr, HVEC& data_to, DType*& x_ptr, 
    int Nsize, DType radius, int scal_pow) {
  for ( int i = 0; i < Nsize; i++ ) {
    DType loc_r = x_ptr[i] * radius;
    data_to[i] = data_fr[i] / pow(loc_r, scal_pow);
  }
  return;
}


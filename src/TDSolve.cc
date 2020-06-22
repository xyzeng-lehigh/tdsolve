#include "TDSolve.h"
#include "Alloc_Utils.h"
#include "TDInput.h"
#include "Alg_Solve.h"
#include "RK_Coef.h"
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

  // Radau RK
  matBlk = NULL;
  matFlt = NULL;
  matSlm = NULL;
  matRK  = NULL;
  prcRK  = NULL;
  solRK  = NULL;
  rhsRK  = NULL;
  resRK  = NULL;
  linSolveRK = NULL;

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

  data_g = NULL;
  data_h = NULL;
  data_n = NULL;
  data_a = NULL;
  data_v = NULL;

  data_r  = NULL;
  data_rp = NULL;
  Ndata   = 3;

  data_cell = NULL;
  data_node = NULL;

  diff_g  = NULL;
  diff_h  = NULL;
  diff_n  = NULL;
  diff_v  = NULL;
  diff_r  = NULL;

  data_p1 = NULL;
  data_p2 = NULL;
  data_dv = NULL;
  data_na = NULL;

  resd_g = NULL;
  resd_h = NULL;
  resd_n = NULL;
  resd_a = NULL;

  resd_p1 = NULL;
  resd_p2 = NULL;

  use_resd = false;
}

TDSolve::~TDSolve(void) {
  erase_array( data_g, Ndata );
  erase_array( data_h, Ndata );
  erase_array( data_n, Ndata );
  erase_array( data_a, Ndata );
  erase_array( data_v, Ndata );

  erase_array( data_cell, 1 );
  erase_array( data_node, 1 );
  // OPN Model
  erase_array( data_p1, Ndata );
  erase_array( data_p2, Ndata );

  erase_array( resd_g, Ndata );
  erase_array( resd_h, Ndata );
  erase_array( resd_n, Ndata );
  erase_array( resd_a, Ndata );
  // OPN Model
  erase_array( resd_p1, Ndata );
  erase_array( resd_p2, Ndata );

  if (data_r!=NULL)  delete[] data_r;
  if (data_rp!=NULL) delete[] data_rp;

  if (diff_g!=NULL)  delete[] diff_g;
  if (diff_h!=NULL)  delete[] diff_h;
  if (diff_n!=NULL)  delete[] diff_n;
  if (diff_v!=NULL)  delete[] diff_v;
  if (diff_r!=NULL)  delete[] diff_r;

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
  
  // Radau RK
  if (matBlk!=NULL) delete matBlk;
  if (matFlt!=NULL) delete matFlt;
  if (matSlm!=NULL) delete matSlm;
  if (matRK!=NULL)  delete matRK;
  if (prcRK!=NULL)  delete prcRK;
  if (rhsRK!=NULL)  delete rhsRK;
  if (resRK!=NULL)  delete resRK;
  if (solRK!=NULL)  delete solRK;
  if (linSolveRK!=NULL) delete linSolveRK;

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

    downscale_var( *data_a[fr_1], *data_a[fr_1], x_a_c, Na, data_r[fr_1] );
    downscale_var( *data_a[fr_2], *data_a[fr_2], x_a_c, Na, data_r[fr_2] );
  }

  *data_g[to] = coef_1 * (*data_g[fr_1]) + coef_2 * (*data_g[fr_2]);
  *data_h[to] = coef_1 * (*data_h[fr_1]) + coef_2 * (*data_h[fr_2]);
  *data_n[to] = coef_1 * (*data_n[fr_1]) + coef_2 * (*data_n[fr_2]);
  *data_a[to] = coef_1 * (*data_a[fr_1]) + coef_2 * (*data_a[fr_2]);
  
  // OPN model -- only implemented for FVM or FVM_Mod
  if ( tdSetup->prob_id == 1 && 
       ( tdSetup->frame == 1 || tdSetup->frame == 2 ) ) {
    downscale_var( *data_p1[fr_1], *data_p1[fr_1], x_c, Nx, data_r[fr_1] );
    downscale_var( *data_p1[fr_2], *data_p1[fr_2], x_c, Nx, data_r[fr_2] );
    downscale_var( *data_p2[fr_1], *data_p2[fr_1], x_c, Nx, data_r[fr_1] );
    downscale_var( *data_p2[fr_2], *data_p2[fr_2], x_c, Nx, data_r[fr_2] );

    *data_p1[to] = coef_1 * (*data_p1[fr_1]) + coef_2 * (*data_p1[fr_2]);
    *data_p2[to] = coef_1 * (*data_p2[fr_1]) + coef_2 * (*data_p2[fr_2]);

    upscale_var( *data_p1[fr_1], *data_p1[fr_1], x_c, Nx, data_r[fr_1] );
    upscale_var( *data_p1[fr_2], *data_p1[fr_2], x_c, Nx, data_r[fr_2] );
    upscale_var( *data_p1[to],   *data_p1[to],   x_c, Nx, data_r[to]   );
    upscale_var( *data_p2[fr_1], *data_p2[fr_1], x_c, Nx, data_r[fr_1] );
    upscale_var( *data_p2[fr_2], *data_p2[fr_2], x_c, Nx, data_r[fr_2] );
    upscale_var( *data_p2[to],   *data_p2[to],   x_c, Nx, data_r[to]   );
  }

  data_r[to]  = coef_1 * data_r[fr_1] + coef_2 * data_r[fr_2];

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

    upscale_var( *data_a[fr_1], *data_a[fr_1], x_a_c, Na, data_r[fr_1] );
    upscale_var( *data_a[fr_2], *data_a[fr_2], x_a_c, Na, data_r[fr_2] );
    upscale_var( *data_a[to],   *data_a[to],   x_a_c, Na, data_r[to]   );
  }

  *data_v[to] = coef_1 * (*data_v[fr_1]) + coef_2 * (*data_v[fr_2]);

  return;
}

void TDSolve::diagnostic(const int nt) {
  if ( tdSetup->silence_screen ) return;
  fprintf( stdout, "Iter %3d: Time %6.8f; Time step %6.8f; CFL %2.4f; Radius %6.8f; Front-velocity %6.8f.\n",
      iter, tsol, dtn, cur_cfl, data_r[nt], data_rp[0] );
  return;
}

void TDSolve::update_rad(const int fr, const int to, const int v_nt) {
  int vnt = (v_nt>=0) ? v_nt : fr;
  data_rp[0] = (*data_v[vnt])[Nx];            // dr_opt == 0
  data_r[to] = data_r[fr] + dtn * data_rp[0]; //  r_opt == 0
  if ( tdSetup->vCoup ) {
    if ( tdSetup->dr_opt == 1 )
      data_rp[0] = data_rp[0] * ( 1.0 - 0.5 * dtn * data_rp[0] / data_r[to] );
    if ( tdSetup->dr_opt == 2 )
      data_rp[0] = data_rp[0] / ( 1.0 - 0.25 * dx * dx );

    if ( tdSetup->r_opt == 1 )
      data_r[to] = std::sqrt( data_r[fr]*data_r[fr] + dtn*dtn*data_rp[0]*data_rp[0]) + dtn*data_rp[0];
  }
  else {
    if ( tdSetup->dr_opt == 1 )
      data_rp[0] = data_rp[0] * ( 1.0 + 0.5 * dtn * data_rp[0] / data_r[fr] );
    if ( tdSetup->dr_opt == 2 )
      data_rp[0] = data_rp[0] / ( 1.0 - 0.25 * dx * dx );

    if ( tdSetup->r_opt == 1 )
      data_r[to] = std::sqrt( data_r[fr]*data_r[fr] + 2.0*dtn*data_r[fr]*data_rp[0] );
  }
  return;
}

void TDSolve::update_rad_radau(RKButcher& rkButcher, const bool diagnos) {
  int Ns = rkButcher.getNstage();
  for ( int i = 1; i <= Ns; i++ )
    downscale_sol(i,i); // Get the primitive Stage i (0)

  // It is supposed that data_r[1] ~ data_r[Ns] are allocated for
  // the 1st ~ Ns-th stages, respectively.

  // Step 1 - compute the tumor velocity for each stage
  DType discount = 1.0/(1.0-0.25*dx*dx);
  for ( int i = 1; i <= Ns; i++ )
    data_rp[i] = discount * (*data_v[i])[Nx];

  // Step 2 - solve the nonlinear system for data_r:
  // r[i]^2 - data_r[0]^2 = 2.0*dtn*sum_{j=1~Ns}a_{ij}rp[j]*r[j]
  // i = 1~Ns
  HVEC sol(Ns,0.0);  // Solution iterates
  for ( int i = 0; i < Ns; i++ )
    sol[i] = data_r[i+1];
  HVEC fval(Ns,0.0); // Function evaluates
  HMAT hess(Ns,0.0); // Hessian
  // Newton-Raphson
  int c_iter = 0;
  eval_rad_radau_fval(fval, rkButcher, data_rp+1, sol);
  while ( fval.norm(-1) > 1e-12 ) {
    if (diagnos) {
      cout << "TDSolve::update_rad_radau -- iteration " << c_iter << std::endl;
      cout << "\t"; sol.print();
      cout << "\tResidual: " << fval.norm(-1) << std::endl;
    }
    eval_rad_radau_hess(hess, rkButcher, data_rp+1, sol);
    sol -= hess.inv() * fval;
    eval_rad_radau_fval(fval, rkButcher, data_rp+1, sol);
    c_iter++;
  }

  // Step 3 - assign the converged solution to data_r
  for ( int i = 0; i < Ns; i++ ) {
    // The difference between new and old data_r
    diff_r[i+1] = std::fabs(data_r[i+1]-sol[i]); 
    data_r[i+1] = sol[i];
  }

  for ( int i = 1; i <= Ns; i++ )
    upscale_sol(i,i); // Get the conservative Stage i using the re-calculated data_r[i]

  return;
}

void TDSolve::eval_rad_radau_fval(HVEC& fval,
    const RKButcher& rkButcher, const DType* rp, const HVEC& sol) {
  int Ns = rkButcher.getNstage();
  fval = 0.0;
  for ( int i = 0; i < Ns; i++ ) {
    for ( int j = 0; j < Ns; j++ )
      fval[i] -= rkButcher.get(i,j) * rp[j] * sol[j];
    fval[i] *= 2.0 * dtn;
    fval[i] += sol[i] * sol[i] - data_r[0] * data_r[0];
  }
  return;
}

void TDSolve::eval_rad_radau_hess(HMAT& hess,
    const RKButcher& rkButcher, const DType* rp, const HVEC& sol) {
  int Ns = rkButcher.getNstage();
  hess = 0.0;
  for ( int i = 0; i < Ns; i++ ) 
    for ( int j = 0; j < Ns; j++ )
      hess[i][j] -= rkButcher.get(i,j) * rp[j];
  hess *= 2.0 * dtn;
  for ( int i = 0; i < Ns; i++ )
    hess[i][i] += 2.0 * sol[i];
  return;
}

void TDSolve::check_conv_radau(const RKButcher& rkButcher, int& n_iter, 
    bool& sol_conv, bool& red_dt, const bool diagnos) {
  DType res_s = 0.0;
  int Ns = rkButcher.getNstage();
  for ( int i = 0; i < Ns; i++ ) {
    // data_r
    res_s = std::max( res_s, diff_r[i+1] );
    // data_v
    res_s = std::max( res_s, diff_v[i+1] );
    if (!tdSetup->vOnly) {
      // data_g
      res_s = std::max( res_s, diff_g[i+1] );
      // data_h
      res_s = std::max( res_s, diff_h[i+1] );
      // data_n
      res_s = std::max( res_s, diff_n[i+1] );
    }
  }
  n_iter++;
  if (n_iter == 0) init_err = 1.0;   // In the first iteration, use absolute error
  sol_conv = res_s <= (tdSetup->vTol * init_err);
  red_dt   = (tdSetup->vIter >= 0) && (n_iter == tdSetup->vIter) && (res_s > (tdSetup->vTol * init_err) );
  if (n_iter == 0) init_err = res_s; // In the next iterations, use relative error
  if (diagnos) {
    std::cout << "Radau PMC iteration: " << n_iter << "; Maximum residue: " << res_s << std::endl;
    std::cout << "  radius: ";
    for ( int i = 0; i < Ns; i++ )
      std::cout << "stage " << i+1 << " - " << std::fabs( diff_r[i+1] ) << "; ";
    std::cout << std::endl;
    std::cout << "  velocity: ";
    for ( int i = 0; i < Ns; i++ )
      std::cout << "stage " << i+1 << " - " << std::fabs( diff_v[i+1] ) << "; ";
    std::cout << std::endl;
    std::cout << "  G: ";
    for ( int i = 0; i < Ns; i++ )
      std::cout << "stage " << i+1 << " - " << std::fabs( diff_g[i+1] ) << "; ";
    std::cout << std::endl;
    std::cout << "  H: ";
    for ( int i = 0; i < Ns; i++ )
      std::cout << "stage " << i+1 << " - " << std::fabs( diff_h[i+1] ) << "; ";
    std::cout << std::endl;
    std::cout << "  N: ";
    for ( int i = 0; i < Ns; i++ )
      std::cout << "stage " << i+1 << " - " << std::fabs( diff_n[i+1] ) << "; ";
    std::cout << std::endl;
  }
  if (sol_conv)
    std::cout << "Radau PMC converged in iteration: " << n_iter << "; Residue: " << res_s << std::endl;
  if (red_dt)
    std::cout << "Radau PMC iteration: " << n_iter << "; Maximum residue: " << res_s << "; Will reduce dt." << std::endl;
  return;
}

void TDSolve::check_conv_radius(const RKButcher& rkButcher, int& n_iter, 
    bool& sol_conv, bool& red_dt, const bool diagnos) {
  DType res_s = 0.0;
  int Ns = rkButcher.getNstage();
  for ( int i = 0; i < Ns; i++ ) 
    res_s = std::max( res_s, diff_r[i+1] );
  n_iter++;
  if (n_iter == 0) init_err = 1.0;   // In the first iteration, use absolute error
  sol_conv = res_s <= (tdSetup->vTol * init_err);
  red_dt   = (tdSetup->vIter >= 0) && (n_iter == tdSetup->vIter) && (res_s > (tdSetup->vTol * init_err) );
  if (n_iter == 0) init_err = res_s; // In the next iterations, use relative error
  if (diagnos) {
    std::cout << "Radau PMC iteration: " << n_iter << "; Maximum diff in R: " << res_s << std::endl;
    std::cout << "  radius: ";
    for ( int i = 0; i < Ns; i++ )
      std::cout << "stage " << i+1 << " - " << std::fabs( diff_r[i+1] ) << "; ";
    std::cout << std::endl;
  }
  if (sol_conv)
    std::cout << "Radau PMC converged in iteration: " << n_iter << "; diff in R: " << res_s << std::endl;
  if (red_dt)
    std::cout << "Radau PMC iteration: " << n_iter << "; Maximum diff in R: " << res_s << "; Will reduce dt." << std::endl;
  return;
}

void TDSolve::reduce_dtn(const DType scale) {
  scale_dtn( scale );
  dtn_old = dtn;
  cur_cfl *= scale;
  reach_term_next = false;
  reach_output_next = false;
}

bool TDSolve::check_Rterm(const int nt) {
  return (data_r[nt] >= tdSetup->Rterm-tdSetup->vTol);
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

void TDSolve::downscale_sol(const int s_nt, const int r_nt) {
  downscale_var( *data_g[s_nt], *data_g[s_nt], x_c, Nx, data_r[r_nt] );
  downscale_var( *data_h[s_nt], *data_h[s_nt], x_c, Nx, data_r[r_nt] );
  downscale_var( *data_n[s_nt], *data_n[s_nt], x_c, Nx, data_r[r_nt] );
}

void TDSolve::upscale_sol(const int s_nt, const int r_nt) {
  upscale_var( *data_g[s_nt], *data_g[s_nt], x_c, Nx, data_r[r_nt] );
  upscale_var( *data_h[s_nt], *data_h[s_nt], x_c, Nx, data_r[r_nt] );
  upscale_var( *data_n[s_nt], *data_n[s_nt], x_c, Nx, data_r[r_nt] );
}

void TDSolve::print_warning(const char* w_message) {
  std::cerr << "Warning: " << w_message << std::endl;
}

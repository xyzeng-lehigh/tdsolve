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

  dtn     = 1.0;
  tsol    = tdSetup->Tinit;

  reach_term = tdSetup->Tterm <= tsol;

  alloc_array_full( data_g,  Nx+1, 3 );
  alloc_array_full( data_h,  Nx+1, 3 );
  alloc_array_full( data_n,  Nx+1, 3 );
  alloc_array_full( data_a,  Nx+1, 3 );
  alloc_array_full( data_v,  Nx+1, 3 );
  alloc_array_full( data_na, Nx+1, 3 );

  FILE* ot_file = fopen( tdInput->outputFile, "a" );
  fprintf( ot_file, "Time, time step size, num. point / x, G, H, N, A, V\n");
  fclose(ot_file);

  mat = NULL;
  prc = NULL;
  sol = NULL;
  rhs = NULL;
  res = NULL;
  linSolve = NULL;

  iter = 0;
}

TDSolve::~TDSolve(void) {
  erase_array( data_g,  3 );
  erase_array( data_h,  3 );
  erase_array( data_n,  3 );
  erase_array( data_a,  3 );
  erase_array( data_v,  3 );
  erase_array( data_na, 3 );

  if (mat!=NULL) delete mat;
  if (prc!=NULL) delete prc;
  if (rhs!=NULL) delete rhs;
  if (res!=NULL) delete res;
  if (sol!=NULL) delete sol;
  if (linSolve!=NULL) delete linSolve;
}

void TDSolve::setup(void) {
  mat = new LinSqSparseMat;
  mat->setBlckDim(1);
  mat->setEqualSize(Nx+1);

  prc = new LinSqSparseMat;
  prc->setBlckDim(1);
  prc->setEqualSize(Nx+1);

  sol = new LinVec;
  sol->setBlckDim(1);
  sol->setSize(Nx+1);

  rhs = new LinVec;
  rhs->setBlckDim(1);
  rhs->setSize(Nx+1);

  res = new LinVec;
  res->setBlckDim(1);
  res->setSize(Nx+1);

  linSolve = new LinSolve;
  linSolve->setBlckDim(1);
  linSolve->setSize(Nx+1);

  linSolve->linkMat(mat);
  linSolve->linkPrc(prc);
  linSolve->linkSol(sol);
  linSolve->linkRhs(rhs);
  linSolve->linkRes(res);

  linSolve->setPrintStatus(tdLin->Print);
  linSolve->setSolvePrec(tdLin->Prec,tdLin->PrecPar);
  switch (tdLin->Type) {
    case LSYS_TYPE_DIR:
      linSolve->setSolveDir();
      break;
    case LSYS_TYPE_ITER:
      linSolve->setMaxIter(tdLin->MaxIter);
      linSolve->setTargetRes(tdLin->TargetRes);
      linSolve->setSolveIter(tdLin->Iter);
      break;
    case LSYS_TYPE_KRYLOV:
      linSolve->setMaxIter(tdLin->MaxIter);
      linSolve->setTargetRes(tdLin->TargetRes);
      linSolve->setSolveKrylov(tdLin->Krylov,tdLin->KrylovCap);
      break;
  }

  return;
}


void TDSolve::set_ic(void) {
  switch ( tdSetup->ic ) {
    case 0: 
      {
        data_r[0] = 2.0;
        (*data_g[0]) = 0.84 * tdCoef->theta;
        (*data_h[0]) = 0.10 * tdCoef->theta;
        (*data_n[0]) = 0.06 * tdCoef->theta;
        DType chem_A = 8.0, chem_a = 1.0;
        for ( int i = 0; i < Nx+1; i++ ) {
          DType tvar = 2.0 * x_f[i] / chem_a;
          (*data_a[0])[i] = chem_A * std::exp( -tvar*tvar);
        }
      }
      break;
    default:
      std::cout << "Error: No initial condition is defined ... Exit ...\n" << std::endl;
      exit(1);
      break;
  }

  update_na(0);
  update_vel(0);
  return;
}

void TDSolve::update_sol(const int fr, const int to) {

  update_rad( fr, to );

  switch( tdSetup->scell ) {
    case 0:
      update_cells_upw( fr, to );
      break;
    default:
      std::cout << "Warning: The spatial scheme for cell not defined, "
        << "using upwind instead\n" << std::endl;
      update_cells_upw( fr, to );
      break;
  }

  switch( tdSetup->schem ) {
    case 0:
      update_chem_upw_be( fr, to );
      break;
    default:
      std::cout << "Warning: The spatial scheme for chem not defined, "
        << "using upw/backward instead\n" << std::endl;
      update_chem_upw_be( fr, to );
      break;
  }

  update_na( to );
  update_vel( to );
  return;
}

void TDSolve::average_sol(const DType coef_1, const int fr_1,
    const DType coef_2, const int fr_2, const int to) {
  *data_g[to] = coef_1 * (*data_g[fr_1]) + coef_2 * (*data_g[fr_2]);
  *data_h[to] = coef_1 * (*data_h[fr_1]) + coef_2 * (*data_h[fr_2]);
  *data_n[to] = coef_1 * (*data_n[fr_1]) + coef_2 * (*data_n[fr_2]);
  *data_a[to] = coef_1 * (*data_a[fr_1]) + coef_2 * (*data_a[fr_2]);
  *data_v[to] = coef_1 * (*data_v[fr_1]) + coef_2 * (*data_v[fr_2]);
  data_r[to]  = coef_1 * data_r[fr_1]    + coef_2 * data_r[fr_2];
  return;
}

void TDSolve::update_rad(const int fr, const int to) {
  data_r[to] = data_r[fr] + dtn * (*data_v[fr])[Nx]; 
  return;
}

void TDSolve::update_vel(const int nt) {
  switch ( tdSetup->vel ) {
    case 0:
      update_vel_cd(nt);
      break;
    default:
      std::cout << "Warning: The velocity scheme not defined, "
        << "using central-difference instead\n" << std::endl;
      update_vel_cd(nt);
      break;
  }
  return;
}

void TDSolve::update_vel_cd(const int nt) {

  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0

  for ( int i = 0; i < Nx; i++ ) { // data_v[nt][i+1]
    DType tvar = x_f[i]*x_f[i]*(*data_v[nt])[i];
    tvar += dx*x_c[i]*x_c[i]*data_r[nt]*
      0.5*( tdCoef->lambda*((*data_g[nt])[i+1]+(*data_g[nt])[i])-
            tdCoef->delta*((*data_h[nt])[i+1]+(*data_h[nt])[i]) )
      /(tdCoef->theta);
    tvar -= (tdCoef->alpha/tdCoef->theta/data_r[nt]) *
      ( (*data_na[nt])[i+1] - (*data_na[nt])[i] );
    (*data_v[nt])[i+1] = tvar / x_f[i+1] / x_f[i+1];
  }

  return;
}

void TDSolve::update_cells_upw(const int fr, const int to) {
  DType lambda = tdCoef->lambda;
  DType delta  = tdCoef->delta;
  DType mu     = tdCoef->mu;
  DType alpha  = tdCoef->alpha;
  DType theta  = tdCoef->theta;

  HVEC rhs_srce(Nx+1);
  
  DType a0 = (*data_a[fr])[0];
  DType a1 = (*data_a[fr])[1];
  DType a2 = (*data_a[fr])[2];
  DType na_over_x2 = 2.0*(*data_n[fr])[0]*(a1-a0)/dx/x_c[0] +
    (*data_n[fr])[0] *(a2-2.0*a1+a0)/dx/dx;

  for ( int i = 0; i < Nx+1; i++ )
    rhs_srce[i] = lambda * (*data_g[fr])[i] - delta * (*data_h[fr])[i];
  rhs_srce[0] -= alpha/(data_r[fr]*data_r[fr]) * na_over_x2;
  for ( int i = 1; i < Nx+1; i++ )
    rhs_srce[i] -= alpha/(x_f[i]*x_f[i]*data_r[fr]*data_r[fr])*(*data_na[fr])[i];

  // Source terms
  for ( int i = 0; i < Nx+1; i++ ) {
    (*data_g[to])[i] = (*data_g[fr])[i] + dtn * (
        lambda * (*data_g[fr])[i] - mu * (*data_g[fr])[i] - 
        rhs_srce[i] * (*data_g[fr])[i] / theta );
    (*data_h[to])[i] = (*data_h[fr])[i] + dtn * (
        mu * (*data_g[fr])[i] - delta * (*data_h[fr])[i] - 
        rhs_srce[i] * (*data_h[fr])[i] / theta );
    (*data_n[to])[i] = (*data_n[fr])[i] + dtn * (
        -rhs_srce[i] * (*data_n[fr])[i] / theta - 
        ( (i==0) ? alpha*na_over_x2/(data_r[fr]*data_r[fr]) :
        alpha/(x_f[i]*x_f[i]*data_r[fr]*data_r[fr])*(*data_na[fr])[i] ) );
  }
  // Advection term, only needed for interior part 
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_vel = 
      ((*data_v[fr])[i]-x_f[i]*(*data_v[fr])[Nx])/data_r[fr];
    (*data_g[to])[i] -= (loc_vel*dtn/dx) * (
        (loc_vel>=0.0) ? ((*data_g[fr])[i]-(*data_g[fr])[i-1]) : 
        ((*data_g[fr])[i+1]-(*data_g[fr])[i]) );
    (*data_h[to])[i] -= (loc_vel*dtn/dx) * (
        (loc_vel>=0.0) ? ((*data_h[fr])[i]-(*data_h[fr])[i-1]) : 
        ((*data_h[fr])[i+1]-(*data_h[fr])[i]) );
    (*data_n[to])[i] -= (loc_vel*dtn/dx) * (
        (loc_vel>=0.0) ? ((*data_n[fr])[i]-(*data_n[fr])[i-1]) : 
        ((*data_n[fr])[i+1]-(*data_n[fr])[i]) );
  }

  return;
}

void TDSolve::update_chem_upw_be(const int fr, const int to) {
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  DType m     = tdCoef->m;
  DType D     = tdCoef->m;
  DType beta  = tdCoef->beta;
  DType alpha = tdCoef->alpha;
  DType gamma = tdCoef->gamma;

  // RHS
  VComp comp(1,0.0);
  for ( int i = 1; i < Nx; i++ ) {
    comp = (*data_a[fr])[i] + 
      dtn * m * (*data_g[to])[i] / (beta+(*data_g[to])[i]);
    rhs->set(i,comp);
  }
  comp = -((*data_n[to])[Nx] - (*data_n[to])[Nx-1])/alpha;
  rhs->set(Nx,comp);

  // Matrix
  Entry elem(1,0.0); 

  // Left boundary
  elem = 1.0;
  mat->set(0,0,elem);
  elem = -1.0;
  mat->set(0,1,elem);

  // Interior node
  for ( int i = 1; i < Nx; i++ ) {
    DType l_vel = x_f[i] * (*data_v[fr])[Nx] / data_r[fr];
    DType coef_adv_l = ( l_vel >= 0.0 ) ? -1.0 :  0.0;
    DType coef_adv_c = ( l_vel >= 0.0 ) ?  1.0 : -1.0;
    DType coef_adv_r = ( l_vel >= 0.0 ) ?  0.0 :  1.0;

    DType coef_diff = 
      D*dtn/(data_r[fr]*data_r[fr]*dx*dx);
    DType x_rat_l = (x_c[i-1]/x_f[i]) * (x_c[i-1]/x_f[i]);
    DType x_rat_r = (x_c[i]/x_f[i]) * (x_c[i]/x_f[i]);

    elem = - coef_diff*x_rat_l - coef_adv_l*(dtn*l_vel/dx);
    mat->set(i,i-1,elem);

    elem = 1.0 + dtn*gamma + coef_diff*(x_rat_l+x_rat_r) -
      coef_adv_c*(dtn*l_vel/dx);
    mat->set(i,i,elem);

    elem = - coef_diff*x_rat_r - coef_adv_r*(dtn*l_vel/dx);
    mat->set(i,i+1,elem);
  }

  // Right boundary
  elem = -1.0;
  mat->set(Nx,Nx-1,elem);
  elem = 1.0;
  mat->set(Nx,Nx,elem);

  // Initialize solution
  for ( int i = 0; i < Nx+1; i++ ) {
    comp = (*data_a[fr])[i];
    sol->set(i,comp);
  }

  linSolve->solve();

  // Fetch solution
  for ( int i = 0; i < Nx+1; i++ ) {
    comp = sol->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_a[to])[i] = vsol[0];
  }

  return;
}

void TDSolve::copy_data(const int nt_fr, const int nt_to) {
  *data_g[nt_to]  = *data_g[nt_fr];
  *data_h[nt_to]  = *data_h[nt_fr];
  *data_n[nt_to]  = *data_n[nt_fr];
  *data_a[nt_to]  = *data_a[nt_fr];
  *data_v[nt_to]  = *data_v[nt_fr];
  *data_na[nt_to] = *data_na[nt_fr];
   
  data_r[nt_to]  =  data_r[nt_fr];
  return;
}

void TDSolve::update_na(const int nt) {
  HVEC tvar(Nx);
  for ( int i = 0; i < Nx; i++ )
    tvar[i] = 0.5*x_c[i]*x_c[i]*((*data_n[nt])[i]+(*data_n[nt])[i+1])
      *((*data_a[nt])[i+1]-(*data_a[nt])[i])/dx;

  (*data_na[nt])[0] = 2.0*tvar[0]/dx;
  for ( int i = 1; i < Nx; i++ )
    (*data_na[nt])[i] = (tvar[i]-tvar[i-1])/dx;

  DType tvar_r = x_f[Nx]*x_f[Nx]*(*data_n[nt])[Nx]*((*data_a[nt])[Nx]-(*data_a[nt])[Nx-1])/dx;
  (*data_na[nt])[Nx] = 2.0*(tvar_r-tvar[Nx-1])/dx;

  return;
}

void TDSolve::calc_dt(const int nt) {
  DType dtn_trial = DTYPE_MAX;
  for ( int i = 0; i < Nx+1; i++ ) {
    DType loc_vel = (*data_v[nt])[i];
    DType m_vel = x_f[i]*(*data_v[nt])[Nx]/data_r[nt];
    DType v_max = max( fabs(m_vel), fabs(loc_vel/data_r[nt]-m_vel) );
    dtn_trial = min( dtn_trial, tdSetup->cfl * dx / v_max );
  }
  dtn = dtn_trial;

  if ( tdSetup->Tterm - tsol < dtn ) {
    dtn = tdSetup->Tterm - tsol;
    reach_term = true;
  }

  tsol += dtn;
  return;
}

void TDSolve::diagnostic(const int nt) {
  fprintf( stdout, "Iter %3d: Time %6.8f; Time step %6.8f; Radius %6.8f; Front-velocity %6.8f.\n",
      iter, tsol, dtn, data_r[nt], (*data_v[nt])[Nx] );
  return;
}

void TDSolve::output(const int nt) {
  FILE* ot_file = fopen( tdInput->outputFile, "a" );
  fprintf( ot_file, "%10.20f\t%10.20f\t%d\n", tsol, dtn, Nx+1 );
  for ( int i = 0; i < Nx+1; i++ )
    fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\n",
        x_f[i] * data_r[nt], (*data_g[nt])[i], (*data_h[nt])[i],
        (*data_n[nt])[i], (*data_a[nt])[i], (*data_v[nt])[i] );
  fclose(ot_file);
}


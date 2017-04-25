#include "TDSolveFDM.h"
#include "Alloc_Utils.h"
#include "TDInput.h"
#include "Alg_Solve.h"
using namespace std;

TDSolveFDM::TDSolveFDM(TDInput* tdInput_ptr) : 
  TDSolve(tdInput_ptr) {

  alloc_array_full( data_g,  Nx+1, 3 );
  alloc_array_full( data_h,  Nx+1, 3 );
  alloc_array_full( data_n,  Nx+1, 3 );
  alloc_array_full( data_v,  Nx+1, 3 );
  alloc_array_full( data_dv, Nx+1, 3 );
  alloc_array_full( data_na, Nx+1, 3 );
  
  alloc_array_full( data_a,  Na+1, 3 );

  FILE* ot_file = fopen( tdInput->outputFile, "w" );
  fprintf( ot_file, "Time, time step size, num. point / x, G, H, N, A, V\n");
  fclose(ot_file);
}

TDSolveFDM::~TDSolveFDM(void) {
  erase_array( data_dv, 3 );
  erase_array( data_na, 3 );
}

void TDSolveFDM::setup(void) {

  setupLinSolve( mat, prc, sol, rhs, res, linSolve, Nx+1 );
  setupLinSolve( mat_a, prc_a, sol_a, rhs_a, res_a, linSolve_a, Na+1 );

  return;
}


void TDSolveFDM::set_ic(void) {
  data_r[0] = tdCoef->R_0;
  switch ( tdSetup->ic ) {
    case 0: 
      {
        (*data_g[0]) = 0.84;
        (*data_h[0]) = 0.16-tdCoef->n_bc;
        (*data_n[0]) = tdCoef->n_bc;
        DType a_coef = 1.0e+3;
        for ( int i = 0; i < Na+1; i++ ) {
          DType tvar = x_a_f[i] * data_r[0];
          (*data_a[0])[i] = a_coef * std::exp( -tvar*tvar );
        }
      }
      break;
    default:
      std::cout << "Error: No initial condition is defined ... Exit ...\n" << std::endl;
      exit(1);
      break;
  }

  update_vel(0);

  FILE* ot_file = fopen( tdInput->radiusFile, "w" );
  fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\n", tsol, data_r[0], (*data_v[0])[Nx]);
  fclose(ot_file);

  return;
}

void TDSolveFDM::copy_data(const int nt_fr, const int nt_to) {
  *data_g[nt_to]  = *data_g[nt_fr];
  *data_h[nt_to]  = *data_h[nt_fr];
  *data_n[nt_to]  = *data_n[nt_fr];
  *data_a[nt_to]  = *data_a[nt_fr];
  *data_v[nt_to]  = *data_v[nt_fr];
  *data_dv[nt_to] = *data_dv[nt_fr];
  *data_na[nt_to] = *data_na[nt_fr];
   
  data_r[nt_to]  =  data_r[nt_fr];
  return;
}

void TDSolveFDM::update_sol(const int fr, const int to) {

  update_rad( fr, to );

  switch( tdSetup->scell ) {
    case 0: // Explicit first-order upwind
      update_cells_upw( fr, to );
      break;
    case 1: // Implicit first-order upwind
      update_cells_upw_be( fr, to );
      break;
    case 2: // Implicit fourth-order Pade method 
      //update_cells_pade_be( fr, to );
      break;
    default:
      std::cout << "Warning: The spatial scheme for cell not defined, "
        << "using upwind instead\n" << std::endl;
      update_cells_upw( fr, to );
      break;
  }

  switch( tdSetup->schem ) {
    case 0: // Implicit first-order upwind for advection and cd for diffusion
      update_chem_upw_be( fr, to );
      break;
    case 1: // Pade for advection and cd for diffusion
      //update_chem_pade_be( fr, to );
      break;
    default:
      std::cout << "Warning: The spatial scheme for chem not defined, "
        << "using upw/backward instead\n" << std::endl;
      update_chem_upw_be( fr, to );
      break;
  }

  update_vel( to );
  return;
}

void TDSolveFDM::update_vel(const int nt) {
  switch ( tdSetup->vel ) {
    case 0: // 2nd-order integral form
      update_vel_int(nt);
      break;
    case 1: // Upwind for advection
      //update_vel_upw(nt);
      break;
    case 2: // BDF2 for advection
      //update_vel_bdf2(nt);
      break;
    case 3: // Compact difference for advection
      //update_vel_cfd(nt);
      break;
    default:
      std::cout << "Warning: The velocity scheme not defined, "
        << "using integral form instead\n" << std::endl;
      update_vel_int(nt);
      break;
  }
  return;
}

void TDSolveFDM::update_vel_upw(const int nt) {

  DType lambda = tdCoef->lambda;
  DType delta  = tdCoef->delta;
  DType alpha  = tdCoef->alpha;

  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0

  for ( int i = 0; i < Nx+1; i++ )
    (*data_dv[nt])[i] = lambda * (*data_g[nt])[i] - delta * (*data_h[nt])[i] - 
      alpha * (*data_na[nt])[i];

  for ( int i = 0; i < Nx; i++ ) { // data_v[nt][i+1]
    DType x_rat = x_f[i]/x_f[i+1];
    (*data_v[nt])[i+1] = x_rat * x_rat * (*data_v[nt])[i] + 
      dx * data_r[nt] * (*data_dv[nt])[i+1];
  }

  return;
}

void TDSolveFDM::update_vel_bdf2(const int nt) {

  DType lambda = tdCoef->lambda;
  DType delta  = tdCoef->delta;
  DType alpha  = tdCoef->alpha;

  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0

  for ( int i = 0; i < Nx+1; i++ )
    (*data_dv[nt])[i] = lambda * (*data_g[nt])[i] - delta * (*data_h[nt])[i] - 
      alpha * (*data_na[nt])[i];

  // data_v[nt][1]
  (*data_v[nt])[1] = dx * data_r[nt] * (*data_dv[nt])[1];

  for ( int i = 1; i < Nx; i++ ) { // data_v[nt][i+1]
    DType x_rat_m1 = x_f[i]/x_f[i+1];
    DType x_rat_m2 = x_f[i-1]/x_f[i+1];
    (*data_v[nt])[i+1] = (4.0/3.0) * x_rat_m1 * x_rat_m1 * (*data_v[nt])[i] - 
                         (1.0/3.0) * x_rat_m2 * x_rat_m2 * (*data_v[nt])[i-1] + 
                         (2.0/3.0) * dx * data_r[nt] * (*data_dv[nt])[i+1];
  }

  return;
}

void TDSolveFDM::update_vel_int(const int nt) {

  DType lambda = tdCoef->lambda;
  DType delta  = tdCoef->delta;
  DType alpha  = tdCoef->alpha;

  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0

  // First do the integration of x^2R(lambda G - delta H)
  for ( int i = 0; i < Nx; i++ ) {
    int j = i+1;
    (*data_v[nt])[j] = (*data_v[nt])[i] + data_r[nt] * 0.5 * dx *
      ( x_f[i] * x_f[i] * ( lambda * (*data_g[nt])[i] - delta * (*data_h[nt])[i] ) +
        x_f[j] * x_f[j] * ( lambda * (*data_g[nt])[j] - delta * (*data_h[nt])[j] ) );
  }

  // Next divide by x^2 and remove the contribution from N and dA
  for ( int i = 1; i < Nx+1; i++ ) {
    (*data_v[nt])[i] = (*data_v[nt])[i]/x_f[i]/x_f[i] - (alpha/data_r[nt]) * 
      (*data_n[nt])[i] * ((*data_a[nt])[i+1]-(*data_a[nt])[i-1]) / (2.0*dx);
  }

  // Lastly store d(x^2V)/dx/(x^2R) in data_dv
  // Interior nodes
  for ( int i = 1; i < Nx; i++ ) 
    (*data_dv[nt])[i] = ((*data_v[nt])[i+1]-(*data_v[nt])[i-1])/(2.0*dx*data_r[nt]) +
                        2.0 * (*data_v[nt])[i] / x_f[i] / data_r[nt];
  // Left boundary
  (*data_dv[nt])[0] = (3.0/data_r[nt]) * ((*data_v[nt])[1]-(*data_v[nt])[0])/dx;
  
  // Right boundary
  (*data_dv[nt])[Nx] = (1.0/data_r[nt]) * ((*data_v[nt])[Nx]-(*data_v[nt])[Nx-1]) / dx +
    2.0 * (*data_v[nt])[Nx] / x_f[Nx] / data_r[nt];

  return;
}

void TDSolveFDM::update_cells_upw(const int fr, const int to) {
  DType lambda = tdCoef->lambda;
  DType delta  = tdCoef->delta;
  DType mu     = tdCoef->mu;
  DType alpha  = tdCoef->alpha;

  // Source terms -- start from i=1, and enforce hom. Neumann b.c. later 
  for ( int i = 1; i < Nx+1; i++ ) {
    DType react_coef = (*data_dv[fr])[i];
    (*data_g[to])[i] = (*data_g[fr])[i] + dtn * (
        lambda * (*data_g[fr])[i] - mu * (*data_g[fr])[i] - 
        react_coef * (*data_g[fr])[i]);
    (*data_h[to])[i] = (*data_h[fr])[i] + dtn * (
        mu * (*data_g[fr])[i] - delta * (*data_h[fr])[i] - 
        react_coef * (*data_h[fr])[i]);
    DType react_coef_n = 
      ( x_a_c[i]*x_a_c[i]*((*data_a[fr])[i+1]-(*data_a[fr])[i]) - 
        x_a_c[i-1]*x_a_c[i-1]*((*data_a[fr])[i]-(*data_a[fr])[i-1]) ) /dx/dx;
    react_coef_n *= alpha/(x_f[i]*data_r[fr])/(x_f[i]*data_r[fr]);
    (*data_n[to])[i] = (*data_n[fr])[i] + dtn * (
        - (react_coef + react_coef_n) * (*data_n[fr])[i] );
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
    loc_vel += (alpha/data_r[fr]/data_r[fr]) * 
      ((*data_a[fr])[i+1]-(*data_a[fr])[i-1])/(2.0*dx);
    (*data_n[to])[i] -= (loc_vel*dtn/dx) * (
        (loc_vel>=0.0) ? ((*data_n[fr])[i]-(*data_n[fr])[i-1]) : 
        ((*data_n[fr])[i+1]-(*data_n[fr])[i]) );
  }
  
  // Advection for N at the right boundary (no flux for the other two)
  {
    DType loc_vel = (alpha/data_r[fr]/data_r[fr]) * ((*data_a[fr])[Nx+1]-(*data_a[fr])[Nx-1])/(2.0*dx);
    (*data_n[to])[Nx] -= (loc_vel*dtn/dx) * (
        (loc_vel>=0.0) ? ((*data_n[fr])[Nx]-(*data_n[fr])[Nx-1]) :
        ( tdCoef->n_bc - (*data_n[fr])[Nx] ) );
  }

  if ( tdSetup->lbc == 0 ) {
    // hom. Neumann b.c. at the origin
    (*data_g[to])[0] = (*data_g[to])[1];
    (*data_h[to])[0] = (*data_h[to])[1];
    (*data_n[to])[0] = (*data_n[to])[1];
  }
  else {
    DType react_coef = (*data_dv[fr])[0];
    (*data_g[to])[0] = (*data_g[fr])[0] + dtn * (
        lambda * (*data_g[fr])[0] - mu * (*data_g[fr])[0] - 
        react_coef * (*data_g[fr])[0]);
    (*data_h[to])[0] = (*data_h[fr])[0] + dtn * (
        mu * (*data_g[fr])[0] - delta * (*data_h[fr])[0] - 
        react_coef * (*data_h[fr])[0]);
    DType react_coef_n = 3.0 * ((*data_a[fr])[1] - (*data_a[fr])[0]) /dx/dx;
    react_coef_n *= alpha/(data_r[fr]*data_r[fr]);
    (*data_n[to])[0] = (*data_n[fr])[0] + dtn * (
        - (react_coef + react_coef_n) * (*data_n[fr])[0] );
  }

  return;
}

void TDSolveFDM::update_cells_upw_be(const int fr, const int to) {
  DType lambda = tdCoef->lambda;
  DType delta  = tdCoef->delta;
  DType mu     = tdCoef->mu;
  DType alpha  = tdCoef->alpha;

  VComp comp(1,0.0);
  Entry elem(1,0.0);

  // First data_g
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  // RHS
  // Left boundary -- no right hand side due to the no flux b.c.
  if ( tdSetup->lbc != 0 ) {
    comp = (*data_g[fr])[0];
    rhs->set(0,comp);
  }
  // Interior nodes and the right boundary
  for ( int i = 1; i < Nx+1; i++ ) {
    comp = (*data_g[fr])[i];
    rhs->set(i,comp);
  }

  // Matrix
  if ( tdSetup->lbc == 0 ) { // Left boundary -- no flux b.c.
    elem = 1.0;
    mat->set(0,0,elem);
    elem = -1.0;
    mat->set(0,1,elem);
  }
  else { // Left boundary -- ode b.c.
    elem = 1.0 - dtn * ( lambda - mu - (*data_dv[fr])[0] );
    mat->set(0,0,elem);
  }
  // Interior node
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_vel = ((*data_v[fr])[i]-x_f[i]*(*data_v[fr])[Nx])/data_r[fr];
    DType adv_coef = loc_vel * dtn / dx;

    if ( loc_vel >= 0.0 ) {
      elem = 1.0 + adv_coef - dtn * ( lambda - mu - (*data_dv[fr])[i] );
      mat->set(i,i,elem);
      elem = -adv_coef;
      mat->set(i,i-1,elem);
    }
    else {
      elem = 1.0 - adv_coef - dtn * ( lambda - mu - (*data_dv[fr])[i] );
      mat->set(i,i,elem);
      elem = adv_coef;
      mat->set(i,i+1,elem);
    }
  }
  // Right boundary
  elem = 1.0 - dtn * ( lambda - mu - (*data_dv[fr])[Nx] );
  mat->set(Nx,Nx,elem);

  // Initialize solution
  for ( int i = 0; i < Nx+1; i++ ) {
    comp = (*data_g[fr])[i];
    sol->set(i,comp);
  }

  linSolve->solve();

  // Fetch solution
  for ( int i = 0; i < Nx+1; i++ ) {
    comp = sol->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_g[to])[i] = vsol[0];
  }


  // Second data_h
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  // RHS
  // Left boundary -- no right hand side due to no flux b.c.
  if ( tdSetup->lbc != 0 ) {
    comp = (*data_h[fr])[0] + dtn * mu * (*data_g[to])[0];
    rhs->set(0,comp);
  }
  // Interior nodes
  for ( int i = 1; i < Nx+1; i++ ) {
    comp = (*data_h[fr])[i] + dtn * mu * (*data_g[to])[i];
    rhs->set(i,comp);
  }

  // Matrix
  if ( tdSetup->lbc == 0 ) { // Left boundary -- no flux b.c.
    elem = 1.0;
    mat->set(0,0,elem);
    elem = -1.0;
    mat->set(0,1,elem);
  }
  else { // Left boundary -- ode b.c.
    elem = 1.0 - dtn * ( lambda - mu - (*data_dv[fr])[0] );
    mat->set(0,0,elem);
  }
  // Interior node
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_vel = 
      ((*data_v[fr])[i]-x_f[i]*(*data_v[fr])[Nx])/data_r[fr];
    DType adv_coef = loc_vel * dtn / dx;

    if ( loc_vel >= 0.0 ) {
      elem = 1.0 + adv_coef - dtn * ( - delta - (*data_dv[fr])[i] );
      mat->set(i,i,elem);
      elem = - adv_coef;
      mat->set(i,i-1,elem);
    }
    else {
      elem = 1.0 - adv_coef - dtn * ( - delta - (*data_dv[fr])[i] );
      mat->set(i,i,elem);
      elem = adv_coef;
      mat->set(i,i+1,elem);
    }
  }
  // Right boundary
  elem = 1.0 - dtn * ( - delta - (*data_dv[fr])[Nx] );
  mat->set(Nx,Nx,elem);

  // Initialize solution
  for ( int i = 0; i < Nx+1; i++ ) {
    comp = (*data_h[fr])[i];
    sol->set(i,comp);
  }
  
  linSolve->solve();

  // Fetch solution
  for ( int i = 0; i < Nx+1; i++ ) {
    comp = sol->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_h[to])[i] = vsol[0];
  }

  
  // Third data_n
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  // RHS
  // Left boundary -- no right hand side due to the no fux b.c.
  if ( tdSetup->lbc != 0 ) {
    comp = (*data_n[fr])[0];
    rhs->set(0,comp);
  }
  // Interior nodes
  for ( int i = 1; i < Nx; i++ ) {
    comp = (*data_n[fr])[i];
    rhs->set(i,comp);
  }
  // Right boundary -- possible incoming data
  {
    DType loc_vel = (alpha/data_r[fr]/data_r[fr]) * 
      ((*data_a[fr])[Nx+1]-(*data_a[fr])[Nx-1])/(2.0*dx);
    if ( loc_vel >= 0.0 ) {
      comp = (*data_n[fr])[Nx];
      rhs->set(Nx,comp);
    }
    else {
      comp = (*data_n[fr])[Nx] - loc_vel * (dtn/dx) * tdCoef->n_bc;
      rhs->set(Nx,comp);
    }
  }

  // Matrix
  if ( tdSetup->lbc == 0 ) { // Left boundary -- no flux b.c.
    elem = 1.0;
    mat->set(0,0,elem);
    elem = -1.0;
    mat->set(0,1,elem);
  }
  else { // Left boundary -- ode b.c.
    DType react_coef = (*data_dv[fr])[0];
    DType react_coef_n = 3.0 * ((*data_a[fr])[1] - (*data_a[fr])[0]) /dx/dx;
    react_coef_n *= alpha/(data_r[fr]*data_r[fr]);
    elem = 1.0 - dtn * ( - react_coef - react_coef_n );
    mat->set(0,0,elem);
  }
  // Interior node
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_vel = 
      ((*data_v[fr])[i]-x_f[i]*(*data_v[fr])[Nx])/data_r[fr] + 
      (alpha/data_r[fr]/data_r[fr]) * ((*data_a[fr])[i+1]-(*data_a[fr])[i-1])/(2.0*dx);
    DType adv_coef = loc_vel * dtn / dx;

    DType react_coef = (*data_dv[fr])[i];
    DType react_coef_n = 
      ( x_a_c[i]*x_a_c[i]*((*data_a[fr])[i+1]-(*data_a[fr])[i]) - 
        x_a_c[i-1]*x_a_c[i-1]*((*data_a[fr])[i]-(*data_a[fr])[i-1]) ) /dx/dx;
    react_coef_n *= alpha/(x_f[i]*data_r[fr])/(x_f[i]*data_r[fr]);

    if ( loc_vel >= 0.0 ) {
      elem = 1.0 + adv_coef - dtn * ( - react_coef - react_coef_n );
      mat->set(i,i,elem);
      elem = - adv_coef;
      mat->set(i,i-1,elem);
    }
    else {
      elem = 1.0 - adv_coef - dtn * ( -react_coef - react_coef_n );
      mat->set(i,i,elem);
      elem = adv_coef;
      mat->set(i,i+1,elem);
    }
  }
  // Right node
  {
    DType loc_vel = (alpha/data_r[fr]/data_r[fr]) * 
      ((*data_a[fr])[Nx+1]-(*data_a[fr])[Nx-1])/(2.0*dx);
    DType adv_coef = loc_vel * dtn / dx;

    DType react_coef = (*data_dv[fr])[Nx];
    DType react_coef_n = 
      ( x_a_c[Nx]*x_a_c[Nx]*((*data_a[fr])[Nx+1]-(*data_a[fr])[Nx]) - 
        x_a_c[Nx-1]*x_a_c[Nx-1]*((*data_a[fr])[Nx]-(*data_a[fr])[Nx-1]) ) /dx/dx;
    react_coef_n *= alpha/(x_f[Nx]*data_r[fr])/(x_f[Nx]*data_r[fr]);

    if ( loc_vel >= 0.0 ) {
      elem = 1.0 + adv_coef - dtn * ( - react_coef - react_coef_n );
      mat->set(Nx,Nx,elem);
      elem = - adv_coef;
      mat->set(Nx,Nx-1,elem);
    }
    else {
      elem = 1.0 - adv_coef - dtn * ( - react_coef - react_coef_n );
      mat->set(Nx,Nx,elem);
    }
  }

  // Initialize solution
  for ( int i = 0; i < Nx+1; i++ ) {
    comp = (*data_n[fr])[i];
    sol->set(i,comp);
  }

  linSolve->solve();

  // Fetch solution
  for ( int i = 0; i < Nx+1; i++ ) {
    comp = sol->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_n[to])[i] = vsol[0];
  }

  return;
}

void TDSolveFDM::update_cells_pade_be(const int fr, const int to) {
  DType lambda = tdCoef->lambda;
  DType delta  = tdCoef->delta;
  DType mu     = tdCoef->mu;
  DType alpha  = tdCoef->alpha;

  VComp comp(1,0.0);
  Entry elem(1,0.0);
  DType one_sixth = 1.0/6.0;
  DType for_sixth = 4.0/6.0;
  //DType one_third = 1.0/3.0;
  //DType two_third = 2.0/3.0;

  // First data_g
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  // RHS
  // Left boundary -- no right hand side due to the hom. Neumann b.c.
  // Interior nodes
  for ( int i = 1; i < Nx; i++ ) {
    comp = one_sixth * (*data_g[fr])[i-1] +
           for_sixth * (*data_g[fr])[i] +
           one_sixth * (*data_g[fr])[i+1];
    rhs->set(i,comp);
  }
  // Right boundary 
  if ( tdSetup->rbc == 0 ) {
    // Do nothing due to no flux condition for G
  }
  else {
    //comp = one_third * (*data_g[fr])[Nx-1] + 
    //       two_third * (*data_g[fr])[Nx];
    comp = (*data_g[fr])[Nx];
    rhs->set(Nx,comp);
  }

  // Matrix
  // Left boundary -- hom. Neumann b.c.
  elem = 1.0;
  mat->set(0,0,elem);
  elem = -1.0;
  mat->set(0,1,elem);
  // Interior node
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_vel = 
      ((*data_v[fr])[i]-x_f[i]*(*data_v[fr])[Nx])/data_r[fr];
    DType adv_coef = 0.5 * loc_vel * dtn / dx;
    
    elem = one_sixth - adv_coef;
    mat->set(i,i-1,elem);

    elem = for_sixth - dtn * ( lambda - mu - (*data_dv[fr])[i] );
    mat->set(i,i,elem);

    elem = one_sixth + adv_coef;
    mat->set(i,i+1,elem);
  }
  // Right boundary
  if ( tdSetup->rbc == 0 ) { // No flux condition
    elem = -1.0;
    mat->set(Nx,Nx-1,elem);
    elem = 1.0;
    mat->set(Nx,Nx,elem);
  }
  else {
    //elem = one_third;
    //mat->set(Nx,Nx-1,elem);

    //elem = two_third - dtn * ( lambda - mu - (*data_dv[fr])[Nx] );
    elem = 1.0 - dtn * ( lambda - mu - (*data_dv[fr])[Nx] );
    mat->set(Nx,Nx,elem);
  }

  // Initialize solution
  for ( int i = 0; i < Nx+1; i++ ) {
    comp = (*data_g[fr])[i];
    sol->set(i,comp);
  }

  linSolve->solve();

  // Fetch solution
  for ( int i = 0; i < Nx+1; i++ ) {
    comp = sol->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_g[to])[i] = vsol[0];
  }


  // Second data_h
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  // RHS
  // Left boundary -- no right hand side due to the hom. Neumann b.c.
  // Interior nodes
  for ( int i = 1; i < Nx; i++ ) {
    comp = one_sixth * (*data_h[fr])[i-1] +
           for_sixth * (*data_h[fr])[i] +
           one_sixth * (*data_h[fr])[i+1] +
           dtn * mu * (*data_g[to])[i];
    rhs->set(i,comp);
  }
  // Right boundary 
  if ( tdSetup->rbc == 0 ) {
    // Do nothing due to no flux condition for G
  }
  else {
    //comp = one_third * (*data_h[fr])[Nx-1] + 
    //       two_third * (*data_h[fr])[Nx] + 
    //       dtn * mu * (*data_g[to])[Nx];
    comp = (*data_h[fr])[Nx] + 
           dtn * mu * (*data_g[to])[Nx];
    rhs->set(Nx,comp);
  }

  // Matrix
  // Left boundary -- hom. Neumann b.c.
  elem = 1.0;
  mat->set(0,0,elem);
  elem = -1.0;
  mat->set(0,1,elem);
  // Interior node
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_vel = 
      ((*data_v[fr])[i]-x_f[i]*(*data_v[fr])[Nx])/data_r[fr];
    DType adv_coef = 0.5 * loc_vel * dtn / dx;
    
    elem = one_sixth - adv_coef;
    mat->set(i,i-1,elem);

    elem = for_sixth - dtn * ( - delta - (*data_dv[fr])[i] );
    mat->set(i,i,elem);

    elem = one_sixth + adv_coef;
    mat->set(i,i+1,elem);
  }
  // Right boundary
  if ( tdSetup->rbc == 0 ) { // No flux condition
    elem = -1.0;
    mat->set(Nx,Nx-1,elem);
    elem = 1.0;
    mat->set(Nx,Nx,elem);
  }
  else {
    //elem = one_third;
    //mat->set(Nx,Nx-1,elem);

    //elem = two_third - dtn * ( - delta - (*data_dv[fr])[Nx] );
    elem = 1.0 - dtn * ( - delta - (*data_dv[fr])[Nx] );
    mat->set(Nx,Nx,elem);
  }

  // Initialize solution
  for ( int i = 0; i < Nx+1; i++ ) {
    comp = (*data_h[fr])[i];
    sol->set(i,comp);
  }

  linSolve->solve();

  // Fetch solution
  for ( int i = 0; i < Nx+1; i++ ) {
    comp = sol->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_h[to])[i] = vsol[0];
  }


  // Third data_n
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  // RHS
  // Left boundary -- no right hand side due to the hom. Neumann b.c.
  // Interior nodes
  for ( int i = 1; i < Nx; i++ ) {
    comp = one_sixth * (*data_n[fr])[i-1] + 
           for_sixth * (*data_n[fr])[i] + 
           one_sixth * (*data_n[fr])[i+1] -
           dtn * alpha * (*data_na[fr])[i];
    rhs->set(i,comp);
  }
  // Right boundary
  //comp = one_third * (*data_n[fr])[Nx-1] + 
  //       two_third * (*data_n[fr])[Nx] -
  //       dtn * alpha * (*data_na[fr])[Nx];
  comp = (*data_n[fr])[Nx] -
         dtn * alpha * (*data_na[fr])[Nx];

  // Matrix
  // Left boundary -- hom. Neumann b.c.
  elem = 1.0;
  mat->set(0,0,elem);
  elem = -1.0;
  mat->set(0,1,elem);
  // Interior node
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_vel = 
      ((*data_v[fr])[i]-x_f[i]*(*data_v[fr])[Nx])/data_r[fr];
    DType adv_coef = 0.5 * loc_vel * dtn / dx;
    
    elem = one_sixth - adv_coef;
    mat->set(i,i-1,elem);

    elem = for_sixth - dtn * ( - (*data_dv[fr])[i] );
    mat->set(i,i,elem);

    elem = one_sixth + adv_coef;
    mat->set(i,i+1,elem);
  }
  // Right boundary
  //elem = one_third;
  //mat->set(Nx,Nx-1,elem);
  
  //elem = two_third - dtn * ( - (*data_dv[fr])[Nx] );
  elem = 1.0 - dtn * ( - (*data_dv[fr])[Nx] );
  mat->set(Nx,Nx,elem);

  // Initialize solution
  for ( int i = 0; i < Nx+1; i++ ) {
    comp = (*data_n[fr])[i];
    sol->set(i,comp);
  }

  linSolve->solve();

  // Fetch solution
  for ( int i = 0; i < Nx+1; i++ ) {
    comp = sol->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_n[to])[i] = vsol[0];
  }

  return;
}

void TDSolveFDM::update_chem_upw_be(const int fr, const int to) {
  rhs_a->reset( 0.0 );
  mat_a->reset( 0.0 );
  sol_a->reset( 0.0 );

  DType m     = tdCoef->m;
  DType D     = tdCoef->D;
  DType beta  = tdCoef->beta;
  DType gamma = tdCoef->gamma;
  DType theta = tdCoef->theta;
  
  // RHS
  VComp comp(1,0.0);
  for ( int i = 1; i < Na; i++ ) {
    comp = (*data_a[fr])[i];
    if ( i < Nx+1 ) // G is active
      comp += dtn * m * (*data_g[to])[i] / (beta/theta+(*data_g[to])[i]);
    rhs_a->set(i,comp);
  }

  // Matrix
  Entry elem(1,0.0); 

  // Left boundary
  elem = 1.0;
  mat_a->set(0,0,elem);
  elem = -1.0;
  mat_a->set(0,1,elem);

  // Interior node
  for ( int i = 1; i < Na; i++ ) {
    DType l_vel = - x_a_f[i] * (*data_v[fr])[Nx] / data_r[fr];
    DType coef_adv_l = ( l_vel >= 0.0 ) ? -1.0 :  0.0;
    DType coef_adv_c = ( l_vel >= 0.0 ) ?  1.0 : -1.0;
    DType coef_adv_r = ( l_vel >= 0.0 ) ?  0.0 :  1.0;

    DType coef_diff = 
      D*(dtn/dx/dx)/(data_r[fr]*data_r[fr]);
    DType x_rat_l = (x_a_c[i-1]/x_a_f[i]) * (x_a_c[i-1]/x_a_f[i]);
    DType x_rat_r = (x_a_c[i]/x_a_f[i]) * (x_a_c[i]/x_a_f[i]);

    elem = - coef_diff*x_rat_l + coef_adv_l*(dtn*l_vel/dx);
    mat_a->set(i,i-1,elem);

    elem = 1.0 + dtn*gamma + coef_diff*(x_rat_l+x_rat_r) +
      coef_adv_c*(dtn*l_vel/dx);
    mat_a->set(i,i,elem);

    elem = - coef_diff*x_rat_r + coef_adv_r*(dtn*l_vel/dx);
    mat_a->set(i,i+1,elem);
  }

  // Right boundary (a = 0)
  elem = 1.0;
  mat_a->set(Na,Na,elem);

  // Initialize solution
  for ( int i = 0; i < Na+1; i++ ) {
    comp = (*data_a[fr])[i];
    sol_a->set(i,comp);
  }

  linSolve_a->solve();

  // Fetch solution
  for ( int i = 0; i < Na+1; i++ ) {
    comp = sol_a->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_a[to])[i] = vsol[0];
  }

  return;
}

void TDSolveFDM::update_chem_pade_be(const int fr, const int to) {
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  DType m     = tdCoef->m;
  DType D     = tdCoef->D;
  DType beta  = tdCoef->beta;
  DType gamma = tdCoef->gamma;
  DType theta = tdCoef->theta;
  DType sigma = tdCoef->sigma;

  // RHS
  VComp comp(1,0.0);
  for ( int i = 1; i < Nx; i++ ) {
    comp = (*data_a[fr])[i] + 
      dtn * m * (*data_g[to])[i] / (beta/theta+(*data_g[to])[i]);
    rhs->set(i,comp);
  }
  comp = ((*data_n[to])[Nx] - (*data_n[to])[Nx-1])/sigma;
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
    DType l_vel = - x_f[i] * (*data_v[fr])[Nx] / data_r[fr];

    DType coef_diff = 
      D*(dtn/dx/dx)/(data_r[fr]*data_r[fr]);
    DType x_rat_l = (x_c[i-1]/x_f[i]) * (x_c[i-1]/x_f[i]);
    DType x_rat_r = (x_c[i]/x_f[i]) * (x_c[i]/x_f[i]);

    elem = 1.0/6.0 - 0.5 * dtn * l_vel / dx - coef_diff * x_rat_l;
    mat->set(i,i-1,elem);

    elem = 4.0/6.0 + dtn * gamma + coef_diff*(x_rat_l+x_rat_r);
    mat->set(i,i,elem);

    elem = 1.0/6.0 + 0.5 * dtn * l_vel / dx - coef_diff * x_rat_r;
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

/*
void TDSolveFDM::update_na(const int nt) {
  HVEC tvar(Nx);
  for ( int i = 0; i < Nx; i++ )
    tvar[i] = 0.5*x_c[i]*x_c[i]*((*data_n[nt])[i]+(*data_n[nt])[i+1])
      *((*data_a[nt])[i+1]-(*data_a[nt])[i])/dx;

  // Left end point
  (*data_na[nt])[0] = (*data_n[nt])[0]*2.0*((*data_a[nt])[1]-(*data_a[nt])[0]) / 
    ((dx*data_r[nt])*(dx*data_r[nt]));

  // Interior points, central difference
  for ( int i = 1; i < Nx; i++ )
    (*data_na[nt])[i] = (tvar[i]-tvar[i-1])/dx /
      ((x_f[i]*data_r[nt])*(x_f[i]*data_r[nt]));

  // Right end point
  DType tvar_r = x_f[Nx]*x_f[Nx]*(*data_n[nt])[Nx]*((*data_a[nt])[Nx]-(*data_a[nt])[Nx-1])/dx;
  (*data_na[nt])[Nx] = 2.0*(tvar_r-tvar[Nx-1])/dx /
      ((x_f[Nx]*data_r[nt])*(x_f[Nx]*data_r[nt]));

  return;
}
*/

void TDSolveFDM::calc_dt(const int nt) {
  if ( iter != 0 )
    tsol += dtn;

  if ( reach_term_next )
    reach_term_now = true;

  DType dtn_trial = DTYPE_MAX;
  DType alpha = tdCoef->alpha;
  for ( int i = 0; i < Nx+1; i++ ) {
    DType loc_vel = (*data_v[nt])[i]/data_r[nt];
    DType m_vel = x_f[i]*(*data_v[nt])[Nx]/data_r[nt];
    DType a_vel = ( i== 0 ) ? 0.0 : alpha * 
      ((*data_a[nt])[i+1]-(*data_a[nt])[i-1])
      /(2.0*dx)/data_r[nt]/data_r[nt];
    DType v_max = max( fabs(m_vel), fabs(loc_vel-m_vel));
    v_max = max( v_max, fabs(loc_vel-m_vel+a_vel) );
    dtn_trial = min( dtn_trial, tdSetup->cfl * dx / v_max );
  }

  if ( iter == 0 )
    dtn = min( dtn_trial, dtn );
  else
    dtn = min( dtn * tdSetup->dTgrow, dtn_trial );

  if ( tdSetup->Tterm - tsol < dtn ) {
    dtn = tdSetup->Tterm - tsol;
    reach_term_next = true;
  }

  if ( iter >= tdSetup->maxit )
    reach_term_now = true;
  return;
}

void TDSolveFDM::output(const int nt) {

  FILE* ot_file = fopen( tdInput->radiusFile, "a" );
  fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\n", tsol, data_r[nt], (*data_v[nt])[Nx] );
  fclose(ot_file);

  if ( !reach_term_now && iter % tdSetup->freq != 0 ) return;
  ot_file = fopen( tdInput->outputFile, "a" );
  fprintf( ot_file, "%10.20f\t%10.20f\t%d\n", tsol, dtn, Nx+1 );
  for ( int i = 0; i < Nx+1; i++ )
    fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\n",
        x_f[i] * data_r[nt], (*data_g[nt])[i], (*data_h[nt])[i],
        (*data_n[nt])[i], (*data_a[nt])[i], (*data_v[nt])[i] );
  fclose(ot_file);

  return;
}


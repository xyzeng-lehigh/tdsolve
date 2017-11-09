#include "TDSolveFVM.h"
#include "Alloc_Utils.h"
#include "TDInput.h"
#include "Alg_Solve.h"
#include "Alg_Limiters.h"
using namespace std;

TDSolveFVM::TDSolveFVM(TDInput* tdInput_ptr) : 
  TDSolve(tdInput_ptr) {

  alloc_array_full( data_g, Nx, 3 );
  alloc_array_full( data_h, Nx, 3 );
  alloc_array_full( data_n, Nx, 3 ); // M
  alloc_array_full( data_v, Nx+1, 3 );
  
  alloc_array_full( data_a, Na, 3 ); 

  alloc_array_full( flux_g, Nx+1, 1 );
  alloc_array_full( flux_h, Nx+1, 1 );
  alloc_array_full( flux_n, Nx+1, 1 );
  alloc_array_full( flux_a, Na+1, 1 );

  alloc_array_full( slope_g, Nx, 1 );
  alloc_array_full( slope_h, Nx, 1 );
  alloc_array_full( slope_n, Nx, 1 );
  alloc_array_full( slope_a, Na, 1 );

  // OPN model
  if ( tdSetup->prob_id == 1 ) {
    alloc_array_full (data_p1, Nx, 3);
    alloc_array_full (data_p2, Nx, 3);

    alloc_array_full (flux_p1, Nx+1, 1);
    alloc_array_full (flux_p2, Nx+1, 1);

    alloc_array_full(slope_p1, Nx, 1);
    alloc_array_full(slope_p2, Nx, 1);
  }

  alloc_array_full( N_dir, Nx+1 );

  FILE* ot_file = fopen( tdInput->outputFile, "w" );
  fprintf( ot_file, "Time, time step size, num. point / x_c, G, H, N, A /newline x_f, V\n");
  fclose(ot_file);
}

TDSolveFVM::~TDSolveFVM(void) {
  erase_array( flux_g, 1 );
  erase_array( flux_h, 1 );
  erase_array( flux_n, 1 );
  erase_array( flux_a, 1 );
  erase_array( slope_g, 1 );
  erase_array( slope_h, 1 );
  erase_array( slope_n, 1 );
  erase_array( slope_a, 1 );
  erase_array( N_dir );

  // OPN model
  if ( tdSetup->prob_id == 1 ) {
    erase_array( flux_p1, 1);
    erase_array( flux_p2, 1);
    erase_array(slope_p1, 1);
    erase_array(slope_p2, 1);
  }
}

void TDSolveFVM::setup(void) {

  setupLinSolve( mat, prc, sol, rhs, res, linSolve, Nx );
  setupLinSolve( mat_a, prc_a, sol_a, rhs_a, res_a, linSolve_a, Na );

  // OPN model
  if ( tdSetup->prob_id == 1 ) {
    setupLinSolve(mat_p1, prc_p1, sol_p1, rhs_p1, res_p1, linSolve_p1, Nx);
    setupLinSolve(mat_p2, prc_p2, sol_p2, rhs_p2, res_p2, linSolve_p2, Nx);
  }

  return;
}

void TDSolveFVM::set_ic(void) {

  if ( !tdSetup->nonDim ) {
    set_ic_dim();
    return;
  }

  data_r[0] = tdCoef->R_0;
  switch ( tdSetup->ic ) {
    case 0: 
      {
        (*data_g[0]) = 0.84;
        (*data_h[0]) = 0.16-tdCoef->n_bc;
        (*data_n[0]) = tdCoef->n_bc;
        DType a_coef = 1.0e+3;
        for ( int i = 0; i < Na; i++ ) {
          DType tvar = x_a_c[i] * data_r[0];
          (*data_a[0])[i] = a_coef * std::exp( -tvar*tvar );
        }
      }
      break;
    case 1: 
      {
        (*data_g[0]) = 5.0*(1.0-tdCoef->n_bc)/6.0;
        (*data_h[0]) =     (1.0-tdCoef->n_bc)/6.0;
        (*data_n[0]) = tdCoef->n_bc;
        DType a_coef = 1.0e+3;
        for ( int i = 0; i < Na; i++ ) {
          DType tvar = x_a_c[i] * data_r[0];
          (*data_a[0])[i] = a_coef * std::exp( -tvar*tvar );
        }
      }
      break;
    default:
      std::cout << "Error: No initial condition is defined ... Exit ...\n" << std::endl;
      exit(1);
      break;
  }

  upscale_var( *data_g[0], *data_g[0], x_c, Nx, data_r[0]);
  upscale_var( *data_h[0], *data_h[0], x_c, Nx, data_r[0]);
  upscale_var( *data_n[0], *data_n[0], x_c, Nx, data_r[0]);

  update_vel(0);

  N_dir[0] = 0;
  for ( int i = 1; i < Nx+1; i++ )
    N_dir[i] = i-1;

  FILE* ot_file = fopen( tdInput->radiusFile, "w" );
  fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\t%10.20f\n", 
      tsol, data_r[0], (*data_v[0])[Nx], 
      -((*data_a[0])[Nx]-(*data_a[0])[Nx-1])/dx/data_r[0]);
  fclose(ot_file);

  return;
}

void TDSolveFVM::set_ic_dim(void) {

  data_r[0] = tdCoef->R_0;
  DType theta = tdCoef->theta;

  switch ( tdSetup->ic ) {
    case 0: 
      {
        (*data_g[0]) = 0.84 * theta;
        (*data_h[0]) = (0.16-tdCoef->n_bc) * theta;
        (*data_n[0]) = (tdCoef->n_bc) * theta;
        DType a_coef = 1.0e+3;
        for ( int i = 0; i < Na; i++ ) {
          DType tvar = x_a_c[i] * data_r[0];
          (*data_a[0])[i] = a_coef * std::exp( -tvar*tvar );
        }
      }
      break;
    default:
      std::cout << "Error: No initial condition is defined ... Exit ...\n" << std::endl;
      exit(1);
      break;
  }

  tdCoef->n_bc = tdCoef->n_bc * theta;

  upscale_var( *data_g[0], *data_g[0], x_c, Nx, data_r[0]);
  upscale_var( *data_h[0], *data_h[0], x_c, Nx, data_r[0]);
  upscale_var( *data_n[0], *data_n[0], x_c, Nx, data_r[0]);

  update_vel_dim(0);

  FILE* ot_file = fopen( tdInput->radiusFile, "w" );
  fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\t%10.20f\n", 
      tsol, data_r[0], (*data_v[0])[Nx],
      -((*data_a[0])[Nx]-(*data_a[0])[Nx-1])/dx/data_r[0]);
  fclose(ot_file);

  return;
}

void TDSolveFVM::copy_data(const int nt_fr, const int nt_to) {
  *data_g[nt_to]  = *data_g[nt_fr];
  *data_h[nt_to]  = *data_h[nt_fr];
  *data_n[nt_to]  = *data_n[nt_fr];
  *data_a[nt_to]  = *data_a[nt_fr];
  *data_v[nt_to]  = *data_v[nt_fr];
  
  // OPN model
  if ( tdSetup->prob_id == 1 ) {
    *data_p1[nt_to] = *data_p1[nt_fr];
    *data_p2[nt_to] = *data_p2[nt_fr];
  }

  data_r[nt_to]  =  data_r[nt_fr];
  return;
}

void TDSolveFVM::upscale_var(
    HVEC& data_fr, HVEC& data_to, DType*& x_ptr, 
    int Nsize, DType radius) {
  for ( int i = 0; i < Nsize; i++ ) {
    DType loc_r = x_ptr[i] * radius;
    data_to[i] = data_fr[i] * loc_r * loc_r;
  }
  return;
}

void TDSolveFVM::downscale_var(
    HVEC& data_fr, HVEC& data_to, DType*& x_ptr, 
    int Nsize, DType radius) {
  for ( int i = 0; i < Nsize; i++ ) {
    DType loc_r = x_ptr[i] * radius;
    data_to[i] = data_fr[i] / loc_r / loc_r;
  }
  return;
}

void TDSolveFVM::calc_dt(const int nt) {
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
      ((*data_a[nt])[i]-(*data_a[nt])[i-1])
      /dx/data_r[nt]/data_r[nt];
    DType v_max = max( fabs(m_vel), fabs(loc_vel-m_vel));
    v_max = max( v_max, fabs(loc_vel-m_vel+a_vel) );
    dtn_trial = min( dtn_trial, tdSetup->cfl * dx / v_max );
  }

  if ( iter == 0 )
    dtn = min( dtn_trial, dtn );
  else
    dtn = min( dtn_old * tdSetup->dTgrow, dtn_trial );
  dtn_old = dtn;

  if ( use_tfreq ) {
    reach_output_now = reach_output_next;
    if ( t_output - tsol <= dtn ) {
      dtn = t_output - tsol;
      reach_output_next = true;
      t_output += tdSetup->Tfreq;
    }
    else
      reach_output_next = false;
  }

  if ( tdSetup->Tterm - tsol < dtn ) {
    dtn = tdSetup->Tterm - tsol;
    reach_term_next = true;
  }

  if ( iter >= tdSetup->maxit )
    reach_term_now = true;

  return;
}

void TDSolveFVM::output(const int nt) {
  FILE* ot_file = fopen( tdInput->radiusFile, "a" );
  fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\t%10.20f\n", tsol, data_r[nt], (*data_v[nt])[Nx],
      -((*data_a[nt])[Nx]-(*data_a[nt])[Nx-1])/dx/data_r[nt]);
  fclose(ot_file);

  if ( !reach_term_now ) {
    if ( use_tfreq ) {
      if ( !reach_output_now ) return;
    }
    else {
      if ( iter % tdSetup->freq != 0 ) return;
    }
  }
  //if ( !reach_term_now && iter % tdSetup->freq != 0 ) return;
  ot_file = fopen( tdInput->outputFile, "a" );
  fprintf( ot_file, "%10.20f\t%10.20f\t%d\t%d\n", tsol, dtn, Nx, Nx+1 );

  downscale_var( *data_g[nt], *data_g[nt], x_c, Nx, data_r[nt] );
  downscale_var( *data_h[nt], *data_h[nt], x_c, Nx, data_r[nt] );
  downscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );

  for ( int i = 0; i < Nx; i++ )
    fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\n",
        x_c[i] * data_r[nt], (*data_g[nt])[i], (*data_h[nt])[i],
        (*data_n[nt])[i], (*data_a[nt])[i] );

  upscale_var( *data_g[nt], *data_g[nt], x_c, Nx, data_r[nt] );
  upscale_var( *data_h[nt], *data_h[nt], x_c, Nx, data_r[nt] );
  upscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );

  for ( int i = 0; i < Nx+1; i++ )
    fprintf( ot_file, "%10.20f\t%10.20f\n", x_f[i] * data_r[nt], (*data_v[nt])[i] );
  fclose(ot_file);

  return;
}

void TDSolveFVM::update_vel(const int nt) {

  if ( !tdSetup->nonDim ) {
    update_vel_dim( nt );
    return;
  }
  
  switch ( tdSetup->prob_id ) {
    case 0:
      switch ( tdSetup->vel ) {
        case 0: // 2nd-order integral form
          update_vel_int(nt);
          //update_vel_int_test(nt);
          break;
        case 1:
          update_vel_int_mod(nt);
          break;
        default:
          std::cout << "Warning: The velocity scheme not defined, "
            << "using integral form instead\n" << std::endl;
          update_vel_int(nt);
          break;
      }
      break;
    case 1: // OPN model
      update_vel_prbid1(nt);
      break;
    default:
      std::cout << "Error: Problem id not defined.\n" << std::endl;
      exit(1);
      break;
  }

  return;
}

// OPN model -- velocity
void TDSolveFVM::update_vel_prbid1(const int nt) {
  DType lambda = tdCoef->opn_lambda;
  DType delta1 = tdCoef->delta1;
  DType delta2 = tdCoef->delta2;
  DType b      = tdCoef->b;
  DType alpha  = tdCoef->alpha;
      
  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0
  // First do the integration of x^2R(lambda*P2*G/b+P2 - delta1*H - delta2*M)
  for ( int i = 0; i < Nx; i++ ) {
    (*data_v[nt])[i+1] = (*data_v[nt])[i] + dx * ( lambda * (*data_g[nt])[i]*(*data_p2[nt])[i]/(b+(*data_p2[nt])[i]) - delta1 * (*data_h[nt])[i]-delta2*(*data_n[nt])[i] ) / data_r[nt];//G,H and M have not been upscaled
  } 
  downscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );
  // Next divide by x^2 and remove the contribution from N and dA
  for ( int i = 1; i < Nx; i++ ) {
    (*data_v[nt])[i] = (*data_v[nt])[i]/x_f[i]/x_f[i] - (alpha/data_r[nt]) * 
      0.5 * ( (*data_n[nt])[i-1] + (*data_n[nt])[i] ) * 
      ( ( (*data_a[nt])[i]-(*data_a[nt])[i-1]) / dx );
  }
  // Right boundary
  {
    (*data_v[nt])[Nx] = (*data_v[nt])[Nx]/x_f[Nx]/x_f[Nx] - (alpha/data_r[nt]) * 
      0.5 * ( (*data_n[nt])[Nx-1] + tdCoef->n_bc ) * 
      ( ( (*data_a[nt])[Nx]-(*data_a[nt])[Nx-1]) / dx );
  } 
  upscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] ); 
  return;
}


void TDSolveFVM::update_vel_int(const int nt) {

  DType lambda = tdCoef->lambda;
  DType delta  = tdCoef->delta;
  DType alpha  = tdCoef->alpha;

  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0

  // First do the integration of x^2R(lambda G - delta H)
  for ( int i = 0; i < Nx; i++ ) {
    (*data_v[nt])[i+1] = (*data_v[nt])[i] + dx * ( lambda * (*data_g[nt])[i] - delta * (*data_h[nt])[i] ) / data_r[nt];
  }

  downscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );

  // Next divide by x^2 and remove the contribution from N and dA
  for ( int i = 1; i < Nx; i++ ) {
    (*data_v[nt])[i] = (*data_v[nt])[i]/x_f[i]/x_f[i] - (alpha/data_r[nt]) * 
      0.5 * ( (*data_n[nt])[i-1] + (*data_n[nt])[i] ) * 
      ( ( (*data_a[nt])[i]-(*data_a[nt])[i-1]) / dx );
  }
  // Right boundary
  {
    (*data_v[nt])[Nx] = (*data_v[nt])[Nx]/x_f[Nx]/x_f[Nx] - (alpha/data_r[nt]) * 
      0.5 * ( (*data_n[nt])[Nx-1] + tdCoef->n_bc ) *
      ( ( (*data_a[nt])[Nx]-(*data_a[nt])[Nx-1]) / dx );
  }

  upscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );

  return;
}

void TDSolveFVM::update_vel_int_mod(const int nt) {

  DType lambda = tdCoef->lambda;
  DType delta  = tdCoef->delta;
  DType alpha  = tdCoef->alpha;

  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0

  // First do the integration of x^2R(lambda G - delta H)
  for ( int i = 0; i < Nx; i++ ) {
    (*data_v[nt])[i+1] = (*data_v[nt])[i] + dx * ( lambda * (*data_g[nt])[i] - delta * (*data_h[nt])[i] ) / data_r[nt];
  }

  // Next divide by x^2 and remove the contribution from N and dA
  for ( int i = 1; i < Nx; i++ ) {
    (*data_v[nt])[i] = (*data_v[nt])[i] - (alpha/data_r[nt]/data_r[nt]/data_r[nt]) * 
      0.5 * ( (*data_n[nt])[i-1] + (*data_n[nt])[i] ) * 
      ( ( (*data_a[nt])[i]-(*data_a[nt])[i-1]) / dx );
    (*data_v[nt])[i] = (*data_v[nt])[i]/x_f[i]/x_f[i];
  }
  // Right boundary
  {
    (*data_v[nt])[Nx] = (*data_v[nt])[Nx] - (alpha/data_r[nt]/data_r[nt]/data_r[nt]) *
      0.5 * ( (*data_n[nt])[Nx-1] + 
          tdCoef->n_bc * (x_f[Nx]+0.5*dx) * (x_f[Nx]+0.5*dx) * data_r[nt] * data_r[nt] ) *
      ( ( (*data_a[nt])[Nx]-(*data_a[nt])[Nx-1]) / dx );
    (*data_v[nt])[Nx] = (*data_v[nt])[Nx]/x_f[Nx]/x_f[Nx];
  }

  return;
}

void TDSolveFVM::update_vel_int_test(const int nt) {

  DType lambda = tdCoef->lambda;
  DType delta  = tdCoef->delta;
  DType alpha  = tdCoef->alpha;

  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0

  // First do the integration of x^2R(lambda G - delta H)
  for ( int i = 0; i < Nx; i++ ) {
    int j = i+1;
    (*data_v[nt])[j] = (*data_v[nt])[i] + dx * ( lambda * (*data_g[nt])[i] - delta * (*data_h[nt])[i] ) / data_r[nt];
  }

  downscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );

  // Next divide by x^2 and remove the contribution from N and dA
  for ( int i = 1; i < Nx; i++ ) {
    (*data_v[nt])[i] = (*data_v[nt])[i]/x_f[i]/x_f[i] - (alpha/data_r[nt]) * 
      (*data_n[nt])[N_dir[i]] * 
      ( ( (*data_a[nt])[i]-(*data_a[nt])[i-1]) / dx );
  }
  // Right boundary
  {
    (*data_v[nt])[Nx] = (*data_v[nt])[Nx]/x_f[Nx]/x_f[Nx] - (alpha/data_r[nt]) * 
      ( ( N_dir[Nx]<Nx ) ? (*data_n[nt])[Nx-1] : tdCoef->n_bc ) *
      ( ( (*data_a[nt])[Nx]-(*data_a[nt])[Nx-1]) / dx );
  }

  upscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );

  return;
}

void TDSolveFVM::update_sol(const int fr, const int to) {

  if ( !tdSetup->nonDim ) {
    update_sol_dim( fr, to );
    return;
  }

  update_rad( fr, to );

  switch( tdSetup->scell ) {
    case 0: // Explicit first-order upwind
      update_cells_upw( fr, to );
      //update_cells_upw_test( fr, to );
      break;
    case 1: // Implicit first-order upwind
      update_cells_upw_be( fr, to );
      break;
    case 2: // Explicit second-order MUSCL
      update_cells_muscl( fr, to );
      break;
    case 3: // Implicit second-order MUSCL
      update_cells_muscl_be( fr, to );
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
    case 1: // Explicit second-order muscl for advection and implicit cd for diffusion
      update_chem_muscl_imex( fr, to );
      break;
    default:
      std::cout << "Warning: The spatial scheme for chem not defined, "
        << "using upw/backward on the primitive equation instead\n" << std::endl;
      update_chem_upw_be( fr, to );
      break;
  }

  update_vel( to );
  return;
}

void TDSolveFVM::update_sol_dim(const int fr, const int to) {
  update_rad( fr, to );

  update_cells_upw( fr, to );

  update_chem_upw_be_dim( fr, to );

  update_vel_dim( to );

  return;
}

void TDSolveFVM::update_vel_dim(const int nt) {

  DType lambda = tdCoef->lambda;
  DType delta  = tdCoef->delta;
  DType alpha  = tdCoef->alpha;

  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0

  for ( int i = 0; i < Nx; i++ ) {
    DType denom = (*data_g[nt])[i] + (*data_h[nt])[i] + (*data_n[nt])[i];
    (*data_v[nt])[i+1] = (*data_v[nt])[i] + dx * x_c[i] * x_c[i] * data_r[nt] *
      ( lambda * (*data_g[nt])[i] - delta * (*data_h[nt])[i] ) / denom;
    DType val_r = 0.0;
    DType val_l = 0.0;

    if ( i == 0 )
      val_l = 0.0;
    else
      val_l = (alpha/data_r[nt]/data_r[nt]/data_r[nt]) * 
        0.5 * ( (*data_n[nt])[i-1] + (*data_n[nt])[i] ) * 
        ( ( (*data_a[nt])[i]-(*data_a[nt])[i-1]) / dx );

    if ( i == Nx-1 )
      val_r = (alpha/data_r[nt]/data_r[nt]/data_r[nt]) * 
        0.5 * ( (*data_n[nt])[i] + 
            tdCoef->n_bc * (x_f[Nx]+0.5*dx) * (x_f[Nx]+0.5*dx) * data_r[nt] * data_r[nt] ) * 
        ( ( (*data_a[nt])[i+1]-(*data_a[nt])[i]) / dx );
    else
      val_r = (alpha/data_r[nt]/data_r[nt]/data_r[nt]) * 
        0.5 * ( (*data_n[nt])[i] + (*data_n[nt])[i+1] ) * 
        ( ( (*data_a[nt])[i+1]-(*data_a[nt])[i]) / dx );

    (*data_v[nt])[i+1] = (*data_v[nt])[i+1] - ( val_r - val_l ) 
      * x_c[i] * x_c[i] * data_r[nt] * data_r[nt] / denom;
  }

  for ( int i = 1; i < Nx+1; i++ )
    (*data_v[nt])[i] = (*data_v[nt])[i]/x_f[i]/x_f[i];

  return;
}

void TDSolveFVM::update_cells_upw(const int fr, const int to) {
  DType lambda  = tdCoef->lambda;
  DType delta   = tdCoef->delta;
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType alpha   = tdCoef->alpha;
  DType rho_mns = tdCoef->rho_mns;
  DType rho_pls = tdCoef->rho_pls;
  DType n_star  = tdCoef->n_star;

  DType dror = (*data_v[fr])[Nx]/data_r[fr];
  DType aorr = alpha/data_r[fr]/data_r[fr];

  // Source terms
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_mu = mu - mu_slp * x_c[i];
    DType rho = ( (*data_n[fr])[i] < n_star ) ? rho_mns : rho_pls;
    //cout << "rho: " << rho << std::endl;
    (*data_g[to])[i] = (*data_g[fr])[i] + dtn * (
        lambda * (*data_g[fr])[i] - loc_mu * (*data_g[fr])[i]
        - dror * (*data_g[fr])[i] );
    (*data_h[to])[i] = (*data_h[fr])[i] + dtn * (
        loc_mu * (*data_g[fr])[i] - delta * (*data_h[fr])[i]
        - dror * (*data_h[fr])[i] );
    (*data_n[to])[i] = (*data_n[fr])[i] + dtn * (
        - rho * (*data_n[fr])[i] - dror * (*data_n[fr])[i] );
  }

  // Fluxes
  { // Left boundary
    (*flux_g[0])[0] = 0.0;
    (*flux_h[0])[0] = 0.0;
    (*flux_n[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_vel = (*data_v[fr])[i]/data_r[fr] - x_f[i]*dror;
    (*flux_g[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_g[fr])[i-1] : (*data_g[fr])[i] );
    (*flux_h[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_h[fr])[i-1] : (*data_h[fr])[i] );
    loc_vel += aorr * (((*data_a[fr])[i]-(*data_a[fr])[i-1])/dx);
    (*flux_n[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_n[fr])[i-1] : (*data_n[fr])[i] );
  }
  { // Right boundary
    (*flux_g[0])[Nx] = 0.0;
    (*flux_h[0])[Nx] = 0.0;
    DType loc_vel = aorr * (((*data_a[fr])[Nx]-(*data_a[fr])[Nx-1])/dx);
    (*flux_n[0])[Nx] = loc_vel * (
        (loc_vel>=0.0) ? (*data_n[fr])[Nx-1] : 
        tdCoef->n_bc * (x_f[Nx]+0.5*dx) * (x_f[Nx]+0.5*dx) * data_r[fr] * data_r[fr] ) ;
  }

  for ( int i = 0; i < Nx; i++ ) {
    (*data_g[to])[i] -= ( (*flux_g[0])[i+1] - (*flux_g[0])[i] ) * (dtn/dx);
    (*data_h[to])[i] -= ( (*flux_h[0])[i+1] - (*flux_h[0])[i] ) * (dtn/dx);
    (*data_n[to])[i] -= ( (*flux_n[0])[i+1] - (*flux_n[0])[i] ) * (dtn/dx);
  }

  return;
}

void TDSolveFVM::update_cells_upw_test(const int fr, const int to) {
  DType lambda = tdCoef->lambda;
  DType delta  = tdCoef->delta;
  DType mu     = tdCoef->mu;
  DType mu_slp = tdCoef->mu_slp;
  DType alpha  = tdCoef->alpha;

  DType dror = (*data_v[fr])[Nx]/data_r[fr];
  DType aorr = alpha/data_r[fr]/data_r[fr];

  // Source terms
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_mu = mu - mu_slp * x_c[i];
    (*data_g[to])[i] = (*data_g[fr])[i] + dtn * (
        lambda * (*data_g[fr])[i] - loc_mu * (*data_g[fr])[i]
        - dror * (*data_g[fr])[i] );
    (*data_h[to])[i] = (*data_h[fr])[i] + dtn * (
        loc_mu * (*data_g[fr])[i] - delta * (*data_h[fr])[i]
        - dror * (*data_h[fr])[i] );
    (*data_n[to])[i] = (*data_n[fr])[i] + dtn * (
        - dror * (*data_n[fr])[i] );
  }

  // Fluxes
  { // Left boundary
    (*flux_g[0])[0] = 0.0;
    (*flux_h[0])[0] = 0.0;
    (*flux_n[0])[0] = 0.0;
    N_dir[0] = 0;
  }
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_vel = (*data_v[fr])[i]/data_r[fr] - x_f[i]*dror;
    (*flux_g[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_g[fr])[i-1] : (*data_g[fr])[i] );
    (*flux_h[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_h[fr])[i-1] : (*data_h[fr])[i] );
    (*flux_n[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_n[fr])[i-1] : (*data_n[fr])[i] );
    loc_vel = aorr * (((*data_a[fr])[i]-(*data_a[fr])[i-1])/dx);
    //(*flux_n[0])[i] += loc_vel * 0.5 * ((*data_n[fr])[i-1] + (*data_n[fr])[i]);
    (*flux_n[0])[i] += loc_vel * (
        (loc_vel>=0.0) ? (*data_n[fr])[i-1] : (*data_n[fr])[i] );
    N_dir[i] = (loc_vel>=0.0) ? i-1 : i;
  }
  { // Right boundary
    (*flux_g[0])[Nx] = 0.0;
    (*flux_h[0])[Nx] = 0.0;
    (*flux_n[0])[Nx] = 0.0;
    DType loc_vel = aorr * (((*data_a[fr])[Nx]-(*data_a[fr])[Nx-1])/dx);
    //(*flux_n[0])[Nx] += loc_vel * 0.5 * ((*data_n[fr])[Nx-1] +
    //    tdCoef->n_bc * (x_f[Nx]+0.5*dx) * (x_f[Nx]+0.5*dx) * data_r[fr] * data_r[fr] );
    (*flux_n[0])[Nx] += loc_vel * (
        (loc_vel>=0.0) ? (*data_n[fr])[Nx-1] :
        tdCoef->n_bc * (x_f[Nx]+0.5*dx) * (x_f[Nx]+0.5*dx) * data_r[fr] * data_r[fr] );
    N_dir[Nx] = (loc_vel>=0.0) ? Nx-1 : Nx;
  }

  for ( int i = 0; i < Nx; i++ ) {
    (*data_g[to])[i] -= ( (*flux_g[0])[i+1] - (*flux_g[0])[i] ) * (dtn/dx);
    (*data_h[to])[i] -= ( (*flux_h[0])[i+1] - (*flux_h[0])[i] ) * (dtn/dx);
    (*data_n[to])[i] -= ( (*flux_n[0])[i+1] - (*flux_n[0])[i] ) * (dtn/dx);
  }

  return;
}

void TDSolveFVM::update_cells_upw_be(const int fr, const int to) {
  DType lambda = tdCoef->lambda;
  DType delta  = tdCoef->delta;
  DType mu     = tdCoef->mu;
  DType mu_slp = tdCoef->mu_slp;
  DType alpha  = tdCoef->alpha;
  DType rho_mns = tdCoef->rho_mns;
  DType rho_pls = tdCoef->rho_pls;
  DType n_star  = tdCoef->n_star;

  DType dror = (*data_v[fr])[Nx]/data_r[fr];
  DType aorr = alpha/data_r[fr]/data_r[fr];

  VComp comp(1,0.0);
  Entry elem(1,0.0);

  // First data_g
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  // RHS
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_g[fr])[i];
    rhs->set(i,comp);
  }

  // Matrix
  // Left boundary
  {
    DType loc_v_r = (*data_v[fr])[1]/data_r[fr] - x_f[1] * dror;
    DType coef_adv_c = 0.0;
    DType coef_adv_r = 0.0;

    if ( loc_v_r >= 0.0 )
      coef_adv_c += dtn*loc_v_r/dx;
    else 
      coef_adv_r += dtn*loc_v_r/dx;

    DType loc_mu = mu - mu_slp * x_c[0];
    elem = 1.0 + coef_adv_c - dtn*( lambda - loc_mu - dror );
    mat->set(0,0,elem);

    elem = coef_adv_r;
    mat->set(0,1,elem);
  }
  // Interior node
  for ( int i = 1; i < Nx-1; i++ ) {
    DType loc_v_r = (*data_v[fr])[i+1]/data_r[fr] - x_f[i+1] * dror;
    DType loc_v_l = (*data_v[fr])[i]/data_r[fr] - x_f[i] * dror;
    DType coef_adv_c = 0.0;
    DType coef_adv_r = 0.0;
    DType coef_adv_l = 0.0;

    if ( loc_v_r >= 0.0 )
      coef_adv_c += dtn*loc_v_r/dx;
    else
      coef_adv_r += dtn*loc_v_r/dx;

    if ( loc_v_l >= 0.0 )
      coef_adv_l -= dtn*loc_v_l/dx;
    else
      coef_adv_c -= dtn*loc_v_l/dx;

    elem = coef_adv_l;
    mat->set(i,i-1,elem);

    DType loc_mu = mu - mu_slp * x_c[i];
    elem = 1.0 + coef_adv_c - dtn*( lambda - loc_mu - dror );
    mat->set(i,i,elem);

    elem = coef_adv_r;
    mat->set(i,i+1,elem);
  }
  // Right boundary
  {
    DType loc_v_l = (*data_v[fr])[Nx-1]/data_r[fr] - x_f[Nx-1] * dror;
    
    DType coef_adv_c = 0.0;
    DType coef_adv_l = 0.0;

    if ( loc_v_l >= 0.0 )
      coef_adv_l -= dtn*loc_v_l/dx;
    else
      coef_adv_c -= dtn*loc_v_l/dx;

    elem = coef_adv_l;
    mat->set(Nx-1,Nx-2,elem);

    DType loc_mu = mu - mu_slp * x_c[Nx-1];
    elem = 1.0 + coef_adv_c - dtn*( lambda - loc_mu - dror );
    mat->set(Nx-1,Nx-1,elem);
  }

  // Initialize solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_g[fr])[i];
    sol->set(i,comp);
  }

  linSolve->solve();

  // Fetch solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = sol->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_g[to])[i] = vsol[0];
  }


  // First data_h
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  // RHS
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_mu = mu - mu_slp * x_c[i];
    comp = (*data_h[fr])[i] + dtn * loc_mu * (*data_g[to])[i];
    rhs->set(i,comp);
  }

  // Matrix
  // Left boundary
  {
    DType loc_v_r = (*data_v[fr])[1]/data_r[fr] - x_f[1] * dror;
    DType coef_adv_c = 0.0;
    DType coef_adv_r = 0.0;

    if ( loc_v_r >= 0.0 )
      coef_adv_c += dtn*loc_v_r/dx;
    else 
      coef_adv_r += dtn*loc_v_r/dx;

    elem = 1.0 + coef_adv_c - dtn*( - delta - dror );
    mat->set(0,0,elem);

    elem = coef_adv_r;
    mat->set(0,1,elem);
  }
  // Interior node
  for ( int i = 1; i < Nx-1; i++ ) {
    DType loc_v_r = (*data_v[fr])[i+1]/data_r[fr] - x_f[i+1] * dror;
    DType loc_v_l = (*data_v[fr])[i]/data_r[fr] - x_f[i] * dror;
    DType coef_adv_c = 0.0;
    DType coef_adv_r = 0.0;
    DType coef_adv_l = 0.0;

    if ( loc_v_r >= 0.0 )
      coef_adv_c += dtn*loc_v_r/dx;
    else
      coef_adv_r += dtn*loc_v_r/dx;

    if ( loc_v_l >= 0.0 )
      coef_adv_l -= dtn*loc_v_l/dx;
    else
      coef_adv_c -= dtn*loc_v_l/dx;

    elem = coef_adv_l;
    mat->set(i,i-1,elem);

    elem = 1.0 + coef_adv_c - dtn*( - delta - dror );
    mat->set(i,i,elem);

    elem = coef_adv_r;
    mat->set(i,i+1,elem);
  }
  // Right boundary
  {
    DType loc_v_l = (*data_v[fr])[Nx-1]/data_r[fr] - x_f[Nx-1] * dror;
    
    DType coef_adv_c = 0.0;
    DType coef_adv_l = 0.0;

    if ( loc_v_l >= 0.0 )
      coef_adv_l -= dtn*loc_v_l/dx;
    else
      coef_adv_c -= dtn*loc_v_l/dx;

    elem = coef_adv_l;
    mat->set(Nx-1,Nx-2,elem);

    elem = 1.0 + coef_adv_c - dtn*( - delta - dror );
    mat->set(Nx-1,Nx-1,elem);
  }

  // Initialize solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_h[fr])[i];
    sol->set(i,comp);
  }

  linSolve->solve();

  // Fetch solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = sol->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_h[to])[i] = vsol[0];
  }


  // First data_n
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  // RHS
  for ( int i = 0; i < Nx-1; i++ ) {
    comp = (*data_n[fr])[i];
    rhs->set(i,comp);
  }
  // Right boundary
  {
    DType loc_v_r = aorr * (((*data_a[fr])[Nx]-(*data_a[fr])[Nx-1])/dx);
    
    DType coef_adv_r = 0.0;

    if ( loc_v_r < 0.0 )
      coef_adv_r += dtn*loc_v_r/dx;

    comp = (*data_n[fr])[Nx-1] - coef_adv_r * tdCoef->n_bc * 
      (x_f[Nx]+0.5*dx) * (x_f[Nx]+0.5*dx) * data_r[fr] * data_r[fr];
    rhs->set(Nx-1,comp);
  }

  // Matrix
  // Left boundary
  {
    DType loc_v_r = (*data_v[fr])[1]/data_r[fr] - x_f[1] * dror
      + aorr * (((*data_a[fr])[1]-(*data_a[fr])[0])/dx);
    DType coef_adv_c = 0.0;
    DType coef_adv_r = 0.0;

    if ( loc_v_r >= 0.0 )
      coef_adv_c += dtn*loc_v_r/dx;
    else 
      coef_adv_r += dtn*loc_v_r/dx;

    DType rho = ( (*data_n[fr])[0] < n_star ) ? rho_mns : rho_pls;
    elem = 1.0 + coef_adv_c - dtn*( - rho - dror );
    mat->set(0,0,elem);

    elem = coef_adv_r;
    mat->set(0,1,elem);
  }
  // Interior node
  for ( int i = 1; i < Nx-1; i++ ) {
    DType loc_v_r = (*data_v[fr])[i+1]/data_r[fr] - x_f[i+1] * dror
      + aorr * (((*data_a[fr])[i+1]-(*data_a[fr])[i])/dx);
    DType loc_v_l = (*data_v[fr])[i]/data_r[fr] - x_f[i] * dror
      + aorr * (((*data_a[fr])[i]-(*data_a[fr])[i-1])/dx);
    DType coef_adv_c = 0.0;
    DType coef_adv_r = 0.0;
    DType coef_adv_l = 0.0;

    if ( loc_v_r >= 0.0 )
      coef_adv_c += dtn*loc_v_r/dx;
    else
      coef_adv_r += dtn*loc_v_r/dx;

    if ( loc_v_l >= 0.0 )
      coef_adv_l -= dtn*loc_v_l/dx;
    else
      coef_adv_c -= dtn*loc_v_l/dx;

    elem = coef_adv_l;
    mat->set(i,i-1,elem);

    DType rho = ( (*data_n[fr])[i] < n_star ) ? rho_mns : rho_pls;
    elem = 1.0 + coef_adv_c - dtn*( - rho - dror );
    mat->set(i,i,elem);

    elem = coef_adv_r;
    mat->set(i,i+1,elem);
  }
  // Right boundary
  {
    DType loc_v_l = (*data_v[fr])[Nx-1]/data_r[fr] - x_f[Nx-1] * dror
      + aorr * (((*data_a[fr])[Nx-1]-(*data_a[fr])[Nx-2])/dx);
    DType loc_v_r = aorr * (((*data_a[fr])[Nx]-(*data_a[fr])[Nx-1])/dx);
    
    DType coef_adv_c = 0.0;
    DType coef_adv_l = 0.0;
    DType coef_adv_r = 0.0;

    if ( loc_v_r >= 0.0 )
      coef_adv_c += dtn*loc_v_r/dx;
    else
      coef_adv_r += dtn*loc_v_r/dx;

    if ( loc_v_l >= 0.0 )
      coef_adv_l -= dtn*loc_v_l/dx;
    else
      coef_adv_c -= dtn*loc_v_l/dx;

    elem = coef_adv_l;
    mat->set(Nx-1,Nx-2,elem);

    DType rho = ( (*data_n[fr])[Nx-1] < n_star ) ? rho_mns : rho_pls;
    elem = 1.0 + coef_adv_c - dtn*( - rho - dror );
    mat->set(Nx-1,Nx-1,elem);
  }

  // Initialize solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_n[fr])[i];
    sol->set(i,comp);
  }

  linSolve->solve();

  // Fetch solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = sol->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_n[to])[i] = vsol[0];
  }

  return;
}

void TDSolveFVM::update_cells_muscl(const int fr, const int to) {
  DType lambda  = tdCoef->lambda;
  DType delta   = tdCoef->delta;
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType alpha   = tdCoef->alpha;
  DType rho_mns = tdCoef->rho_mns;
  DType rho_pls = tdCoef->rho_pls;
  DType n_star  = tdCoef->n_star;

  DType dror = (*data_v[fr])[Nx]/data_r[fr];
  DType aorr = alpha/data_r[fr]/data_r[fr];

  // Source terms
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_mu = mu - mu_slp * x_c[i];
    DType rho = ( (*data_n[fr])[i] < n_star ) ? rho_mns : rho_pls;
    //cout << "rho: " << rho << std::endl;
    (*data_g[to])[i] = (*data_g[fr])[i] + dtn * (
        lambda * (*data_g[fr])[i] - loc_mu * (*data_g[fr])[i]
        - dror * (*data_g[fr])[i] );
    (*data_h[to])[i] = (*data_h[fr])[i] + dtn * (
        loc_mu * (*data_g[fr])[i] - delta * (*data_h[fr])[i]
        - dror * (*data_h[fr])[i] );
    (*data_n[to])[i] = (*data_n[fr])[i] + dtn * (
        - rho * (*data_n[fr])[i] - dror * (*data_n[fr])[i] );
  }

  // Slopes
  { // Left boundary
    (*slope_g[0])[0] = 0.0;
    (*slope_h[0])[0] = 0.0;
    (*slope_n[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx-1; i++ ) {
    (*slope_g[0])[i] = Alg_LimSlope_Minmod( (*data_g[fr])[i], (*data_g[fr])[i-1], (*data_g[fr])[i+1], dx );
    (*slope_h[0])[i] = Alg_LimSlope_Minmod( (*data_h[fr])[i], (*data_h[fr])[i-1], (*data_h[fr])[i+1], dx );
    (*slope_n[0])[i] = Alg_LimSlope_Minmod( (*data_n[fr])[i], (*data_n[fr])[i-1], (*data_n[fr])[i+1], dx );
  }
  { // Right boundary
    (*slope_g[0])[Nx-1] = 0.0;
    (*slope_h[0])[Nx-1] = 0.0;
    if ( tdSetup->bc_lim )
      (*slope_n[0])[Nx-1] = Alg_LimSlope_Minmod( (*data_n[fr])[Nx-1], (*data_n[fr])[Nx-2], 
          tdCoef->n_bc * (x_f[Nx]+0.5*dx) * (x_f[Nx]+0.5*dx) * data_r[fr] * data_r[fr] , dx );
    else
      (*slope_n[0])[Nx-1] = 0.0; 
  }

  // Fluxes
  { // Left boundary
    (*flux_g[0])[0] = 0.0;
    (*flux_h[0])[0] = 0.0;
    (*flux_n[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_vel = (*data_v[fr])[i]/data_r[fr] - x_f[i]*dror;
    (*flux_g[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_g[fr])[i-1]+0.5*dx*(*slope_g[0])[i-1] : 
                         (*data_g[fr])[i]-0.5*dx*(*slope_g[0])[i] );
    (*flux_h[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_h[fr])[i-1]+0.5*dx*(*slope_h[0])[i-1] : 
                         (*data_h[fr])[i]-0.5*dx*(*slope_h[0])[i] );
    loc_vel += aorr * (((*data_a[fr])[i]-(*data_a[fr])[i-1])/dx);
    (*flux_n[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_n[fr])[i-1]+0.5*dx*(*slope_n[0])[i-1] : 
                         (*data_n[fr])[i]-0.5*dx*(*slope_n[0])[i] );
  }
  { // Right boundary
    (*flux_g[0])[Nx] = 0.0;
    (*flux_h[0])[Nx] = 0.0;
    DType loc_vel = aorr * (((*data_a[fr])[Nx]-(*data_a[fr])[Nx-1])/dx);
    (*flux_n[0])[Nx] = loc_vel * (
        (loc_vel>=0.0) ? (*data_n[fr])[Nx-1]+0.5*dx*(*slope_n[0])[Nx-1] : 
        tdCoef->n_bc * x_f[Nx] * x_f[Nx] * data_r[fr] * data_r[fr] ) ;
  }

  for ( int i = 0; i < Nx; i++ ) {
    (*data_g[to])[i] -= ( (*flux_g[0])[i+1] - (*flux_g[0])[i] ) * (dtn/dx);
    (*data_h[to])[i] -= ( (*flux_h[0])[i+1] - (*flux_h[0])[i] ) * (dtn/dx);
    (*data_n[to])[i] -= ( (*flux_n[0])[i+1] - (*flux_n[0])[i] ) * (dtn/dx);
  }

  return;
}

void TDSolveFVM::update_cells_muscl_be(const int fr, const int to) {
}
 
void TDSolveFVM::update_chem_upw_be(const int fr, const int to) {
  rhs_a->reset( 0.0 );
  mat_a->reset( 0.0 );
  sol_a->reset( 0.0 );

  DType m     = tdCoef->m;
  DType D     = tdCoef->D;
  DType beta  = tdCoef->beta;
  DType gamma = tdCoef->gamma;
  DType theta = tdCoef->theta;

  DType dror = (*data_v[fr])[Nx]/data_r[fr];

  // RHS
  VComp comp(1,0.0);
  for ( int i = 0; i < Na; i++ ) {
    DType loc_r_fr = x_a_c[i] * data_r[fr];
    DType loc_r_to = x_a_c[i] * data_r[to];
    comp = loc_r_fr * loc_r_fr * (*data_a[fr])[i];
    if ( i < Nx ) { // G is active 
      comp += dtn * m * (*data_g[to])[i] / 
        (beta/theta+(*data_g[to])[i]/loc_r_to/loc_r_to);
    }
    comp /= loc_r_to*loc_r_to;
    rhs_a->set(i,comp);
  }

  // Matrix
  Entry elem(1,0.0);

  // Left boundary
  {
    DType loc_v_r = - x_a_f[1] * dror;
    DType coef_adv_0 = 0.0;
    DType coef_adv_1 = 0.0;

    DType loc_r_0 = x_a_c[0] * data_r[to];
    DType loc_r_1 = x_a_c[1] * data_r[to];

    if ( loc_v_r >= 0.0 )
      coef_adv_0 += (dtn*loc_v_r/dx)*loc_r_0*loc_r_0;
    else
      coef_adv_1 += (dtn*loc_v_r/dx)*loc_r_1*loc_r_1;

    DType coef_diff_1 = dtn * D * x_a_f[1] * x_a_f[1] / dx / dx;

    elem = 1.0 + dtn * ( gamma + dror ) + 
      (coef_diff_1 + coef_adv_0) / loc_r_0 / loc_r_0;
    mat_a->set(0,0,elem);

    elem = (- coef_diff_1 + coef_adv_1) / loc_r_0 / loc_r_0;
    mat_a->set(0,1,elem);
  }
  // Interior cells
  for ( int i = 1; i < Na-1; i++ ) {
    DType loc_v_r = - x_a_f[i+1] * dror;
    DType loc_v_l = - x_a_f[i] * dror;
    DType coef_adv_c = 0.0;
    DType coef_adv_r = 0.0;
    DType coef_adv_l = 0.0;

    DType loc_r_r = x_a_c[i+1] * data_r[to];
    DType loc_r_c = x_a_c[i] * data_r[to];
    DType loc_r_l = x_a_c[i-1] * data_r[to];

    if ( loc_v_r >= 0.0 ) 
      coef_adv_c += (dtn*loc_v_r/dx)*loc_r_c*loc_r_c;
    else
      coef_adv_r += (dtn*loc_v_r/dx)*loc_r_r*loc_r_r;

    if ( loc_v_l >= 0.0 )
      coef_adv_l -= (dtn*loc_v_l/dx)*loc_r_l*loc_r_l;
    else
      coef_adv_c -= (dtn*loc_v_l/dx)*loc_r_c*loc_r_c;

    DType coef_diff_r = dtn * D * x_a_f[i+1] * x_a_f[i+1] / dx / dx;
    DType coef_diff_l = dtn * D * x_a_f[i] * x_a_f[i] / dx / dx;

    elem = (- coef_diff_l + coef_adv_l) / loc_r_c / loc_r_c;
    mat_a->set(i,i-1,elem);

    elem = 1.0 + dtn * ( gamma + dror ) + 
      (coef_diff_l + coef_diff_r + coef_adv_c) / loc_r_c / loc_r_c;
    mat_a->set(i,i,elem);

    elem = (- coef_diff_r + coef_adv_r) / loc_r_c / loc_r_c;
    mat_a->set(i,i+1,elem);
  }
  // Right boundary
  {
    DType loc_v_l = - x_a_f[Na-1] * dror;
    DType coef_adv_l = 0.0;
    DType coef_adv_c = 0.0;

    DType loc_r_l = x_a_c[Na-2] * data_r[to];
    DType loc_r_c = x_a_c[Na-1] * data_r[to];

    if ( loc_v_l >= 0.0 )
      coef_adv_l -= (dtn*loc_v_l/dx)*loc_r_l*loc_r_l;
    else
      coef_adv_c -= (dtn*loc_v_l/dx)*loc_r_c*loc_r_c;

    DType coef_diff_r = dtn * D * x_a_f[Na] * x_a_f[Na] / dx / dx;
    DType coef_diff_l = dtn * D * x_a_f[Na-1] * x_a_f[Na-1] / dx / dx;

    elem = (- coef_diff_l + coef_adv_l) / loc_r_c / loc_r_c;
    mat_a->set(Na-1,Na-2,elem);

    elem = 1.0 + dtn * ( gamma + dror ) + 
      (coef_diff_l + coef_diff_r + coef_adv_c) / loc_r_c / loc_r_c;
    mat_a->set(Na-1,Na-1,elem);
  }

  // Initialize solution
  for ( int i = 0; i < Na; i++ ) {
    comp = (*data_a[fr])[i];
    sol_a->set(i,comp);
  }

  linSolve_a->solve();

  // Fetch solution
  for ( int i = 0; i < Na; i++ ) {
    comp = sol_a->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_a[to])[i] = vsol[0];
  }

  return;
}

void TDSolveFVM::update_chem_muscl_imex(const int fr, const int to) {
  rhs_a->reset( 0.0 );
  mat_a->reset( 0.0 );
  sol_a->reset( 0.0 );

  DType m     = tdCoef->m;
  DType D     = tdCoef->D;
  DType beta  = tdCoef->beta;
  DType gamma = tdCoef->gamma;
  DType theta = tdCoef->theta;

  DType dror = (*data_v[fr])[Nx]/data_r[fr];

  // Slopes
  (*slope_a[0])[0] = 0.0;
  for ( int i = 1; i < Na-1; i++ ) {
    DType loc_r_l = x_a_c[i-1] * data_r[fr];
    DType loc_r_c = x_a_c[i] * data_r[fr];
    DType loc_r_r = x_a_c[i+1] * data_r[fr];
    (*slope_a[0])[i] = Alg_LimSlope_Minmod(
        loc_r_c * loc_r_c * (*data_a[fr])[i], loc_r_l * loc_r_l * (*data_a[fr])[i-1],
        loc_r_r * loc_r_r * (*data_a[fr])[i+1], dx );
  }
  (*slope_a[0])[Na-1] = 0.0;

  // Fluxes
  (*flux_a[0])[0] = 0.0;
  for ( int i = 1; i < Na; i++ ) {
    DType loc_r_l = x_a_c[i-1] * data_r[fr];
    DType loc_r_r = x_a_c[i]   * data_r[fr];
    DType loc_vel = - x_a_f[i] * dror;
    (*flux_a[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? loc_r_l*loc_r_l*(*data_a[fr])[i-1]+0.5*dx*(*slope_a[0])[i-1] :
                         loc_r_r*loc_r_r*(*data_a[fr])[i]-0.5*dx*(*slope_a[0])[i] );
  }
  (*flux_a[0])[Na] = 0.0; // This is an approximation, truly valid only when R is always increasing

  // RHS
  VComp comp(1,0.0);
  for ( int i = 0; i < Na; i++ ) {
    DType loc_r_fr = x_a_c[i] * data_r[fr];
    DType loc_r_to = x_a_c[i] * data_r[to];
    comp = loc_r_fr * loc_r_fr * (*data_a[fr])[i];
    if ( i < Nx ) { // G is active 
      comp += dtn * m * (*data_g[to])[i] / 
        (beta/theta+(*data_g[to])[i]/loc_r_to/loc_r_to);
    }
    // Advection part
    comp -= ((*flux_a[0])[i+1] - (*flux_a[0])[i]) * (dtn/dx);
    comp /= loc_r_to*loc_r_to;
    rhs_a->set(i,comp);
  }

  // Matrix
  Entry elem(1,0.0);

  // Left boundary
  {
    DType loc_r_0 = x_a_c[0] * data_r[to];

    DType coef_diff_1 = dtn * D * x_a_f[1] * x_a_f[1] / dx / dx;

    elem = 1.0 + dtn * ( gamma + dror ) + coef_diff_1 / loc_r_0 / loc_r_0;
    mat_a->set(0,0,elem);

    elem = - coef_diff_1 / loc_r_0 / loc_r_0;
    mat_a->set(0,1,elem);
  }
  // Interior cells
  for ( int i = 1; i < Na-1; i++ ) {
    DType loc_r_c = x_a_c[i] * data_r[to];

    DType coef_diff_r = dtn * D * x_a_f[i+1] * x_a_f[i+1] / dx / dx;
    DType coef_diff_l = dtn * D * x_a_f[i] * x_a_f[i] / dx / dx;

    elem = - coef_diff_l / loc_r_c / loc_r_c;
    mat_a->set(i,i-1,elem);

    elem = 1.0 + dtn * ( gamma + dror ) + 
      (coef_diff_l + coef_diff_r ) / loc_r_c / loc_r_c;
    mat_a->set(i,i,elem);

    elem = - coef_diff_r / loc_r_c / loc_r_c;
    mat_a->set(i,i+1,elem);
  }
  // Right boundary
  {
    DType loc_r_c = x_a_c[Na-1] * data_r[to];

    DType coef_diff_r = dtn * D * x_a_f[Na] * x_a_f[Na] / dx / dx;
    DType coef_diff_l = dtn * D * x_a_f[Na-1] * x_a_f[Na-1] / dx / dx;

    elem = - coef_diff_l / loc_r_c / loc_r_c;
    mat_a->set(Na-1,Na-2,elem);

    elem = 1.0 + dtn * ( gamma + dror ) + 
      (coef_diff_l + coef_diff_r) / loc_r_c / loc_r_c;
    mat_a->set(Na-1,Na-1,elem);
  }

  // Initialize solution
  for ( int i = 0; i < Na; i++ ) {
    comp = (*data_a[fr])[i];
    sol_a->set(i,comp);
  }

  linSolve_a->solve();

  // Fetch solution
  for ( int i = 0; i < Na; i++ ) {
    comp = sol_a->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_a[to])[i] = vsol[0];
  }

  return;
}

void TDSolveFVM::update_chem_upw_be_dim(const int fr, const int to) {
  rhs_a->reset( 0.0 );
  mat_a->reset( 0.0 );
  sol_a->reset( 0.0 );

  DType m     = tdCoef->m;
  DType D     = tdCoef->D;
  DType beta  = tdCoef->beta;
  DType gamma = tdCoef->gamma;

  DType dror = (*data_v[fr])[Nx]/data_r[fr];

  // RHS
  VComp comp(1,0.0);
  for ( int i = 0; i < Na; i++ ) {
    DType loc_r_fr = x_a_c[i] * data_r[fr];
    DType loc_r_to = x_a_c[i] * data_r[to];
    comp = loc_r_fr * loc_r_fr * (*data_a[fr])[i];
    if ( i < Nx ) { // G is active 
      comp += dtn * m * (*data_g[to])[i] / 
        (beta+(*data_g[to])[i]/loc_r_to/loc_r_to);
    }
    comp /= loc_r_to*loc_r_to;
    rhs_a->set(i,comp);
  }

  // Matrix
  Entry elem(1,0.0);

  // Left boundary
  {
    DType loc_v_r = - x_a_f[1] * dror;
    DType coef_adv_0 = 0.0;
    DType coef_adv_1 = 0.0;

    DType loc_r_0 = x_a_c[0] * data_r[to];
    DType loc_r_1 = x_a_c[1] * data_r[to];

    if ( loc_v_r >= 0.0 )
      coef_adv_0 += (dtn*loc_v_r/dx)*loc_r_0*loc_r_0;
    else
      coef_adv_1 += (dtn*loc_v_r/dx)*loc_r_1*loc_r_1;

    DType coef_diff_1 = dtn * D * x_a_f[1] * x_a_f[1] / dx / dx;

    elem = 1.0 + dtn * ( gamma + dror ) + 
      (coef_diff_1 + coef_adv_0) / loc_r_0 / loc_r_0;
    mat_a->set(0,0,elem);

    elem = (- coef_diff_1 + coef_adv_1) / loc_r_0 / loc_r_0;
    mat_a->set(0,1,elem);
  }
  // Interior cells
  for ( int i = 1; i < Na-1; i++ ) {
    DType loc_v_r = - x_a_f[i+1] * dror;
    DType loc_v_l = - x_a_f[i] * dror;
    DType coef_adv_c = 0.0;
    DType coef_adv_r = 0.0;
    DType coef_adv_l = 0.0;

    DType loc_r_r = x_a_c[i+1] * data_r[to];
    DType loc_r_c = x_a_c[i] * data_r[to];
    DType loc_r_l = x_a_c[i-1] * data_r[to];

    if ( loc_v_r >= 0.0 ) 
      coef_adv_c += (dtn*loc_v_r/dx)*loc_r_c*loc_r_c;
    else
      coef_adv_r += (dtn*loc_v_r/dx)*loc_r_r*loc_r_r;

    if ( loc_v_l >= 0.0 )
      coef_adv_l -= (dtn*loc_v_l/dx)*loc_r_l*loc_r_l;
    else
      coef_adv_c -= (dtn*loc_v_l/dx)*loc_r_c*loc_r_c;

    DType coef_diff_r = dtn * D * x_a_f[i+1] * x_a_f[i+1] / dx / dx;
    DType coef_diff_l = dtn * D * x_a_f[i] * x_a_f[i] / dx / dx;

    elem = (- coef_diff_l + coef_adv_l) / loc_r_c / loc_r_c;
    mat_a->set(i,i-1,elem);

    elem = 1.0 + dtn * ( gamma + dror ) + 
      (coef_diff_l + coef_diff_r + coef_adv_c) / loc_r_c / loc_r_c;
    mat_a->set(i,i,elem);

    elem = (- coef_diff_r + coef_adv_r) / loc_r_c / loc_r_c;
    mat_a->set(i,i+1,elem);
  }
  // Right boundary
  {
    DType loc_v_l = - x_a_f[Na-1] * dror;
    DType coef_adv_l = 0.0;
    DType coef_adv_c = 0.0;

    DType loc_r_l = x_a_c[Na-2] * data_r[to];
    DType loc_r_c = x_a_c[Na-1] * data_r[to];

    if ( loc_v_l >= 0.0 )
      coef_adv_l -= (dtn*loc_v_l/dx)*loc_r_l*loc_r_l;
    else
      coef_adv_c -= (dtn*loc_v_l/dx)*loc_r_c*loc_r_c;

    DType coef_diff_r = dtn * D * x_a_f[Na] * x_a_f[Na] / dx / dx;
    DType coef_diff_l = dtn * D * x_a_f[Na-1] * x_a_f[Na-1] / dx / dx;

    elem = (- coef_diff_l + coef_adv_l) / loc_r_c / loc_r_c;
    mat_a->set(Na-1,Na-2,elem);

    elem = 1.0 + dtn * ( gamma + dror ) + 
      (coef_diff_l + coef_diff_r + coef_adv_c) / loc_r_c / loc_r_c;
    mat_a->set(Na-1,Na-1,elem);
  }

  // Initialize solution
  for ( int i = 0; i < Na; i++ ) {
    comp = (*data_a[fr])[i];
    sol_a->set(i,comp);
  }

  linSolve_a->solve();

  // Fetch solution
  for ( int i = 0; i < Na; i++ ) {
    comp = sol_a->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_a[to])[i] = vsol[0];
  }

  return;
}


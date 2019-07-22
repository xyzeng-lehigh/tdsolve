#include <cmath>
#include "TDSolveFVM.h"
#include "Alloc_Utils.h"
#include "TDInput.h"
#include "Alg_Solve.h"
#include "Alg_Limiters.h"
using namespace std;

TDSolveFVM::TDSolveFVM(TDInput* tdInput_ptr) : 
  TDSolve(tdInput_ptr) {

  slope_g = NULL;
  slope_h = NULL;
  slope_n = NULL;
  flux_g  = NULL;
  flux_h  = NULL;
  flux_n  = NULL;
  flux_a  = NULL;
  slope_a = NULL;
  flux_p1  = NULL;
  slope_p1 = NULL;
  flux_p2  = NULL;
  slope_p2 = NULL;
  N_dir    = NULL;

  if ( !(tdSetup->silence) ) {
    FILE* ot_file = fopen( tdInput->outputFile, "w" );
    if ( tdSetup->prob_id < 0 ) // This is only relevant for FVMMod
      fprintf( ot_file, "Time, time step size, radius, num. point / x_c, G, H, N, 1-G-H-N /newline x_f, V\n");
    if ( tdSetup->prob_id == 0 )
      fprintf( ot_file, "Time, time step size, radius, num. point / x_c, G, H, N, A /newline x_f, V\n");
    if ( tdSetup->prob_id == 1 )
      fprintf( ot_file, "Time, time step size, radius, num. point / x_c, G, H, M, A, P1, P2 /newline x_f, V\n");
    fclose(ot_file);
  }
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
    erase_array( flux_p1, 1 );
    erase_array( flux_p2, 1 );
    erase_array( slope_p1, 1 );
    erase_array( slope_p2, 1 );
  }
}

void TDSolveFVM::setup(void) {

  setupLinSolve( mat, prc, sol, rhs, res, linSolve, Nx );
  setupLinSolve( mat_a, prc_a, sol_a, rhs_a, res_a, linSolve_a, Na );

  // OPN model
  if ( tdSetup->prob_id == 1 ) 
    setupLinSolve( mat_p, prc_p, sol_p, rhs_p, res_p, linSolve_p, Nx );

  return;
}

void TDSolveFVM::set_ic(void) {

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
    case 2:
      {
        (*data_g[0]) = 1.0 - tdCoef->n_bc;
        (*data_h[0]) = 0.0;
        (*data_n[0]) = tdCoef->n_bc;
        DType r0 = tdCoef->R_0;
        DType lambda = tdCoef->lambda;
        DType alpha  = tdCoef->alpha;
        DType A0     = 10.0 * r0 * r0 * lambda / 6.0 / alpha;
        for ( int i = 0; i < Nx; i++ ) {
          DType loc_r = x_a_c[i] * data_r[0];
          (*data_a[0])[i] = A0 - loc_r * loc_r * lambda / 6.0 / alpha;
        }
        DType AR     = A0 - r0 * r0 * lambda / 6.0 / alpha;
        DType slope  = r0 * lambda / 3.0 / alpha / AR;
        for ( int i = Nx; i< Na; i++ ) {
          DType loc_r = x_a_c[i] * data_r[0];
          (*data_a[0])[i] = AR * exp(-slope*(loc_r-r0));
        }
      }
      break;
    default:
      std::cout << "Error: No initial condition is defined ... Exit ...\n" << std::endl;
      exit(1);
      break;
  }

  upscale_var( *data_g[0], *data_g[0], x_c, Nx, data_r[0] );
  upscale_var( *data_h[0], *data_h[0], x_c, Nx, data_r[0] );
  upscale_var( *data_n[0], *data_n[0], x_c, Nx, data_r[0] );
  upscale_var( *data_a[0], *data_a[0], x_a_c, Na, data_r[0] );

  if ( tdSetup->prob_id == 1 ) {
    (*data_p1[0]) = 0.0;
    (*data_p2[0]) = 0.0;
  }

  this->update_vel(0, 1, 0.0);

  N_dir[0] = 0;
  for ( int i = 1; i < Nx+1; i++ )
    N_dir[i] = i-1;

  if ( tdSetup->frame == 1 )
  {
    rprime = (*data_v[0])[Nx];
  }

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

void TDSolveFVM::calc_dt(const int nt) {
  if ( iter != 0 )
    tsol += dtn;

  if ( data_r[nt] <= tdSetup->Rdead || data_r[nt] >= tdSetup->Rterm )
    reach_term_now = true;

  if ( reach_term_next )
    reach_term_now = true;

  downscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );

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

  upscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );

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
  
  if ( tdSetup->silence ) return;

  switch ( tdSetup->prob_id ) {
    case 0:
    case -1: // Test
    case -2: // Test
    case -3: // Test
    case -4: // Test
      output_pdgf( nt );
      break;
    case 1:
      output_opn( nt );
      break;
    default:
      std::cout << "Error: Problem id not defined.\n" << std::endl;
      exit(1);
      break;
  }

  return;
}

void TDSolveFVM::output_pdgf(const int nt) {
  if ( tdSetup->silence ) return;

  downscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );

  FILE* ot_file = fopen( tdInput->radiusFile, "a" );
  if ( tdSetup->frame == 2 && tdSetup->prob_id < 0 )
    this->print_custom_vel_u( ot_file, nt );
  else
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
  fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\t%d\t%d\n", tsol, dtn, data_r[nt], Nx, Nx+1 );

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
  upscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );

  for ( int i = 0; i < Nx+1; i++ )
    fprintf( ot_file, "%10.20f\t%10.20f\n", x_f[i] * data_r[nt], (*data_v[nt])[i] );
  fclose(ot_file);

  return;
}

void TDSolveFVM::output_opn(const int nt) {
  if ( tdSetup->silence ) return;

  downscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );

  FILE* ot_file = fopen( tdInput->radiusFile, "a" );
  if ( tdSetup->frame == 2 && tdSetup->prob_id < 0 )
    this->print_custom_vel_u( ot_file, nt );
  else
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
  fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\t%d\t%d\n", tsol, dtn, data_r[nt], Nx, Nx+1 );

  downscale_var( *data_g[nt], *data_g[nt], x_c, Nx, data_r[nt] );
  downscale_var( *data_h[nt], *data_h[nt], x_c, Nx, data_r[nt] );
  downscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );
  downscale_var( *data_p1[nt], *data_p1[nt], x_c, Nx, data_r[nt] );
  downscale_var( *data_p2[nt], *data_p2[nt], x_c, Nx, data_r[nt] );

  for ( int i = 0; i < Nx; i++ )
    fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\n",
        x_c[i] * data_r[nt], (*data_g[nt])[i], (*data_h[nt])[i],
        (*data_n[nt])[i], (*data_a[nt])[i], (*data_p1[nt])[i], (*data_p2[nt])[i] );

  upscale_var( *data_g[nt], *data_g[nt], x_c, Nx, data_r[nt] );
  upscale_var( *data_h[nt], *data_h[nt], x_c, Nx, data_r[nt] );
  upscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );
  upscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );
  upscale_var( *data_p1[nt], *data_p1[nt], x_c, Nx, data_r[nt] );
  upscale_var( *data_p2[nt], *data_p2[nt], x_c, Nx, data_r[nt] );

  for ( int i = 0; i < Nx+1; i++ )
    fprintf( ot_file, "%10.20f\t%10.20f\n", x_f[i] * data_r[nt], (*data_v[nt])[i] );
  fclose(ot_file);

  return;
}

void TDSolveFVM::alloc(void) {

  switch (tdSetup->tmeth) {
    case 0:
      Ndata = 2;
      break;
    case 1:
      Ndata = 3;
      break;
    case 2:
      Ndata = 5;
      break;
  }

  alloc_array_full( data_g, Nx, Ndata );
  alloc_array_full( data_h, Nx, Ndata );
  alloc_array_full( data_n, Nx, Ndata ); // M
  alloc_array_full( data_v, Nx+1, Ndata );

  data_r = new DType[Ndata];
  
  alloc_array_full( data_a, Na, Ndata ); 

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
    alloc_array_full (data_p1, Nx, Ndata );
    alloc_array_full (data_p2, Nx, Ndata );

    alloc_array_full (flux_p1, Nx+1, 1);
    alloc_array_full (flux_p2, Nx+1, 1);

    alloc_array_full(slope_p1, Nx, 1);
    alloc_array_full(slope_p2, Nx, 1);
  }
  else {
    data_p1 = NULL;
    data_p2 = NULL;
    flux_p1 = NULL;
    flux_p2 = NULL;
    slope_p1 = NULL;
    slope_p2 = NULL;
  }

  alloc_array_full( N_dir, Nx+1 );
  return;
}

void TDSolveFVM::update_vel(const int nt, const int nt_dummy, const DType shift_dt, const bool last_iter) {

  switch ( tdSetup->prob_id ) {
    case 0:
      switch ( tdSetup->vel ) {
        case 0: // 2nd-order integral form
          update_vel_int(nt);
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
  }

  return;
}

DType TDSolveFVM::get_nperc(const int nt) {
  DType tot_g = 0.0;
  DType tot_h = 0.0;
  DType tot_n = 0.0;

  for ( int i = 0; i < Nx; i++ ) {
    DType r_loc = data_r[nt] * x_c[i];
    tot_g += 4.0 * PI * r_loc * r_loc * (*data_g[nt])[i];
    tot_h += 4.0 * PI * r_loc * r_loc * (*data_h[nt])[i];
    tot_n += 4.0 * PI * r_loc * r_loc * (*data_n[nt])[i];
  }

  return tot_n / ( tot_g + tot_h + tot_n );
}

// OPN model -- velocity
void TDSolveFVM::update_vel_prbid1(const int nt) {
  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;
  DType b       = tdCoef->b;
  DType b_prime = tdCoef->b_prime;
  DType alpha   = tdCoef->alpha;
  DType omega   = tdCoef->omega;
      
  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0
  // First do the integration of
  // x^2R[lambda(1+b'P2/(b+P2))*G-delta_h*H-delta_n*N]
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    DType loc_r = x_c[i] * data_r[nt];
    DType logistic_term = 1.0; //= ( tdSetup->logistic ) ? 1.0 : 0.0;
    if ( tdSetup->logistic == 1 )
      logistic_term = ( 1.0 - (*data_g[nt])[i] / loc_r / loc_r );
    if ( tdSetup->logistic == 2 )
      logistic_term = max( 0.0, ( 1.0 - (*data_g[nt])[i] / loc_r / loc_r ) );
    (*data_v[nt])[i+1] = (*data_v[nt])[i] + dx * (
        loc_lambda * ( 1.0 + (b_prime * (*data_p2[nt])[i])/(b + (*data_p2[nt])[i]) ) * (*data_g[nt])[i] * 
        logistic_term - 
        //( 1.0 - logistic_term * (*data_g[nt])[i] / loc_r / loc_r ) - 
        delta_h * (*data_h[nt])[i] - delta_n * (*data_n[nt])[i] -
        omega * (*data_g[nt])[i] * (*data_n[nt])[i] / loc_r / loc_r ) / data_r[nt]; // G, H, and M have not been upscaled
  }

  downscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );
  downscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );

  // Next divide by x^2 and remove the contribution from N and dA
  for ( int i = 1; i < Nx; i++ ) {
    (*data_v[nt])[i] = (*data_v[nt])[i]/x_f[i]/x_f[i] - (alpha/data_r[nt]) * 
      0.5 * ( (*data_n[nt])[i-1] + (*data_n[nt])[i] ) * 
      ( ( (*data_a[nt])[i]-(*data_a[nt])[i-1]) / dx ); // x = eta
  }
  // Right boundary
  {
    (*data_v[nt])[Nx] = (*data_v[nt])[Nx]/x_f[Nx]/x_f[Nx] - (alpha/data_r[nt]) * 
      0.5 * ( (*data_n[nt])[Nx-1] + tdCoef->n_bc ) * 
      ( ( (*data_a[nt])[Nx]-(*data_a[nt])[Nx-1]) / dx );
  } 

  upscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] ); 
  upscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );

  return;
}

void TDSolveFVM::update_vel_int(const int nt) {

  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;
  DType alpha   = tdCoef->alpha;

  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0

  // First do the integration of x^2R(lambda G - delta_h H)
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    (*data_v[nt])[i+1] = (*data_v[nt])[i] + dx * ( loc_lambda * (*data_g[nt])[i] - delta_h * (*data_h[nt])[i] - delta_n * (*data_n[nt])[i] ) / data_r[nt];
  }

  downscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );
  downscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );

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
  upscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );

  return;
}

void TDSolveFVM::update_vel_int_mod(const int nt) {

  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType alpha   = tdCoef->alpha;

  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0

  // First do the integration of x^2R(lambda G - delta_h H)
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    (*data_v[nt])[i+1] = (*data_v[nt])[i] + dx * ( loc_lambda * (*data_g[nt])[i] - delta_h * (*data_h[nt])[i] ) / data_r[nt];
  }

  downscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );

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

  upscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );

  return;
}

bool TDSolveFVM::update_sol(const int fr, const int to, 
    const DType shift_dt_fr, const DType shift_dt_to) {

  update_vel( fr, 0, shift_dt_fr );
  update_rad( fr, to );

  // Both original and OPN
  switch ( tdSetup->scell ) {
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

  if ( tdSetup->prob_id == 1 ) // OPN
    update_cells_opn_add( fr, to );

  switch( tdSetup->prob_id ) {
    case 0:
      switch ( tdSetup->schem ) {
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
      break;
    case 1:
      switch ( tdSetup->schem ) {
        case 0: // Implicit first-order upwind for advection and cd for diffusion
          // Update A
          update_chem_upw_be( fr, to ); 
          // Update P1 and P2
          update_chem_opn_upw_be( tsol+shift_dt_fr*dtn, fr, to, tdCoef->D_1, tdCoef->alpha_1, 
              tdCoef->beta_1, tdCoef->gamma_1, data_n, data_p1, tdCoef->treatP1, tdCoef->treatP1BC );
          update_chem_opn_upw_be( tsol+shift_dt_fr*dtn, fr, to, tdCoef->D_2, tdCoef->alpha_2, 
              tdCoef->beta_2, tdCoef->gamma_2, data_g, data_p2 );
          break;
        case 1: // Explicit second-order muscl for advection and implicit cd for diffusion
          // Update A
          update_chem_muscl_imex( fr, to ); 
          // Update P1 and P2
          update_chem_opn_muscl_imex( tsol+shift_dt_fr*dtn, fr, to, tdCoef->D_1, tdCoef->alpha_1, 
              tdCoef->beta_1, tdCoef->gamma_1, data_n, data_p1, flux_p1, slope_p1, 
              tdCoef->treatP1, tdCoef->treatP1BC );
          update_chem_opn_muscl_imex( tsol+shift_dt_fr*dtn, fr, to, tdCoef->D_2, tdCoef->alpha_2, 
              tdCoef->beta_2, tdCoef->gamma_2, data_g, data_p2, flux_p2, slope_p2 );
          break;
        default:
          std::cout << "Warning: The spatial scheme for chem not defined, "
            << "using upw/backward on the primitive equation instead\n" << std::endl;
          // Update A
          update_chem_upw_be( fr, to );
          // Update P1 and P2
          update_chem_opn_upw_be( tsol+shift_dt_fr*dtn, fr, to, tdCoef->D_1, tdCoef->alpha_1, 
              tdCoef->beta_1, tdCoef->gamma_1, data_n, data_p1, tdCoef->treatP1, tdCoef->treatP1BC );
          update_chem_opn_upw_be( tsol+shift_dt_fr*dtn, fr, to, tdCoef->D_2, tdCoef->alpha_2, 
              tdCoef->beta_2, tdCoef->gamma_2, data_g, data_p2 );
          break;
      }
      break;
  }

  update_vel( to, 0, shift_dt_to );
  return true;
}

void TDSolveFVM::update_cells_upw(const int fr, const int to) {
  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType alpha   = tdCoef->alpha;
  DType delta_n = tdCoef->delta_n;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

  DType dror = (*data_v[fr])[Nx]/data_r[fr];
  DType aorr = alpha/data_r[fr]/data_r[fr];

  // Source terms
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[i],p) ); 
    (*data_g[to])[i] = (*data_g[fr])[i] + dtn * (
        loc_lambda * (*data_g[fr])[i] - loc_mu * (*data_g[fr])[i]
        - dror * (*data_g[fr])[i] );
    (*data_h[to])[i] = (*data_h[fr])[i] + dtn * (
        loc_mu * (*data_g[fr])[i] - delta_h * (*data_h[fr])[i]
        - dror * (*data_h[fr])[i] );
    (*data_n[to])[i] = (*data_n[fr])[i] + dtn * (
        - delta_n * (*data_n[fr])[i] - dror * (*data_n[fr])[i] );
  }
  
  downscale_var( *data_a[fr], *data_a[fr], x_a_c, Na, data_r[fr] );

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

  upscale_var( *data_a[fr], *data_a[fr], x_a_c, Na, data_r[fr] );

  for ( int i = 0; i < Nx; i++ ) {
    (*data_g[to])[i] -= ( (*flux_g[0])[i+1] - (*flux_g[0])[i] ) * (dtn/dx);
    (*data_h[to])[i] -= ( (*flux_h[0])[i+1] - (*flux_h[0])[i] ) * (dtn/dx);
    (*data_n[to])[i] -= ( (*flux_n[0])[i+1] - (*flux_n[0])[i] ) * (dtn/dx);
  }

  return;
}

void TDSolveFVM::update_cells_upw_test(const int fr, const int to) {
  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType alpha   = tdCoef->alpha;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

  DType dror = (*data_v[fr])[Nx]/data_r[fr];
  DType aorr = alpha/data_r[fr]/data_r[fr];

  // Source terms
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[i],p) ); 
    (*data_g[to])[i] = (*data_g[fr])[i] + dtn * (
        loc_lambda * (*data_g[fr])[i] - loc_mu * (*data_g[fr])[i]
        - dror * (*data_g[fr])[i] );
    (*data_h[to])[i] = (*data_h[fr])[i] + dtn * (
        loc_mu * (*data_g[fr])[i] - delta_h * (*data_h[fr])[i]
        - dror * (*data_h[fr])[i] );
    (*data_n[to])[i] = (*data_n[fr])[i] + dtn * (
        - dror * (*data_n[fr])[i] );
  }

  downscale_var( *data_a[fr], *data_a[fr], x_a_c, Na, data_r[fr] );

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

  upscale_var( *data_a[fr], *data_a[fr], x_a_c, Na, data_r[fr] );

  for ( int i = 0; i < Nx; i++ ) {
    (*data_g[to])[i] -= ( (*flux_g[0])[i+1] - (*flux_g[0])[i] ) * (dtn/dx);
    (*data_h[to])[i] -= ( (*flux_h[0])[i+1] - (*flux_h[0])[i] ) * (dtn/dx);
    (*data_n[to])[i] -= ( (*flux_n[0])[i+1] - (*flux_n[0])[i] ) * (dtn/dx);
  }

  return;
}

void TDSolveFVM::update_cells_upw_be(const int fr, const int to) {
  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType alpha   = tdCoef->alpha;
  DType delta_n = tdCoef->delta_n;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

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

    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[0],q);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[0],p) ); 
    elem = 1.0 + coef_adv_c - dtn*( loc_lambda - loc_mu - dror );
    mat->set(0,0,elem);

    elem = coef_adv_r;
    mat->set(0,1,elem);
  }
  // Interior node
  for ( int i = 1; i < Nx-1; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
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

    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[i],p) ); 
    elem = 1.0 + coef_adv_c - dtn*( loc_lambda - loc_mu - dror );
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

    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[Nx-1],q);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[Nx-1],p) ); 
    elem = 1.0 + coef_adv_c - dtn*( loc_lambda - loc_mu - dror );
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


  // Second data_h
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  // RHS
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[i],p) ); 
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

    elem = 1.0 + coef_adv_c - dtn*( - delta_h - dror );
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

    elem = 1.0 + coef_adv_c - dtn*( - delta_h - dror );
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

    elem = 1.0 + coef_adv_c - dtn*( - delta_h - dror );
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


  // Last data_n
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  downscale_var( *data_a[fr], *data_a[fr], x_a_c, Na, data_r[fr] );

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

    elem = 1.0 + coef_adv_c - dtn*( - delta_n - dror );
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

    elem = 1.0 + coef_adv_c - dtn*( - delta_n - dror );
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

    elem = 1.0 + coef_adv_c - dtn*( - delta_n - dror );
    mat->set(Nx-1,Nx-1,elem);
  }

  upscale_var( *data_a[fr], *data_a[fr], x_a_c, Na, data_r[fr] );

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
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType alpha   = tdCoef->alpha;
  DType delta_n = tdCoef->delta_n;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

  DType dror = (*data_v[fr])[Nx]/data_r[fr];
  DType aorr = alpha/data_r[fr]/data_r[fr];

  // Source terms
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[i],p) ); 
    (*data_g[to])[i] = (*data_g[fr])[i] + dtn * (
        loc_lambda * (*data_g[fr])[i] - loc_mu * (*data_g[fr])[i]
        - dror * (*data_g[fr])[i] );
    (*data_h[to])[i] = (*data_h[fr])[i] + dtn * (
        loc_mu * (*data_g[fr])[i] - delta_h * (*data_h[fr])[i]
        - dror * (*data_h[fr])[i] );
    (*data_n[to])[i] = (*data_n[fr])[i] + dtn * (
        - delta_n * (*data_n[fr])[i] - dror * (*data_n[fr])[i] );
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
    (*slope_n[0])[Nx-1] = 0.0; 
  }

  downscale_var( *data_a[fr], *data_a[fr], x_a_c, Na, data_r[fr] );

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

  upscale_var( *data_a[fr], *data_a[fr], x_a_c, Na, data_r[fr] );

  for ( int i = 0; i < Nx; i++ ) {
    (*data_g[to])[i] -= ( (*flux_g[0])[i+1] - (*flux_g[0])[i] ) * (dtn/dx);
    (*data_h[to])[i] -= ( (*flux_h[0])[i+1] - (*flux_h[0])[i] ) * (dtn/dx);
    (*data_n[to])[i] -= ( (*flux_n[0])[i+1] - (*flux_n[0])[i] ) * (dtn/dx);
  }

  return;
}

void TDSolveFVM::update_cells_muscl_be(const int fr, const int to) {
}

// OPN model

void TDSolveFVM::update_cells_opn_add(const int fr, const int to,
    const bool only_res ) {
  DType a       = tdCoef->a;
  DType b       = tdCoef->b;
  DType a_prime = tdCoef->a_prime;
  DType b_prime = tdCoef->b_prime;
  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType omega   = tdCoef->omega;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

  DType coef_space = only_res ? 1.0 : dtn;

  downscale_var( *data_p1[fr], *data_p1[fr], x_c, Nx, data_r[fr] );
  downscale_var( *data_p2[fr], *data_p2[fr], x_c, Nx, data_r[fr] );

  // Additional source terms
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[i],p) ); 
    (*data_g[to])[i] += coef_space * ( 
        loc_lambda * b_prime * (*data_p2[fr])[i] / ( b + (*data_p2[fr])[i] ) * (*data_g[fr])[i] - 
        loc_mu * a_prime * (*data_p1[fr])[i] / ( a + (*data_p1[fr])[i] ) * (*data_g[fr])[i] );
    (*data_h[to])[i] += coef_space * (
        loc_mu * a_prime * (*data_p1[fr])[i] / ( a + (*data_p1[fr])[i] ) * (*data_g[fr])[i] );
  }

  if ( omega != 0.0 )
    for ( int i = 0; i < Nx; i++ ) {
      DType loc_r = x_c[i] * data_r[fr];
      (*data_g[to])[i] += coef_space * (
          - omega * (*data_g[fr])[i] * (*data_n[fr])[i] / loc_r / loc_r );
    }

  if ( tdSetup->logistic == 1 ) 
    for ( int i = 0; i < Nx; i++ ) {
      DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
      DType loc_r = x_c[i] * data_r[fr];
      (*data_g[to])[i] += coef_space * loc_lambda * ( 1.0 + b_prime * (*data_p2[fr])[i] / ( b + (*data_p2[fr])[i] ) ) * 
        (*data_g[fr])[i] * ( - (*data_g[fr])[i] / loc_r / loc_r );
    }

  if ( tdSetup->logistic == 2 )
    for ( int i = 0; i < Nx; i++ ) {
      DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
      DType loc_r = x_c[i] * data_r[fr];
      DType logistic_term = min( 1.0, (*data_g[fr])[i] / loc_r / loc_r );
      (*data_g[to])[i] += coef_space * loc_lambda * ( 1.0 + b_prime * (*data_p2[fr])[i] / ( b + (*data_p2[fr])[i] ) ) * 
        (*data_g[fr])[i] * ( - logistic_term ); //( - (*data_g[fr])[i] / loc_r / loc_r );
    }

  upscale_var( *data_p1[fr], *data_p1[fr], x_c, Nx, data_r[fr] );
  upscale_var( *data_p2[fr], *data_p2[fr], x_c, Nx, data_r[fr] );

  return;
}

DType TDSolveFVM::calc_treatment(const DType tsol_loc) {
  std::vector<DType> t0 = tdCoef->treatP1_t0;
  DType T  = tdCoef->treatP1_T;
  DType h  = tdCoef->treatP1_h;

  if ( t0.size() == 0 )
    return 0.0;

  int pass_count = 0;

  for ( unsigned i = 0; i < t0.size(); i++ )
    if ( tsol_loc <= t0[i] )
      break;
    else
      pass_count++;

  if ( pass_count%2==0 )
    return 0.0;
  else 
    return T * (1.0+cos(2.0*PI*tsol_loc/h));
}
 
void TDSolveFVM::update_chem_upw_be(const int fr, const int to) {
  rhs_a->reset( 0.0 );
  mat_a->reset( 0.0 );
  sol_a->reset( 0.0 );

  DType alpha_a = tdCoef->alpha_a;
  DType D       = tdCoef->D;
  DType beta_a  = tdCoef->beta_a;
  DType gamma_a = tdCoef->gamma_a;
  DType theta   = tdCoef->theta;

  DType dror = (*data_v[fr])[Nx]/data_r[fr];

  downscale_var( *data_a[fr], *data_a[fr], x_a_c, Na, data_r[fr] );

  // RHS
  VComp comp(1,0.0);
  for ( int i = 0; i < Na; i++ ) {
    DType loc_r_fr = x_a_c[i] * data_r[fr];
    DType loc_r_to = x_a_c[i] * data_r[to];
    comp = loc_r_fr * loc_r_fr * (*data_a[fr])[i];
    if ( i < Nx ) { // G is active 
      comp += dtn * alpha_a * (*data_g[to])[i] / 
        (beta_a/theta+(*data_g[to])[i]/loc_r_to/loc_r_to);
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

    elem = 1.0 + dtn * ( gamma_a + dror ) + 
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

    elem = 1.0 + dtn * ( gamma_a + dror ) + 
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

    elem = 1.0 + dtn * ( gamma_a + dror ) + 
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

  upscale_var( *data_a[fr], *data_a[fr], x_a_c, Na, data_r[fr] );
  upscale_var( *data_a[to], *data_a[to], x_a_c, Na, data_r[to] );

  return;
}

// Update req accordingly for OPN
void TDSolveFVM::update_chem_muscl_imex(const int fr, const int to) {
  rhs_a->reset( 0.0 );
  mat_a->reset( 0.0 );
  sol_a->reset( 0.0 );

  DType alpha_a = tdCoef->alpha_a;
  DType D       = tdCoef->D;
  DType beta_a  = tdCoef->beta_a;
  DType gamma_a = tdCoef->gamma_a;
  DType theta   = tdCoef->theta;

  DType dror = (*data_v[fr])[Nx]/data_r[fr];

  downscale_var( *data_a[fr], *data_a[fr], x_a_c, Na, data_r[fr] );

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
      comp += dtn * alpha_a * (*data_g[to])[i] / 
        (beta_a/theta+(*data_g[to])[i]/loc_r_to/loc_r_to);
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

    elem = 1.0 + dtn * ( gamma_a + dror ) + coef_diff_1 / loc_r_0 / loc_r_0;
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

    elem = 1.0 + dtn * ( gamma_a + dror ) + 
      (coef_diff_l + coef_diff_r) / loc_r_c / loc_r_c;
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

    elem = 1.0 + dtn * ( gamma_a + dror ) + 
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

  upscale_var( *data_a[fr], *data_a[fr], x_a_c, Na, data_r[fr] );
  upscale_var( *data_a[to], *data_a[to], x_a_c, Na, data_r[to] );

  return;
}

void TDSolveFVM::update_chem_opn_upw_be(
    const DType tsol_loc, const int fr, const int to, 
    const DType coef_D, const DType coef_alpha, const DType coef_beta, 
    const DType coef_gamma, HVEC**& data_source, HVEC**& data_P,
    const bool treatment, const bool treatment_bc) {
  rhs_p->reset( 0.0 );
  mat_p->reset( 0.0 );
  sol_p->reset( 0.0 );

  DType theta = tdCoef->theta;

  DType dror = (*data_v[fr])[Nx]/data_r[fr];

  DType f_treat = 0.0;
  if ( treatment ) // Backward Euler
    f_treat = calc_treatment(tsol_loc+dtn);
  DType bc_treat = 0.0;
  if ( treatment_bc ) // Backward Euler
    bc_treat = tdCoef->treatP1_bc * calc_treatment(tsol_loc+dtn);

  downscale_var( *data_P[fr], *data_P[fr], x_c, Nx, data_r[fr] );

  // RHS
  VComp comp(1,0.0);
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_r_fr = x_c[i] * data_r[fr];
    DType loc_r_to = x_c[i] * data_r[to];
    comp = loc_r_fr * loc_r_fr * (*data_P[fr])[i];
    comp += dtn * coef_alpha * (*data_source[to])[i] / 
      (coef_beta/theta+(*data_source[to])[i]/loc_r_to/loc_r_to);
    comp += dtn * loc_r_to * loc_r_to * f_treat;
    comp /= loc_r_to * loc_r_to;
    rhs_p->set(i,comp);
  }
  if ( bc_treat != 0.0 ) { // Fix the non-hom Neumann condition
    comp = rhs_p->get(Nx-1);
    DType loc_r_to = x_c[Nx-1] * data_r[to];
    comp *= loc_r_to * loc_r_to;

    // Advection part
    DType loc_v_r = -x_f[Nx] * dror;
    if ( loc_v_r<0.0 ) {
      DType loc_r_r = (2.0*x_c[Nx-1]-x_c[Nx-2]) * data_r[to];
      comp += dtn * loc_v_r * loc_r_r * loc_r_r * data_r[to] * bc_treat;
    }

    // Diffusion part
    comp += dtn * coef_D * x_f[Nx] * x_f[Nx] * 
      ( - data_r[to] * bc_treat ) / dx;

    comp /= loc_r_to * loc_r_to;
    rhs_p->set(Nx-1,comp);
  }

  // Matrix
  Entry elem(1,0.0);

  // Left boundary dP/dx = 0 at x = 0
  {
    DType loc_v_r = - x_f[1] * dror;
    DType coef_adv_0 = 0.0;
    DType coef_adv_1 = 0.0;

    DType loc_r_0 = x_c[0] * data_r[to];
    DType loc_r_1 = x_c[1] * data_r[to];

    if ( loc_v_r >= 0.0 )
      coef_adv_0 += (dtn*loc_v_r/dx)*loc_r_0*loc_r_0;
    else
      coef_adv_1 += (dtn*loc_v_r/dx)*loc_r_1*loc_r_1;

    DType coef_diff_1 = dtn * coef_D * x_f[1] * x_f[1] / dx / dx;

    elem = 1.0 + dtn * ( coef_gamma + dror ) + 
      (coef_diff_1 + coef_adv_0) / loc_r_0 / loc_r_0;
    mat_p->set(0,0,elem);

    elem = (- coef_diff_1 + coef_adv_1) / loc_r_0 / loc_r_0;
    mat_p->set(0,1,elem);
  }
  // Interior cells
  for ( int i = 1; i < Nx-1; i++ ) {
    DType loc_v_r = - x_f[i+1] * dror;
    DType loc_v_l = - x_f[i] * dror;
    DType coef_adv_c = 0.0;
    DType coef_adv_r = 0.0;
    DType coef_adv_l = 0.0;

    DType loc_r_r = x_c[i+1] * data_r[to];
    DType loc_r_c = x_c[i] * data_r[to];
    DType loc_r_l = x_c[i-1] * data_r[to];

    if ( loc_v_r >= 0.0 )
      coef_adv_c += (dtn*loc_v_r/dx)*loc_r_c*loc_r_c;
    else
      coef_adv_r += (dtn*loc_v_r/dx)*loc_r_r*loc_r_r;

    if ( loc_v_l >= 0.0 )
      coef_adv_l -= (dtn*loc_v_l/dx)*loc_r_l*loc_r_l;
    else
      coef_adv_c -= (dtn*loc_v_l/dx)*loc_r_c*loc_r_c;

    DType coef_diff_r = dtn * coef_D * x_f[i+1] * x_f[i+1] / dx / dx;
    DType coef_diff_l = dtn * coef_D * x_f[i] * x_f[i] / dx / dx;


    elem = (- coef_diff_l + coef_adv_l ) / loc_r_c / loc_r_c;
    mat_p->set(i,i-1,elem);

    elem = 1.0 + dtn * ( coef_gamma + dror ) + 
      (coef_diff_l + coef_diff_r + coef_adv_c) / loc_r_c / loc_r_c;
    mat_p->set(i,i,elem);

    elem = (- coef_diff_r + coef_adv_r ) / loc_r_c / loc_r_c;
    mat_p->set(i,i+1,elem);
  }
  // Right boundary dP/dx = 0 at x = 1
  {
    DType loc_v_l = -x_f[Nx-1] * dror;
    DType loc_v_r = -x_f[Nx] * dror;
    DType coef_adv_l = 0.0;
    DType coef_adv_c = 0.0;

    DType loc_r_l = x_c[Nx-2] * data_r[to];
    DType loc_r_c = x_c[Nx-1] * data_r[to];
    DType loc_r_r = (2.0*x_c[Nx-1]-x_c[Nx-2]) * data_r[to];

    if ( loc_v_r >= 0.0 )
      coef_adv_c += (dtn*loc_v_r/dx)*loc_r_c*loc_r_c;
    else
      coef_adv_c += (dtn*loc_v_r/dx)*loc_r_r*loc_r_r;

    if ( loc_v_l >= 0.0 )
      coef_adv_l -= (dtn*loc_v_l/dx)*loc_r_l*loc_r_l;
    else
      coef_adv_c -= (dtn*loc_v_l/dx)*loc_r_c*loc_r_c;

    DType coef_diff_l = dtn * coef_D * x_f[Nx-1] * x_f[Nx-1] / dx / dx;

    elem = (- coef_diff_l + coef_adv_l) / loc_r_c / loc_r_c;
    mat_p->set(Nx-1,Nx-2,elem);

    elem = 1.0 + dtn * ( coef_gamma + dror ) + 
      ( coef_diff_l +coef_adv_c) / loc_r_c / loc_r_c;
    mat_p->set(Nx-1,Nx-1,elem);
  }

  // Initialize solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_P[fr])[i];
    sol_p->set(i,comp);
  }

  upscale_var( *data_P[fr], *data_P[fr], x_c, Nx, data_r[fr] );

  linSolve_p->solve();

  // Fetch solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = sol_p->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_P[to])[i] = vsol[0];
  }

  upscale_var( *data_P[to], *data_P[to], x_c, Nx, data_r[to] );

  return;
}

void TDSolveFVM::update_chem_opn_muscl_imex(
    const DType tsol_loc, const int fr, const int to, 
    const DType coef_D, const DType coef_alpha, const DType coef_beta, 
    const DType coef_gamma, HVEC**& data_source, 
    HVEC**& data_P, HVEC**& flux_P, HVEC**& slope_P,
    const bool treatment, const bool treatment_bc) {
  rhs_p->reset( 0.0 );
  mat_p->reset( 0.0 );
  sol_p->reset( 0.0 );

  DType theta = tdCoef->theta;

  DType dror = (*data_v[fr])[Nx]/data_r[fr];

  DType f_treat = 0.0;
  if ( treatment ) // Backward Euler
    f_treat = calc_treatment(tsol_loc+dtn);
  DType bc_treat = 0.0;
  if ( treatment_bc ) // Forward Euler Euler
    bc_treat = tdCoef->treatP1_bc * calc_treatment(tsol_loc);

  downscale_var( *data_P[fr], *data_P[fr], x_c, Nx, data_r[fr] );

  // Slopes
  (*slope_P[0])[0] = 0.0;
  for ( int i = 1; i < Nx-1; i++ ) {
    DType loc_r_l = x_c[i-1] * data_r[fr];
    DType loc_r_c = x_c[i] * data_r[fr];
    DType loc_r_r = x_c[i+1] * data_r[fr];
    (*slope_P[0])[i] = Alg_LimSlope_Minmod(
        loc_r_c * loc_r_c * (*data_P[fr])[i], loc_r_l * loc_r_l * (*data_P[fr])[i-1],
        loc_r_r * loc_r_r * (*data_P[fr])[i+1], dx );
  }
  // Right boundary
  {
    DType loc_r_l = x_c[Nx-2] * data_r[fr];
    DType loc_r_c = x_c[Nx-1] * data_r[fr];
    DType loc_r_r = (2.0*x_c[Nx-1]-x_c[Nx-2]) * data_r[fr];
    (*slope_P[0])[Nx-1] = Alg_LimSlope_Minmod(
        loc_r_c * loc_r_c * (*data_P[fr])[Nx-1], loc_r_l * loc_r_l * (*data_P[fr])[Nx-2],
        loc_r_r * loc_r_r * ((*data_P[fr])[Nx-1] - dx*data_r[fr]*bc_treat), dx );
  }

  // Fluxes
  (*flux_P[0])[0] = 0.0;
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_r_l = x_c[i-1] * data_r[fr];
    DType loc_r_r = x_c[i] * data_r[fr];
    DType loc_vel = - x_f[i] * dror;
    (*flux_P[0])[i] = loc_vel * ( 
        (loc_vel>=0.0) ? loc_r_l*loc_r_l*(*data_P[fr])[i-1] + 0.5*dx*(*slope_P[0])[i-1] :
                         loc_r_r*loc_r_r*(*data_P[fr])[i] - 0.5*dx*(*slope_P[0])[i] );
  }
  // Right boundary
  {
    DType loc_r_l = x_c[Nx-1] * data_r[fr];
    DType loc_r_r = (2.0*x_c[Nx-1]-x_c[Nx-2]) * data_r[fr];
    DType loc_vel = - x_f[Nx] * dror;
    (*flux_P[0])[Nx] = loc_vel * (
      (loc_vel>=0.0) ? loc_r_l*loc_r_l*(*data_P[fr])[Nx-1] + 0.5*dx*(*slope_P[0])[Nx-1] : 
                       loc_r_r*loc_r_r*((*data_P[fr])[Nx-1] - dx*data_r[fr]*bc_treat) );
  }

  // RHS
  VComp comp(1,0.0);
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_r_fr = x_c[i] * data_r[fr];
    DType loc_r_to = x_c[i] * data_r[to];
    comp = loc_r_fr * loc_r_fr * (*data_P[fr])[i];
    comp += dtn * coef_alpha * (*data_source[to])[i] / 
      (coef_beta/theta+(*data_source[to])[i]/loc_r_to/loc_r_to);
    comp -= ((*flux_P[0])[i+1] - (*flux_P[0])[i]) * (dtn/dx);
    comp += dtn * loc_r_to * loc_r_to * f_treat;
    comp /= loc_r_to * loc_r_to;
    rhs_p->set(i,comp);
  }
  if ( bc_treat != 0.0 ) { // Fix the non-hom Neumann condition
    comp = rhs_p->get(Nx-1);
    DType loc_r_to = x_c[Nx-1] * data_r[to];
    comp *= loc_r_to * loc_r_to;

    comp += dtn * coef_D * x_f[Nx] * x_f[Nx] * 
      ( - data_r[to] * bc_treat ) / dx;

    comp /= loc_r_to * loc_r_to;
    rhs_p->set(Nx-1,comp);
  }

  // Matrix
  Entry elem(1,0.0);

  // Left boundary
  {
    DType loc_r_0 = x_c[0] * data_r[to];

    DType coef_diff_1 = dtn * coef_D * x_f[1] * x_f[1] / dx / dx;

    elem = 1.0 + dtn * ( coef_gamma + dror ) + coef_diff_1 / loc_r_0 / loc_r_0;
    mat_p->set(0,0,elem);

    elem = - coef_diff_1 / loc_r_0 / loc_r_0;
    mat_p->set(0,1,elem);
  }
  // Interior cells
  for ( int i = 1; i < Nx-1; i++ ) {
    DType loc_r_c = x_c[i] * data_r[to];
    DType coef_diff_r = dtn * coef_D * x_f[i+1] * x_f[i+1] / dx / dx;
    DType coef_diff_l = dtn * coef_D * x_f[i] * x_f[i] / dx / dx;
    
    elem = - coef_diff_l / loc_r_c / loc_r_c;
    mat_p->set(i,i-1,elem);

    elem = 1.0 + dtn * ( coef_gamma + dror ) + 
      (coef_diff_l + coef_diff_r) / loc_r_c / loc_r_c;
    mat_p->set(i,i,elem);

    elem = - coef_diff_r / loc_r_c / loc_r_c;
    mat_p->set(i,i+1,elem);
  }
  // Right boundary
  {
    DType loc_r_c = x_c[Nx-1] * data_r[to];

    DType coef_diff_l = dtn * coef_D * x_f[Nx-1] * x_f[Nx-1] / dx / dx;

    elem = - coef_diff_l / loc_r_c / loc_r_c;
    mat_p->set(Nx-1,Nx-2,elem);

    elem = 1.0 + dtn * ( coef_gamma + dror ) + coef_diff_l / loc_r_c / loc_r_c;
    mat_p->set(Nx-1,Nx-1,elem);
  }

  // Initialize solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_P[fr])[i];
    sol_p->set(i,comp);
  }

  upscale_var( *data_P[fr], *data_P[fr], x_c, Nx, data_r[fr] );

  linSolve_p->solve();

  // Fetch solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = sol_p->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_P[to])[i] = vsol[0];
  }

  upscale_var( *data_P[to], *data_P[to], x_c, Nx, data_r[to] );

  return;
}


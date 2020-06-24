#include <cmath>
#include <stdexcept>
#include "TDSolveSacado.h"
#include "Alloc_Utils.h"
#include "TDInput.h"
#include "Alg_Solve.h"
#include "Alg_Limiters.h"

TDSolveSacado::TDSolveSacado(TDInput* tdInput_ptr) : 
  TDSolve(tdInput_ptr) {
  data_inc = NULL;
  if ( !(tdSetup->silence_output) ) {
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

TDSolveSacado::~TDSolveSacado(void) {
  erase_array( data_inc, 1 );
}

void TDSolveSacado::setup(void) {
  setupLinSolve( mat, prc, sol, rhs, res, linSolve, Nx );
  
  // Non test problems
  if ( tdSetup->prob_id >= 0 )
    setupLinSolve( mat_a, prc_a, sol_a, rhs_a, res_a, linSolve_a, Na );

  // OPN model
  if ( tdSetup->prob_id == 1 ) 
    setupLinSolve( mat_p, prc_p, sol_p, rhs_p, res_p, linSolve_p, Nx );

  return;
}

void TDSolveSacado::set_ic(void) {
  data_r[0] = tdCoef->R_0;
  bool invalid_ic = false;
  if ( tdSetup->prob_id >= 0 ) { // Tumor model
    switch ( tdSetup->ic ) {
      case 0: 
        setic_unif_cell(0.84,0.16-tdCoef->n_bc,tdCoef->n_bc);
        setic_expd_chem(1.0e+3,1.0,Na);
        break;
      case 1: 
        setic_unif_cell(5.0*(1.0-tdCoef->n_bc)/6.0,(1.0-tdCoef->n_bc)/6.0,tdCoef->n_bc);
        setic_expd_chem(1.0e+3,1.0,Na);
        break;
      default:
        invalid_ic = true;
        break;
    }
  }
  else { // Simpler model
    switch ( tdSetup->ic ) {
      case 0:
        setic_unif_cell(0.0,1.0-tdCoef->n_bc,tdCoef->n_bc);
        break;
      case 1:
        setic_unif_cell(1.0-tdCoef->n_bc,0.0,tdCoef->n_bc);
        break;
      case 2:
        setic_unif_cell(0.5-0.5*tdCoef->n_bc,0.5-0.5*tdCoef->n_bc,tdCoef->n_bc);
        break;
      default:
        invalid_ic = true;
        break;
    }
  }

  if ( invalid_ic ) {
    std::cout << "Error: No initial condition is defined ... Exit ...\n" << std::endl;
    exit(1);
  }

  upscale_var( *data_g[0], *data_g[0], x_c, Nx, data_r[0] );
  upscale_var( *data_h[0], *data_h[0], x_c, Nx, data_r[0] );
  upscale_var( *data_n[0], *data_n[0], x_c, Nx, data_r[0] );
  upscale_var( *data_a[0], *data_a[0], x_a_c, Na, data_r[0] );

  calc_total_inc_viol(0);

  if ( tdSetup->prob_id == 1 ) {
    (*data_p1[0]) = 0.0;
    (*data_p2[0]) = 0.0;
  }

  update_vel(0, 1, 0.0);

  return;
}

void TDSolveSacado::calc_dt(const int nt) {
  return;
}

void TDSolveSacado::output(const int nt) {
  return;
}

void TDSolveSacado::alloc(void) {
  Ndata = 3;
  
  alloc_array_full( data_g, Nx, Ndata );
  alloc_array_full( data_h, Nx, Ndata );
  alloc_array_full( data_n, Nx, Ndata ); // M
  alloc_array_full( data_v, Nx+1, Ndata );

  alloc_array_full( data_cell, Nx, 1 );
  alloc_array_full( data_node, Nx+1, 1 );
  
  data_r  = new DType[Ndata];

  alloc_array_full( data_a, Na, Ndata ); 

  if ( tdSetup->prob_id == 1 ) {
    alloc_array_full( data_p1, Nx, Ndata );
    alloc_array_full( data_p2, Nx, Ndata );
  }

  alloc_sacado_ind();

  alloc_sacado_dep();

  return;
}

bool TDSolveSacado::update_sol(const int fr, const int to, const DType shift_dt_fr, const DType shift_dt_to) {
  return true;
}

void TDSolveSacado::copy_data(const int fr, const int to) {
  return;
}

void TDSolveSacado::update_vel(const int fr, const int to, const DType shift_dt, const bool last_iter) {
  // Compute the velocity u using the gradient
  compute_vel_u( fr, tsol+shift_dt*dtn );
  return;
}

void TDSolveSacado::setic_unif_cell(DType val_g, DType val_h, DType val_n, int nt) {
  (*data_g[nt]) = val_g;
  (*data_h[nt]) = val_h;
  (*data_n[nt]) = val_n;
  return;
}

void TDSolveSacado::setic_expd_chem(DType multip, DType decay, int Na, int nt) {
  for ( int i = 0; i < Na; i++ ) {
    DType tvar = x_a_c[i] * data_r[nt];
    (*data_a[nt])[i] = multip * std::exp( - decay * tvar * tvar );
  }
  return;
}

void TDSolveSacado::alloc_sacado_ind(void) {
  return;
}

void TDSolveSacado::alloc_sacado_ind_array(std::string name, const int num) {
  std::map<std::string, std::vector<Sacado::Fad::DFad<DType>>>::iterator it = S_data_var.find( name );
  if ( it != S_data_var.end() ) {
    std::string error_message = "TDSolveSacado::alloc_sacado_ind_array - the independent array " + name + " is already allocated.\n";
    throw std::runtime_error( error_message );
  }
  // Get the number of existing independent variables
  int cur_num = 0;
  for ( it = S_data_var.begin(); it != S_data_var.end(); it++ ) 
    cur_num += it->second.size();
  S_data_offset.insert( std::pair<std::string, int>(name, cur_num) );
  return;
}

void TDSolveSacado::alloc_sacado_dep(void) {
  return;
}

void TDSolveSacado::calc_total_inc_viol(const int nt) {
  downscale_var( *data_g[nt], *data_g[nt], x_c, Nx, data_r[nt]);
  downscale_var( *data_h[nt], *data_h[nt], x_c, Nx, data_r[nt]);
  downscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt]);

  for ( int i = 0; i < Nx; i++ )
    (*data_inc[0])[i] = 1.0 - (*data_g[nt])[i] - (*data_h[nt])[i] - (*data_n[nt])[i];

  upscale_var( *data_g[nt], *data_g[nt], x_c, Nx, data_r[nt] );
  upscale_var( *data_h[nt], *data_h[nt], x_c, Nx, data_r[nt] );
  upscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );
  upscale_var( *data_inc[0], *data_inc[0], x_c, Nx, data_r[nt] );

  return;
}

void TDSolveSacado::compute_vel_u(const int nt, const DType tsol_loc) {
  if ( tdSetup->prob_id < 0 )
    compute_vel_u_pres( nt, tsol_loc );
  else
    compute_vel_u_chem( nt );

  return;
}

void TDSolveSacado::compute_vel_u_pres(const int nt, const DType tsol_loc) {

  switch ( tdSetup->prob_id ) {
    case -1: // Test 1, u is only a function of x
      break;
    case -2: // Test 2, space homogeneous solutions
      {
        //DType v_par = tdSetup->init_u_par;
        //(*data_u[0])[0] = 0.0;
        //for ( int i = 1; i <= Nx; i++ )
        //  (*data_u[0])[i] = - x_f[i] * data_r[nt] * v_par / ( 1.0 - std::exp(-3.0*v_par*tsol_loc) * (1.0-tdCoef->n_bc));
      }
      break;
    case -3: // Test 3, linear growth
      {
        //DType v_par = tdSetup->init_u_par;
        //(*data_u[0])[0] = 0.0;
        //for ( int i = 1; i <= Nx; i++ )
        //  (*data_u[0])[i] = -  v_par * std::sin( 0.5 * PI * x_f[i] ) / tdCoef->n_bc;
      }
      break;
    case -4: // Test 4, non-monotone radius change
      {
        //DType v_par = tdSetup->init_u_par;
        //(*data_u[0])[0] = 0.0;
        //for ( int i = 1; i <= Nx; i++ )
        //  (*data_u[0])[i] = v_par * std::sin( 4.0 * x_f[i] * data_r[nt] * ( 1.0 + tsol_loc ) );
      }
      break;
  }

  return;
}

void TDSolveSacado::compute_vel_u_chem(const int nt) {
  //*data_u[0] = 0.0;

  //DType alpha = tdCoef->alpha;
  //DType aor  = alpha/data_r[nt];

  //downscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );

  //for ( int i = 1; i <= Nx; i++ )
  //  (*data_u[0])[i] = aor * (((*data_a[nt])[i]-(*data_a[nt])[i-1]) / dx);

  //upscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );

  return;
}

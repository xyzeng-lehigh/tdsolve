#include <cmath>
#include <stdexcept>
#include "TDSolveSacado.h"
#include "Alloc_Utils.h"
#include "TDInput.h"
#include "Alg_Solve.h"
#include "Alg_Limiters.h"

TDSolveSacado::TDSolveSacado(TDInput* tdInput_ptr) : 
  TDSolve(tdInput_ptr) {
  data_inc   = NULL;
  S_data_num = 0;
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

  calc_total_inc_viol(0);

  if ( tdSetup->prob_id == 1 ) {
    (*data_p1[0]) = 0.0;
    (*data_p2[0]) = 0.0;
  }

  set_sacado_ind(0);

  update_vel(0, 1, 0.0, false);

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

  S_data_num += Nx;   // _g
  S_data_num += Nx;   // _h
  S_data_num += Nx;   // _n

  alloc_array_full( data_cell, Nx, 1 );
  alloc_array_full( data_node, Nx+1, 1 );
  
  data_r  = new DType[Ndata];
  S_data_num += 1;    // _r

  alloc_array_full( data_a, Na, Ndata ); 

  S_data_num += Na;   // _a

  if ( tdSetup->prob_id == 1 ) {
    alloc_array_full( data_p1, Nx, Ndata );
    alloc_array_full( data_p2, Nx, Ndata );
    S_data_num += Nx;  // _p1
    S_data_num += Nx;  // _p2
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
  compute_vel_u( tsol+shift_dt*dtn );
  // Compute the flux associated with U
  // XYZ
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
  alloc_sacado_ind_array("G", Nx);
  alloc_sacado_ind_array("H", Nx);
  alloc_sacado_ind_array("N", Nx);
  alloc_sacado_ind_array("R", 1);
  alloc_sacado_ind_array("A", Na);
  
  if ( tdSetup->prob_id == 1 ) {
    alloc_sacado_ind_array("P1", Nx);
    alloc_sacado_ind_array("P2", Nx);
  }
  return;
}

void TDSolveSacado::alloc_sacado_ind_array(std::string name, const int num) {
  std::map<std::string, SacadoIndVarArray>::iterator it = S_data_var.find( name );
  if ( it != S_data_var.end() ) {
    std::string error_message = "TDSolveSacado::alloc_sacado_ind_array - the independent array " + name + " is already allocated.\n";
    throw std::runtime_error( error_message );
  }
  // Get the number of existing independent variables
  int cur_num = 0;
  for ( it = S_data_var.begin(); it != S_data_var.end(); it++ ) 
    cur_num += it->second.size();
  if ( cur_num + num > S_data_num ) {
    std::string error_message = "TDSolveSacado::alloc_sacado_ind_array - requesting " + std::to_string(num) + " Sacado independent variables but there are only " + std::to_string(S_data_num-cur_num) + " opening left (" + std::to_string(cur_num) + "/" + std::to_string(S_data_num) + ").\n";
    throw std::runtime_error( error_message );
  }
  SacadoIndVarArray S_array;
  for ( int n = 0; n < num; n++ ) {
    Sacado::Fad::DFad<DType> value(S_data_num,cur_num+n,0.0);
    S_array.push_back(value);
  }
  S_data_var.insert( std::pair<std::string, SacadoIndVarArray>(name, S_array) );
  S_data_offset.insert( std::pair<std::string, int>(name, cur_num) );
  S_data_size.insert( std::pair<std::string, int>(name, num) );
  return;
}

void TDSolveSacado::alloc_sacado_dep(void) {
  std::vector<std::string>::iterator it_name;

  // Infiltration velocity
  std::vector<std::string> velnames = { "V", "U" };
  for ( it_name = velnames.begin(); it_name != velnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, Nx+1);
    S_data_vel.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }

  // Total flux and segregated fluxes for cell numbers
  std::vector<std::string> cellnames   = { "G", "H", "N" };
  std::vector<std::string> cellnames_u = { "N" };
  // Total flux
  for ( it_name = cellnames.begin(); it_name != cellnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, Nx+1);
    S_flux_cell.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  // Flux associated with v
  for ( it_name = cellnames.begin(); it_name != cellnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, Nx+1);
    S_flux_cell_v.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  // Flux associated with r
  for ( it_name = cellnames.begin(); it_name != cellnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, Nx+1);
    S_flux_cell_r.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  // Flux associated with u
  for ( it_name = cellnames_u.begin(); it_name != cellnames_u.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, Nx+1);
    S_flux_cell_u.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  // Limiter for PPM flux, two per interval
  for ( it_name = cellnames.begin(); it_name != cellnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, Nx);
    S_lim_ppm_cell[0].insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  for ( it_name = cellnames.begin(); it_name != cellnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, Nx);
    S_lim_ppm_cell[1].insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  // Synced PPM limiters
  alloc_sacado_dep_array(S_lim_ppm_sync[0],Nx);
  alloc_sacado_dep_array(S_lim_ppm_sync[1],Nx);

  // Chemoattractant/OPN
  std::vector<std::string> chemnames = { "A" };
  if ( tdSetup->prob_id == 1 ) {
    chemnames.push_back( "P1" );
    chemnames.push_back( "P2" );
  }
  // Fluxes for chems
  for ( it_name = chemnames.begin(); it_name != chemnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, S_data_size[*it_name]);
    S_flux_chem.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  return;
}

void TDSolveSacado::alloc_sacado_dep_array(SacadoDepVarArray& S_array, const int num) {
  S_array.clear();
  for ( int n = 0; n < num; n++ ) {
    SacadoDepVar value(0.0);
    S_array.push_back(value);
  }
  return;
}

void TDSolveSacado::set_sacado_ind(const int nt) {
  set_sacado_ind_array("G", *data_g[nt]);
  set_sacado_ind_array("H", *data_h[nt]);
  set_sacado_ind_array("N", *data_n[nt]);
  set_sacado_ind_array("A", *data_a[nt]);
  set_sacado_ind_array("P1", *data_p1[nt]);
  set_sacado_ind_array("P2", *data_p2[nt]);
  set_sacado_ind_array("R", data_r[nt]);
  return;
}

void TDSolveSacado::set_sacado_ind_array(std::string name, HVEC& data) {
  std::map<std::string, SacadoIndVarArray>::iterator it = S_data_var.find( name );
  if ( it == S_data_var.end() ) return;
  SacadoIndVarArray& s_var = it->second;
  int num = s_var.size();
  for ( int i = 0; i < num; i++ ) {
    s_var[i] = data[i];
  }
  return;
}

void TDSolveSacado::set_sacado_ind_array(std::string name, DType data) {
  std::map<std::string, SacadoIndVarArray>::iterator it = S_data_var.find( name );
  if ( it == S_data_var.end() ) return; 
  SacadoIndVarArray& s_var = it->second;
  s_var[0] = data;
  return;
}

void TDSolveSacado::calc_total_inc_viol(const int nt) {
  for ( int i = 0; i < Nx; i++ )
    (*data_inc[0])[i] = 1.0 - (*data_g[nt])[i] - (*data_h[nt])[i] - (*data_n[nt])[i];
  return;
}

void TDSolveSacado::compute_vel_u(const DType tsol_loc) {
  if ( tdSetup->prob_id < 0 )
    compute_vel_u_pres(tsol_loc );
  else
    compute_vel_u_chem();
  return;
}

void TDSolveSacado::compute_vel_u_pres(const DType tsol_loc) {
  SacadoDepVarArray& s_u = S_data_vel["U"];
  SacadoIndVarArray& s_r = S_data_var["R"];

  switch ( tdSetup->prob_id ) {
    case -1: // Test 1, u is only a function of x
      break;
    case -2: // Test 2, space homogeneous solutions
      {
        DType v_par = tdSetup->init_u_par;
        s_u[0] = 0.0;
        for ( int i = 1; i <= Nx; i++ )
          s_u[i] = - x_f[i] * s_r[0] * v_par / ( 1.0 - std::exp(-3.0*v_par*tsol_loc) * (1.0-tdCoef->n_bc));
      }
      break;
    case -3: // Test 3, linear growth
      {
        DType v_par = tdSetup->init_u_par;
        s_u[0] = 0.0;
        for ( int i = 1; i <= Nx; i++ )
          s_u[i] = -  v_par * std::sin( 0.5 * PI * x_f[i] ) / tdCoef->n_bc;
      }
      break;
    case -4: // Test 4, non-monotone radius change
      {
        DType v_par = tdSetup->init_u_par;
        s_u[0] = 0.0;
        for ( int i = 1; i <= Nx; i++ )
          s_u[i] = v_par * std::sin( 4.0 * x_f[i] * s_r[0] * ( 1.0 + tsol_loc ) );
      }
      break;
  }

  return;
}

void TDSolveSacado::compute_vel_u_chem(void) {
  SacadoDepVarArray& s_u = S_data_vel["U"];
  SacadoIndVarArray& s_r = S_data_var["R"];
  SacadoIndVarArray& s_a = S_data_var["A"];

  DType alpha = tdCoef->alpha;
  SacadoDepVar aor  = alpha/s_r[0];

  for ( int i = 1; i <= Nx; i++ )
    s_u[i] = aor * ( (s_a[i] - s_a[i-1]) / dx );

  return;
}

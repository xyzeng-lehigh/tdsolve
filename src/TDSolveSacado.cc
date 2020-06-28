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
  cellnames.clear();
  cellnames_u.clear();
  velnames.clear();
  chemnames.clear();
  varsizes.clear();

  cellnames.push_back( "G" ); 
  cellnames.push_back( "H" ); 
  cellnames.push_back( "N" ); 

  varsizes.insert(std::pair<std::string,int>("G",Nx));
  varsizes.insert(std::pair<std::string,int>("H",Nx));
  varsizes.insert(std::pair<std::string,int>("N",Nx));

  varcoord.insert(std::pair<std::string,DType*>("G",x_c));
  varcoord.insert(std::pair<std::string,DType*>("H",x_c));
  varcoord.insert(std::pair<std::string,DType*>("N",x_c));

  cellnames_u.push_back( "N" );

  velnames.push_back( "V" ); 
  velnames.push_back( "U" ); 

  varsizes.insert(std::pair<std::string,int>("V",Nx+1));
  varsizes.insert(std::pair<std::string,int>("U",Nx+1));

  varcoord.insert(std::pair<std::string,DType*>("V",x_f));
  varcoord.insert(std::pair<std::string,DType*>("U",x_f));

  if ( tdSetup->prob_id == 0 ) { // PDGF
    chemnames.push_back( "A" ); 
    varsizes.insert(std::pair<std::string,int>("A",Na));
    varcoord.insert(std::pair<std::string,DType*>("A",x_a_c));
  }
  if ( tdSetup->prob_id == 1 ) { // OPN
    chemnames.push_back( "A" );  
    chemnames.push_back( "P1" ); 
    chemnames.push_back( "P2" ); 

    varsizes.insert(std::pair<std::string,int>("A", Na));
    varsizes.insert(std::pair<std::string,int>("P1",Nx));
    varsizes.insert(std::pair<std::string,int>("P2",Nx));

    varcoord.insert(std::pair<std::string,DType*>("A", x_a_c));
    varcoord.insert(std::pair<std::string,DType*>("P1",x_c));
    varcoord.insert(std::pair<std::string,DType*>("P2",x_c));
  }

  if ( !(tdSetup->silence_output) ) {
    FILE* ot_file = fopen( tdInput->outputFile, "w" );
    fprintf( ot_file, "Time Time step size Radius Num. point /newline X_c ");
    for ( std::vector<std::string>::iterator it = cellnames.begin(); 
          it != cellnames.end(); it++ ) 
      fprintf( ot_file, "%s", *it + " " );
    for ( std::vector<std::string>::iterator it = chemnames.begin();
          it != chemnames.end(); it++ )
      fprintf( ot_file, "%s", *it + " " );
    fprintf( ot_file, "/newline X_f" );
    for ( std::vector<std::string>::iterator it = velnames.begin();
          it != velnames.end(); it++ )
      fprintf( ot_file, "%s", " " + *it );
    fprintf( ot_file, "\n" );
    fclose( ot_file );
  }
}

TDSolveSacado::~TDSolveSacado(void) {
  erase_array( data_inc, 1 );
  for ( std::map<std::string,HVEC**>::iterator it = H_data_res.begin();
        it != H_data_res.end(); it++ ) {
    erase_array( it->second, 1 );
  }
}

void TDSolveSacado::setup(void) {
  setupLinSolve( mat, prc, sol, rhs, res, linSolve, Nx );
  
  // PDGF problems
  if ( tdSetup->prob_id == 0 )
    setupLinSolve( mat_a, prc_a, sol_a, rhs_a, res_a, linSolve_a, Na );

  // OPN model
  if ( tdSetup->prob_id == 1 ) {
    setupLinSolve( mat_a, prc_a, sol_a, rhs_a, res_a, linSolve_a, Na );
    setupLinSolve( mat_p, prc_p, sol_p, rhs_p, res_p, linSolve_p, Nx );
  }

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

  // This is only to compute the initial time step size
  // All Sacado arrays, velocity, etc., will be recomputed 
  // at the beginning of update_sol
  {
    set_sacado_ind(0);
    set_sacado_r(0);
    update_vel(0, 1, 0.0, false);
    set_from_sacado_v(0);
  }

  return;
}

void TDSolveSacado::calc_dt(const int nt) {
  if ( iter != 0 )
    tsol += dtn;

  if ( data_r[nt] <= tdSetup->Rdead || data_r[nt] >= tdSetup->Rterm )
    reach_term_now = true;

  if ( reach_term_next )
    reach_term_now = true;

  DType dtn_trial = DTYPE_MAX;
  DType dtn_cfl1  = DTYPE_MAX; // dtn corresponding to CFL 1

  SacadoDepVarArray& s_v = S_data_vel["V"];
  SacadoDepVarArray& s_u = S_data_vel["U"];

  for ( int i = 0; i < Nx+1; i++ ) {
    DType v_max = std::fabs( x_f[i] * S_data_rp.val() / data_r[nt] );
    for ( SacadoDepVarMap::iterator it = S_data_vel.begin();
          it != S_data_vel.end(); it++ ) {
      v_max += std::fabs( it->second[i].val() / data_r[nt] );
    }
    dtn_trial = min( dtn_trial, tdSetup->cfl * dx / v_max );
    dtn_cfl1  = min( dtn_cfl1, dx / v_max );
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

  cur_cfl = dtn/dtn_cfl1;

  return;
}

void TDSolveSacado::output(const int nt) {
  if ( tdSetup->silence_output ) return;

  FILE* ot_file = fopen( tdInput->radiusFile, "a" );
  fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\t%10.20f\n", tsol, data_r[nt], 
      S_data_rp.val(), S_data_vel["U"][varsizes["U"]-1].val() );
  fclose( ot_file );

  if ( !reach_term_now ) {
    if ( use_tfreq ) {
      if ( !reach_output_now ) return;
    }
    else {
      if ( iter % tdSetup->freq != 0 ) return;
    }
  }

  ot_file = fopen( tdInput->outputFile, "a" );
  fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\t%d\t%d\n", tsol, dtn, data_r[nt], Nx, Nx+1 );

  for ( int i = 0; i < Nx; i++ ) {
    fprintf( ot_file, "%10.20f\t", x_c[i] * data_r[nt] );
    for ( std::vector<std::string>::iterator it = cellnames.begin();
          it != cellnames.end(); it++ )
      fprintf( ot_file, "%10.20f\t", (*H_data_var[*it][nt])[i] );
    for ( std::vector<std::string>::iterator it = chemnames.begin();
          it != chemnames.end(); it++ )
      fprintf( ot_file, "%10.20f\t", (*H_data_var[*it][nt])[i] );
  }
  fprintf( ot_file, "\n" );

  for ( int i = 0; i < Nx+1; i++ ) {
    fprintf( ot_file, "%10.20f\t", x_f[i] * data_r[nt] );
    for ( std::vector<std::string>::iterator it = velnames.begin();
          it != velnames.end(); it++ )
      fprintf( ot_file, "%10.20f\t", S_data_vel[*it][i].val() );
  }
  fprintf( ot_file, "\n" );

  fclose( ot_file );

  return;
}

void TDSolveSacado::alloc(void) {
  Ndata = 3;
  
  alloc_array_full( data_g, Nx, Ndata );
  alloc_array_full( data_h, Nx, Ndata );
  alloc_array_full( data_n, Nx, Ndata ); // M
  alloc_array_full( data_v, Nx+1, Ndata );

  H_data_var.insert(std::pair<std::string,HVEC**>("G",data_g));
  H_data_var.insert(std::pair<std::string,HVEC**>("H",data_h));
  H_data_var.insert(std::pair<std::string,HVEC**>("N",data_g));

  S_data_num += Nx;   // _g
  S_data_num += Nx;   // _h
  S_data_num += Nx;   // _n

  alloc_array_full( data_cell, Nx, 1 );
  alloc_array_full( data_node, Nx+1, 1 );
  
  data_r  = new DType[Ndata];

  alloc_array_full( data_a, Na, Ndata ); 

  H_data_var.insert(std::pair<std::string,HVEC**>("A",data_a));

  S_data_num += Na;   // _a

  if ( tdSetup->prob_id == 1 ) {
    alloc_array_full( data_p1, Nx, Ndata );
    alloc_array_full( data_p2, Nx, Ndata );
    H_data_var.insert(std::pair<std::string,HVEC**>("P1",data_p1));
    H_data_var.insert(std::pair<std::string,HVEC**>("P2",data_p2));
    S_data_num += Nx;  // _p1
    S_data_num += Nx;  // _p2
  }

  alloc_hvec_res();

  alloc_sacado_ind();

  alloc_sacado_dep();

  alloc_functions();

  return;
}

bool TDSolveSacado::update_sol(const int fr, const int to, const DType shift_dt_fr, const DType shift_dt_to) {
  // fr - solution at t_n
  // to - solution at t_{n+1}/Newton iterates
  copy_data( fr, to );

  // Setup updates schemes
  int  numIter  = 1;
  bool assembleJacobian = false;
  switch ( tdSetup->assemb ) {
    case 0: // explicit (no Jacobian)
      numIter = 1;
      assembleJacobian = false;
      break;
    case 1: // loosely implicit
      numIter = 1;
      assembleJacobian = true;
      break;
    case 2: // fully impliict
      numIter = tdSetup->ntIter;
      assembleJacobian = true;
      break;
    case 3: // PMC - just add for fun, not really using
      numIter = tdSetup->ntIter;
      assembleJacobian = false;
      break;
    default:
      break;
  }

  for ( int iter = 0; iter < numIter; iter++ ) {
    // Step 0. Setup the Sacado variables
    reset_sacado_ind_array();
    reset_sacado_dep_array();
    // Set the independent variables iterates to [to] 
    set_sacado_ind( to );

    // Step 1. Compute S_data_r and S_data_rp.
    //   assemb == 0   : S_data_r = data_r[fr], S_data_rp ~ data_v[fr]
    //   assemb == 1,2 : S_data_r = data_r[t0], S_data_rp ~ data_v[to]
    // Note that because [fr] is copied to [to] in iter 0,
    // data_v[to] is always used in this function to compute S_data_rp
    compute_rad( fr, to );

    // Step 2. Compute velocity S_data_vel["V"], which requires S_data_r
    //         but not S_data_rp
    // This will also compute the following in sequel:
    //   S_data_vel["U"]
    //   S_flux_cell_v["U"]["N"]
    //   S_srce_cell["X"] (but excluding the geometric term -eta^2R'RX
    //   S_data_vel["V"]
    update_vel( fr, to, shift_dt_to, false );  

    // Step 3. Compute the source terms for chemicals.
    //         Then complete the source terms for both cells and chems
    //         with geometric terms.
    compute_source_chem_prob( tsol+shift_dt_to*dtn );
    update_source_geom();

  }

  return true;
}

void TDSolveSacado::copy_data(const int fr, const int to) {
  *data_g[to]  = *data_g[fr];
  *data_h[to]  = *data_h[fr];
  *data_n[to]  = *data_n[fr];
  *data_v[to]  = *data_v[fr];
  
  // PDGF
  if ( tdSetup->prob_id == 0 )
    *data_a[to]  = *data_a[fr];
  
  // OPN
  if ( tdSetup->prob_id == 1 ) {
    *data_a[to]  = *data_a[fr];
    *data_p1[to] = *data_p1[fr];
    *data_p2[to] = *data_p2[fr];
  }

  data_r[to]  =  data_r[fr];

  return;
}

void TDSolveSacado::update_vel(const int fr, const int to, const DType shift_dt, const bool last_iter) {
  // Compute the velocity u using the gradient
  compute_vel_u( tsol+shift_dt*dtn );
  // Compute the flux associated with U
  compute_cells_varflux( "U", "N", tsol+shift_dt*dtn, true );
  // Compute the right-hand sides (source terms)
  compute_source_cell_prob( tsol+shift_dt*dtn );
  // Compute the velocity (integral form)
  compute_vel_v();
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

void TDSolveSacado::alloc_hvec_res(void) {
  H_data_res.clear();
  for ( std::map<std::string,HVEC**>::iterator it = H_data_var.begin();
        it != H_data_var.end(); it++ ) {
    HVEC** loc_res;
    alloc_array_full( loc_res, varsizes[it->first], 1 );
    H_data_res.insert(std::pair<std::string,HVEC**>(it->first,loc_res));
  }
  return;
}

void TDSolveSacado::alloc_sacado_ind(void) {
  for ( std::vector<std::string>::iterator it = cellnames.begin();
        it != cellnames.end(); it++ )
    alloc_sacado_ind_array(*it,varsizes[*it]);
  for ( std::vector<std::string>::iterator it = chemnames.begin();
        it != chemnames.end(); it++ )
    alloc_sacado_ind_array(*it,varsizes[*it]);
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
  return;
}

void TDSolveSacado::alloc_sacado_dep(void) {
  std::vector<std::string>::iterator it_name;

  // Variable velocities
  for ( it_name = velnames.begin(); it_name != velnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, varsizes[*it_name]);
    S_data_vel.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }

  // Total flux and segregated fluxes for cell numbers
  // Total flux
  for ( it_name = cellnames.begin(); it_name != cellnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, Nx+1);
    S_flux_cell.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  // Flux associated with r
  for ( it_name = cellnames.begin(); it_name != cellnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, varsizes[*it_name]+1);
    S_flux_cell_r.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  // Flux associated with v
  SacadoDepVarMap S_array_v; S_array_v.clear();
  for ( it_name = cellnames.begin(); it_name != cellnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, varsizes[*it_name]+1);
    S_array_v.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  S_flux_cell_v.insert( std::pair<std::string, SacadoDepVarMap>("V", S_array_v) );
  // Flux associated with u
  SacadoDepVarMap S_array_u; S_array_u.clear();
  for ( it_name = cellnames_u.begin(); it_name != cellnames_u.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, varsizes[*it_name]+1);
    S_array_u.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  S_flux_cell_v.insert( std::pair<std::string, SacadoDepVarMap>("U", S_array_u) );
  // Limiter for PPM flux, two per interval
  for ( it_name = cellnames.begin(); it_name != cellnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, varsizes[*it_name]);
    S_lim_ppm_cell[0].insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  for ( it_name = cellnames.begin(); it_name != cellnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, varsizes[*it_name]);
    S_lim_ppm_cell[1].insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  // Synced PPM limiters
  alloc_sacado_dep_array(S_lim_ppm_sync[0],Nx);
  alloc_sacado_dep_array(S_lim_ppm_sync[1],Nx);
  // Source terms
  for ( it_name = cellnames.begin(); it_name != cellnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, varsizes[*it_name]);
    S_srce_cell.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }

  // Chemoattractant/OPN
  // Source terms
  for ( it_name = chemnames.begin(); it_name != chemnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, varsizes[*it_name]);
    S_srce_chem.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  // Fluxes for chems
  for ( it_name = chemnames.begin(); it_name != chemnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, varsizes[*it_name]);
    S_flux_chem.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  return;
}

void TDSolveSacado::alloc_functions(void) {
  F_calcSrceMap.clear();
  if ( tdSetup->prob_id == 0 ) { // PDGF
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("G", &TDSolveSacado::calc_srce_g_pdgf) );
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("H", &TDSolveSacado::calc_srce_h_pdgf) );
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("N", &TDSolveSacado::calc_srce_n_pdgf) );

    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("A", &TDSolveSacado::calc_srce_a_pdgf) );
  }
  else if ( tdSetup->prob_id == 1 ) { // OPN
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("G", &TDSolveSacado::calc_srce_g_opn) );
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("H", &TDSolveSacado::calc_srce_h_opn) );
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("N", &TDSolveSacado::calc_srce_n_opn) );

    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("A", &TDSolveSacado::calc_srce_a_pdgf) );
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("P1", &TDSolveSacado::calc_srce_p1_opn) );
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("P2", &TDSolveSacado::calc_srce_p2_opn) );
  }
  else {
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("G", &TDSolveSacado::calc_srce_g_base) );
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("H", &TDSolveSacado::calc_srce_h_base) );
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("N", &TDSolveSacado::calc_srce_n_base) );
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

void TDSolveSacado::set_sacado_r(const int nt) {
  S_data_r = data_r[nt];
  return;
}

void TDSolveSacado::set_sacado_v(const int nt) {
  SacadoDepVarArray& s_v = S_data_vel["V"];
  int num = s_v.size();
  for ( int i = 0; i < num; i++ ) {
    s_v[i] = (*data_v[nt])[i];
  }
  return;
}

void TDSolveSacado::set_from_sacado_ind(const int nt) {
  set_from_sacado_ind_array("G", *data_g[nt]);
  set_from_sacado_ind_array("H", *data_h[nt]);
  set_from_sacado_ind_array("N", *data_n[nt]);
  set_from_sacado_ind_array("A", *data_a[nt]);
  set_from_sacado_ind_array("P1", *data_p1[nt]);
  set_from_sacado_ind_array("P2", *data_p2[nt]);
  return;
}

void TDSolveSacado::set_from_sacado_ind_array(std::string name, HVEC& data) {
  std::map<std::string, SacadoIndVarArray>::iterator it = S_data_var.find( name );
  if ( it == S_data_var.end() ) return;
  SacadoIndVarArray& s_var = it->second;
  int num = s_var.size();
  for ( int i = 0; i < num; i++ ) {
    data[i] = s_var[i].val();
  }
  return;
}

void TDSolveSacado::set_from_sacado_r(const int nt) {
  data_r[nt] = S_data_r.val();
  return;
}

void TDSolveSacado::set_from_sacado_v(const int nt) {
  SacadoDepVarArray& s_v = S_data_vel["V"];
  int num = s_v.size();
  for ( int i = 0; i < num; i++ ) {
    (*data_v[nt])[i] = s_v[i].val();
  }
  return;
}

void TDSolveSacado::reset_sacado_ind_array(void) {
  for ( std::map<std::string, SacadoIndVarArray>::iterator it = S_data_var.begin();
        it != S_data_var.end(); it++ ) {
    reset_sacado_ind_array( it->second );
  }
  return;
}

void TDSolveSacado::reset_sacado_ind_array(SacadoIndVarArray& s_var) {
  for ( int i = 0; i < s_var.size(); i++ ) {
    s_var[i] = 0.0;
  }
  return;
}

void TDSolveSacado::reset_sacado_dep_array(void) {
  S_data_r  = 0.0;
  S_data_rp = 0.0;
  reset_sacado_dep_array( S_data_vel );
  reset_sacado_dep_array( S_flux_cell );
  reset_sacado_dep_array( S_flux_cell_r );
  reset_sacado_dep_array( S_flux_cell_v );
  reset_sacado_dep_array( S_srce_cell );
  for ( int i = 0; i < 2; i++ ) {
    reset_sacado_dep_array( S_lim_ppm_cell[i] );
    reset_sacado_dep_array( S_lim_ppm_sync[i] );
  }
  reset_sacado_dep_array( S_flux_chem );
  return;
}

void TDSolveSacado::reset_sacado_dep_array(SacadoDepVarMapMap& s_maps) {
  for ( SacadoDepVarMapMap::iterator it = s_maps.begin();
        it != s_maps.end(); it++ ) {
    reset_sacado_dep_array( it->second );
  }
  return;
}

void TDSolveSacado::reset_sacado_dep_array(SacadoDepVarMap& s_map) {
  for ( SacadoDepVarMap::iterator it = s_map.begin();
        it != s_map.end(); it++ ) {
    reset_sacado_dep_array( it->second );
  }
  return;
}

void TDSolveSacado::reset_sacado_dep_array(SacadoDepVarArray& s_var) {
  for ( int i = 0; i < s_var.size(); i++ ) {
    s_var[i] = 0.0;
  }
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

void TDSolveSacado::compute_vel_v(void) {
  SacadoDepVarArray& s_v = S_data_vel["V"];
  SacadoDepVarArray& s_flux_u = S_flux_cell_v["U"]["N"];
  
  s_v[0] = 0.0;

  for ( int i = 0; i < Nx; i++ ) {
    s_v[i+1] = s_v[i];
    for ( SacadoDepVarMap::iterator it = S_srce_cell.begin();
          it != S_srce_cell.end(); it++ ) {
      s_v[i+1] += (it->second)[i];
    }
  }

  for ( int i = 1; i < Nx+1; i++ ) {
    s_v[i] = (s_v[i] - s_flux_u[i])/x_f[i]/x_f[i]/S_data_r;
  }

  return;
}

void TDSolveSacado::compute_cells_varflux(std::string velname, std::string cellname, const DType tsol_loc, const bool rbc) {
  SacadoDepVarMap&   s_fluxes = S_flux_cell_v[velname];
  SacadoDepVarArray& s_flux   = s_fluxes[cellname];
  SacadoDepVarArray& s_vel    = S_data_vel[velname];
  SacadoIndVarArray& s_var    = S_data_var[cellname];
  DType* rdata = NULL;
  if (rbc) {
    rdata = new DType(calc_n_bc(tsol_loc));
  }

  switch ( tdSetup->scell ) {
    case 0: // First-order upwind
      compute_cells_varflux_upw( s_vel, s_var, s_flux, rdata );
      break;
    default:
      throw std::invalid_argument( "TDSolveSacado::compute_cells_fluxes_u -- infiltration flux method tag " + std::to_string(tdSetup->scell) + " not defined.\n" );
      break;
  }

  if (rdata!=NULL)
    delete rdata;

  return;
}

void TDSolveSacado::compute_cells_varflux_upw(SacadoDepVarArray& s_vel, SacadoIndVarArray& s_var, SacadoDepVarArray& s_flux, const DType* rdata) {
  reset_sacado_dep_array(s_flux);
  SacadoIndVarArray& s_r = S_data_var["R"];
  // Left boundary
  s_flux[0] = 0.0;
  // Interior points
  for ( int i = 1; i < Nx; i++ ) {
    SacadoDepVar s_loc_u  = x_f[i] * x_f[i] * s_r[0] * s_vel[i];
    if ( s_loc_u >= 0.0 )
      s_flux[i] = s_loc_u * s_var[i-1];
    else
      s_flux[i] = s_loc_u * s_var[i];
  }
  // Right boundary
  if ( rdata == NULL )
    s_flux[Nx] = 0.0;
  else {
    SacadoDepVar s_loc_u = x_f[Nx] * x_f[Nx] * s_r[0] * s_vel[Nx];
    if ( s_loc_u >= 0.0 )
      s_flux[Nx] = s_loc_u * s_var[Nx-1];
    else
      s_flux[Nx] = s_loc_u * (*rdata);
  }
  return;
}

void TDSolveSacado::compute_rad(const int fr, const int to) {
  if ( tdSetup->dr_opt == 2 ) // For legacy issue, there is no opt 1
    S_data_rp = (*data_v[to])[Nx]/(1.0-0.25*dx*dx);
  else
    S_data_rp = (*data_v[to])[Nx];

  if ( tdSetup->r_opt == 1 ) {
    if ( tdSetup->assemb == 0 ) { // Explicit
      data_r[to] = std::sqrt( data_r[fr]*data_r[fr] + 
          2.0*dtn*data_r[fr]*S_data_rp.val() );
    }
    else {
      data_r[to] = std::sqrt( data_r[fr]*data_r[fr] + dtn*dtn*
          S_data_rp.val()*S_data_rp.val()) + dtn*S_data_rp.val();
    }
  }
  else {
    data_r[to] = data_r[fr] + dtn * S_data_rp.val();
  }

  if ( tdSetup->assemb == 1 || tdSetup->assemb == 2 )
    S_data_r = data_r[to];
  else
    S_data_r = data_r[fr];

  return;
}

void TDSolveSacado::update_source_geom(void) {
  for ( std::vector<std::string>::iterator it = cellnames.begin();
        it != cellnames.end(); it++ ) {
    SacadoIndVarArray& s_var  = S_data_var[*it];
    SacadoDepVarArray& s_srce = S_srce_cell[*it];
    calc_srce_geom( s_var, s_srce, varcoord[*it] );
  }
  for ( std::vector<std::string>::iterator it = chemnames.begin();
        it != chemnames.end(); it++ ) {
    SacadoIndVarArray& s_var  = S_data_var[*it];
    SacadoDepVarArray& s_srce = S_srce_chem[*it];
    calc_srce_geom( s_var, s_srce, varcoord[*it] );
  }
  return;
}

void TDSolveSacado::calc_srce_geom(SacadoIndVarArray& s_var, SacadoIndVarArray& s_srce, DType* coord) {
  int num = s_var.size(); // This could be either Nx or Na
  for ( int i = 0; i < num; i++ ) {
    s_srce[i] -= coord[i] * coord[i] * S_data_rp * S_data_r * s_var[i]; 
  }
  return;
}

void TDSolveSacado::compute_source_cell_prob(const DType tsol_loc) {
  std::vector<std::string>::iterator it_name;
  for ( it_name = cellnames.begin(); it_name != cellnames.end(); it_name++ ) {
    SacadoDepVarArray& s_srce = S_srce_cell[*it_name];
    (this->*F_calcSrceMap[*it_name])(s_srce, tsol_loc);
  }
  return;
}

void TDSolveSacado::calc_srce_g_base(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  reset_sacado_dep_array(s_srce);
  return;
}

void TDSolveSacado::calc_srce_h_base(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  reset_sacado_dep_array(s_srce);
  return;
}

void TDSolveSacado::calc_srce_n_base(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  reset_sacado_dep_array(s_srce);
  return;
}

void TDSolveSacado::calc_srce_g_pdgf(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  calc_srce_g_base( s_srce, tsol_loc );

  DType lambda = tdCoef->lambda;
  DType mu     = tdCoef->mu;
  SacadoIndVarArray& s_g = S_data_var["G"];

  for ( int i = 0; i < Nx; i++ ) {
    s_srce[i] += x_c[i] * x_c[i] * S_data_r * S_data_r * ((lambda-mu) * s_g[i]);
  }

  return;
}

void TDSolveSacado::calc_srce_h_pdgf(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  calc_srce_h_base( s_srce, tsol_loc );

  DType mu      = tdCoef->mu;
  DType delta_h = tdCoef->delta_h;
  SacadoIndVarArray& s_g = S_data_var["G"];
  SacadoIndVarArray& s_h = S_data_var["H"];

  for ( int i = 0; i < Nx; i++ ) {
    s_srce[i] += (x_c[i] * x_c[i] * S_data_r * S_data_r) * ( mu * s_g[i] - delta_h * s_h[i] );
  }

  return;
}

void TDSolveSacado::calc_srce_n_pdgf(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  calc_srce_n_base( s_srce, tsol_loc );

  DType delta_n = tdCoef->delta_h;
  SacadoIndVarArray& s_n = S_data_var["N"];

  for ( int i = 0; i < Nx; i++ ) {
    s_srce[i] += (x_c[i] * x_c[i] * S_data_r * S_data_r) * ( - delta_n * s_n[i] );
  }

  return;
}

void TDSolveSacado::calc_srce_g_opn(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  calc_srce_g_base( s_srce, tsol_loc );
  return;
}

void TDSolveSacado::calc_srce_h_opn(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  calc_srce_h_base( s_srce, tsol_loc );
  return;
}

void TDSolveSacado::calc_srce_n_opn(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  calc_srce_n_base( s_srce, tsol_loc );
  return;
}

void TDSolveSacado::compute_source_chem_prob(const DType tsol_loc) {
  std::vector<std::string>::iterator it_name;
  for ( it_name = chemnames.begin(); it_name != chemnames.end(); it_name++ ) {
    SacadoDepVarArray& s_srce = S_srce_chem[*it_name];
    (this->*F_calcSrceMap[*it_name])(s_srce, tsol_loc);
  }
  return;
}

void TDSolveSacado::calc_srce_a_pdgf(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  reset_sacado_dep_array(s_srce);

  // To implement

  return;
}

void TDSolveSacado::calc_srce_p1_opn(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  reset_sacado_dep_array(s_srce);

  // To implement

  return;
}

void TDSolveSacado::calc_srce_p2_opn(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  reset_sacado_dep_array(s_srce);

  // To implement

  return;
}

DType TDSolveSacado::calc_n_bc( const DType tsol_loc ) {
  DType n_bc = tdCoef->n_bc;
  if ( tdSetup->prob_id == -2 ) {
    DType v_par = tdSetup->init_u_par;
    n_bc = 1.0 - std::exp(-3.0*v_par*tsol_loc) * ( 1.0-n_bc );
  }
  return n_bc;
}

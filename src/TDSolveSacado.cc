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
  varcoord.clear();
  varflxcoord.clear();
  chempars.clear();
  varlspacks.clear();
  dimlspacks.clear();
  H_cell_offset.clear();

  // Common for all existing prob_id:
  //  -1 : model problem
  //   0 : pdgf
  //   1 : opn
  cellnames.push_back( "G" ); 
  cellnames.push_back( "N" ); 
  varsizes.insert(std::pair<std::string,int>("G",Nx));
  varsizes.insert(std::pair<std::string,int>("N",Nx));
  varcoord.insert(std::pair<std::string,DType*>("G",x_c));
  varcoord.insert(std::pair<std::string,DType*>("N",x_c));
  varflxcoord.insert(std::pair<std::string,DType*>("G",x_f));
  varflxcoord.insert(std::pair<std::string,DType*>("N",x_f));

  cellnames_u.push_back( "N" );

  velnames.push_back( "V" ); 
  velnames.push_back( "U" ); 
  varsizes.insert(std::pair<std::string,int>("V",Nx+1));
  varsizes.insert(std::pair<std::string,int>("U",Nx+1));
  varcoord.insert(std::pair<std::string,DType*>("V",x_f));
  varcoord.insert(std::pair<std::string,DType*>("U",x_f));

  if ( tdSetup->prob_id == 0 || tdSetup->prob_id == 1 ) { // PDGF or OPN
    cellnames.push_back( "H" ); 
    varsizes.insert(std::pair<std::string,int>("H",Nx));
    varcoord.insert(std::pair<std::string,DType*>("H",x_c));
    varflxcoord.insert(std::pair<std::string,DType*>("H",x_f));

    chemnames.push_back( "A" ); 
    varsizes.insert(std::pair<std::string,int>("A",Na));
    varcoord.insert(std::pair<std::string,DType*>("A",x_a_c));
    varflxcoord.insert(std::pair<std::string,DType*>("A",x_a_f));
  }

  if ( tdSetup->prob_id == 1 ) { // OPN
    chemnames.push_back( "P1" ); 
    chemnames.push_back( "P2" ); 
    varsizes.insert(std::pair<std::string,int>("P1",Nx));
    varsizes.insert(std::pair<std::string,int>("P2",Nx));
    varcoord.insert(std::pair<std::string,DType*>("P1",x_c));
    varcoord.insert(std::pair<std::string,DType*>("P2",x_c));
    varflxcoord.insert(std::pair<std::string,DType*>("P1",x_f));
    varflxcoord.insert(std::pair<std::string,DType*>("P2",x_f));
  }

  if ( !(tdSetup->silence_output) ) {
    FILE* ot_file = fopen( tdInput->outputFile, "w" );
    fprintf( ot_file, "Time Time step size Radius Num. point /newline X_c ");
    for ( std::vector<std::string>::iterator it = cellnames.begin(); 
          it != cellnames.end(); it++ ) 
      fprintf( ot_file, "%s ", it->c_str() );
    fprintf( ot_file, "Inc-Viol "); // Incompressibility violation
    for ( std::vector<std::string>::iterator it = chemnames.begin();
          it != chemnames.end(); it++ )
      fprintf( ot_file, "%s ", it->c_str() );
    fprintf( ot_file, "/newline X_f" );
    for ( std::vector<std::string>::iterator it = velnames.begin();
          it != velnames.end(); it++ )
      fprintf( ot_file, " %s", it->c_str() );
    fprintf( ot_file, "\n" );
    fclose( ot_file );
  }
}

TDSolveSacado::~TDSolveSacado(void) {
  erase_array( data_inc, 1 );
  erase_array( data_u, Ndata );
  for ( std::map<std::string,HVEC**>::iterator it = H_data_res.begin();
        it != H_data_res.end(); it++ ) {
    erase_array( it->second, 1 );
  }
  for ( std::map<int,LSPack>::iterator it = dimlspacks.begin();
        it != dimlspacks.end(); it++ ) {
    delete it->second.mat;
    delete it->second.prc;
    delete it->second.rhs;
    delete it->second.res;
    delete it->second.sol;
    delete it->second.linSolve;
  }
}

void TDSolveSacado::setup(void) {
  std::vector<std::string>::iterator it;

  // Setup linear solver for cell numbers if an implicit method is used
  if ( tdSetup->assemb == 1 ) {
    //for ( it = cellnames.begin(); it != cellnames.end(); it++ ) {
    //  setup_new_lin_solve( varsizes[*it] );
    //  std::map<int,LSPack>::iterator it_ls = dimlspacks.find( varsizes[*it] );
    //  varlspacks.insert( std::pair<std::string,LSPack>( *it, it_ls->second ) );
    //}
    int size = 0;
    for ( it = cellnames.begin(); it != cellnames.end(); it++ ) {
      H_cell_offset.insert( std::pair<std::string,int>( *it, size ) );
      size += varsizes[*it];
    }
    setup_new_lin_solve( size );
    std::map<int,LSPack>::iterator it_ls = dimlspacks.find( size );
    cell_lspack = it_ls->second;
  }
  
  // Setup linear solver for chemicals
  for ( it = chemnames.begin(); it != chemnames.end(); it++ ) {
    setup_new_lin_solve( varsizes[*it] );
    std::map<int,LSPack>::iterator it_ls = dimlspacks.find( varsizes[*it] );
    varlspacks.insert( std::pair<std::string,LSPack>( *it, it_ls->second ) );
  }

  // OPN or PDGF problems
  if ( tdSetup->prob_id == 0 || tdSetup->prob_id == 1 ) {
    setup_chem_par( "A", tdCoef->gamma_a, tdCoef->D, "G", tdCoef->alpha_a, tdCoef->beta_a/tdCoef->theta, true );
  }

  // OPN model
  if ( tdSetup->prob_id == 1 ) {
    setup_chem_par( "P1", tdCoef->gamma_1, tdCoef->D_1, "N", tdCoef->alpha_1, tdCoef->beta_1/tdCoef->theta, false );
    setup_chem_par( "P2", tdCoef->gamma_2, tdCoef->D_2, "G", tdCoef->alpha_2, tdCoef->beta_2/tdCoef->theta, false );
  }

  return;
}

void TDSolveSacado::setup_new_lin_solve(const int num) {
  std::map<int,LSPack>::iterator it = dimlspacks.find( num );
  if ( it != dimlspacks.end() ) // A linear solver of size num already exist
    return;
  LSPack lsPack;
  setupLinSolve( lsPack.mat, lsPack.prc, lsPack.sol, lsPack.rhs, lsPack.res, lsPack.linSolve, num );
  dimlspacks.insert( std::pair<int,LSPack>(num,lsPack) );
  return;
}

void TDSolveSacado::set_ic(void) {
  data_r[0] = tdCoef->R_0;
  bool invalid_ic = false;
  if ( tdSetup->prob_id >= 0 ) { // Tumor model
    switch ( tdSetup->ic ) {
      case 0: 
        setic_unif_cell("G",0.84);
        setic_unif_cell("H",0.16-tdCoef->n_bc);
        setic_unif_cell("N",tdCoef->n_bc);
        setic_expd_chem(1.0e+3,1.0);
        break;
      case 1: 
        setic_unif_cell("G",5.0*(1.0-tdCoef->n_bc)/6.0);
        setic_unif_cell("H",(1.0-tdCoef->n_bc)/6.0);
        setic_unif_cell("N",tdCoef->n_bc);
        setic_expd_chem(1.0e+3,1.0);
        break;
      case 2:
        setic_unif_cell("G",1.0-tdCoef->n_bc);
        setic_unif_cell("H",0.0);
        setic_unif_cell("N",tdCoef->n_bc);
        setic_polyexp_chem(5.0/3.0,1.0/6.0);
        break;
      default:
        invalid_ic = true;
        break;
    }
  }
  else if ( tdSetup->prob_id == -1 ) { // Model problem, only G and N
    switch ( tdSetup->ic ) {
      case 0:
        setic_unif_cell("G",1.0-tdCoef->n_bc);
        setic_unif_cell("N",tdCoef->n_bc);
        break;
      default:
        invalid_ic = true;
        break;
    }
  }
  else { // Other models
    invalid_ic = true;
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
    // Initialize S_data_var with H_data_var["X"][0]
    set_sacado_ind(0);
    // Set S_data_r to data_r[0]
    set_sacado_r(0);
    // First two arguments are not used at all.
    // Compute S_data_vel["V"] from S_data_var and S_data_r
    // Hence in this case they corresponds to nt=0
    update_vel(0, 1, 0.0);
    // Copy S_data_vel["V"] to data_v[0]
    set_from_sacado_vel("V",0);
    set_from_sacado_vel("U",0);
    // This will predict a S_data_rp that is used for the
    // first time step size computation.
    // Using 0, 0, only S_data_rp will be computed.
    compute_rad(0, 0);
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
  // Predicting R'. Cannot use S_data_rp, as some of calc_dt
  // is called when dtn is not valid, hence also S_data_rp
  DType dr = (*H_data_vel["V"][nt])[varsizes["V"]-1];

  for ( int i = 0; i < Nx+1; i++ ) {
    DType v_max = std::fabs( x_f[i] * dr / data_r[nt] );
    for ( std::vector<std::string>::iterator it = velnames.begin();
          it != velnames.end(); it++ ) {
      v_max += std::fabs( (*H_data_vel[*it][nt])[i] / data_r[nt] );
    }
    dtn_trial = min( dtn_trial, tdSetup->cfl * dx / v_max );
    dtn_trial = min( dtn_trial, tdSetup->dTmax );
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

  calc_total_inc_viol( nt );

  if ( tdSetup->silence_output ) return;

  FILE* ot_file = fopen( tdInput->radiusFile, "a" );
  fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\t%10.20f\n", tsol, data_r[nt], 
      (*H_data_vel["V"][nt])[varsizes["V"]-1], (*H_data_vel["U"][nt])[varsizes["U"]-1] );
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
    fprintf( ot_file, "%10.20f\t", (*data_inc[0])[i] );
    for ( std::vector<std::string>::iterator it = chemnames.begin();
          it != chemnames.end(); it++ )
      fprintf( ot_file, "%10.20f\t", (*H_data_var[*it][nt])[i] );
    fprintf( ot_file, "\n" );
  }

  for ( int i = 0; i < Nx+1; i++ ) {
    fprintf( ot_file, "%10.20f\t", x_f[i] * data_r[nt] );
    for ( std::vector<std::string>::iterator it = velnames.begin();
          it != velnames.end(); it++ )
      fprintf( ot_file, "%10.20f\t", (*H_data_vel[*it][nt])[i] );
    fprintf( ot_file, "\n" );
  }

  fclose( ot_file );

  return;
}

void TDSolveSacado::diagnostic(const int nt) {
  if ( tdSetup->silence_screen ) return;
  fprintf( stdout, "Iter %3d: Time %6.8f; Time step %6.8f; CFL %2.4f; Radius %6.8f; Front-velocity %6.8f.\n",
      iter, tsol, dtn, cur_cfl, data_r[nt], S_data_rp.val() );
  return;
}

void TDSolveSacado::alloc(void) {
  Ndata = 3;
  
  alloc_array_full( data_g, Nx, Ndata );
  alloc_array_full( data_n, Nx, Ndata ); // M
  alloc_array_full( data_v, Nx+1, Ndata );
  alloc_array_full( data_u, Nx+1, Ndata );

  alloc_array_full( data_inc, Nx, 1 );

  H_data_var.insert(std::pair<std::string,HVEC**>("G",data_g));
  H_data_var.insert(std::pair<std::string,HVEC**>("N",data_n));

  H_data_vel.insert(std::pair<std::string,HVEC**>("V",data_v));
  H_data_vel.insert(std::pair<std::string,HVEC**>("U",data_u));

  S_data_num += Nx;   // _g
  S_data_num += Nx;   // _n

  data_r  = new DType[Ndata];

  if ( tdSetup->prob_id == 0 || tdSetup->prob_id == 1 ) {
    alloc_array_full( data_h, Nx, Ndata );
    H_data_var.insert(std::pair<std::string,HVEC**>("H",data_h));
    S_data_num += Nx;   // _h

    alloc_array_full( data_a, Na, Ndata ); 
    H_data_var.insert(std::pair<std::string,HVEC**>("A",data_a));
    S_data_num += Na;   // _a
  }

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
  // Time argument for prescribed functions
  DType shift_dt_space = shift_dt_to;
  switch ( tdSetup->assemb ) {
    case 0: // explicit (no Jacobian)
      numIter = 1;
      assembleJacobian = false;
      shift_dt_space = shift_dt_fr;
      break;
    case 1: // implicit
      numIter = tdSetup->ntIter;
      assembleJacobian = true;
      break;
    case 2: // PMC - just add for fun, not really using
      numIter = tdSetup->ntIter;
      assembleJacobian = false;
      shift_dt_space = shift_dt_fr;
      break;
    default:
      break;
  }

  // ! All fluxes are on the left-hand side, including the diffusion ones
  // ! All source terms are on the right-hand side.
  for ( int iter = 0; iter < numIter; iter++ ) {
    // Step 0. Setup the Sacado variables
    reset_sacado_ind_array();
    reset_sacado_dep_array();
    // Set the independent variables iterates to [to] 
    set_sacado_ind( to );

    // Step 1. Compute S_data_r and S_data_rp.
    //   assemb == 0 : S_data_r = data_r[fr], S_data_rp ~ data_v[fr]
    //   assemb == 1 : S_data_r = data_r[to], S_data_rp ~ data_v[to]
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
    update_vel( fr, to, shift_dt_space, false );  

    // Step 3. Compute the geometric terms for both cells and chems
    // This will complete
    //   S_srce_cell["X"]
    // and update
    //   S_srce_chem["X"]
    // Note this one does not use S_data_r, as it may point to either fr or to.
    // Since implicit method is always used for chemicals, we use data_r[to] here
    //
    // Problem-dependent source for chems are computed after cells are updated
    update_source_geom( to );

    // Step 4. Compute all remaining fluxes for cells. In order:
    // 4.1 - S_flux_cell_v["X"]
    compute_cells_vflux();
    // 4.2 - S_lim_ppm_cell[i]["X"], i=0,1
    //       S_lim_ppm_sync[i], i=0,1
    //       S_flux_cell_r["X"] 
    compute_cells_rflux();
    // 4.3 - S_flux_cell
    compute_cells_flux();

    // Step 5. Compute all residuals for cells
    // The results are stored in H_data_res
    // Note that dt is multuplied to the spatial terms
    reset_residual();
    compute_cells_residual( fr, to );

    // Step 6. Update the solution/iterates for cells
    if ( assembleJacobian ) {
      // Assemble the Jacobian, solve for one Newton iteration
      // 6.1 - Implicit updates of cell numbers 
      update_cells_by_res_imp( to );
    }
    else {
      // This will use data_r instead of S_data_r
      // 6.1 - Explicit updates of cell numbers
      update_cells_by_res( to );
    }

    // Step 7. Now that cells are updated, compute the
    //         problem-dependent source for chems
    // This will complete
    //   S_srce_chem["X"]
    // using updated cell numbers
    set_sacado_ind( cellnames, to );
    compute_source_chem_prob( to, tsol+shift_dt_to*dtn );

    // Step 8. Compute the fluxes for chems. In order:
    // 8.1 - S_flux_chem_d["X"]
    //     - it is independent of R, hence no passing of nt
    compute_chems_dflux();
    // 8.2 - S_flux_chem_r["X"]
    compute_chems_rflux( to );
    // 8.2 - S_flux_chem["X"]
    compute_chems_flux();

    // Step 9. Compute all residuals for chems
    // The results are stored in H_data_res
    // Note that dt is multuplied to the spatial terms
    compute_chems_residual( fr, to );

    // Step X. Update the solution/iterates for chems
    // X-1. Implicit updates of chemicals
    update_chems_by_res_imp( to );
    // X-2. Make consistent the velocity for dtn
    reset_sacado_ind_array();
    reset_sacado_dep_array();
    set_sacado_ind( to );
    set_sacado_r( to );
    update_vel( to, to, shift_dt_to );
    set_from_sacado_vel( "V", to );
    set_from_sacado_vel( "U", to );
    compute_rad( to, to );
  }

  return true;
}

void TDSolveSacado::copy_data(const int fr, const int to) {
  for ( std::map<std::string,HVEC**>::iterator it = H_data_var.begin();
        it != H_data_var.end(); it++ ) {
    *it->second[to] = *it->second[fr];
  }

  for ( std::map<std::string,HVEC**>::iterator it = H_data_vel.begin();
        it != H_data_vel.end(); it++ ) {
    *it->second[to] = *it->second[fr];
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

void TDSolveSacado::setup_chem_par( std::string chemname, 
    const DType decay, const DType diffusion, std::string cellname, 
    const DType mm_ceil, const DType mm_half, const bool dbc ) {
  std::map<std::string, ChemPar>::iterator it = chempars.find( chemname );
  if ( it != chempars.end() ) {
    std::string error_message = "TDSolveSacado::setup_chem_par - chemoattractant/opn parameter for " + chemname + " has already been setup.\n";
    throw std::runtime_error( error_message );
  }
  ChemPar chemPar;
  chemPar.decay        = decay;
  chemPar.diffusion    = diffusion;
  chemPar.mm_ceil      = mm_ceil;
  chemPar.mm_half      = mm_half;
  chemPar.driver       = cellname;
  chemPar.dirichlet_bc = dbc;
  chempars.insert( std::pair<std::string, ChemPar>( chemname, chemPar ) );
  return;
}

void TDSolveSacado::setic_unif_cell(std::string name, DType val, int nt) {
  *H_data_var[name][nt] = val;
  return;
}

void TDSolveSacado::setic_expd_chem(DType multip, DType decay, int nt) {
  for ( int i = 0; i < Na; i++ ) {
    DType tvar = x_a_c[i] * data_r[nt];
    (*data_a[nt])[i] = multip * std::exp( - decay * tvar * tvar );
  }
  return;
}

void TDSolveSacado::setic_polyexp_chem(DType poly_0, DType poly_2, int nt) {
  for ( int i = 0; i < Nx; i++ ) {
    (*data_a[nt])[i] = poly_0 - poly_2 * x_a_c[i] * x_a_c[i];
  }
  DType exp_coef = poly_0 - poly_2;
  DType exp_rate = 2*poly_2/exp_coef;
  for ( int i = Nx; i < Na; i++ ) {
    (*data_a[nt])[i] = exp_coef * exp(-exp_rate*(x_a_c[i]-1.0));
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
    alloc_sacado_dep_array(S_array, varsizes[*it_name]+1);
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
  // Reconstructed interval face value for PPM flux
  for ( it_name = cellnames.begin(); it_name != cellnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, varsizes[*it_name]+1);
    S_rec_ppm_cell.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
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
    alloc_sacado_dep_array(S_array, varsizes[*it_name]+1);
    S_flux_chem.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  for ( it_name = chemnames.begin(); it_name != chemnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, varsizes[*it_name]+1);
    S_flux_chem_r.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  for ( it_name = chemnames.begin(); it_name != chemnames.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, varsizes[*it_name]+1);
    S_flux_chem_d.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );
  }
  return;
}

void TDSolveSacado::alloc_functions(void) {
  F_calcSrceMap.clear();
  if ( tdSetup->prob_id == 0 ) { // PDGF
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("G", &TDSolveSacado::calc_srce_g_pdgf) );
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("H", &TDSolveSacado::calc_srce_h_pdgf) );
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("N", &TDSolveSacado::calc_srce_n_pdgf) );
  }
  else if ( tdSetup->prob_id == 1 ) { // OPN
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("G", &TDSolveSacado::calc_srce_g_opn) );
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("H", &TDSolveSacado::calc_srce_h_opn) );
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("N", &TDSolveSacado::calc_srce_n_opn) );
  }
  else if ( tdSetup->prob_id == -1 ) { // Model problem
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("G", &TDSolveSacado::calc_srce_g_model) );
    F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("N", &TDSolveSacado::calc_srce_n_model) );
  }
  else {
    for ( std::vector<std::string>::iterator it = cellnames.begin(); 
          it != cellnames.end(); it++ ) {
      F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>(*it, &TDSolveSacado::calc_srce_cell_base) );
    }
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
  //std::vector<std::string>::iterator it;
  //for ( it = cellnames.begin(); it != cellnames.end(); it++ ) {
  //  set_sacado_ind_array(*it, nt);
  //}
  //for ( it = chemnames.begin(); it != chemnames.end(); it++ ) {
  //  set_sacado_ind_array(*it, nt);
  //}
  set_sacado_ind( cellnames, nt );
  set_sacado_ind( chemnames, nt );
  return;
}

void TDSolveSacado::set_sacado_ind(std::vector<std::string>& names, const int nt) {
  std::vector<std::string>::iterator it;
  for ( it = names.begin(); it != names.end(); it++ ) {
    set_sacado_ind_array(*it, nt);
  }
  return;
}

void TDSolveSacado::set_sacado_ind_array(std::string name, const int nt) {
  std::map<std::string, SacadoIndVarArray>::iterator it = S_data_var.find( name );
  if ( it == S_data_var.end() ) return;
  SacadoIndVarArray& s_var = it->second;
  int s_offset = S_data_offset[it->first];
  HVEC& data = *H_data_var[it->first][nt];
  int num = s_var.size();
  for ( int i = 0; i < num; i++ ) {
    s_var[i] = Sacado::Fad::DFad<DType>(S_data_num,s_offset+i,data[i]);
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
  std::vector<std::string>::iterator it;
  for ( it = cellnames.begin(); it != cellnames.end(); it++ ) {
    set_from_sacado_ind_array( *it, *H_data_var[*it][nt] );
  }
  for ( it = chemnames.begin(); it != chemnames.end(); it++ ) {
    set_from_sacado_ind_array( *it, *H_data_var[*it][nt] );
  }
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

void TDSolveSacado::set_from_sacado_vel( std::string name, const int nt) {
  SacadoDepVarArray& s_vel = S_data_vel[name];
  int num = varsizes[name];
  for ( int i = 0; i < num; i++ ) {
    (*H_data_vel[name][nt])[i] = s_vel[i].val();
  }
  return;
}

void TDSolveSacado::reset_sacado_ind_array(void) {
  for ( std::map<std::string, SacadoIndVarArray>::iterator it = S_data_var.begin();
        it != S_data_var.end(); it++ ) {
    reset_sacado_ind_array( it->second, S_data_offset[it->first], 0.0 );
  }
  return;
}

void TDSolveSacado::reset_sacado_ind_array(SacadoIndVarArray& s_var, const int s_offset, DType val) {
  for ( int i = 0; i < s_var.size(); i++ ) {
    s_var[i] = Sacado::Fad::DFad<DType>(S_data_num, s_offset+i, val);
  }
  return;
}

void TDSolveSacado::reset_sacado_dep_array(void) {
  S_data_r  = 0.0;
  S_data_rp = 0.0;
  reset_sacado_dep_array( S_data_vel, 0.0 );
  reset_sacado_dep_array( S_flux_cell, 0.0 );
  reset_sacado_dep_array( S_flux_cell_r, 0.0 );
  reset_sacado_dep_array( S_flux_cell_v, 0.0 );
  reset_sacado_dep_array( S_srce_cell, 0.0 );
  reset_sacado_dep_array( S_rec_ppm_cell, 0.0 );
  for ( int i = 0; i < 2; i++ ) {
    reset_sacado_dep_array( S_lim_ppm_cell[i], 1.0 );
    reset_sacado_dep_array( S_lim_ppm_sync[i], 1.0 );
  }
  reset_sacado_dep_array( S_srce_chem, 0.0 );
  reset_sacado_dep_array( S_flux_chem, 0.0 );
  reset_sacado_dep_array( S_flux_chem_r, 0.0 );
  reset_sacado_dep_array( S_flux_chem_d, 0.0 );
  return;
}

void TDSolveSacado::reset_sacado_dep_array(SacadoDepVarMapMap& s_maps, DType val) {
  for ( SacadoDepVarMapMap::iterator it = s_maps.begin();
        it != s_maps.end(); it++ ) {
    reset_sacado_dep_array( it->second, val );
  }
  return;
}

void TDSolveSacado::reset_sacado_dep_array(SacadoDepVarMap& s_map, DType val) {
  for ( SacadoDepVarMap::iterator it = s_map.begin();
        it != s_map.end(); it++ ) {
    reset_sacado_dep_array( it->second, val );
  }
  return;
}

void TDSolveSacado::reset_sacado_dep_array(SacadoDepVarArray& s_var, DType val) {
  for ( int i = 0; i < s_var.size(); i++ ) {
    s_var[i] = val;
  }
  return;
}

void TDSolveSacado::sum_in_sacado_dep_array(SacadoDepVarArray& s_tot, DType multip, SacadoDepVarArray& s_var) {
  int num = s_tot.size();
  for ( int i = 0; i < num; i++ ) {
    s_tot[i] = s_tot[i] + multip * s_var[i];
  }
  return;
}

void TDSolveSacado::print_sacado_ind_array(SacadoIndVarArray& s_var) {
  int num = s_var.size();
  cout << "Begin a new independent Sacado array:" << std::endl;
  for ( int i = 0; i < num; i++ ) {
    cout << s_var[i] << std::endl;
  }
  return;
}

void TDSolveSacado::print_sacado_dep_array(SacadoDepVarArray& s_var) {
  int num = s_var.size();
  cout << "Begin a new dependent Sacado array:" << std::endl;
  for ( int i = 0; i < num; i++ ) {
    cout << s_var[i] << std::endl;
  }
  return;
}

void TDSolveSacado::calc_total_inc_viol(const int nt) {
  *data_inc[0] = 1.0;
  for ( std::vector<std::string>::iterator it = cellnames.begin();
        it != cellnames.end(); it++ ) {
    *data_inc[0] -= *H_data_var[*it][nt];
  }
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
  int num = varsizes["U"];
  DType* l_x_f = varcoord["U"];

  reset_sacado_dep_array(s_u, 0.0);

  switch ( tdSetup->infil_u ) {
    case 0: // Test 0, u is always zero
      break;
    case 1: // Test 1, u is only a function of x
      break;
    case 2: // Test 2, space homogeneous solutions
      {
        DType v_par = tdSetup->init_u_par;
        s_u[0] = 0.0;
        for ( int i = 1; i < num; i++ )
          s_u[i] = - l_x_f[i] * S_data_r * v_par / ( 1.0 - std::exp(-3.0*v_par*tsol_loc) * (1.0-tdCoef->n_bc));
      }
      break;
    case 3: // Test 3, linear growth
      {
        DType v_par = tdSetup->init_u_par;
        s_u[0] = 0.0;
        for ( int i = 1; i < num; i++ )
          s_u[i] = -  v_par * std::sin( 0.5 * PI * l_x_f[i] ) / tdCoef->n_bc;
      }
      break;
    case 4: // Test 4, non-monotone radius change
      {
        DType v_par = tdSetup->init_u_par;
        s_u[0] = 0.0;
        for ( int i = 1; i < num; i++ )
          s_u[i] = v_par * std::sin( 4.0 * l_x_f[i] * S_data_r * ( 1.0 + tsol_loc ) );
      }
      break;
    case 5:
      {
        DType v_par = tdSetup->init_u_par;
        s_u[0] = 0.0;
        for ( int i = 0; i < num; i++ ) {
          SacadoDepVar s_r = l_x_f[i] * S_data_r;
          s_u[i] = - 2.0 * v_par * s_r * std::exp( - v_par * s_r * s_r );
        }
      }
      break;
  }

  return;
}

void TDSolveSacado::compute_vel_u_chem(void) {
  SacadoDepVarArray& s_u = S_data_vel["U"];
  SacadoIndVarArray& s_a = S_data_var["A"];

  DType alpha = tdCoef->alpha;
  SacadoDepVar aor  = alpha/S_data_r;

  s_u[0] = 0.0;
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
      s_v[i+1] += dx * (it->second)[i];
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
  SacadoDepVar s_loc_u;
  // Left boundary
  s_flux[0] = 0.0;
  // Interior points
  for ( int i = 1; i < Nx; i++ ) {
    if ( tdSetup->jac_vel ) s_loc_u = x_f[i] * x_f[i] * S_data_r * s_vel[i]; 
    else s_loc_u = x_f[i] * x_f[i] * S_data_r * s_vel[i].val();
    if ( s_loc_u >= 0.0 )
      s_flux[i] = s_loc_u * s_var[i-1];
    else
      s_flux[i] = s_loc_u * s_var[i];
  }
  // Right boundary
  if ( rdata == NULL )
    s_flux[Nx] = 0.0;
  else {
    if ( tdSetup->jac_vel ) s_loc_u = x_f[Nx] * x_f[Nx] * S_data_r * s_vel[Nx];
    else s_loc_u = x_f[Nx] * x_f[Nx] * S_data_r * s_vel[Nx].val();
    if ( s_loc_u >= 0.0 )
      s_flux[Nx] = s_loc_u * s_var[Nx-1];
    else
      s_flux[Nx] = s_loc_u * (*rdata);
  }
  return;
}

void TDSolveSacado::compute_cells_vflux(void) {
  for ( std::vector<std::string>::iterator it = cellnames.begin(); 
        it != cellnames.end(); it++ ) {
    // The third argument is not used at all
    compute_cells_varflux( "V", *it, 0.0, false );
  }
  return;
}

void TDSolveSacado::LimFun_PPM_BC(DType dval_c, DType dval_r, SacadoDepVar& s_lim) {
  if ( std::fabs(dval_c)<TOL_ZERO && std::fabs(dval_r)<TOL_ZERO )
    s_lim = 1.0;
  else if ( dval_c * dval_r <= 0.0 || 3.0*std::fabs(dval_r) <= 8.0*std::fabs(dval_c) )
    s_lim = 0.0;
  else
    s_lim = min(1.0, 5.0/(std::fabs(dval_r)/std::fabs(dval_c)-1.0));
  return;
}

void TDSolveSacado::LimFun_PPM_BC(SacadoDepVar& s_dval_c, SacadoDepVar& s_dval_r, SacadoDepVar& s_lim) {
  if ( Sacado::Fad::fabs(s_dval_c)<TOL_ZERO && Sacado::Fad::fabs(s_dval_r)<TOL_ZERO )
    s_lim = 1.0;
  else if ( s_dval_c * s_dval_r <= 0.0 || 3.0*Sacado::Fad::fabs(s_dval_r) <= 8.0*Sacado::Fad::fabs(s_dval_c) )
    s_lim = 0.0;
  else
    s_lim = Sacado::Fad::min(1.0, 5.0/(Sacado::Fad::fabs(s_dval_r)/Sacado::Fad::fabs(s_dval_c)-1.0));
  return;
}

void TDSolveSacado::LimFun_PPM_INT(DType alpha1, DType alpha2, SacadoDepVar& s_lim_l, SacadoDepVar& s_lim_r) {
  if ( alpha2 > alpha1 + TOL_ZERO ) {
    s_lim_l = 0.0; s_lim_r = 0.0;
  } else if ( alpha1 < 1.0 ) {
    s_lim_l = 1.0; s_lim_r = alpha1;
  } else if ( alpha2 > 1.0 ) {
    s_lim_l = 1.0/alpha2; s_lim_r = 1.0;
  } else {
    s_lim_l = 1.0; s_lim_r = 1.0;
  }
  return;
}

void TDSolveSacado::LimFun_PPM_INT(SacadoDepVar& s_alpha1, SacadoDepVar& s_alpha2, SacadoDepVar& s_lim_l, SacadoDepVar& s_lim_r) {
  if ( s_alpha2 > s_alpha1 + TOL_ZERO ) {
    s_lim_l = 0.0; s_lim_r = 0.0;
  } else if ( s_alpha1 < 1.0 ) {
    s_lim_l = 1.0; s_lim_r = s_alpha1;
  } else if ( s_alpha2 > 1.0 ) {
    s_lim_l = 1.0/s_alpha2; s_lim_r = 1.0;
  } else {
    s_lim_l = 1.0; s_lim_r = 1.0;
  }
  return;
}

void TDSolveSacado::compute_cells_rflux(void) {
  // Preparation for the synced limiter
  std::vector<bool> set_zero(S_lim_ppm_sync[0].size(), false);
  SacadoDepVarArray s_ppm_alpha[2];
  alloc_sacado_dep_array( s_ppm_alpha[0], S_lim_ppm_sync[0].size() );
  alloc_sacado_dep_array( s_ppm_alpha[1], S_lim_ppm_sync[1].size() );
  reset_sacado_dep_array( s_ppm_alpha[0], DTYPE_MAX );
  reset_sacado_dep_array( s_ppm_alpha[1], 0.0 );
  
  // Compute individual limiter for cell species
  for ( std::vector<std::string>::iterator it = cellnames.begin();
        it != cellnames.end(); it++ ) {
    SacadoIndVarArray& s_var = S_data_var[*it];
    SacadoDepVarArray& s_rec = S_rec_ppm_cell[*it];
    SacadoDepVarArray& s_lim_l = S_lim_ppm_cell[0][*it];
    SacadoDepVarArray& s_lim_r = S_lim_ppm_cell[1][*it];
    // Passing -1 as nt, so that S_data_r is used
    calc_ppm_lim_var( s_var, s_rec, s_lim_l, s_lim_r, varcoord[*it], -1, true, &set_zero, s_ppm_alpha);
  }

  // Compute the synced limiter
  calc_ppm_lim_sync( S_lim_ppm_sync, S_lim_ppm_cell, set_zero, s_ppm_alpha );

  // Compute the fluxes
  for ( std::vector<std::string>::iterator it = cellnames.begin();
        it != cellnames.end(); it++ ) {
    SacadoIndVarArray& s_var = S_data_var[*it];
    SacadoDepVarArray& s_rec = S_rec_ppm_cell[*it];
    SacadoDepVarArray& s_flx = S_flux_cell_r[*it];
    calc_ppm_flux_var( s_var, s_rec, s_flx, S_lim_ppm_sync, varcoord[*it] );
  }
  return;
}

void TDSolveSacado::compute_cells_flux(void) {
  std::vector<std::string>::iterator it;
  // Flux associated with R'
  for ( it = cellnames.begin(); it != cellnames.end(); it++ ) {
    sum_in_sacado_dep_array( S_flux_cell[*it], 1.0, S_flux_cell_r[*it] );
  }
  // Flux associated with velocities
  for ( std::vector<std::string>::iterator it_vel = velnames.begin();
        it_vel != velnames.end(); it_vel++ ) {
    for ( it = cellnames.begin(); it != cellnames.end(); it++ ) {
      SacadoDepVarMap::iterator it_map = S_flux_cell_v[*it_vel].find(*it);
      if ( it_map == S_flux_cell_v[*it_vel].end() ) // Not all cells are associated with all velocities
        continue;
      sum_in_sacado_dep_array( S_flux_cell[*it], 1.0, S_flux_cell_v[*it_vel][*it] );
    }
  }
  return;
}

void TDSolveSacado::calc_ppm_lim_var(SacadoIndVarArray& s_var, SacadoDepVarArray& s_rec, SacadoDepVarArray& s_lim_l, SacadoDepVarArray& s_lim_r, DType* l_x_c, const int nt, const bool tosync, std::vector<bool>* set_zero, SacadoDepVarArray* s_alpha) {
  reset_sacado_dep_array( s_lim_l, 1.0 );
  reset_sacado_dep_array( s_lim_r, 1.0 );
  const DType sot = 7.0/12.0;
  const DType oot = 1.0/12.0;

  int num = s_var.size();
  SacadoDepVarArray s_bar;
  alloc_sacado_dep_array( s_bar, num );
  if ( nt == - 1 ) {
    for ( int i = 0; i < num; i++ ) {
      s_bar[i] = l_x_c[i] * l_x_c[i] * l_x_c[i] * S_data_r * S_data_r * s_var[i];
    }
  }
  else {
    for ( int i = 0; i < num; i++ ) {
      s_bar[i] = l_x_c[i] * l_x_c[i] * l_x_c[i] * data_r[nt] * data_r[nt] * s_var[i];
    }
  }

  // Reconstruct interface face value - near left boundary
  {
    s_rec[0] = 0.0;
    s_rec[1] = sot*(s_bar[0]+s_bar[1]) - oot*(s_bar[2]-s_bar[0]);
  }
  // - interior points
  for ( int i = 2; i < num-1; i++ )
    s_rec[i] = sot*(s_bar[i-1]+s_bar[i]) - oot*(s_bar[i-2]+s_bar[i+1]);
  // - near right boundary
  {
    s_rec[num-1] = oot*(3.0*s_bar[num-1]+13.0*s_bar[num-2]-5.0*s_bar[num-3]+s_bar[num-4]);
    s_rec[num] = oot*(25.0*s_bar[num-1]-23.0*s_bar[num-2]+13.0*s_bar[num-3]-3.0*s_bar[num-4]);
  }

  // Computing the limiters - left boundary 
  {
    s_lim_l[0] = 0.0;
    if ( tdSetup->jac_lim ) LimFun_PPM_BC( s_bar[0], s_rec[1], s_lim_r[0] );
    else LimFun_PPM_BC( s_bar[0].val(), s_rec[1].val(), s_lim_r[0] );
  }
  // - all but left boundary
  for ( int i = 1; i < num; i++ ) {
    SacadoDepVar s_diff_l = s_rec[i]   - s_bar[i];
    SacadoDepVar s_diff_r = s_rec[i+1] - s_bar[i];
    SacadoDepVar s_alpha1, s_alpha2;
    bool zero_lim = false;
    bool no_lim = std::fabs(s_diff_l.val()) + std::fabs(s_diff_r.val()) < 2.0*TOL_ZERO;
    if ( !no_lim ) {
      if ( std::fabs(s_diff_l.val())<=TOL_ZERO || std::fabs(s_diff_r.val())<=TOL_ZERO || s_diff_l.val()*s_diff_r.val() > TOL_ZERO ) {
        zero_lim = true;
        s_lim_l[i] = 0.0;
        s_lim_r[i] = 0.0;
      }
      else {
        s_alpha1 = 2.0*Sacado::Fad::fabs(s_diff_l)/Sacado::Fad::fabs(s_diff_r);
        s_alpha2 = 0.25 * s_alpha1;
        if ( tdSetup->jac_lim ) LimFun_PPM_INT( s_alpha1, s_alpha2, s_lim_l[i], s_lim_r[i] );
        else LimFun_PPM_INT( s_alpha1.val(), s_alpha2.val(), s_lim_l[i], s_lim_r[i] );
      }
    }
    else {
      s_lim_l[i] = 1.0;
      s_lim_r[i] = 1.0;
      s_alpha1 = DTYPE_MAX;
      s_alpha2 = 0.0;
    }
    if ( tosync ) {
      if ( zero_lim ) {
        (*set_zero)[i] = true;
      }
      else {
        s_alpha[0][i] = Sacado::Fad::min( s_alpha[0][i], s_alpha1 );
        s_alpha[1][i] = Sacado::Fad::max( s_alpha[1][i], s_alpha2 );
      }
    }
  }
  return;
}

void TDSolveSacado::calc_ppm_lim_sync(SacadoDepVarArray* s_sync, SacadoDepVarMap* s_lims, std::vector<bool>& set_zero, SacadoDepVarArray* s_alpha) {
  int num = s_sync[0].size();
  reset_sacado_dep_array( s_sync[0], 1.0 );
  reset_sacado_dep_array( s_sync[1], 1.0 );

  // - left boundary
  {
    s_sync[0][0] = 0.0;
    for ( SacadoDepVarMap::iterator it = s_lims[1].begin(); it != s_lims[1].end(); it++ ) {
      s_sync[1][0] = Sacado::Fad::min( s_sync[1][0], it->second[0] );
    }
  }
  // - all but left boundary
  for ( int i = 1; i < num; i++ ) {
    if ( set_zero[i] ) {
      s_sync[0][i] = 0.0;
      s_sync[1][i] = 0.0;
    }
    else {
      if ( tdSetup->jac_lim ) LimFun_PPM_INT( s_alpha[0][i], s_alpha[1][i], s_sync[0][i], s_sync[1][i] );
      else LimFun_PPM_INT( s_alpha[0][i].val(), s_alpha[1][i].val(), s_sync[0][i], s_sync[1][i] );
    }
  }

  return;
}

void TDSolveSacado::calc_ppm_flux_var(SacadoIndVarArray& s_var, SacadoDepVarArray& s_rec, SacadoDepVarArray& s_flx, SacadoDepVarArray* s_lim, DType* l_x_c, const int nt) {
  int num = s_var.size();

  SacadoDepVarArray s_bar;
  alloc_sacado_dep_array( s_bar, num );
  if ( nt == -1 ) {
    for ( int i = 0; i < num; i++ ) {
      s_bar[i] = l_x_c[i] * l_x_c[i] * l_x_c[i] * S_data_r * S_data_r * s_var[i];
    }
  }
  else {
    for ( int i = 0; i < num; i++ ) {
      s_bar[i] = l_x_c[i] * l_x_c[i] * l_x_c[i] * data_r[nt] * data_r[nt] * s_var[i];
    }
  }

  // - left boundary
  {
    s_flx[0] = 0.0;
  }
  // - interior points
  SacadoDepVar loc_vel;
  if ( nt == -1 )
    loc_vel = -S_data_rp/S_data_r;
  else
    loc_vel = -S_data_rp/data_r[nt];
  for ( int i = 1; i < num; i++ ) {
    SacadoDepVar lim_l = s_lim[1][i-1];
    SacadoDepVar lim_r = s_lim[0][i];
    if ( loc_vel >= 0.0 ) 
      s_flx[i] = loc_vel * ( (1.0-lim_l) * s_bar[i-1] + lim_l * s_rec[i] );
    else
      s_flx[i] = loc_vel * ( (1.0-lim_r) * s_bar[i] + lim_r * s_rec[i] );
  }
  // - right boundary
  {
    s_flx[num] = 0.0;
  }

  return;
}

void TDSolveSacado::compute_chems_dflux(void) {
  for ( std::vector<std::string>::iterator it = chemnames.begin();
        it != chemnames.end(); it++ ) {
    calc_dflux_chem_var( *it );
  }
  return;
}

void TDSolveSacado::calc_dflux_chem_var(std::string chemname) {
  SacadoIndVarArray& s_var = S_data_var[chemname];
  SacadoDepVarArray& s_flx = S_flux_chem_d[chemname];
  reset_sacado_dep_array(s_flx);
  DType* l_x_f = varflxcoord[chemname];
  int   num   = s_var.size();
  DType dcoef = chempars[chemname].diffusion;
  bool  dbc   = chempars[chemname].dirichlet_bc;

  // - left boundary
  {
    s_flx[0] = 0.0;
  }
  // - interior points
  for ( int i = 1; i < num; i++ ) {
    s_flx[i] = - dcoef * l_x_f[i] * l_x_f[i] * (s_var[i] - s_var[i-1]) / dx;
  }
  // - right boundary
  {
    if ( dbc ) // Homogeneous D. b.c.
      s_flx[num] = - dcoef * l_x_f[num] * l_x_f[num] * (0.0 - s_var[num-1]) / dx;
    else // Homogeneous N. b.c.
      s_flx[num] = 0.0;
  }

  return;
}

void TDSolveSacado::compute_chems_rflux(const int nt) {
  for ( std::vector<std::string>::iterator it = chemnames.begin();
        it != chemnames.end(); it++ ) {
    if ( tdSetup->chmppm ) calc_rflux_chem_var_ppm( *it, nt );
    else calc_rflux_chem_var_upw( *it, nt );
  }
  return;
}

void TDSolveSacado::calc_rflux_chem_var_upw(std::string chemname, const int nt) {
  SacadoIndVarArray& s_var = S_data_var[chemname];
  SacadoDepVarArray& s_flx = S_flux_chem_r[chemname];
  reset_sacado_dep_array(s_flx);
  DType* l_x_c = varcoord[chemname];
  DType* l_x_f = varflxcoord[chemname];
  int  num = s_var.size();
  bool dbc = chempars[chemname].dirichlet_bc;

  // - left boundary
  {
    s_flx[0] = 0.0;
  }
  // - interior points
  for ( int i = 1; i < num; i++ ) {
    SacadoDepVar loc_vel_r = - l_x_f[i] * S_data_rp * data_r[nt];
    int offset = ( loc_vel_r >= 0.0 ) ? -1 : 0;
    int idx = i + offset;
    s_flx[i] = loc_vel_r * l_x_c[idx] * l_x_c[idx] * s_var[idx];
  }
  // - right boundary
  {
    SacadoDepVar loc_vel_r = - l_x_f[num] * S_data_rp * data_r[nt];
    if ( dbc ) { // Homogeneous D. b.c.
      if ( loc_vel_r >= 0.0 )
        s_flx[num] = loc_vel_r * l_x_c[num-1] * l_x_c[num-1] * s_var[num-1];
      else
        s_flx[num] = 0.0;
    }
    else { // Homogeneous N. b.c.
      DType loc_x = ( loc_vel_r >= 0.0 ) ? l_x_c[num-1] : 2.0*l_x_c[num-1]-l_x_c[num-2];
      s_flx[num] = loc_vel_r * loc_x * loc_x * s_var[num-1];
    }
  }

  return;
}

void TDSolveSacado::calc_rflux_chem_var_ppm(std::string chemname, const int nt) {
  SacadoIndVarArray& s_var = S_data_var[chemname];
  SacadoDepVarArray& s_flx = S_flux_chem_r[chemname];

  SacadoDepVarArray  s_lim[2];
  SacadoDepVarArray  s_rec;
  alloc_sacado_dep_array( s_lim[0], s_var.size() );
  alloc_sacado_dep_array( s_lim[1], s_var.size() );
  alloc_sacado_dep_array( s_rec, s_var.size()+1 );
  calc_ppm_lim_var( s_var, s_rec, s_lim[0], s_lim[1], varcoord[chemname], nt, false );

  calc_ppm_flux_var( s_var, s_rec, s_flx, s_lim, varcoord[chemname], nt );

  return;
}

void TDSolveSacado::compute_chems_flux(void) {
  std::vector<std::string>::iterator it;
  // Diffusive flux
  for ( it = chemnames.begin(); it != chemnames.end(); it++ ) {
    sum_in_sacado_dep_array( S_flux_chem[*it], 1.0, S_flux_chem_d[*it] );
  }
  // Flux associated with R'
  for ( it = chemnames.begin(); it != chemnames.end(); it++ ) {
    sum_in_sacado_dep_array( S_flux_chem[*it], 1.0, S_flux_chem_r[*it] );
  }
  return;
}

void TDSolveSacado::reset_residual(void) {
  for ( std::map<std::string,HVEC**>::iterator it = H_data_res.begin();
        it != H_data_res.end(); it++ ) {
    *it->second[0] = 0.0;
  }
  return;
}

void TDSolveSacado::compute_cells_residual(const int fr, const int to) {
  std::vector<std::string>::iterator it;
  for ( it = cellnames.begin(); it != cellnames.end(); it++ ) {
    HVEC& h_var = *H_data_var[*it][fr];
    HVEC& h_res = *H_data_res[*it][0];
    SacadoIndVarArray& s_var = S_data_var[*it];
    SacadoDepVarArray& s_flx = S_flux_cell[*it];
    SacadoDepVarArray& s_src = S_srce_cell[*it];
    h_res = 0.0;
    assemble_inertia_var( h_res, data_r[fr], h_var, data_r[to], s_var, varcoord[*it] );
    assemble_flux_var( h_res, s_flx );
    assemble_srce_var( h_res, s_src );
  }
  return;
}

void TDSolveSacado::assemble_inertia_var(HVEC& h_res, DType r_fr, HVEC& h_var, DType r_to, SacadoIndVarArray s_var, DType* l_x_c) {
  int num = s_var.size();
  DType r_fr2 = r_fr * r_fr;
  DType r_to2 = r_to * r_to;
  for ( int i = 0; i < num; i++ ) {
    h_res[i] += l_x_c[i] * l_x_c[i] * ( r_to2 * s_var[i].val() - r_fr2 * h_var[i] );
  }
  return;
}

void TDSolveSacado::assemble_flux_var(HVEC& h_res, SacadoDepVarArray& s_flx) {
  int num = s_flx.size() - 1;
  for ( int i = 0; i < num; i++ ) {
    h_res[i] += (dtn/dx) * (s_flx[i+1].val() - s_flx[i].val());
  }
  return;
}

void TDSolveSacado::assemble_srce_var(HVEC& h_res, SacadoDepVarArray& s_src) {
  int num = s_src.size();
  for ( int i = 0; i < num; i++ ) {
    h_res[i] -= dtn * s_src[i].val();
  }
  return;
}

void TDSolveSacado::compute_chems_residual(const int fr, const int to) {
  std::vector<std::string>::iterator it;
  for ( it = chemnames.begin(); it != chemnames.end(); it++ ) {
    HVEC& h_var = *H_data_var[*it][fr];
    HVEC& h_res = *H_data_res[*it][0];
    SacadoIndVarArray& s_var = S_data_var[*it];
    SacadoDepVarArray& s_flx = S_flux_chem[*it];
    SacadoDepVarArray& s_src = S_srce_chem[*it];
    h_res = 0.0;
    assemble_inertia_var( h_res, data_r[fr], h_var, data_r[to], s_var, varcoord[*it] );
    assemble_flux_var( h_res, s_flx );
    assemble_srce_var( h_res, s_src );
  }
  return;
}

void TDSolveSacado::update_cells_by_res(const int nt) {
  std::vector<std::string>::iterator it;
  for ( it = cellnames.begin(); it != cellnames.end(); it++ ) {
    HVEC& h_var = *H_data_var[*it][nt];
    HVEC& h_res = *H_data_res[*it][0];
    DType h_r   = data_r[nt];
    update_sol_by_res_var_exp( h_var, h_res, h_r, varcoord[*it], varsizes[*it] );
  }
  return;
}

void TDSolveSacado::update_sol_by_res_var_exp(HVEC& h_var, HVEC& h_res, DType h_r, DType* l_x_c, const int num) {
  DType h_r2 = h_r * h_r;
  for ( int i = 0; i < num; i++ ) {
    h_var[i] -= h_res[i]/h_r2/l_x_c[i]/l_x_c[i];
  }
  return;
}

void TDSolveSacado::update_cells_by_res_imp(const int nt) {
  cell_lspack.mat->reset( 0.0 );
  cell_lspack.rhs->reset( 0.0 );
  cell_lspack.sol->reset( 0.0 );

  std::vector<std::string>::iterator it;

  for ( it = cellnames.begin(); it != cellnames.end(); it++ ) {
    HVEC& h_res = *H_data_res[*it][0];
    DType h_r   = data_r[nt];
    int s_offset = S_data_offset[*it];
    int l_offset = H_cell_offset[*it];
    update_cell_linsys_by_res_var( h_res, h_r, varcoord[*it], S_flux_cell[*it], S_srce_cell[*it], varsizes[*it], s_offset, l_offset, cell_lspack );
  }

  cell_lspack.linSolve->solve();

  for ( it = cellnames.begin(); it != cellnames.end(); it++ ) {
    HVEC& h_var = *H_data_var[*it][nt];
    int l_offset = H_cell_offset[*it];
    int num = varsizes[*it];
    update_hsol_by_linsol( h_var, l_offset, num, cell_lspack );
  }

  return;
}

void TDSolveSacado::update_cell_linsys_by_res_var(HVEC& h_res, DType h_r, DType* l_x_c, SacadoDepVarArray& s_flx, SacadoDepVarArray& s_src, const int num, const int s_offset, const int l_offset, LSPack& lsPack) {
  VComp comp(1,0.0);
  Entry elem(1,0.0);

  // RHS
  DType s_r_to2 = h_r * h_r;
  for ( int i = 0; i < num; i++ ) {
    comp = h_res[i]/l_x_c[i]/l_x_c[i]/s_r_to2;
    lsPack.rhs->set(l_offset+i,comp);
  }

  // Matrix
  // - mass matrix
  for ( int i = 0; i < num; i++ ) {
    elem = 1.0;
    lsPack.mat->set(l_offset+i,l_offset+i,elem);
  }
  // - source terms - assuming no nearby data is needed!
  for ( int i = 0; i < num; i++ ) {
    for ( std::vector<std::string>::iterator it = cellnames.begin();
          it != cellnames.end(); it++ ) {
      int l_offset_cell = H_cell_offset[*it];
      int s_offset_cell = S_data_offset[*it];
      elem = lsPack.mat->get(l_offset+i,l_offset_cell+i);
      elem -= dtn * s_src[i].dx( s_offset_cell+i )/l_x_c[i]/l_x_c[i]/s_r_to2;
      lsPack.mat->set(l_offset+i,l_offset_cell+i,elem);
    }
  }
  // - fluxes 
  if ( tdSetup->jac_vel ) { 
    for ( int i = 0; i < num; i++ ) {
      DType denom = l_x_c[i]*l_x_c[i]*s_r_to2;
      for ( std::vector<std::string>::iterator it = cellnames.begin();
            it != cellnames.end(); it++ ) {
        int l_offset_cell = H_cell_offset[*it];
        int s_offset_cell = S_data_offset[*it];
        for ( int j = 0; j < num; j++ ) {
          elem = lsPack.mat->get(l_offset+i,l_offset_cell+j);
          elem += (dtn/dx) * (s_flx[i+1].dx( s_offset_cell+j ) - s_flx[i].dx( s_offset_cell+j ))/denom;
          lsPack.mat->set(l_offset+i,l_offset_cell+j,elem);
        }
      }
    }
  }
  else {
    // assuming the stencil within a centered stencil of:
    // - 3 intervals if upwind flux is used
    // - 7 intervals if ppm flux is used 
    int sten = ( tdSetup->chmppm ) ? 3 : 1;
    for ( int i = 0; i < num; i++ ) {
      int i_min = max(0,i-sten);
      DType denom = l_x_c[i]*l_x_c[i]*s_r_to2;
      int i_max = min(num-1,i+sten);
      for ( int j = i_min; j <= i_max; j++ ) {
        elem = lsPack.mat->get(l_offset+i,l_offset+j);
        elem += (dtn/dx) * (s_flx[i+1].dx( s_offset+j )-s_flx[i].dx( s_offset+j ))/denom;
        lsPack.mat->set(l_offset+i,l_offset+j,elem);
      }
    }
  }

  return;
}

void TDSolveSacado::update_hsol_by_linsol(HVEC& h_var, const int l_offset, const int num, LSPack& lsPack) {
  VComp comp(1,0.0);

  for ( int i = 0; i < num; i++ ) {
    comp = lsPack.sol->get(l_offset+i);
    HVEC vsol(1);
    vsol = comp.getVal();
    h_var[i] -= vsol[0];
  }

  return;
}

void TDSolveSacado::update_chems_by_res_imp(const int nt) {
  std::vector<std::string>::iterator it;
  for ( it = chemnames.begin(); it != chemnames.end(); it++ ) {
    HVEC& h_var = *H_data_var[*it][nt];
    HVEC& h_res = *H_data_res[*it][0];
    DType h_r   = data_r[nt];
    int s_offset = S_data_offset[*it];
    update_chem_by_res_var_imp( h_var, h_res, h_r, varcoord[*it], S_flux_chem[*it], S_srce_chem[*it], varsizes[*it], s_offset, varlspacks[*it] );
  }
  return;
}

void TDSolveSacado::update_chem_by_res_var_imp(HVEC& h_var, HVEC& h_res, DType h_r, DType* l_x_c, SacadoDepVarArray& s_flx, SacadoDepVarArray& s_src, const int num, const int s_offset, LSPack& lsPack) {
  lsPack.mat->reset( 0.0 );
  lsPack.rhs->reset( 0.0 );
  lsPack.sol->reset( 0.0 );

  VComp comp(1,0.0);
  Entry elem(1,0.0);

  // RHS
  DType s_r_to2 = h_r * h_r;
  for ( int i = 0; i < num; i++ ) {
    comp = h_res[i]/l_x_c[i]/l_x_c[i]/s_r_to2;
    lsPack.rhs->set(i,comp);
  }

  // Matrix
  // - mass matrix
  for ( int i = 0; i < num; i++ ) {
    elem = 1.0;
    lsPack.mat->set(i,i,elem);
  }
  // - source terms - assuming no nearby data is needed!
  for ( int i = 0; i < num; i++ ) {
    elem = lsPack.mat->get(i,i);
    elem -= dtn*s_src[i].dx( s_offset+i )/l_x_c[i]/l_x_c[i]/s_r_to2;
    lsPack.mat->set(i,i,elem);
  }
  // - fluxes - assuming the stencil within a centered stencil of:
  // - 3 intervals if upwind flux is used
  // - 7 intervals if ppm flux is used 
  int sten = ( tdSetup->chmppm ) ? 3 : 1;
  for ( int i = 0; i < num; i++ ) {
    int i_min = max(0,i-sten);
    DType denom = l_x_c[i]*l_x_c[i]*s_r_to2;
    int i_max = min(num-1,i+sten);
    for ( int j = i_min; j <= i_max; j++ ) {
      elem = lsPack.mat->get(i,j);
      elem += (dtn/dx) * (s_flx[i+1].dx( s_offset+j )-s_flx[i].dx( s_offset+j ))/denom;
      lsPack.mat->set(i,j,elem);
    }
  }

  // Solve the linear system
  lsPack.linSolve->solve();

  // Fetch and update solution
  for ( int i = 0; i < num; i++ ) {
    comp = lsPack.sol->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    h_var[i] -= vsol[0];
  }

  return;
}

void TDSolveSacado::compute_rad(const int fr, const int to) {
  if ( tdSetup->dr_opt == 2 ) // For legacy issue, there is no opt 1
    S_data_rp = (*data_v[to])[Nx]/(1.0-0.25*dx*dx);
  else
    S_data_rp = (*data_v[to])[Nx];
  
  // Only S_data_rp is required
  if ( fr == to )
    return;

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

  if ( tdSetup->assemb == 1 )
    S_data_r = data_r[to];
  else
    S_data_r = data_r[fr];

  return;
}

void TDSolveSacado::update_source_geom(int nt) {
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
    calc_srce_geom( s_var, s_srce, varcoord[*it], nt );
  }
  return;
}

void TDSolveSacado::calc_srce_geom(SacadoIndVarArray& s_var, SacadoIndVarArray& s_srce, DType* l_x_c) {
  int num = s_var.size(); // This could be either Nx or Na
  for ( int i = 0; i < num; i++ ) {
    s_srce[i] -= l_x_c[i] * l_x_c[i] * S_data_rp * S_data_r * s_var[i]; 
  }
  return;
}

void TDSolveSacado::calc_srce_geom(SacadoIndVarArray& s_var, SacadoIndVarArray& s_srce, DType* l_x_c, const int nt) {
  int num = s_var.size(); // This could be either Nx or Na
  for ( int i = 0; i < num; i++ ) {
    s_srce[i] -= l_x_c[i] * l_x_c[i] * S_data_rp * data_r[nt] * s_var[i]; 
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

void TDSolveSacado::calc_srce_cell_base(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  reset_sacado_dep_array(s_srce);
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

  DType delta_n = tdCoef->delta_n;
  SacadoIndVarArray& s_n = S_data_var["N"];

  for ( int i = 0; i < Nx; i++ ) {
    s_srce[i] += (x_c[i] * x_c[i] * S_data_r * S_data_r) * ( - delta_n * s_n[i] );
  }

  return;
}

void TDSolveSacado::calc_srce_g_opn(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  calc_srce_g_base( s_srce, tsol_loc );

  DType lambda  = tdCoef->lambda;
  DType omega   = tdCoef->omega;
  DType a       = tdCoef->a;
  DType b       = tdCoef->b;
  DType a_prime = tdCoef->a_prime;
  DType b_prime = tdCoef->b_prime;
  DType mu      = tdCoef->mu;
  DType mu_coef = tdCoef->mu_coef;
  DType mu_slp  = tdCoef->mu_slp;
  DType p       = tdCoef->p;
  SacadoIndVarArray& s_g  = S_data_var["G"];
  SacadoIndVarArray& s_n  = S_data_var["N"];
  SacadoIndVarArray& s_p1 = S_data_var["P1"];
  SacadoIndVarArray& s_p2 = S_data_var["P2"];

  SacadoDepVar s_lambda;
  SacadoDepVar s_mu;

  for ( int i = 0; i < Nx; i++ ) {
    s_lambda = icalc_scal_logist_mm(lambda,s_g[i],1.0,b_prime,b,s_p2[i]);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * std::pow(x_c[i],p) );
    s_mu = icalc_scal_logist_mm(loc_mu,0.0,1.0,a_prime,a,s_p1[i]);
    s_srce[i] += (x_c[i] * x_c[i] * S_data_r * S_data_r) * ( (s_lambda-s_mu-omega*s_n[i]) * s_g[i] );
  }

  return;
}

void TDSolveSacado::calc_srce_h_opn(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  calc_srce_h_base( s_srce, tsol_loc );

  DType a       = tdCoef->a;
  DType a_prime = tdCoef->a_prime;
  DType mu      = tdCoef->mu;
  DType mu_coef = tdCoef->mu_coef;
  DType mu_slp  = tdCoef->mu_slp;
  DType p       = tdCoef->p;
  DType delta_h = tdCoef->delta_h;
  SacadoIndVarArray& s_g  = S_data_var["G"];
  SacadoIndVarArray& s_h  = S_data_var["H"];
  SacadoIndVarArray& s_p1 = S_data_var["P1"];

  SacadoDepVar s_mu;

  for ( int i = 0; i < Nx; i++ ) {
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * std::pow(x_c[i],p) );
    s_mu = icalc_scal_logist_mm(loc_mu,0.0,1.0,a_prime,a,s_p1[i]);
    s_srce[i] += (x_c[i] * x_c[i] * S_data_r * S_data_r) * ( s_mu * s_g[i] - delta_h * s_h[i] );
  }

  return;
}

void TDSolveSacado::calc_srce_n_opn(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  // The same as in the pdgf model
  calc_srce_n_pdgf( s_srce, tsol_loc );
  return;
}

void TDSolveSacado::calc_srce_g_model(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  calc_srce_g_base( s_srce, tsol_loc );

  DType lambda = tdCoef->lambda;
  SacadoIndVarArray& s_g = S_data_var["G"];

  for ( int i = 0; i < Nx; i++ ) {
    s_srce[i] += x_c[i] * x_c[i] * S_data_r * S_data_r * ( lambda * s_g[i] );
  }

  return;
}

void TDSolveSacado::calc_srce_n_model(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  calc_srce_n_base( s_srce, tsol_loc );

  DType mu      = tdCoef->mu;
  DType delta_n = tdCoef->delta_n;
  SacadoIndVarArray& s_g = S_data_var["G"];
  SacadoIndVarArray& s_n = S_data_var["N"];

  for ( int i = 0; i < Nx; i++ ) {
    s_srce[i] += (x_c[i] * x_c[i] * S_data_r * S_data_r) * ( mu * s_g[i] - delta_n * s_n[i] );
  }

  return;
}

void TDSolveSacado::compute_source_chem_prob(const int nt, const DType tsol_loc) {
  std::vector<std::string>::iterator it_name;
  for ( it_name = chemnames.begin(); it_name != chemnames.end(); it_name++ ) {
    calc_srce_chem_var( *it_name, nt, tsol_loc );
  }
  return;
}

void TDSolveSacado::calc_srce_chem_var(std::string chemname, const int nt, const DType tsol_loc) {
  // May need to incorporate the additional source such as injection in ChemPar
  SacadoDepVarArray& s_srce  = S_srce_chem[chemname];
  DType*             l_x_c   = varcoord[chemname];
  ChemPar&           chemPar = chempars[chemname];
  SacadoIndVarArray& s_chem  = S_data_var[chemname];
  SacadoIndVarArray& s_cell  = S_data_var[chemPar.driver];

  int num_chem = s_chem.size();
  int num_cell = s_cell.size(); // The driver may have a shorter range than the chemicals

  for ( int i = 0; i < num_chem; i++ ) {
    s_srce[i] += l_x_c[i] * l_x_c[i] * data_r[nt] * data_r[nt] * ( - chemPar.decay * s_chem[i] );
  }

  for ( int i = 0; i < num_cell; i++ ) {
    s_srce[i] += l_x_c[i] * l_x_c[i] * data_r[nt] * data_r[nt] * 
      ( ( chemPar.mm_ceil * s_cell[i] ) / ( chemPar.mm_half + s_cell[i] ) );
  }

  return;
}

DType TDSolveSacado::calc_n_bc( const DType tsol_loc ) {
  DType n_bc = tdCoef->n_bc;
  if ( tdSetup->prob_id == -1 ) {
    switch ( tdSetup->infil_u ) {
      case 2:
        DType v_par = tdSetup->init_u_par;
        n_bc = 1.0 - std::exp(-3.0*v_par*tsol_loc) * ( 1.0-n_bc );
        break;
    }
  }
  return n_bc;
}


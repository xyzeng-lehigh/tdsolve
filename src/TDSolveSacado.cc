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
  prob_name = "";
  force_nextIter = false;
  cellnames.clear();
  cellnames_u.clear();
  velnames.clear();
  chemnames.clear();
  ale_chemnames.clear();
  varsizes.clear();
  varcoord.clear();
  varflxcoord.clear();
  chempars.clear();
  varlspacks.clear();
  dimlspacks.clear();
  H_cell_offset.clear();
  cellvels.clear();
  setupVarMap.clear();
  setupICMap.clear();
  setupSrceMap.clear();
  setupVelUFlxMap.clear();
  setupVelUMap.clear();
  setupExtraVFlxMap.clear();
  setupAuxMap.clear();

  setupVarMap.insert( std::pair<int, fSetupFun>(-1, &TDSolveSacado::setup_var_model) );
  setupVarMap.insert( std::pair<int, fSetupFun>( 0, &TDSolveSacado::setup_var_pdgf) );
  setupVarMap.insert( std::pair<int, fSetupFun>( 1, &TDSolveSacado::setup_var_opn) );
  setupVarMap.insert( std::pair<int, fSetupFun>( 2, &TDSolveSacado::setup_var_oncolytic) );
  setupVarMap.insert( std::pair<int, fSetupFun>( 3, &TDSolveSacado::setup_var_potassium) );

  // Setup variables
  (this->*setupVarMap[tdSetup->prob_id])();

  setupICMap.insert( std::pair<int, fSetupFun>(-1, &TDSolveSacado::setup_ic_model) );
  setupICMap.insert( std::pair<int, fSetupFun>( 0, &TDSolveSacado::setup_ic_pdgf) );
  setupICMap.insert( std::pair<int, fSetupFun>( 1, &TDSolveSacado::setup_ic_opn) );
  setupICMap.insert( std::pair<int, fSetupFun>( 2, &TDSolveSacado::setup_ic_oncolytic) );
  setupICMap.insert( std::pair<int, fSetupFun>( 3, &TDSolveSacado::setup_ic_potassium) );

  setupSrceMap.insert( std::pair<int, fSetupFun>(-1, &TDSolveSacado::setup_srce_model) );
  setupSrceMap.insert( std::pair<int, fSetupFun>( 0, &TDSolveSacado::setup_srce_pdgf) );
  setupSrceMap.insert( std::pair<int, fSetupFun>( 1, &TDSolveSacado::setup_srce_opn) );
  setupSrceMap.insert( std::pair<int, fSetupFun>( 2, &TDSolveSacado::setup_srce_oncolytic) );
  setupSrceMap.insert( std::pair<int, fSetupFun>( 3, &TDSolveSacado::setup_srce_potassium) );

  setupVelUFlxMap.insert( std::pair<int, fSetupFun>(-1, &TDSolveSacado::add_vupdate_uflx) );
  setupVelUFlxMap.insert( std::pair<int, fSetupFun>( 0, &TDSolveSacado::add_vupdate_uflx) );
  setupVelUFlxMap.insert( std::pair<int, fSetupFun>( 1, &TDSolveSacado::add_vupdate_uflx) );

  setupVelUMap.insert( std::pair<int, fSetupTFun>(-1, &TDSolveSacado::compute_vel_u_pres) );
  setupVelUMap.insert( std::pair<int, fSetupTFun>( 0, &TDSolveSacado::compute_vel_u_chem) );
  setupVelUMap.insert( std::pair<int, fSetupTFun>( 1, &TDSolveSacado::compute_vel_u_chem) );

  setupExtraVFlxMap.insert( std::pair<int, fSetupTFun>(-1, &TDSolveSacado::compute_extra_vflx_infiltration) );
  setupExtraVFlxMap.insert( std::pair<int, fSetupTFun>( 0, &TDSolveSacado::compute_extra_vflx_infiltration) );
  setupExtraVFlxMap.insert( std::pair<int, fSetupTFun>( 1, &TDSolveSacado::compute_extra_vflx_infiltration) );

  setupAuxMap.insert( std::pair<int, fSetupTFun>(2, &TDSolveSacado::setup_aux_oncolytic) );

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

  srce_max = NULL;
}

TDSolveSacado::~TDSolveSacado(void) {
  erase_array( data_inc, 1 );
  for ( std::map<std::string,HVEC**>::iterator it = H_data_var.begin();
        it != H_data_var.end(); it++ ) {
    erase_array( it->second, Ndata );
  }
  for ( std::map<std::string,HVEC**>::iterator it = H_data_vel.begin();
        it != H_data_vel.end(); it++ ) {
    erase_array( it->second, Ndata );
  }
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
  erase_array( srce_max, 1 );
  if ( tdSetup->chemo_opt == 1 ) {
    for ( int i = 0; i < 2; i++ ) {
      delete[] x_ale_c[i];
      delete[] x_ale_f[i];
      delete[] dx_ale[i];
    }
  }
}

void TDSolveSacado::setup(void) {
  std::vector<std::string>::iterator it;

  // Setup linear solver for cell numbers if an implicit method is used
  if ( tdSetup->assemb == 1 ) {
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
  bool invalid_ic = (this->*setupICMap[tdSetup->prob_id])();

  if ( invalid_ic ) {
    std::cout << "Error: No initial condition is defined ... Exit ...\n" << std::endl;
    exit(1);
  }

  calc_total_inc_viol(0);

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
    // Copy S_data_vel["V"] to H_data_vel["V"][0]
    set_from_sacado_vel(0);
    // This will predict a S_data_rp that is used for the
    // first time step size computation.
    // Using 0, 0, only S_data_rp will be computed.
    compute_rad(0, 0);
  }

  return;
}

void TDSolveSacado::calc_dt(const int nt) {
  force_nextIter = false;
  if ( iter != 0 )
    tsol += dtn;

  if ( !tdSetup->fix_domain ) {
    if ( data_r[nt] <= tdSetup->Rdead || data_r[nt] >= tdSetup->Rterm )
      reach_term_now = true;
  }

  if ( reach_term_next )
    reach_term_now = true;

  DType dtn_trial = DTYPE_MAX;
  DType dtn_cfl1  = DTYPE_MAX; // dtn corresponding to CFL 1

  if ( tdSetup->cellmotion ) {
    // Predicting R'. Cannot use S_data_rp, as some of calc_dt
    // is called when dtn is not valid, hence also S_data_rp
    DType dr = (*H_data_vel["V"][nt])[varsizes["V"]-1];

    for ( int i = 0; i < Nx+1; i++ ) {
      DType v_max = std::fabs( x_f[i] * dr / data_r[nt] );
      for ( std::vector<std::string>::iterator it = velnames.begin();
            it != velnames.end(); it++ ) {
        v_max += std::fabs( (*H_data_vel[*it][nt])[i] / data_r[nt] );
      }
      if ( tdSetup->pos_fix != 0 ) {
        DType lip_srce = 0.0;
        if ( i > 0 )
          lip_srce = max( lip_srce, (*srce_max[0])[i-1] );
        if ( i < Nx )
          lip_srce = max( lip_srce, (*srce_max[0])[i] );
        dtn_trial = min( dtn_trial, tdSetup->cfl*(1.0/(lip_srce+v_max/dx)) );
      }
      else {
        dtn_trial = min( dtn_trial, tdSetup->cfl * dx / v_max );
      }
      dtn_trial = min( dtn_trial, tdSetup->dTmax );
      dtn_cfl1  = min( dtn_cfl1, dx / v_max );
    }
  }

  if ( iter == 0 )
    dtn = min( dtn_trial, dtn );
  else
    dtn = min( dtn_old * tdSetup->dTgrow, dtn_trial );
  dtn = min( dtn, tdSetup->dTmax );
  if ( dtn < tdSetup->dTmin ) {
    fprintf( stdout, "    --- Time step size too small, switch to minimum value.\n" );
    dtn = tdSetup->dTmin;
    force_nextIter = true;
  }
  dtn_old = dtn;

  if ( use_tfreq ) {
    if ( iter == 0 )
      reach_output_now = true;
    else
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
  fprintf( ot_file, "%10.20f\t%10.20f", tsol, data_r[nt] );
  for ( std::vector<std::string>::iterator it = velnames.begin();
        it != velnames.end(); it++ )
    fprintf( ot_file, "\t%10.20f", (*H_data_vel[*it][nt])[varsizes[*it]-1] );
  fprintf( ot_file, "\n" );
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
  if ( tdSetup->cellmotion )
    fprintf( stdout, "Iter %3d: Time %6.8f; Time step %6.8f; CFL %2.4f; Radius %6.8f; Front-velocity %6.8f.\n",
        iter, tsol, dtn, cur_cfl, data_r[nt], S_data_rp.val() );
  else 
    fprintf( stdout, "Iter %3d: Time %6.8f; Time step %6.8f.\n", iter, tsol, dtn );
  return;
}

void TDSolveSacado::alloc(void) {
  Ndata = 3;

  for ( std::vector<std::string>::iterator it = cellnames.begin();
        it != cellnames.end(); it++ ) {
    HVEC** data = NULL;
    alloc_array_full( data, varsizes[*it], Ndata );
    H_data_var.insert(std::pair<std::string,HVEC**>(*it,data));
    S_data_num += varsizes[*it];
  }

  for ( std::vector<std::string>::iterator it = chemnames.begin();
        it != chemnames.end(); it++ ) {
    HVEC** data = NULL;
    alloc_array_full( data, varsizes[*it], Ndata );
    H_data_var.insert(std::pair<std::string,HVEC**>(*it,data));
    S_data_num += varsizes[*it];
  }

  for ( std::vector<std::string>::iterator it = velnames.begin();
        it != velnames.end(); it++ ) {
    HVEC** data = NULL;
    alloc_array_full( data, varsizes[*it], Ndata );
    H_data_vel.insert(std::pair<std::string,HVEC**>(*it,data));
  }

  alloc_array_full( data_inc, Nx, 1 );

  data_r  = new DType[Ndata];

  alloc_hvec_res();

  alloc_sacado_ind();

  alloc_sacado_dep();

  alloc_functions();

  if ( tdSetup->pos_fix != 0 )
    alloc_array_full( srce_max, Nx, 1 );

  if ( tdSetup->chemo_opt == 1 ) {
    for ( int i = 0; i < 2; i++ ) {
      x_ale_f[i] = new DType[Na+1];
      x_ale_c[i] = new DType[Na];
      dx_ale[i]  = new DType[Na];
    }
  }

  return;
}

bool TDSolveSacado::update_sol(const int fr, const int to, const DType shift_dt_fr, const DType shift_dt_to) {
  // fr - solution at t_n
  // to - solution at t_{n+1}/Newton iterates
  copy_data( fr, to );

  // Setup updates schemes
  int  numIter  = 1;
  bool assembleJacobian = false;
  DType ntTol, cur_res=DTYPE_MAX;
  bool displayNewton = false;
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
      ntTol = tdSetup->ntTol;
      displayNewton = true;
      break;
    case 2: // PMC - just add for fun, not really using
      numIter = tdSetup->ntIter;
      assembleJacobian = false;
      ntTol = -1.0;
      shift_dt_space = shift_dt_fr;
      break;
    default:
      break;
  }

  // ! All fluxes are on the left-hand side, including the diffusion ones
  // ! All source terms are on the right-hand side.
  int iter = 0;
  while ( iter < numIter && cur_res > ntTol ) {
    // Step 0. Setup the Sacado variables and initialize local inf-res
    DType loc_res = 0.0;
    reset_sacado_ind_array();
    reset_sacado_dep_array();
    // Set the independent variables iterates to [to] 
    set_sacado_ind( to );

    // Step 1. Compute S_data_r and S_data_rp.
    //   assemb == 0 : S_data_r = data_r[fr], S_data_rp ~ H_data_vel["V"][fr]
    //   assemb == 1 : S_data_r = data_r[to], S_data_rp ~ H_data_vel["V"][to]
    // Note that because [fr] is copied to [to] in iter 0,
    // H_data_vel["V"][to] is always used in this function to compute S_data_rp
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
      update_cells_by_res_imp( to, loc_res );
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

    // Step 8*. Compute the fluxes for chemo on fixed range
    if ( tdSetup->chemo_opt == 1 ) {
      update_ale_grids( fr ); // For now, suppose fr = 0
      update_ale_grids( to ); // For now, suppose to = 0
      compute_ale_source_chem_prob( to, tsol+shift_dt_to*dtn );
      compute_ale_chems_dflux( to );
      // S_flux_chem_r["A"] is treated as the difference in flux here
      // hence it is cell-centered and the last entry is not used at all
      compute_ale_chems_rflux( fr, to  );
      compute_ale_chems_residual( fr, to );
      update_ale_chems_by_res_imp( to );
    }

    // Step 8. Compute the fluxes for chems. In order:
    // 8.1 - S_flux_chem_d["X"]
    //     - it is independent of R, hence no passing of nt
    compute_chems_dflux( to, tsol+shift_dt_to*dtn );
    // 8.2 - S_flux_chem_r["X"]
    compute_chems_rflux( to, tsol+shift_dt_to*dtn );
    // 8.2 - S_flux_chem["X"]
    compute_chems_flux();

    // Step 9. Compute all residuals for chems
    // The results are stored in H_data_res
    // Note that dt is multuplied to the spatial terms
    compute_chems_residual( fr, to );

    // Step X. Update the solution/iterates for chems
    // X-1. Implicit updates of chemicals
    update_chems_by_res_imp( to, loc_res );
    // X-2. Make consistent the velocity for dtn
    reset_sacado_ind_array();
    reset_sacado_dep_array();
    set_sacado_ind( to );
    set_sacado_r( to );
    update_vel( to, to, shift_dt_to );
    set_from_sacado_vel(to);
    compute_rad( to, to );
    
    cur_res = loc_res;
    iter++;
  }

  if ( displayNewton ) {
    fprintf( stdout, "    --- Newton method: %3d iterations, Linf residual %e.\n", iter, cur_res );
  }

  // Aux variables
  std::map<int,fSetupTFun>::iterator it = setupAuxMap.find( tdSetup->prob_id );
  if ( it != setupAuxMap.end() )
    (this->*setupAuxMap[tdSetup->prob_id])(tsol+shift_dt_to*dtn);

  if ( !tdSetup->auto_red_dt || force_nextIter ) // Ignore convergence of Newton solve
    return true;
  else if ( iter == numIter && cur_res > ntTol ) { // Non-convergent Newton solve
    fprintf( stdout, "    --- Newton method not converging, reduce dt by half.\n" );
    return false;
  }
  else // Successful Newton solve
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
  if ( tdSetup->cellmotion ) {
    std::map<int,fSetupTFun>::iterator it = setupExtraVFlxMap.find( tdSetup->prob_id );
    if ( it != setupExtraVFlxMap.end() )
      (this->*setupExtraVFlxMap[tdSetup->prob_id])(tsol+shift_dt*dtn);
  }
  // Compute the right-hand sides (source terms)
  compute_source_cell_prob( tsol+shift_dt*dtn );
  if ( tdSetup->cellmotion ) {
    // Compute the velocity (integral form)
    compute_vel_v();
  }
  return;
}

DType TDSolveSacado::get_gperc(const int nt) {
  return get_cell_perc("G",nt);
}

DType TDSolveSacado::get_nperc(const int nt) {
  return get_cell_perc("N",nt);
}

DType TDSolveSacado::get_argmax(const std::string name, const int nt) {
  std::vector<std::string>::iterator it = std::find( cellnames.begin(), cellnames.end(), name );
  if ( it == cellnames.end() ) {
    std::cout << " Warning: the cell variable name " << name << " is not registered." << std::endl;
    return 0.0;
  }
  int num = varsizes[name];
  DType max_val = -DTYPE_MAX;
  int i_max = -1;
  DType* l_x_c = varcoord[name];
  for ( int i = 0; i < num; i++ ) {
    if ( max_val < (*H_data_var[name][nt])[i] ) {
      max_val = (*H_data_var[name][nt])[i];
      i_max = i;
    }
  }
  return data_r[nt] * l_x_c[i_max];
}

DType TDSolveSacado::get_total(const std::string name, const int nt) {
  std::vector<std::string>::iterator it = std::find( cellnames.begin(), cellnames.end(), name );
  if ( it == cellnames.end() ) {
    std::cout << " Warning: the cell variable name " << name << " is not registered." << std::endl;
    return 0.0;
  }
  int num = varsizes[name];
  DType* l_x_c = varcoord[name];
  DType tot_var = 0.0;
  DType dr = data_r[nt] * dx;
  for ( int i = 0; i < num; i++ ) {
    DType r_loc = data_r[nt] * l_x_c[i];
    tot_var += dr * tdSetup->sweight * std::pow(r_loc,tdSetup->dim-1) * (*H_data_var[name][nt])[i];
  }
  return tot_var;
}

void TDSolveSacado::add_cellular_var(std::vector<std::string>& names, std::string name, const int size, DType* x_var, DType* x_flx) {
  names.push_back( name );
  varsizes.insert(std::pair<std::string,int>(name,size));
  varcoord.insert(std::pair<std::string,DType*>(name,x_var));
  varflxcoord.insert(std::pair<std::string,DType*>(name,x_flx));
  return;
}

void TDSolveSacado::add_nodal_var(std::vector<std::string>& names, std::string name, const int size, DType* x_var) {
  names.push_back( name );
  varsizes.insert(std::pair<std::string,int>(name,size));
  varcoord.insert(std::pair<std::string,DType*>(name,x_var));
  return;
}

bool TDSolveSacado::setup_var_model(void) {
  prob_name = "model";
  add_cellular_var( cellnames, "G", Nx, x_c, x_f );
  add_cellular_var( cellnames, "N", Nx, x_c, x_f );
  add_nodal_var( velnames, "V", Nx+1, x_f );
  add_nodal_var( velnames, "U", Nx+1, x_f );
  cellnames_u.push_back( "N" );
  return true;
}

bool TDSolveSacado::setup_var_pdgf(void) {
  prob_name = "pdgf";
  add_cellular_var( cellnames, "G", Nx, x_c, x_f );
  add_cellular_var( cellnames, "N", Nx, x_c, x_f );
  add_cellular_var( cellnames, "H", Nx, x_c, x_f );
  add_nodal_var( velnames, "V", Nx+1, x_f );
  add_nodal_var( velnames, "U", Nx+1, x_f );
  cellnames_u.push_back( "N" );

  add_cellular_var( chemnames, "A", Na, x_a_c, x_a_f );
  if ( tdSetup->chemo_opt == 1 )
    ale_chemnames.push_back( "A" );

  setup_chem_par( "A", tdCoef->gamma_a, tdCoef->D, "G", tdCoef->alpha_a, tdCoef->beta_a/tdCoef->theta, true );
  // XYZ - testing for a mild infiltration velocity only
  //setup_chem_par( "A", tdCoef->gamma_a, tdCoef->D, "G", tdCoef->alpha_a, tdCoef->beta_a/tdCoef->theta, true, 1 );

  return true;
}

bool TDSolveSacado::setup_var_opn(void) {
  prob_name = "opn";
  add_cellular_var( cellnames, "G", Nx, x_c, x_f );
  add_cellular_var( cellnames, "N", Nx, x_c, x_f );
  add_cellular_var( cellnames, "H", Nx, x_c, x_f );
  add_nodal_var( velnames, "V", Nx+1, x_f );
  add_nodal_var( velnames, "U", Nx+1, x_f );
  cellnames_u.push_back( "N" );

  add_cellular_var( chemnames, "A", Na, x_a_c, x_a_f );
  if ( tdSetup->chemo_opt == 1 )
    ale_chemnames.push_back( "A" );
  add_cellular_var( chemnames, "P1", Nx, x_c, x_f );
  add_cellular_var( chemnames, "P2", Nx, x_c, x_f );

  setup_chem_par( "A", tdCoef->gamma_a, tdCoef->D, "G", tdCoef->alpha_a, tdCoef->beta_a/tdCoef->theta, true );
  //setup_chem_par( "P1", tdCoef->gamma_1, tdCoef->D_1, "N", tdCoef->alpha_1, tdCoef->beta_1/tdCoef->theta, false );
  setup_chem_par( "P1", tdCoef->gamma_1, tdCoef->D_1, "N", tdCoef->alpha_1, tdCoef->beta_1/tdCoef->theta, tdCoef->use_dbc_P1 );
  setup_chem_par( "P2", tdCoef->gamma_2, tdCoef->D_2, "G", tdCoef->alpha_2, tdCoef->beta_2/tdCoef->theta, false );
  if ( tdCoef->treatP1BC )
    chempars["P1"].F_chemBC = &TDSolveSacado::calc_treatment_P1;

  return true;
}

bool TDSolveSacado::setup_var_oncolytic(void) {
  prob_name = "oncolytic";
  add_cellular_var( cellnames, "B",  Nx, x_c, x_f );
  add_cellular_var( chemnames, "Vi", Nx, x_c, x_f );
  add_cellular_var( cellnames, "H",  Nx, x_c, x_f );
  add_cellular_var( cellnames, "N",  Nx, x_c, x_f );
  add_cellular_var( cellnames, "G",  Nx, x_c, x_f );
  add_cellular_var( cellnames, "Mv", Nx, x_c, x_f );
  add_cellular_var( cellnames, "Mb", Nx, x_c, x_f );
  add_cellular_var( chemnames, "A",  Nx, x_c, x_f );
  add_nodal_var( velnames, "V", Nx+1, x_f );

  setup_chem_par_gen( "Vi", tdCoef->probVarParVals[prob_name]["Vi"]["grow_neg"], 
      tdCoef->probVarParVals[prob_name]["Vi"]["diffusion"], 
      false, &TDSolveSacado::calc_chem_srce_vi_oncolytic );
  setup_chem_par_gen( "A",  0.0, tdCoef->probVarParVals[prob_name]["A"]["diffusion"],
      false, &TDSolveSacado::calc_chem_srce_a_oncolytic );
  return true;
}

bool TDSolveSacado::setup_var_potassium(void) {
  prob_name = "potassium";
  assert( tdSetup->chemo_opt == 0 );
  tdSetup->fix_domain = true;  // data_r will not change
  tdSetup->cellmotion = false; // no velnames
  add_cellular_var( cellnames, "B1", Nx, x_c, x_f );
  add_cellular_var( cellnames, "B2", Nx, x_c, x_f );
  add_cellular_var( cellnames, "B3", Nx, x_c, x_f );
  add_cellular_var( chemnames, "P", Nx, x_c, x_f );
  add_cellular_var( chemnames, "V", Nx, x_c, x_f );
  
  if ( tdCoef->probVarParVals[prob_name]["P"]["slope"] == 0.0 ) {
    std::cout << " Warning: Setting the P.slope to zero will turn off the blocking of potassium entirely instead of using the smoothed Heaviside function." << std::endl;
  }

  setup_chem_par_gen( "P", tdCoef->probVarParVals[prob_name]["P"]["grow_neg"],
      tdCoef->probVarParVals[prob_name]["P"]["diffusion"], 
      false, &TDSolveSacado::calc_chem_srce_p_potassium ); // false -> do not use Dirichlet bc
  setup_chem_par_gen( "V", tdCoef->probVarParVals[prob_name]["V"]["grow_neg"], 
      tdCoef->probVarParVals[prob_name]["V"]["diffusion"], 
      false, &TDSolveSacado::calc_chem_srce_v_potassium ); // false -> do not use Dirichlet bc
  return true;
}

bool TDSolveSacado::setup_ic_model(void) {
  bool invalid_ic = false;
  switch ( tdSetup->ic ) {
    case 0:
      setic_unif_cell("G",1.0-tdCoef->n_bc);
      setic_unif_cell("N",tdCoef->n_bc);
      break;
    case 1:
      setic_sine_cell("G",0.5,0.5,tdCoef->ic_sine_freq);
      setic_sine_cell("N",0.5,-0.5,tdCoef->ic_sine_freq);
      break;
    case 2:
      setic_step_cell("G",0.5,tdCoef->n_bc,1.0-tdCoef->n_bc);
      setic_step_cell("N",0.5,1.0-tdCoef->n_bc,tdCoef->n_bc);
      break;
    default:
      invalid_ic = true;
      break;
  }
  return invalid_ic;
}

bool TDSolveSacado::setup_ic_pdgf(void) {
  bool invalid_ic = false;
  switch ( tdSetup->ic ) {
    case 0: 
      setic_unif_cell("G",0.84);
      setic_unif_cell("H",0.16-tdCoef->n_bc);
      setic_unif_cell("N",tdCoef->n_bc);
      setic_expd_chem(1.0e+3,1.0);
      break;
    case 10:
      setic_unif_cell("G",0.84);
      setic_unif_cell("H",0.16-tdCoef->n_bc);
      setic_unif_cell("N",tdCoef->n_bc);
      (*(H_data_var["A"])[0]) = 0.0;
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
  return invalid_ic;
}

bool TDSolveSacado::setup_ic_opn(void) {
  bool invalid_ic = false;
  switch ( tdSetup->ic ) {
    case 0: 
      setic_unif_cell("G",0.84);
      setic_unif_cell("H",0.16-tdCoef->n_bc);
      setic_unif_cell("N",tdCoef->n_bc);
      setic_expd_chem(1.0e+3,1.0);
      break;
    case 10:
      setic_unif_cell("G",0.84);
      setic_unif_cell("H",0.16-tdCoef->n_bc);
      setic_unif_cell("N",tdCoef->n_bc);
      (*(H_data_var["A"])[0]) = 0.0;
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

  (*(H_data_var["P1"])[0]) = 0.0;
  (*(H_data_var["P2"])[0]) = 0.0;

  return invalid_ic;
}

bool TDSolveSacado::setup_ic_oncolytic(void) {
  in_oncolytic_stage[0] = true;  // The first stage 
  in_oncolytic_stage[1] = false; // The second stage 
  in_oncolytic_stage[2] = false; // The second stage 

  bool invalid_ic = false;
  DType theta = tdCoef->theta;

  switch ( tdSetup->ic ) {
    case 0:
      for ( std::map<std::string,HVEC**>::iterator it = H_data_var.begin();
            it != H_data_var.end(); it++ )
        *(it->second)[0] = tdCoef->probVarParVals[prob_name][it->first]["init"];
      break;
    default:
      invalid_ic = true;
  }

  return invalid_ic;
}

bool TDSolveSacado::setup_ic_potassium(void) {
  bool invalid_ic = false;
  setic_unif_cell("B1",1.0);
  setic_unif_cell("B2",0.0);
  setic_unif_cell("B3",0.0);
  setic_unif_cell("P", 0.0);
  setic_unif_cell("V", 0.0);
  (*(H_data_var["V"])[0])[0] = tdCoef->probVarParVals[prob_name]["V"]["init_num"] / tdSetup->vweight / std::pow(tdCoef->R_0*dx,tdSetup->dim);
  return invalid_ic;
}

bool TDSolveSacado::setup_srce_model(void) {
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("G", &TDSolveSacado::calc_srce_g_model) );
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("N", &TDSolveSacado::calc_srce_n_model) );
  return true;
}

bool TDSolveSacado::setup_srce_pdgf(void) {
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("G", &TDSolveSacado::calc_srce_g_pdgf) );
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("H", &TDSolveSacado::calc_srce_h_pdgf) );
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("N", &TDSolveSacado::calc_srce_n_pdgf) );
  return true;
}

bool TDSolveSacado::setup_srce_opn(void) {
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("G", &TDSolveSacado::calc_srce_g_opn) );
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("H", &TDSolveSacado::calc_srce_h_opn) );
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("N", &TDSolveSacado::calc_srce_n_opn) );
  return true;
}

bool TDSolveSacado::setup_srce_oncolytic(void) {
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("B",  &TDSolveSacado::calc_srce_b_oncolytic) );
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("H",  &TDSolveSacado::calc_srce_h_oncolytic) );
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("N",  &TDSolveSacado::calc_srce_n_oncolytic) );
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("G",  &TDSolveSacado::calc_srce_g_oncolytic) );
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("Mv", &TDSolveSacado::calc_srce_mv_oncolytic) );
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("Mb", &TDSolveSacado::calc_srce_mb_oncolytic) );
  return true;
}

bool TDSolveSacado::setup_srce_potassium(void) {
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("B1",  &TDSolveSacado::calc_srce_b1_potassium) );
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("B2",  &TDSolveSacado::calc_srce_b2_potassium) );
  F_calcSrceMap.insert( std::pair<std::string, fCalcSrce>("B3",  &TDSolveSacado::calc_srce_b3_potassium) );
  return true;
}

bool TDSolveSacado::add_vupdate_uflx(void) {
  SacadoDepVarArray& s_v = S_data_vel["V"];
  SacadoDepVarArray& s_flux_u = S_flux_cell_v["U"]["N"];
  for ( int i = 1; i < Nx+1; i++ ) {
    s_v[i] = s_v[i] - s_flux_u[i];
  }
  return true;
}

bool TDSolveSacado::compute_extra_vflx_infiltration(const DType tsol_loc) {
  // Compute the velocity u using the gradient
  compute_vel_u( tsol_loc );
  // Compute the flux associated with U
  compute_cells_varflux( "U", "N", tsol_loc, true );
  return true;
}

bool TDSolveSacado::setup_aux_oncolytic(const DType tsol_loc) {
  return true;
}

void TDSolveSacado::add_vel_to_cellvels(std::string cname, std::string vname) {
  std::map<std::string,std::set<std::string>>::iterator it = cellvels.find( cname );
  if ( it == cellvels.end() ) {
    std::set<std::string> vnames;
    vnames.clear();
    cellvels.insert( std::pair<std::string,std::set<std::string>>( cname, vnames ) );
  }
  it = cellvels.find( cname );
  std::set<std::string>::iterator vit = it->second.find( vname );
  if ( vit == it->second.end() )
    it->second.insert( vname );
  return;
}

void TDSolveSacado::setup_chem_par( std::string chemname, 
    const DType decay, const DType diffusion, std::string cellname, 
    const DType mm_ceil, const DType mm_half, const bool dbc, const int test_opt ) {
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
  chemPar.F_chemBC     = NULL;
  chemPar.test_opt     = test_opt;

  chemPar.srceFun      = &TDSolveSacado::calc_srce_chem_var;

  chempars.insert( std::pair<std::string, ChemPar>( chemname, chemPar ) );
  return;
}

bool TDSolveSacado::setup_chem_par_gen(std::string chemname, const DType decay,
    const DType diffusion, const bool dbc, fCalcTSrce srceFun ) {
  std::map<std::string, ChemPar>::iterator it = chempars.find( chemname );
  if ( it != chempars.end() ) {
    std::string error_message = "TDSolveSacado::setup_chem_par - chemoattractant/opn parameter for " + chemname + " has already been setup.\n";
    throw std::runtime_error( error_message );
  }
  ChemPar chemPar;
  chemPar.decay        = decay;
  chemPar.diffusion    = diffusion;
  chemPar.dirichlet_bc = dbc;
  chemPar.F_chemBC     = NULL;
  chemPar.srceFun      = srceFun;

  chempars.insert( std::pair<std::string, ChemPar>( chemname, chemPar ) );
  return true;
}

void TDSolveSacado::setic_unif_cell(std::string name, DType val, int nt) {
  *H_data_var[name][nt] = val;
  return;
}

void TDSolveSacado::setic_sine_cell(std::string name, DType c0, DType c1, DType cs, int nt) {
  DType* l_x_c = varcoord[name];
  int    num   = varsizes[name];
  for ( int i = 0; i < num; i++ ) {
    (*H_data_var[name][nt])[i] = c0 + c1 * std::sin( cs * PI * l_x_c[i] );
  }
  return;
}

void TDSolveSacado::setic_step_cell(std::string name, DType x_sep, DType val_l, DType val_r, int nt) {
  DType* l_x_c = varcoord[name];
  int    num   = varsizes[name];
  for ( int i = 0; i < num; i++ ) {
    if ( l_x_c[i] < x_sep )
      (*H_data_var[name][nt])[i] = val_l;
    else
      (*H_data_var[name][nt])[i] = val_r;
  }
  return;
}

void TDSolveSacado::setic_expd_chem(DType multip, DType decay, int nt) {
  HVEC** data_chem = H_data_var["A"];
  if ( tdSetup->chemo_opt == 1 ) {
    update_ale_grids( nt );
    for ( int i = 0; i < Na; i++ ) {
      (*data_chem[nt])[i] = multip * std::exp( - decay * x_ale_c[nt][i] * x_ale_c[nt][i] );
    }
  }
  else {
    for ( int i = 0; i < Na; i++ ) {
      DType tvar = x_a_c[i] * data_r[nt];
      (*data_chem[nt])[i] = multip * std::exp( - decay * tvar * tvar );
    }
  }
  return;
}

void TDSolveSacado::setic_polyexp_chem(DType poly_0, DType poly_2, int nt) {
  HVEC** data_chem = H_data_var["A"];
  if ( tdSetup->chemo_opt == 1 ) {
    for ( int i = 0; i < Nx; i++ ) {
      (*data_chem[nt])[i] = poly_0 - poly_2 * x_ale_c[nt][i] * x_ale_c[nt][i] / data_r[nt] / data_r[nt];
    }
    DType exp_coef = poly_0 - poly_2;
    DType exp_rate = 2*poly_2/exp_coef;
    for ( int i = Nx; i < Na; i++ ) {
      (*data_chem[nt])[i] = exp_coef * exp(-exp_rate*(x_ale_c[nt][i]/data_r[nt]-1.0));
    }
  }
  else {
    for ( int i = 0; i < Nx; i++ ) {
      (*data_chem[nt])[i] = poly_0 - poly_2 * x_a_c[i] * x_a_c[i];
    }
    DType exp_coef = poly_0 - poly_2;
    DType exp_rate = 2*poly_2/exp_coef;
    for ( int i = Nx; i < Na; i++ ) {
      (*data_chem[nt])[i] = exp_coef * exp(-exp_rate*(x_a_c[i]-1.0));
    }
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

    add_vel_to_cellvels( *it_name, "V" );
  }
  S_flux_cell_v.insert( std::pair<std::string, SacadoDepVarMap>("V", S_array_v) );
  // Flux associated with u
  SacadoDepVarMap S_array_u; S_array_u.clear();
  for ( it_name = cellnames_u.begin(); it_name != cellnames_u.end(); it_name++ ) {
    SacadoDepVarArray S_array;
    alloc_sacado_dep_array(S_array, varsizes[*it_name]+1);
    S_array_u.insert( std::pair<std::string, SacadoDepVarArray>(*it_name, S_array) );

    add_vel_to_cellvels( *it_name, "U" );
  }
  S_flux_cell_v.insert( std::pair<std::string, SacadoDepVarMap>("U", S_array_u) );
  if ( tdSetup->cell_rflux == 3 ) {
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
  }
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
  (this->*setupSrceMap[tdSetup->prob_id])();
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
  HVEC** data_vel = H_data_vel["V"];
  int num = s_v.size();
  for ( int i = 0; i < num; i++ ) {
    s_v[i] = (*data_vel[nt])[i];
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

void TDSolveSacado::set_from_sacado_vel(const int nt) {
  if ( tdSetup->cellmotion ) {
    for ( std::vector<std::string>::iterator it = velnames.begin();
          it != velnames.end(); it++ )
      set_from_sacado_vel_ind( *it, nt );
  }
  return;
}

void TDSolveSacado::set_from_sacado_vel_ind( std::string name, const int nt) {
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
  if ( tdSetup->cell_rflux == 3 ) {
    reset_sacado_dep_array( S_rec_ppm_cell, 0.0 );
    for ( int i = 0; i < 2; i++ ) {
      reset_sacado_dep_array( S_lim_ppm_cell[i], 1.0 );
      reset_sacado_dep_array( S_lim_ppm_sync[i], 1.0 );
    }
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

void TDSolveSacado::reset_sacado_dep_array(SacadoDepVarArray& s_var, DType val, int num) {
  int num_max = s_var.size();
  int num_loc = std::min(num,num_max);
  for ( int i = 0; i < num_loc; i++ ) {
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
  std::map<int,fSetupTFun>::iterator it = setupVelUMap.find( tdSetup->prob_id );
  if ( it != setupVelUMap.end() )
    (this->*setupVelUMap[tdSetup->prob_id])(tsol_loc);
  return;
}

bool TDSolveSacado::compute_vel_u_pres(const DType tsol_loc) {
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
    case 5: // Exponential growth
      {
        DType v_par = tdSetup->init_u_par;
        s_u[0] = 0.0;
        for ( int i = 0; i < num; i++ ) {
          //s_u[i] = - v_par * l_x_f[i] * S_data_r / tdCoef->n_bc;
          s_u[i] = - 0.5 * v_par * tsol_loc * tsol_loc * l_x_f[i] * S_data_r / tdCoef->n_bc;
        }
      }
      break;
  }

  return true;
}

bool TDSolveSacado::compute_vel_u_chem(const DType tsol_loc) {
  SacadoDepVarArray& s_u = S_data_vel["U"];
  SacadoIndVarArray& s_a = S_data_var["A"];

  DType alpha = tdCoef->alpha;
  SacadoDepVar aor  = alpha/S_data_r;

  s_u[0] = 0.0;
  for ( int i = 1; i <= Nx; i++ )
    s_u[i] = aor * ( (s_a[i] - s_a[i-1]) / dx );

  return true;
}

void TDSolveSacado::compute_vel_v(void) {
  SacadoDepVarArray& s_v = S_data_vel["V"];
  
  s_v[0] = 0.0;

  for ( int i = 0; i < Nx; i++ ) {
    s_v[i+1] = s_v[i];
    for ( SacadoDepVarMap::iterator it = S_srce_cell.begin();
          it != S_srce_cell.end(); it++ ) {
      s_v[i+1] += dx * (it->second)[i];
    }
  }

  std::map<int,fSetupFun>::iterator it = setupVelUFlxMap.find( tdSetup->prob_id );
  if ( it != setupVelUFlxMap.end() )
    (this->*setupVelUFlxMap[tdSetup->prob_id])();

  for ( int i = 1; i < Nx+1; i++ ) {
    s_v[i] = s_v[i]/x_f[i]/x_f[i]/S_data_r;
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
  if ( tdSetup->cellmotion ) {
    for ( std::vector<std::string>::iterator it = cellnames.begin(); 
          it != cellnames.end(); it++ ) {
      // The third argument is not used at all
      compute_cells_varflux( "V", *it, 0.0, false );
    }
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
  if ( std::fabs(s_dval_c)<TOL_ZERO && std::fabs(s_dval_r)<TOL_ZERO )
    s_lim = 1.0;
  else if ( s_dval_c * s_dval_r <= 0.0 || 3.0*std::fabs(s_dval_r) <= 8.0*std::fabs(s_dval_c) )
    s_lim = 0.0;
  else
    s_lim = std::min(1.0, 5.0/(std::fabs(s_dval_r)/std::fabs(s_dval_c)-1.0));
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
  if ( tdSetup->fix_domain ) {
    for ( std::vector<string>::iterator it = cellnames.begin();
          it != cellnames.end(); it++ ) {
      reset_sacado_dep_array( S_flux_cell_r[*it] );
    }
    return;
  }

  if ( tdSetup->cell_rflux == 0 ) { // upwind
    for ( std::vector<string>::iterator it = cellnames.begin();
          it != cellnames.end(); it++ ) {
      calc_rflux_cell_var_upw( *it );
    }
  }
  else if ( tdSetup->cell_rflux == 3 ) { // synced ppm
    compute_cells_rflux_ppm();
  }
  else { // weighted upwind
    for ( std::vector<string>::iterator it = cellnames.begin();
          it != cellnames.end(); it++ ) {
      calc_rflux_cell_var_wupw( *it );
    }
  }
  return;
}

void TDSolveSacado::calc_rflux_cell_var_upw(std::string cellname) {
  SacadoIndVarArray& s_var = S_data_var[cellname];
  SacadoDepVarArray& s_flx = S_flux_cell_r[cellname];
  reset_sacado_dep_array(s_flx);
  DType* l_x_c = varcoord[cellname];
  DType* l_x_f = varflxcoord[cellname];
  int num = s_var.size();

  // - left boundary
  {
    s_flx[0] = 0.0;
  }
  // - interior points
  for ( int i = 1; i < num; i++ ) {
    SacadoDepVar loc_vel_r = - l_x_f[i] * S_data_rp * S_data_r;
    int offset = ( loc_vel_r >= 0.0 ) ? -1 : 0;
    int idx = i + offset;
    s_flx[i] = loc_vel_r * l_x_c[idx] * l_x_c[idx] * s_var[idx];
  }
  // - right boundary
  {
    s_flx[num] = 0.0;
  }
  return;
}

void TDSolveSacado::calc_rflux_cell_var_wupw(std::string cellname) {
  SacadoIndVarArray& s_var = S_data_var[cellname];
  SacadoDepVarArray& s_flx = S_flux_cell_r[cellname];
  reset_sacado_dep_array(s_flx);
  DType* l_x_f = varflxcoord[cellname];
  int num = s_var.size();

  SacadoDepVar loc_vel_r = - S_data_rp * S_data_r;
  int offset = ( loc_vel_r >= 0.0 ) ? -1 : 0;

  // - left boundary
  {
    s_flx[0] = 0.0;
  }
  // - interior points
  for ( int i = 1; i < num; i++ ) {
    int idx = i + offset;
    s_flx[i] = loc_vel_r * calc_wupw( l_x_f[i], dx ) * s_var[idx];
  }
  // - right boundary
  {
    s_flx[num] = 0.0;
  }
  return;
}

void TDSolveSacado::compute_cells_rflux_ppm(void) {
  assert( tdSetup->dim == 3 );
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

  if ( tdSetup->pos_fix == 2 ) {
    adjust_pos_ppm_lim( s_ppm_alpha );
  }

  // Compute the synced limiter
  calc_ppm_lim_sync( S_lim_ppm_sync, S_lim_ppm_cell, set_zero, s_ppm_alpha );

  // Turn off the PPM limiter in the first a few intervals!
  if ( tdSetup->ppm_reset_num < 0 ) {
    reset_sacado_dep_array( S_lim_ppm_sync[0], 1.0 );
    reset_sacado_dep_array( S_lim_ppm_sync[1], 1.0 );
  }
  if ( tdSetup->ppm_reset_num > 0 ) {
    reset_sacado_dep_array( S_lim_ppm_sync[0], 1.0, tdSetup->ppm_reset_num );
    reset_sacado_dep_array( S_lim_ppm_sync[1], 1.0, tdSetup->ppm_reset_num );
  }

  // Compute the fluxes
  for ( std::vector<std::string>::iterator it = cellnames.begin();
        it != cellnames.end(); it++ ) {
    SacadoIndVarArray& s_var = S_data_var[*it];
    SacadoDepVarArray& s_rec = S_rec_ppm_cell[*it];
    SacadoDepVarArray& s_flx = S_flux_cell_r[*it];
    calc_ppm_flux_var( s_var, s_rec, s_flx, S_lim_ppm_sync, varcoord[*it] );
  }
}

void TDSolveSacado::compute_cells_flux(void) {
  std::vector<std::string>::iterator it;
  // Flux associated with R'
  for ( it = cellnames.begin(); it != cellnames.end(); it++ ) {
    sum_in_sacado_dep_array( S_flux_cell[*it], 1.0, S_flux_cell_r[*it] );
  }
  if ( tdSetup->cellmotion ) {
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
        s_alpha1 = 2.0*std::fabs(s_diff_l)/std::fabs(s_diff_r);
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
        s_alpha[0][i] = std::min( s_alpha[0][i], s_alpha1 );
        s_alpha[1][i] = std::max( s_alpha[1][i], s_alpha2 );
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
      s_sync[1][0] = std::min( s_sync[1][0], it->second[0] );
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

void TDSolveSacado::compute_chems_dflux(const int nt, const DType tsol_loc) {
  for ( std::vector<std::string>::iterator it = chemnames.begin();
        it != chemnames.end(); it++ ) {
    if ( check_chem_skip( *it ) )
      continue;
    calc_dflux_chem_var( *it, nt, tsol_loc );
  }
  return;
}

void TDSolveSacado::calc_dflux_chem_var(std::string chemname, const int nt, const DType tsol_loc) {
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
    s_flx[i] = - dcoef * std::pow( l_x_f[i] * data_r[nt], tdSetup->dim-1 ) / data_r[nt] / data_r[nt] * (s_var[i] - s_var[i-1]) / dx;
    //s_flx[i] = - dcoef * l_x_f[i] * l_x_f[i] * (s_var[i] - s_var[i-1]) / dx;
  }
  // - right boundary
  {
    DType chem_bc = 0.0;
    if ( chempars[chemname].F_chemBC != NULL )
      chem_bc = (this->*(chempars[chemname].F_chemBC))(tsol_loc);
    if ( dbc ) { // Dirichlet b.c. 
      s_flx[num] = - dcoef * std::pow( l_x_f[num] * data_r[nt], tdSetup->dim-1 ) / data_r[nt] / data_r[nt] * (chem_bc - s_var[num-1]) / dx;
      //s_flx[num] = - dcoef * l_x_f[num] * l_x_f[num] * (chem_bc - s_var[num-1]) / dx;
    }
    else { // Neumann b.c.
      s_flx[num] = - dcoef * std::pow( l_x_f[num] * data_r[nt], tdSetup->dim-1 ) / data_r[nt] * chem_bc;
      //s_flx[num] = - dcoef * l_x_f[num] * l_x_f[num] * data_r[nt] * chem_bc;
    }
  }

  return;
}

void TDSolveSacado::compute_chems_rflux(const int nt, const DType tsol_loc) {
  if ( tdSetup->fix_domain ) {
    for ( std::vector<std::string>::iterator it = chemnames.begin();
          it != chemnames.end(); it++ ) {
      reset_sacado_dep_array( S_flux_chem_r[*it] );
    }
  }

  for ( std::vector<std::string>::iterator it = chemnames.begin();
        it != chemnames.end(); it++ ) {
    if ( check_chem_skip( *it ) )
      continue;
    if ( tdSetup->chem_rflux == 0 )
      calc_rflux_chem_var_upw( *it, nt );
    else if ( tdSetup->chem_rflux == 3 )
      calc_rflux_chem_var_ppm( *it, nt );
    else // Including == 1
      calc_rflux_chem_var_wupw( *it, nt, tsol_loc );
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
    s_flx[i] = loc_vel_r * std::pow( l_x_c[idx], tdSetup->dim-1 ) * s_var[idx];
    //s_flx[i] = loc_vel_r * l_x_c[idx] * l_x_c[idx] * s_var[idx];
  }
  // - right boundary
  {
    SacadoDepVar loc_vel_r = - l_x_f[num] * S_data_rp * data_r[nt];
    if ( dbc ) { // Homogeneous D. b.c.
      if ( loc_vel_r >= 0.0 ) {
        s_flx[num] = loc_vel_r * std::pow( l_x_c[num-1], tdSetup->dim-1 ) * s_var[num-1];
        //s_flx[num] = loc_vel_r * l_x_c[num-1] * l_x_c[num-1] * s_var[num-1];
      }
      else
        s_flx[num] = 0.0;
    }
    else { // Homogeneous N. b.c.
      DType loc_x = ( loc_vel_r >= 0.0 ) ? l_x_c[num-1] : 2.0*l_x_c[num-1]-l_x_c[num-2];
      s_flx[num] = loc_vel_r * std::pow( loc_x, tdSetup->dim-1 ) * s_var[num-1];
      //s_flx[num] = loc_vel_r * loc_x * loc_x * s_var[num-1];
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

void TDSolveSacado::calc_rflux_chem_var_wupw(std::string chemname, const int nt, const DType tsol_loc) {
  SacadoIndVarArray& s_var = S_data_var[chemname];
  SacadoDepVarArray& s_flx = S_flux_chem_r[chemname];
  reset_sacado_dep_array(s_flx);
  DType* l_x_f = varflxcoord[chemname];
  int  num = s_var.size();
  bool dbc = chempars[chemname].dirichlet_bc;

  SacadoDepVar loc_vel_r = - S_data_rp * data_r[nt];
  int offset = ( loc_vel_r >= 0.0 ) ? -1 : 0;

  // - left boundary
  {
    s_flx[0] = 0.0;
  }
  // - interior points
  for ( int i = 1; i < num; i++ ) {
    int idx = i + offset;
    s_flx[i] = loc_vel_r * calc_wupw( l_x_f[i], dx ) * s_var[idx];
  }
  // - right boundary
  {
    DType chem_bc = 0.0;
    if ( chempars[chemname].F_chemBC != NULL )
      chem_bc = (this->*(chempars[chemname].F_chemBC))(tsol_loc);
    if ( dbc ) { // Dirichlet b.c.
      if ( loc_vel_r >= 0.0 )
        s_flx[num] = loc_vel_r * calc_wupw( l_x_f[num], dx ) * s_var[num-1];
      else
        s_flx[num] = loc_vel_r * calc_wupw( l_x_f[num], dx ) * chem_bc;
    }
    else { // Neumann b.c.
      if ( loc_vel_r >= 0.0 )
        s_flx[num] = loc_vel_r * calc_wupw( l_x_f[num], dx ) * s_var[num-1];
      else // XYZ - pay attention whether to use 0.5 here
        s_flx[num] = loc_vel_r * calc_wupw( l_x_f[num], dx ) * (s_var[num-1] + dx * data_r[nt] * chem_bc);
    }
  }
  return;
}

void TDSolveSacado::compute_chems_flux(void) {
  std::vector<std::string>::iterator it;
  // Diffusive flux
  for ( it = chemnames.begin(); it != chemnames.end(); it++ ) {
    if ( check_chem_skip( *it ) )
      continue;
    sum_in_sacado_dep_array( S_flux_chem[*it], 1.0, S_flux_chem_d[*it] );
  }
  // Flux associated with R'
  for ( it = chemnames.begin(); it != chemnames.end(); it++ ) {
    if ( check_chem_skip( *it ) )
      continue;
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
  //DType r_fr2 = r_fr * r_fr;
  //DType r_to2 = r_to * r_to;
  for ( int i = 0; i < num; i++ ) {
    h_res[i] += std::pow( l_x_c[i] * r_to, tdSetup->dim-1 ) * s_var[i].val() -
                std::pow( l_x_c[i] * r_fr, tdSetup->dim-1 ) * h_var[i];
    //h_res[i] += l_x_c[i] * l_x_c[i] * ( r_to2 * s_var[i].val() - r_fr2 * h_var[i] );
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
    if ( check_chem_skip( *it ) )
      continue;
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
  //DType h_r2 = h_r * h_r;
  for ( int i = 0; i < num; i++ ) {
    h_var[i] -= h_res[i]/std::pow( l_x_c[i]*h_r, tdSetup->dim-1 );
    //h_var[i] -= h_res[i]/h_r2/l_x_c[i]/l_x_c[i];
  }
  return;
}

void TDSolveSacado::update_cells_by_res_imp(const int nt, DType& loc_res) {
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
    update_hsol_by_linsol( h_var, l_offset, num, cell_lspack, loc_res );
  }

  return;
}

void TDSolveSacado::update_cell_linsys_by_res_var(HVEC& h_res, DType h_r, DType* l_x_c, SacadoDepVarArray& s_flx, SacadoDepVarArray& s_src, const int num, const int s_offset, const int l_offset, LSPack& lsPack) {
  VComp comp(1,0.0);
  Entry elem(1,0.0);

  // RHS
  //DType s_r_to2 = h_r * h_r;
  for ( int i = 0; i < num; i++ ) {
    comp = h_res[i]/std::pow(l_x_c[i]*h_r, tdSetup->dim-1);
    //comp = h_res[i]/l_x_c[i]/l_x_c[i]/s_r_to2;
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
      elem -= dtn * s_src[i].dx( s_offset_cell+i )/std::pow(l_x_c[i]*h_r, tdSetup->dim-1);
      //elem -= dtn * s_src[i].dx( s_offset_cell+i )/l_x_c[i]/l_x_c[i]/s_r_to2;
      lsPack.mat->set(l_offset+i,l_offset_cell+i,elem);
    }
  }
  // - fluxes 
  if ( tdSetup->jac_vel ) { 
    for ( int i = 0; i < num; i++ ) {
      DType denom = std::pow( l_x_c[i]*h_r, tdSetup->dim-1 ); 
      //DType denom = l_x_c[i]*l_x_c[i]*s_r_to2;
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
    int sten = ( tdSetup->chem_rflux == 3 ) ? 3 : 1;
    for ( int i = 0; i < num; i++ ) {
      int i_min = max(0,i-sten);
      DType denom = std::pow( l_x_c[i]*h_r, tdSetup->dim-1 ); 
      //DType denom = l_x_c[i]*l_x_c[i]*s_r_to2;
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

void TDSolveSacado::update_hsol_by_linsol(HVEC& h_var, const int l_offset, const int num, LSPack& lsPack, DType& loc_res) {
  VComp comp(1,0.0);

  loc_res = std::max(loc_res, lsPack.sol->getInfNorm() / std::max( h_var.norm(-1), 1.0 ) );

  for ( int i = 0; i < num; i++ ) {
    comp = lsPack.sol->get(l_offset+i);
    HVEC vsol(1);
    vsol = comp.getVal();
    h_var[i] -= vsol[0];
  }

  return;
}

void TDSolveSacado::update_chems_by_res_imp(const int nt, DType& loc_res) {
  std::vector<std::string>::iterator it;
  for ( it = chemnames.begin(); it != chemnames.end(); it++ ) {
    if ( check_chem_skip( *it ) )
      continue;
    HVEC& h_var = *H_data_var[*it][nt];
    HVEC& h_res = *H_data_res[*it][0];
    DType h_r   = data_r[nt];
    int s_offset = S_data_offset[*it];
    update_chem_by_res_var_imp( h_var, h_res, h_r, varcoord[*it], S_flux_chem[*it], S_srce_chem[*it], varsizes[*it], s_offset, varlspacks[*it], loc_res );
  }
  return;
}

void TDSolveSacado::update_chem_by_res_var_imp(HVEC& h_var, HVEC& h_res, DType h_r, DType* l_x_c, SacadoDepVarArray& s_flx, SacadoDepVarArray& s_src, const int num, const int s_offset, LSPack& lsPack, DType& loc_res) {
  lsPack.mat->reset( 0.0 );
  lsPack.rhs->reset( 0.0 );
  lsPack.sol->reset( 0.0 );

  VComp comp(1,0.0);
  Entry elem(1,0.0);

  // RHS
  //DType s_r_to2 = h_r * h_r;
  for ( int i = 0; i < num; i++ ) {
    comp = h_res[i]/std::pow( l_x_c[i]*h_r, tdSetup->dim-1 );
    //comp = h_res[i]/l_x_c[i]/l_x_c[i]/s_r_to2;
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
    elem -= dtn*s_src[i].dx( s_offset+i )/std::pow( l_x_c[i]*h_r, tdSetup->dim-1 );
    //elem -= dtn*s_src[i].dx( s_offset+i )/l_x_c[i]/l_x_c[i]/s_r_to2;
    lsPack.mat->set(i,i,elem);
  }
  // - fluxes - assuming the stencil within a centered stencil of:
  // - 3 intervals if upwind flux is used
  // - 7 intervals if ppm flux is used 
  int sten = ( tdSetup->chem_rflux == 3 ) ? 3 : 1;
  for ( int i = 0; i < num; i++ ) {
    int i_min = max(0,i-sten);
    DType denom = std::pow( l_x_c[i]*h_r, tdSetup->dim-1 );
    //DType denom = l_x_c[i]*l_x_c[i]*s_r_to2;
    int i_max = min(num-1,i+sten);
    for ( int j = i_min; j <= i_max; j++ ) {
      elem = lsPack.mat->get(i,j);
      elem += (dtn/dx) * (s_flx[i+1].dx( s_offset+j )-s_flx[i].dx( s_offset+j ))/denom;
      lsPack.mat->set(i,j,elem);
    }
  }

  // Solve the linear system
  lsPack.linSolve->solve();

  loc_res = std::max(loc_res, lsPack.sol->getInfNorm() / std::max( h_var.norm(-1), 1.0 ) );

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
  if ( tdSetup->cellmotion ) {
    HVEC** data_vel = H_data_vel["V"];
    if ( tdSetup->dr_opt == 2 ) // For legacy issue, there is no opt 1
      S_data_rp = (*data_vel[to])[Nx]/(1.0-0.25*dx*dx);
    else
      S_data_rp = (*data_vel[to])[Nx];
    
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
  }
  else {
    S_data_rp = 0.0;
    data_r[to] = data_r[fr];
  }

  if ( tdSetup->assemb == 1 )
    S_data_r = data_r[to];
  else
    S_data_r = data_r[fr];

  return;
}

void TDSolveSacado::update_source_geom(int nt) {
  if ( tdSetup->fix_domain ) 
    return;

  if ( tdSetup->pos_fix != 0 )
    adjust_srce_max_geom();
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
  if ( tdSetup->pos_fix != 0 )
    reset_srce_max();
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

void TDSolveSacado::calc_srce_g_model(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  calc_srce_g_base( s_srce, tsol_loc );

  DType lambda = tdCoef->lambda;
  SacadoIndVarArray& s_g = S_data_var["G"];

  for ( int i = 0; i < Nx; i++ ) {
    s_srce[i] += x_c[i] * x_c[i] * S_data_r * S_data_r * ( lambda * s_g[i] );
  }

  if ( tdSetup->pos_fix != 0 ) {
    for ( int i = 0; i < Nx; i++ )
      adjust_srce_max( i, lambda );
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

  if ( tdSetup->pos_fix != 0 ) {
    for ( int i = 0; i < Nx; i++ )
      adjust_srce_max( i, -delta_n );
  }

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

  if ( tdSetup->pos_fix != 0 ) {
    for ( int i = 0; i < Nx; i++ )
      adjust_srce_max( i, lambda - mu );
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

  if ( tdSetup->pos_fix != 0 ) {
    for ( int i = 0; i < Nx; i++ )
      adjust_srce_max( i, -delta_h );
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

  if ( tdSetup->pos_fix != 0 ) {
    for ( int i = 0; i < Nx; i++ )
      adjust_srce_max( i, -delta_n );
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
    if ( tdSetup->logistic == 0 )
      s_lambda = icalc_scal_logist_mm(lambda,s_g[i],1.0,b_prime,b,s_p2[i]);
    else
      s_lambda = icalc_scal_logist_mm(lambda,0.0,1.0,b_prime,b,s_p2[i]);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * std::pow(x_c[i],p) );
    s_mu = icalc_scal_logist_mm(loc_mu,0.0,1.0,a_prime,a,s_p1[i]);
    s_srce[i] += (x_c[i] * x_c[i] * S_data_r * S_data_r) * ( (s_lambda-s_mu-omega*s_n[i]) * s_g[i] );
    if ( tdSetup->pos_fix != 0 )
      adjust_srce_max( i, s_lambda.val()-s_mu.val()-omega*s_n[i].val() );
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

  if ( tdSetup->pos_fix != 0 ) {
    for ( int i = 0; i < Nx; i++ )
      adjust_srce_max( i, -delta_h );
  }

  return;
}

void TDSolveSacado::calc_srce_n_opn(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  // The same as in the pdgf model
  calc_srce_n_pdgf( s_srce, tsol_loc );
  return;
}

void TDSolveSacado::calc_srce_b_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  reset_sacado_dep_array(s_srce);
  if (!in_oncolytic_stage[0]) return;

  DType theta   = tdCoef->theta;
  DType mm_ceil = tdCoef->probVarParVals[prob_name]["B"]["mm_ceil_grow_pos"];
  DType mm_half = tdCoef->probVarParVals[prob_name]["B"]["mm_half"];
  DType coef_vi = tdCoef->probVarParVals[prob_name]["B"]["coef_vi_neg"];
  DType coef_h  = tdCoef->probVarParVals[prob_name]["B"]["coef_h_neg"] * theta;
  DType coef_mb = tdCoef->probVarParVals[prob_name]["B"]["mm_ceil_coef_mb_neg"] * theta;
  SacadoIndVarArray& s_b  = S_data_var["B"];
  SacadoIndVarArray& s_vi = S_data_var["Vi"];
  SacadoIndVarArray& s_h  = S_data_var["H"];
  SacadoIndVarArray& s_mb = S_data_var["Mb"];
  SacadoIndVarArray& s_a  = S_data_var["A"];

  int    num   = varsizes["B"];
  DType* l_x_c = varcoord["B"];

  SacadoDepVar coef_b(0.0);

  for ( int i = 0; i < num; i++ ) {
    coef_b = (mm_ceil*s_a[i])/(mm_half+s_a[i]) - coef_vi*s_vi[i] - coef_h*s_h[i] - 
        ((coef_mb*s_a[i])/(mm_half+s_a[i]))*s_mb[i];
    s_srce[i] += (l_x_c[i] * l_x_c[i] * S_data_r * S_data_r) * coef_b * s_b[i];
    if ( tdSetup->pos_fix != 0 )
      adjust_srce_max( i, coef_b.val() );
  }

  return;
}

void TDSolveSacado::calc_srce_h_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  reset_sacado_dep_array(s_srce);
  if (!in_oncolytic_stage[1]) return;

  DType theta   = tdCoef->theta;
  DType mm_ceil   = tdCoef->probVarParVals[prob_name]["H"]["mm_ceil_grow_pos"];
  DType mm_half   = tdCoef->probVarParVals[prob_name]["H"]["mm_half"];
  DType coef_mb   = tdCoef->probVarParVals[prob_name]["H"]["mm_ceil_coef_mb_neg"] * theta;
  DType coef_mv   = tdCoef->probVarParVals[prob_name]["H"]["mm_ceil_coef_mv_neg"] * theta;
  DType decay     = tdCoef->probVarParVals[prob_name]["H"]["grow_neg"];
  DType prod_vi_b = tdCoef->probVarParVals[prob_name]["H"]["prod_vi_b_pos"];
  DType coef_b    = tdCoef->probVarParVals[prob_name]["H"]["coef_b_pos"] * theta;
  SacadoIndVarArray& s_h  = S_data_var["H"];
  SacadoIndVarArray& s_b  = S_data_var["B"];
  SacadoIndVarArray& s_vi = S_data_var["Vi"];
  SacadoIndVarArray& s_mb = S_data_var["Mb"];
  SacadoIndVarArray& s_mv = S_data_var["Mv"];
  SacadoIndVarArray& s_a  = S_data_var["A"];

  int    num   = varsizes["H"];
  DType* l_x_c = varcoord["H"];

  SacadoDepVar coef_h(0.0);

  for ( int i = 0; i < num; i++ ) {
    coef_h = (mm_ceil*s_a[i])/(mm_half+s_a[i]) + coef_b*s_b[i] - ((coef_mb*s_a[i])/(mm_half+s_a[i]))*s_mb[i] - 
      ((coef_mv*s_a[i])/(mm_half+s_a[i]))*s_mv[i] - decay;
    s_srce[i] += (l_x_c[i] * l_x_c[i] * S_data_r * S_data_r) * (coef_h * s_h[i] + prod_vi_b * s_vi[i] * s_b[i] );
    if ( tdSetup->pos_fix != 0 )
      adjust_srce_max( i, coef_h.val() );
  }

  return;
}

void TDSolveSacado::calc_srce_n_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  reset_sacado_dep_array(s_srce);
  if (!in_oncolytic_stage[1]) return;

  DType theta     = tdCoef->theta;
  DType mm_half   = tdCoef->probVarParVals[prob_name]["N"]["mm_half"];
  DType prod_mb_b = tdCoef->probVarParVals[prob_name]["N"]["mm_ceil_prod_mb_b_pos"] * theta;
  DType prod_mb_h = tdCoef->probVarParVals[prob_name]["N"]["mm_ceil_prod_mb_h_pos"] * theta;
  DType prod_mv_h = tdCoef->probVarParVals[prob_name]["N"]["mm_ceil_prod_mv_h_pos"] * theta;
  DType prod_h    = tdCoef->probVarParVals[prob_name]["N"]["prod_h_pos"];
  SacadoIndVarArray& s_n  = S_data_var["N"];
  SacadoIndVarArray& s_h  = S_data_var["H"];
  SacadoIndVarArray& s_b  = S_data_var["B"];
  SacadoIndVarArray& s_mb = S_data_var["Mb"];
  SacadoIndVarArray& s_mv = S_data_var["Mv"];
  SacadoIndVarArray& s_a  = S_data_var["A"];

  int    num   = varsizes["N"];
  DType* l_x_c = varcoord["N"];

  for ( int i = 0; i < num; i++ ) {
    s_srce[i] += (l_x_c[i] * l_x_c[i] * S_data_r * S_data_r) *
      ( ((prod_mb_b*s_a[i])/(mm_half+s_a[i])) * s_mb[i] * s_b[i] + 
        ((prod_mb_h*s_a[i])/(mm_half+s_a[i])) * s_mb[i] * s_h[i] + 
        ((prod_mv_h*s_a[i])/(mm_half+s_a[i])) * s_mv[i] * s_h[i] + 
        prod_h * s_h[i] );
  }

  return;
}

void TDSolveSacado::calc_srce_g_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  reset_sacado_dep_array(s_srce);

  DType theta  = tdCoef->theta;
  DType coef_b = tdCoef->probVarParVals[prob_name]["G"]["coef_b_pos"] * theta;
  DType decay  = tdCoef->probVarParVals[prob_name]["G"]["grow_neg"];
  SacadoIndVarArray& s_g = S_data_var["G"];
  SacadoIndVarArray& s_b = S_data_var["B"];

  int    num   = varsizes["G"];
  DType* l_x_c = varcoord["G"];

  SacadoDepVar coef_g(0.0);

  for ( int i = 0; i < num; i++ ) {
    coef_g = coef_b * s_b[i] - decay;
    s_srce[i] += (l_x_c[i] * l_x_c[i] * S_data_r * S_data_r) * ( coef_g * s_g[i] );
    if ( tdSetup->pos_fix != 0 )
      adjust_srce_max( i, coef_g.val() );
  }

  return;
}

void TDSolveSacado::calc_srce_mv_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  reset_sacado_dep_array(s_srce);
  if (!in_oncolytic_stage[2]) return;

  DType theta   = tdCoef->theta;
  DType mm_ceil = tdCoef->probVarParVals[prob_name]["Mv"]["mm_ceil_grow_pos"];
  DType mm_half = tdCoef->probVarParVals[prob_name]["Mv"]["mm_half"];
  DType decay   = tdCoef->probVarParVals[prob_name]["Mv"]["grow_neg"];
  DType coef_h  = tdCoef->probVarParVals[prob_name]["Mv"]["coef_h_pos"] * theta;
  SacadoIndVarArray& s_mv = S_data_var["Mv"];
  SacadoIndVarArray& s_h  = S_data_var["H"];
  SacadoIndVarArray& s_a  = S_data_var["A"];

  int    num   = varsizes["Mv"];
  DType* l_x_c = varcoord["Mv"];

  SacadoDepVar coef_mv(0.0);

  for ( int i = 0; i < num; i++ ) {
    coef_mv = coef_h*s_h[i] + (mm_ceil*s_a[i])/(mm_half+s_a[i]) - decay;
    s_srce[i] += (l_x_c[i] * l_x_c[i] * S_data_r * S_data_r) * (coef_mv * s_mv[i]);
    if ( tdSetup->pos_fix != 0 )
      adjust_srce_max( i, coef_mv.val() );
  }

  return;
}

void TDSolveSacado::calc_srce_mb_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  reset_sacado_dep_array(s_srce);
  if (!in_oncolytic_stage[2]) return;

  DType theta   = tdCoef->theta;
  DType mm_ceil = tdCoef->probVarParVals[prob_name]["Mb"]["mm_ceil_grow_pos"];
  DType mm_half = tdCoef->probVarParVals[prob_name]["Mb"]["mm_half"];
  DType decay   = tdCoef->probVarParVals[prob_name]["Mb"]["grow_neg"];
  DType coef_b  = tdCoef->probVarParVals[prob_name]["Mb"]["coef_b_pos"] * theta;
  DType coef_h  = tdCoef->probVarParVals[prob_name]["Mb"]["coef_h_pos"] * theta;
  SacadoIndVarArray& s_mb = S_data_var["Mb"];
  SacadoIndVarArray& s_h  = S_data_var["H"];
  SacadoIndVarArray& s_b  = S_data_var["B"];
  SacadoIndVarArray& s_a  = S_data_var["A"];

  int    num   = varsizes["Mb"];
  DType* l_x_c = varcoord["Mb"];

  SacadoDepVar coef_mb(0.0);

  for ( int i = 0; i < num; i++ ) {
    coef_mb = coef_b*s_b[i] + coef_h*s_h[i] + (mm_ceil*s_a[i])/(mm_half+s_a[i]) - decay;
    s_srce[i] += (l_x_c[i] * l_x_c[i] * S_data_r * S_data_r) * (coef_mb * s_mb[i]);
    if ( tdSetup->pos_fix != 0 )
      adjust_srce_max( i, coef_mb.val() );
  }

  return;
}

void TDSolveSacado::calc_srce_b1_potassium(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  reset_sacado_dep_array(s_srce);

  DType decay = tdCoef->probVarParVals[prob_name]["B1"]["grow_neg_virus"];
  DType k     = tdCoef->probVarParVals[prob_name]["P"]["slope"];
  DType Ps    = tdCoef->probVarParVals[prob_name]["P"]["threshold"];
  SacadoIndVarArray& s_b1 = S_data_var["B1"];
  SacadoIndVarArray& s_p  = S_data_var["P"];
  SacadoIndVarArray& s_v  = S_data_var["V"];

  int    num   = varsizes["B1"];
  DType* l_x_c = varcoord["B1"];

  for ( int i = 0; i < num; i++ ) {
    s_srce[i] += std::pow( l_x_c[i] * S_data_r, tdSetup->dim-1 ) * 
      ( - decay * icalc_scal_block_potassium( k, Ps, s_p[i] ) * s_b1[i] * s_v[i] );
  }

  return;
}

void TDSolveSacado::calc_srce_b2_potassium(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  reset_sacado_dep_array(s_srce);

  DType decay_b1 = tdCoef->probVarParVals[prob_name]["B1"]["grow_neg_virus"];
  DType k        = tdCoef->probVarParVals[prob_name]["P"]["slope"];
  DType Ps       = tdCoef->probVarParVals[prob_name]["P"]["threshold"];
  DType decay_b2 = tdCoef->probVarParVals[prob_name]["B2"]["grow_neg"];
  DType p_block = ( tdCoef->probVarParVals[prob_name]["B2"].find( "grow_block" ) == 
                    tdCoef->probVarParVals[prob_name]["B2"].end() ) ? 0.0 : 
                    tdCoef->probVarParVals[prob_name]["B2"]["grow_block"];
  SacadoIndVarArray& s_b1 = S_data_var["B1"];
  SacadoIndVarArray& s_b2 = S_data_var["B2"];
  SacadoIndVarArray& s_p  = S_data_var["P"];
  SacadoIndVarArray& s_v  = S_data_var["V"];

  int    num   = varsizes["B2"];
  DType* l_x_c = varcoord["B2"];

  for ( int i = 0; i < num; i++ ) {
    s_srce[i] += std::pow( l_x_c[i] * S_data_r, tdSetup->dim-1 ) * 
      ( decay_b1 * icalc_scal_block_potassium( k, Ps, s_p[i] ) * s_b1[i] * s_v[i] - 
        decay_b2 * icalc_scal_block_potassium( k, Ps, s_p[i], p_block ) * s_b2[i] );
  }

  return;
}

void TDSolveSacado::calc_srce_b3_potassium(SacadoDepVarArray& s_srce, const DType tsol_loc) {
  reset_sacado_dep_array(s_srce);

  DType k        = tdCoef->probVarParVals[prob_name]["P"]["slope"];
  DType Ps       = tdCoef->probVarParVals[prob_name]["P"]["threshold"];
  DType decay_b2 = tdCoef->probVarParVals[prob_name]["B2"]["grow_neg"];
  DType p_block = ( tdCoef->probVarParVals[prob_name]["B2"].find( "grow_block" ) == 
                    tdCoef->probVarParVals[prob_name]["B2"].end() ) ? 0.0 : 
                    tdCoef->probVarParVals[prob_name]["B2"]["grow_block"];
  SacadoIndVarArray& s_b2 = S_data_var["B2"];
  SacadoIndVarArray& s_p  = S_data_var["P"];

  int    num   = varsizes["B3"];
  DType* l_x_c = varcoord["B3"];

  for ( int i = 0; i < num; i++ ) {
    s_srce[i] += std::pow( l_x_c[i] * S_data_r, tdSetup->dim-1 ) * ( decay_b2 * 
        icalc_scal_block_potassium( k, Ps, s_p[i], p_block ) * s_b2[i] );
  }

  return;
}

void TDSolveSacado::compute_source_chem_prob(const int nt, const DType tsol_loc) {
  std::vector<std::string>::iterator it_name;
  for ( it_name = chemnames.begin(); it_name != chemnames.end(); it_name++ ) {
    if ( check_chem_skip( *it_name ) )
      continue;
    (this->*chempars[*it_name].srceFun)( *it_name, nt, tsol_loc );
    //calc_srce_chem_var( *it_name, nt, tsol_loc );
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

void TDSolveSacado::calc_chem_srce_vi_oncolytic(std::string chemname, const int nt, const DType tsol_loc) {
  // assert( chemname.compare("Vi") == 0 );
  SacadoDepVarArray& s_srce  = S_srce_chem[chemname];
  DType*             l_x_c   = varcoord[chemname];
  ChemPar&           chemPar = chempars[chemname];
  SacadoIndVarArray& s_chem  = S_data_var[chemname];
  SacadoIndVarArray& s_h     = S_data_var["H"];

  DType theta  = tdCoef->theta;
  DType prod_h = tdCoef->probVarParVals[prob_name][chemname]["prod_growth_h_pos"] * 
                 tdCoef->probVarParVals[prob_name][chemname]["growth_h"] * theta;
  DType decay  = chemPar.decay;

  int num = s_chem.size();

  for ( int i = 0; i < num; i++ ) {
    s_srce[i] += (l_x_c[i] * l_x_c[i] * data_r[nt] * data_r[nt]) * 
      ( prod_h * s_h[i] - decay * s_chem[i] );
  }

  return;
}

void TDSolveSacado::calc_chem_srce_a_oncolytic(std::string chemname, const int nt, const DType tsol_loc) {
  // assert( chemname.compare("A") == 0 );
  SacadoDepVarArray& s_srce  = S_srce_chem[chemname];
  DType*             l_x_c   = varcoord[chemname];
  ChemPar&           chemPar = chempars[chemname];
  SacadoIndVarArray& s_chem  = S_data_var[chemname];
  SacadoIndVarArray& s_b     = S_data_var["B"];
  SacadoIndVarArray& s_h     = S_data_var["H"];
  SacadoIndVarArray& s_mv    = S_data_var["Mv"];
  SacadoIndVarArray& s_mb    = S_data_var["Mb"];
  SacadoIndVarArray& s_g     = S_data_var["G"];

  DType theta   = tdCoef->theta;
  DType mm_half = tdCoef->probVarParVals[prob_name][chemname]["mm_half"];
  DType coef_b  = tdCoef->probVarParVals[prob_name][chemname]["coef_mm_b_neg"] * 
                  tdCoef->probVarParVals[prob_name][chemname]["mm_ceil_b"] * theta;
  DType coef_h  = tdCoef->probVarParVals[prob_name][chemname]["coef_mm_h_neg"] * 
                  tdCoef->probVarParVals[prob_name][chemname]["mm_ceil_h"] * theta;
  DType coef_mv = tdCoef->probVarParVals[prob_name][chemname]["coef_mm_mv_neg"] * 
                  tdCoef->probVarParVals[prob_name][chemname]["mm_ceil_mv"] * theta;
  DType coef_mb = tdCoef->probVarParVals[prob_name][chemname]["coef_mm_mb_neg"] * 
                  tdCoef->probVarParVals[prob_name][chemname]["mm_ceil_mb"] * theta;
  DType prod_g  = tdCoef->probVarParVals[prob_name][chemname]["prod_g_neg"] * theta;
  DType prod    = tdCoef->probVarParVals[prob_name][chemname]["prod_constant"];

  int num = s_chem.size();

  for ( int i = 0; i < num; i++ ) {
    s_srce[i] += (l_x_c[i] * l_x_c[i] * data_r[nt] * data_r[nt]) *
      ( prod - prod_g * s_g[i] - (s_chem[i]/(mm_half+s_chem[i])) * 
        (coef_b * s_b[i] + coef_h * s_h[i] + coef_mv * s_mv[i] + coef_mb * s_mb[i]) );
  }

  return;
}

void TDSolveSacado::calc_chem_srce_p_potassium(std::string chemname, const int nt, const DType tsol_loc) {
  // assert( chemname.compare("P") == 0 );
  SacadoDepVarArray& s_srce  = S_srce_chem[chemname];
  DType*             l_x_c   = varcoord[chemname];
  ChemPar&           chemPar = chempars[chemname];
  SacadoIndVarArray& s_chem  = S_data_var[chemname]; // This is just "P"
  SacadoIndVarArray& s_b2    = S_data_var["B2"];

  DType theta   = tdCoef->theta;
  DType k       = tdCoef->probVarParVals[prob_name]["P"]["slope"];
  DType Ps      = tdCoef->probVarParVals[prob_name]["P"]["threshold"];
  DType prod_b2 = tdCoef->probVarParVals[prob_name][chemname]["burst"] *
                  tdCoef->probVarParVals[prob_name]["B2"]["grow_neg"] * theta;
  DType p_block = ( tdCoef->probVarParVals[prob_name]["B2"].find( "grow_block" ) == 
                    tdCoef->probVarParVals[prob_name]["B2"].end() ) ? 0.0 : 
                    tdCoef->probVarParVals[prob_name]["B2"]["grow_block"];
  DType decay   = chemPar.decay;

  int num = s_chem.size();
  
  for ( int i = 0; i < num; i++ ) {
    s_srce[i] += std::pow( l_x_c[i] * data_r[nt], tdSetup->dim-1 ) *
      ( prod_b2 * icalc_scal_block_potassium( k, Ps, s_chem[i], p_block ) * s_b2[i] - decay * s_chem[i] );
  }

  return;
}

void TDSolveSacado::calc_chem_srce_v_potassium(std::string chemname, const int nt, const DType tsol_loc) {
  // assert( chemname.compare("V") == 0 );
  SacadoDepVarArray& s_srce  = S_srce_chem[chemname];
  DType*             l_x_c   = varcoord[chemname];
  ChemPar&           chemPar = chempars[chemname];
  SacadoIndVarArray& s_chem  = S_data_var[chemname];
  SacadoIndVarArray& s_p     = S_data_var["P"];
  SacadoIndVarArray& s_b1    = S_data_var["B1"];
  SacadoIndVarArray& s_b2    = S_data_var["B2"];

  DType theta   = tdCoef->theta;
  DType prod_b2 = tdCoef->probVarParVals[prob_name][chemname]["burst"] *
                  tdCoef->probVarParVals[prob_name]["B2"]["grow_neg"] * theta;
  DType p_block = ( tdCoef->probVarParVals[prob_name]["B2"].find( "grow_block" ) == 
                    tdCoef->probVarParVals[prob_name]["B2"].end() ) ? 0.0 : 
                    tdCoef->probVarParVals[prob_name]["B2"]["grow_block"];
  DType decay_b = tdCoef->probVarParVals[prob_name]["B1"]["grow_neg_virus"] * theta;
  DType decay   = chemPar.decay;
  DType k       = tdCoef->probVarParVals[prob_name]["P"]["slope"];
  DType Ps      = tdCoef->probVarParVals[prob_name]["P"]["threshold"];

  int num = s_chem.size();

  for ( int i = 0; i < num; i++ ) {
    s_srce[i] += std::pow( l_x_c[i] * data_r[nt], tdSetup->dim-1 ) *
      ( icalc_scal_block_potassium( k, Ps, s_p[i] ) * ( prod_b2 * 
        icalc_scal_block_potassium( k, Ps, s_p[i], p_block ) * s_b2[i] - decay_b * s_b1[i] * s_chem[i] ) 
        - decay * s_chem[i] );
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

DType TDSolveSacado::get_cell_perc(std::string cellname, const int nt) {
  DType tot_cell = 0.0;
  DType tot_all  = 0.0;

  //DType pi4 = 4.0 * PI;
  
  for ( std::vector<std::string>::iterator it = cellnames.begin();
        it != cellnames.end(); it++ ) {
    int num = varsizes[*it];
    DType* l_x_c = varcoord[*it];
    DType tot_loc = 0.0;
    DType dr = data_r[nt] * dx;
    for ( int i = 0; i < num; i++ ) {
      DType r_loc = data_r[nt] * l_x_c[i];
      tot_loc += dr * tdSetup->sweight * std::pow(r_loc,tdSetup->dim-1) * (*H_data_var[*it][nt])[i];
      //tot_loc += dr * pi4 * r_loc * r_loc * (*H_data_var[*it][nt])[i];
    }
    tot_all += tot_loc;
    if ( *it == cellname )
      tot_cell = tot_loc;
  }

  return tot_cell / tot_all;
}

DType TDSolveSacado::calc_treatment_P1(const DType tsol_loc) {
  std::vector<DType> t0 = tdCoef->treatP1_t0;
  DType T = tdCoef->treatP1_T;
  DType h = tdCoef->treatP1_h;

  if ( t0.size() == 0 )
    return 0.0;

  int pass_count = 0;
  for ( unsigned i = 0; i < t0.size(); i++ ) 
    if ( tsol_loc <= t0[i] )
      break;
    else
      pass_count++;

  if ( pass_count%2 == 0 )
    return 0.0;
  else {
    if ( h <= 0.0 )
      return T;
    else
      return T * ( 1.0 + std::cos(2.0*PI*tsol_loc/h) );
  }
}

void TDSolveSacado::reset_srce_max(void) {
  *srce_max[0] = 0.0;
  return;
}

void TDSolveSacado::adjust_srce_max(const int i, DType lip_loc) {
  if ( lip_loc < 0.0 ) 
    (*srce_max[0])[i] = std::max( (*srce_max[0])[i], -lip_loc );
  return;
}

void TDSolveSacado::adjust_srce_max_geom(void) {
  if (S_data_rp.val() < 0.0)
    return;
  else {
    for ( int i = 0; i < Nx; i++ ) {
      (*srce_max[0])[i] += S_data_rp.val()/S_data_r.val();
    }
  }
  return;
}

void TDSolveSacado::adjust_pos_ppm_lim( SacadoDepVarArray* s_alpha ) {
  if ( S_data_rp.val() < 0 )
    for ( int i = 0; i < Nx-1; i++ ) {
      DType coef = - (x_f[i+1]/x_c[i+1])*(x_f[i+1]/x_c[i+1])/fabs(S_data_rp.val())/x_c[i+1];
      DType loc_x_r = x_c[i+1]*S_data_r.val();
      DType loc_x_r3 = loc_x_r*loc_x_r*loc_x_r;
      DType loc_x = x_c[i]*S_data_r.val();
      DType loc_x3 = loc_x*loc_x*loc_x;
      for ( std::vector<std::string>::iterator it = cellnames.begin(); 
            it != cellnames.end(); it++ ) {
        SacadoIndVarArray& s_var = S_data_var[*it];
        SacadoDepVarArray& s_rec = S_rec_ppm_cell[*it];
        DType loc_v = 0.0;
        for ( std::set<std::string>::iterator vit = cellvels[*it].begin();
              vit != cellvels[*it].end(); vit++ ) {
          SacadoDepVarArray& s_v = S_data_vel[*vit];
          loc_v += min(0.0,s_v[i+1].val());
        }
        s_alpha[0][i] = min( s_alpha[0][i].val(), 
            coef * loc_v * std::fabs(loc_x_r3*s_var[i+1].val()-loc_x3*s_var[i].val()) /
            max( std::fabs(s_rec[i+1].val()-loc_x3*s_var[i].val()), 1.0e-6 ) );
      }
    }
  if ( S_data_rp.val() > 0 )
    for ( int i = 1; i < Nx; i++ ) {
      DType coef = (x_f[i]/x_c[i-1])*(x_f[i]/x_c[i-1])/fabs(S_data_rp.val())/x_c[i-1];
      DType loc_x_l = x_c[i-1]*S_data_r.val();
      DType loc_x_l3 = loc_x_l*loc_x_l*loc_x_l;
      DType loc_x = x_c[i]*S_data_r.val();
      DType loc_x3 = loc_x*loc_x*loc_x;
      for ( std::vector<std::string>::iterator it = cellnames.begin(); 
            it != cellnames.end(); it++ ) {
        SacadoIndVarArray& s_var = S_data_var[*it];
        SacadoDepVarArray& s_rec = S_rec_ppm_cell[*it];
        DType loc_v = 0.0;
        for ( std::set<std::string>::iterator vit = cellvels[*it].begin();
              vit != cellvels[*it].end(); vit++ ) {
          SacadoDepVarArray& s_v = S_data_vel[*vit];
          loc_v += max(0.0,s_v[i].val());
        }
        s_alpha[1][i] = max( s_alpha[0][i].val(), 
            (1.0/coef/max(loc_v,1.0e-6)) * std::fabs(loc_x3*s_var[i].val() - s_rec[i].val()) /
            max( std::fabs(loc_x3*s_var[i].val() - loc_x_l3*s_var[i-1].val()), 1.0e-6) );
      }
    }
  return;
}

void TDSolveSacado::update_ale_grids(const int nt) {
  //if ( tdSetup->chemo_L < data_r[nt] * (1.0*Na/Nx) ) {
  if ( tdSetup->chemo_L < data_r[nt] * 2 ) {
    throw std::runtime_error( "TDSolveSacado::update_ale_grids - range is too short!\n" );
  }

  DType dr = dx * data_r[nt];

  // First part, r from 0 to R, uniform grid
  for ( int i = 0; i <= Nx; i++ ) {
    x_ale_f[nt][i] = i * dr;
  }

  // Second part, r from R to chemo_L, geometric grid
  // (g^{Na-Nx}-1)/(g-1) = (L-R)/(R/Nx)
  DType g;
  if ( !calc_ale_growth(g,Na-Nx,(tdSetup->chemo_L-data_r[nt])/dr) ) {
    throw std::runtime_error( "TDSolveSacado::calc_ale_growth - bisection fail.\n" );
  }
  DType g_scal = 1.0;
  for ( int i = Nx+1; i < Na; i++ ) {
    x_ale_f[nt][i] = x_ale_f[nt][i-1] + dr * g_scal;
    g_scal *= g;
  }
  x_ale_f[nt][Na] = tdSetup->chemo_L;

  for ( int i = 0; i < Na; i++ ) {
    x_ale_c[nt][i] = 0.5 * (x_ale_f[nt][i] + x_ale_f[nt][i+1]);
  }

  for ( int i = 0; i < Na; i++ ) {
    dx_ale[nt][i] = x_ale_f[nt][i+1] - x_ale_f[nt][i];
  }
  return;
}

bool TDSolveSacado::calc_ale_growth(DType& g, const int num, const DType val) {
  // Bisection
  const int   iter_max = 100;
  const DType loc_eps = 1e-10;
  DType g_l = std::pow(val/num,1.0/(num-1.0));
  DType g_r = std::pow(val,1.0/(num-1.0));
  g = 0.5 * (g_l+g_r);
  DType f_l = std::pow(g_l,num)-1.0-val*(g_l-1.0);
  DType f_r = std::pow(g_r,num)-1.0-val*(g_r-1.0);

  int niter = 0;

  if ( f_l * f_r > 0.0 )
    return false;

  while ( std::fabs(g_r-g_l) > loc_eps && niter < iter_max ) {
    niter++;
    DType g_m = 0.5 * ( g_l + g_r );
    DType f_m = std::pow(g_m,num)-1.0-val*(g_m-1.0); 
    if ( std::fabs( f_m ) < loc_eps )
      break;
    else if ( f_m * f_l < 0.0 )
      g_r = g_m;
    else
      g_l = g_m;
  }

  g = 0.5 * ( g_l + g_r );

  return niter <  iter_max;
}

void TDSolveSacado::print_data_array(std::string message, DType* loc_d, const int num) {
  std::cout << message << std::endl;
  for ( int i = 0; i < num; i++ ) {
    std::cout << "\t" << i << ": " << loc_d[i] << std::endl;
  }
  return;
}

bool TDSolveSacado::check_chem_skip(std::string name) {
  for ( std::vector<std::string>::iterator it = ale_chemnames.begin();
        it != ale_chemnames.end(); it++ ) {
    if ( *it == name )
      return true;
  }
  return false;
}

void TDSolveSacado::compute_ale_source_chem_prob(const int nt, const DType tsol_loc) {
  std::vector<std::string>::iterator it_name;
  for ( it_name = ale_chemnames.begin(); it_name != ale_chemnames.end(); it_name++ ) {
    calc_ale_srce_chem_var( *it_name, nt, tsol_loc );
  }
  return;
}

void TDSolveSacado::calc_ale_srce_chem_var(std::string chemname, const int nt, const DType tsol_loc) {
  // May need to incorporate the additional source such as injection in ChemPar
  SacadoDepVarArray& s_srce  = S_srce_chem[chemname];
  ChemPar&           chemPar = chempars[chemname];
  SacadoIndVarArray& s_chem  = S_data_var[chemname];
  SacadoIndVarArray& s_cell  = S_data_var[chemPar.driver];

  for ( int i = 0; i < Na; i++ ) {
    s_srce[i] += ( - chemPar.decay * s_chem[i] );
  }

  if ( chemPar.test_opt == 0 ) { // Regular method
    for ( int i = 0; i < Nx; i++ ) {
      s_srce[i] += ( chemPar.mm_ceil * s_cell[i] ) / ( chemPar.mm_half + s_cell[i] );
    }
  }
  else if ( chemPar.test_opt == 1 ) {
    Sacado::Fad::DFad<DType> tot_cell = 0.0;
    Sacado::Fad::DFad<DType> tot_vlmn = 0.0;
    for ( int i = 0; i < Nx; i++ ) {
      tot_cell = x_c[i] * x_c[i] * s_cell[i];
      tot_vlmn = x_c[i] * x_c[i];
    }
    tot_cell = tot_cell/tot_vlmn;
    Sacado::Fad::DFad<DType> growth = ( chemPar.mm_ceil * tot_cell ) / ( chemPar.mm_half + tot_cell );
    for ( int i = 0; i < Nx; i++ ) {
      s_srce[i] += growth * s_chem[i];
    }
  }
  return;
}

void TDSolveSacado::compute_ale_chems_dflux(const int nt) {
  for ( std::vector<std::string>::iterator it = ale_chemnames.begin();
        it != ale_chemnames.end(); it++ ) {
    calc_ale_dflux_chem_var( *it, nt );
  }
  return;
}

void TDSolveSacado::calc_ale_dflux_chem_var(std::string chemname, const int nt) {
  SacadoIndVarArray& s_var = S_data_var[chemname];
  SacadoDepVarArray& s_flx = S_flux_chem_d[chemname];
  reset_sacado_dep_array(s_flx);
  DType dcoef = chempars[chemname].diffusion;
  bool  dbc   = chempars[chemname].dirichlet_bc;

  // - left boundary
  {
    s_flx[0] = 0.0;
  }
  // - interior points
  for ( int i = 1; i < Na; i++ ) {
    s_flx[i] = - dcoef * x_ale_f[nt][i] * x_ale_f[nt][i] * 
      ( (s_var[i] - s_var[i-1])/(x_ale_c[nt][i] - x_ale_c[nt][i-1]) );
  }
  // - right boundary
  {
    if ( dbc ) // Homogeneous D. b.c.
      s_flx[Na] = - dcoef * x_ale_f[nt][Na] * x_ale_f[nt][Na] * (0.0 - s_var[Na-1]) / (0.5*dx_ale[nt][Na-1]);
    else // Homogeneous N. b.c.
      s_flx[Na] = 0.0;
  }

  return;
}

void TDSolveSacado::compute_ale_chems_rflux(const int fr, const int to ) {
  for ( std::vector<std::string>::iterator it = ale_chemnames.begin();
        it != ale_chemnames.end(); it++ ) {
    calc_ale_rflux_chem_var( *it, fr, to );
  }
  return;
}

void TDSolveSacado::calc_ale_rflux_chem_var(std::string chemname, const int fr, const int to ) {
  SacadoIndVarArray& s_var = S_data_var[chemname];
  SacadoDepVarArray& s_flx = S_flux_chem_r[chemname];
  reset_sacado_dep_array(s_flx);
  bool dbc = chempars[chemname].dirichlet_bc;

  // - first cell
  {
    DType loc_v = -(x_ale_c[to][0]-x_ale_c[fr][0])/dtn;
    if (loc_v>=0.0)
      s_flx[0] = 0.0;
    else
      s_flx[0] = loc_v * (s_var[1]-s_var[0])/(x_ale_c[to][1]-x_ale_c[to][0]);
  }
  // - interior cells
  for ( int i = 1; i < Na-1; i++ ) {
    DType loc_v = -(x_ale_c[to][i]-x_ale_c[fr][i])/dtn;
    if (loc_v>=0.0)
      s_flx[i] = loc_v * (s_var[i]-s_var[i-1])/(x_ale_c[to][i]-x_ale_c[to][i-1]);
    else
      s_flx[i] = loc_v * (s_var[i+1]-s_var[i])/(x_ale_c[to][i+1]-x_ale_c[to][i]);
  }
  // - last cell
  {
    DType loc_v = -(x_ale_c[to][Na-1]-x_ale_c[fr][Na-1])/dtn;
    if (loc_v>=0.0)
      s_flx[Na-1] = loc_v * (s_var[Na-1]-s_var[Na-2])/(x_ale_c[to][Na-1]-x_ale_c[to][Na-2]);
    else {
      if ( dbc ) { // Homogeneous D. b.c.
        s_flx[Na-1] = loc_v * (0.0-s_var[Na-1])/(0.5*dx_ale[to][Na-1]);
      }
      else { // Homogenous N. b.c.
        s_flx[Na-1] = 0.0;
      }
    }
  }

  return;
}

void TDSolveSacado::compute_ale_chems_residual(const int fr, const int to) {
  for ( std::vector<std::string>::iterator it = ale_chemnames.begin();
        it != ale_chemnames.end(); it++ ) {
    HVEC& h_var = *H_data_var[*it][fr];
    HVEC& h_res = *H_data_res[*it][0];
    SacadoIndVarArray& s_var   = S_data_var[*it];
    SacadoDepVarArray& s_flx_r = S_flux_chem_r[*it];
    SacadoDepVarArray& s_flx_d = S_flux_chem_d[*it];
    SacadoDepVarArray& s_src   = S_srce_chem[*it];
    h_res = 0.0;
    assemble_ale_residual_var( h_res, fr, to, h_var, s_var, s_flx_r, s_flx_d, s_src );
  }
}

void TDSolveSacado::assemble_ale_residual_var(HVEC& h_res, const int fr, const int to, HVEC& h_var, SacadoIndVarArray& s_var, SacadoDepVarArray& s_flx_r, SacadoDepVarArray& s_flx_d, SacadoDepVarArray& s_src) {
  for ( int i = 0; i < Na; i++ ) {
    h_res[i] += s_var[i].val() - h_var[i]; // Inertia term
    h_res[i] -= dtn * s_src[i].val(); // Source term
    h_res[i] += (dtn/dx_ale[to][i]) * (s_flx_d[i+1].val() - s_flx_d[i].val())/(x_ale_c[to][i]*x_ale_c[to][i]); // Diffusion
    h_res[i] += dtn * s_flx_r[i].val();
  }
  return;
}

void TDSolveSacado::update_ale_chems_by_res_imp(const int nt) {
  std::vector<std::string>::iterator it;
  for ( it = ale_chemnames.begin(); it != ale_chemnames.end(); it++ ) {
    HVEC& h_var = *H_data_var[*it][nt];
    HVEC& h_res = *H_data_res[*it][0];
    SacadoDepVarArray& s_flx_r = S_flux_chem_r[*it];
    SacadoDepVarArray& s_flx_d = S_flux_chem_d[*it];
    SacadoDepVarArray& s_src   = S_srce_chem[*it];
    int s_offset = S_data_offset[*it];
    update_ale_chem_by_res_var_imp( h_var, h_res, s_flx_r, s_flx_d, s_src, s_offset, varlspacks[*it], nt );
  }
  return;
}

void TDSolveSacado::update_ale_chem_by_res_var_imp(HVEC& h_var, HVEC& h_res, SacadoDepVarArray& s_flx_r, SacadoDepVarArray& s_flx_d, SacadoDepVarArray& s_src, const int s_offset, LSPack& lsPack, const int nt) {
  lsPack.mat->reset( 0.0 );
  lsPack.rhs->reset( 0.0 );
  lsPack.sol->reset( 0.0 );

  VComp comp(1,0.0);
  Entry elem(1,0.0);

  // RHS
  for ( int i = 0; i < Na; i++ ) {
    comp = h_res[i];
    lsPack.rhs->set(i,comp);
  }

  // Matrix
  // - mass matrix
  for ( int i = 0; i < Na; i++ ) {
    elem = 1.0;
    lsPack.mat->set(i,i,elem);
  }
  // - source terms - assuming no nearby data is needed!
  for ( int i = 0; i < Na; i++ ) {
    elem = lsPack.mat->get(i,i);
    elem -= dtn * s_src[i].dx( s_offset+i );
    lsPack.mat->set(i,i,elem);
  }
  // - fluxes -
  int sten = 1;
  for ( int i = 0; i < Na; i++ ) {
    int i_min = max(0,i-sten);
    int i_max = min(Na-1,i+sten);
    DType denom = x_ale_c[nt][i] * x_ale_c[nt][i] * dx_ale[nt][i];
    for ( int j = i_min; j <= i_max; j++ ) {
      elem = lsPack.mat->get(i,j);
      elem += (dtn/denom) * (s_flx_d[i+1].dx(s_offset+j)-s_flx_d[i].dx(s_offset+j)); 
      elem += dtn*s_flx_r[i].dx(s_offset+j);
      lsPack.mat->set(i,j,elem);
    }
  }

  // Solve the linear system
  lsPack.linSolve->solve();

  // Fetch and update solution
  for ( int i = 0; i < Na; i++ ) {
    comp = lsPack.sol->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    h_var[i] -= vsol[0];
  }

  return;
}


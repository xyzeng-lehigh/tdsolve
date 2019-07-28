#include "TDSolveFVMMod.h"
#include "Alloc_Utils.h"
#include "TDInput.h"
#include "Alg_Solve.h"
#include "Alg_Limiters.h"
using namespace std;

TDSolveFVMMod::TDSolveFVMMod(TDInput* tdInput_ptr) :
  TDSolveFVM(tdInput_ptr) {

  lim_g     = NULL;
  lim_h     = NULL;
  lim_n     = NULL;
  lim_s     = NULL;
  lim_muscl = NULL;

  flux_g_v  = NULL;
  flux_g_r  = NULL;
  flux_h_v  = NULL;
  flux_h_r  = NULL;
  flux_n_v  = NULL;
  flux_n_r  = NULL;
  flux_n_u  = NULL;

  data_u    = NULL;
  vel_tmp   = NULL;

}

TDSolveFVMMod::~TDSolveFVMMod(void) {
  erase_array( lim_g, 2 );
  erase_array( lim_h, 2 );
  erase_array( lim_n, 2 );
  erase_array( lim_s, 2 );
  erase_array( lim_muscl, 1 );
  erase_array( flux_g_v, 1 );
  erase_array( flux_g_r, 1 );
  erase_array( flux_h_v, 1 );
  erase_array( flux_h_r, 1 );
  erase_array( flux_n_v, 1 ); // or M
  erase_array( flux_n_r, 1 ); // or M
  erase_array( flux_n_u, Ndata ); // or M
  erase_array( data_u,   Ndata );
  erase_array( vel_tmp,  1 );
}

void TDSolveFVMMod::set_ic(void) {
  if ( tdSetup->prob_id >= 0 ) {
    TDSolveFVM::set_ic();
    return;
  }

  // Testing only data_u is prescribed rather than computed
  data_r[0] = tdCoef->R_0;
  switch ( tdSetup->ic ) {
    case 0:
      {
        *data_g[0] = 0.0;
        *data_h[0] = 1.0-tdCoef->n_bc;
        *data_n[0] = tdCoef->n_bc;
      }
      break;
    case 1:
      {
        *data_g[0] = 1.0-tdCoef->n_bc;
        *data_h[0] = 0.0;
        *data_n[0] = tdCoef->n_bc;
      }
      break;
    case 2:
      {
        *data_g[0] = 0.5-0.5*tdCoef->n_bc;
        *data_h[0] = 0.5-0.5*tdCoef->n_bc;
        *data_n[0] = tdCoef->n_bc;
      }
      break;
  }

  upscale_var( *data_g[0], *data_g[0], x_c, Nx, data_r[0]);
  upscale_var( *data_h[0], *data_h[0], x_c, Nx, data_r[0]);
  upscale_var( *data_n[0], *data_n[0], x_c, Nx, data_r[0]);

  calc_total_diff_model(0);

  *data_u[0] = 0.0;
  switch ( tdSetup->prob_id ) {
    case -1:
      {
        DType u_pow = tdSetup->init_u_par;
        if ( u_pow > 0 )
          for ( int i = 1; i <= Nx; i++ )
            (*data_u[0])[i] = -std::pow(i*dx,u_pow);
      }
      break;
    case -2:
      {
        DType v_par = tdSetup->init_u_par;
        for ( int i = 1; i <= Nx; i++ )
          (*data_u[0])[i] = - x_f[i] * data_r[0] * v_par / tdCoef->n_bc;
      }
      break;
    case -3:
      {
        DType v_par = tdSetup->init_u_par;
        for ( int i = 1; i <= Nx; i++ )
          (*data_u[0])[i] =  - v_par * std::sin( 0.5 * PI * x_f[i] ) / tdCoef->n_bc;
      }
      break;
    case -4:
      {
        DType v_par = tdSetup->init_u_par;
        for ( int i = 1; i <= Nx; i++ )
          (*data_u[0])[i] = v_par * std::sin( 4.0 * x_f[i] * data_r[0] );
      }
      break;
  }

  // Set all the right hand sides to zero:
  tdCoef->lambda  = 0.0;
  tdCoef->mu      = 0.0;
  tdCoef->delta_h = 0.0;
  tdCoef->delta_n = 0.0;
  tdCoef->omega   = 0.0;
  tdSetup->logistic = 0;

  this->update_vel( 0, 1, 0.0 );

  return;
}

void TDSolveFVMMod::calc_dt(const int nt) {
  // Compute dt using *data_u[0]
  if ( iter != 0 )
    tsol += dtn;

  if ( data_r[nt] <= tdSetup->Rdead || data_r[nt] >= tdSetup->Rterm )
    reach_term_now = true;

  if ( reach_term_next )
    reach_term_now = true;

  DType dtn_trial = DTYPE_MAX;
  DType dr = (*data_v[nt])[Nx]; // A predictor

  compute_vel_u( nt, tsol );

  if ( tdSetup->scell < 0 ) {
    for ( int i = 0; i < Nx+1; i++ ) {
      DType loc_vel = (*data_v[nt])[i]/data_r[nt];
      DType m_vel = x_f[i]*dr/data_r[nt];
      DType u_vel = (*data_u[0])[i]/data_r[nt];
      DType v_max = max( fabs(loc_vel-m_vel), fabs(loc_vel-m_vel+u_vel) );
      v_max = max( fabs(m_vel), v_max );
      dtn_trial = min( dtn_trial, tdSetup->cfl * dx / v_max );
    }
  }
  else {
    for ( int i = 0; i < Nx+1; i++ ) {
      DType loc_vel = (*data_v[nt])[i]/data_r[nt];
      DType m_vel = x_f[i]*dr/data_r[nt];
      DType u_vel = (*data_u[0])[i]/data_r[nt];
      DType v_max = fabs( loc_vel ) + fabs( m_vel ) + fabs( u_vel );
      dtn_trial = min( dtn_trial, tdSetup->cfl * dx / v_max );
    }
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

void TDSolveFVMMod::print_custom_vel_u( FILE* ot_file, const int nt ){
  if ( tdSetup->silence ) return;
  fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\t%10.20f\n", tsol, data_r[nt], (*data_v[nt])[Nx], -(*data_u[0])[Nx] );
  return;
}

void TDSolveFVMMod::alloc(void) {

  TDSolveFVM::alloc();

  alloc_array_full( lim_g, Nx, 2 );
  alloc_array_full( lim_h, Nx, 2 );
  alloc_array_full( lim_n, Nx, 2 ); // M
  alloc_array_full( lim_s, Nx, 2 ); // Sync
  alloc_array_full( lim_muscl, Nx, 1 );

  alloc_array_full( flux_g_v, Nx+1, 1 );
  alloc_array_full( flux_g_r, Nx+1, 1 );
  alloc_array_full( flux_h_v, Nx+1, 1 );
  alloc_array_full( flux_h_r, Nx+1, 1 );
  alloc_array_full( flux_n_v, Nx+1, 1 ); // or M
  alloc_array_full( flux_n_r, Nx+1, 1 ); // or M
  alloc_array_full( flux_n_u, Nx+1, Ndata ); // or M

  alloc_array_full( data_u, Nx+1, Ndata );
  if ( tdSetup->vCoup )
    alloc_array_full( vel_tmp, Nx+1, 1 );

  return;
}

bool TDSolveFVMMod::update_sol(const int fr, const int to,
    const DType shift_dt_fr, const DType shift_dt_to) {

  // At this point dtn and tsol are known

  if ( tdSetup->vCoup ) {
    // In the sub-iteration, always use v_fr to compute fluxes
    // Do not use copy_data for now, which will overwrite data_r
    *data_g[to] = *data_g[fr];
    *data_h[to] = *data_h[fr];
    *data_n[to] = *data_n[fr];
    *data_v[to] = *data_v[fr];
    data_r[to] = data_r[fr] + (shift_dt_to - shift_dt_fr) * dtn * (*data_v[fr])[Nx];
    int   cur_iter = 0;
    int   max_iter = tdSetup->vIter;
    DType max_tol  = tdSetup->vTol;
    DType cur_res;
    while ( max_iter < 0 || cur_iter < max_iter ) {
      // Use *data_v[to], data_r[to], and data_rp[0] to compute
      // the cells implicitly
      switch ( tdSetup->scell ) {
        case 0: // Implicit first-order upwind
          update_cells_imp_upw( fr, to, tsol+dtn*shift_dt_to);
          break;
        case 2: // Implicit second-order MUSCL 
          update_cells_imp_muscl( fr, to, tsol+dtn*shift_dt_to);
          break;
        case -1: // Implicit first-order upwind, origin
          update_cells_orig_imp_upw( fr, to, tsol+dtn*shift_dt_to);
          break;
        case -2: // Implicit second-order MUSCL, origin
          update_cells_orig_imp_muscl( fr, to, tsol+dtn*shift_dt_to);
          break;
        default:
          std::cout << "Warning: The flux scheme for cell not defined, "
            << "using implicit upwind instead\n" << std::endl;
          update_cells_imp_upw( fr, to, tsol+dtn*shift_dt_to);
          break;
      }

      // Then to work on the chems
      if ( tdSetup->prob_id == 0 ) // PDGF
        *data_a[to] = *data_a[fr];
      if ( tdSetup->prob_id == 1 ) { // OPN
        *data_a[to]  = *data_a[fr];
        *data_p1[to] = *data_p1[fr];
        *data_p2[to] = *data_p2[fr];
      }
      update_chem_imp( fr, to, shift_dt_fr, shift_dt_to );

      *vel_tmp[0] = *data_v[to];
      update_vel( to, 0, shift_dt_to, true ); // Update data_v[to] no updating on radius
      update_rad( fr, to, to );           // Update the radius data_r[to] and data_rp[0] 
      *vel_tmp[0] -= *data_v[to];
      cur_res = (*vel_tmp[0]).norm(-1);
      if ( (*vel_tmp[0]).norm(-1) < max_tol )
        break;

      cur_iter++;
    }

    if ( cur_iter == max_iter && cur_res > max_tol ) {
      cout << "Coupling iteration: " << cur_iter << "; Residue: " << cur_res << "; Will reduce dt." << std::endl;
      return false;
    }

    if ( tdSetup->diagn_vcoup )
      cout << "Coupling iteration: " << cur_iter << "; Residue: " << cur_res << std::endl;
    if ( tdSetup->prob_id < 0 )
      calc_total_diff_model( to );
    return true;
  }

  // Will also compute fluxes for the current time step
  // and update_rad
  update_vel( fr, to, shift_dt_fr );

  // The fluxes are already computed in previous iteration
  if ( tdSetup->scell >= -50 ) // Non-auto-inc methods
    update_cells_common( fr, to );
  else
    update_cells_common_auto( fr, to );

  if ( tdSetup->prob_id == 1 ) // OPN
    update_cells_opn_add( fr, to );

  // Simple model
  if ( tdSetup->prob_id < 0 )
    calc_total_diff_model( to );

  update_chem_imp( fr, to, shift_dt_fr, shift_dt_to );
  
  update_vel( to, to, shift_dt_to, true ); // Update data_v[to] no updating on radius

  return true;
}

void TDSolveFVMMod::compute_rhs(const DType alpha, const int nt) {
  // At this moment, flux_n_u[nt] and slope_*[nt] is computed
  compute_cells_fluxes_v( nt, nt ); // Will compute flux_*_v[0]
  compute_cells_fluxes_r( nt, nt ); // Will compute flux_*_r[0] and flux_*[0]
  compute_rhs_common( nt );     // Stores the residuals to resd_*[nt]
  return;
}

void TDSolveFVMMod::compute_rhs_common(const int nt) {
  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

  DType dror = data_rp[nt]/data_r[nt];

  // Source terms
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[i],p) ); 
    (*resd_g[nt])[i] = loc_lambda * (*data_g[nt])[i] - loc_mu * (*data_g[nt])[i] - dror * (*data_g[nt])[i];
    (*resd_h[nt])[i] = loc_mu * (*data_g[nt])[i] - delta_h * (*data_h[nt])[i] - dror * (*data_h[nt])[i];
    (*resd_n[nt])[i] = - delta_n * (*data_n[nt])[i] - dror * (*data_n[nt])[i]; 
  }

  for ( int i = 0; i < Nx; i++ ) {
    (*resd_g[nt])[i] -= ( (*flux_g[0])[i+1] - (*flux_g[0])[i] ) / dx;
    (*resd_h[nt])[i] -= ( (*flux_h[0])[i+1] - (*flux_h[0])[i] ) / dx;
    (*resd_n[nt])[i] -= ( (*flux_n[0])[i+1] - (*flux_n[0])[i] ) / dx;
  }

  return;
}

void TDSolveFVMMod::update_sol(const int rkns, RKButcher& rkButcher, const int fr, const int to) {
  if ( tdSetup->vCoup ) {
    // Implicit update
    switch ( tdSetup->scell ) {
      case 0: // Implicit first-order upwind
        update_sol_radau_imp_upw( rkns, rkButcher, fr, to );
        break;
      case 2: // Implicit second-order MUSCL
        update_sol_radau_imp_muscl( rkns, rkButcher, fr, to );
        break;
      default: // Implicit first-order upwind
        update_sol_radau_imp_upw( rkns, rkButcher, fr, to );
        break;
    }
  }
  else {
    // Explicit update -- should only use resd_*
    update_sol_radau_exp( rkns, rkButcher, fr, to );
  }

  return;
}

void TDSolveFVMMod::update_vel(const DType alpha, const int nt) {
  // Compute the velocity u data_u[nt] using the gradient
  compute_vel_u( nt, tsol+alpha*dtn, nt );
  // Compute fluxes - flux_n_u[nt]
  // and (if MUSCL) - slope_*[nt]
  compute_cells_fluxes_u( nt, tsol+alpha*dtn, nt );
  switch ( tdSetup->prob_id ) {
    case 0:  // PDGF
    case -1: // Test
    case -2: // Test
    case -3: // Test
    case -4: // Test
      update_vel_int( nt, nt ); // Will use flux_n_u[0]
      break;
    case 1: // OPN
      update_vel_int_opn( nt, nt ); // Will use flux_n_u[0]
      break;
  }

  return;
}

void TDSolveFVMMod::copy_aux(const int fr, const int to) {
  *data_u[to]   = *data_u[fr];
  *flux_n_u[to] = *flux_n_u[fr];
  *slope_g[to]  = *slope_g[fr];
  *slope_h[to]  = *slope_h[fr];
  *slope_n[to]  = *slope_n[fr];
}

void TDSolveFVMMod::calc_total_diff_model(const int nt) {
  downscale_var( *data_g[nt], *data_g[nt], x_c, Nx, data_r[nt]);
  downscale_var( *data_h[nt], *data_h[nt], x_c, Nx, data_r[nt]);
  downscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt]);

  for ( int i = 0; i < Nx; i++ )
    (*data_a[nt])[i] = 1.0 - (*data_g[nt])[i] - (*data_h[nt])[i] - (*data_n[nt])[i];
  
  upscale_var( *data_g[nt], *data_g[nt], x_c, Nx, data_r[nt] );
  upscale_var( *data_h[nt], *data_h[nt], x_c, Nx, data_r[nt] );
  upscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );
  upscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );
  return;
}

void TDSolveFVMMod::update_vel(const int fr, const int to, const DType shift_dt, const bool last_iter) {

  // Compute the velocity u using the gradient
  compute_vel_u( fr, tsol+shift_dt*dtn );

  if ( tdSetup->scell >= 0 ) { // Segregated fluxes
    // The fluxes that does not need V
    compute_cells_fluxes_u( fr, tsol+shift_dt*dtn );

    switch ( tdSetup->prob_id ) {
      case 0:  // PDGF
      case -1: // Test
      case -2: // Test
      case -3: // Test
      case -4: // Test
        update_vel_int(fr);
        break;
      case 1: // OPN
        update_vel_int_opn(fr);
        break;
    }
  }
  else { // Original fluxes
    switch ( tdSetup->prob_id ) {
      case 0:  // PDGF
        update_vel_int_orig( fr );
        break;
      case 1:  // OPN
        update_vel_int_opn_orig( fr );
        break;
    }
  }

  if ( last_iter ) return;

  // Collocation between fr and to should be dtn
  if ( fr != to )
    update_rad( fr, to );

  if ( tdSetup->scell >= 0 ) { // Segregated fluxes
    // The fluxes that need V
    compute_cells_fluxes_v( fr );

    // The fluxes associated with R'
    compute_cells_fluxes_r( fr );
  }
  else {
    compute_cells_fluxes_orig( fr, tsol+shift_dt*dtn );
  }

  return;
}

void TDSolveFVMMod::compute_vel_u(const int nt, const DType tsol_loc, const int u_nt) {
  switch ( tdSetup->prob_id ) {
    case -1: // Test 1, u is only a function of x
      break;
    case -2: // Test 2, space homogeneous solutions
      {
        DType v_par = tdSetup->init_u_par;
        (*data_u[u_nt])[0] = 0.0;
        for ( int i = 1; i <= Nx; i++ )
          (*data_u[u_nt])[i] = - x_f[i] * data_r[nt] * v_par / ( 1.0 - std::exp(-3.0*v_par*tsol_loc) * (1.0-tdCoef->n_bc));
      }
      break;
    case -3: // Test 3, linear growth
      {
        DType v_par = tdSetup->init_u_par;
        (*data_u[u_nt])[0] = 0.0;
        for ( int i = 1; i <= Nx; i++ )
          (*data_u[u_nt])[i] = -  v_par * std::sin( 0.5 * PI * x_f[i] ) / tdCoef->n_bc;
      }
      break;
    case -4: // Test 4, non-monotone radius change
      {
        DType v_par = tdSetup->init_u_par;
        (*data_u[u_nt])[0] = 0.0;
        for ( int i = 1; i <= Nx; i++ )
          (*data_u[u_nt])[i] = v_par * std::sin( 4.0 * x_f[i] * data_r[nt] * ( 1.0 + tsol_loc ) );
      }
      break;
    case 0: // PDGF
    case 1: // OPN
      // u = (alpha/R^2)(dA/dx)
      compute_vel_u_chem( nt, u_nt );
      break;
  }
  return;
}

void TDSolveFVMMod::compute_vel_u_chem(const int nt, const int u_nt) {
  *data_u[u_nt] = 0.0;

  DType alpha = tdCoef->alpha;
  DType aor  = alpha/data_r[nt];

  downscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );

  for ( int i = 1; i <= Nx; i++ )
    (*data_u[u_nt])[i] = aor * (((*data_a[nt])[i]-(*data_a[nt])[i-1]) / dx);

  upscale_var( *data_a[nt], *data_a[nt], x_a_c, Na, data_r[nt] );

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_u(const int nt, const DType tsol_loc, const int u_nt) {
  switch ( tdSetup->scell ) {
    case 0: // Explicit first-order upwind
      compute_cells_fluxes_u_upw( nt, tsol_loc, u_nt );
      break;
    case 2: // Explicit second-order MUSCL
      compute_cells_fluxes_u_muscl( nt, tsol_loc, u_nt );
      break;
    default:
      std::cout << "Warning: The flux scheme for cell not defined, "
        << "using upwind instead\n" << std::endl;
      compute_cells_fluxes_u_upw( nt, tsol_loc, u_nt );
      break;
  }

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_v(const int nt, const int aux_nt) {
  switch ( tdSetup->scell ) {
    case 0: // Explicit first-order upwind
      compute_cells_fluxes_v_upw( nt );
      break;
    case 2: // Explicit second-order MUSCL
      compute_cells_fluxes_v_muscl( nt, aux_nt );
      break;
    default:
      std::cout << "Warning: The flux scheme for cell not defined, "
        << "using upwind instead\n" << std::endl;
      compute_cells_fluxes_v_upw( nt );
      break;
  }

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_u_upw(const int nt, const DType tsol_loc, const int u_nt) {

  downscale_var(*data_n[nt],*data_n[nt],x_c,Nx,data_r[nt]);

  DType n_bc = calc_n_bc( tsol_loc );

  // Fluxes
  { // Left boundary
    (*flux_n_u[u_nt])[0] = 0.0;
  }
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_r  = x_f[i] * data_r[nt];
    DType loc_r2 = loc_r * loc_r;

    DType loc_vel_u = (*data_u[u_nt])[i]/data_r[nt];
    int i_idx_u = (loc_vel_u>=0.0) ? i-1 : i;
    (*flux_n_u[u_nt])[i] = loc_r2 * loc_vel_u * (*data_n[nt])[i_idx_u];
  }
  { // Right boundary
    DType loc_vel_u = (*data_u[u_nt])[Nx]/data_r[nt];
    DType loc_r = x_f[Nx]*data_r[nt];
    (*flux_n_u[u_nt])[Nx] = loc_r * loc_r * loc_vel_u * 
      ((loc_vel_u>=0.0) ? (*data_n[nt])[Nx-1] : n_bc);
  }

  upscale_var(*data_n[nt],*data_n[nt],x_c,Nx,data_r[nt]);

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_u_muscl(const int nt, const DType tsol_loc, const int u_nt) {

  downscale_var(*data_g[nt],*data_g[nt],x_c,Nx,data_r[nt]);
  downscale_var(*data_h[nt],*data_h[nt],x_c,Nx,data_r[nt]);
  downscale_var(*data_n[nt],*data_n[nt],x_c,Nx,data_r[nt]);

  DType n_bc = calc_n_bc( tsol_loc );

  // Compute the limiters based on the primitive variables
  compute_flux_muscl_lims( nt );

  // Slopes
  { // Left boundary
    (*slope_g[u_nt])[0] = 0.0;
    (*slope_h[u_nt])[0] = 0.0;
    (*slope_n[u_nt])[0] = 0.0;
  }
  for ( int i = 1; i < Nx-1; i++ ) {
    (*slope_g[u_nt])[i] = (*lim_muscl[0])[i] * 
      ((*data_g[nt])[i+1]-(*data_g[nt])[i])/dx;
    (*slope_h[u_nt])[i] = (*lim_muscl[0])[i] * 
      ((*data_h[nt])[i+1]-(*data_h[nt])[i])/dx;
    (*slope_n[u_nt])[i] = (*lim_muscl[0])[i] * 
      ((*data_n[nt])[i+1]-(*data_n[nt])[i])/dx;
  }
  { // Right boundary
    (*slope_g[u_nt])[Nx-1] = 0.0;
    (*slope_h[u_nt])[Nx-1] = 0.0;
    (*slope_n[u_nt])[Nx-1] = 0.0;
  }

  // Fluxes
  { // Left boundary
    (*flux_n_u[u_nt])[0] = 0.0;
  }
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_r  = x_f[i] * data_r[nt];
    DType loc_r2 = loc_r * loc_r;

    DType loc_vel_u = (*data_u[u_nt])[i]/data_r[nt];
    if ( loc_vel_u >= 0.0 )
      (*flux_n_u[u_nt])[i] = loc_r2 * loc_vel_u * (
          (*data_n[nt])[i-1]+0.5*dx*(*slope_n[u_nt])[i-1] );
    else
      (*flux_n_u[u_nt])[i] = loc_r2 * loc_vel_u * (
          (*data_n[nt])[i]-0.5*dx*(*slope_n[u_nt])[i] );
  }
  { // Right boundary
    DType loc_vel_u = (*data_u[u_nt])[Nx]/data_r[nt];
    DType loc_r = x_f[Nx]*data_r[nt];
    (*flux_n_u[u_nt])[Nx] = loc_r * loc_r * loc_vel_u * 
      ((loc_vel_u>=0.0) ? (*data_n[nt])[Nx-1]+0.5*dx*(*slope_n[u_nt])[Nx-1] : n_bc );
  }

  upscale_var(*data_g[nt],*data_g[nt],x_c,Nx,data_r[nt]);
  upscale_var(*data_h[nt],*data_h[nt],x_c,Nx,data_r[nt]);
  upscale_var(*data_n[nt],*data_n[nt],x_c,Nx,data_r[nt]);

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_v_upw(const int nt) {

  downscale_var(*data_g[nt],*data_g[nt],x_c,Nx,data_r[nt]);
  downscale_var(*data_h[nt],*data_h[nt],x_c,Nx,data_r[nt]);
  downscale_var(*data_n[nt],*data_n[nt],x_c,Nx,data_r[nt]);

  // Fluxes
  { // Left boundary
    (*flux_g_v[0])[0] = 0.0;
    (*flux_h_v[0])[0] = 0.0;
    (*flux_n_v[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_r  = x_f[i] * data_r[nt];
    DType loc_r2 = loc_r * loc_r;

    DType loc_vel_v = (*data_v[nt])[i]/data_r[nt];
    int i_idx_v = (loc_vel_v>=0.0) ? i-1 : i;
    (*flux_g_v[0])[i] = loc_r2 * loc_vel_v * (*data_g[nt])[i_idx_v];
    (*flux_h_v[0])[i] = loc_r2 * loc_vel_v * (*data_h[nt])[i_idx_v];
    (*flux_n_v[0])[i] = loc_r2 * loc_vel_v * (*data_n[nt])[i_idx_v];
  }
  { // Right boundary
    (*flux_g_v[0])[Nx] = 0.0;
    (*flux_h_v[0])[Nx] = 0.0;
    (*flux_n_v[0])[Nx] = 0.0;
  }

  upscale_var(*data_g[nt],*data_g[nt],x_c,Nx,data_r[nt]);
  upscale_var(*data_h[nt],*data_h[nt],x_c,Nx,data_r[nt]);
  upscale_var(*data_n[nt],*data_n[nt],x_c,Nx,data_r[nt]);

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_v_muscl(const int nt, const int slp_nt) {

  downscale_var(*data_g[nt],*data_g[nt],x_c,Nx,data_r[nt]);
  downscale_var(*data_h[nt],*data_h[nt],x_c,Nx,data_r[nt]);
  downscale_var(*data_n[nt],*data_n[nt],x_c,Nx,data_r[nt]);

  // Fluxes
  { // Left boundary
    (*flux_g_v[0])[0] = 0.0;
    (*flux_h_v[0])[0] = 0.0;
    (*flux_n_v[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_r  = x_f[i] * data_r[nt];
    DType loc_r2 = loc_r * loc_r;

    DType loc_vel_v = (*data_v[nt])[i]/data_r[nt];
    if ( loc_vel_v >= 0.0 ) {
      (*flux_g_v[0])[i] = loc_r2 * loc_vel_v * (
          (*data_g[nt])[i-1]+0.5*dx*(*slope_g[slp_nt])[i-1] );
      (*flux_h_v[0])[i] = loc_r2 * loc_vel_v * (
          (*data_h[nt])[i-1]+0.5*dx*(*slope_h[slp_nt])[i-1] );
      (*flux_n_v[0])[i] = loc_r2 * loc_vel_v * (
          (*data_n[nt])[i-1]+0.5*dx*(*slope_n[slp_nt])[i-1] );
    }
    else {
      (*flux_g_v[0])[i] = loc_r2 * loc_vel_v * (
          (*data_g[nt])[i]-0.5*dx*(*slope_g[slp_nt])[i] );
      (*flux_h_v[0])[i] = loc_r2 * loc_vel_v * (
          (*data_h[nt])[i]-0.5*dx*(*slope_h[slp_nt])[i] );
      (*flux_n_v[0])[i] = loc_r2 * loc_vel_v * (
          (*data_n[nt])[i]-0.5*dx*(*slope_n[slp_nt])[i] );
    }
  }
  { // Right boundary
    (*flux_g_v[0])[Nx] = 0.0;
    (*flux_h_v[0])[Nx] = 0.0;
    (*flux_n_v[0])[Nx] = 0.0;
  }

  upscale_var(*data_g[nt],*data_g[nt],x_c,Nx,data_r[nt]);
  upscale_var(*data_h[nt],*data_h[nt],x_c,Nx,data_r[nt]);
  upscale_var(*data_n[nt],*data_n[nt],x_c,Nx,data_r[nt]);

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_r(const int nt, const int rp_nt) {
  DType dror = data_rp[rp_nt]/data_r[nt];

  const DType sot = 7.0/12.0;
  const DType oot = 1.0/12.0;

  // Compute the limiters and store them in lim_?, as well as 
  // unlimited interpolated values, and store them in flux_?
  compute_flux_r_ppm_lims( nt );

  // x^2R^2N -> x^3R^2N, etc.
  upscale_var(*data_g[nt],*data_g[nt],x_c,Nx,1.0,1);
  upscale_var(*data_h[nt],*data_h[nt],x_c,Nx,1.0,1);
  upscale_var(*data_n[nt],*data_n[nt],x_c,Nx,1.0,1);

  // Step 3 -- compute the fluxes
  // Left boundary
  {
    (*flux_g_r[0])[0] = 0.0;
    (*flux_h_r[0])[0] = 0.0;
    (*flux_n_r[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_vel_r = - dror;
    DType lim_l = (*lim_s[1])[i-1];
    DType lim_r = (*lim_s[0])[i];
    //DType lim_l = 1.0;
    //DType lim_r = 1.0;
    if ( loc_vel_r >= 0.0 ) {
      (*flux_g_r[0])[i] = loc_vel_r * ( (1.0-lim_l) * (*data_g[nt])[i-1] + lim_l*(*flux_g[0])[i] );
      (*flux_h_r[0])[i] = loc_vel_r * ( (1.0-lim_l) * (*data_h[nt])[i-1] + lim_l*(*flux_h[0])[i] );
      (*flux_n_r[0])[i] = loc_vel_r * ( (1.0-lim_l) * (*data_n[nt])[i-1] + lim_l*(*flux_n[0])[i] );
    }
    else {
      (*flux_g_r[0])[i] = loc_vel_r * ( (1.0-lim_r) * (*data_g[nt])[i] + lim_r*(*flux_g[0])[i] );
      (*flux_h_r[0])[i] = loc_vel_r * ( (1.0-lim_r) * (*data_h[nt])[i] + lim_r*(*flux_h[0])[i] );
      (*flux_n_r[0])[i] = loc_vel_r * ( (1.0-lim_r) * (*data_n[nt])[i] + lim_r*(*flux_n[0])[i] );
    }
  }
  // Right boundary
  {
    (*flux_g_r[0])[Nx] = 0.0;
    (*flux_h_r[0])[Nx] = 0.0;
    (*flux_n_r[0])[Nx] = 0.0;
  }

  // x^3R^2N -> x^2R^2N, etc.
  downscale_var(*data_g[nt],*data_g[nt],x_c,Nx,1.0,1);
  downscale_var(*data_h[nt],*data_h[nt],x_c,Nx,1.0,1);
  downscale_var(*data_n[nt],*data_n[nt],x_c,Nx,1.0,1);

  for ( int i = 0; i < Nx+1; i++ ) {
    (*flux_g[0])[i] = (*flux_g_v[0])[i]+(*flux_g_r[0])[i];
    (*flux_h[0])[i] = (*flux_h_v[0])[i]+(*flux_h_r[0])[i];
    (*flux_n[0])[i] = (*flux_n_v[0])[i]+(*flux_n_r[0])[i]+(*flux_n_u[rp_nt])[i];
  }

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_orig(const int nt, const DType tsol_loc) {
  switch ( tdSetup->scell ) {
    case -1: // Explicit first-order upwind
      compute_cells_fluxes_orig_upw( nt, tsol_loc );
      break;
    case -2: // Explicit second-order MUSCL
      compute_cells_fluxes_orig_muscl( nt, tsol_loc );
      break;
    case -101: // Explicit first-order upwind, auto g = 1 - h - n
      compute_cells_fluxes_orig_upw_auto_g( nt, tsol_loc );
      break;
    case -102: // Explicit first-order upwind, auto h = 1 - g - n
      compute_cells_fluxes_orig_upw_auto_h( nt, tsol_loc );
      break;
    case -103: // Explicit first-order upwind, auto n = 1 - g - h
      compute_cells_fluxes_orig_upw_auto_n( nt, tsol_loc );
      break;
    default:
      std::cout << "Warning: The flux scheme for cell not defined, "
        << "using upwind instead\n" << std::endl;
      compute_cells_fluxes_orig_upw( nt, tsol_loc );
      break;
  }

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_orig_upw(const int nt, const DType tsol_loc) {

  DType dror = data_rp[0]/data_r[nt];

  DType n_bc = calc_n_bc( tsol_loc );

  // Fluxes
  { // Left boundary
    (*flux_g[0])[0] = 0.0;
    (*flux_h[0])[0] = 0.0;
    (*flux_n[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_vel = (*data_v[nt])[i]/data_r[nt] - x_f[i]*dror;
    (*flux_g[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_g[nt])[i-1] : (*data_g[nt])[i] );
    (*flux_h[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_h[nt])[i-1] : (*data_h[nt])[i] );
    loc_vel += (*data_u[0])[i]/data_r[nt];
    (*flux_n[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_n[nt])[i-1] : (*data_n[nt])[i] );
  }
  { // Right boundary
    (*flux_g[0])[Nx] = 0.0;
    (*flux_h[0])[Nx] = 0.0;
    DType loc_vel = (*data_u[0])[Nx]/data_r[nt];
    (*flux_n[0])[Nx] = loc_vel * (
        (loc_vel>=0.0) ? (*data_n[nt])[Nx-1] : 
        n_bc * (x_f[Nx]+0.5*dx)*(x_f[Nx]+0.5*dx)*data_r[nt]*data_r[nt] );
  }

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_orig_upw_auto_g(const int nt, const DType tsol_loc) {

  DType dror = data_rp[0]/data_r[nt];

  DType n_bc = calc_n_bc( tsol_loc );

  // Fluxes
  { // Left boundary
    (*flux_h[0])[0] = 0.0;
    (*flux_n[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_vel = (*data_v[nt])[i]/data_r[nt] - x_f[i]*dror;
    (*flux_h[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_h[nt])[i-1] : (*data_h[nt])[i] );
    loc_vel += (*data_u[0])[i]/data_r[nt];
    (*flux_n[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_n[nt])[i-1] : (*data_n[nt])[i] );
  }
  { // Right boundary
    (*flux_h[0])[Nx] = 0.0;
    DType loc_vel = (*data_u[0])[Nx]/data_r[nt];
    (*flux_n[0])[Nx] = loc_vel * (
        (loc_vel>=0.0) ? (*data_n[nt])[Nx-1] : 
        n_bc * (x_f[Nx]+0.5*dx)*(x_f[Nx]+0.5*dx)*data_r[nt]*data_r[nt] );
  }

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_orig_upw_auto_h(const int nt, const DType tsol_loc) {

  DType dror = data_rp[0]/data_r[nt];

  DType n_bc = calc_n_bc( tsol_loc );

  // Fluxes
  { // Left boundary
    (*flux_g[0])[0] = 0.0;
    (*flux_n[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_vel = (*data_v[nt])[i]/data_r[nt] - x_f[i]*dror;
    (*flux_g[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_g[nt])[i-1] : (*data_g[nt])[i] );
    loc_vel += (*data_u[0])[i]/data_r[nt];
    (*flux_n[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_n[nt])[i-1] : (*data_n[nt])[i] );
  }
  { // Right boundary
    (*flux_g[0])[Nx] = 0.0;
    DType loc_vel = (*data_u[0])[Nx]/data_r[nt];
    (*flux_n[0])[Nx] = loc_vel * (
        (loc_vel>=0.0) ? (*data_n[nt])[Nx-1] : 
        n_bc * (x_f[Nx]+0.5*dx)*(x_f[Nx]+0.5*dx)*data_r[nt]*data_r[nt] );
  }

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_orig_upw_auto_n(const int nt, const DType tsol_loc) {

  DType dror = data_rp[0]/data_r[nt];

  // Fluxes
  { // Left boundary
    (*flux_g[0])[0] = 0.0;
    (*flux_h[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_vel = (*data_v[nt])[i]/data_r[nt] - x_f[i]*dror;
    (*flux_g[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_g[nt])[i-1] : (*data_g[nt])[i] );
    (*flux_h[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_h[nt])[i-1] : (*data_h[nt])[i] );
  }
  { // Right boundary
    (*flux_g[0])[Nx] = 0.0;
    (*flux_h[0])[Nx] = 0.0;
  }

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_orig_muscl(const int nt, const DType tsol_loc) {

  DType dror = data_rp[0]/data_r[nt];

  DType n_bc = calc_n_bc( tsol_loc );

  // Slopes
  { // Left boundary
    (*slope_g[0])[0] = 0.0;
    (*slope_h[0])[0] = 0.0;
    (*slope_n[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx-1; i++ ) {
    (*slope_g[0])[i] = Alg_LimSlope_Minmod( (*data_g[nt])[i], (*data_g[nt])[i-1], (*data_g[nt])[i+1], dx );
    (*slope_h[0])[i] = Alg_LimSlope_Minmod( (*data_h[nt])[i], (*data_h[nt])[i-1], (*data_h[nt])[i+1], dx );
    (*slope_n[0])[i] = Alg_LimSlope_Minmod( (*data_n[nt])[i], (*data_n[nt])[i-1], (*data_n[nt])[i+1], dx );
  }
  { // Right boundary
    (*slope_g[0])[Nx-1] = 0.0;
    (*slope_h[0])[Nx-1] = 0.0;
    (*slope_n[0])[Nx-1] = 0.0;
  }

  // Fluxes
  { // Left boundary
    (*flux_g[0])[0] = 0.0;
    (*flux_h[0])[0] = 0.0;
    (*flux_n[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_vel = (*data_v[nt])[i]/data_r[nt] - x_f[i]*dror;
    (*flux_g[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_g[nt])[i-1]+0.5*dx*(*slope_g[0])[i-1] : 
                         (*data_g[nt])[i]-0.5*dx*(*slope_g[0])[i] );
    (*flux_h[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_h[nt])[i-1]+0.5*dx*(*slope_h[0])[i-1] : 
                         (*data_h[nt])[i]-0.5*dx*(*slope_h[0])[i] );
    loc_vel += (*data_u[0])[i]/data_r[nt];
    (*flux_n[0])[i] = loc_vel * (
        (loc_vel>=0.0) ? (*data_n[nt])[i-1]+0.5*dx*(*slope_n[0])[i-1] : 
                         (*data_n[nt])[i]-0.5*dx*(*slope_n[0])[i] );
  }
  { // Right boundary
    (*flux_g[0])[Nx] = 0.0;
    (*flux_h[0])[Nx] = 0.0;
    DType loc_vel = (*data_u[0])[Nx]/data_r[nt];
    (*flux_n[0])[Nx] = loc_vel * (
        (loc_vel>=0.0) ? (*data_n[nt])[Nx-1]+0.5*dx*(*slope_n[0])[Nx-1] : 
        n_bc * x_f[Nx] * x_f[Nx] * data_r[nt] * data_r[nt] ) ;
  }

  return;
}

void TDSolveFVMMod::update_cells_orig_imp_upw(const int fr, const int to, const DType tsol_loc) {
  DType lambda = tdCoef->lambda;
  DType q      = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

  DType dror = data_rp[0]/data_r[to];

  DType n_bc = calc_n_bc( tsol_loc );

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
    DType loc_v_r = (*data_v[to])[1]/data_r[to] - x_f[1] * dror;
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
    DType loc_v_r = (*data_v[to])[i+1]/data_r[to] - x_f[i+1] * dror;
    DType loc_v_l = (*data_v[to])[i]/data_r[to] - x_f[i] * dror;
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
    DType loc_v_l = (*data_v[to])[Nx-1]/data_r[to] - x_f[Nx-1] * dror;
    
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
    DType loc_v_r = (*data_v[to])[1]/data_r[to] - x_f[1] * dror;
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
    DType loc_v_r = (*data_v[to])[i+1]/data_r[to] - x_f[i+1] * dror;
    DType loc_v_l = (*data_v[to])[i]/data_r[to] - x_f[i] * dror;
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
    DType loc_v_l = (*data_v[to])[Nx-1]/data_r[to] - x_f[Nx-1] * dror;
    
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

  // RHS
  for ( int i = 0; i < Nx-1; i++ ) {
    comp = (*data_n[fr])[i];
    rhs->set(i,comp);
  }
  // Right boundary
  {
    DType loc_v_r = (*data_u[0])[Nx]/data_r[to];
    
    DType coef_adv_r = 0.0;

    if ( loc_v_r < 0.0 )
      coef_adv_r += dtn*loc_v_r/dx;

    comp = (*data_n[fr])[Nx-1] - coef_adv_r * n_bc * 
      (x_f[Nx]+0.5*dx) * (x_f[Nx]+0.5*dx) * data_r[to] * data_r[to];
    rhs->set(Nx-1,comp);
  }

  // Matrix
  // Left boundary
  {
    DType loc_v_r = (*data_v[to])[1]/data_r[to] - x_f[1] * dror + (*data_u[0])[1]/data_r[to];
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
    DType loc_v_r = (*data_v[to])[i+1]/data_r[to] - x_f[i+1] * dror + (*data_u[0])[i+1]/data_r[to];
    DType loc_v_l = (*data_v[to])[i]/data_r[to] - x_f[i] * dror + (*data_u[0])[i]/data_r[to];
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
    DType loc_v_l = (*data_v[to])[Nx-1]/data_r[to] - x_f[Nx-1] * dror + (*data_u[0])[Nx-1]/data_r[to];
    DType loc_v_r = (*data_u[0])[Nx]/data_r[to];
    
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

void TDSolveFVMMod::update_cells_orig_imp_muscl(const int fr, const int to, const DType tsol_loc) {
  DType lambda = tdCoef->lambda;
  DType q      = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

  DType dror = data_rp[0]/data_r[to];

  DType n_bc = calc_n_bc( tsol_loc );

  compute_flux_muscl_lims( to );

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
  // Source term
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[i],p) ); 
    elem = 1.0 - dtn * ( loc_lambda - loc_mu - dror );
    mat->set(i,i,elem);
  }

  // Temporary nodal velocity for data_g and data_h
  for ( int i = 0; i < Nx; i++ )
    (*vel_tmp[0])[i] = (*data_v[to])[i]/data_r[to] - x_f[i] * dror;

  update_mat_orig_muscl_v( to, vel_tmp, 0, lim_g, false );

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
  for ( int i = 0; i < Nx; i++ ) {
    elem = 1.0 - dtn * ( -delta_h - dror );
    mat->set(i,i,elem);
  }
 
  // No need to modify vel_tmp
  update_mat_orig_muscl_v( to, vel_tmp, 0, lim_h, false );

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

  // RHS
  for ( int i = 0; i < Nx-1; i++ ) {
    comp = (*data_n[fr])[i];
    rhs->set(i,comp);
  }
  // Right boundary
  {
    DType loc_u = (*data_u[0])[Nx]/data_r[to];
    DType coef_u_r = 0.0;
    if ( loc_u < 0.0 )
      coef_u_r += (dtn/dx) * loc_u;
    comp = (*data_n[fr])[Nx-1] - coef_u_r * n_bc * x_f[Nx] * x_f[Nx] * data_r[to] * data_r[to];
    rhs->set(Nx-1,comp);
  }

  // Matrix
  for ( int i = 0; i < Nx; i++ ) {
    elem = 1.0 - dtn * ( - delta_n - dror );
    mat->set(i,i,elem);
  }

  // Nodal velocity
  for ( int i = 0; i < Nx; i++ )
    (*vel_tmp[0])[i] = (*data_v[to])[i]/data_r[to] - x_f[i] * dror + (*data_u[0])[i]/data_r[to];

  update_mat_orig_muscl_v( to, vel_tmp, 0, lim_n, true );

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

void TDSolveFVMMod::compute_flux_muscl_lims(const int nt){

  // Limiters for the velocity associated with u and with v
  { // Left boundary
    (*lim_g[0])[0] = 0.0;
    (*lim_h[0])[0] = 0.0;
    (*lim_n[0])[0] = 0.0;
    (*lim_muscl[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx-1; i++ ) {
    DType loc_lim_g = LimFun_Minmod(
        (*data_g[nt])[i] - (*data_g[nt])[i-1], (*data_g[nt])[i+1] - (*data_g[nt])[i] );
    DType loc_lim_h = LimFun_Minmod(
        (*data_h[nt])[i] - (*data_h[nt])[i-1], (*data_h[nt])[i+1] - (*data_h[nt])[i] );
    DType loc_lim_n = LimFun_Minmod(
        (*data_n[nt])[i] - (*data_n[nt])[i-1], (*data_n[nt])[i+1] - (*data_n[nt])[i] );
    DType loc_lim_theta = 1.0;
    if ( tdSetup->slim_tot ) {
      DType loc_theta[3];
      loc_theta[0] = (*data_g[nt])[i-1] + (*data_h[nt])[i-1] + (*data_n[nt])[i-1];
      loc_theta[1] = (*data_g[nt])[i] + (*data_h[nt])[i] + (*data_n[nt])[i];
      loc_theta[2] = (*data_g[nt])[i+1] + (*data_h[nt])[i+1] + (*data_n[nt])[i+1];
      loc_lim_theta = LimFun_Minmod(loc_theta[1]-loc_theta[0], loc_theta[2]-loc_theta[1]);
    }

    (*lim_g[0])[i] = loc_lim_g;
    (*lim_h[0])[i] = loc_lim_h;
    (*lim_n[0])[i] = loc_lim_n;
    (*lim_muscl[0])[i] = min( min(loc_lim_g,loc_lim_h), min(loc_lim_n,loc_lim_theta) );
  }
  { // Right boundary
    (*lim_g[0])[Nx-1] = 0.0;
    (*lim_h[0])[Nx-1] = 0.0;
    (*lim_n[0])[Nx-1] = 0.0;
    (*lim_muscl[0])[Nx-1] = 0.0;
  }

  return;
}

void TDSolveFVMMod::compute_flux_r_ppm_lims(const int nt) {

  const DType sot = 7.0/12.0;
  const DType oot = 1.0/12.0;

  // x^2R^2N -> x^3R^2N, etc.
  upscale_var(*data_g[nt],*data_g[nt],x_c,Nx,1.0,1);
  upscale_var(*data_h[nt],*data_h[nt],x_c,Nx,1.0,1);
  upscale_var(*data_n[nt],*data_n[nt],x_c,Nx,1.0,1);

  // Step one -- use flux_n, etc, to store the surface values
  // Near left boundary
  {
    (*flux_g[0])[0] = 0.0;
    (*flux_h[0])[0] = 0.0;
    (*flux_n[0])[0] = 0.0;
   
    (*flux_g[0])[1] = sot*((*data_g[nt])[0]+(*data_g[nt])[1]) - oot*((*data_g[nt])[2]-(*data_g[nt])[0]);
    (*flux_h[0])[1] = sot*((*data_h[nt])[0]+(*data_h[nt])[1]) - oot*((*data_h[nt])[2]-(*data_h[nt])[0]);
    (*flux_n[0])[1] = sot*((*data_n[nt])[0]+(*data_n[nt])[1]) - oot*((*data_n[nt])[2]-(*data_n[nt])[0]);
  }
  for ( int i = 2; i < Nx-1; i++ )
  {
    (*flux_g[0])[i] = sot*((*data_g[nt])[i-1]+(*data_g[nt])[i]) - oot*((*data_g[nt])[i-2]+(*data_g[nt])[i+1]);
    (*flux_h[0])[i] = sot*((*data_h[nt])[i-1]+(*data_h[nt])[i]) - oot*((*data_h[nt])[i-2]+(*data_h[nt])[i+1]);
    (*flux_n[0])[i] = sot*((*data_n[nt])[i-1]+(*data_n[nt])[i]) - oot*((*data_n[nt])[i-2]+(*data_n[nt])[i+1]);
  }
  // Near right boundary (biased stencils)
  {
    (*flux_g[0])[Nx-1] = oot*(3.0*(*data_g[nt])[Nx-1]+13.0*(*data_g[nt])[Nx-2]-5.0*(*data_g[nt])[Nx-3]+(*data_g[nt])[Nx-4]);
    (*flux_h[0])[Nx-1] = oot*(3.0*(*data_h[nt])[Nx-1]+13.0*(*data_h[nt])[Nx-2]-5.0*(*data_h[nt])[Nx-3]+(*data_h[nt])[Nx-4]);
    (*flux_n[0])[Nx-1] = oot*(3.0*(*data_n[nt])[Nx-1]+13.0*(*data_n[nt])[Nx-2]-5.0*(*data_n[nt])[Nx-3]+(*data_n[nt])[Nx-4]);

    (*flux_g[0])[Nx] = oot*(25.0*(*data_g[nt])[Nx-1]-23.0*(*data_g[nt])[Nx-2]+13.0*(*data_g[nt])[Nx-3]-3.0*(*data_g[nt])[Nx-4]);
    (*flux_h[0])[Nx] = oot*(25.0*(*data_h[nt])[Nx-1]-23.0*(*data_h[nt])[Nx-2]+13.0*(*data_h[nt])[Nx-3]-3.0*(*data_h[nt])[Nx-4]);
    (*flux_n[0])[Nx] = oot*(25.0*(*data_n[nt])[Nx-1]-23.0*(*data_n[nt])[Nx-2]+13.0*(*data_n[nt])[Nx-3]-3.0*(*data_n[nt])[Nx-4]);
  }

  // Step two -- compute the left limiters in lim_?[0] and the right ones
  //             in lim_?[1]
  // Left boundary
  {
    (*lim_g[0])[0] = 0.0;
    (*lim_h[0])[0] = 0.0;
    (*lim_n[0])[0] = 0.0;
    (*lim_s[0])[0] = 0.0;

    (*lim_g[1])[0] = LimFun_PPM_BC( (*data_g[nt])[0], (*flux_g[0])[1] );
    (*lim_h[1])[0] = LimFun_PPM_BC( (*data_h[nt])[0], (*flux_h[0])[1] );
    (*lim_n[1])[0] = LimFun_PPM_BC( (*data_n[nt])[0], (*flux_n[0])[1] );

    (*lim_s[1])[0] = min( min((*lim_g[1])[0],(*lim_h[1])[0]), (*lim_n[1])[0] );

    if ( tdSetup->slim_tot ) {
      DType loc_lim_s = LimFun_PPM_BC( (*data_g[nt])[0] + (*data_h[nt])[0] + (*data_n[nt])[0],
          (*flux_g[0])[1] + (*flux_h[0])[1] + (*flux_n[0])[1] );
      (*lim_s[1])[0] = min( (*lim_s[1])[0], loc_lim_s );
    }
  }
  for ( int i = 1; i < Nx; i++ ) {
    DType diff_g_l = (*flux_g[0])[i]   - (*data_g[nt])[i];
    DType diff_g_r = (*flux_g[0])[i+1] - (*data_g[nt])[i];
    DType diff_h_l = (*flux_h[0])[i]   - (*data_h[nt])[i];
    DType diff_h_r = (*flux_h[0])[i+1] - (*data_h[nt])[i];
    DType diff_n_l = (*flux_n[0])[i]   - (*data_n[nt])[i];
    DType diff_n_r = (*flux_n[0])[i+1] - (*data_n[nt])[i];
    DType diff_s_l, diff_s_r;

    DType alpha1_g, alpha1_h, alpha1_n, alpha1_s;
    DType alpha2_g, alpha2_h, alpha2_n, alpha2_s;

    bool zero_lim_g = false;
    bool zero_lim_h = false;
    bool zero_lim_n = false;
    bool no_lim_g = fabs(diff_g_l)+fabs(diff_g_r)<2.0*TOL_ZERO;
    bool no_lim_h = fabs(diff_h_l)+fabs(diff_h_r)<2.0*TOL_ZERO;
    bool no_lim_n = fabs(diff_n_l)+fabs(diff_n_r)<2.0*TOL_ZERO;

    if (!no_lim_g) {
      if ( fabs(diff_g_l)<=TOL_ZERO || fabs(diff_g_r)<=TOL_ZERO || diff_g_l*diff_g_r > TOL_ZERO ) {
        zero_lim_g = true;
        (*lim_g[0])[i] = 0.0; (*lim_g[1])[i] = 0.0;
      }
      else {
        alpha1_g = 2.0*fabs(diff_g_l)/fabs(diff_g_r);
        alpha2_g = 0.25 * alpha1_g;
        LimFun_PPM_INT(alpha1_g, alpha2_g, (*lim_g[0])[i], (*lim_g[1])[i]);
      }
    }
    else {
      (*lim_g[0])[i] = 1.0; (*lim_g[1])[i] = 1.0;
      alpha1_g = DTYPE_MAX; alpha2_g = 0.0;
    }

    if (!no_lim_h) {
      if ( fabs(diff_h_l)<=TOL_ZERO || fabs(diff_h_r)<=TOL_ZERO || diff_h_l*diff_h_r > TOL_ZERO ) {
        zero_lim_h = true;
        (*lim_h[0])[i] = 0.0; (*lim_h[1])[i] = 0.0;
      }
      else {
        alpha1_h = 2.0*fabs(diff_h_l)/fabs(diff_h_r);
        alpha2_h = 0.25 * alpha1_h;
        LimFun_PPM_INT(alpha1_h, alpha2_h, (*lim_h[0])[i], (*lim_h[1])[i]);
      }
    }
    else {
      (*lim_h[0])[i] = 1.0; (*lim_h[1])[i] = 1.0;
      alpha1_h = DTYPE_MAX; alpha2_h = 0.0;
    }

    if (!no_lim_n) {
      if ( fabs(diff_n_l)<=TOL_ZERO || fabs(diff_n_r)<=TOL_ZERO || diff_n_l*diff_n_r > TOL_ZERO ) {
        zero_lim_n = true;
        (*lim_n[0])[i] = 0.0; (*lim_n[1])[i] = 0.0;
      }
      else {
        alpha1_n = 2.0*fabs(diff_n_l)/fabs(diff_n_r);
        alpha2_n = 0.25 * alpha1_n;
        LimFun_PPM_INT(alpha1_n, alpha2_n, (*lim_n[0])[i], (*lim_n[1])[i]);
      }
    }
    else {
      (*lim_n[0])[i] = 1.0; (*lim_n[1])[i] = 1.0;
      alpha1_n = DTYPE_MAX; alpha2_n = 0.0;
    }

    bool zero_lim_s = false;
    bool no_lim_s = false; 

    if ( tdSetup->slim_tot ) {
      diff_s_l = diff_g_l + diff_h_l + diff_n_l;
      diff_s_r = diff_g_r + diff_h_r + diff_n_r;
      no_lim_s = fabs(diff_g_l)+fabs(diff_g_r)<2.0*TOL_ZERO;
      if (!no_lim_s) {
        if ( fabs(diff_s_l)<=TOL_ZERO || fabs(diff_s_r)<=TOL_ZERO || diff_s_l*diff_s_r > TOL_ZERO ) {
          zero_lim_s = true;
          (*lim_s[0])[i] = 0.0; (*lim_s[1])[i] = 0.0;
        }
        else {
          alpha1_s = 2.0*fabs(diff_s_l)/fabs(diff_s_r);
          alpha2_s = 0.25 * alpha1_s;
          LimFun_PPM_INT(alpha1_s, alpha2_s, (*lim_s[0])[i], (*lim_s[1])[i]);
        }
      }
      else {
        (*lim_s[0])[i] = 1.0; (*lim_s[1])[i] = 1.0;
        alpha1_s = DTYPE_MAX; alpha2_s = 0.0;
      }
    }

    if ( tdSetup->slim_tot ) {
      if ( zero_lim_g || zero_lim_h || zero_lim_n || zero_lim_s ) {
        (*lim_s[0])[i] = 0.0; (*lim_s[1])[i] = 0.0;
      }
      else {
        DType alpha1 = DTYPE_MAX;
        DType alpha2 = 0.0;
        if ( !no_lim_g ) { alpha1 = min(alpha1,alpha1_g); alpha2 = max(alpha2,alpha2_g); }
        if ( !no_lim_h ) { alpha1 = min(alpha1,alpha1_h); alpha2 = max(alpha2,alpha2_h); }
        if ( !no_lim_n ) { alpha1 = min(alpha1,alpha1_n); alpha2 = max(alpha2,alpha2_n); }
        if ( !no_lim_s ) { alpha1 = min(alpha1,alpha1_s); alpha2 = max(alpha2,alpha2_s); }
        LimFun_PPM_INT(alpha1, alpha2, (*lim_s[0])[i], (*lim_s[1])[i]);
      }
    }
    else {
      if ( zero_lim_g || zero_lim_h || zero_lim_n ) {
        (*lim_s[0])[i] = 0.0; (*lim_s[1])[i] = 0.0;
      }
      else {
        DType alpha1 = DTYPE_MAX;
        DType alpha2 = 0.0;
        if ( !no_lim_g ) { alpha1 = min(alpha1,alpha1_g); alpha2 = max(alpha2,alpha2_g); }
        if ( !no_lim_h ) { alpha1 = min(alpha1,alpha1_h); alpha2 = max(alpha2,alpha2_h); }
        if ( !no_lim_n ) { alpha1 = min(alpha1,alpha1_n); alpha2 = max(alpha2,alpha2_n); }
        LimFun_PPM_INT(alpha1, alpha2, (*lim_s[0])[i], (*lim_s[1])[i]);
      }
    }
  }

  downscale_var(*data_g[nt],*data_g[nt],x_c,Nx,1.0,1);
  downscale_var(*data_h[nt],*data_h[nt],x_c,Nx,1.0,1);
  downscale_var(*data_n[nt],*data_n[nt],x_c,Nx,1.0,1);

  return;
}

void TDSolveFVMMod::update_mat_upw_v(const DType multip, const int nt, HVEC**& vel, 
    const int v_nt, const bool rbc, const int rp_nt) {

  DType dror = data_rp[rp_nt]/data_r[nt];
  DType r2nt = data_r[nt]*data_r[nt];
  
  Entry elem(1,0.0);

  // Cell face at x_f[1]
  {
    DType xc2_l = x_c[0]*x_c[0];
    DType xc2_r = x_c[1]*x_c[1];
    DType xf2   = x_f[1]*x_f[1];

    DType loc_v = (*vel[v_nt])[1]/data_r[nt];
    DType coef_v_l = 0.0;
    DType coef_v_r = 0.0;

    if ( loc_v >= 0.0 )
      coef_v_l = multip * (dtn/dx) * xf2 * r2nt * loc_v;
    else
      coef_v_r = multip * (dtn/dx) * xf2 * r2nt * loc_v;

    elem = mat->get(0,0); elem += coef_v_l / xc2_l; mat->set(0,0,elem);
    elem = mat->get(0,1); elem += coef_v_r / xc2_l; mat->set(0,1,elem);

    elem = mat->get(1,0); elem -= coef_v_l / xc2_r; mat->set(1,0,elem);
    elem = mat->get(1,1); elem -= coef_v_r / xc2_r; mat->set(1,1,elem);
  }
  // Interior cell faces
  for ( int i = 2; i < Nx; i++ ) {
    DType xc2_l = x_c[i-1]*x_c[i-1];
    DType xc2_r = x_c[i]*x_c[i];
    DType xf2   = x_f[i]*x_f[i];

    DType loc_v = (*vel[v_nt])[i]/data_r[nt];
    DType coef_v_l = 0.0;
    DType coef_v_r = 0.0;

    if ( loc_v >= 0.0 )
      coef_v_l = multip * (dtn/dx) * xf2 * r2nt * loc_v;
    else
      coef_v_r = multip * (dtn/dx) * xf2 * r2nt * loc_v;

    elem = mat->get(i-1,i-1); elem += coef_v_l / xc2_l; mat->set(i-1,i-1,elem);
    elem = mat->get(i-1,i);   elem += coef_v_r / xc2_l; mat->set(i-1,i,elem);

    elem = mat->get(i,i-1);   elem -= coef_v_l / xc2_r; mat->set(i,i-1,elem);
    elem = mat->get(i,i);     elem -= coef_v_r / xc2_r; mat->set(i,i,elem);
  }

  // Cell face x_f[Nx]
  if (rbc) {
    DType xc2 = x_c[Nx-1]*x_c[Nx-1];
    DType xf2 = x_f[Nx]*x_f[Nx];

    DType loc_v = (*vel[v_nt])[Nx]/data_r[nt];
    DType coef_v_l = 0.0;

    if ( loc_v >= 0.0 )
      coef_v_l = multip * (dtn/dx) * xf2 * r2nt * loc_v;

    elem = mat->get(Nx-1,Nx-1); elem += coef_v_l / xc2; mat->set(Nx-1,Nx-1,elem);
  }

  return;
}

void TDSolveFVMMod::update_mat_orig_muscl_v(const int nt, HVEC**& vel, const int v_nt, HVEC**& lim, const bool rbc) {

  Entry elem(1,0.0);

  // Cell face at x_f[1]
  {
    DType loc_v = (*vel[v_nt])[1];
    DType coef_v_l  = 0.0;
    DType coef_v_r  = 0.0;
    DType coef_v_rr = 0.0;
    DType lim_r = (*lim[0])[1];

    if ( loc_v >= 0.0 )
      coef_v_l  = (dtn/dx) * loc_v;
    else {
      coef_v_r  = (dtn/dx) * loc_v * ( 1.0 + 0.5*lim_r );
      coef_v_rr = (dtn/dx) * loc_v * ( - 0.5*lim_r );
    }

    elem = mat->get(0,0); elem += coef_v_l;  mat->set(0,0,elem);
    elem = mat->get(0,1); elem += coef_v_r;  mat->set(0,1,elem);
    elem = mat->get(0,2); elem += coef_v_rr; mat->set(0,2,elem);
    
    elem = mat->get(1,0); elem -= coef_v_l;  mat->set(1,0,elem);
    elem = mat->get(1,1); elem -= coef_v_r;  mat->set(1,1,elem);
    elem = mat->get(1,2); elem -= coef_v_rr; mat->set(1,2,elem);
  }
  // Interior cell faces
  for ( int i = 2; i < Nx-1; i++ ) {
    DType loc_v = (*vel[v_nt])[i];
    DType coef_v_l  = 0.0;
    DType coef_v_r  = 0.0;
    DType coef_v_rr = 0.0;
    DType lim_l = (*lim[0])[i-1];
    DType lim_r = (*lim[0])[i];

    if ( loc_v >= 0.0 ) {
      coef_v_l  = (dtn/dx) * loc_v * ( 1.0 - 0.5*lim_l );
      coef_v_r  = (dtn/dx) * loc_v * ( 0.5*lim_l );
    }
    else {
      coef_v_r  = (dtn/dx) * loc_v * ( 1.0 + 0.5*lim_r );
      coef_v_rr = (dtn/dx) * loc_v * ( - 0.5*lim_r );
    }

    elem = mat->get(i-1,i-1); elem += coef_v_l;  mat->set(i-1,i-1,elem);
    elem = mat->get(i-1,i);   elem += coef_v_r;  mat->set(i-1,i,  elem);
    elem = mat->get(i-1,i+1); elem += coef_v_rr; mat->set(i-1,i+1,elem);

    elem = mat->get(i,i-1); elem -= coef_v_l;  mat->set(i,i-1,elem);
    elem = mat->get(i,i);   elem -= coef_v_r;  mat->set(i,i,  elem);
    elem = mat->get(i,i+1); elem -= coef_v_rr; mat->set(i,i+1,elem);
  }
  // Cell face x_f[Nx-1]
  {
    DType loc_v = (*vel[v_nt])[Nx-1];
    DType coef_v_l = 0.0;
    DType coef_v_r = 0.0;
    DType lim_l = (*lim[0])[Nx-2];

    if ( loc_v >= 0.0 ) {
      coef_v_l = (dtn/dx) * loc_v * ( 1.0 - 0.5*lim_l );
      coef_v_r = (dtn/dx) * loc_v * ( 0.5*lim_l );
    }
    else
      coef_v_r = (dtn/dx) * loc_v;

    elem = mat->get(Nx-2,Nx-2); elem += coef_v_l; mat->set(Nx-2,Nx-2,elem);
    elem = mat->get(Nx-2,Nx-1); elem += coef_v_r; mat->set(Nx-2,Nx-1,elem);

    elem = mat->get(Nx-1,Nx-2); elem -= coef_v_l; mat->set(Nx-1,Nx-2,elem);
    elem = mat->get(Nx-1,Nx-1); elem -= coef_v_r; mat->set(Nx-1,Nx-1,elem);
  }

  // Cell face x_f[Nx]
  if (rbc) {
    DType loc_v = (*vel[v_nt])[Nx];
    DType coef_v_l = 0.0;

    if ( loc_v >= 0.0 )
      coef_v_l = (dtn/dx) * loc_v;

    elem = mat->get(Nx-1,Nx-1); elem += coef_v_l; mat->set(Nx-1,Nx-1,elem);
  }

  return;
}

void TDSolveFVMMod::update_mat_muscl_v(const DType multip, 
    const int nt, HVEC**& vel, const int v_nt, const bool rbc, const int rp_nt) {

  DType dror = data_rp[rp_nt]/data_r[nt];
  DType r2nt = data_r[nt]*data_r[nt];

  Entry elem(1,0.0);

  // Cell face at x_f[1]
  {
    DType xc2_l = x_c[0]*x_c[0];
    DType xc2_r = x_c[1]*x_c[1];
    DType xf2   = x_f[1]*x_f[1];

    DType loc_v = (*vel[v_nt])[1]/data_r[nt];
    DType coef_v_l  = 0.0;
    DType coef_v_r  = 0.0;
    DType coef_v_rr = 0.0;
    DType lim_r = (*lim_muscl[0])[1];

    if ( loc_v >= 0.0 )
      coef_v_l  = (dtn/dx) * xf2 * r2nt * loc_v;
    else {
      coef_v_r  = (dtn/dx) * xf2 * r2nt * loc_v * ( 1.0 + 0.5*lim_r );
      coef_v_rr = (dtn/dx) * xf2 * r2nt * loc_v * ( - 0.5*lim_r );
    }

    elem = mat->get(0,0); elem += coef_v_l  / xc2_l;  mat->set(0,0,elem);
    elem = mat->get(0,1); elem += coef_v_r  / xc2_l;  mat->set(0,1,elem);
    elem = mat->get(0,2); elem += coef_v_rr / xc2_l; mat->set(0,2,elem);
    
    elem = mat->get(1,0); elem -= coef_v_l  / xc2_r;  mat->set(1,0,elem);
    elem = mat->get(1,1); elem -= coef_v_r  / xc2_r;  mat->set(1,1,elem);
    elem = mat->get(1,2); elem -= coef_v_rr / xc2_r; mat->set(1,2,elem);
  }
  // Interior cell faces
  for ( int i = 2; i < Nx-1; i++ ) {
    DType xc2_l = x_c[i-1]*x_c[i-1];
    DType xc2_r = x_c[i]*x_c[i];
    DType xf2   = x_f[i]*x_f[i];

    DType loc_v = (*vel[v_nt])[i]/data_r[nt];
    DType coef_v_l  = 0.0;
    DType coef_v_r  = 0.0;
    DType coef_v_rr = 0.0;
    DType lim_l = (*lim_muscl[0])[i-1];
    DType lim_r = (*lim_muscl[0])[i];

    if ( loc_v >= 0.0 ) {
      coef_v_l  = (dtn/dx) * xf2 * r2nt * loc_v * ( 1.0 - 0.5*lim_l );
      coef_v_r  = (dtn/dx) * xf2 * r2nt * loc_v * ( 0.5*lim_l );
    }
    else {
      coef_v_r  = (dtn/dx) * xf2 * r2nt * loc_v * ( 1.0 + 0.5*lim_r );
      coef_v_rr = (dtn/dx) * xf2 * r2nt * loc_v * ( - 0.5*lim_r );
    }

    elem = mat->get(i-1,i-1); elem += coef_v_l  / xc2_l;  mat->set(i-1,i-1,elem);
    elem = mat->get(i-1,i);   elem += coef_v_r  / xc2_l;  mat->set(i-1,i,  elem);
    elem = mat->get(i-1,i+1); elem += coef_v_rr / xc2_l; mat->set(i-1,i+1,elem);

    elem = mat->get(i,i-1); elem -= coef_v_l  / xc2_r;  mat->set(i,i-1,elem);
    elem = mat->get(i,i);   elem -= coef_v_r  / xc2_r;  mat->set(i,i,  elem);
    elem = mat->get(i,i+1); elem -= coef_v_rr / xc2_r; mat->set(i,i+1,elem);
  }
  // Cell face x_f[Nx-1]
  {
    DType xc2_l = x_c[Nx-2]*x_c[Nx-2];
    DType xc2_r = x_c[Nx-1]*x_c[Nx-1];
    DType xf2   = x_f[Nx-1]*x_f[Nx-1];

    DType loc_v = (*vel[v_nt])[Nx-1]/data_r[nt];
    DType coef_v_l = 0.0;
    DType coef_v_r = 0.0;
    DType lim_l = (*lim_muscl[0])[Nx-2];

    if ( loc_v >= 0.0 ) {
      coef_v_l = (dtn/dx) * xf2 * r2nt * loc_v * ( 1.0 - 0.5*lim_l );
      coef_v_r = (dtn/dx) * xf2 * r2nt * loc_v * ( 0.5*lim_l );
    }
    else
      coef_v_r = (dtn/dx) * xf2 * r2nt * loc_v;

    elem = mat->get(Nx-2,Nx-2); elem += coef_v_l / xc2_l; mat->set(Nx-2,Nx-2,elem);
    elem = mat->get(Nx-2,Nx-1); elem += coef_v_r / xc2_l; mat->set(Nx-2,Nx-1,elem);

    elem = mat->get(Nx-1,Nx-2); elem -= coef_v_l / xc2_r; mat->set(Nx-1,Nx-2,elem);
    elem = mat->get(Nx-1,Nx-1); elem -= coef_v_r / xc2_r; mat->set(Nx-1,Nx-1,elem);
  }

  // Cell face x_f[Nx]
  if (rbc) {
    DType xc2 = x_c[Nx-1]*x_c[Nx-1];
    DType xf2 = x_f[Nx]*x_f[Nx];

    DType loc_v = (*vel[v_nt])[Nx]/data_r[nt];
    DType coef_v_l = 0.0;

    if ( loc_v >= 0.0 )
      coef_v_l = (dtn/dx) * xf2 * r2nt * loc_v;

    elem = mat->get(Nx-1,Nx-1); elem += coef_v_l / xc2; mat->set(Nx-1,Nx-1,elem);
  }

  return;
}

void TDSolveFVMMod::update_mat_ppm_r(const DType multip, const int nt, const int rp_nt) {

  DType dror = data_rp[rp_nt]/data_r[nt];
  DType r2nt = data_r[nt]*data_r[nt];

  Entry elem(1,0.0);

  DType sot = 7.0/12.0;
  DType oot = 1.0/12.0;

  // Cell face at x_f[1]
  {
    DType xc2_l = x_c[0]*x_c[0];
    DType xc2_r = x_c[1]*x_c[1];
    DType xf2   = x_f[1]*x_f[1];

    DType xc3_l  = x_c[0]*x_c[0]*x_c[0];
    DType xc3_r  = x_c[1]*x_c[1]*x_c[1];
    DType xc3_rr = x_c[2]*x_c[2]*x_c[2];

    DType loc_r = - dror;
    DType lim_r_l = (*lim_s[1])[0];
    DType lim_r_r = (*lim_s[0])[1];
    DType coef_r_l  = 0.0;
    DType coef_r_r  = 0.0;
    DType coef_r_rr = 0.0;

    if ( loc_r >= 0.0 ) {
      coef_r_l  = multip * (dtn/dx) * r2nt * loc_r * xc3_l  * ( (1.0-lim_r_l) + lim_r_l*(sot+oot) );
      coef_r_r  = multip * (dtn/dx) * r2nt * loc_r * xc3_r  * ( lim_r_l*sot );
      coef_r_rr = multip * (dtn/dx) * r2nt * loc_r * xc3_rr * ( -lim_r_l*oot );
    }
    else {
      coef_r_l  = multip * (dtn/dx) * r2nt * loc_r * xc3_l  * ( lim_r_r*(sot+oot) );
      coef_r_r  = multip * (dtn/dx) * r2nt * loc_r * xc3_r  * ( (1.0-lim_r_r) + lim_r_r*sot );
      coef_r_rr = multip * (dtn/dx) * r2nt * loc_r * xc3_rr * ( -lim_r_r*oot );
    }

    elem = mat->get(0,0); elem += coef_r_l  / xc2_l; mat->set(0,0,elem);
    elem = mat->get(0,1); elem += coef_r_r  / xc2_l; mat->set(0,1,elem);
    elem = mat->get(0,2); elem += coef_r_rr / xc2_l; mat->set(0,2,elem);

    elem = mat->get(1,0); elem -= coef_r_l  / xc2_r; mat->set(1,0,elem);
    elem = mat->get(1,1); elem -= coef_r_r  / xc2_r; mat->set(1,1,elem);
    elem = mat->get(1,2); elem -= coef_r_rr / xc2_r; mat->set(1,2,elem);
  }
  // Interior cell faces
  for ( int i = 2; i < Nx-1; i++ ) {
    DType xc2_l = x_c[i-1]*x_c[i-1];
    DType xc2_r = x_c[i]*x_c[i];
    DType xf2   = x_f[i]*x_f[i];

    DType xc3_ll = x_c[i-2]*x_c[i-2]*x_c[i-2];
    DType xc3_l  = x_c[i-1]*x_c[i-1]*x_c[i-1];
    DType xc3_r  = x_c[i]*x_c[i]*x_c[i];
    DType xc3_rr = x_c[i+1]*x_c[i+1]*x_c[i+1];

    DType loc_r = -dror;
    DType lim_r_l = (*lim_s[1])[i-1];
    DType lim_r_r = (*lim_s[0])[i];
    DType coef_r_ll = 0.0;
    DType coef_r_l  = 0.0;
    DType coef_r_r  = 0.0;
    DType coef_r_rr = 0.0;

    if ( loc_r >= 0.0 ) {
      coef_r_ll = multip * (dtn/dx) * r2nt * loc_r * xc3_ll * ( -lim_r_l*oot );
      coef_r_l  = multip * (dtn/dx) * r2nt * loc_r * xc3_l  * ( (1.0-lim_r_l) + lim_r_l*sot );
      coef_r_r  = multip * (dtn/dx) * r2nt * loc_r * xc3_r  * ( lim_r_l*sot );
      coef_r_rr = multip * (dtn/dx) * r2nt * loc_r * xc3_rr * ( -lim_r_l*oot );
    }
    else {
      coef_r_ll = multip * (dtn/dx) * r2nt * loc_r * xc3_ll * ( -lim_r_r*oot );
      coef_r_l  = multip * (dtn/dx) * r2nt * loc_r * xc3_l  * ( lim_r_r*sot );
      coef_r_r  = multip * (dtn/dx) * r2nt * loc_r * xc3_r  * ( (1.0-lim_r_r) + lim_r_r*sot );
      coef_r_rr = multip * (dtn/dx) * r2nt * loc_r * xc3_rr * ( -lim_r_r*oot );
    }

    elem = mat->get(i-1,i-2); elem += coef_r_ll / xc2_l; mat->set(i-1,i-2,elem);
    elem = mat->get(i-1,i-1); elem += coef_r_l  / xc2_l; mat->set(i-1,i-1,elem);
    elem = mat->get(i-1,i);   elem += coef_r_r  / xc2_l; mat->set(i-1,i,elem);
    elem = mat->get(i-1,i+1); elem += coef_r_rr / xc2_l; mat->set(i-1,i+1,elem);

    elem = mat->get(i,i-2);   elem -= coef_r_ll / xc2_r; mat->set(i,i-2,elem);
    elem = mat->get(i,i-1);   elem -= coef_r_l  / xc2_r; mat->set(i,i-1,elem);
    elem = mat->get(i,i);     elem -= coef_r_r  / xc2_r; mat->set(i,i,elem);
    elem = mat->get(i,i+1);   elem -= coef_r_rr / xc2_r; mat->set(i,i+1,elem);
  }
  // Cell face x_f[Nx-1]
  {
    DType xc2_l = x_c[Nx-2]*x_c[Nx-2];
    DType xc2_r = x_c[Nx-1]*x_c[Nx-1];
    DType xf2   = x_f[Nx-1]*x_f[Nx-1];

    DType xc3_l3 = x_c[Nx-4]*x_c[Nx-4]*x_c[Nx-4];
    DType xc3_ll = x_c[Nx-3]*x_c[Nx-3]*x_c[Nx-3];
    DType xc3_l  = x_c[Nx-2]*x_c[Nx-2]*x_c[Nx-2];
    DType xc3_r  = x_c[Nx-1]*x_c[Nx-1]*x_c[Nx-1];

    DType loc_r = -dror;
    DType lim_r_l = (*lim_s[1])[Nx-2];
    DType lim_r_r = (*lim_s[0])[Nx-1];
    DType coef_r_l3 = 0.0;
    DType coef_r_ll = 0.0;
    DType coef_r_l  = 0.0;
    DType coef_r_r  = 0.0;

    if ( loc_r >= 0.0 ) {
      coef_r_l3 = multip * (dtn/dx) * r2nt * loc_r * xc3_l3 * ( lim_r_l*oot );
      coef_r_ll = multip * (dtn/dx) * r2nt * loc_r * xc3_ll * ( -lim_r_l*5.0*oot );
      coef_r_l  = multip * (dtn/dx) * r2nt * loc_r * xc3_l  * ( (1.0-lim_r_l) + lim_r_l*13.0*oot );
      coef_r_r  = multip * (dtn/dx) * r2nt * loc_r * xc3_r  * ( lim_r_l*3.0*oot );
    }
    else {
      coef_r_l3 = multip * (dtn/dx) * r2nt * loc_r * xc3_l3 * ( lim_r_r*oot );
      coef_r_ll = multip * (dtn/dx) * r2nt * loc_r * xc3_ll * ( -lim_r_r*5.0*oot );
      coef_r_l  = multip * (dtn/dx) * r2nt * loc_r * xc3_l  * ( lim_r_r*13.0*oot );
      coef_r_r  = multip * (dtn/dx) * r2nt * loc_r * xc3_r  * ( (1.0-lim_r_r) + lim_r_r*3.0*oot );
    }

    elem = mat->get(Nx-2,Nx-4); elem += coef_r_l3 / xc2_l; mat->set(Nx-2,Nx-4,elem);
    elem = mat->get(Nx-2,Nx-3); elem += coef_r_ll / xc2_l; mat->set(Nx-2,Nx-3,elem);
    elem = mat->get(Nx-2,Nx-2); elem += coef_r_l  / xc2_l; mat->set(Nx-2,Nx-2,elem);
    elem = mat->get(Nx-2,Nx-1); elem += coef_r_r  / xc2_l; mat->set(Nx-2,Nx-1,elem);

    elem = mat->get(Nx-1,Nx-4); elem -= coef_r_l3 / xc2_r; mat->set(Nx-1,Nx-4,elem);
    elem = mat->get(Nx-1,Nx-3); elem -= coef_r_ll / xc2_r; mat->set(Nx-1,Nx-3,elem);
    elem = mat->get(Nx-1,Nx-2); elem -= coef_r_l  / xc2_r; mat->set(Nx-1,Nx-2,elem);
    elem = mat->get(Nx-1,Nx-1); elem -= coef_r_r  / xc2_r; mat->set(Nx-1,Nx-1,elem);
  }

  return;
}

void TDSolveFVMMod::update_cells_imp_upw(const int fr, const int to, const DType tsol_loc) {
  DType lambda = tdCoef->lambda;
  DType q      = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

  DType dror = data_rp[0]/data_r[to];
  DType r2fr = data_r[fr]*data_r[fr];
  DType r2to = data_r[to]*data_r[to];

  DType n_bc = calc_n_bc( tsol_loc );

  // Compute the limiters and store them in lim_?, as well as 
  // unlimited interpolated values, and store them in flux_?
  // We won't use flux_? here, though
  compute_flux_r_ppm_lims( to );

  downscale_var(*data_g[fr],*data_g[fr],x_c,Nx,data_r[fr]);
  downscale_var(*data_h[fr],*data_h[fr],x_c,Nx,data_r[fr]);
  downscale_var(*data_n[fr],*data_n[fr],x_c,Nx,data_r[fr]);

  VComp comp(1,0.0);
  Entry elem(1,0.0);

  // First data_g
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  // RHS
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_g[fr])[i] * r2fr;
    rhs->set(i,comp);
  }

  // Matrix
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[i],p) ); 
    elem = r2to - dtn * r2to * ( loc_lambda - loc_mu - dror );
    mat->set(i,i,elem);
  }

  update_mat_upw_v( 1.0, to, data_v, to );
  update_mat_ppm_r( 1.0, to );

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
    comp = (*data_h[fr])[i] * r2fr + dtn * loc_mu * (*data_g[to])[i] * r2to;
    rhs->set(i,comp);
  }

  // Matrix
  for ( int i = 0; i < Nx; i++ ) {
    elem = r2to - dtn * r2to * ( - delta_h - dror );
    mat->set(i,i,elem);
  }

  update_mat_upw_v( 1.0, to, data_v, to );
  update_mat_ppm_r( 1.0, to );

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

  // RHS
  for ( int i = 0; i < Nx-1; i++ ) {
    comp = (*data_n[fr])[i] * r2fr;
    rhs->set(i,comp);
  }
  // Right boundary
  {
    DType xc2 = x_c[Nx-1]*x_c[Nx-1];
    DType xf2 = x_f[Nx]*x_f[Nx];

    DType loc_u = (*data_u[0])[Nx]/data_r[to];
    DType coef_u_r = 0.0;
    if ( loc_u < 0.0 )
      coef_u_r += (dtn/dx) * xf2 * r2to * loc_u;
    comp = (*data_n[fr])[Nx-1] * r2fr - coef_u_r * n_bc / xc2;
    rhs->set(Nx-1,comp);
  }

  // Matrix
  for ( int i = 0; i < Nx; i++ ) {
    elem = r2to - dtn * r2to * ( - delta_n - dror );
    mat->set(i,i,elem);
  }

  update_mat_upw_v( 1.0, to, data_v, to );
  update_mat_ppm_r( 1.0, to );
  update_mat_upw_v( 1.0, to, data_u, 0, true );

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

  upscale_var(*data_g[fr],*data_g[fr],x_c,Nx,data_r[fr]);
  upscale_var(*data_h[fr],*data_h[fr],x_c,Nx,data_r[fr]);
  upscale_var(*data_n[fr],*data_n[fr],x_c,Nx,data_r[fr]);

  upscale_var(*data_g[to],*data_g[to],x_c,Nx,data_r[to]);
  upscale_var(*data_h[to],*data_h[to],x_c,Nx,data_r[to]);
  upscale_var(*data_n[to],*data_n[to],x_c,Nx,data_r[to]);
  
  return;
}

void TDSolveFVMMod::update_cells_imp_muscl(const int fr, const int to, const DType tsol_loc) {
  DType lambda = tdCoef->lambda;
  DType q      = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

  DType dror = data_rp[0]/data_r[to];
  DType r2fr = data_r[fr]*data_r[fr];
  DType r2to = data_r[to]*data_r[to];

  DType n_bc = calc_n_bc( tsol_loc );

  compute_flux_r_ppm_lims( to );

  downscale_var(*data_g[fr],*data_g[fr],x_c,Nx,data_r[fr]);
  downscale_var(*data_h[fr],*data_h[fr],x_c,Nx,data_r[fr]);
  downscale_var(*data_n[fr],*data_n[fr],x_c,Nx,data_r[fr]);

  // XYZ -- this seems to be a bug in the previous implementation
  downscale_var(*data_g[to],*data_g[to],x_c,Nx,data_r[to]);
  downscale_var(*data_h[to],*data_h[to],x_c,Nx,data_r[to]);
  downscale_var(*data_n[to],*data_n[to],x_c,Nx,data_r[to]);

  // For the modified FVM, use primitive variables to compute muscl limiters
  compute_flux_muscl_lims( to );
  
  VComp comp(1,0.0);
  Entry elem(1,0.0);

  // First data_g
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  // RHS
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_g[fr])[i] * r2fr;
    rhs->set(i,comp);
  }

  // Matrix
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[i],p) ); 
    elem = r2to - dtn * r2to * ( loc_lambda - loc_mu - dror );
    mat->set(i,i,elem);
  }

  update_mat_muscl_v( 1.0, to, data_v, to );
  update_mat_ppm_r( 1.0, to );

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
    comp = (*data_h[fr])[i] * r2fr + dtn * loc_mu * (*data_g[to])[i] * r2to;
    rhs->set(i,comp);
  }

  // Matrix
  for ( int i = 0; i < Nx; i++ ) {
    elem = r2to - dtn * r2to * ( - delta_h - dror );
    mat->set(i,i,elem);
  }

  update_mat_muscl_v( 1.0, to, data_v, to );
  update_mat_ppm_r( 1.0, to );

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

  // RHS
  for ( int i = 0; i < Nx-1; i++ ) {
    comp = (*data_n[fr])[i] * r2fr;
    rhs->set(i,comp);
  }
  // Right boundary
  {
    DType xc2 = x_c[Nx-1]*x_c[Nx-1];
    DType xf2 = x_f[Nx]*x_f[Nx];

    DType loc_u = (*data_u[0])[Nx]/data_r[to];
    DType coef_u_r = 0.0;
    if ( loc_u < 0.0 )
      coef_u_r += (dtn/dx) * xf2 * r2to * loc_u;
    comp = (*data_n[fr])[Nx-1] * r2fr - coef_u_r * n_bc / xc2;
    rhs->set(Nx-1,comp);
  }

  // Matrix
  for ( int i = 0; i < Nx; i++ ) {
    elem = r2to - dtn * r2to * ( - delta_n - dror );
    mat->set(i,i,elem);
  }

  update_mat_muscl_v( 1.0, to, data_v, to );
  update_mat_ppm_r( 1.0, to );
  update_mat_muscl_v( 1.0, to, data_u, 0, true );

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

  upscale_var(*data_g[fr],*data_g[fr],x_c,Nx,data_r[fr]);
  upscale_var(*data_h[fr],*data_h[fr],x_c,Nx,data_r[fr]);
  upscale_var(*data_n[fr],*data_n[fr],x_c,Nx,data_r[fr]);

  upscale_var(*data_g[to],*data_g[to],x_c,Nx,data_r[to]);
  upscale_var(*data_h[to],*data_h[to],x_c,Nx,data_r[to]);
  upscale_var(*data_n[to],*data_n[to],x_c,Nx,data_r[to]);

  return;
}

void TDSolveFVMMod::update_vel_int(const int nt, const int u_nt) {

  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;

  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0
  
  // First do the integration of x^2R(lambda G - delta_h H)
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    (*data_v[nt])[i+1] = (*data_v[nt])[i] + dx * ( loc_lambda * (*data_g[nt])[i] - delta_h * (*data_h[nt])[i] - delta_n * (*data_n[nt])[i] ) / data_r[nt];
  }

  // Remove the advection part due to chemical gradient
  for ( int i = 1; i < Nx+1; i++ )
    (*data_v[nt])[i] = ((*data_v[nt])[i] - (*flux_n_u[u_nt])[i]/data_r[nt])/x_f[i]/x_f[i];

  return;
}

void TDSolveFVMMod::update_chem_imp( const int fr, const int to, 
    const DType shift_dt_fr, const DType shift_dt_to ) {
  // Chemicals -- prob_id should be non-negative
  switch ( tdSetup->prob_id ) {
    case 0: // PDGF
      switch ( tdSetup->schem ) {
        case 0: // Implicit first-order upwind for advection and cd for diffusion
          update_chem_upw_be( fr, to );
          break;
        case 1: // Explicit second-order muscl for advection and implicit cd for diffusion
          update_chem_muscl_imex( fr, to );
          break;
        case 2: // Implicit second-order muscl for advection and implicit cd for diffusion
          update_chem_muscl_be( fr, to );
          break;
        default:
          std::cout << "Warning: The spatial scheme for chem not defined, "
            << "using upw/backward on the primitive equation instead\n" << std::endl;
          update_chem_upw_be( fr, to );
          break;
      }
      break;
    case 1: // OPN
      switch ( tdSetup->schem ) {
        case 0: // Implicit first-order upwind for advection and cd for diffusion
          // Update A
          update_chem_upw_be( fr, to );
          // Update P1 and P2
          update_chem_opn_upw_be( tsol+shift_dt_fr*dtn, fr, to, tdCoef->D_1, tdCoef->alpha_1, 
              tdCoef->beta_1, tdCoef->gamma_1, data_n, data_p1, tdCoef->treatP1, tdCoef->treatP1BC );
          update_chem_opn_upw_be( tsol+shift_dt_fr*dtn, fr, to, tdCoef->D_2, 
              tdCoef->alpha_2, tdCoef->beta_2, tdCoef->gamma_2, data_g, data_p2 );
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
        case 2: // Implicit second-order muscl for advection and implicit cd for diffusion
          // Update A
          update_chem_muscl_be( fr, to ); 
          // Update P1 and P2
          update_chem_opn_muscl_be( tsol+shift_dt_fr*dtn, fr, to, tdCoef->D_1, tdCoef->alpha_1, 
              tdCoef->beta_1, tdCoef->gamma_1, data_n, data_p1, flux_p1, slope_p1, 
              tdCoef->treatP1, tdCoef->treatP1BC );
          update_chem_opn_muscl_be( tsol+shift_dt_fr*dtn, fr, to, tdCoef->D_2, tdCoef->alpha_2, 
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
  return;
}

void TDSolveFVMMod::compute_flux_chem_muscl_lims(const int nt, HVEC**& data, HVEC**& lim, const int N_loc, const DType* xc_loc, const DType* xf_loc, const int rbc, const DType bc_val) {
  DType r2nt = data_r[nt] * data_r[nt];
  // Left boundary
  (*lim[0])[0] = 0.0;
  // Interior cells
  for ( int i = 1; i < N_loc-1; i++ ) {
    DType xc2_l = xc_loc[i-1]*xc_loc[i-1];
    DType xc2_c = xc_loc[i]*xc_loc[i];
    DType xc2_r = xc_loc[i+1]*xc_loc[i+1];

    DType data_loc_l = (*data[nt])[i-1]*xc2_l*r2nt;
    DType data_loc_c = (*data[nt])[i]  *xc2_c*r2nt;
    DType data_loc_r = (*data[nt])[i+1]*xc2_r*r2nt;

    (*lim[0])[i] = LimFun_Minmod( data_loc_c-data_loc_l, data_loc_r-data_loc_c );
  }
  // Right boundary
  if ( rbc == 1 ) { // Flux boundary
    DType xc2_l = xc_loc[N_loc-2]*xc_loc[N_loc-2];
    DType xc2_c = xc_loc[N_loc-1]*xc_loc[N_loc-1];
    DType xc2_r = (2.0*xc_loc[N_loc-1]-xc_loc[N_loc-2])*(2.0*xc_loc[N_loc-1]-xc_loc[N_loc-2]);

    DType data_loc_l = (*data[nt])[N_loc-2]*xc2_l*r2nt;
    DType data_loc_c = (*data[nt])[N_loc-1]*xc2_c*r2nt;
    DType data_loc_r = ((*data[nt])[N_loc-1]-dx*data_r[nt]*bc_val)*xc2_r*r2nt;
    // XYZ -- not finished?
  }
  else { // Dirichlet b.c., simply set limiter to zero
    (*lim[0])[N_loc-1] = 0.0;
  }

  return;
}

void TDSolveFVMMod::update_mat_chem_muscl_r(const int nt, LinMat*& mat_loc, HVEC**& lim, const int N_loc, const DType* xc_loc, const DType* xf_loc, const int rbc, const int rp_nt) {

  DType dror = data_rp[rp_nt]/data_r[nt];
  DType r2nt = data_r[nt]*data_r[nt];

  Entry elem(1,0.0);

  // Cell face at xf_loc[1]
  {
    DType xc2_l  = xc_loc[0]*xc_loc[0];
    DType xc2_r  = xc_loc[1]*xc_loc[1];
    DType xc2_rr = xc_loc[2]*xc_loc[2];

    DType loc_v = - xf_loc[1] * dror;
    DType coef_v_l  = 0.0;
    DType coef_v_r  = 0.0;
    DType coef_v_rr = 0.0;
    DType lim_r = (*lim[0])[1];

    if ( loc_v >= 0.0 )
      coef_v_l  = (dtn/dx) * xc2_l  * r2nt * loc_v;
    else {
      coef_v_r  = (dtn/dx) * xc2_r  * r2nt * loc_v * ( 1.0 + 0.5*lim_r );
      coef_v_rr = (dtn/dx) * xc2_rr * r2nt * loc_v * ( - 0.5*lim_r );
    }

    elem = mat_loc->get(0,0); elem += coef_v_l  / xc2_l / r2nt; mat_loc->set(0,0,elem);
    elem = mat_loc->get(0,1); elem += coef_v_r  / xc2_l / r2nt; mat_loc->set(0,1,elem);
    elem = mat_loc->get(0,2); elem += coef_v_rr / xc2_l / r2nt; mat_loc->set(0,2,elem);

    elem = mat_loc->get(1,0); elem -= coef_v_l  / xc2_r / r2nt; mat_loc->set(1,0,elem);
    elem = mat_loc->get(1,1); elem -= coef_v_r  / xc2_r / r2nt; mat_loc->set(1,1,elem);
    elem = mat_loc->get(1,2); elem -= coef_v_rr / xc2_r / r2nt; mat_loc->set(1,2,elem);
  }
  // Interior cell faces
  for ( int i = 2; i < N_loc-1; i++ ) {
    DType xc2_l  = xc_loc[i-1]*xc_loc[i-1];
    DType xc2_r  = xc_loc[i]*xc_loc[i];
    DType xc2_rr = xc_loc[i+1]*xc_loc[i+1];

    DType loc_v = - xf_loc[i] * dror;
    DType coef_v_l  = 0.0;
    DType coef_v_r  = 0.0;
    DType coef_v_rr = 0.0;
    DType lim_l = (*lim[0])[i-1];
    DType lim_r = (*lim[0])[i];

    if ( loc_v >= 0.0 ) {
      coef_v_l  = (dtn/dx) * xc2_l  * r2nt * loc_v * ( 1.0 - 0.5*lim_l );
      coef_v_r  = (dtn/dx) * xc2_r  * r2nt * loc_v * ( 0.5*lim_l );
    }
    else {
      coef_v_r  = (dtn/dx) * xc2_r  * r2nt * loc_v * ( 1.0 + 0.5*lim_r );
      coef_v_rr = (dtn/dx) * xc2_rr * r2nt * loc_v * ( - 0.5*lim_r );
    }

    elem = mat_loc->get(i-1,i-1); elem += coef_v_l  / xc2_l / r2nt; mat_loc->set(i-1,i-1,elem);
    elem = mat_loc->get(i-1,i);   elem += coef_v_r  / xc2_l / r2nt; mat_loc->set(i-1,i,  elem);
    elem = mat_loc->get(i-1,i+1); elem += coef_v_rr / xc2_l / r2nt; mat_loc->set(i-1,i+1,elem);

    elem = mat_loc->get(i,i-1); elem -= coef_v_l  / xc2_r / r2nt; mat_loc->set(i,i-1,elem);
    elem = mat_loc->get(i,i);   elem -= coef_v_r  / xc2_r / r2nt; mat_loc->set(i,i,  elem);
    elem = mat_loc->get(i,i+1); elem -= coef_v_rr / xc2_r / r2nt; mat_loc->set(i,i+1,elem);
  }
  // Right boundary
  {
    DType xc2_l = xc_loc[N_loc-2]*xc_loc[N_loc-2];
    DType xc2_r = xc_loc[N_loc-1]*xc_loc[N_loc-1];

    DType loc_v = - xc_loc[N_loc-1] * dror;
    DType coef_v_l = 0.0;
    DType coef_v_r = 0.0;
    DType lim_l = (*lim[0])[N_loc-2];

    if ( loc_v >= 0.0 ) {
      coef_v_l = (dtn/dx) * xc2_l * r2nt * loc_v * ( 1.0 - 0.5*lim_l );
      coef_v_r = (dtn/dx) * xc2_r * r2nt * loc_v * ( 0.5*lim_l );
    }
    else
      coef_v_r = (dtn/dx) * xc2_r * r2nt * loc_v;

    elem = mat_loc->get(N_loc-2,N_loc-2); elem += coef_v_l / xc2_l / r2nt; mat_loc->set(N_loc-2,N_loc-2,elem);
    elem = mat_loc->get(N_loc-2,N_loc-1); elem += coef_v_r / xc2_l / r2nt; mat_loc->set(N_loc-2,N_loc-1,elem);

    elem = mat_loc->get(N_loc-1,N_loc-2); elem -= coef_v_l / xc2_r / r2nt; mat_loc->set(N_loc-1,N_loc-2,elem);
    elem = mat_loc->get(N_loc-1,N_loc-1); elem -= coef_v_r / xc2_r / r2nt; mat_loc->set(N_loc-1,N_loc-1,elem);
  }

  // if rbc==0, homogeneous Dirichlet b.c., no need to do anything
  // if rbc==1, Neumann b.c.
  if ( rbc==1 ) {
  }

  return;
}

void TDSolveFVMMod::update_mat_chem_diff(const int nt, LinMat*& mat_loc, const int N_loc, const DType* xc_loc, const DType* xf_loc, const DType diff_coef, const int rbc ) {
  DType r2nt = data_r[nt] * data_r[nt];
  DType dx2  = dx * dx;
  
  Entry elem(1,0.0);

  // Left boundary
  {
    DType xc2 = xc_loc[0]*xc_loc[0];
    DType xf2 = xf_loc[1]*xf_loc[1];
    DType coef_diff_r = dtn * diff_coef * xf2 / dx2;
    elem = mat_loc->get(0,0); elem += coef_diff_r / xc2 / r2nt; mat_loc->set(0,0,elem);
    elem = mat_loc->get(0,1); elem -= coef_diff_r / xc2 / r2nt; mat_loc->set(0,1,elem);
  }
  // Interior points
  for ( int i = 1; i < N_loc-1; i++ ) {
    DType xc2   = xc_loc[i]*xc_loc[i];
    DType xf2_l = xf_loc[i]*xf_loc[i];
    DType xf2_r = xf_loc[i+1]*xf_loc[i+1];

    DType coef_diff_r = dtn * diff_coef * xf2_r / dx2;
    DType coef_diff_l = dtn * diff_coef * xf2_l / dx2;

    elem = mat_loc->get(i,i-1); elem -=  coef_diff_l / xc2 / r2nt;              mat_loc->set(i,i-1,elem);
    elem = mat_loc->get(i,i);   elem += (coef_diff_l+coef_diff_r) / xc2 / r2nt; mat_loc->set(i,i,  elem);
    elem = mat_loc->get(i,i+1); elem -=  coef_diff_r / xc2 / r2nt;              mat_loc->set(i,i+1,elem);
  }
  // Right boundary
  if ( rbc == 0 ) { // Homogeneous Dirichlet b.c.
    DType xf2_l = xf_loc[N_loc-1]*xf_loc[N_loc-1];
    DType xf2_r = xf_loc[N_loc]*xf_loc[N_loc];
    DType xc2   = xc_loc[N_loc-1]*xc_loc[N_loc-1];

    DType coef_diff_r = dtn * diff_coef * xf2_r / dx2;
    DType coef_diff_l = dtn * diff_coef * xf2_l / dx2;
    
    elem = mat_loc->get(N_loc-1,N_loc-2); elem -=  coef_diff_l / xc2 / r2nt;                mat_loc->set(N_loc-1,N_loc-2,elem);
    elem = mat_loc->get(N_loc-1,N_loc-1); elem += (coef_diff_l + coef_diff_r) / xc2 / r2nt; mat_loc->set(N_loc-1,N_loc-1,elem);
  }

  return;
}

// Implicit MUSCL
void TDSolveFVMMod::update_chem_muscl_be(const int fr, const int to, const int rp_nt) {
  rhs_a->reset( 0.0 );
  mat_a->reset( 0.0 );
  sol_a->reset( 0.0 );

  DType alpha_a = tdCoef->alpha_a;
  DType D       = tdCoef->D;
  DType beta_a  = tdCoef->beta_a;
  DType gamma_a = tdCoef->gamma_a;
  DType theta   = tdCoef->theta;

  DType dror = data_rp[rp_nt]/data_r[to];
  DType r2fr = data_r[fr]*data_r[fr];
  DType r2to = data_r[to]*data_r[to];

  downscale_var( *data_a[to], *data_a[to], x_a_c, Na, data_r[to] );
  compute_flux_chem_muscl_lims( to, data_a, slope_a, Na, x_a_c, x_a_f );
  upscale_var( *data_a[to], *data_a[to], x_a_c, Na, data_r[to] );

  VComp comp(1,0.0);
  Entry elem(1,0.0);

  downscale_var( *data_a[fr], *data_a[fr], x_a_c, Na, data_r[fr] );

  // RHS
  for ( int i = 0; i < Na; i++ ) {
    DType xc2 = x_a_c[i] * x_a_c[i];
    comp = xc2 * r2fr * (*data_a[fr])[i];
    if ( i < Nx ) { // G is active 
      comp += dtn * alpha_a * (*data_g[to])[i] / 
        (beta_a/theta+(*data_g[to])[i]/xc2/r2to);
    }
    comp /= xc2*r2to;
    rhs_a->set(i,comp);
  }

  // Matrix
  for ( int i = 0; i < Na; i++ ) {
    elem = 1.0 + dtn * ( gamma_a + dror );
    mat_a->set(i,i,elem);
  }

  update_mat_chem_muscl_r( to, mat_a, slope_a, Na, x_a_c, x_a_f );
  update_mat_chem_diff( to, mat_a, Na, x_a_c, x_a_f, D, 0 );

  // Initialize solution
  for ( int i = 0; i < Na; i++ ) {
    comp = (*data_a[fr])[i];
    sol_a->set(i,comp);
  }

  upscale_var( *data_a[fr], *data_a[fr], x_a_c, Na, data_r[fr] );

  linSolve_a->solve();

  // Fetch solution
  for ( int i = 0; i < Na; i++ ) {
    comp = sol_a->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_a[to])[i] = vsol[0];
  }

  upscale_var( *data_a[to], *data_a[to], x_a_c, Na, data_r[to] );

  return;
}

void TDSolveFVMMod::update_chem_opn_muscl_be(
    const DType tsol_loc, const int fr, const int to, 
    const DType coef_D, const DType coef_alpha, const DType coef_beta, 
    const DType coef_gamma, HVEC**& data_source, 
    HVEC**& data_P, HVEC**& flux_P, HVEC**& slope_P,
    const bool treatment, const bool treatment_bc, const int rp_nt) {
  rhs_p->reset( 0.0 );
  mat_p->reset( 0.0 );
  sol_p->reset( 0.0 );

  DType theta = tdCoef->theta;

  DType dror = data_rp[rp_nt]/data_r[to];
  DType r2fr = data_r[fr]*data_r[fr];
  DType r2to = data_r[to]*data_r[to];

  DType f_treat = 0.0;
  if ( treatment ) // Backward Euler
    f_treat = calc_treatment(tsol_loc+dtn);
  DType bc_treat = 0.0;
  if ( treatment_bc ) // Forward Euler Euler
    bc_treat = tdCoef->treatP1_bc * calc_treatment(tsol_loc);

  downscale_var( *data_P[fr], *data_P[fr], x_c, Nx, data_r[fr] );
  downscale_var( *data_P[to], *data_P[to], x_c, Nx, data_r[to] );

  compute_flux_chem_muscl_lims( to, data_P, slope_P, Nx, x_c, x_f, 1, bc_treat );

  VComp comp(1,0.0);
  Entry elem(1,0.0);

  // RHS
  for ( int i = 0; i < Nx; i++ ) {
    DType xc2 = x_c[i] * x_c[i];
    comp = xc2 * r2fr * (*data_P[fr])[i] + dtn * coef_alpha * (*data_source[to])[i] / 
      (coef_beta/theta+(*data_source[to])[i]/xc2/r2to) + dtn * xc2 * r2to * f_treat;
    comp /= xc2*r2to;
    rhs_p->set(i,comp);
  }
  if ( bc_treat != 0.0 ) { // Fix the non-hom Neumann condition
    comp = rhs_p->get(Nx-1);
    DType xc2 = x_c[Nx-1]*x_c[Nx-1];
    comp *= xc2 * r2to;

    comp += dtn * coef_D * x_f[Nx] * x_f[Nx] * ( - data_r[to] * bc_treat ) / dx;
    comp /= xc2 * r2to;
    rhs_p->set(Nx-1,comp);
  }

  // Matrix
  for ( int i = 0; i < Nx; i++ ) {
    elem = 1.0 + dtn * ( coef_gamma + dror );
    mat_p->set(i,i,elem);
  }

  update_mat_chem_diff( to, mat_p, Nx, x_c, x_f, coef_D );

  // Initialize solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_P[fr])[i];
    sol_p->set(i,comp);
  }

  linSolve_p->solve();

  // Fetch solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = sol_p->get(i);
    HVEC vsol(1);
    vsol = comp.getVal();
    (*data_P[to])[i] = vsol[0];
  }
  
  upscale_var( *data_P[fr], *data_P[fr], x_c, Nx, data_r[fr] );
  upscale_var( *data_P[to], *data_P[to], x_c, Nx, data_r[to] );

  return;
}

void TDSolveFVMMod::update_vel_int_opn(const int nt, const int u_nt) {

  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;
  DType b       = tdCoef->b;
  DType b_prime = tdCoef->b_prime;
  DType omega   = tdCoef->omega;

  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0

  downscale_var( *data_p2[nt], *data_p2[nt], x_c, Nx, data_r[nt] );

  // First do the integration of
  // x^2R[lambda(1+b'P2/(b+P2))*G-delta_h*H-delta_n*N]
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);

    DType loc_r = x_c[i] * data_r[nt];
    DType logistic_term = 1.0; 
    if ( tdSetup->logistic == 1 )
      logistic_term = ( 1.0 - (*data_g[nt])[i] / loc_r / loc_r );
    if ( tdSetup->logistic == 2 )
      logistic_term = max( 0.0, ( 1.0 - (*data_g[nt])[i] / loc_r / loc_r ) );

    (*data_v[nt])[i+1] = (*data_v[nt])[i] + dx * (
        loc_lambda * ( 1.0 + (b_prime * (*data_p2[nt])[i])/(b + (*data_p2[nt])[i]) ) * (*data_g[nt])[i] * logistic_term - 
        delta_h * (*data_h[nt])[i] - delta_n * (*data_n[nt])[i] -
        omega * (*data_g[nt])[i] * (*data_n[nt])[i] / loc_r / loc_r ) / data_r[nt]; // G, H, and M have not been upscaled
  }

  upscale_var( *data_p2[nt], *data_p2[nt], x_c, Nx, data_r[nt] );

  // Remove the advection part due to chemical gradient
  for ( int i = 1; i < Nx+1; i++ )
    (*data_v[nt])[i] = ((*data_v[nt])[i] - (*flux_n_u[u_nt])[i]/data_r[nt])/x_f[i]/x_f[i];

  return;
}

void TDSolveFVMMod::update_vel_int_orig(const int nt) {

  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;

  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0
  
  // First do the integration of x^2R(lambda G - delta_h H)
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    (*data_v[nt])[i+1] = (*data_v[nt])[i] + dx * ( loc_lambda * (*data_g[nt])[i] - delta_h * (*data_h[nt])[i] - delta_n * (*data_n[nt])[i] ) / data_r[nt];
  }

  downscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );

  // Remove the advection part due to chemical gradient
  for ( int i = 1; i < Nx; i++ )
    (*data_v[nt])[i] = ((*data_v[nt])[i] - 0.5*x_f[i]*x_f[i]*(*data_u[0])[i]*((*data_n[nt])[i-1]+(*data_n[nt])[i]))/x_f[i]/x_f[i];
  (*data_v[nt])[Nx] = ((*data_v[nt])[Nx] - 0.5*x_f[Nx]*x_f[Nx]*(*data_u[0])[Nx]*((*data_n[nt])[Nx-1]+tdCoef->n_bc))/x_f[Nx]/x_f[Nx];

  upscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );

  return;
}

void TDSolveFVMMod::update_vel_int_opn_orig(const int nt) {

  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;
  DType b       = tdCoef->b;
  DType b_prime = tdCoef->b_prime;
  DType alpha   = tdCoef->alpha;
  DType omega   = tdCoef->omega;

  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0

  downscale_var( *data_p2[nt], *data_p2[nt], x_c, Nx, data_r[nt] );
  
  // First do the integration of
  // x^2R[lambda(1+b'P2/(b+P2))*G-delta_h*H-delta_n*N]
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    DType loc_r = x_c[i] * data_r[nt];
    DType logistic_term = 1.0;
    if ( tdSetup->logistic == 1 )
      logistic_term = ( 1.0 - (*data_g[nt])[i] / loc_r / loc_r );
    if ( tdSetup->logistic == 2 )
      logistic_term = max( 0.0, ( 1.0 - (*data_g[nt])[i] / loc_r / loc_r ) );
    (*data_v[nt])[i+1] = (*data_v[nt])[i] + dx * ( 
        loc_lambda * ( 1.0 + (b_prime*(*data_p2[nt])[i])/(b+(*data_p2[nt])[i]) ) *  (*data_g[nt])[i] * logistic_term - 
        delta_h * (*data_h[nt])[i] - delta_n * (*data_n[nt])[i] - 
        omega * (*data_g[nt])[i] * (*data_n[nt])[i] / loc_r / loc_r ) / data_r[nt];
  }

  upscale_var( *data_p2[nt], *data_p2[nt], x_c, Nx, data_r[nt] );

  downscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );

  // Remove the advection part due to chemical gradient
  for ( int i = 1; i < Nx; i++ )
    (*data_v[nt])[i] = ((*data_v[nt])[i] - 0.5*x_f[i]*x_f[i]*(*data_u[0])[i]*((*data_n[nt])[i-1]+(*data_n[nt])[i]))/x_f[i]/x_f[i];
  (*data_v[nt])[Nx] = ((*data_v[nt])[Nx] - 0.5*x_f[Nx]*x_f[Nx]*(*data_u[0])[Nx]*((*data_n[nt])[Nx-1]+tdCoef->n_bc))/x_f[Nx]/x_f[Nx];

  upscale_var( *data_n[nt], *data_n[nt], x_c, Nx, data_r[nt] );

  return;
}

void TDSolveFVMMod::update_cells_common(const int fr, const int to, 
    const bool only_res ) {
  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

  DType dror = data_rp[0]/data_r[fr];

  DType coef_time  = only_res ? 0.0 : 1.0;
  DType coef_space = only_res ? 1.0 : dtn;

  // Source terms
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[i],p) ); 
    (*data_g[to])[i] = coef_time * (*data_g[fr])[i] + coef_space * (
        loc_lambda * (*data_g[fr])[i] - loc_mu * (*data_g[fr])[i]
        - dror * (*data_g[fr])[i] );
    (*data_h[to])[i] = coef_time * (*data_h[fr])[i] + coef_space * (
        loc_mu * (*data_g[fr])[i] - delta_h * (*data_h[fr])[i]
        - dror * (*data_h[fr])[i] );
    (*data_n[to])[i] = coef_time * (*data_n[fr])[i] + coef_space * (
        - delta_n * (*data_n[fr])[i] - dror * (*data_n[fr])[i] );
  }

  for ( int i = 0; i < Nx; i++ ) {
    (*data_g[to])[i] -= ( (*flux_g[0])[i+1] - (*flux_g[0])[i] ) * (coef_space/dx);
    (*data_h[to])[i] -= ( (*flux_h[0])[i+1] - (*flux_h[0])[i] ) * (coef_space/dx);
    (*data_n[to])[i] -= ( (*flux_n[0])[i+1] - (*flux_n[0])[i] ) * (coef_space/dx);
  }

  return;
}

void TDSolveFVMMod::update_cells_common_auto(const int fr, const int to) {
  switch ( tdSetup->scell ) {
    case -101: // Explicit first-order upwind, auto g = 1 - h - n
      update_cells_common_auto_g( fr, to );
      break;
    case -102: // Explicit first-order upwind, auto h = 1 - g - n
      update_cells_common_auto_h( fr, to );
      break;
    case -103: // Explicit first-order upwind, auto n = 1 - g - h
      update_cells_common_auto_n( fr, to );
      break;
    default:
      update_cells_common( fr, to );
      break;
  }
  return;
}

void TDSolveFVMMod::update_cells_common_auto_g(const int fr, const int to) {
  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

  DType dror = data_rp[0]/data_r[fr];

  // Source terms
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[i],p) ); 
    (*data_h[to])[i] = (*data_h[fr])[i] + dtn * (
        loc_mu * (*data_g[fr])[i] - delta_h * (*data_h[fr])[i]
        - dror * (*data_h[fr])[i] );
    (*data_n[to])[i] = (*data_n[fr])[i] + dtn * (
        - delta_n * (*data_n[fr])[i] - dror * (*data_n[fr])[i] );
  }

  for ( int i = 0; i < Nx; i++ ) {
    (*data_h[to])[i] -= ( (*flux_h[0])[i+1] - (*flux_h[0])[i] ) * (dtn/dx);
    (*data_n[to])[i] -= ( (*flux_n[0])[i+1] - (*flux_n[0])[i] ) * (dtn/dx);
  }

  downscale_var( *data_h[to], *data_h[to], x_c, Nx, data_r[to] );
  downscale_var( *data_n[to], *data_n[to], x_c, Nx, data_r[to] );

  for ( int i = 0; i < Nx; i++ )
    (*data_g[to])[i] = 1.0 - (*data_h[to])[i] - (*data_n[to])[i];

  upscale_var( *data_g[to], *data_g[to], x_c, Nx, data_r[to]);
  upscale_var( *data_h[to], *data_h[to], x_c, Nx, data_r[to]);
  upscale_var( *data_n[to], *data_n[to], x_c, Nx, data_r[to]);

  return;
}

void TDSolveFVMMod::update_cells_common_auto_h(const int fr, const int to) {
  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

  DType dror = data_rp[0]/data_r[fr];

  // Source terms
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[i],p) ); 
    (*data_g[to])[i] = (*data_g[fr])[i] + dtn * (
        loc_lambda * (*data_g[fr])[i] - loc_mu * (*data_g[fr])[i]
        - dror * (*data_g[fr])[i] );
    (*data_n[to])[i] = (*data_n[fr])[i] + dtn * (
        - delta_n * (*data_n[fr])[i] - dror * (*data_n[fr])[i] );
  }

  for ( int i = 0; i < Nx; i++ ) {
    (*data_g[to])[i] -= ( (*flux_g[0])[i+1] - (*flux_g[0])[i] ) * (dtn/dx);
    (*data_n[to])[i] -= ( (*flux_n[0])[i+1] - (*flux_n[0])[i] ) * (dtn/dx);
  }

  downscale_var( *data_g[to], *data_g[to], x_c, Nx, data_r[to] );
  downscale_var( *data_n[to], *data_n[to], x_c, Nx, data_r[to] );

  for ( int i = 0; i < Nx; i++ )
    (*data_h[to])[i] = 1.0 - (*data_g[to])[i] - (*data_n[to])[i];

  upscale_var( *data_g[to], *data_g[to], x_c, Nx, data_r[to]);
  upscale_var( *data_h[to], *data_h[to], x_c, Nx, data_r[to]);
  upscale_var( *data_n[to], *data_n[to], x_c, Nx, data_r[to]);

  return;
}

void TDSolveFVMMod::update_cells_common_auto_n(const int fr, const int to) {
  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

  DType dror = data_rp[0]/data_r[fr];

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
  }

  for ( int i = 0; i < Nx; i++ ) {
    (*data_g[to])[i] -= ( (*flux_g[0])[i+1] - (*flux_g[0])[i] ) * (dtn/dx);
    (*data_h[to])[i] -= ( (*flux_h[0])[i+1] - (*flux_h[0])[i] ) * (dtn/dx);
  }

  downscale_var( *data_g[to], *data_g[to], x_c, Nx, data_r[to] );
  downscale_var( *data_h[to], *data_h[to], x_c, Nx, data_r[to] );

  for ( int i = 0; i < Nx; i++ )
    (*data_n[to])[i] = 1.0 - (*data_g[to])[i] - (*data_h[to])[i];

  upscale_var( *data_g[to], *data_g[to], x_c, Nx, data_r[to]);
  upscale_var( *data_h[to], *data_h[to], x_c, Nx, data_r[to]);
  upscale_var( *data_n[to], *data_n[to], x_c, Nx, data_r[to]);

  return;
}

DType TDSolveFVMMod::calc_n_bc( const DType tsol_loc ) {
  DType n_bc = tdCoef->n_bc;
  if ( tdSetup->prob_id == -2 ) {
    DType v_par = tdSetup->init_u_par;
    n_bc = 1.0 - std::exp(-3.0*v_par*tsol_loc) * ( 1.0-n_bc );
  }
  return n_bc;
}

DType TDSolveFVMMod::LimFun_PPM_BC(DType dval_c, DType dval_r) {
  if ( fabs(dval_c)<TOL_ZERO && fabs(dval_r)<TOL_ZERO )
    return 1.0;
  if ( dval_c * dval_r <= 0.0 || 3.0*fabs(dval_r) <= 8.0*fabs(dval_c) )
    return 0.0;
  else
    return min(1.0, 5.0/(fabs(dval_r)/fabs(dval_c)-1.0));
}

void TDSolveFVMMod::LimFun_PPM_INT(DType alpha1, DType alpha2, DType& lim_l, DType& lim_r) {
  if ( alpha2 > alpha1 + TOL_ZERO ) {
    lim_l = 0.0; lim_r = 0.0;
  } else if ( alpha1 < 1.0 ) {
    lim_l = 1.0; lim_r = alpha1;
  } else if ( alpha2 > 1.0 ) {
    lim_l = 1.0/alpha2; lim_r = 1.0;
  } else {
    lim_l = 1.0; lim_r = 1.0;
  }
  return;
}

void TDSolveFVMMod::update_chem_rhs( const int fr, const int to, const DType shift_dt ) {
  switch ( tdSetup->prob_id ) {
    case 0: // PDGF
      switch ( tdSetup->schem ) {
        case 0:
          update_chem_upw_rhs( fr, to );
          break;
        case 1:
        case 2:
          update_chem_muscl_rhs( fr, to );
          break;
      }
      break;
    case 1: // OPN
      switch ( tdSetup->schem ) {
        case 0:
          update_chem_upw_rhs( fr, to );
          //update_chem_opn_upw_rhs( tsol+shift_dt_fr*dtn, fr, to, tdCoef->D_1, tdCoef->alpha_1, 
          //    tdCoef->beta_1, tdCoef->gamma_1, data_n, data_p1, tdCoef->treatP1, tdCoef->treatP1BC );
          //update_chem_opn_upw_rhs( tsol+shift_dt_fr*dtn, fr, to, tdCoef->D_2, 
          //    tdCoef->alpha_2, tdCoef->beta_2, tdCoef->gamma_2, data_g, data_p2 );
          break;
        case 1:
        case 2:
          update_chem_muscl_rhs( fr, to );
          //update_chem_opn_muscl_rhs( tsol+shift_dt_fr*dtn, fr, to, tdCoef->D_1, tdCoef->alpha_1, 
          //    tdCoef->beta_1, tdCoef->gamma_1, data_n, data_p1, tdCoef->treatP1, tdCoef->treatP1BC );
          //update_chem_opn_muscl_rhs( tsol+shift_dt_fr*dtn, fr, to, tdCoef->D_2, 
          //    tdCoef->alpha_2, tdCoef->beta_2, tdCoef->gamma_2, data_g, data_p2 );
          break;
      }
  }
  return;
}

void TDSolveFVMMod::update_chem_upw_rhs( const int fr, const int to ) {
  DType alpha_a = tdCoef->alpha_a;
  DType D       = tdCoef->D;
  DType beta_a  = tdCoef->beta_a;
  DType gamma_a = tdCoef->gamma_a;
  DType theta   = tdCoef->theta;

  DType dror = (*data_v[fr])[Nx]/data_r[fr];

  for ( int i = 0; i < Na; i++ ) {
    DType loc_r_fr = x_a_c[i] * data_r[fr];
    DType loc_r_to = x_a_c[i] * data_r[to];
    //comp = loc_r_fr * loc_r_fr * (*data_a[fr])[i];
    //if ( i < Nx ) { // G is active 
    //  comp += dtn * alpha_a * (*data_g[to])[i] / 
    //    (beta_a/theta+(*data_g[to])[i]/loc_r_to/loc_r_to);
    //}
    //comp /= loc_r_to*loc_r_to;
    //rhs_a->set(i,comp);
  }
}

void TDSolveFVMMod::update_chem_muscl_rhs( const int fr, const int to ) {
}

void TDSolveFVMMod::update_sol_radau_exp( const int rkns, 
    RKButcher& rkButcher, const int fr, const int to) {
  assert( to == rkns+1 );
  *data_g[to] = *data_g[fr];
  *data_h[to] = *data_h[fr];
  *data_n[to] = *data_n[fr];

  int Ns = rkButcher.getNstage();
  for ( int j = 0; j < Ns; j++ ) {
    DType rk_coef = rkButcher.get(rkns,j) * dtn;
    *data_g[to] += rk_coef * (*resd_g[j+1]);
    *data_h[to] += rk_coef * (*resd_h[j+1]);
    *data_n[to] += rk_coef * (*resd_n[j+1]);
  }

  return;
}

void TDSolveFVMMod::update_sol_radau_imp_upw( const int rkns, 
    RKButcher& rkButcher, const int fr, const int to) {
  assert( to == rkns+1 );
  int Ns = rkButcher.getNstage();
  //DType rk_coef = rkButcher.get(rkns,rkns);
  DType rk_coef = rkButcher.stage(rkns);

  DType lambda = tdCoef->lambda;
  DType q      = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

  DType dror = data_rp[to]/data_r[to];
  DType r2fr = data_r[fr]*data_r[fr];
  DType r2to = data_r[to]*data_r[to];

  DType n_bc = calc_n_bc( tsol+rkButcher.stage(rkns)*dtn );

  // Compute the limiters and store them in lim_?, as well as 
  // unlimited interpolated values, and store them in flux_?
  // We won't use flux_? here, though
  compute_flux_r_ppm_lims( to );

  downscale_var(*data_g[fr],*data_g[fr],x_c,Nx,data_r[fr]);
  downscale_var(*data_h[fr],*data_h[fr],x_c,Nx,data_r[fr]);
  downscale_var(*data_n[fr],*data_n[fr],x_c,Nx,data_r[fr]);

  downscale_var(*data_g[to],*data_g[to],x_c,Nx,data_r[to]);
  downscale_var(*data_h[to],*data_h[to],x_c,Nx,data_r[to]);
  downscale_var(*data_n[to],*data_n[to],x_c,Nx,data_r[to]);

  VComp comp(1,0.0);
  Entry elem(1,0.0);

  // First data_g
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  // RHS
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_g[fr])[i] * r2fr;
    DType xc2 = x_c[i] * x_c[i];
    DType loc_res = (*resd_g[rkns+1])[i];
    for ( int j = 0; j < Ns; j++ )
      if ( j != rkns )
        comp += rkButcher.get(rkns,j) * dtn * ((*resd_g[j+1])[i]-loc_res) / xc2;
    rhs->set(i,comp);
  }

  // Matrix
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[i],p) ); 
    elem = r2to - rk_coef * dtn * r2to * ( loc_lambda - loc_mu - dror );
    mat->set(i,i,elem);
  }

  update_mat_upw_v( rk_coef, to, data_v, to, false, to );
  update_mat_ppm_r( rk_coef, to, to );

  // Initialize solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_g[to])[i];
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
    comp = (*data_h[fr])[i] * r2fr + rk_coef * dtn * loc_mu * (*data_g[to])[i] * r2to;
    DType xc2 = x_c[i] * x_c[i];
    DType loc_res = (*resd_h[rkns+1])[i];
    for ( int j = 0; j < Ns; j++ )
      if ( j!= rkns )
        comp += rkButcher.get(rkns,j) * dtn * ((*resd_h[j+1])[i]-loc_res) / xc2;
    rhs->set(i,comp);
  }

  // Matrix
  for ( int i = 0; i < Nx; i++ ) {
    elem = r2to - rk_coef * dtn * r2to * ( - delta_h - dror );
    mat->set(i,i,elem);
  }

  update_mat_upw_v( rk_coef, to, data_v, to, false, to );
  update_mat_ppm_r( rk_coef, to, to );

  // Initialize solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_h[to])[i];
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

  // RHS
  for ( int i = 0; i < Nx-1; i++ ) {
    comp = (*data_n[fr])[i] * r2fr;
    DType xc2 = x_c[i] * x_c[i];
    DType loc_res = (*resd_n[rkns+1])[i];
    for ( int j = 0; j < Ns; j++ )
      if ( j != rkns )
        comp += rkButcher.get(rkns,j) * dtn * ((*resd_n[j+1])[i]-loc_res) / xc2;
    rhs->set(i,comp);
  }
  // Right boundary
  {
    DType xc2 = x_c[Nx-1]*x_c[Nx-1];
    DType xf2 = x_f[Nx]*x_f[Nx];

    DType loc_u = (*data_u[0])[Nx]/data_r[to];
    DType coef_u_r = 0.0;
    if ( loc_u < 0.0 )
      coef_u_r += rk_coef * (dtn/dx) * xf2 * r2to * loc_u;
    comp = (*data_n[fr])[Nx-1] * r2fr - coef_u_r * n_bc / xc2;
    DType loc_res = (*resd_n[rkns+1])[Nx-1];
    for ( int j = 0; j < Ns; j++ )
      if ( j != rkns )
        comp += rkButcher.get(rkns,j) * dtn * ((*resd_n[j+1])[Nx-1]-loc_res) / xc2;
    rhs->set(Nx-1,comp);
  }

  // Matrix
  for ( int i = 0; i < Nx; i++ ) {
    elem = r2to - rk_coef * dtn * r2to * ( - delta_n - dror );
    mat->set(i,i,elem);
  }

  update_mat_upw_v( rk_coef, to, data_v, to, false, to );
  update_mat_ppm_r( rk_coef, to, to );
  update_mat_upw_v( rk_coef, to, data_u, to, true, to );

  // Initialize solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_n[to])[i];
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

  upscale_var(*data_g[fr],*data_g[fr],x_c,Nx,data_r[fr]);
  upscale_var(*data_h[fr],*data_h[fr],x_c,Nx,data_r[fr]);
  upscale_var(*data_n[fr],*data_n[fr],x_c,Nx,data_r[fr]);

  upscale_var(*data_g[to],*data_g[to],x_c,Nx,data_r[to]);
  upscale_var(*data_h[to],*data_h[to],x_c,Nx,data_r[to]);
  upscale_var(*data_n[to],*data_n[to],x_c,Nx,data_r[to]);

  // // Debug -- recompute as explicit
  // *data_g[to] = *data_g[fr];
  // *data_n[to] = *data_n[fr];

  // Ns = rkButcher.getNstage();
  // for ( int j = 0; j < Ns; j++ ) {
  //   DType rk_coef = rkButcher.get(rkns,j) * dtn;
  //   *data_g[to] += rk_coef * (*resd_g[j+1]);
  //   *data_n[to] += rk_coef * (*resd_n[j+1]);
  // }

  return;
}

void TDSolveFVMMod::update_sol_radau_imp_muscl( const int rkns, 
    RKButcher& rkButcher, const int fr, const int to) {
  assert( to == rkns+1 );
  int Ns = rkButcher.getNstage();
  //DType rk_coef = rkButcher.get(rkns,rkns);
  DType rk_coef = rkButcher.stage(rkns);

  DType lambda = tdCoef->lambda;
  DType q      = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

  DType dror = data_rp[to]/data_r[to];
  DType r2fr = data_r[fr]*data_r[fr];
  DType r2to = data_r[to]*data_r[to];

  DType n_bc = calc_n_bc( tsol+rkButcher.stage(rkns)*dtn );

  compute_flux_r_ppm_lims( to );

  downscale_var(*data_g[fr],*data_g[fr],x_c,Nx,data_r[fr]);
  downscale_var(*data_h[fr],*data_h[fr],x_c,Nx,data_r[fr]);
  downscale_var(*data_n[fr],*data_n[fr],x_c,Nx,data_r[fr]);

  downscale_var(*data_g[to],*data_g[to],x_c,Nx,data_r[to]);
  downscale_var(*data_h[to],*data_h[to],x_c,Nx,data_r[to]);
  downscale_var(*data_n[to],*data_n[to],x_c,Nx,data_r[to]);

  // For the modified FVM, use primitive variables to compute muscl limiters
  compute_flux_muscl_lims( to );
  
  VComp comp(1,0.0);
  Entry elem(1,0.0);

  // First data_g
  rhs->reset( 0.0 );
  mat->reset( 0.0 );
  sol->reset( 0.0 );

  // RHS
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_g[fr])[i] * r2fr;
    DType xc2 = x_c[i] * x_c[i];
    DType loc_res = (*resd_g[rkns+1])[i];
    for ( int j = 0; j < Ns; j++ )
      if ( j != rkns )
        comp += rkButcher.get(rkns,j) * dtn * ((*resd_g[j+1])[i]-loc_res) / xc2;
    rhs->set(i,comp);
  }

  // Matrix
  for ( int i = 0; i < Nx; i++ ) {
    DType loc_lambda = ( q == 0.0 ) ? lambda : lambda * pow(x_c[i],q);
    DType loc_mu = ( p == 0.0 ) ? mu * ( mu_coef - mu_slp ) : mu * ( mu_coef - mu_slp * pow(x_c[i],p) ); 
    elem = r2to - rk_coef * dtn * r2to * ( loc_lambda - loc_mu - dror );
    mat->set(i,i,elem);
  }

  update_mat_muscl_v( rk_coef, to, data_v, to, false, to );
  update_mat_ppm_r( rk_coef, to, to );

  // Initialize solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_g[to])[i];
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
    comp = (*data_h[fr])[i] * r2fr + rk_coef * dtn * loc_mu * (*data_g[to])[i] * r2to;
    DType xc2 = x_c[i] * x_c[i];
    DType loc_res = (*resd_h[rkns+1])[i];
    for ( int j = 0; j < Ns; j++ )
      if ( j!= rkns )
        comp += rkButcher.get(rkns,j) * dtn * ((*resd_h[j+1])[i]-loc_res) / xc2;
    rhs->set(i,comp);
  }

  // Matrix
  for ( int i = 0; i < Nx; i++ ) {
    elem = r2to - rk_coef * dtn * r2to * ( - delta_h - dror );
    mat->set(i,i,elem);
  }

  update_mat_muscl_v( rk_coef, to, data_v, to, false, to );
  update_mat_ppm_r( rk_coef, to, to );

  // Initialize solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_h[to])[i];
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

  // RHS
  for ( int i = 0; i < Nx-1; i++ ) {
    comp = (*data_n[fr])[i] * r2fr;
    DType xc2 = x_c[i] * x_c[i];
    DType loc_res = (*resd_n[rkns+1])[i];
    for ( int j = 0; j < Ns; j++ )
      if ( j != rkns )
        comp += rkButcher.get(rkns,j) * dtn * ((*resd_n[j+1])[i]-loc_res) / xc2;
    rhs->set(i,comp);
  }
  // Right boundary
  {
    DType xc2 = x_c[Nx-1]*x_c[Nx-1];
    DType xf2 = x_f[Nx]*x_f[Nx];

    DType loc_u = (*data_u[0])[Nx]/data_r[to];
    DType coef_u_r = 0.0;
    if ( loc_u < 0.0 )
      coef_u_r += rk_coef * (dtn/dx) * xf2 * r2to * loc_u;
    comp = (*data_n[fr])[Nx-1] * r2fr - coef_u_r * n_bc / xc2;
    DType loc_res = (*resd_n[rkns+1])[Nx-1];
    for ( int j = 0; j < Ns; j++ )
      if ( j != rkns )
        comp += rkButcher.get(rkns,j) * dtn * ((*resd_n[j+1])[Nx-1]-loc_res) / xc2;
    rhs->set(Nx-1,comp);
  }

  // Matrix
  for ( int i = 0; i < Nx; i++ ) {
    elem = r2to - rk_coef * dtn * r2to * ( - delta_n - dror );
    mat->set(i,i,elem);
  }

  update_mat_muscl_v( rk_coef, to, data_v, to, false, to );
  update_mat_ppm_r( rk_coef, to, to );
  update_mat_muscl_v( rk_coef, to, data_u, to, true, to );

  // Initialize solution
  for ( int i = 0; i < Nx; i++ ) {
    comp = (*data_n[to])[i];
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

  upscale_var(*data_g[fr],*data_g[fr],x_c,Nx,data_r[fr]);
  upscale_var(*data_h[fr],*data_h[fr],x_c,Nx,data_r[fr]);
  upscale_var(*data_n[fr],*data_n[fr],x_c,Nx,data_r[fr]);

  upscale_var(*data_g[to],*data_g[to],x_c,Nx,data_r[to]);
  upscale_var(*data_h[to],*data_h[to],x_c,Nx,data_r[to]);
  upscale_var(*data_n[to],*data_n[to],x_c,Nx,data_r[to]);

  return;
}


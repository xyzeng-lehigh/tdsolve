#include "TDSolveFVMMod.h"
#include "Alloc_Utils.h"
#include "TDInput.h"
#include "Alg_Solve.h"
#include "Alg_Limiters.h"
using namespace std;

TDSolveFVMMod::TDSolveFVMMod(TDInput* tdInput_ptr) :
  TDSolveFVM(tdInput_ptr) {

  alloc_array_full( lim_g, Nx, 2 );
  alloc_array_full( lim_h, Nx, 2 );
  alloc_array_full( lim_n, Nx, 2 ); // M
  alloc_array_full( lim_s, Nx, 2 ); // Sync

  alloc_array_full( flux_g_v, Nx+1, 1 );
  alloc_array_full( flux_g_r, Nx+1, 1 );
  alloc_array_full( flux_h_v, Nx+1, 1 );
  alloc_array_full( flux_h_r, Nx+1, 1 );
  alloc_array_full( flux_n_v, Nx+1, 1 ); // or M
  alloc_array_full( flux_n_r, Nx+1, 1 ); // or M
  alloc_array_full( flux_n_u, Nx+1, 1 ); // or M

  alloc_array_full( data_u, Nx+1, 1 );

}

TDSolveFVMMod::~TDSolveFVMMod(void) {
  erase_array( lim_g, 2 );
  erase_array( lim_h, 2 );
  erase_array( lim_n, 2 );
  erase_array( lim_s, 2 );
  erase_array( flux_g_v, 1 );
  erase_array( flux_g_r, 1 );
  erase_array( flux_h_v, 1 );
  erase_array( flux_h_r, 1 );
  erase_array( flux_n_v, 1 ); // or M
  erase_array( flux_n_r, 1 ); // or M
  erase_array( flux_n_u, 1 ); // or M
  erase_array( data_u, 1 );
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

  *data_u[0] = 0.0;
  int u_pow = tdSetup->init_u_pow;
  if ( u_pow != 0 )
    for ( int i = 1; i <= Nx; i++ )
      (*data_u[0])[i] = -std::pow(i*dx,u_pow);

  // Set all the right hand sides to zero:
  tdCoef->lambda  = 0.0;
  tdCoef->mu      = 0.0;
  tdCoef->delta_h = 0.0;
  tdCoef->delta_n = 0.0;
  tdCoef->omega   = 0.0;
  tdSetup->logistic = 0;

  this->update_vel( 0 );

  FILE* ot_file = fopen( tdInput->radiusFile, "w" );
  fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\t%10.20f\n", 
      tsol, data_r[0], (*data_v[0])[Nx], -(*data_u[0])[Nx] );
  fclose(ot_file);

  return;
}

void TDSolveFVMMod::print_custom_vel_u( FILE* ot_file, const int nt ){
  fprintf( ot_file, "%10.20f\t%10.20f\t%10.20f\t%10.20f\n", tsol, data_r[nt], (*data_v[nt])[Nx], -(*data_u[0])[Nx] );
}

void TDSolveFVMMod::update_sol(const int fr, const int to, const int shift_dt) {

  // Use data_v[fr], which is computed in previous iteration
  update_rad( fr, to );

  // The fluxes are already computed in previous iteration
  update_cells_common( fr, to );

  if ( tdSetup->prob_id == 1 ) // OPN
    update_cells_opn_add( fr, to );

  // Chemicals
  switch ( tdSetup->prob_id ) {
    case 0: // PDGF
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
    case 1: // OPN
      switch ( tdSetup->schem ) {
        case 0: // Implicit first-order upwind for advection and cd for diffusion
          // Update A
          update_chem_upw_be( fr, to );
          // Update P1 and P2
          update_chem_opn_upw_be( tsol+shift_dt*dtn, fr, to, tdCoef->D_1, tdCoef->alpha_1, 
              tdCoef->beta_1, tdCoef->gamma_1, data_n, data_p1, tdCoef->treatP1 );
          update_chem_opn_upw_be( tsol+shift_dt*dtn, fr, to, tdCoef->D_2, 
              tdCoef->alpha_2, tdCoef->beta_2, tdCoef->gamma_2, data_g, data_p2 );
          break;
        case 1: // Explicit second-order muscl for advection and implicit cd for diffusion
          // Update A
          update_chem_muscl_imex( fr, to ); 
          // Update P1 and P2
          update_chem_opn_muscl_imex( tsol+shift_dt*dtn, fr, to, tdCoef->D_1, tdCoef->alpha_1, 
              tdCoef->beta_1, tdCoef->gamma_1, data_n, data_p1, flux_p1, slope_p1, tdCoef->treatP1 );
          update_chem_opn_muscl_imex( tsol+shift_dt*dtn, fr, to, tdCoef->D_2, tdCoef->alpha_2, 
              tdCoef->beta_2, tdCoef->gamma_2, data_g, data_p2, flux_p2, slope_p2 );
          break;
        default:
          std::cout << "Warning: The spatial scheme for chem not defined, "
            << "using upw/backward on the primitive equation instead\n" << std::endl;
          // Update A
          update_chem_upw_be( fr, to );
          // Update P1 and P2
          update_chem_opn_upw_be( tsol+shift_dt*dtn, fr, to, tdCoef->D_1, tdCoef->alpha_1, 
              tdCoef->beta_1, tdCoef->gamma_1, data_n, data_p1, tdCoef->treatP1 );
          update_chem_opn_upw_be( tsol+shift_dt*dtn, fr, to, tdCoef->D_2, tdCoef->alpha_2, 
              tdCoef->beta_2, tdCoef->gamma_2, data_g, data_p2 );
          break;
      }
      break;
  }

  // Will also compute fluxes for the current time step
  update_vel( to );

  return;
}

void TDSolveFVMMod::update_vel(const int nt) {

  // Compute the velocity u using the gradient
  compute_vel_u( nt );

  // The fluxes that does not need V
  compute_cells_fluxes_u( nt );

  switch ( tdSetup->prob_id ) {
    case 0:  // PDGF
    case -1: // Test
      update_vel_int(nt);
      break;
    case 1: // OPN
      update_vel_int_opn(nt);
      break;
  }

  // The fluxes that need V
  compute_cells_fluxes_v( nt );

  // The fluxes associated with R'
  compute_cells_fluxes_r( nt );

  return;
}

void TDSolveFVMMod::compute_vel_u(const int nt) {
  switch ( tdSetup->prob_id ) {
    case -1: // Test, u is independent of x
      break;
    case 0: // PDGF
    case 1: // OPN
      // u = (alpha/R^2)(dA/dx)
      compute_vel_u_chem( nt );
      break;
  }
  return;
}

void TDSolveFVMMod::compute_vel_u_chem(const int nt) {
  *data_u[0] = 0.0;

  DType alpha = tdCoef->alpha;
  DType aor  = alpha/data_r[nt]/data_r[nt];

  for ( int i = 1; i <= Nx; i++ )
    (*data_u[0])[i] = aor * (((*data_a[nt])[i]-(*data_a[nt])[i-1]) / dx);

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_u(const int nt) {
  switch ( tdSetup->scell ) {
    case 0: // Explicit first-order upwind
      compute_cells_fluxes_u_upw( nt );
      break;
    case 2: // Explicit second-order MUSCL
      compute_cells_fluxes_u_muscl( nt );
      break;
    default:
      std::cout << "Warning: The flux scheme for cell not defined, "
        << "using upwind instead\n" << std::endl;
      compute_cells_fluxes_u_upw( nt );
      break;
  }

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_v(const int nt) {
  switch ( tdSetup->scell ) {
    case 0: // Explicit first-order upwind
      compute_cells_fluxes_v_upw( nt );
      break;
    case 2: // Explicit second-order MUSCL
      compute_cells_fluxes_v_muscl( nt );
      break;
    default:
      std::cout << "Warning: The flux scheme for cell not defined, "
        << "using upwind instead\n" << std::endl;
      compute_cells_fluxes_v_upw( nt );
      break;
  }

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_u_upw(const int nt) {

  downscale_var(*data_n[nt],*data_n[nt],x_c,Nx,data_r[nt]);

  // Fluxes
  { // Left boundary
    (*flux_n_u[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_r  = x_f[i] * data_r[nt];
    DType loc_r2 = loc_r * loc_r;

    DType loc_vel_u = (*data_u[0])[i]/data_r[nt];
    int i_idx_u = (loc_vel_u>=0.0) ? i-1 : i;
    (*flux_n_u[0])[i] = loc_r2 * loc_vel_u * (*data_n[nt])[i_idx_u];
  }
  { // Right boundary
    DType loc_vel_u = (*data_u[0])[Nx]/data_r[nt];
    DType loc_r = x_f[Nx]*data_r[nt];
    (*flux_n_u[0])[Nx] = loc_r * loc_r * loc_vel_u * 
      ((loc_vel_u>=0.0) ? (*data_n[nt])[Nx-1] : tdCoef->n_bc);
  }

  upscale_var(*data_n[nt],*data_n[nt],x_c,Nx,data_r[nt]);

  return;
}

void TDSolveFVMMod::compute_cells_fluxes_u_muscl(const int nt) {

  downscale_var(*data_g[nt],*data_g[nt],x_c,Nx,data_r[nt]);
  downscale_var(*data_h[nt],*data_h[nt],x_c,Nx,data_r[nt]);
  downscale_var(*data_n[nt],*data_n[nt],x_c,Nx,data_r[nt]);

  // Limiters for the velocity associated with u and with v
  { // Left boundary
    (*lim_g[0])[0] = 0.0;
    (*lim_h[0])[0] = 0.0;
    (*lim_n[0])[0] = 0.0;
    (*lim_s[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx-1; i++ ) {
    DType loc_lim_g = LimFun_Minmod(
        (*data_g[nt])[i] - (*data_g[nt])[i-1],
        (*data_g[nt])[i+1] - (*data_g[nt])[i] );
    DType loc_lim_h = LimFun_Minmod(
        (*data_h[nt])[i] - (*data_h[nt])[i-1],
        (*data_h[nt])[i+1] - (*data_h[nt])[i] );
    DType loc_lim_n = LimFun_Minmod(
        (*data_n[nt])[i] - (*data_n[nt])[i-1],
        (*data_n[nt])[i+1] - (*data_n[nt])[i] );
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
    (*lim_s[0])[i] = min( min(loc_lim_g,loc_lim_h),
                          min(loc_lim_n,loc_lim_theta) );
  }
  { // Right boundary
    (*lim_g[0])[Nx-1] = 0.0;
    (*lim_h[0])[Nx-1] = 0.0;
    (*lim_n[0])[Nx-1] = 0.0;
    (*lim_s[0])[Nx-1] = 0.0;
  }

  // Slopes
  { // Left boundary
    (*slope_g[0])[0] = 0.0;
    (*slope_h[0])[0] = 0.0;
    (*slope_n[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx-1; i++ ) {
    (*slope_g[0])[i] = (*lim_s[0])[i] * 
      ((*data_g[nt])[i+1]-(*data_g[nt])[i])/dx;
    (*slope_h[0])[i] = (*lim_s[0])[i] * 
      ((*data_h[nt])[i+1]-(*data_h[nt])[i])/dx;
    (*slope_n[0])[i] = (*lim_s[0])[i] * 
      ((*data_n[nt])[i+1]-(*data_n[nt])[i])/dx;
  }
  { // Right boundary
    (*slope_g[0])[Nx-1] = 0.0;
    (*slope_h[0])[Nx-1] = 0.0;
    (*slope_n[0])[Nx-1] = 0.0;
  }

  // Fluxes
  { // Left boundary
    (*flux_n_u[0])[0] = 0.0;
  }
  for ( int i = 1; i < Nx; i++ ) {
    DType loc_r  = x_f[i] * data_r[nt];
    DType loc_r2 = loc_r * loc_r;

    DType loc_vel_u = (*data_u[0])[i]/data_r[nt];
    if ( loc_vel_u >= 0.0 )
      (*flux_n_u[0])[i] = loc_r2 * loc_vel_u * (
          (*data_n[nt])[i-1]+0.5*dx*(*slope_n[0])[i-1] );
    else
      (*flux_n_u[0])[i] = loc_r2 * loc_vel_u * (
          (*data_n[nt])[i]-0.5*dx*(*slope_n[0])[i] );
  }
  { // Right boundary
    DType loc_vel_u = (*data_u[0])[Nx]/data_r[nt];
    DType loc_r = x_f[Nx]*data_r[nt];
    (*flux_n_u[0])[Nx] = loc_r * loc_r * loc_vel_u * 
      ((loc_vel_u>=0.0) ? (*data_n[nt])[Nx-1]+0.5*dx*(*slope_n[0])[Nx-1] : tdCoef->n_bc );
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

void TDSolveFVMMod::compute_cells_fluxes_v_muscl(const int nt) {

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
          (*data_g[nt])[i-1]+0.5*dx*(*slope_g[0])[i-1] );
      (*flux_h_v[0])[i] = loc_r2 * loc_vel_v * (
          (*data_h[nt])[i-1]+0.5*dx*(*slope_h[0])[i-1] );
      (*flux_n_v[0])[i] = loc_r2 * loc_vel_v * (
          (*data_n[nt])[i-1]+0.5*dx*(*slope_n[0])[i-1] );
    }
    else {
      (*flux_g_v[0])[i] = loc_r2 * loc_vel_v * (
          (*data_g[nt])[i]-0.5*dx*(*slope_g[0])[i] );
      (*flux_h_v[0])[i] = loc_r2 * loc_vel_v * (
          (*data_h[nt])[i]-0.5*dx*(*slope_h[0])[i] );
      (*flux_n_v[0])[i] = loc_r2 * loc_vel_v * (
          (*data_n[nt])[i]-0.5*dx*(*slope_n[0])[i] );
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

void TDSolveFVMMod::compute_cells_fluxes_r(const int nt) {
  DType dror = (*data_v[nt])[Nx]/data_r[nt];
  dror = dror * ( 1.0 + 0.5 * dtn * dror );

  const DType sot = 7.0/12.0;
  const DType oot = 1.0/12.0;

  // x^2R^2N -> x^3R^2N, etc.
  upscale_var(*data_g[nt],*data_g[nt],x_c,Nx,1.0,1);
  upscale_var(*data_h[nt],*data_h[nt],x_c,Nx,1.0,1);
  upscale_var(*data_n[nt],*data_n[nt],x_c,Nx,1.0,1);

  // Step one -- use flux_n, etc, to store the surface values
  DType rbc_coef   = pow((x_c[Nx-1]+dx)/x_c[Nx-1],3);
  DType rbc_coef_2 = pow((x_c[Nx-1]+2.0*dx)/x_c[Nx-1],3);
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
  // Near right boundary (constant extrapolation in the primitive variables)
  {
    (*flux_g[0])[Nx-1] = sot*((*data_g[nt])[Nx-2]+(*data_g[nt])[Nx-1]) 
      - oot*((*data_g[nt])[Nx-3]+rbc_coef*(*data_g[nt])[Nx-1]);
    (*flux_h[0])[Nx-1] = sot*((*data_h[nt])[Nx-2]+(*data_h[nt])[Nx-1]) 
      - oot*((*data_h[nt])[Nx-3]+rbc_coef*(*data_h[nt])[Nx-1]);
    (*flux_n[0])[Nx-1] = sot*((*data_n[nt])[Nx-2]+(*data_n[nt])[Nx-1]) 
      - oot*((*data_n[nt])[Nx-3]+rbc_coef*(*data_n[nt])[Nx-1]);

    (*flux_g[0])[Nx] = sot*((*data_g[nt])[Nx-1]+rbc_coef*(*data_g[nt])[Nx-1])
      - oot*((*data_g[nt])[Nx-2]+rbc_coef_2*(*data_g[nt])[Nx-1]);
    (*flux_h[0])[Nx] = sot*((*data_h[nt])[Nx-1]+rbc_coef*(*data_h[nt])[Nx-1])
      - oot*((*data_h[nt])[Nx-2]+rbc_coef_2*(*data_h[nt])[Nx-1]);
    (*flux_n[0])[Nx] = sot*((*data_n[nt])[Nx-1]+rbc_coef*(*data_n[nt])[Nx-1])
      - oot*((*data_n[nt])[Nx-2]+rbc_coef_2*(*data_n[nt])[Nx-1]);
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
        LimFun_PPM_INT(alpha1, alpha2, (*lim_s[0])[i], (*lim_n[1])[i]);
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
        LimFun_PPM_INT(alpha1, alpha2, (*lim_s[0])[i], (*lim_n[1])[i]);
      }
    }
  }

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
    (*flux_n[0])[i] = (*flux_n_v[0])[i]+(*flux_n_r[0])[i]+(*flux_n_u[0])[i];
  }

  return;
}

void TDSolveFVMMod::update_vel_int(const int nt) {

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
    (*data_v[nt])[i] = ((*data_v[nt])[i] - (*flux_n_u[0])[i]/data_r[nt])/x_f[i]/x_f[i];

  return;
}

void TDSolveFVMMod::update_vel_int_opn(const int nt) {

  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;
  DType b       = tdCoef->b;
  DType b_prime = tdCoef->b_prime;
  DType omega   = tdCoef->omega;

  (*data_v[nt])[0] = 0.0; // b.c. for velocity at x=0

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

  // Remove the advection part due to chemical gradient
  for ( int i = 1; i < Nx+1; i++ )
    (*data_v[nt])[i] = ((*data_v[nt])[i] - (*flux_n_u[0])[i]/data_r[nt])/x_f[i]/x_f[i];

  return;
}

void TDSolveFVMMod::update_cells_common(const int fr, const int to) {
  DType lambda  = tdCoef->lambda;
  DType q       = tdCoef->q;
  DType delta_h = tdCoef->delta_h;
  DType delta_n = tdCoef->delta_n;

  // mu = mu(mu_coef-mu_slp*x^p)
  DType mu      = tdCoef->mu;
  DType mu_slp  = tdCoef->mu_slp;
  DType mu_coef = tdCoef->mu_coef;
  DType p       = tdCoef->p;

  DType dror = (*data_v[fr])[Nx]/data_r[fr];
  dror = dror * ( 1.0 + 0.5 * dtn * dror );

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

  for ( int i = 0; i < Nx; i++ ) {
    (*data_g[to])[i] -= ( (*flux_g[0])[i+1] - (*flux_g[0])[i] ) * (dtn/dx);
    (*data_h[to])[i] -= ( (*flux_h[0])[i+1] - (*flux_h[0])[i] ) * (dtn/dx);
    (*data_n[to])[i] -= ( (*flux_n[0])[i+1] - (*flux_n[0])[i] ) * (dtn/dx);
  }

  return;
}

DType TDSolveFVMMod::LimFun_PPM_BC(DType dval_c, DType dval_r) {
  if ( fabs(dval_c)<TOL_ZERO && fabs(dval_r)<TOL_ZERO )
    return 1.0;
  if ( dval_c * dval_r <= 0.0 ||
       3.0*fabs(dval_r) <= 8.0*fabs(dval_c) )
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

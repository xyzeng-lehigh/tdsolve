#include <iostream>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>

#include "Data_Type.h"
#include "Comput_Choices.h"
#include "TDInput.h"
#include "TDSolve.h"
#include "TDSolveFDM.h"
#include "TDSolveFVM.h"
#include "TDSolveFVMMod.h"
#include "TDSolveSacado.h"
#include "RK_Coef.h"
using namespace std;

char inputFile[40];
char* program_name;

void printlogo(void);
void parse(int*, char**);
void usage(void);
void print_help(void);

int main(int argc, char** argv) {
  printlogo();
  parse(&argc,argv);

  TDInput tdInput;
  tdInput.setInputFile((const char*)inputFile);
  tdInput.parse();

  return 0;

  TDSolve* tdSolve;
  bool classical_tmeth = false;
  switch ( tdInput.tdSetup.frame ) {
    case 0:
      tdSolve = new TDSolveFDM(&tdInput);
      break;
    case 1:
      tdSolve = new TDSolveFVM(&tdInput);
      break;
    case 2:
      tdSolve = new TDSolveFVMMod(&tdInput);
      break;
    case 3:
      tdSolve = new TDSolveSacado(&tdInput);
      classical_tmeth = true;
      break;
    default:
      tdSolve = new TDSolveSacado(&tdInput);
      classical_tmeth = true;
      break;
  }

  //tdSolve->getTDInput()->display_parameter();
  //tdSolve->getTDInput()->display_option();
  //tdSolve->getTDInput()->display_toggle();
  
  tdSolve->alloc();

  tdSolve->setup();
  tdSolve->set_ic();
  tdSolve->calc_dt(0);
  tdSolve->diagnostic(0);

  if (classical_tmeth) { // A cleaner version
    while ( tdSolve->nextIter() ) {
      switch ( tdInput.tdSetup.tmeth ) {
        case 0:  // First-order explicit/implicit Euler
          tdSolve->update_sol(0,1);
          break;
        default:
          tdSolve->update_sol(0,1);
          break;
      }
      tdSolve->output(0);
      tdSolve->copy_data(1,0);
      tdSolve->calc_dt(0);
      tdSolve->diagnostic(0);
    }
    tdSolve->update_vel(0,1,0.0,true);
    tdSolve->output(0);
  }
  else { // A little messy with the old code
    while ( tdSolve->nextIter() ) {
      switch ( tdInput.tdSetup.tmeth ) {
        case 0: // First-order segregated
          {
            bool vcoup_conv;
            while (1) {
              vcoup_conv = tdSolve->update_sol(0,1); // BE or FE
              if ( !vcoup_conv && tdInput.tdSetup.auto_red_dt ) {
                tdSolve->reduce_dtn( 0.5 );
                tdSolve->diagnostic(0);
                continue;
              }
              else
                break;
            }
          }
          break;
        case 1: // Second-order segregated
          if ( tdInput.tdSetup.vCoup == 0 ) { // TVD RK2
            tdSolve->update_sol(0,1);
            tdSolve->update_sol(1,2,1.0,2.0);
            tdSolve->average_sol(0.5,0,0.5,2,1);
          }
          else { // DIRK2
            // Butcher table: (alpha = 1-sqrt(2)/2)
            //
            // alpha   |   alpha      0
            // 1       |   1-alpha    alpha
            // ----------------------------
            //         |   1-alpha    alpha
            //
            //=============================== Standard Formulation =============
            //
            // Stage 1:
            //   S_1 ~ t_n + alpha * dtn
            //   S_1 = S_n + alpha * dtn * F(S_1)
            // Stage 2:
            //   S_2 ~ t_n + dtn
            //   S_2 = S_n + (1-alpha) * dtn * F(S_1) + alpha * dtn * F(S_1)
            // Finally:
            //   S_{n+1} = S_2
            //
            //=============================== Modified Formulation =============
            //
            // Stage 1:
            //   S_1' = S_n + alpha * dtn * F(S_1') = S_1
            // Stage 2:
            //   omega = (1-alpha)/alpha
            //   S_2' = (1-omega) * S_n + omega * S_1' ~ t_n + (1-alpha) * dtn
            //   S_2' = S_n + (1-alpha) * dtn * F(S_1')
            //   S_3' ~ t_n + (1-alpha) * dtn + alpha * dtn
            //   S_3' = S_2' + alpha * dtn * F(S_3')
            // Finally:
            //   S_{n+1} = S_3'
            
            bool vcoup_conv;
            while (1) {
              DType alpha = 1.0-0.5*std::sqrt(2.0);
              tdSolve->scale_dtn( alpha );
              vcoup_conv = tdSolve->update_sol(0,1);
              if ( !vcoup_conv && tdInput.tdSetup.auto_red_dt ) {
                tdSolve->scale_dtn( 1.0/alpha );
                tdSolve->reduce_dtn( 0.5 );
                tdSolve->diagnostic(0);
                continue;
              }
              tdSolve->average_sol(2.0-1.0/alpha,0,1.0/alpha-1.0,1,2);
              vcoup_conv = tdSolve->update_sol(2,1,1.0/alpha-1.0,1.0/alpha);
              if ( !vcoup_conv && tdInput.tdSetup.auto_red_dt ) {
                tdSolve->scale_dtn( 1.0/alpha );
                tdSolve->reduce_dtn( 0.5 );
                tdSolve->diagnostic(0);
                continue;
              }
              else {
                tdSolve->scale_dtn( 1.0/alpha );
                break;
              }
            }
          }
          break;
        case 2: // Second-order Radau collocation RK
          {
            // Butcher table
            //
            // 1/3    |    5/12   -1/12
            // 1      |    3/4    1/4
            // ------------------------
            //        |    3/4    1/4
            //
            //=============================== Standard Formulation =============
            //
            //   S_1 ~ t_1 = t_n + 1/3 * dtn
            //   S_2 ~ t_2 = t_n + dtn
            //
            //   S_1 = S_n + 5/12 * dtn * F(S_1) - 1/12 * dtn * F(S_2)
            //   S_2 = S_n + 3/4  * dtn * F(S_1) + 1/4  * dtn * F(S_2)
            // 
            // Finally
            //   S_{n+1} = S_2
            //
            //=============================== Modified Formulation =============
            //
            //   Predictor S_2^(0) ( FE, BE, etc )
            //
            // Given S_2^(k-1), k >= 1:
            //   Solve for S_1^(k):
            //     S_1^(k) = S_n + 5/12 * dtn * F(S_1^(k)) 
            //                   - 1/12 * dtn * F(S_2^(k-1))
            //   Solve for S_2^(k):
            //     S_2^(k) = S_n + 3/4 * dtn * F(S_1^(k-1))
            //                   + 1/4 * dtn * F(S_2^(k))
            //
            //   Until converge
            //   
            //   Data usage:
            //   data_*[0] - solutions at t_n
            //   data_*[1] - solutions at t_1
            //   data_*[2] - solutions at t_2

            int Ns = 2;
            RKButcher rkButcher(Ns);
            rkButcher.get(0,0) = 5.0/12.0;
            rkButcher.get(0,1) = -1.0/12.0;
            rkButcher.get(1,0) = 3.0/4.0;
            rkButcher.get(1,1) = 1.0/4.0;

            bool sol_conv, red_dt;

            if ( !(tdInput.tdSetup.use_pmc) ) { // Direct solve for a larger system, freeze limiters though
              // This will setup the matrices if not already done so.
              tdSolve->setup_rk_lin( rkButcher );  

              // Predictor corrector only on the radius!
              while (1) {
                for ( int i = 1; i <= Ns; i++ )
                  tdSolve->copy_data(0,i); // Only data_v and data_r are needed
                // Predict data_rp and data_r at all stages
                tdSolve->update_rad_radau(rkButcher); 
                int n_iter = -1;
                tdSolve->check_conv_radius(rkButcher, n_iter, sol_conv, red_dt, true);

                while (!sol_conv && !red_dt) {
                  // Freeze limiters, and solve for cell species given data_rp and data_r
                  if (tdSolve->update_rk_sol(rkButcher)) { // Newton iteration converged
                    // Update data_rp and data_r according to the new data_v
                    tdSolve->update_rad_radau(rkButcher);
                    // Check if the data_r has converged
                    tdSolve->check_conv_radius(rkButcher, n_iter, sol_conv, red_dt, true);
                  }
                  else {
                    sol_conv = false; red_dt = true;
                  }
                }

                // Either converged or need to reduce dt
                if (sol_conv) // Converged
                  break;
                else if (tdInput.tdSetup.auto_red_dt) { // Reduce dt
                  tdSolve->reduce_dtn( 0.5 );
                  tdSolve->diagnostic( 0 );
                  continue;
                }
                else { // Not converged but no reduce of dt
                  tdSolve->print_warning("Not converged");
                  break;
                }
              }

              // Legacy issue, the solution at t_{n+1} needs to be 
              // stored at data_*[1]
              tdSolve->copy_data(2,1);
            }
            else { // Predictor/multi-corrector
              // Each loop corresponds to a dt (could be reduced!)
              // data_*[0], slope_*[0], and flux_n_u[0] should NOT be modified!
              while (1) {
                // Copy data_g, data_h, data_n, and data_v
                for ( int i = 1; i <= Ns; i++ ) {
                  tdSolve->copy_data(0,i); // Stage i (0) at data_*[i]
                  // flux_n_u and slope_*[nt] should have been computed 
                  // in set_ic or previous iteration
                  tdSolve->copy_aux(0,i);
                }

                // Step 3 - compute data_r[1] and data_r[2]
                // Compute data_rp and data_r at all stages
                tdSolve->update_rad_radau(rkButcher);
                // Step 4 - compute the residuals for S_1 and S_2
                for ( int i = 1; i <= Ns; i++ )
                  tdSolve->compute_rhs(rkButcher.stage(i-1), i); // Stored in resd_*[i]
                // Step 5 - evaluate full residual to check convergence
                int  n_iter = -1;
                tdSolve->check_conv_radau(rkButcher, n_iter, sol_conv, red_dt);

                while (!sol_conv && !red_dt) {
                  // Step 1 - Implicit solve for new iterates
                  for ( int i = 1; i <= Ns; i++ )
                    tdSolve->update_sol(i-1,rkButcher,0,i); // Update S_i
                  // Step 2 - update the velocities
                  //        - and flux_n_u and lim_muscl
                  for ( int i = 1; i <= Ns; i++ )
                    tdSolve->update_vel(rkButcher.stage(i-1), i); // V_i
                  // Step 3 - compute data_r[1] and data_r[2]
                  tdSolve->update_rad_radau(rkButcher);
                  // Step 4 - compute the residuals for S_1 and S_2
                  for ( int i = 1; i <= Ns; i++ )
                    tdSolve->compute_rhs(rkButcher.stage(i-1), i); // S_i
                  // Step 5 - compute non-cell species if needed
                  // Step 6 - evaluate full residual to check convergence
                  tdSolve->check_conv_radau(rkButcher, n_iter, sol_conv, red_dt);
                }

                // Either converged or need to reduce dt
                if (sol_conv) // Converged
                  break;
                else if ( tdInput.tdSetup.auto_red_dt ) { // Reduce dt
                  tdSolve->reduce_dtn( 0.5 );
                  tdSolve->diagnostic( 0 );
                  continue;
                }
                else { // Not converged but no reduce of dt
                  tdSolve->print_warning("Not converged");
                  break;
                }
              }
              // Legacy issue, the solution at t_{n+1} needs to be 
              // stored at data_*[1]
              tdSolve->copy_data(2,1);
            }

          }
          break;
        case 100: // Test first-order fully backward Euler
          {
            int Ns = 1;
            RKButcher rkButcher(Ns);
            rkButcher.get(0,0) = 1.0;
            
            bool sol_conv, red_dt;

            tdSolve->setup_rk_lin( rkButcher );
            
            while (1) {
              tdSolve->copy_data(0,1);
              tdSolve->update_rad_radau(rkButcher);
              int n_iter = -1;
              tdSolve->check_conv_radius(rkButcher,n_iter,sol_conv,red_dt);
              while (!sol_conv && !red_dt) {
                if (tdSolve->update_rk_sol(rkButcher)) {
                  tdSolve->update_rad_radau(rkButcher);
                  tdSolve->check_conv_radius(rkButcher,n_iter,sol_conv,red_dt);
                }
                else {
                  sol_conv = false; red_dt = true;
                }
              }
              if (sol_conv)
                break;
              else if (tdInput.tdSetup.auto_red_dt) {
                tdSolve->reduce_dtn( 0.5 );
                tdSolve->diagnostic( 0 );
                continue;
              }
              else {
                tdSolve->print_warning("Not converged");
                break;
              }
            }
          }
          break;
        default: // First-order
          tdSolve->update_sol(0,1);
          break;
      }
      tdSolve->output(0); // This is the solution before updating
      tdSolve->copy_data(1,0);
      tdSolve->calc_dt(0);
      tdSolve->diagnostic(0);
    }
    tdSolve->update_vel(0,1,0.0,true);
    tdSolve->output(0);
  }

  return 1;
}

void printlogo(void) {
  printf("\n");
  printf("*************************\n");
  printf("* Tumor Dynamics Solver *\n");
  printf("*************************\n");
  printf("\n");
  return;
}

void parse(int* argc, char** argv) {
  program_name = argv[0];
  while ((*argc>1)&&(argv[1][0]=='-')) {
    switch (argv[1][1]) {
      case 'i':
        strcpy(inputFile,argv[2]);
        break;
      case 'h':
        print_help();
        break;
      default:
        printf("Bad option %s\n",argv[1]);
        usage();
        break;
    }
    ++argv; ++argv;
    --(*argc); --(*argc);
  }
  if ((*argc)>1) {
    printf("Unknown or too many parameters.\n");
    usage();
  }
  return;
}

void usage(void) {
  printf("Usage is %s [-flag] [value]\n",program_name);
  printf("Options\n");
  printf("  -i    Specify the input file name\n");
  printf("  -h    Print a brief documentation of options available\n");
  printf("\n");
  exit(1);
}

void print_help(void) {
  printf("In the input file: [name] [{default value}] [description]\n");
  printf(" Nx           {100}    Number of cells (number of grid points is Nx+1)\n");
  printf(" A_Multiplier {2}      Multiplier of cells for chemicals\n");
  printf(" cfl          {0.9}    The Courant number\n");
  printf(" ic           {0}      An integer denoting specific initial values (problems)\n");
  printf("                       Examplar input files are in the folder \"probs\"\n");
  printf("                       0 -- A simple test to check the code\n");
  printf(" lbc          {0}      The boundary condition for cells at the left boundary\n");
  printf("                       0 -- No-flux boundary\n");
  printf("                       1 -- ODE-based boundary condition\n");
  printf(" space_cell   {0}      Scheme to discretize the equations for cell numbers in space\n");
  printf("                       0 -- Explicit 1st-order upwind\n");
  printf("                       1 -- Implicit 1st-order upwind\n");
  printf("                       2 -- Partially implicit 4th-order Pade for each cell number\n");
  printf(" space_chem   {0}      Scheme to discretize the equation for chemicals in space\n");
  printf("                       0 -- Implicit 1st-order upwind\n");
  printf("                       1 -- Implicit 4th-order Pade\n"); 
  printf(" vel          {0}      Scheme to update the velocity\n");
  printf("                       0 -- 2nd-order integration form\n");
  printf("                       1 -- 1st-order upwind\n");
  printf("                       2 -- 2nd-order BDF\n");
  printf("                       3 -- 4th-order Compact differencing\n");
  printf(" time_methd   {1}      Scheme for time-integration\n");
  printf("                       0 -- FE/BE\n");
  printf("                       1 -- RK2\n");
  printf("                       2 -- Radau RK2\n");
  //printf(" rbc          {0}      Boundary condition at the tumor front\n");
  //printf("                       0 -- Homogeneous Neumann b.c. for G and H\n");
  //printf("                       1 -- ODE based (no b.c. at the front)\n");
}

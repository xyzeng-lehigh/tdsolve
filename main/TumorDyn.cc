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

  TDSolve* tdSolve;
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
    default:
      tdSolve = new TDSolveFVM(&tdInput);
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

  while ( tdSolve->nextIter() ) {
    switch ( tdInput.tdSetup.tmeth ) {
      case 0: // First-order segregated
        {
          bool vcoup_conv;
          while (1) {
            vcoup_conv = tdSolve->update_sol(0,1); // BE or FE
            if ( !vcoup_conv ) {
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
          bool vcoup_conv;
          while (1) {
            DType alpha = 1.0-0.5*std::sqrt(2.0);
            tdSolve->scale_dtn( alpha );
            vcoup_conv = tdSolve->update_sol(0,1);
            if ( !vcoup_conv ) {
              tdSolve->scale_dtn( 1.0/alpha );
              tdSolve->reduce_dtn( 0.5 );
              tdSolve->diagnostic(0);
              continue;
            }
            tdSolve->average_sol(2.0-1.0/alpha,0,1.0/alpha-1.0,1,2);
            vcoup_conv = tdSolve->update_sol(2,1,1.0/alpha-1.0,1.0/alpha);
            if ( !vcoup_conv ) {
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
  //printf(" rbc          {0}      Boundary condition at the tumor front\n");
  //printf("                       0 -- Homogeneous Neumann b.c. for G and H\n");
  //printf("                       1 -- ODE based (no b.c. at the front)\n");
}

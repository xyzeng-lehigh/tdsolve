#include <iostream>
#include <cstdio>
#include <cstring>
#include <ctime>

#include "Data_Type.h"
#include "Comput_Choices.h"
#include "TDInput.h"
#include "TDSolve.h"
using namespace std;

char inputFile[40];
char* program_name;

void printlogo(void);
void parse(int*, char**);
void usage(void);

int main(int argc, char** argv) {
  printlogo();
  parse(&argc,argv);

  TDInput tdInput;
  tdInput.setInputFile((const char*)inputFile);
  tdInput.parse();

  TDSolve tdSolve(&tdInput);

  tdSolve.setup();
  tdSolve.set_ic();
  tdSolve.calc_dt(0);
  tdSolve.output(0);
  tdSolve.diagnostic(0);

  while ( tdSolve.nextIter() ) {
    tdSolve.update_sol(0,1);
    tdSolve.update_sol(1,2);
    tdSolve.average_sol(0.5,0,0.5,2,1);
    tdSolve.copy_data(1,0);
    tdSolve.calc_dt(0);
    tdSolve.diagnostic(0);
    tdSolve.output(0);
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
  printf("\n");
  exit(1);
}


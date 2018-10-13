#include <iostream>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <fstream>
#include <vector>
#include <map>
#include <set>

#include "Data_Type.h"
#include "Comput_Choices.h"
#include "TDInput.h"
#include "TDSolve_API.h"
#include "Fun_Estimator.h"
using namespace std;

char  inputFileBase[40];
char  inputFilePar[40];
char  inputFileData[40];
char  inputFileRefs[40];
bool  use_output_ref;
char  outputFile[40];
char  fitMethod[40];
char  fitType[40];
char  probName[40];
char* program_name;

void printlogo(void);
void parse(int*, char**);
void usage(void);

char delims[40];
char delims_var[40];
char delims_nospace[40];
void parse_pars(void);
void parse_orefs(void);
void create_factories(std::vector<TDFactory>&, const TDInput&);
bool increment_var_array(int*&,const int*,const int);
void add_output( TDFactory& );
void add_error( TDFactory& );
std::map< std::string, std::set<DType> > pars;
std::set< std::string >                  keys;
std::map< std::string, DType >           o_refs;
std::set< std::string >                  okeys;

int main(int argc, char** argv) {
  printlogo();
  parse(&argc,argv);

  pars.clear();
  keys.clear();
  o_refs.clear();
  okeys.clear();

  // Step 1: Create the input decks for all computations
  
  // Setup the baseline input file
  TDInput tdInputBase;
  tdInputBase.setInputFile((const char*)inputFileBase);
  tdInputBase.parse();

  // Create an array of input decks
  std::strcpy(delims_var,    " ,=\t");
  std::strcpy(delims,        (const char*)tdInputBase.delims);
  std::strcpy(delims_nospace,(const char*)tdInputBase.delims_nospace);
  parse_pars();
  parse_orefs();
  std::vector<TDFactory> factories;
  factories.clear();
  create_factories(factories, tdInputBase);

  // Step 2: Run the computations in parallel and append the new results
  //         to the end of the provided data file
#pragma omp parallel for
  for ( int n = 0; n < factories.size(); n++ ) {
    bool valid = factories[n].run();
#pragma omp critical
    if ( valid )
      add_output( factories[n] );
    else
      add_error( factories[n] );
  }

  // Step 3: Parametric study
  FunEstimator* funEst = new FunEstimator(fitMethod,fitType);
  funEst->setNumInput(keys.size());
  funEst->setNumOutput(okeys.size());
  funEst->getReferenceInput(tdInputBase,keys);
  funEst->getReferenceOutput(o_refs,okeys);
  funEst->readData(inputFileData);
  funEst->run();
  funEst->outputPar(outputFile);
  delete funEst;

  return 0;

}

void printlogo(void) {
  printf("\n");
  printf("*****************************\n");
  printf("* Parametric Analysis Suite *\n");
  printf("*****************************\n");
  printf("\n");
  return;
}

void parse(int* argc, char** argv) {
  program_name = argv[0];
  use_output_ref = false;
  while ((*argc>1)&&(argv[1][0]=='-')) {
    switch (argv[1][1]) {
      case 'd':
        strcpy(inputFileData,argv[2]);
        break;
      case 'i':
        strcpy(inputFileBase,argv[2]);
        break;
      case 'm':
        strcpy(fitMethod,argv[2]);
        break;
      case 'n':
        strcpy(probName,argv[2]);
        break;
      case 'o':
        strcpy(outputFile,argv[2]);
        break;
      case 'p':
        strcpy(inputFilePar,argv[2]);
        break;
      case 't':
        strcpy(fitType,argv[2]);
        break;
      case 'r':
        use_output_ref = true;
        strcpy(inputFileRefs,argv[2]);
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
  printf("  -d    Specify the data set file name\n");
  printf("  -i    Specify the baseline input file name\n");
  printf("  -j    Specify the number of threads\n");
  printf("  -m    Fitting method (LS0,LS1,LS2,LS3,LS4,LS5)\n");
  printf("  -n    Problem name -- decide which solver to use and the outputs\n");
  printf("          PDGF - PDGF/Tterm, Nperc\n");
  printf("  -o    Specify the output file to store fitting information\n");
  printf("  -p    Specify the parametric input file name\n");
  printf("  -r    Specify the output reference values\n");
  printf("  -t    Fitting type (Poly,ExpPoly)\n");
  printf("\n");
  exit(1);
}

void parse_pars(void) {
  std::ifstream fin;
  std::string file_line_s;
  fin.open(inputFilePar);
  char* name;
  char* value;

  while (getline(fin,file_line_s)) {
    // This ignore a short line
    if (file_line_s.length()<2)
      continue;
    // This indicates the comment lines
    if ((file_line_s[0]=='/')&&(file_line_s[1]=='/'))
      continue;
    // This ignore a non-comment line with no separation tokens
    if (file_line_s.find_first_of((const char*)delims_nospace)==string::npos)
      continue;
    name = strtok((char*)file_line_s.c_str(),delims);
    keys.insert(name);
    std::map< std::string, std::set<DType> >::iterator it = pars.find( name );
    if ( it == pars.end() ) {
      std::set<DType> empty_var;
      empty_var.clear();
      pars[name] = empty_var;
      it = pars.find( name );
    }
    std::set<DType>& cur_var = it->second;
    while ( (value = strtok(NULL,delims_var)) )
      cur_var.insert( atof(value) );
  }

  return;
}

void parse_orefs(void) {
  if ( !use_output_ref ) return;

  std::ifstream fin;
  std::string file_line_s;
  fin.open(inputFileRefs);
  char* name;
  char* value;

  while (getline(fin,file_line_s)) {
    // This ignore a short line
    if (file_line_s.length()<2)
      continue;
    // This indicates the comment lines
    if ((file_line_s[0]=='/')&&(file_line_s[1]=='/'))
      continue;
    // This ignore a non-comment line with no separation tokens
    if (file_line_s.find_first_of((const char*)delims_nospace)==string::npos)
      continue;
    name  = strtok((char*)file_line_s.c_str(),delims);
    value = strtok(NULL,delims);
    o_refs[name] = atof(value);
  }

  return;
}

void create_factories(std::vector<TDFactory>& factories, const TDInput& input_base) {
  int* cur_idx = new int[keys.size()];
  int* max_idx = new int[keys.size()];

  int cur_i = 0;
  for ( std::set<std::string>::iterator it = keys.begin(); it != keys.end(); it++, cur_i++ ) {
    cur_idx[cur_i] = 0;
    max_idx[cur_i] = pars[*it].size();
  }

  do {
    TDInput tdInput = input_base;
    std::map< std::string, DType > fac_pars;
    fac_pars.clear();
    int cur_i = 0;
    for ( std::set<std::string>::iterator it = keys.begin(); it != keys.end(); it++, cur_i++ ) {
      std::map< std::string, std::set<DType> >::iterator it_par = pars.find(*it);
      if ( it_par == pars.end() ) {
        fprintf(stderr, "Error: There is a key %s that is not mapped in pars.\n", (*it).c_str() );
        exit(1);
      }
      std::set<DType>& loc_pars = it_par->second;
      std::set<DType>::iterator it_val = loc_pars.begin();
      std::advance(it_val,cur_idx[cur_i]);
      tdInput.set_parameter(*it, *it_val);
      fac_pars[*it] = *it_val;
    }

    TDFactory tdFactory(tdInput,probName);
    tdFactory.pars = fac_pars;
    if ( okeys.size() != 0 )
      okeys = tdFactory.okeys;
    factories.push_back(tdFactory);
  }
  while ( increment_var_array(cur_idx,max_idx,keys.size()) );

  delete[] cur_idx;
  delete[] max_idx;
  return;
}

bool increment_var_array(int*& cur_idx, const int* max_idx, const int size) {
  int last_inc=-1;
  for ( int i = size-1; i>=0; i-- ) {
    if ( cur_idx[i] < max_idx[i]-1 ) {
      last_inc = i;
      break;
    }
  }

  if ( last_inc == -1 )
    return false;
  else {
    cur_idx[last_inc]++;
    for ( int i = last_inc+1; i < size; i++ )
      cur_idx[i] = 0;
    return true;
  }

}

void add_output( TDFactory& tdFac ) {
  FILE* data_file = fopen(inputFileData,"a");
  for ( std::set<std::string>::iterator it = keys.begin(); it != keys.end(); it++ )
    fprintf( data_file, "%lf\t", tdFac.pars[*it] );
  fprintf( data_file, "\n" );
  for ( std::set<std::string>::iterator it = okeys.begin(); it != okeys.end(); it++ )
    fprintf( data_file, "%lf\t", tdFac.sols[*it] );
  fprintf( data_file, "\n" );
  fclose(data_file);
  return;
}

void add_error( TDFactory& tdFac ) {
  FILE* data_file = fopen(inputFileData,"a");
  for ( std::set<std::string>::iterator it = keys.begin(); it != keys.end(); it++ )
    fprintf( data_file, "%lf\t", tdFac.pars[*it] );
  fprintf( data_file, "\n" );
  fprintf( data_file, "Unsuccessful run!\n" );
  fclose(data_file);
  return;
}



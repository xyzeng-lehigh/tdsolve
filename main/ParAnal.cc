#include <iostream>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <cmath>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <filesystem>

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
bool  perform_fitting;
bool  use_output_ref;
bool  use_input_ref;
char  outputFile[40];
char  fitMethod[40];
char  probName[40];
char* program_name;

bool  simulDataKeep;
char  simulDataDrct[40];
char  simulDataBase[40];

void printlogo(void);
void parse(int*, char**);
void usage(void);

char delims[40];
char delims_var[40];
char delims_nospace[40];
void parse_pars(void);
void parse_orefs(void);
void create_factories(std::vector<TDFactory*>&, const TDInput&);
bool increment_var_array(int*&,const int*,const int);
void add_output( TDFactory& );
void add_error( TDFactory& );
std::map< std::string, std::set<DType> > pars;
std::set< std::string >                  keys;
std::set< std::string >                  okeys;
std::map< std::string, DType >           o_refs;
std::map< std::string, std::string >     o_funs;

int main(int argc, char** argv) {
  printlogo();
  parse(&argc,argv);

  pars.clear();
  keys.clear();
  okeys.clear();
  o_refs.clear();
  o_funs.clear();

  // Step 1: Create the input decks for all computations
  
  // Setup the baseline input file
  TDInput tdInputBase;
  tdInputBase.setInputFile((const char*)inputFileBase);
  tdInputBase.parse();
  tdInputBase.tdSetup.silence_screen = true;
  tdInputBase.tdSetup.silence_output = !simulDataKeep;

  // Create an array of input decks
  std::strcpy(delims_var,    " ,=\t");
  std::strcpy(delims,        (const char*)tdInputBase.delims);
  std::strcpy(delims_nospace,(const char*)tdInputBase.delims_nospace);
  parse_pars();
  parse_orefs();
  std::vector<TDFactory*> factories;
  factories.clear();
  create_factories(factories, tdInputBase);

  // Step 2: Run the computations in parallel and append the new results
  //         to the end of the provided data file
#pragma omp parallel for
  for ( int n = 0; n < factories.size(); n++ ) {
    factories[n]->alloc();
    bool valid = factories[n]->run();
#pragma omp critical
    {
      cout << " - Complete simulation of factory " << n+1 << " out of " << factories.size() 
           << ", writing results ... ";
      if ( valid )
        add_output( *(factories[n]) );
      else
        add_error( *(factories[n]) );
      cout << "Completed." << std::endl;
      delete factories[n]->getTDInput();
      delete factories[n];
    }
  }

  if ( !perform_fitting ) return 0;

  // Step 3: Parametric study
  FunEstimator* funEst = new FunEstimator(fitMethod);
  funEst->setInputKeys(keys);
  funEst->setOutputKeys(okeys);
  funEst->getReferenceInput(use_input_ref, tdInputBase);
  funEst->getReferenceOutput(o_refs);
  funEst->getFittingMethod(o_funs);
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
  perform_fitting = false;
  use_output_ref = false;
  use_input_ref  = true;
  std::strcpy(fitMethod,"LS2");
  std::strcpy(probName,"PDGF");
  std::strcpy(inputFileData,"regression.dat");
  simulDataKeep  = false;
  std::strcpy(simulDataDrct,"data");
  std::strcpy(simulDataBase,"simul");
  while ((*argc>1)&&(argv[1][0]=='-')) {
    bool bad_option = false;
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
        perform_fitting = true;
        strcpy(outputFile,argv[2]);
        break;
      case 'p':
        strcpy(inputFilePar,argv[2]);
        break;
      case 'r':
        use_output_ref = true;
        strcpy(inputFileRefs,argv[2]);
        break;
      case 'R':
        use_input_ref = ((std::atoi(argv[2])) != 0);
        break;
      case 'S':
        switch (argv[1][2]) {
          case 'b':
            strcpy(simulDataBase,argv[2]);
            break;
          case 'd':
            strcpy(simulDataDrct,argv[2]);
            break;
          case 'k':
            simulDataKeep = ((std::atoi(argv[2])) != 0);
            break;
          default:
            bad_option = true;
            break;
        }
        break;
      default:
        bad_option = true;
        break;
    }
    if (bad_option) {
      printf("Bad option %s\n",argv[1]);
      usage();
    }
    ++argv; ++argv;
    --(*argc); --(*argc);
  }
  if ((*argc)>1) {
    printf("Unknown or too many parameters.\n");
    usage();
  }
  if (simulDataKeep) {
    if (std::filesystem::exists(simulDataDrct)) {
      std::filesystem::remove_all(simulDataDrct);
    }
    std::filesystem::create_directory(simulDataDrct);
  }
  return;
}

void usage(void) {
  printf("Usage is %s [-flag] [value]\n",program_name);
  printf("Options\n");
  printf("  -d    Specify the data set file name\n");
  printf("  -i    Specify the baseline input file name\n");
  printf("  -m    Fitting method (LS0D,LS1D,LS2D,LS3D,LS4D,LS5D) - degree of poly\n");
  printf("                       (LS0T,LS1T,LS2T,LS3T,LS4T,LS5T) - degree in each var\n");
  printf("  -n    Problem name -- decide which solver to use and the outputs\n");
  printf("          PDGF - PDGF/Tterm, Gperc, Nperc\n");
  printf("          OPN  - OPN/Tterm, Gperc, Nperc\n");
  printf("          Potassium - Potassium/B2peak, B2vol, B3vol, Pvol, Vvol\n");
  printf("  -o    Specify the output file to store fitting information\n");
  printf("  -p    Specify the parametric input file name\n");
  printf("  -r    Specify the output reference values\n");
  printf("  -R    Toggle of using reference value for inputs\n");
  printf("  -Sb   Base name for sub-folder of single-simulation files\n");
  printf("  -Sd   Directory to put simulation files\n");
  printf("  -Sk   Option to keep simulation files\n");
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

  fin.close();

  return;
}

void parse_orefs(void) {
  if ( !use_output_ref ) return;

  std::ifstream fin;
  std::string file_line_s;
  fin.open(inputFileRefs);
  char* name;
  char* value;

  // Format
  //  name = numeric.value
  // or
  //  name = string
  // Here numeric.value cannot be zero, otherwise it will
  // be interpreted as a string

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
    DType trial_value = atof(value);
    if ( trial_value == 0.0 )
      o_funs[name] = value;
    else
      o_refs[name] = trial_value;
  }

  fin.close();

  return;
}

void create_factories(std::vector<TDFactory*>& factories, const TDInput& input_base) {
  int* cur_idx = new int[keys.size()];
  int* max_idx = new int[keys.size()];

  int cur_i = 0;
  int tot_count = 1;
  for ( std::set<std::string>::iterator it = keys.begin(); it != keys.end(); it++, cur_i++ ) {
    cur_idx[cur_i] = 0;
    max_idx[cur_i] = pars[*it].size();
    tot_count *= pars[*it].size();
  }

  if ( tot_count == 0 ) {
    TDInput* tdInput = new TDInput(input_base);
    TDFactory tdFactory(tdInput,probName);
    okeys = tdFactory.okeys;
    delete tdInput;
    delete[] cur_idx;
    delete[] max_idx;
    return;
  }

  cout << "Start creating factories ..." << std::endl;

  int f_count = 0;
  do {
    TDInput* tdInput = new TDInput(input_base);
    string subfolder = string(simulDataDrct) + "/" + simulDataBase;
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
      tdInput->set_parameter(*it, *it_val);
      fac_pars[*it] = *it_val;

      if ( std::strcmp( probName, "PDGF" ) == 0 )
        tdInput->set_option("prob_id",0);
      if ( std::strcmp( probName, "OPN" ) == 0 )
        tdInput->set_option("prob_id",1);

      std::stringstream stream;
      if ( *it_val > 1000 )
        stream << std::scientific << std::setprecision(2) << *it_val;
      else
        stream << std::fixed << std::setprecision(2) << *it_val;
      subfolder += "_" + *it + stream.str();
    }

    if (simulDataKeep) {
      std::filesystem::create_directory(subfolder);
      string simul_sol_file = subfolder+"/sol.dat";
      string simul_rad_file = subfolder+"/rad.dat";
      std::strcpy(tdInput->outputFile,simul_sol_file.c_str());
      std::strcpy(tdInput->radiusFile,simul_rad_file.c_str());
    }
    
    TDFactory* tdFactory = new TDFactory(tdInput,probName);
    tdFactory->pars = fac_pars;
    //tdFactory->getTDInput()->display_parameter();
    //tdFactory->getTDInput()->display_option();
    //tdFactory->getTDInput()->display_toggle();
    if ( okeys.size() == 0 )
      okeys = tdFactory->okeys;
    factories.push_back(tdFactory);
    f_count++;
    cout << " - Complete factory creation " << f_count << " out of " << tot_count << std::endl;
  }
  while ( increment_var_array(cur_idx,max_idx,keys.size()) );

  cout << "Complete creating factories, start simulations ..." << std::endl;

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



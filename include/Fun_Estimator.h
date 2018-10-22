#ifndef FUN_ESTIMATOR_H
#define FUN_ESTIMATOR_H

#include <vector>
#include <utility>
#include <cstring>
#include <set>
#include <map>
#include "Data_Type.h"
#include "Alloc_Utils.h"

#ifndef FLOAT_ALG
extern "C" void dposv_(char *UPLO, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, int *INFO);
#else
extern "C" void sposv_(char *UPLO, int *N, int *NRHS, float *A,  int *LDA, double *B, int *LDB, int *INFO);
#endif

class TDInput;

class FunEstimator {
  public:
    FunEstimator(char*);

    enum FitFuns {
      FUN_POLY = 0, FUN_EXP_POLY, FUN_NUM
    };

    void readData(const char*);

    void run(void) {
      preprocess();
      (this->*estimate)(); 
    };

    void outputPar(const char* ofile) {
      write(ofile); };

    void setInputKeys(const std::set<std::string> keys_in ) {
      ikeys  = keys_in;
      //for ( std::set<std::string>::iterator it = ikeys.begin();
      //      it != ikeys.end(); it++ )
      //  cout << " (" << *it << ") ";
      //cout << std::endl;
      return;
    }

    void setOutputKeys(const std::set<std::string> keys_in ) {
      okeys   = keys_in;
      //for ( std::set<std::string>::iterator it = okeys.begin();
      //      it != okeys.end(); it++ )
      //  cout << " (" << *it << ") ";
      //cout << std::endl;
      return;
    }

    void getReferenceInput(const bool,TDInput&);
    void getReferenceOutput(std::map<std::string,DType>&);
    void getFittingMethod(std::map<std::string,std::string>&);

  private:

    // A vector of ikeys.size() unsigned powers
    typedef std::vector<unsigned> PolyTerm;
    std::vector<PolyTerm> basis;

    // The first one is input, the second one is output
    typedef std::pair< std::vector<DType>, std::vector<DType> > VarPair;
    std::vector<VarPair> data;

    // Reference values for the parameters
    // Used for normalizing the inputs
    std::set<std::string>        ikeys;
    std::map<std::string, DType> irefs;
    // Reference values for the outputs, could be zero (no reference)
    std::set<std::string>          okeys;
    std::map<std::string, DType>   orefs;
    std::map<std::string, FitFuns> ofuns;
    std::map<std::string, std::vector<DType> > fcoef;

    typedef void (FunEstimator::*funRun)(void);

    funRun estimate;

    // Preprocessor
    void preprocess(void);

    // Estimation
    unsigned ls_degree;
    void est_ls(const unsigned n);
    void est_ls0(void) {est_ls(0); return;};
    void est_ls1(void) {est_ls(1); return;};
    void est_ls2(void) {est_ls(2); return;};
    void est_ls3(void) {est_ls(3); return;};
    void est_ls4(void) {est_ls(4); return;};
    void est_ls5(void) {est_ls(5); return;};

    // Writer
    void write(const char*);

    void write_header(FILE*);
    void write_variable(FILE*, const std::string);

    // Polynomial least squares
    void setup_ls_basis(const unsigned p);
    std::vector<DType> eval_poly_basis(const unsigned p, std::vector<DType>&);


};

#endif

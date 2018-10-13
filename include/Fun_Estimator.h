#ifndef FUN_ESTIMATOR_H
#define FUN_ESTIMATOR_H

#include <vector>
#include <utility>
#include <cstring>
#include <set>
#include <map>
#include "Data_Type.h"

class TDInput;

class FunEstimator {
  public:
    FunEstimator(char*, char*);

    void readData(const char*);

    void run(void) {
      (this->*preprocess)();
      (this->*estimate)(); };

    void outputPar(const char* ofile) {
      (this->*write)(ofile); };

    void setNumInput(const int N_in) {
      Ninput = N_in;
      return;
    }

    void setNumOutput(const int N_out) {
      Noutput = N_out;
      return;
    }

    void getReferenceInput(TDInput&, std::set<std::string>&);
    void getReferenceOutput(std::map< std::string, DType >&, std::set<std::string>&);

  private:
    int Ninput, Noutput;

    // The first one is input, the second one is output
    typedef std::pair< std::vector<DType>, std::vector<DType> > VarPair;
    std::vector<VarPair> data;

    // Reference values for the parameters
    // Used for normalizing the inputs
    std::vector<DType> refs;
    // Reference values for the outputs, could be zero (no reference)
    std::vector<DType> o_refs;

    typedef void (FunEstimator::*funPrec)(void);
    typedef void (FunEstimator::*funRun)(void);
    typedef void (FunEstimator::*funOutput)(const char*);

    funPrec   preprocess;
    funRun    estimate;
    funOutput write;

    // Preprocessor
    void prec_none(void) {}; // Do nothing
    void prec_log(void);

    // Estimation
    void est_ls(const unsigned n);
    void est_ls0(void) {est_ls(0); return;};
    void est_ls1(void) {est_ls(1); return;};
    void est_ls2(void) {est_ls(2); return;};
    void est_ls3(void) {est_ls(3); return;};
    void est_ls4(void) {est_ls(4); return;};
    void est_ls5(void) {est_ls(5); return;};

    // Writer
    void write_poly(const char*);
    void write_exp_poly(const char*);

};

#endif

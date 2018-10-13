#ifndef TDSOLVE_API_H
#define TDSOLVE_API_H

#include "Data_Type.h"
#include "TDSolve.h"
#include "TDInput.h"
#include <cstring>
#include <map>
#include <set>

class TDInput;
class TDSolve;

class TDFactory {
  public:
    TDFactory(TDInput, const char*);
    ~TDFactory(void) {
      if (tdSolve) delete tdSolve;
    }

    bool run(void);

    // Pair of input parameters
    std::map< std::string, DType > pars;
    // Pair of output parameters, indexed by oKeys;
    std::map< std::string, DType > sols;
    // Names of output variables
    std::set< std::string >        okeys;

  private:
    char     probName[40];
    TDInput  tdInput;
    TDSolve* tdSolve;

    bool run_pdgf(DType& Tterm, DType& Nperc);
};

#endif

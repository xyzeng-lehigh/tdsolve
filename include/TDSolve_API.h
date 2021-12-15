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
    TDFactory(TDInput*, const char*);
    ~TDFactory(void) {
      if (tdSolve) delete tdSolve;
    }

    void alloc(void);
    bool run(void);

    // Pair of input parameters
    std::map< std::string, DType > pars;
    // Pair of output parameters, indexed by oKeys;
    std::map< std::string, DType > sols;
    // Names of output variables
    std::set< std::string >        okeys;

    TDInput*& getTDInput(void) {
      return tdInput;
    }

  private:
    char     probName[40];
    bool     classical_tmeth;
    TDInput* tdInput;
    TDSolve* tdSolve;

    bool run_pdgf(DType& Tterm, DType& Gperc, DType& Nperc);
    bool run_opn(DType& Tterm, DType& Gperc, DType& Nperc);
    bool run_potassium(DType& B2peak, DType& B2vol, DType& B3vol, DType& Pvol, DType& Vvol);
};

#endif

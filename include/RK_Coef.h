#ifndef RK_COEF_H
#define RK_COEF_H

#include "Data_Type.h"

struct RKButcher {
  int    Nstage; // This will ends in a Nstage x Nstage table
  DType* coef;   // Butcher table
  DType& get(int, int);
  DType get(int, int) const;
  DType  stage(int);
  int getNstage(void) const {
    return Nstage;}

  RKButcher(int); 
};

#endif

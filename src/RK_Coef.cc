#include "RK_Coef.h"
#include <cassert>

RKButcher::RKButcher(int Nstage_in) {
  assert( Nstage_in > 0 );
  Nstage = Nstage_in;
  coef = new DType[Nstage*Nstage];
  memset( coef, 0.0, Nstage*Nstage*sizeof(DType) );
}

DType& RKButcher::get(int idx_i, int idx_j) {
  assert( idx_i >= 0 && idx_i < Nstage );
  assert( idx_j >= 0 && idx_j < Nstage );
  return coef[idx_i*Nstage+idx_j];
}

DType RKButcher::get(int idx_i, int idx_j) const {
  assert( idx_i >= 0 && idx_i < Nstage );
  assert( idx_j >= 0 && idx_j < Nstage );
  return coef[idx_i*Nstage+idx_j];
}

DType RKButcher::stage(int idx) {
  assert( idx >= 0 && idx < Nstage );
  DType row_sum = 0.0;
  int offset = idx*Nstage;
  for ( int i = 0; i < Nstage; i++ )
    row_sum += coef[offset+i];
  return row_sum;
}

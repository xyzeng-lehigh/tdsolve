#include "TDSolve.h"
#include "Alloc_Utils.h"
#include "TDInput.h"
#include "Alg_Solve.h"
using namespace std;

TDSolve::TDSolve(TDInput* tdInput_ptr) {
  tdInput = tdInput_ptr;
  tdCoef  = &(tdInput->tdCoef);
  tdMesh  = &(tdInput->tdMesh);
  tdSetup = &(tdInput->tdSetup);
  tdLin   = &(tdInput->tdLin);
  Nx      = tdMesh->Nx;
  dx      = tdMesh->dx;
  x_f     = tdMesh->x_f;
  x_c     = tdMesh->x_c;
  Na_m    = tdMesh->Na_m;
  Na      = tdMesh->Na;
 

  x_a_f   = tdMesh->x_a_f;
  x_a_c   = tdMesh->x_a_c;

  dtn     = tdSetup->dTinit;
  dtn_old = tdSetup->dTinit;
  tsol    = tdSetup->Tinit;

  reach_term_now  = tdSetup->Tterm <= tsol;
  reach_term_next = false;

  mat = NULL;
  prc = NULL;
  sol = NULL;
  rhs = NULL;
  res = NULL;
  linSolve = NULL;

//OPN Model
  mat_p1 = NULL;
  prc_p1 = NULL;
  sol_p1 = NULL;
  rhs_p1 = NULL;
  res_p1 = NULL;
  linSolve_p1 = NULL;


      mat_p2 = NULL;
      prc_p2 = NULL;
      sol_p2 = NULL;
      rhs_p2 = NULL;
      res_p2 = NULL;
      linSolve_p2 = NULL;
//
            mat_a = NULL;
            prc_a = NULL;
            sol_a = NULL;
            rhs_a = NULL;
	    res_a = NULL;
	    linSolve_a = NULL;
  iter = 0;

  use_tfreq = tdSetup->Tfreq > 0.0;
  reach_output_now  = false;
  reach_output_next = false;
  if ( use_tfreq ) 
    t_output = tdSetup->Tfreq;
}

TDSolve::~TDSolve(void) {
  erase_array( data_g, 3 );
  erase_array( data_h, 3 );
  erase_array( data_n, 3 );
  erase_array( data_a, 3 );
  //OPN Model
  erase_array(data_p1, 3);
  erase_array(data_p2, 3);
  //
  erase_array( data_v, 3 );

  if (mat!=NULL) delete mat;
  if (prc!=NULL) delete prc;
  if (rhs!=NULL) delete rhs;
  if (res!=NULL) delete res;
  if (sol!=NULL) delete sol;
  if (linSolve!=NULL) delete linSolve;

  if (mat_a!=NULL) delete mat_a;
  if (prc_a!=NULL) delete prc_a;
  if (rhs_a!=NULL) delete rhs_a;
  if (res_a!=NULL) delete res_a;
  if (sol_a!=NULL) delete sol_a;
  if (linSolve_a!=NULL) delete linSolve_a;

//OPN Model
if (mat_p1!=NULL) delete mat_p1;
if (prc_p1!=NULL) delete prc_p1;
if (rhs_p1!=NULL) delete rhs_p1;
if (res_p1!=NULL) delete res_p1;
if (sol_p1!=NULL) delete sol_p1;
if (linSolve_p1!=NULL) delete linSolve_p1;



  if (mat_p2!=NULL) delete mat_p2;
  if (prc_p2!=NULL) delete prc_p2;
  if (rhs_p2!=NULL) delete rhs_p2;
  if (res_p2!=NULL) delete res_p2;
  if (sol_p2!=NULL) delete sol_p2;
  if (linSolve_p2!=NULL) delete linSolve_p2;
}


void TDSolve::setupLinSolve( 
    LinMat*& mat_in, LinMat*& prc_in, LinVec*& sol_in, LinVec*& rhs_in, 
    LinVec*& res_in, LinSolve*& linSolve_in, int Nsize ) {
  mat_in = new LinSqSparseMat;
  mat_in->setBlckDim(1);
  mat_in->setEqualSize(Nsize);

  prc_in = new LinSqSparseMat;
  prc_in->setBlckDim(1);
  prc_in->setEqualSize(Nsize);

  sol_in = new LinVec;
  sol_in->setBlckDim(1);
  sol_in->setSize(Nsize);

  rhs_in = new LinVec;
  rhs_in->setBlckDim(1);
  rhs_in->setSize(Nsize);

  res_in = new LinVec;
  res_in->setBlckDim(1);
  res_in->setSize(Nsize);

  linSolve_in = new LinSolve;
  linSolve_in->setBlckDim(1);
  linSolve_in->setSize(Nsize);

  linSolve_in->linkMat(mat_in);
  linSolve_in->linkPrc(prc_in);
  linSolve_in->linkSol(sol_in);
  linSolve_in->linkRhs(rhs_in);
  linSolve_in->linkRes(res_in);

  linSolve_in->setPrintStatus(tdLin->Print);
  linSolve_in->setSolvePrec(tdLin->Prec,tdLin->PrecPar);
  switch (tdLin->Type) {
    case LSYS_TYPE_DIR:
      linSolve_in->setSolveDir();
      break;
    case LSYS_TYPE_ITER:
      linSolve_in->setMaxIter(tdLin->MaxIter);
      linSolve_in->setTargetRes(tdLin->TargetRes);
      linSolve_in->setSolveIter(tdLin->Iter);
      break;
    case LSYS_TYPE_KRYLOV:
      linSolve_in->setMaxIter(tdLin->MaxIter);
      linSolve_in->setTargetRes(tdLin->TargetRes);
      linSolve_in->setSolveKrylov(tdLin->Krylov,tdLin->KrylovCap);
      break;
  }

  return;
}

void TDSolve::average_sol(const DType coef_1, const int fr_1,
    const DType coef_2, const int fr_2, const int to) {
  *data_g[to] = coef_1 * (*data_g[fr_1]) + coef_2 * (*data_g[fr_2]);
  *data_h[to] = coef_1 * (*data_h[fr_1]) + coef_2 * (*data_h[fr_2]);
  *data_n[to] = coef_1 * (*data_n[fr_1]) + coef_2 * (*data_n[fr_2]);
  *data_a[to] = coef_1 * (*data_a[fr_1]) + coef_2 * (*data_a[fr_2]);
//OPN model
  *data_p1[to] = coef_1 * (*data_p1[fr_1]) + coef_2 * (*data_p1[fr_2]);
  *data_p2[to] = coef_1 * (*data_p2[fr_1]) + coef_2 * (*data_p2[fr_2]);

  data_r[to]  = coef_1 * data_r[fr_1]    + coef_2 * data_r[fr_2];

  this->update_vel(to);
  return;
}

void TDSolve::diagnostic(const int nt) {
  fprintf( stdout, "Iter %3d: Time %6.8f; Time step %6.8f; Radius %6.8f; Front-velocity %6.8f.\n",
      iter, tsol, dtn, data_r[nt], (*data_v[nt])[Nx] );
  return;
}

void TDSolve::update_rad(const int fr, const int to) {
  data_r[to] = data_r[fr] + dtn * (*data_v[fr])[Nx]; 
  return;
}


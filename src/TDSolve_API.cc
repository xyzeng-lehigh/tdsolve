#include <cmath>
#include "TDSolve_API.h"
#include "TDInput.h"
#include "TDSolve.h"
#include "TDSolveFDM.h"
#include "TDSolveFVM.h"
#include "TDSolveFVMMod.h"

TDFactory::TDFactory(TDInput* tdInput_in, const char* probName_in)  {
  tdInput = tdInput_in;
  tdSolve = NULL;
  pars.clear();
  sols.clear();
  okeys.clear();

  std::strcpy( probName, probName_in );

  switch ( tdInput->tdSetup.frame ) {
    case 0:
      tdSolve = new TDSolveFDM(tdInput);
      break;
    case 1:
      tdSolve = new TDSolveFVM(tdInput);
      break;
    case 2:
      tdSolve = new TDSolveFVMMod(tdInput);
      break;
    default:
      tdSolve = new TDSolveFVMMod(tdInput);
      break;
  }

  if ( std::strcmp( probName, "PDGF" ) == 0 ||
       std::strcmp( probName, "OPN" ) == 0 ) {
    okeys.insert("Tterm");
    okeys.insert("Nperc");
  }

  for ( std::set<std::string>::iterator it = okeys.begin();
        it != okeys.end(); it++ )
    sols[*it] = 0.0;

}

void TDFactory::alloc(void) {
  tdSolve->alloc();
  return;
}

bool TDFactory::run(void) {
  bool success = false;

  if ( std::strcmp( probName, "PDGF" ) == 0 ) {
    DType Tterm, Nperc;
    success = run_pdgf(Tterm, Nperc);
    sols["Tterm"] = Tterm;
    sols["Nperc"] = Nperc;
  }

  if ( std::strcmp( probName, "OPN" ) == 0 ) {
    DType Tterm, Nperc;
    success = run_opn(Tterm, Nperc);
    sols["Tterm"] = Tterm;
    sols["Nperc"] = Nperc;
  }

  return success;
}

bool TDFactory::run_pdgf(DType& Tterm, DType& Nperc) {
  tdSolve->setup();
  tdSolve->set_ic();
  tdSolve->calc_dt(0);
  tdSolve->diagnostic(0);

  while ( tdSolve->nextIter() ) {
    switch ( tdInput->tdSetup.tmeth ) {
      case 0: // First-order
        {
          bool vcoup_conv;
          while (1) {
            vcoup_conv = tdSolve->update_sol(0,1); // BE or FE
            if ( !vcoup_conv ) {
              tdSolve->reduce_dtn(0.5);
              tdSolve->diagnostic(0);
              continue;
            }
            else
              break;
          }
        }
        break;
      case 1: // Second-order
        if ( tdInput->tdSetup.vCoup == 0 ) { // TVD RK2
          tdSolve->update_sol(0,1);
          tdSolve->update_sol(1,2,1.0,2.0);
          tdSolve->average_sol(0.5,0,0.5,2,1);
        }
        else { // DIRK2
          bool vcoup_conv;
          while (1) {
            DType alpha = 1.0-0.5*std::sqrt(2.0);
            tdSolve->scale_dtn( alpha );
            vcoup_conv = tdSolve->update_sol(0,1);
            if ( !vcoup_conv ) {
              tdSolve->scale_dtn( 1.0/alpha );
              tdSolve->reduce_dtn( 0.5 );
              tdSolve->diagnostic(0);
              continue;
            }
            tdSolve->average_sol(2.0-1.0/alpha,0,1.0/alpha-1.0,1,2);
            vcoup_conv = tdSolve->update_sol(2,1,1.0/alpha-1.0,1.0/alpha);
            if ( !vcoup_conv ) {
              tdSolve->scale_dtn( 1.0/alpha );
              tdSolve->reduce_dtn( 0.5 );
              tdSolve->diagnostic(0);
              continue;
            }
            else {
              tdSolve->scale_dtn( 1.0/alpha );
              break;
            }
          }
        }
        break;
      default: // First-order
        tdSolve->update_sol(0,1);
        break;
    }
    tdSolve->output(0);
    tdSolve->copy_data(1,0);
    tdSolve->calc_dt(0);
    tdSolve->diagnostic(0);
  }

  tdSolve->update_vel(0,1,0.0,true);
  tdSolve->output(0);

  bool run_success = tdSolve->check_Rterm(0);
  if ( run_success ) {
    Tterm = tdSolve->get_tsol();
    Nperc = tdSolve->get_nperc(0);
  }

  return run_success;
}

bool TDFactory::run_opn(DType& Tterm, DType& Nperc) {
  tdSolve->setup();
  tdSolve->set_ic();
  tdSolve->calc_dt(0);
  tdSolve->diagnostic(0);

  while ( tdSolve->nextIter() ) {
    switch ( tdInput->tdSetup.tmeth ) {
      case 0: // First-order
        {
          bool vcoup_conv;
          while (1) {
            vcoup_conv = tdSolve->update_sol(0,1);
            if ( !vcoup_conv ) {
              tdSolve->reduce_dtn(0.5);
              tdSolve->diagnostic(0);
              continue;
            }
            else
              break;
          }
        }
        break;
      default: // First-order
        tdSolve->update_sol(0,1);
        break;
    }
    tdSolve->output(0);
    tdSolve->copy_data(1,0);
    tdSolve->calc_dt(0);
    tdSolve->diagnostic(0);
  }

  tdSolve->update_vel(0,1,0.0,true);
  tdSolve->output(0);

  bool run_success = tdSolve->check_Rterm(0);
  if ( run_success ) {
    Tterm = tdSolve->get_tsol();
    Nperc = tdSolve->get_nperc(0);
  }

  return run_success;
}


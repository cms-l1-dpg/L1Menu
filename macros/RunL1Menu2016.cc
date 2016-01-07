// ===========================================================================
// 
//       Filename:  RunL1Menu2015.cc
// 
//    Description:  An executable for running the L1Menu2015
// 
//        Version:  1.0
//        Created:  12/22/2015 16:36:48
//       Compiler:  g++ -std=c++11
// 
//         Author:  Zhenbin Wu (benwu)
//          Email:  zhenbin.wu@gmail.com
//        Company:  UIC, CMS@LPC, CDF@FNAL
// 
// ===========================================================================


#include <cstdlib>
#include <iostream>  
#include <string>
#include <vector>

#include "BasicRatePlots.C"
#include "L1Menu2016_minbias_cross_section.C"

#include "TSystem.h"
#include <FWCore/FWLite/interface/FWLiteEnabler.h>
// ===  FUNCTION  ============================================================
//         Name:  main
//  Description:  
// ===========================================================================
int main ( int argc, char *argv[] )
{
  gSystem->Load("libFWCoreFWLite");
  FWLiteEnabler::enable();

  //goRatePlots("RUN256843_Stage2",0,0);
  RunL1(true, true, 5);
  return EXIT_SUCCESS;
}				// ----------  end of function main  ----------

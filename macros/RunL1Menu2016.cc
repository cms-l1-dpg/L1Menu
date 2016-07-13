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
#include "getDistributions.C"
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

  //RunL1(true, true, 5);
  goRatePlots("run80xMC_zerobias_l1tint_v63p1_muonFix_PU28",0,12000000);
  //goRatePlots("run274241_zerobias_l1tint_v61p1_test_v2",0,100000);
  //goRatePlots("run20PUMC_stage2_l1tint_v58p1",0,10000000);
  //goGetDistributions(1000000);
  return EXIT_SUCCESS;
}				// ----------  end of function main  ----------

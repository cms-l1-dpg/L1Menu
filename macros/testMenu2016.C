// ===========================================================================
// 
//       Filename:  testMenu2016.C
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  01/14/2016 10:04:16
//       Compiler:  g++ -std=c++11
// 
//         Author:  Zhenbin Wu (benwu)
//          Email:  zhenbin.wu@gmail.com
//        Company:  UIC, CMS@LPC, CDF@FNAL
// 
// ===========================================================================
//
#include "L1Menu2016.h"

#include <cstdlib>
#include <iostream>  
#include <string>
#include <vector>

// ===  FUNCTION  ============================================================
//         Name:  main
//  Description:  
// ===========================================================================
int main ( int argc, char *argv[] )
{

  //// Declare the supported options.
  //po::options_description desc("Allowed options");
  //desc.add_options()
    //("help", "produce help message")
    //("compression", po::value<int>(), "set compression level")
    //;

  //po::variables_map vm;
  //po::store(po::parse_command_line(ac, av, desc), vm);
  //po::notify(vm);    

  //if (vm.count("help")) {
    //cout << desc << "\n";
    //return 1;
  //}

  //if (vm.count("compression")) {
    //cout << "Compression level was set to " 
      //<< vm["compression"].as<int>() << ".\n";
  //} else {
    //cout << "Compression level was not set.\n";
  //}


  L1Menu2016 men;
  men.ReadMenu("Menu_256843_Tune.txt");
  men.OpenWithList("./ntuples_256843_stage2.list");
  men.PrintConfig();
  men.Loop();
  return EXIT_SUCCESS;
}				// ----------  end of function main  ----------

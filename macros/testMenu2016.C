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

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// ===  FUNCTION  ============================================================
//         Name:  main
//  Description:  
// ===========================================================================
int main ( int argc, char *argv[] )
{
  namespace po = boost::program_options;
  // Declare the supported options.
  boost::program_options::options_description desc("Allowed options");
  const std::string defaultMenu = "menu/Menu_256843_Stage2.txt";
  const std::string defaultntuple = "ntuple/Run256843_stage2_Len.list";
  desc.add_options()
    ("help,h", "produce help message")
    ("menufile,m", po::value<std::string>()->default_value(defaultMenu), "set the input menu")
    ("filelist,l", po::value<std::string>()->default_value(defaultntuple), "set the input ntuple list")
    ("writetext,t", po::value<bool>()->default_value(true), "write rate to output")
    ("writecsv,c", po::value<bool>()->default_value(true), "write rate to output in CSV format")
    ("writeplot,p", po::value<bool>()->default_value(true), "write plot to output")
    ("outfilename,o", po::value<std::string>()->default_value("Auto"), "set output file name")
    ("outputdir,d", po::value<std::string>()->default_value("results"), "set output directory")
    ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Create output directory ~~~~~
  if (!boost::filesystem::is_directory(vm["outputdir"].as<std::string>() ))
  {
    boost::filesystem::create_directory(vm["outputdir"].as<std::string>());
  }

  L1Menu2016 men(vm["menufile"].as<std::string>(), vm["filelist"].as<std::string>());

  men.ConfigOutput(vm["writetext"].as<bool>(), 
      vm["writecsv"].as<bool>(), 
      vm["writeplot"].as<bool>(), 
      vm["outputdir"].as<std::string>(),
      vm["outfilename"].as<std::string>());

  men.PreLoop();
  men.Loop();
  men.PostLoop();

  //men.PrintConfig();
  //men.Loop();
  return EXIT_SUCCESS;
}				// ----------  end of function main  ----------

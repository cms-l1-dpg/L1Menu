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
  const std::string defaultMenu = "menu/Menu_259721_Stage2_Tune.txt";
  const std::string defaultntuple = "ntuple/Run259721_stage2_Len_New.list";
  //const std::string defaultntuple = "ntuple/Run259721_stage2_Len.list";
  desc.add_options()
    ("help,h", "produce help message")
    ("menufile,m",    po::value<std::string>()->default_value(defaultMenu),   "set the input menu")
    ("filelist,l",    po::value<std::string>()->default_value(defaultntuple), "set the input ntuple list")
    ("outfilename,o", po::value<std::string>()->default_value("Auto"),        "set output file name")
    ("outputdir,d",   po::value<std::string>()->default_value("results"),     "set output directory")
    ("writetext,t",   po::value<bool>()->default_value(true),                 "write rate to output")
    ("writecsv,c",    po::value<bool>()->default_value(true),                 "write rate to output in CSV format")
    ("writeplot,p",   po::value<bool>()->default_value(true),                 "write plot to output")
    ("doPlotRate",    po::bool_switch(),                                      "save rate plot to output")
    ("doPlotEff",     po::bool_switch(),                                      "save efficiency plot to output")
    ("doPrintLS",     po::bool_switch(),                                      "print out rate per LS to file")
    ("doPrintPU",     po::bool_switch(),                                      "print out rate per PU to file")
    ("maxEvent,n",    po::value<int>()->default_value(-1),                    "run number of events; -1 for all")
    ("nBunches,b",    po::value<int>(),                                       "set number of bunches")
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


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Get the running time config ~~~~~
  std::map<std::string, float> L1Config;
  for (const auto& it : vm) {
    if (!vm.count(it.first)) continue;
    auto& value = it.second.value();
    if (auto v = boost::any_cast<bool>(&value))
      L1Config[it.first] = *v;
    else if (auto v = boost::any_cast<int>(&value))
      L1Config[it.first] = *v;
    else if (auto v = boost::any_cast<float>(&value))
      L1Config[it.first] = *v;
  }
  men.PreLoop(L1Config);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Start to Run ~~~~~
  men.Loop();
  men.PostLoop();

  return EXIT_SUCCESS;
}				// ----------  end of function main  ----------

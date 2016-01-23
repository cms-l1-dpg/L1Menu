// ===========================================================================
// 
//       Filename:  L1Menu2016.h
// 
//    Description:
// 
//        Version:  1.0
//        Created:  01/13/2016 18:39:36
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Zhenbin Wu (benwu), zhenbin.wu@gmail.com
//        Company:  UIC, CMS@LPC, CDF@FNAL
// 
// ===========================================================================

#ifndef  MY_L1MENU2016_INC
#define  MY_L1MENU2016_INC


#include <map>
#include <regex>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>

// ROOT
#include "TH1F.h"
#include "TH2F.h"

// boost
#include "boost/tokenizer.hpp"
#include "boost/bind.hpp"

// Local
#include "L1AlgoFactory.h"
struct L1Seed
{
  std::string name;
  int bit;
  unsigned int prescale;
  std::vector<std::string> POG;
  std::vector<std::string> PAG;
  unsigned int ncounts;
  L1Seed()
  {
    name = "";
    bit = -1;
    prescale = 0;
    POG.clear();
    PAG.clear();
    ncounts = 0;
  }
};

struct StructL1Event
{
  float JetPt;
  float JetCenPt;
  float TauPt;
  float EGPt;
  float EGerPt;
  float MuPt;
  float MuerPt;
  float IsoEGPt;
  float IsoEGerPt;
  float HTT;
  float ETM;
  float ETT;

  StructL1Event()
  {
    JetPt     = -10;
    JetCenPt  = -10;
    TauPt     = -10;
    EGPt      = -10;
    EGerPt    = -10;
    IsoEGPt   = -10;
    IsoEGerPt = -10;
    MuPt      = -10;
    MuerPt    = -10;
  }
  
};

typedef bool(*FUNCPTR)();
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

// ===========================================================================
//        Class:  L1Menu2016
//  Description:  A class to handle the L1Menu
// ===========================================================================
class L1Menu2016 : public L1AlgoFactory 
{
  public:

    // ====================  LIFECYCLE     ===============================
    L1Menu2016 ();                             // constructor
    L1Menu2016 ( const L1Menu2016 &other );   // copy constructor
    ~L1Menu2016 ();                            // destructor

    // ====================  ACCESSORS     ===============================
    bool BindAlgo();
    bool InitConfig();
    bool OpenWithList(std::string filelist);
    bool ParseConfig(const std::string line);
    bool PrintConfig() const;
    bool PostLoop();
    bool Loop();
    bool CheckL1Seed(const std::string L1Seed);
    bool L1SeedFunc() const;
    bool PreLoop();
    bool BookHistogram();
    bool ReadMenu(std::string MenuName);
    bool BuildRelation();

    // ====================  MUTATORS      ===============================

    // ====================  OPERATORS     ===============================

    L1Menu2016& operator = ( const L1Menu2016 &other ); // assignment operator

    // ====================  DATA MEMBERS  ===============================
    bool writefiles;
    bool writecsv;
    bool drawplots;

  protected:
    // ====================  METHODS       ===============================

    bool GetL1Event();
    void CorrectScale(TH1F* h, Float_t scal);
    bool InsertInMenu(std::string L1name, bool value);
    Bool_t EGamma();
    bool ParseSingleObject(const std::string SeedName);
    // ====================  DATA MEMBERS  ===============================

  private:
    // ====================  METHODS       ===============================

    // ====================  DATA MEMBERS  ===============================
    StructL1Event L1Event;
    std::map<std::string, float*> L1ObjectMap;
    bool RunMenu();
    std::map<std::string, double> L1Config;
    std::map<std::string, L1Seed> mL1Seed;
    std::map<std::string, TH1F*> HistMap;

	std::map<std::string, std::function<bool()>> L1SeedFun;
    std::map<std::string, std::map<int, int> > L1LSCount; // counting lumi section



}; // -----  end of class L1Menu2016  -----

#endif   // ----- #ifndef MY_L1MENU2016_INC  -----

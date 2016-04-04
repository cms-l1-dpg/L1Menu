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
#include <iomanip>
#include <fstream>
#include <sstream>
#include <functional>

// ROOT
#include "TH1F.h"
#include "TH2F.h"

// boost
#include "boost/tokenizer.hpp"
#include "boost/filesystem.hpp"

// Local
#include "L1AlgoFactory.h"
#include "L1Plot.h"
#include "L1Struct.h"

#define INFTY 262139
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;


// ===  FUNCTION  ============================================================
//         Name:  SingleObjPt
//  Description:  /* cursor */
// ===========================================================================
inline bool SingleObjPt(float* obj, double pt)
{
  return *obj >=pt ;
}       // -----  end of function SingleObjPt  -----
// ===========================================================================
//        Class:  L1Menu2016
//  Description:  A class to handle the L1Menu
// ===========================================================================
class L1Menu2016 : public L1AlgoFactory 
{
  public:

    // ====================  LIFECYCLE     ===============================
    L1Menu2016 (std::string MenuName, std::string filelist);                             // constructor
    ~L1Menu2016 ();                            // destructor

    // ====================  ACCESSORS     ===============================
    bool BindAlgo();
    bool InitConfig();
    bool PrintRates(std::ostream &out);
    bool OpenWithList(std::string filelist);
    bool ParseConfig(const std::string line);
    bool PrintConfig() const;
    bool PostLoop();
    bool Loop();
    bool CheckL1Seed(const std::string L1Seed);
    bool L1SeedFunc();
    bool PreLoop(std::map<std::string, float> &config);
    bool BookHistogram();
    bool ReadMenu();
    bool ReadDataPU();
    bool BuildRelation();
    bool PrintPUCSV();
    bool WriteHistogram();
    bool GetRunConfig(std::map<std::string, float> &config);
    bool InitOutput();


    bool ConfigOutput(bool writetext_, bool writecsv_, bool writeplot_, 
        std::string outputdir_, std::string outputname_);
    std::string SetOutputName() const;
    // ====================  MUTATORS      ===============================

    bool ParseL1Seed(const std::string SeedName);
    bool ParseSingleObject(const std::string SeedName);
    bool ParseTripleJetVBF(const std::string& SeedName);
    bool ParseDoubleTau(const std::string& SeedName);
    bool ParseDoubleJet(const std::string& SeedName);
    bool ParseQuadJet(const std::string& SeedName);
    bool ParseDoubleEG(const std::string& SeedName);
    bool ParseTripleEG(const std::string& SeedName);
    bool ParseEGSum(const std::string& SeedName);

    bool ParseCrossMu(const std::string& SeedName);
    std::function<bool()> ParseBptx(const std::string Seedtoken);
    bool ParseMultiEGMass(const std::string& SeedName);
    bool ParseMuEG(const std::string& SeedName);
    bool ParseMuerTauer(const std::string& SeedName);
    bool ParseMuSum(const std::string& SeedName);
    // ====================  OPERATORS     ===============================

    L1Menu2016& operator = ( const L1Menu2016 &other ); // assignment operator

    // ====================  DATA MEMBERS  ===============================
    bool writefiles;
    bool writecsv;
    bool writeplots;
    std::string  outputdir;
    std::string  outputname;
    std::fstream *outfile;
    std::fstream *outcsv;
    TFile        *outrootfile;


  protected:
    // ====================  METHODS       ===============================

    bool GetL1Event();
    void CorrectScale(TH1F* h, Float_t scal);
    bool InsertInMenu(std::string L1name, bool value);
    bool FillLumiSection(int currentLumi);
    bool FillPileUpSec();
    bool PrintCSV(std::ostream &out);
    bool CheckPureFire();
    bool CheckPhysFire();
    // ====================  DATA MEMBERS  ===============================

  private:
    // ====================  METHODS       ===============================
    double CalScale(int nEvents_ = 0, int nBunches_ = 0, bool print=false);
    bool RunMenu();
    bool FillDefHist1D();
    void CalLocalHT(float &HTTcut);
    void CalLocalHTM(float &HTMcut);

    // ====================  DATA MEMBERS  ===============================
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Configurations ~~~~~
    std::string outfilename;
    std::string outfiledir;
    std::string menufilename;
    std::string tuplefilename;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Rate variables ~~~~~
    double scale;
    double nFireevents;
    unsigned int nZeroBiasevents;
    unsigned int nLumi;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ L1Seeds ~~~~~
    StructL1Event L1Event;
    L1Plot *l1Plot;

    std::map<std::string, float*> L1ObjectMap;
    std::map<std::string, float> L1Config;
    std::map<std::string, L1Seed> mL1Seed;
    std::map<std::string, TH1F*> HistMap;
    std::map<std::string, TH2F*> Hist2D;
	std::map<std::string, std::function<bool()>> L1SeedFun;
    std::map<int, std::string> BitMap;

    //Relationship
    std::set<std::string> FireSeed;
    std::map<std::string, std::map<int, int> > L1LSCount; // counting lumi section

    // Seed, PU, count
    std::map<std::string, std::map<double, int> > L1PUCount; // counting lumi section
    std::map<unsigned, std::map<unsigned, double> > DataLSPU; // mapping of PU for data

    std::map<std::string, std::vector<int> > POGMap;
    std::map<std::string, std::vector<int> > PAGMap;
    std::map<std::string, int > PhyCounts;
    std::map<std::string, int > PhyPureCounts;
    std::set<std::string> FiredPhy;


}; // -----  end of class L1Menu2016  -----

#endif   // ----- #ifndef MY_L1MENU2016_INC  -----

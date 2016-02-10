// ===========================================================================
// 
//       Filename:  L1Plot.h
// 
//    Description: G
// 
//        Version:  1.0
//        Created:  01/31/2016 03:55:43 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Zhenbin Wu (benwu), zhenbin.wu@gmail.com
//        Company:  UIC, CMS@LPC, CDF@FNAL
// 
// ===========================================================================

#ifndef  MY_L1PLOT_INC
#define  MY_L1PLOT_INC

#include <map>
#include <iostream>


#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

#include "L1Trigger/L1TNtuples/interface/L1AnalysisEventDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoJetDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMetDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoElectronDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMuon2DataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoTauDataFormat.h"

#include "L1Struct.h"

// ===========================================================================
//        Class:  L1Plot
//  Description:  
// ===========================================================================
class L1Plot
{
  public:

    // ====================  LIFECYCLE     ===============================
    L1Plot (
        TFile* outrootfile_,
        L1Analysis::L1AnalysisEventDataFormat        *event__ = nullptr,
        L1Analysis::L1AnalysisL1UpgradeDataFormat    *upgrade__ = nullptr,
        L1Analysis::L1AnalysisRecoJetDataFormat      *recoJet__ = nullptr,
        L1Analysis::L1AnalysisRecoMetDataFormat      *recoSum__ = nullptr,
        L1Analysis::L1AnalysisRecoElectronDataFormat *recoEle__ = nullptr,
        L1Analysis::L1AnalysisRecoMuon2DataFormat    *recoMuon__ = nullptr,
        L1Analysis::L1AnalysisRecoTauDataFormat      *recoTau__ = nullptr
        );

    L1Plot ( const L1Plot &other );   // copy constructor
    ~L1Plot ();                            // destructor

    // ====================  ACCESSORS     ===============================

    // ====================  MUTATORS      ===============================
    bool RunPlot( StructL1Event *L1Event_, std::map<std::string, L1Seed> *mL1Seed_);
    bool PostRun(double scale);
    bool PreRun();

    // ====================  OPERATORS     ===============================

    L1Plot& operator = ( const L1Plot &other ); // assignment operator

    // ====================  DATA MEMBERS  ===============================

  protected:
    // ====================  METHODS       ===============================
    bool BookRateHistogram();
    bool FillWriteRateHistogram();
    bool WriteRateHistogram(double scale) const;

    // ====================  DATA MEMBERS  ===============================

  private:
    // ====================  METHODS       ===============================

    // ====================  DATA MEMBERS  ===============================
    TFile        *outfile;
    L1Analysis::L1AnalysisEventDataFormat        *event_;
    L1Analysis::L1AnalysisL1UpgradeDataFormat    *upgrade_;
    L1Analysis::L1AnalysisRecoJetDataFormat      *recoJet_;
    L1Analysis::L1AnalysisRecoMetDataFormat      *recoSum_;
    L1Analysis::L1AnalysisRecoElectronDataFormat *recoEle_;
    L1Analysis::L1AnalysisRecoMuon2DataFormat     *recoMuon_;
    L1Analysis::L1AnalysisRecoTauDataFormat      *recoTau_;



    StructL1Event *L1Event;
    std::map<std::string, L1Seed> *mL1Seed;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Plots ~~~~~
    std::map<std::string,TH1F*> hRate1F;
    std::map<std::string,TH2F*> hRate2F;
}; // -----  end of class L1Plot  -----


#endif   // ----- #ifndef MY_L1PLOT_INC  -----

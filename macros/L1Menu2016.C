// ===========================================================================
// 
//       Filename:  L1Menu2016.C
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  01/13/2016 18:39:20
//       Compiler:  g++ -std=c++11
// 
//         Author:  Zhenbin Wu (benwu)
//          Email:  zhenbin.wu@gmail.com
//        Company:  UIC, CMS@LPC, CDF@FNAL
// 
// ===========================================================================

#include "L1Menu2016.h"

//----------------------------------------------------------------------------
//       Class:  L1Menu2016
//      Method:  L1Menu2016
// Description:  constructor
//----------------------------------------------------------------------------
L1Menu2016::L1Menu2016 (std::string MenuName, std::string filelist):
  writefiles(true),writecsv(false),writeplots(true),
  menufilename(MenuName), 
  tuplefilename(filelist),
  scale(0),
  l1Plot(nullptr),
  l1TnP(nullptr),
  l1uGT(nullptr),
  l1unpackuGT(nullptr)
{
  InitConfig();
}  // -----  end of method L1Menu2016::L1Menu2016  (constructor)  -----

//----------------------------------------------------------------------------
//       Class:  L1Menu2016
//      Method:  ~L1Menu2016
// Description:  destructor
//----------------------------------------------------------------------------
L1Menu2016::~L1Menu2016 ()
{
  outfile->close();
  outcsv->close();
  outrootfile->Close();
  fChain->Reset();
  //delete outfile;
  //delete outcsv;
  //delete outrootfile;
}  // -----  end of method L1Menu2016::-L1Menu2016  (destructor)  -----

//----------------------------------------------------------------------------
//       Class:  L1Menu2016
//      Method:  operator =
// Description:  assignment operator
//----------------------------------------------------------------------------
  L1Menu2016&
L1Menu2016::operator = ( const L1Menu2016 &other )
{
  if ( this != &other ) {
  }
  return *this;
}  // -----  end of method L1Menu2016::operator =  (assignment operator)  ---

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ConfigOutput
//  Description:  
// ===========================================================================
bool L1Menu2016::ConfigOutput(bool writetext_, bool writecsv_, bool writeplot_, 
    std::string outputdir_, std::string outputname_)
{
  writefiles = writetext_;
  writecsv = writecsv_;
  writeplots = writeplot_;
  outputdir = outputdir_;
  if (outputname_ == "Auto")
  {
    outputname = SetOutputName();
  }
  else
    outputname = outputname_;

  if (writefiles)
    outfile = new std::fstream( outputdir + "/" + outputname +".txt", std::fstream::out  );
  if (writecsv)
    outcsv = new std::fstream( outputdir + "/" + outputname +".csv", std::fstream::out  );
  if (writeplots)
  {
    std::string rootfilename = outputdir + "/" + outputname +".root";
    outrootfile =  new TFile( rootfilename.c_str(), "RECREATE");
  }
  return true;
}       // -----  end of function L1Menu2016::ConfigOutput  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::InitConfig
//  Description:  
// ===========================================================================
bool L1Menu2016::InitConfig()
{
  L1Config["SumJetET"]       = 0;
  L1Config["SumJetEta"]      = 999;
  L1Config["nBunches"]       = 2592; //default for 2017 nBunches
  L1Config["doPlotRate"]     = 0;
  L1Config["doPlotEff"]      = 0;
  L1Config["doPlotTest"]     = 0;
  L1Config["doPlotuGt"]      = 0;
  L1Config["doTnPMuon"]      = 0;
  L1Config["doPlotLS"]       = 0;
  L1Config["doPrintPU"]      = 0;
  L1Config["doCompuGT"]      = 0;
  L1Config["maxEvent"]       = -1;
  L1Config["SetMuonER"]      = -1;
  L1Config["SetNoPrescale"]  = 0;
  L1Config["IgnorePrescale"] = 0;
  L1Config["UseUpgradeLyr1"] = -1;
  L1Config["UseL1CaloTower"] = -1;
  L1Config["SelectRun"]      = -1;
  L1Config["SelectEvent"]    = -1;
  L1Config["UsePFMETNoMuon"] = 0;
  L1Config["UseuGTDecision"] = 0;
  L1Config["UseUnpackTree"]  = 0;
  L1Config["doScanLS"]       = 0;
  
  L1ConfigStr["SelectLS"] = "";
  L1ConfigStr["SelectBX"] = "";
  L1ConfigStr["Lumilist"] = "";
  L1ConfigStr["SelectCol"] = "";

  L1ObjectMap["Jet"]     = &L1Event.JetPt;
  L1ObjectMap["JetC"]    = &L1Event.JetCenPt;
  L1ObjectMap["Jeter3p0"]    = &L1Event.JetCenPt;
  L1ObjectMap["Tau"]     = &L1Event.TauPt;
  L1ObjectMap["Tauer"]   = &L1Event.TauerPt;
  L1ObjectMap["Tauer2p1"]   = &L1Event.TauerPt;
  L1ObjectMap["IsoTau"]  = &L1Event.IsoTauPt;
  L1ObjectMap["EG"]      = &L1Event.EGPt;
  L1ObjectMap["EGer"]    = &L1Event.EGerPt;
  L1ObjectMap["EGer2p1"]    = &L1Event.EGerPt;
  L1ObjectMap["IsoEG"]   = &L1Event.IsoEGPt;
  L1ObjectMap["IsoEGer"] = &L1Event.IsoEGerPt;
  L1ObjectMap["IsoEGer2p1"] = &L1Event.IsoEGerPt;
  L1ObjectMap["Mu"]      = &L1Event.MuPt;
  L1ObjectMap["MuOpen"]  = &L1Event.MuOpenPt;
  L1ObjectMap["Muer"]    = &L1Event.MuerPt;
  L1ObjectMap["Muer2p1"]    = &L1Event.MuerPt;
  L1ObjectMap["HTT"]     = &L1Event.HTT;
  L1ObjectMap["HTTer"]     = &L1Event.HTT;
  L1ObjectMap["HTM"]     = &L1Event.HTM;
  L1ObjectMap["ETM"]     = &L1Event.ETM;
  L1ObjectMap["ETT"]     = &L1Event.ETT;
  L1ObjectMap["ETMHF"]     = &L1Event.ETMHF;
  L1ObjectMap["HTTHF"]     = &L1Event.HTTHF;


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Map to old func for now. ~~~~~
  // MutliJets
  L1SeedFun["L1_QuadJetC36_TauJet52"] = std::bind(&L1AlgoFactory::QuadJetCentral_TauJet, this, 36.,52., false, false);
  L1SeedFun["L1_QuadJetC36_Tau52"] = std::bind(&L1AlgoFactory::QuadJetCentral_TauJet, this, 36.,52., false, false);
  L1SeedFun["L1_QuadJet36er3p0_Tau52"] = std::bind(&L1AlgoFactory::QuadJetCentral_TauJet, this, 36.,52., false, false);

  // MultiMuon
  L1SeedFun["L1_DoubleMu0er1p6_dEta_Max1p8_OS"] = std::bind(&L1AlgoFactory::Onia2015, this, 0.,0.,true,true,18, 1.6);
  L1SeedFun["L1_DoubleMu0er1p6_dEta_Max1p8"] = std::bind(&L1AlgoFactory::Onia2015, this, 0.,0.,true,false,18, 1.6);
  L1SeedFun["L1_DoubleMu0er1p25_dEta_Max1p8_OS"] = std::bind(&L1AlgoFactory::Onia2015, this, 0.,0.,true,true,18, 1.25);
  L1SeedFun["L1_DoubleMu0er1p0_dEta_Max1p8_OS"] = std::bind(&L1AlgoFactory::Onia2015, this, 0.,0.,true,true,18, 1.25);
  L1SeedFun["L1_DoubleMu0er1p4_dEta_Max1p8_OS"] = std::bind(&L1AlgoFactory::Onia2015, this, 0.,0.,true,true,18, 1.4);
  L1SeedFun["L1_DoubleMu_10_0_dEta_Max1p8"] = std::bind(&L1AlgoFactory::Onia2015, this, 10.,0.,false,false,18, 1.6);
  L1SeedFun["L1_DoubleMuOpen"] = std::bind(&L1AlgoFactory::DoubleMuOpen, this, 0.);
  L1SeedFun["L1_DoubleMu_10_Open"] = std::bind(&L1AlgoFactory::DoubleMuXOpen, this, 10.);
  L1SeedFun["L1_DoubleMu0"]       = std::bind(&L1AlgoFactory::DoubleMu , this , 0.  , 0.  , 1 , false);
  L1SeedFun["L1_DoubleMu_10_0"]   = std::bind(&L1AlgoFactory::DoubleMu , this , 10. , 0.  , 1 , false);
  L1SeedFun["L1_DoubleMu_10_3p5"] = std::bind(&L1AlgoFactory::DoubleMu , this , 10. , 3.5 , 1 , false);
  L1SeedFun["L1_DoubleMu_12_5"]   = std::bind(&L1AlgoFactory::DoubleMu , this , 12. , 5.  , 1 , false);
  L1SeedFun["L1_DoubleMu_11_4"]   = std::bind(&L1AlgoFactory::DoubleMu , this , 11. , 4.  , 1 , false);
  L1SeedFun["L1_DoubleMu_12_8"]   = std::bind(&L1AlgoFactory::DoubleMu , this , 12. , 8.  , 1 , false);
  L1SeedFun["L1_DoubleMu_13_6"]   = std::bind(&L1AlgoFactory::DoubleMu , this , 13. , 6.  , 1 , false);
  L1SeedFun["L1_DoubleMu_15_5"]   = std::bind(&L1AlgoFactory::DoubleMu , this , 15. , 5.  , 1 , false);
  L1SeedFun["L1_TripleMu0"] = std::bind(&L1AlgoFactory::TripleMu, this, 0.,0.,0.,1);
  L1SeedFun["L1_TripleMuOpen"] = std::bind(&L1AlgoFactory::TripleMu, this, 0.,0.,0.,0);
  L1SeedFun["L1_TripleMu_5_5_3"] = std::bind(&L1AlgoFactory::TripleMu, this, 5.,5.,3.,1);
  L1SeedFun["L1_TripleMu_5_0_0"] = std::bind(&L1AlgoFactory::TripleMu, this, 5.,0.,0.,1);
  L1SeedFun["L1_TripleMu_3_0_0"] = std::bind(&L1AlgoFactory::TripleMu, this, 3.,0.,0.,1);
  L1SeedFun["L1_QuadMu0"] = std::bind(&L1AlgoFactory::QuadMu, this, 0.,0.,0.,0.,1);

  //Cross
  L1SeedFun["L1_IsoEG20er_Tau20er_dEta_Min0p2"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 20.,20.,false);
  L1SeedFun["L1_IsoEG22er_Tau20er_dEta_Min0p2"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 22.,20.,false);
  L1SeedFun["L1_IsoEG20er_Tau24er_dEta_Min0p2"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 20.,24.,false);
  L1SeedFun["L1_IsoEG22er_IsoTau26er_dEta_Min0p2"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 22.,26.,true);
  L1SeedFun["L1_IsoEG20er_IsoTau25er_dEta_Min0p2"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 20.,25.,true);
  L1SeedFun["L1_IsoEG18er_IsoTau23er_dEta_Min0p2"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 18.,23.,true);
  L1SeedFun["L1_IsoEG18er_IsoTau25er_dEta_Min0p2"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 18.,25.,true);
  L1SeedFun["L1_IsoEG18er_IsoTau24er_dEta_Min0p2"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 18.,24.,true);
  L1SeedFun["L1_DoubleMu6_EG6"] = std::bind(&L1AlgoFactory::DoubleMu_EG, this, 6.,6.,1, false);
  L1SeedFun["L1_DoubleMu6_EG16"] = std::bind(&L1AlgoFactory::DoubleMu_EG, this, 6.,16.,1, false); // l1t-tsg-v3:  L1_DoubleMu6_EG6
  L1SeedFun["L1_DoubleMu7_EG7"] = std::bind(&L1AlgoFactory::DoubleMu_EG, this, 7,7.,1, false);
  L1SeedFun["L1_DoubleMu7_EG14"] = std::bind(&L1AlgoFactory::DoubleMu_EG, this, 7,14.,1, false);  // l1t-tsg-v3:  L1_DoubleMu7_EG7
  L1SeedFun["L1_Mu5_DoubleEG5"] = std::bind(&L1AlgoFactory::Mu_DoubleEG, this, 5., 5.);
  L1SeedFun["L1_Mu6_DoubleEG10"] = std::bind(&L1AlgoFactory::Mu_DoubleEG, this, 6., 10.);
  L1SeedFun["L1_Mu6_DoubleEG17"] = std::bind(&L1AlgoFactory::Mu_DoubleEG, this, 6., 17.); // l1t-tsg-v3:  L1_Mu6_DoubleEG10

  //MultiCross
  L1SeedFun["L1_Mu3_JetC16_dEta_Max0p4_dPhi_Max0p4"] = std::bind(&L1AlgoFactory::Mu_JetCentral_delta, this, 3.,16.);
  L1SeedFun["L1_Mu3_JetC52_dEta_Max0p4_dPhi_Max0p4"] = std::bind(&L1AlgoFactory::Mu_JetCentral_delta, this, 3.,52.);
  L1SeedFun["L1_Mu3_JetC60_dEta_Max0p4_dPhi_Max0p4"] = std::bind(&L1AlgoFactory::Mu_JetCentral_delta, this, 3.,60.);
  L1SeedFun["L1_Mu3_JetC120_dEta_Max0p4_dPhi_Max0p4"] = std::bind(&L1AlgoFactory::Mu_JetCentral_delta, this, 3.,120.);
  L1SeedFun["L1_DoubleJetC56_ETM60"] = std::bind(&L1AlgoFactory::DoubleJetCentral_ETM, this, 56.,56.,60.);
  L1SeedFun["L1_DoubleJetC60_ETM60"] = std::bind(&L1AlgoFactory::DoubleJetCentral_ETM, this, 60.,60.,60.);
  L1SeedFun["L1_DoubleJetC60_ETM70"] = std::bind(&L1AlgoFactory::DoubleJetCentral_ETM, this, 60.,60.,70.);
  L1SeedFun["L1_DoubleEG6_HTT150"] = std::bind(&L1AlgoFactory::DoubleEG_HT, this, 6., 150.);
  L1SeedFun["L1_DoubleEG6_HTT255"] = std::bind(&L1AlgoFactory::DoubleEG_HT, this, 6., 255.); // l1t-tsg-v3:  L1_DoubleEG6_HTT150
  L1SeedFun["L1_Jet32_DoubleMuOpen_Mu10_dPhi_Jet_Mu0_Max1p05_dPhi_Mu_Mu_Min1p0"] = std::bind(&L1AlgoFactory::Jet_MuOpen_Mu_dPhiMuMu1, this, 32.,10., 0);
  L1SeedFun["L1_Jet32_DoubleMuOpen_Mu10_dPhi_Jet_Mu0_Max0p4_dPhi_Mu_Mu_Min1p0"] = std::bind(&L1AlgoFactory::Jet_MuOpen_Mu_dPhiMuMu1, this, 32.,10., 0);
  L1SeedFun["L1_Jet32_DoubleMu_10_0_dPhi_Jet_Mu0_Max0p4_dPhi_Mu_Mu_Min1p0"] = std::bind(&L1AlgoFactory::Jet_MuOpen_Mu_dPhiMuMu1, this, 32.,10., 1 );
  L1SeedFun["L1_Jet32_MuOpen_EG10_dPhi_Jet_Mu_Max1p05_dPhi_Mu_EG_Min1p05"] = std::bind(&L1AlgoFactory::Jet_MuOpen_EG_dPhiMuEG1, this, 32.,10., 0);
  L1SeedFun["L1_Jet32_MuOpen_EG10_dPhi_Jet_Mu_Max0p4_dPhi_Mu_EG_Min1p0"] = std::bind(&L1AlgoFactory::Jet_MuOpen_EG_dPhiMuEG1, this, 32.,10.,0);
  L1SeedFun["L1_Jet32_Mu0_EG10_dPhi_Jet_Mu_Max0p4_dPhi_Mu_EG_Min1p0"] = std::bind(&L1AlgoFactory::Jet_MuOpen_EG_dPhiMuEG1, this, 32.,10.,1);
  L1SeedFun["L1_Jet32MuOpen_EG17_dPhiMu_EG1"] = std::bind(&L1AlgoFactory::Jet_MuOpen_EG_dPhiMuEG1, this, 32.,17., 0); // l1t-tsg-v3:  L1_Jet32MuOpen_EG10_dPhiMu_EG1
  L1SeedFun["L1_DoubleMu0_ETM40"] = std::bind(&L1AlgoFactory::DoubleMu_ETM, this, 0, 0, 40, false); 
  L1SeedFun["L1_DoubleMu0_ETM50"] = std::bind(&L1AlgoFactory::DoubleMu_ETM, this, 0, 0, 50, false); 
  L1SeedFun["L1_DoubleMu0_ETM55"] = std::bind(&L1AlgoFactory::DoubleMu_ETM, this, 0, 0, 55, false); 
  L1SeedFun["L1_DoubleMu0_ETM60"] = std::bind(&L1AlgoFactory::DoubleMu_ETM, this, 0, 0, 60, false); 
  L1SeedFun["L1_DoubleMu0_ETM65"] = std::bind(&L1AlgoFactory::DoubleMu_ETM, this, 0, 0, 65, false); 
  L1SeedFun["L1_DoubleMu0_ETM70"] = std::bind(&L1AlgoFactory::DoubleMu_ETM, this, 0, 0, 70, false); 
  L1SeedFun["L1_DoubleMu0_ETM75"] = std::bind(&L1AlgoFactory::DoubleMu_ETM, this, 0, 0, 75, false); 

  L1SeedFun["L1_DoubleJet8_ForwardBackward"] = std::bind(&L1AlgoFactory::DoubleJet_ForwardBackward, this, 8., 8.); 
  L1SeedFun["L1_DoubleJet12_ForwardBackward"] = std::bind(&L1AlgoFactory::DoubleJet_ForwardBackward, this, 12., 12.); 
  L1SeedFun["L1_DoubleJet16_ForwardBackward"] = std::bind(&L1AlgoFactory::DoubleJet_ForwardBackward, this, 16., 16.); 
  L1SeedFun["L1_ETM60_Jet60_dPhi_Min0p4"] = std::bind(&L1AlgoFactory::ETM_Jet, this, 60., 60., false); 
  L1SeedFun["L1_ETM70_Jet60_dPhi_Min0p4"] = std::bind(&L1AlgoFactory::ETM_Jet, this, 70., 60., false); 
  L1SeedFun["L1_ETM75_Jet60_dPhi_Min0p4"] = std::bind(&L1AlgoFactory::ETM_Jet, this, 75., 60., false); 
  L1SeedFun["L1_ETM85_Jet60_dPhi_Min0p4"] = std::bind(&L1AlgoFactory::ETM_Jet, this, 85., 60., false); 
  L1SeedFun["L1_HTM60_HTT260"] = std::bind(&L1AlgoFactory::HTM_HTT, this, 60., 260.); 
  L1SeedFun["L1_HTM80_HTT220"] = std::bind(&L1AlgoFactory::HTM_HTT, this, 80., 220.); 
  L1SeedFun["L1_Mu3_JetC35"] = std::bind(&L1AlgoFactory::Mu_Jet, this, 3., 35., false, true); 
  L1SeedFun["L1_Mu3_JetC16"] = std::bind(&L1AlgoFactory::Mu_Jet, this, 3., 16., false, true); 
  L1SeedFun["L1_Mu3_JetC60"] = std::bind(&L1AlgoFactory::Mu_Jet, this, 3., 60., false, true); 
  L1SeedFun["L1_Mu3_JetC120"] = std::bind(&L1AlgoFactory::Mu_Jet, this, 3., 120., false, true); 
  L1SeedFun["L1_MU20_EG15"] = std::bind(&L1AlgoFactory::Mu_EG, this,20, 15, false, 2);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Mass Trigger ~~~~~
  L1SeedFun["L1_DoubleJet_90_30_Mj30j30_580"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 90, 30., false, 30, 30, false, 580); 
  L1SeedFun["L1_DoubleJet_90_50_Mj30j30_580"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 90, 50., false, 30, 30, false, 580); 
  L1SeedFun["L1_DoubleJet_90_60_Mj30j30_580"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 90, 60., false, 30, 30, false, 580); 
  L1SeedFun["L1_DoubleJet_90_70_Mj30j30_580"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 90, 70., false, 30, 30, false, 580); 
  L1SeedFun["L1_DoubleJet_90_80_Mj30j30_580"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 90, 80., false, 30, 30, false, 580); 
  L1SeedFun["L1_DoubleJet_90_90_Mj30j30_580"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 90, 90., false, 30, 30, false, 580); 
  L1SeedFun["L1_DoubleJet_90_30_Mj30j30_580"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 90, 30., false, 30, 30, false, 580); 
  L1SeedFun["L1_DoubleJet_90_30_Mj30j30_610"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 90, 30., false, 30, 30, false, 610); 
  L1SeedFun["L1_DoubleJet_90_30_Mj30j30_640"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 90, 30., false, 30, 30, false, 640); 
  L1SeedFun["L1_DoubleJet_90_30_Mj30j30_670"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 90, 30., false, 30, 30, false, 670); 
  L1SeedFun["L1_DoubleJet_90_30_Mj30j30_700"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 90, 30., false, 30, 30, false, 700); 
  L1SeedFun["L1_DoubleJet_90_30_Mj30j30_730"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 90, 30., false, 30, 30, false, 730); 
  L1SeedFun["L1_DoubleJet_90_30_Mj30j30_760"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 90, 30., false, 30, 30, false, 760); 


  L1SeedFun["L1_DoubleJet30_Mj30j30_360_Mu6"] = std::bind(&L1AlgoFactory::DoubleJetMass_Mu, this, 30, 30., false, 30, 30, false, 360, 6, false, 2); 
  L1SeedFun["L1_DoubleJet30_Mj30j30_360_Mu10"] = std::bind(&L1AlgoFactory::DoubleJetMass_Mu, this, 30, 30., false, 30, 30, false, 360, 10, false, 2); 
  L1SeedFun["L1_DoubleJet_40_30_Mj40j30_540_IsoEG12"] = std::bind(&L1AlgoFactory::DoubleJetMass_Mu, this, 40, 30., false, 40, 30, false, 540, 12, true, false); 

//**************************************************************************//
//                           2017 L1Menu Proposals                          //
//**************************************************************************//

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Top:Halil ~~~~~
  L1SeedFun["L1_IsoEG24er2p1_Jet26er3p0_dR_Min0p3"] = std::bind(&L1AlgoFactory::EG_Jet, this, 24., 26., 0.3, true, true, true);
  L1SeedFun["L1_IsoEG26er2p1_Jet34er3p0_dR_Min0p3"] = std::bind(&L1AlgoFactory::EG_Jet, this, 26,  30,  0.3, true, true, true);
  L1SeedFun["L1_IsoEG28er2p1_Jet34er3p0_dR_Min0p3"] = std::bind(&L1AlgoFactory::EG_Jet, this, 28,  34,  0.3, true, true, true);
  L1SeedFun["L1_IsoEG30er2p1_Jet34er3p0_dR_Min0p3"] = std::bind(&L1AlgoFactory::EG_Jet, this, 30,  34,  0.3, true, true, true);
  L1SeedFun["L1_IsoEG24er2p1_HTT100er"] = std::bind(&L1AlgoFactory::EG_HTT, this, 24, 100, true, true);
  L1SeedFun["L1_IsoEG26er2p1_HTT100er"] = std::bind(&L1AlgoFactory::EG_HTT, this, 26, 100, true, true);
  L1SeedFun["L1_IsoEG28er2p1_HTT100er"] = std::bind(&L1AlgoFactory::EG_HTT, this, 28, 100, true, true);
  L1SeedFun["L1_IsoEG30er2p1_HTT100er"] = std::bind(&L1AlgoFactory::EG_HTT, this, 30, 100, true, true);
  //L1SeedFun["L1_IsoEG24er_TripleJetC26"] = std::bind(&L1AlgoFactory::EGer_TripleJetCentral, this, 24, 26, true, true);
  L1SeedFun["L1_Mu18_HTT100er"] = std::bind(&L1AlgoFactory::Mu_HTT, this, 18, 100);
  L1SeedFun["L1_Mu18_Jet24er3p0"] = std::bind(&L1AlgoFactory::Mu_Jet, this, 18, 24, false, true);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Top:Marina ~~~~~
  L1SeedFun["L1_HTT250er_QuadJet_70_55_40_35_er2p5"] = std::bind(&L1AlgoFactory::HTT_QuadJet, this, 250, 70, 55, 40, 35, 2.5);
  L1SeedFun["L1_HTT280er_QuadJet_70_55_40_35_er2p5"] = std::bind(&L1AlgoFactory::HTT_QuadJet, this, 280, 70, 55, 40, 35, 2.5);
  L1SeedFun["L1_HTT300er_QuadJet_70_55_40_35_er2p5"] = std::bind(&L1AlgoFactory::HTT_QuadJet, this, 300, 70, 55, 40, 35, 2.5);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Top:JIRA38 ~~~~~
  L1SeedFun["L1_DoubleJet100er2p3_dEta_Max1p6"] = std::bind(&L1AlgoFactory::DoubleJet_EtaRes_deltaEta, this, 100, 100, 2.3, 1.6);
  L1SeedFun["L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6"] = std::bind(&L1AlgoFactory::Mu_DoubleJet_Cor, this, 10, 32, 2.3, 0.4, 1.6);
  L1SeedFun["L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6"] = std::bind(&L1AlgoFactory::Mu_DoubleJet_Cor, this, 12, 40, 2.3, 0.4, 1.6);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EXO:JIRA48 ~~~~~
  L1SeedFun["L1_IsoTau40er_ETM80"] = std::bind(&L1AlgoFactory::Tau_ETM, this, 40, 80,  true, true, false);
  L1SeedFun["L1_IsoTau40er_ETM85"]  = std::bind(&L1AlgoFactory::Tau_ETM, this, 40, 85,  true, true, false);
  L1SeedFun["L1_IsoTau40er_ETM90"]  = std::bind(&L1AlgoFactory::Tau_ETM, this, 40, 90,  true, true, false);
  L1SeedFun["L1_IsoTau40er_ETM95"]  = std::bind(&L1AlgoFactory::Tau_ETM, this, 40, 95,  true, true, false);
  L1SeedFun["L1_IsoTau40er_ETM100"] = std::bind(&L1AlgoFactory::Tau_ETM, this, 40, 100, true, true, false);
  L1SeedFun["L1_IsoTau40er_ETMHF80"] = std::bind(&L1AlgoFactory::Tau_ETM, this, 40, 80,  true, true, true);
  L1SeedFun["L1_IsoTau40er_ETMHF90"] = std::bind(&L1AlgoFactory::Tau_ETM, this, 40, 90,  true, true, true);
  L1SeedFun["L1_QuadJet36er3p0_IsoTau52er2p1"] = std::bind(&L1AlgoFactory::QuadJetCentral_TauJet, this, 36, 52, true, true);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Higgs:Chiara ~~~~~
  L1SeedFun["L1_DoubleJet_90_30_DoubleJet30_Mass_Min620"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 90, 30., false, 30, 30, false, 620); 
  L1SeedFun["L1_DoubleJet_100_30_DoubleJet30_Mass_Min620"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 100, 30., false, 30, 30, false, 620); 
  L1SeedFun["L1_DoubleJet_100_35_DoubleJet35_Mass_Min620"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 100, 35., false, 35, 35, false, 620); 
  L1SeedFun["L1_DoubleJet_110_35_DoubleJet35_Mass_Min620"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 110, 35., false, 35, 35, false, 620); 
  L1SeedFun["L1_DoubleJet_110_40_DoubleJet40_Mass_Min620"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 110, 40., false, 35, 35, false, 620); 
  L1SeedFun["L1_DoubleJet_115_35_DoubleJet35_Mass_Min620"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 115, 35., false, 35, 35, false, 620); 
  L1SeedFun["L1_DoubleJet_115_40_DoubleJet40_Mass_Min620"] = std::bind(&L1AlgoFactory::DoubleJetMass, this, 115, 40., false, 40, 40, false, 620); 
  L1SeedFun["L1_DoubleJet30_Mass_Min400_Mu10"] = std::bind(&L1AlgoFactory::DoubleJetMass_Mu, this, 30, 30., false, 30, 30, false, 400, 10, false, 2); 
  L1SeedFun["L1_DoubleJet30_Mass_Min400_Mu6"] = std::bind(&L1AlgoFactory::DoubleJetMass_Mu, this, 30, 30., false, 30, 30, false, 400, 6, false, 2); 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BPH:Sergey ~~~~~
//
  L1SeedFun["L1_DoubleMu4_OS_EG12"] = std::bind(&L1AlgoFactory::DoubleMu_EG, this, 4.,12.,1, true);
  L1SeedFun["L1_DoubleMu5_OS_EG12"] = std::bind(&L1AlgoFactory::DoubleMu_EG, this, 5.,12.,1, true);
  L1SeedFun["L1_DoubleMu6_HighQ_OS"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 6.  , 6. , 999, 2 , true, 999, 999);
  //
  L1SeedFun["L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"] = std::bind(&L1AlgoFactory::DoubleMudRMax , this , 0   , 0   , 1.5, 2 , true, 1.4);
  L1SeedFun["L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4"] = std::bind(&L1AlgoFactory::DoubleMudRMax , this , 0   , 0   , 1.4, 2 , true, 1.4);
  L1SeedFun["L1_DoubleMu4_SQ_OS_dR_Max1p2"]      = std::bind(&L1AlgoFactory::DoubleMudRMax , this , 4   , 4   , 999, 2 , true, 1.2);
  L1SeedFun["L1_DoubleMu4p5_SQ_OS_dR_Max1p2"] = std::bind(&L1AlgoFactory::DoubleMudRMax , this , 4.5   , 4.5   , 999, 2 , true, 1.2);

  L1SeedFun["L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18"]  = std::bind(&L1AlgoFactory::DoubleMuMass  , this , 4.5 , 4.5 , 2.0, 2 , true, 7,    18);
  L1SeedFun["L1_DoubleMu5_SQ_OS_Mass7to18"]  = std::bind(&L1AlgoFactory::DoubleMuMass  , this , 5 , 5 , 999, 2 , true, 7,    18);
  L1SeedFun["L1_DoubleMu_20_2_SQ_Mass_Max20"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 20.  , 2. , 999, 2 , false, 0, 20);

  L1SeedFun["L1_TripleMu_5_0_0_DoubleMu_5_0_OS_Mass_Max17"] = std::bind(&L1AlgoFactory::TripleMu_DoubleMuMass, this, 5, 0, 0, 1, 5, 0, true, 0, 17);
  L1SeedFun["L1_TripleMu_5_3_0_DoubleMu_5_3_OS_Mass_Max17"] = std::bind(&L1AlgoFactory::TripleMu_DoubleMuMass, this, 5, 3, 0, 1, 5, 3, true, 0, 17);
  L1SeedFun["L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17"] = std::bind(&L1AlgoFactory::TripleMu_DoubleMuMass, this, 5, 3.5, 2.5, 1, 5, 2.5, true, 5, 17);
  L1SeedFun["L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17"] = std::bind(&L1AlgoFactory::TripleMu_DoubleMuMass, this, 5, 4, 2.5, 1, 5, 2.5, true, 5, 17);
  

  L1SeedFun["L1_DoubleMu0_SQ"]    = std::bind(&L1AlgoFactory::DoubleMuMass , this , 0. , 0. , 999 , 2, false, 999, 999);
  L1SeedFun["L1_DoubleMu0_SQ_OS"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 0. , 0. , 999 , 2, true,  999, 999);
  L1SeedFun["L1_DoubleMu0er1p5_SQ_OS"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 0. , 0. , 1.5061875 , 2, true,  999, 999);
  L1SeedFun["L1_DoubleMu4_SQ_OS"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 4. , 4. , 999 , 2, true,  999, 999);
  L1SeedFun["L1_DoubleMu4p5_SQ"]    = std::bind(&L1AlgoFactory::DoubleMuMass , this , 4.5 , 4.5 , 999 , 2, false, 999, 999);
  L1SeedFun["L1_DoubleMu4p5_SQ_OS"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 4.5 , 4.5 , 999 , 2, true,  999, 999);
  L1SeedFun["L1_DoubleMu4p5er2p0_SQ_OS"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 4.5 , 4.5 ,2.0 , 2, true,  999, 999);
  L1SeedFun["L1_DoubleMu5_SQ_OS"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 5 , 5 , 999 , 2, true,  999, 999);
  L1SeedFun["L1_DoubleMu6_SQ_OS"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 6 , 6 , 999 , 2, true,  999, 999);
  L1SeedFun["L1_DoubleMu8_SQ"]    = std::bind(&L1AlgoFactory::DoubleMuMass , this , 8. , 8. , 999 , 2, false, 999, 999);

  L1SeedFun["L1_TripleMu0_OQ"] = std::bind(&L1AlgoFactory::TripleMu, this, 0.,0.,0.,0);
  L1SeedFun["L1_TripleMu_5_3p5_2p5"] = std::bind(&L1AlgoFactory::TripleMu, this, 5.,3.5,2.5,1);
  L1SeedFun["L1_TripleMu_5_5_5_OQ_OS"] = std::bind(&L1AlgoFactory::TripleMuOS, this, 5., 5, 5, 0, 5, 5, true);
  L1SeedFun["L1_DoubleMu_12_8_SQ"]   = std::bind(&L1AlgoFactory::DoubleMu , this , 12. , 8.  , 2 , false);
  L1SeedFun["L1_DoubleMu_15_5_SQ"]   = std::bind(&L1AlgoFactory::DoubleMu , this , 15. , 5.  , 2 , false);

  //L1SeedFun["L1_DoubleMu_12_5_OS"]      = std::bind(&L1AlgoFactory::DoubleMuMass , this , 12. , 5. , 999, 1 , true, 999, 999);
  //L1SeedFun["L1_DoubleMu_9_5_HighQ_OS"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 9.  , 5. , 999, 2 , true, 999, 999);
  //L1SeedFun["L1_DoubleMu6_OS_M0to9"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 6.  , 6. , 999, 1 , true, 0, 9);
  //L1SeedFun["L1_DoubleMu5_HighQ_OS_M0to9"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 5.  , 5. , 999, 2 , true, 0, 9);
  //L1SeedFun["L1_DoubleMu6_HighQ_OS_M0to16"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 6.  , 6. , 999, 2 , true, 0, 16);

  //L1SeedFun["L1_DoubleMu0er1p4_HighQ_OS_M0to9"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 0.  , 0. , 1.4, 2 , true, 0, 9);
  //L1SeedFun["L1_DoubleMu0er1p5_HighQ_OS_M0to9"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 0.  , 0. , 1.5, 2 , true, 0, 9);
  //L1SeedFun["L1_DoubleMu3er1p6_HighQ_OS_M6to20"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 3.  , 3. , 1.6, 2 , true, 6, 20);
  //L1SeedFun["L1_DoubleMu3er1p8_HighQ_OS_M6to20"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 3.  , 3. , 1.8, 2 , true, 6, 20);
  //L1SeedFun["L1_DoubleMu4er1p8_HighQ_OS_M6to20"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 4.  , 4. , 1.8, 2 , true, 6, 20);

  //L1SeedFun["L1_DoubleMu4p5_HighQ_OS_M0to9"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 4.5  , 4.5 , 999, 2 , true, 0, 9);
  //L1SeedFun["L1_DoubleMu6_HighQ_OS_M0to20"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 6  , 6 , 999, 2 , true, 0, 20);
  //L1SeedFun["L1_DoubleMu_7_4_HighQ_OS_M0to9"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 7  , 4 , 999, 2 , true, 0, 9);
  //L1SeedFun["L1_DoubleMu_7_5_HighQ_OS_M0to20"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 7  , 5 , 999, 2 , true, 0, 20);
  //L1SeedFun["L1_DoubleMu_9_2_HighQ_OS_M0to9"] = std::bind(&L1AlgoFactory::DoubleMuMass , this , 9  , 2 , 999, 2 , true, 0, 9);

  //L1SeedFun["L1_DoubleMu0er1p4_HighQ_OS_dR_Max1p3"] = std::bind(&L1AlgoFactory::DoubleMudRMax , this , 0  , 0 , 1.4, 2 , true, 1.3);
  //L1SeedFun["L1_DoubleMu3er1p6_HighQ_OS_dR_Max2p0"] = std::bind(&L1AlgoFactory::DoubleMudRMax , this , 3  , 3 , 1.6, 2 , true, 2.0);
  //L1SeedFun["L1_DoubleMu3er1p8_HighQ_OS_dR_Max2p0"] = std::bind(&L1AlgoFactory::DoubleMudRMax , this , 3  , 3 , 1.8, 2 , true, 2.0);
  //L1SeedFun["L1_DoubleMu4er1p8_HighQ_OS_dR_Max2p0"] = std::bind(&L1AlgoFactory::DoubleMudRMax , this , 4  , 4 , 1.8, 2 , true, 2.0);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SUSY:Laurent ~~~~~
  L1SeedFun["L1_DoubleMu_15_7"]       = std::bind(&L1AlgoFactory::DoubleMu , this , 15. , 7. , 1 , false);
  L1SeedFun["L1_DoubleMu_15_7_SQ"] = std::bind(&L1AlgoFactory::DoubleMu , this , 15. , 7. , 2 , false);
  L1SeedFun["L1_DoubleMu_15_7_SQ_Mass_Min4"] = std::bind(&L1AlgoFactory::DoubleMuMass, this , 15. , 7. , 999, 2 , false, 4, 999);
  L1SeedFun["L1_TripleMu3"] = std::bind(&L1AlgoFactory::TripleMu, this, 3.,3.,3.,1);
  L1SeedFun["L1_TripleMu3_SQ"] = std::bind(&L1AlgoFactory::TripleMu, this, 3.,3.,3.,2);
  L1SeedFun["L1_TripleMu_5_3_3"] = std::bind(&L1AlgoFactory::TripleMu, this, 5.,3.,3.,1);
  L1SeedFun["L1_TripleMu_4_4_4"] = std::bind(&L1AlgoFactory::TripleMu, this, 4.,4.,4.,1);
  L1SeedFun["L1_DoubleMu7_SQ_EG7"] = std::bind(&L1AlgoFactory::DoubleMu_EG, this, 7,7.,2, false);

  L1SeedFun["L1_DoubleEG6_HTT240er"] = std::bind(&L1AlgoFactory::DoubleEG_HT, this, 6., 240.);
  L1SeedFun["L1_DoubleEG6_HTT250er"] = std::bind(&L1AlgoFactory::DoubleEG_HT, this, 6., 250.);
  L1SeedFun["L1_DoubleEG6_HTT255er"] = std::bind(&L1AlgoFactory::DoubleEG_HT, this, 6., 255.);
  L1SeedFun["L1_DoubleEG6_HTT270er"] = std::bind(&L1AlgoFactory::DoubleEG_HT, this, 6., 270.);
  L1SeedFun["L1_DoubleEG6_HTT300er"] = std::bind(&L1AlgoFactory::DoubleEG_HT, this, 6., 300.);

  L1SeedFun["L1_DoubleMu0_ETMHF40_Jet60_OR_DoubleJet30"] = std::bind(&L1AlgoFactory::DoubleMu_ETMHF_Jets, this, 0, 0, 1, false, 40, true, 60, 30);
  L1SeedFun["L1_DoubleMu0_ETMHF50_Jet60_OR_DoubleJet30"] = std::bind(&L1AlgoFactory::DoubleMu_ETMHF_Jets, this, 0, 0, 1, false, 50, true, 60, 30);
  L1SeedFun["L1_DoubleMu0_ETMHF60_Jet60_OR_DoubleJet30"] = std::bind(&L1AlgoFactory::DoubleMu_ETMHF_Jets, this, 0, 0, 1, false, 60, true, 60, 30);
  L1SeedFun["L1_DoubleMu0_ETMHF70_Jet60_OR_DoubleJet30"] = std::bind(&L1AlgoFactory::DoubleMu_ETMHF_Jets, this, 0, 0, 1, false, 70, true, 60, 30);
  L1SeedFun["L1_DoubleMu0_ETMHF80_Jet60_OR_DoubleJet30"] = std::bind(&L1AlgoFactory::DoubleMu_ETMHF_Jets, this, 0, 0, 1, false, 80, true, 60, 30);
  L1SeedFun["L1_DoubleMu3_SQ_ETMHF40_Jet60_OR_DoubleJet30"] = std::bind(&L1AlgoFactory::DoubleMu_ETMHF_Jets, this, 3, 3, 2, false, 40, true, 60, 30);
  L1SeedFun["L1_DoubleMu3_SQ_ETMHF50_Jet60_OR_DoubleJet30"] = std::bind(&L1AlgoFactory::DoubleMu_ETMHF_Jets, this, 3, 3, 2, false, 50, true, 60, 30);
  L1SeedFun["L1_DoubleMu3_SQ_ETMHF60_Jet60_OR_DoubleJet30"] = std::bind(&L1AlgoFactory::DoubleMu_ETMHF_Jets, this, 3, 3, 2, false, 60, true, 60, 30);
  L1SeedFun["L1_DoubleMu3_SQ_ETMHF70_Jet60_OR_DoubleJet30"] = std::bind(&L1AlgoFactory::DoubleMu_ETMHF_Jets, this, 3, 3, 2, false, 70, true, 60, 30);
  L1SeedFun["L1_DoubleMu3_SQ_ETMHF80_Jet60_OR_DoubleJet30"] = std::bind(&L1AlgoFactory::DoubleMu_ETMHF_Jets, this, 3, 3, 2, false, 80, true, 60, 30);


  L1SeedFun["L1_DoubleMu3_SQ_HTT100er"] = std::bind(&L1AlgoFactory::DoubleMu_HTT, this, 3, 3, 2, false, 100);
  L1SeedFun["L1_DoubleMu3_SQ_HTT200er"] = std::bind(&L1AlgoFactory::DoubleMu_HTT, this, 3, 3, 2, false, 200);
  L1SeedFun["L1_DoubleMu3_SQ_HTT220er"] = std::bind(&L1AlgoFactory::DoubleMu_HTT, this, 3, 3, 2, false, 220);
  L1SeedFun["L1_DoubleMu3_SQ_HTT240er"] = std::bind(&L1AlgoFactory::DoubleMu_HTT, this, 3, 3, 2, false, 240);
  //L1SeedFun["L1_DoubleMu3_SQ_HTT100er"] = std::bind(&L1AlgoFactory::DoubleMu_HTT, this, 0, 0, 2, false, 100);

  L1SeedFun["L1_ETMHF70_Jet60_OR_DoubleJet30"]  = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 70 , true, 60, 30, 0, false);
  L1SeedFun["L1_ETMHF75_Jet60_OR_DoubleJet30"]  = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 75 , true, 60, 30, 0, false);
  L1SeedFun["L1_ETMHF80_Jet60_OR_DoubleJet30"]  = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 80 , true, 60, 30, 0, false);
  L1SeedFun["L1_ETMHF85_Jet60_OR_DoubleJet30"]  = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 85 , true, 60, 30, 0, false);
  L1SeedFun["L1_ETMHF90_Jet60_OR_DoubleJet30"]  = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 90 , true, 60, 30, 0, false);
  L1SeedFun["L1_ETMHF95_Jet60_OR_DoubleJet30"]  = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 95 , true, 60, 30, 0, false);
  L1SeedFun["L1_ETMHF100_Jet60_OR_DoubleJet30"] = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 100, true, 60, 30, 0, false);
  L1SeedFun["L1_ETMHF105_Jet60_OR_DoubleJet30"] = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 105, true, 60, 30, 0, false);
  L1SeedFun["L1_ETMHF110_Jet60_OR_DoubleJet30"] = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 110, true, 60, 30, 0, false);
  L1SeedFun["L1_ETMHF115_Jet60_OR_DoubleJet30"] = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 115, true, 60, 30, 0, false);
  L1SeedFun["L1_ETMHF120_Jet60_OR_DoubleJet30"] = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 120, true, 60, 30, 0, false);

  L1SeedFun["L1_ETMHF70_Jet90_OR_DoubleJet45_OR_TripleJet30"] = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 70, true, 90, 45, 30, false);
  L1SeedFun["L1_ETMHF80_Jet90_OR_DoubleJet45_OR_TripleJet30"] = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 80, true, 90, 45, 30, false);
  L1SeedFun["L1_ETMHF90_Jet90_OR_DoubleJet45_OR_TripleJet30"] = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 90, true, 90, 45, 30, false);
  L1SeedFun["L1_ETMHF100_Jet90_OR_DoubleJet45_OR_TripleJet30"] = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 100, true, 90, 45, 30, false);
  L1SeedFun["L1_ETMHF110_Jet90_OR_DoubleJet45_OR_TripleJet30"] = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 110, true, 90, 45, 30, false);
  L1SeedFun["L1_ETMHF120_Jet90_OR_DoubleJet45_OR_TripleJet30"] = std::bind(&L1AlgoFactory::ETM_JetsComb, this, 120, true, 90, 45, 30, false);

  L1SeedFun["L1_ETMHF100_HTTer60"] = std::bind(&L1AlgoFactory::ETM_HTT, this, 100, 60, true);
  L1SeedFun["L1_DoubleJet30_Mass_Min300_dEta_Max1p5"] = std::bind(&L1AlgoFactory::DoubleJet_dEtaMass, this, 30, 30, 999, false, 1.5, 300);
  L1SeedFun["L1_DoubleJet30_Mass_Min320_dEta_Max1p5"] = std::bind(&L1AlgoFactory::DoubleJet_dEtaMass, this, 30, 30, 999, false, 1.5, 320);
  L1SeedFun["L1_DoubleJet30_Mass_Min340_dEta_Max1p5"] = std::bind(&L1AlgoFactory::DoubleJet_dEtaMass, this, 30, 30, 999, false, 1.5, 340);
  L1SeedFun["L1_DoubleJet30_Mass_Min360_dEta_Max1p5"] = std::bind(&L1AlgoFactory::DoubleJet_dEtaMass, this, 30, 30, 999, false, 1.5, 360);
  L1SeedFun["L1_DoubleJet30_Mass_Min380_dEta_Max1p5"] = std::bind(&L1AlgoFactory::DoubleJet_dEtaMass, this, 30, 30, 999, false, 1.5, 380);
  L1SeedFun["L1_DoubleJet30_Mass_Min400_dEta_Max1p5"] = std::bind(&L1AlgoFactory::DoubleJet_dEtaMass, this, 30, 30, 999, false, 1.5, 400);
  L1SeedFun["L1_DoubleLeadJet30_Mass_Min300_dEta_Max1p5"] = std::bind(&L1AlgoFactory::DoubleJet_dEtaMass, this, 30, 30, 999, true, 1.5, 300);
  L1SeedFun["L1_DoubleLeadJet30_Mass_Min320_dEta_Max1p5"] = std::bind(&L1AlgoFactory::DoubleJet_dEtaMass, this, 30, 30, 999, true, 1.5, 320);
  L1SeedFun["L1_DoubleLeadJet30_Mass_Min340_dEta_Max1p5"] = std::bind(&L1AlgoFactory::DoubleJet_dEtaMass, this, 30, 30, 999, true, 1.5, 340);
  L1SeedFun["L1_DoubleLeadJet30_Mass_Min360_dEta_Max1p5"] = std::bind(&L1AlgoFactory::DoubleJet_dEtaMass, this, 30, 30, 999, true, 1.5, 360);
  L1SeedFun["L1_DoubleLeadJet30_Mass_Min380_dEta_Max1p5"] = std::bind(&L1AlgoFactory::DoubleJet_dEtaMass, this, 30, 30, 999, true, 1.5, 380);
  L1SeedFun["L1_DoubleLeadJet30_Mass_Min400_dEta_Max1p5"] = std::bind(&L1AlgoFactory::DoubleJet_dEtaMass, this, 30, 30, 999, true, 1.5, 400);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HIG:Javier ~~~~~
  L1SeedFun["L1_ETM80_Jet60_dPhi_Min0p4"] = std::bind(&L1AlgoFactory::ETM_Jet, this, 80., 60., false); 
  L1SeedFun["L1_ETM90_Jet60_dPhi_Min0p4"] = std::bind(&L1AlgoFactory::ETM_Jet, this, 90., 60., false); 
  L1SeedFun["L1_ETM100_Jet60_dPhi_Min0p4"] = std::bind(&L1AlgoFactory::ETM_Jet, this, 100., 60., false); 
  L1SeedFun["L1_ETM110_Jet60_dPhi_Min0p4"] = std::bind(&L1AlgoFactory::ETM_Jet, this, 110., 60., false); 
  L1SeedFun["L1_DoubleJet60er3p0_ETM60"] = std::bind(&L1AlgoFactory::DoubleJetCentral_ETM, this, 60.,60.,60.);
  L1SeedFun["L1_DoubleJet60er3p0_ETM70"] = std::bind(&L1AlgoFactory::DoubleJetCentral_ETM, this, 60.,60.,70.);
  L1SeedFun["L1_DoubleJet60er3p0_ETM80"] = std::bind(&L1AlgoFactory::DoubleJetCentral_ETM, this, 60.,60.,80.);
  L1SeedFun["L1_DoubleJet60er3p0_ETM90"] = std::bind(&L1AlgoFactory::DoubleJetCentral_ETM, this, 60.,60.,90.);
  L1SeedFun["L1_DoubleJet60er3p0_ETM100"] = std::bind(&L1AlgoFactory::DoubleJetCentral_ETM, this, 60.,60.,100.);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HIG:Sumit ~~~~~
  L1SeedFun["L1_DoubleEG_Iso23_10"] = std::bind(&L1AlgoFactory::DoubleEGIsoPer, this, 23.,10.,true, false, false );
  L1SeedFun["L1_DoubleEG_Iso24_10"] = std::bind(&L1AlgoFactory::DoubleEGIsoPer, this, 24.,10.,true, false, false );
  L1SeedFun["L1_TripleEG_Iso20_10_5"] = std::bind(&L1AlgoFactory::TripleEGIsoPer, this, 20, 10, 5, true, false, false);
  L1SeedFun["L1_DoubleIsoEG22er2p1"] = std::bind(&L1AlgoFactory::DoubleEGIsoPer, this, 22.,22.,true, true, true );
  L1SeedFun["L1_DoubleMu18er2p1"]   = std::bind(&L1AlgoFactory::DoubleMu , this , 18. , 18.  , 1 , true);
  L1SeedFun["L1_DoubleMu22er2p1"]   = std::bind(&L1AlgoFactory::DoubleMu , this , 22. , 22.  , 1 , true);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ECal:Andrea ~~~~~
  L1SeedFun["L1_SingleIsoEG33_Mt40"] = std::bind(&L1AlgoFactory::SingleEGMT, this, 33.,40., 0, true, false);
  L1SeedFun["L1_SingleIsoEG33_Mt43"] = std::bind(&L1AlgoFactory::SingleEGMT, this, 33.,43., 0, true, false);
  L1SeedFun["L1_SingleIsoEG33_Mt47"] = std::bind(&L1AlgoFactory::SingleEGMT, this, 33.,47., 0, true, false);
  L1SeedFun["L1_SingleIsoEG33_Mt48"] = std::bind(&L1AlgoFactory::SingleEGMT, this, 33.,48., 0, true, false);
  L1SeedFun["L1_SingleIsoEG34_Mt44"] = std::bind(&L1AlgoFactory::SingleEGMT, this, 34.,44., 0, true, false);
  L1SeedFun["L1_SingleIsoEG34_Mt46"] = std::bind(&L1AlgoFactory::SingleEGMT, this, 34.,46., 0, true, false);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tracker:L1DPG74 ~~~~~
  L1SeedFun["L1_CDC_3_er1p2_TOP120_DPHI2p618_3p665"] = std::bind(&L1AlgoFactory::MuonCDC_dPhi, this, 3, 0, 1.2, 2.618, 3.665);
  L1SeedFun["L1_CDC_4_er1p2_TOP120_DPHI2p618_3p665"] = std::bind(&L1AlgoFactory::MuonCDC_dPhi, this, 4, 0, 1.2, 2.618, 3.665);

  L1SeedFun["L1_ZeroBias"] = [](){return true;};
  return true;

}       // -----  end of function L1Menu2016::InitConfig  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::BookHistogram
//  Description:  
// ===========================================================================
bool L1Menu2016::BookHistogram()
{
  for(auto col : ColumnMap)
  {
    col.second->BookHistogram();
  }
  return true;
}       // -----  end of function L1Menu2016::BookHistogram  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::WriteHistogram
//  Description:  
// ===========================================================================
bool L1Menu2016::WriteHistogram() 
{
  for(auto col : ColumnMap)
  {
    col.second->WriteHistogram(outrootfile);
  }
  return true;
}       // -----  end of function L1Menu2016::WriteHistogram  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ReadMenu
//  Description:  /* cursor */
// ===========================================================================
bool L1Menu2016::ReadMenu()
{
  mL1Seed.clear();
  vL1Seed.clear();

  //Read the prescales table
  std::ifstream menufile(menufilename);
  if (!menufile)
  {
    std::cout << "MenuFile "<<menufilename<<" is not found !"<<std::endl;
    return false;
  }

  if (writefiles)
    *outfile <<  "---------------------------- Input Menu -------------------------" << std::endl;

  if (menufilename.find_last_of("txt")!= std::string::npos)
    ReadMenuTXT(menufile);
  else if (menufilename.find_last_of("csv")!= std::string::npos)
    ReadMenuCSV(menufile);
  else
  {
    std::cout << "Can not understand MenuFile "<<menufilename<<"! Please provide TXT or CSV format."<<std::endl;
    return false;
  }

  if (writefiles)
    *outfile <<  "---------------------------- Input Menu -------------------------" <<std::endl << std::endl;

  return true;
}       // -----  end of function L1Menu2016::ReadMenu  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ReadMenuTXT
//  Description:  
// ===========================================================================
bool L1Menu2016::ReadMenuTXT(std::ifstream &menufile)
{
  std::string line;

  while (std::getline(menufile, line))
  {
    if (line.empty()) continue;
    if (writefiles)
      *outfile << line <<std::endl;
    if (line.at(0) == '#')
      continue;
    if (line.at(0) == '%')
    {
      std::cout << "Parsing config from menu: " << line  << std::endl;
      ParseConfig(line);
      continue;
    }

    std::size_t commentpos = line.find_first_of("#");
    std::string goodline = "";
    std::string comline = "";

    if (commentpos != std::string::npos)
    {
      goodline = line.substr(0, commentpos);
      comline = line.substr(commentpos, line.length() - commentpos);
    }
    else
      goodline = line;
    std::istringstream iss(goodline);

    std::string seed;
    int bit;
    int prescale;
    std::string pog, pag;

    iss >> seed >> bit >> prescale >> pog >> pag;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Form L1Seed ~~~~~
    L1Seed temp;
    temp.name = seed;
    vL1Seed.push_back(seed);
    temp.bit = bit;
    temp.comment = comline;
    temp.prescale = prescale;

    if (L1Config["doCompuGT"] || L1Config["SetNoPrescale"] )
      temp.prescale = 1;

    if (L1Config["IgnorePrescale"] && temp.prescale > 1 )
      temp.prescale = 0;

    if (pog.length() != 0)
      temp.POG = TokenGroups(pog);
    if (pag.length() != 0)
      temp.PAG = TokenGroups(pag);
    mL1Seed[seed] = temp;
  }
  return true;
}       // -----  end of function L1Menu2016::ReadMenuTXT  -----


// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ReadMenuCSV
//  Description:  To be finished
// ===========================================================================
bool L1Menu2016::ReadMenuCSV(std::ifstream &menufile)
{
  // Get the first line
  std::string line;
  while (std::getline(menufile, line))
  {
    line.erase( std::remove(line.begin(), line.end(), '\r'), line.end() );
    if (line.empty()) continue;
    if (line.find_first_not_of(", ") == std::string::npos) continue;
    if (line.at(0) == '#')
      continue;
    break; // Get the first line
  }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Getting the header  ~~~~~
  boost::char_separator<char> sep(",", "", boost::keep_empty_tokens);
  tokenizer tokens(line, sep);
  std::map<int, std::string> premap;
  std::map<int, std::string> infomap;
  int j = 0;
  for(auto i : tokens)
  {
    i.erase(std::remove_if(i.begin(), i.end(), ::isspace), i.end());
    try
    {
      boost::lexical_cast<double>(i);
      premap[j] = i;
    }
    catch (const boost::bad_lexical_cast &) {
      infomap[j] = i;
    };
    j++;
  }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Select the prescale column ~~~~~
  int targetcol = -1;
  if (premap.size() == 0)
  {
    std::cout << "No prescale column found in " << menufilename<<". Exiting"  << std::endl;
    exit(1);
  }

  if (premap.size() == 1)
  {
    targetcol = premap.begin()->first;
  }
  else{
    if (premap.size() > 1 && L1ConfigStr["SelectCol"] == "")
    {
      std::cout << "Select prescale columns from: ";
      for(const auto &i : premap)
        std::cout << i.second <<", ";
      std::cout << std::endl;
      exit(1);
    }

    for(const auto &i : premap)
    {
      if (i.second == L1ConfigStr["SelectCol"] )
      {
        targetcol = i.first;
        break;
      }
    }
    if (targetcol == -1)
    {
      std::cout << "Can not find " << L1ConfigStr["SelectCol"] <<" in ";
      for(const auto &i : premap)
        std::cout << i.second <<", ";
      std::cout <<"Exiting" << std::endl;
      exit(1);
    }
  }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Finally reading in the menu ~~~~~
  std::vector<int> outidx;
  std::stringstream ss;
  ss << std::left;
  infomap[targetcol] = "Prescale";
  for(const auto &k : infomap)
  {
    tokenizer::iterator t = tokens.begin();
    std::advance(t, k.first);
    std::string it = std::regex_replace(*t, std::regex("^ +| +$|( ) +"), "$1");
    ss << it <<" "; 
  }
  if (writefiles)
  {
    assert(outfile != nullptr);
    *outfile << ss.str() <<std::endl;
  }

  while (std::getline(menufile, line))
  {
    line.erase( std::remove(line.begin(), line.end(), '\r'), line.end() );
    if (line.empty()) continue;
    if (line.find_first_not_of(", ") == std::string::npos) continue;
    if (line.at(0) == '#')
      continue;
    if (line.at(0) == '%')
    {
      std::cout << "Parsing config from menu: " << line  << std::endl;
      ParseConfig(line);
      continue;
    }

    ss.str("");
    L1Seed temp;

    tokenizer tokens(line, sep);
    for(const auto &k : infomap)
    {
      tokenizer::iterator t = tokens.begin();
      std::advance(t, k.first);
      std::string it = std::regex_replace(*t, std::regex("^ +|\t+$| +$|( ) +"), "$1");
      
      if (k.second == "n")
      {
        ss << std::setw(4) << it <<" "; 
        try
        {
          temp.bit = boost::lexical_cast<int>(it);
        }
        catch (const boost::bad_lexical_cast &)
        {
          std::cout << "Can't cast bit " << it<< " to int type in line: " << line << std::endl;
        }
      }
      if (k.second == "L1AlgoName")
      {
        ss << std::setw(65) << it <<" "; 
        temp.name = it;
        vL1Seed.push_back(it);
      }
      if (k.second == "Prescale")
      {
        ss << std::setw(5) << it <<" "; 
        try
        {
          temp.prescale = boost::lexical_cast<int>(it);
        }
        catch (const boost::bad_lexical_cast &)
        {
          std::cout << "Can't cast prescale " << it<< " to int type in line: " << line <<"; set to disable" << std::endl;
          temp.prescale = 0;
        }
      }
      if (k.second == "Comment")
      {

        ss << it <<" "; 
        temp.comment = it;
      }
      if (k.second == "POG")
      {
        ss << std::setw(15) << it <<" "; 
        temp.POG = TokenGroups(it);
      }
      if (k.second == "PAG")
      {
        ss << std::setw(15) << it <<" "; 
        temp.PAG = TokenGroups(it);
      }

      if (L1Config["doCompuGT"] || L1Config["SetNoPrescale"] )
        temp.prescale = 1;

      if (L1Config["IgnorePrescale"] && temp.prescale > 1 )
        temp.prescale = 0;
    }

    if (writefiles)
      *outfile << ss.str() <<std::endl;
    mL1Seed[temp.name] = temp;
  }

  return true;
}       // -----  end of function L1Menu2016::ReadMenuCSV  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::FormPrescaleColumns
//  Description:  
// ===========================================================================
bool L1Menu2016::FormPrescaleColumns()
{
  std::vector<std::string> VarSeeds;
  for(auto &seed : mL1Seed)
  {
    if (seed.second.prescale < 0)
    {
      VarSeeds.push_back(seed.first);
      std::cout << "Varing " << seed.first << std::endl;
    }
  }
 
  const int varSize = VarSeeds.size();
  int varCounts = 1 << varSize;
  for (int i = 0; i < varCounts; ++i)
  {
    std::map<std::string, L1Seed> tempL1Seed = mL1Seed;
    for (int j = 0; j < varSize; ++j)
    {
      tempL1Seed.at(VarSeeds.at(j)).prescale = (i & 1 << j) > 0 ; 
    }
    ColumnMap[i] = new PreColumn(i, tempL1Seed);
    ColumnMap[i]->PassRelation( vL1Seed, BitMap, POGMap, PAGMap);
  }

  std::cout << "In total ColumnMap "<< ColumnMap.size() << std::endl;
  return true;
}       // -----  end of function L1Menu2016::FormPrescaleColumns  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ReadFilelist
//  Description:  
// ===========================================================================
bool L1Menu2016::OpenWithList(std::string filelist)
{
  L1Ntuple::SelectTree(L1Config["UseUnpackTree"]);
  if (filelist.find(".root") != std::string::npos)
  {
    L1Ntuple::Open(filelist);
    return false;
  }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ OpenNtupleList ~~~~~

  std::ifstream flist(filelist.c_str());
  if (!flist)
  {
    std::cout << "File "<<filelist<<" is not found !"<<std::endl;
    return false;
  }

  std::string line;
  while (std::getline(flist, line))
  {
    if (line.empty()) continue;
    if (line.at(0) == '#')
      continue;
    if (line.at(0) == '%')
    {
      std::cout << "Parsing config from ntuple list: " << line  << std::endl;
      ParseConfig(line);
      continue;
    }
    if (!flist.fail())
    {
      listNtuples.push_back(line);
    }
  }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ OpenNtupleList ~~~~~
  if (!CheckFirstFile())      exit(0);
  if (!OpenWithoutInit())     exit(0);

  std::cout.flush();
  std::cout<<"Going to init the available trees..."<<std::endl;
  std::cout.flush();
  Init();

  return true;
}       // -----  end of function L1Menu2016::ReadFilelist  -----


// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseConfig
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseConfig(const std::string line)
{
  if (line.find('=') == std::string::npos)
  {
    std::cout << "Wrong config: " << line << std::endl;
    return false;
  }

  std::istringstream iss(line);

  std::string key, sign;
  double value;

  iss >> sign >> key >> sign >> value;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Parse for double inputs ~~~~~
  bool parDouble = iss.fail();
  if (!parDouble && L1Config.find(key) != L1Config.end())
  {
    L1Config[key] = value;
    return true;
  }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Parse for string inputs ~~~~~
  bool parString = true;
  std::istringstream isss(line);
  key=sign="";
  std::string valstr="";

  if (parDouble)
  {
    isss >> sign >> std::ws >> key >>std::ws >> sign >>std::ws >> valstr;
    parString = isss.fail();
  }

  if (!parString && L1ConfigStr.find(key) != L1ConfigStr.end())
  {
    L1ConfigStr[key] = valstr;
    return true;
  }

  if (parDouble && parString)
    std::cout<<"\033[0;31mCan't parse config:\033[0m "<<line<< std::endl; 
  else
    std::cout << "Not reconfiguzed config key " << key<< std::endl;
  
  return false;
}       // -----  end of function L1Menu2016::ParseConfig  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::GetRunConfig
//  Description:  Get the running time config from command line
// ===========================================================================
bool L1Menu2016::GetRunConfig(std::map<std::string, float> &config, 
    std::map<std::string, std::string> &configstr)
{

  for(auto c : config)
  {
    if (L1Config.find(c.first) != L1Config.end())
    {
      L1Config[c.first] = c.second;
    }
  }

  for(auto c : configstr)
  {
    if (L1ConfigStr.find(c.first) != L1ConfigStr.end())
    {
      L1ConfigStr[c.first] = c.second;
    }
  }
  return true;
}       // -----  end of function L1Menu2016::GetRunConfig  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::PrintConfig
//  Description:  
// ===========================================================================
bool L1Menu2016::PrintConfig() const
{
  std::cout << "Printing configuration ___________________________ " << std::endl;
  for(auto &x: L1Config)
    std::cout << std::setw(20) <<x.first <<" : " << x.second << std::endl;
  for(auto &x: L1ConfigStr)
    std::cout << std::setw(20) <<x.first <<" : " << x.second << std::endl;
  std::cout << "Printed configuration ============================ " << std::endl;
  return true;
}       // -----  end of function L1Menu2016::PrintConfig  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::PreLoop
//  Description:  
// ===========================================================================
bool L1Menu2016::PreLoop(std::map<std::string, float> &config, std::map<std::string, std::string> &configstr)
{
  GetRunConfig(config, configstr);
  OpenWithList(tuplefilename);

  //Prepare Menu
  std::cout << "Preparing menus __________________________________ " << std::endl;
  ReadMenu();
  BuildRelation();
  L1SeedFunc();
  FormPrescaleColumns();

  PrintConfig();
  BookHistogram();
  
  if (writeplots)
  {
    GlobalAlgBlk *l1uGTsel_ = l1uGT_;
    TChain       *fl1uGTsel = fl1uGT;
    
    if (L1Config["UseUnpackTree"]){
      l1uGTsel_ = l1unpackuGT_;
      fl1uGTsel = fl1unpackuGT;      
    }
    
    l1Plot = new L1Plot(outrootfile, event_, upgrade_, recoJet_,
			recoSum_, recoEle_, recoMuon_, recoTau_, recoFilter_, l1CaloTower_, recoVtx_, l1uGTsel_);
    l1Plot->SetTodo(L1Config);
    l1Plot->PreRun(&L1Event, &mL1Seed, L1Ntuple::GetuGTAlias(fl1uGTsel));
  }

  if (L1Config["doPrintPU"])
  {
    ReadDataPU();
  }
    
  if (L1Config["doTnPMuon"])
  {
    l1TnP = new L1TnP(outrootfile, event_, upgrade_, recoJet_,
		      recoSum_, recoEle_, recoMuon_, recoTau_, recoFilter_, l1CaloTower_, recoVtx_, l1uGT_);
    if (L1Config["doTnPMuon"])
      l1TnP->DoMuonTnP();
  }

  if (l1unpackuGT_ != NULL)
  {
    l1unpackuGT = new L1uGT( outrootfile, event_, l1unpackuGT_, &L1Event, &mL1Seed);
    l1unpackuGT->GetTreeAlias(L1Ntuple::GetuGTAlias(fl1unpackuGT));
  }

  if (L1Config["doCompuGT"] || L1Config["UseuGTDecision"] || L1Config["doPlotuGt"])
  {
    assert(l1uGT_ != NULL);
    l1uGT = new L1uGT( outrootfile, event_, l1uGT_, &L1Event, &mL1Seed);
    l1uGT->GetTreeAlias(L1Ntuple::GetuGTAlias(fl1uGT));
  }


  if (L1Config["SetMuonER"] != -1) SetMuonER(L1Config["SetMuonER"]);
  if (L1Config["UseUpgradeLyr1"] != -1) SetUseUpgradeLyr1(L1Config["UseUpgradeLyr1"]);
  if (L1Config["UseL1CaloTower"] != -1) SetUseL1CaloTower(L1Config["UseL1CaloTower"]);

  if (L1ConfigStr["SelectLS"] != "") 
    ParseRanges("SelectLS", pLS);

  if (L1ConfigStr["SelectBX"] != "") 
    ParseRanges("SelectBX", pBX);
  return true;
}       // -----  end of function L1Menu2016::PreLoop  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::GetL1Event
//  Description:  
// ===========================================================================
bool L1Menu2016::GetL1Event()
{
  L1Event = {};

  //Jet
  L1AlgoFactory::SingleJetPt(L1Event.JetPt,false);
  L1AlgoFactory::SingleJetPt(L1Event.JetCenPt,true);

  //EG
  L1AlgoFactory::SingleEGPt(L1Event.EGPt,false, false);
  L1AlgoFactory::SingleEGPt(L1Event.EGerPt,false, true);
  L1AlgoFactory::SingleEGPt(L1Event.IsoEGPt,true, false);
  L1AlgoFactory::SingleEGPt(L1Event.IsoEGerPt,true, true);

  //Tau
  L1AlgoFactory::SingleTauPt(L1Event.TauPt, false, false); 
  L1AlgoFactory::SingleTauPt(L1Event.TauerPt, true, false);
  L1AlgoFactory::SingleTauPt(L1Event.IsoTauPt, false, true);

  //Mu
  L1AlgoFactory::SingleMuPt(L1Event.MuPt, false, 2);
  L1AlgoFactory::SingleMuPt(L1Event.MuOpenPt, false, 0);
  L1AlgoFactory::SingleMuPt(L1Event.MuerPt, true, 2);

  //Sum
  if (L1Config["SumJetET"] != 0 || L1Config["SumJetEta"] != 999)
  {
    CalLocalHT(L1Event.HTT, false);
    CalLocalHT(L1Event.HTTHF, true);
    CalLocalHTM(L1Event.HTM);
  } else {
    L1AlgoFactory::HTTVal(L1Event.HTT);
    L1AlgoFactory::HTTHFVal(L1Event.HTTHF); // Not in L1Ntuple yet
    L1AlgoFactory::HTMVal(L1Event.HTM);
  }

  if (L1Config["UseL1CaloTower"])
  {
    CalLocalETM(L1Event.ETM, false);
    CalLocalETM(L1Event.ETMHF, true);
  }
  else
  {
    L1AlgoFactory::ETMVal(L1Event.ETM);
    L1AlgoFactory::ETMHFVal(L1Event.ETMHF);
  }

  L1AlgoFactory::ETTVal(L1Event.ETT);

  // Mulit
  float dummy = 0;
  L1AlgoFactory::DoubleMuPt(L1Event.doubleMuPt1,L1Event.doubleMuPt2, true, false);
  L1AlgoFactory::Onia2015Pt(L1Event.oniaMuPt1, L1Event.oniaMuPt2,true, false, 18);
  L1AlgoFactory::DoubleJetPt(L1Event.dijetPt1,L1Event.dijetPt2);
  L1AlgoFactory::DoubleJetPt(L1Event.diCenjetPt1,L1Event.diCenjetPt2,true);
  L1AlgoFactory::DoubleTauJetEta2p17Pt(dummy,L1Event.ditauPt, false);
  L1AlgoFactory::DoubleTauJetEta2p17Pt(dummy,L1Event.diIsotauPt, true);
  L1AlgoFactory::QuadJetPt(dummy,dummy,dummy,L1Event.quadjetPt);
  L1AlgoFactory::QuadJetPt(dummy,dummy,dummy,L1Event.quadjetCPt,true);
  L1AlgoFactory::DoubleEGPt(L1Event.diEG1,L1Event.diEG2,false);
  L1AlgoFactory::DoubleEGPt(L1Event.diIsolEG1,L1Event.diIsolEG2,true);

  return true;
}       // -----  end of function L1Menu2016::GetL1Event  -----
// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::Loop
//  Description:  
// ===========================================================================
bool L1Menu2016::Loop()
{
  ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Initial ~~~~~
  //Int_t nevents = fChain->GetEntriesFast();//GetEntries();
  unsigned int currentLumi(-1);
  nZeroBiasevents = 0.;
  int i = -1;
  nLumi.clear();
  bool skipLS = false;

  while(true)
  {
    i++;
    Long64_t ientry = LoadTree(i); 
    if (ientry < 0) break;
    GetEntry(i);
    if (L1Config["maxEvent"] != -1 && i > L1Config["maxEvent"]) break;

    if (event_ != NULL )
    {
      if (L1Config["SelectRun"] != -1 && event_->run != L1Config["SelectRun"])
        continue;

      if (L1Config["SelectEvent"] != -1 && event_->event != L1Config["SelectEvent"])
        continue;

      if(event_ -> lumi != currentLumi){
        currentLumi = event_ -> lumi;
        skipLS      = CheckLS(currentLumi);
        if (!skipLS)
          nLumi.insert(currentLumi);
      } 

      if (L1ConfigStr["SelectLS"] != "" && skipLS)
        continue;

      if (L1ConfigStr["SelectBX"] != "" && CheckBX(event_->bx) )
        continue;

      if (L1Config["doScanLS"])
        continue;
    }

    if (i % 200000 == 0)
    {
      std::cout << "Processed " << i << " events." << std::endl;
    }

    //Use Final decision by default, unless for PlotLS
    if (l1unpackuGT != NULL && !l1unpackuGT->GetuGTDecision("L1_ZeroBias", L1Config["doPlotLS"])) 
      continue;

    nZeroBiasevents++;

    GetL1Event();
    RunMenu();

    if (L1Config["doPlotLS"])
      FillLumiSection(currentLumi);

    if (L1Config["doPrintPU"] && event_ != NULL)
      FillPileUpSec();

    if (l1Plot != NULL)
      l1Plot->RunPlot();

    if (l1TnP != NULL)
      l1TnP->RunTnP();

    if (L1Config["doCompuGT"])
      l1uGT->CompEvents();
  }

  std::cout << "Total Event: " << i <<" ZeroBias Event: " << nZeroBiasevents << std::endl;
  return true;
}       // -----  end of function L1Menu2016::Loop  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::PostLoop
//  Description:  
// ===========================================================================
bool L1Menu2016::PostLoop()
{
  scale = CalScale(0, 0, true);

  std::cout << "==================================== Summary" << std::endl;

  for(auto col : ColumnMap)
  {
    col.second->CalRate(scale);
    col.second->FillDefHist1D(scale);
    col.second->FillDefHist2D(scale);
  }
  
  if (ColumnMap.size() == 1)
    PrintRates(std::cout);
  else{
    for(auto col : ColumnMap)
      col.second->PrintMenuRate(scale);
  }
  if (writefiles)
    PrintRates(*outfile);
  PrintCSV(*outcsv);

  if (l1Plot != NULL)
  {
    l1Plot->PostRun(scale);
    if (L1Config["doPlotLS"])
      l1Plot->PlotRatePerLS(L1LSCount, L1Config["nBunches"]);
  }

  if (l1TnP != NULL)
    l1TnP->PostRun();

  if (L1Config["doPrintPU"])
  {
    PrintPUCSV();
  }
  WriteHistogram();
  return true;
}       // -----  end of function L1Menu2016::PostLoop  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::PrintPUCSV
//  Description:  
// ===========================================================================
bool L1Menu2016::PrintPUCSV()
{
  //const int nBunches = 2736;
  std::fstream pucsv (outputdir + "/" + outputname+"_PU" +".csv", std::fstream::out );

  // L1Seeds
  std::vector<std::string> csvout;
  for(auto col : ColumnMap)
  {
    col.second->PrintPUCSV(csvout);
  }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Writing out to file ~~~~~
  for(auto i : csvout)
  {
    pucsv << i <<std::endl;
  }
  pucsv.close();

  return true;
}       // -----  end of function L1Menu2016::PrintPUCSV  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::PrintRates
//  Description:  
// ===========================================================================
bool L1Menu2016::PrintRates(std::ostream &out)
{
  for(auto col : ColumnMap)
  {
    col.second->PrintRates(out, scale);
  }
  return true;
}       // -----  end of function L1Menu2016::PrintRates  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::BuildRelation
//  Description:  
// ===========================================================================
bool L1Menu2016::BuildRelation()
{
  int misbit= -1;
  for(auto &l1 : mL1Seed)
  {

    // bit
    if (BitMap.find(l1.second.bit) != BitMap.end())
    {
      std::cout << "Duplicate bit number at " << l1.second.bit <<" for " 
        << BitMap[l1.second.bit] <<" and " << l1.first 
        << "; setting " << l1.first  << " to bit " << misbit << std::endl;
      l1.second.bit = misbit;
      BitMap[l1.second.bit] = l1.first;
      misbit--;
    } else BitMap[l1.second.bit] = l1.first;

    for(auto &pog : l1.second.POG)
    {
      if (POGMap.find(pog) != POGMap.end())
      {
        POGMap[pog] = {};
        assert(POGMap[pog].size() == 0);
      } else POGMap[pog].push_back(l1.second.bit);
    }

    for(auto &pag : l1.second.PAG)
    {
      if (PAGMap.find(pag) != PAGMap.end())
      {
        PAGMap[pag] = {};
        assert(PAGMap[pag].size() == 0);
      } else PAGMap[pag].push_back(l1.second.bit);
    }
  }

  return true;
}       // -----  end of function L1Menu2016::BuildRelation  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::L1SeedFunc
//  Description:  
// ===========================================================================
bool L1Menu2016::L1SeedFunc()
{
#ifdef UTM_MENULIB
  addFuncFromName(L1SeedFun, upgrade_, l1CaloTower_);
  //addFuncFromName(L1SeedFun, upgrade_);
#endif
    
  for(auto &L1Seed : mL1Seed)
  {
    if (L1SeedFun.find(L1Seed.first) != L1SeedFun.end())
      continue;

    if(ParseL1Seed(L1Seed.first))
      continue;

    std::cout << "No function call for " << L1Seed.first <<"; setting to no fire"<< std::endl;
  }

  return true;
}       // -----  end of function L1Menu2016::L1SeedFunc  -----


// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::CheckL1Seed
//  Description:  /* cursor */
// ===========================================================================
bool L1Menu2016::CheckL1Seed(const std::string L1Seed)
{
  if (L1SeedFun.find(L1Seed) != L1SeedFun.end())
  {
    return L1SeedFun[L1Seed]();
  }
  return false;
}       // -----  end of function L1Menu2016::CheckL1Seed  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::RunMenu
//  Description:  
// ===========================================================================
bool L1Menu2016::RunMenu()
{
  for(auto col : ColumnMap)
  {
    col.second->EventReset();
  }

  for(auto& seed: mL1Seed)
  {
    bool IsFired = false;
    if (L1Config["UseuGTDecision"])
    {
      assert(l1uGT != NULL);
      IsFired = l1uGT->GetuGTDecision(seed.first);
    }
    else
      IsFired = CheckL1Seed(seed.first);

    for(auto col : ColumnMap)
    {
      col.second->InsertInMenu(seed.first, IsFired);
    }
  }

  for(auto col : ColumnMap)
  {
    col.second->CheckCorrelation();
  }

  return true;
}       // -----  end of function L1Menu2016::RunMenu  -----


// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::CalScale
//  Description:  
// ===========================================================================
double L1Menu2016::CalScale(int nEvents_, int nBunches_, bool print) 
{
  double scale = 0.0;
  int nEvents = nEvents_ == 0 ? nZeroBiasevents : nEvents_;
  double nBunches = nBunches_ == 0 ?  L1Config["nBunches"] : nBunches_;

  if (L1Config["nBunches"] < 0)
  {
    scale = (-1.*L1Config["nBunches"])/(nLumi.size()*23.31);      
    if (print)
      std::cout << "Scale by "   << "("<< -1. * L1Config["nBunches"] <<")/(nLumi*23.31) with nLumi = " << nLumi.size()      << std::endl;
  } else {
    scale = 11246.; // ZB per bunch in Hz
    //scale /= nZeroBiasevents*1000.; // in kHz
    scale /= nEvents; // in Hz
    scale *= nBunches;
    if (print)
      std::cout << "Scale by "   << " 11246 / nZeroBiasevents * NumberOfBunches, with nZeroBiasevents = " 
        << nEvents    <<" NumberOfBunches = " << nBunches << std::endl;
  }
  return scale;
}       // -----  end of function L1Menu2016::CalScale  -----


// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::SetOutputName
//  Description:  
// ===========================================================================
std::string L1Menu2016::SetOutputName() const
{
  boost::filesystem::path menupath(menufilename);
  boost::filesystem::path flistpath(tuplefilename);
  std::stringstream ss;
  ss << menupath.stem().string() << "-" << flistpath.stem().string();
  return ss.str();
}       // -----  end of function L1Menu2016::SetOutputName  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseL1Seed
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseL1Seed(const std::string SeedName)
{
  if (SeedName.find("L1_") == std::string::npos)
  {
    return false;
  }

  if (ParseSingleObject(SeedName)) return true;

  // Jets
  if (ParseDoubleJet(SeedName)) return true;
  if (ParseTripleJetVBF(SeedName)) return true;
  if (ParseQuadJet(SeedName)) return true;
  // EG
  if (ParseDoubleEG(SeedName)) return true;
  if (ParseTripleEG(SeedName)) return true;
  // Tau
  if (ParseDoubleTau(SeedName)) return true;

  // Mu_Tau
  if (ParseMuerTauer(SeedName)) return true;
  // Mu_EG
  if (ParseMuEG(SeedName)) return true;
  // Mu_Sum
  if (ParseMuSum(SeedName)) return true;

  // EG_Sum
  if (ParseEGSum(SeedName)) return true;

  if (ParseEGStrategy(SeedName)) return true;

  if (ParseETMJetdPhi(SeedName)) return true;

  if (ParseComplexSingleMu(SeedName)) return true;
  // EGMass
  //if (ParseMultiEGMass(SeedName)) return true;

  return false;
}       // -----  end of function L1Menu2016::ParseL1Seed  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseSingleObject
//  Description:  /* cursor */
// ===========================================================================
bool L1Menu2016::ParseSingleObject(const std::string SeedName)
{
  std::string L1object ="";
  std::string postfix = "";
  int pt = -10;

  boost::char_separator<char> sep("_");
  tokenizer tokens(SeedName, sep);
  if (std::distance(tokens.begin(), tokens.end()) < 2) return false;
  boost::tokenizer<boost::char_separator<char> >::iterator tokenit = tokens.begin();
  tokenit++;
  std::string Seedtoken(*tokenit);

  std::smatch base_match;
  std::regex integerobj("Single([^0-9]+)([0-9]+)(er|er2p1|)");
  std::regex integerSum("(ETMHF|ETM|HTT|ETT|HTM)([0-9]+)(er|)");
  if (std::regex_match(Seedtoken, base_match, integerobj))
  {
	// The first sub_match is the whole string; the next
	// sub_match is the first parenthesized expression.
    L1object = base_match[1].str();
    pt = std::stoi(base_match[2].str(), nullptr);
    postfix = base_match[3].str();
  }else if (std::regex_match(Seedtoken, base_match, integerSum))
  {
	// The first sub_match is the whole string; the next
	// sub_match is the first parenthesized expression.
    L1object = base_match[1].str();
    pt = std::stoi(base_match[2].str(), nullptr);
    postfix = base_match[3].str();
  } else if(Seedtoken == "SingleMuOpen")
  {
    L1object = "Mu";
    postfix = "Open";
    pt = 0;
  }
    

  L1object += postfix;
  mL1Seed[SeedName].singleObj = L1object;
  std::vector<std::function<bool()>> funs;

  //std::cout <<  std::distance(tokenit, tokens.end())<< std::endl;
  if (std::distance(tokenit, tokens.end()) > 1) 
  {
    tokenit++;
    Seedtoken = *tokenit;
    if (Seedtoken == "NotBptxOR" || Seedtoken == "BptxAND")
    {
      funs.push_back(ParseBptx(Seedtoken));
    }
    else return false;
  }

  if (L1ObjectMap.find(L1object) != L1ObjectMap.end())
  {
    //funs.push_back(std::bind(&SingleObjPt, L1ObjectMap[L1object], pt));
    // No idea for the funs vector 
    L1SeedFun[SeedName] = std::bind(&SingleObjPt, L1ObjectMap[L1object], pt);
  } else return false;

  return true;
}       // -----  end of function L1Menu2016::ParseSingleObject  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseBptx
//  Description:  /* cursor */
// ===========================================================================
std::function<bool()> L1Menu2016::ParseBptx(const std::string /*Seedtoken*/)
{
  return [](){return true;};
}       // -----  end of function L1Menu2016::ParseBptx  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseDoubleJet
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseDoubleJet(const std::string& SeedName)
{
  std::smatch base_match;
  std::regex integer("L1_DoubleJet(C|)([0-9]+)(er3p0|)");
  if (std::regex_match(SeedName, base_match, integer))
  {
    bool isCentral = base_match.length(1) == 1 || (base_match.length(3) == 5 );
    unsigned int pt = std::stoi(base_match[2].str(), nullptr);
    L1SeedFun[SeedName] = std::bind(&L1AlgoFactory::DoubleJet, this, pt, pt, isCentral);
    return true;
  }
  else return false;
}       // -----  end of function L1Menu2016::ParseDoubleJet  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseDoubleTau
//  Description:  /* cursor */
// ===========================================================================
bool L1Menu2016::ParseDoubleTau(const std::string& SeedName) 
{
  std::smatch base_match;
  std::regex integer("L1_Double(Iso|)Tau([0-9]+)er(2p1|)");
  if (std::regex_match(SeedName, base_match, integer))
  {
    bool isIso = base_match.length(1) == 3;
    unsigned int pt = std::stoi(base_match[2].str(), nullptr);
    L1SeedFun[SeedName] = std::bind(&L1AlgoFactory::DoubleTauJetEta2p17, this, pt, pt, isIso);
    return true;
  }
  else return false;
}       // -----  end of function L1Menu2016::ParseDoubleTau  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseTripleJetVBF
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseTripleJetVBF(const std::string& SeedName)
{
  const int jetclass = 0; 
  std::smatch base_match;
  std::regex integer("L1_TripleJet_([0-9]+)_([0-9]+)_([0-9]+)_VBF");
  if (std::regex_match(SeedName, base_match, integer))
  {
    L1SeedFun[SeedName] = std::bind(&L1AlgoFactory::TripleJet_VBF, this, 
        std::stoi(base_match[1].str(), nullptr),
        std::stoi(base_match[2].str(), nullptr),
        std::stoi(base_match[3].str(), nullptr), jetclass);
    return true;
  }
  else return false;
}       // -----  end of function L1Menu2016::ParseTripleJetVBF  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseQuadJet
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseQuadJet(const std::string& SeedName)
{
  int pt1 = -10;
  int pt2 = -10;
  int pt3 = -10;
  int pt4 = -10;
  bool isCentral = false;

  std::smatch base_match;
  std::regex integer_sys("L1_QuadJet([C]*)([0-9]+)(er3p0|)");
  std::regex integer_asys("L1_QuadJet([C]*)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)");
  if (std::regex_match(SeedName, base_match, integer_sys))
  {
    isCentral = base_match.length(1) == 1 || base_match.length(3) == 5;
    pt1 = std::stoi(base_match[2].str(), nullptr);
    pt2 = pt3 = pt4 = pt1;
  } else if (std::regex_match(SeedName, base_match, integer_asys))
  {
    isCentral = base_match.length(1) == 1;
    pt1 = std::stoi(base_match[2].str(), nullptr);
    pt2 = std::stoi(base_match[3].str(), nullptr);
    pt3 = std::stoi(base_match[4].str(), nullptr);
    pt4 = std::stoi(base_match[5].str(), nullptr);
  }

  if (pt1 != -10 && pt2 != -10 && pt3 != -10 && pt4 != -10)
  {
    L1SeedFun[SeedName] = std::bind(&L1AlgoFactory::QuadJet, this, pt1, pt2, pt3, pt4, isCentral);
    return true;
  }
  else return false;
}       // -----  end of function L1Menu2016::ParseQuadJet  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseDoubleEG
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseDoubleEG(const std::string& SeedName)
{
  int pt1 = -10;
  int pt2 = -10;
  bool isIso = false;

  std::smatch base_match;
  std::regex integer_sys("L1_Double(Iso|)EG([0-9]+)");
  std::regex integer_asys("L1_Double(Iso|)EG_([0-9]+)_([0-9]+)");
  if (std::regex_match(SeedName, base_match, integer_sys))
  {
    isIso = base_match.length(1) == 3;
    pt1 = std::stoi(base_match[2].str(), nullptr);
    pt2 = pt1;
  }else if (std::regex_match(SeedName, base_match, integer_asys))
  {
    isIso = base_match.length(1) == 3;
    pt1 = std::stoi(base_match[2].str(), nullptr);
    pt2 = std::stoi(base_match[3].str(), nullptr);
  }

  if (pt1 != -10 && pt2 != -10)
  {
    L1SeedFun[SeedName] = std::bind(&L1AlgoFactory::DoubleEG, this, pt1, pt2, isIso);
    return true;
  }
  else
    return false;
}       // -----  end of function L1Menu2016::ParseDoubleEG  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseTripleEG
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseTripleEG(const std::string& SeedName)
{
  int pt1 = -10;
  int pt2 = -10;
  int pt3 = -10;

  std::smatch base_match;
  std::regex integer_sys("L1_TripleEG([0-9]+)");
  std::regex integer_asys("L1_TripleEG_([0-9]+)_([0-9]+)_([0-9]+)");
  if (std::regex_match(SeedName, base_match, integer_sys))
  {
    pt1 = std::stoi(base_match[1].str(), nullptr);
    pt2 = pt1;
    pt3 = pt1;
  }else if (std::regex_match(SeedName, base_match, integer_asys))
  {
    pt1 = std::stoi(base_match[1].str(), nullptr);
    pt2 = std::stoi(base_match[2].str(), nullptr);
    pt3 = std::stoi(base_match[3].str(), nullptr);
  }

  if (pt1 != -10 && pt2 != -10 && pt3 != -10)
  {
    L1SeedFun[SeedName] = std::bind(&L1AlgoFactory::TripleEG, this, pt1, pt2, pt3);
    return true;
  }
  else
    return false;
}       // -----  end of function L1Menu2016::ParseTripleEG  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseEGSum
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseEGSum(const std::string& SeedName)
{
  int EGpt = -10;
  int Sumpt = -10;
  bool isIsoEG= false;

  std::smatch base_match;
  std::regex integerEGerHTT("L1_(Iso|)EG([0-9]+)er(2p1|)_HTT([0-9]+)(er|)");
  if (std::regex_match(SeedName, base_match, integerEGerHTT))
  {
    isIsoEG = base_match.length(1) == 3;
    EGpt =  std::stoi(base_match[2].str(), nullptr);
    Sumpt = std::stoi(base_match[4].str(), nullptr);
    L1SeedFun[SeedName] = std::bind(&L1AlgoFactory::SingleEG_Eta2p1_HTT, this, EGpt, Sumpt, isIsoEG);
    return true;
  }

  return false;
}       // -----  end of function L1Menu2016::ParseEGSum  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseEGStrategy
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseEGStrategy(const std::string & SeedName)
{
  int pt = -10;

  std::smatch base_match;
  std::regex EGStrategy("L1_StrategyEG([0-9]+)");
  if (std::regex_match(SeedName, base_match, EGStrategy))
  {
    pt =  std::stoi(base_match[1].str(), nullptr);
    L1SeedFun[SeedName] = 
      boost::bind(&SingleObjPt, L1ObjectMap["EG"], pt) ||
      boost::bind(&SingleObjPt, L1ObjectMap["IsoEG"], pt-2) ||
      boost::bind(&SingleObjPt, L1ObjectMap["IsoEGer"], pt-4);
    return true;
  }

  return false;
}       // -----  end of function L1Menu2016::ParseEGStrategy  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseETMJetdPhi
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseETMJetdPhi(const std::string & SeedName)
{
  int ETMpt = -10;
  int Jetpt = -10;
  bool isCentral = false;

  std::smatch base_match;
  std::regex ETMJetdPhi("L1_ETM([0-9]+)_Jet([C]*)([0-9]+)_dPhi_Min0p4");
  if (std::regex_match(SeedName, base_match, ETMJetdPhi))
  {
    ETMpt =  std::stoi(base_match[1].str(), nullptr);
    isCentral = base_match.length(2) == 1;
    Jetpt =  std::stoi(base_match[3].str(), nullptr);
    L1SeedFun[SeedName] = std::bind(&L1AlgoFactory::ETM_Jet, this, ETMpt, Jetpt, isCentral);
    return true;
  }

  return false;
}       // -----  end of function L1Menu2016::ParseETMJetdPhi  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseMultiEGMass
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseMultiEGMass(const std::string& SeedName)
{
  
  int pts = -10;
  int pt1 = -10;
  int pt2 = -10;
  int pt3 = -10;
  int pt4 = -10;
  int Mcut = -10;
  bool isIso = false;
  bool isER = false;
  int EGcount = 0;

  std::map<std::string, int> CountMap;
  CountMap["Single"] = 1;
  CountMap["Double"] = 2;
  CountMap["Triple"] = 3;
  CountMap["Quad"] = 4;


  std::smatch base_match;
  std::regex integer("L1_(Single|Double|Triple|Quad)(Iso|)EG([0-9]*)(er|)(_{0,1}[0-9]*)(_{0,1}[0-9]*)(_{0,1}[0-9]*)(_{0,1}[0-9]*)_M([0-9]+)");
  if (std::regex_match(SeedName, base_match, integer))
  {
    EGcount = CountMap[base_match[1].str()];
    isIso = base_match.length(2) != 0;
    if (base_match[3].str() != "")
      pts =std::stoi(base_match[3].str(), nullptr);
    isER = base_match.length(4) != 0;
    if (base_match[5].str() != "")
      pt1 = std::stoi(base_match[5].str().erase(0, 1), nullptr);
    if (base_match[6].str() != "")
      pt2 = std::stoi(base_match[6].str().erase(0, 1), nullptr);
    if (base_match[7].str() != "")
      pt3 = std::stoi(base_match[7].str().erase(0, 1), nullptr);
    if (base_match[8].str() != "")
      pt4 = std::stoi(base_match[8].str().erase(0, 1), nullptr);
    Mcut = std::stoi(base_match[9].str(), nullptr);
  }


  std::vector<int> Ptcuts;
  if (pts > 0) // Set to all legs
  {
    pt1 = pt2 = pt3 = pt4 = pts;
    for (int i = 0; i < 4; ++i)
    {
      if (i <= EGcount)
        Ptcuts.push_back(pts);
      else
        Ptcuts.push_back(-10);
    }
  } else{
    Ptcuts.push_back(pt1);
    Ptcuts.push_back(pt2);
    Ptcuts.push_back(pt3);
    Ptcuts.push_back(pt4);
    std::sort(Ptcuts.begin(), Ptcuts.end(), std::greater<int>());
  }

  assert(EGcount == std::count_if(Ptcuts.begin(), Ptcuts.end(), [](int i){return i > 0;}));
  L1SeedFun[SeedName] = std::bind(&L1AlgoFactory::MultiEGMass, this, Ptcuts.at(0), 
      Ptcuts.at(1), Ptcuts.at(2),  Ptcuts.at(3), Mcut, isIso, isER);

  return true;
}       // -----  end of function L1Menu2016::ParseMultiEGMass  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseCrossMu
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseCrossMu(const std::string& /*SeedName*/)
{
  //std::smatch base_match;
  //std::regex integer("L1_QuadJet([C]*)([0-9]+)");
  //if (std::regex_match(SeedName, base_match, integer))
  //{
    //if (base_match.size() != 2) return false;
    //bool isCentral = base_match.length(1) == 1;
    //unsigned int pt = std::stoi(base_match[2].str(), nullptr);
    //L1SeedFun[SeedName] = std::bind(&L1AlgoFactory::QuadJet, this, pt, pt, pt, pt, isCentral);
    //return true;
  //}
  //else return false;
  return false;
}       // -----  end of function L1Menu2016::ParseCrossMu  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseMuerTauer
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseMuerTauer(const std::string& SeedName)
{
  int Mupt = -10;
  int Taupt = -10;
  bool isIsoTau = false;

  std::smatch base_match;
  std::regex integer("L1_Mu([0-9]+)(er|er2p1)_(Iso|)Tau(Jet|)([0-9]+)(er|er2p1)");
  if (std::regex_match(SeedName, base_match, integer))
  {
    Mupt =  std::stoi(base_match[1].str(), nullptr);
    isIsoTau = base_match.length(3) == 3;
    Taupt = std::stoi(base_match[5].str(), nullptr);
    L1SeedFun[SeedName] = std::bind(&L1AlgoFactory::Muer_TauJetEta2p17, this, Mupt, Taupt, isIsoTau);
    return true;
  }

  return false;
}       // -----  end of function L1Menu2016::ParseMuerTauer  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseMuEG
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseMuEG(const std::string& SeedName)
{
  const int muonQual = 2;
  int Mupt = -10;
  int EGpt = -10;
  bool isIsoEG = false;

  std::smatch base_match;
  std::regex integer("L1_Mu([0-9]+)_(Iso|)EG([0-9]+)");
  if (std::regex_match(SeedName, base_match, integer))
  {
    Mupt =  std::stoi(base_match[1].str(), nullptr);
    isIsoEG = base_match.length(2) == 3;
    EGpt = std::stoi(base_match[3].str(), nullptr);
    L1SeedFun[SeedName] = std::bind(&L1AlgoFactory::Mu_EG, this, Mupt, EGpt, isIsoEG, muonQual);
    return true;
  }

  return false;
}       // -----  end of function L1Menu2016::ParseMuEG  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseMuSum
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseMuSum(const std::string& SeedName)
{
  int Mupt = -10;
  int Sumpt = -10;

  std::smatch base_match;
  std::regex integerMuHTT("L1_Mu([0-9]+)_HTT([0-9]+)(er|)");
  std::regex integerMuerETM("L1_Mu([0-9]+)(er|er2p1)_ETM([0-9]+)");
  if (std::regex_match(SeedName, base_match, integerMuHTT))
  {
    Mupt =  std::stoi(base_match[1].str(), nullptr);
    Sumpt = std::stoi(base_match[2].str(), nullptr);
    L1SeedFun[SeedName] = std::bind(&L1AlgoFactory::Mu_HTT, this, Mupt, Sumpt);
    return true;
  }

  if (std::regex_match(SeedName, base_match, integerMuerETM))
  {
    Mupt =  std::stoi(base_match[1].str(), nullptr);
    assert(base_match.length(2) > 0);
    Sumpt = std::stoi(base_match[3].str(), nullptr);
    L1SeedFun[SeedName] = std::bind(&L1AlgoFactory::Muer_ETM, this, Mupt, Sumpt);
    return true;
  }

  return false;
}       // -----  end of function L1Menu2016::ParseMuSum  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseComplexSingleMu
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseComplexSingleMu(const std::string& SeedName) 
{
  // This function overlap wiht ParseSingleObject. The code shouldn't come
  // here for simple singleMu trigger. Check for safety
  if (L1SeedFun.find(SeedName)!= L1SeedFun.end())
    return false;

  std::smatch base_match;
  std::regex integer("L1_Single(Mu)([0-9]+)([^0-9_]*)(_(EMTF|OMTF|BMTF))*(_Bx([-+0-9]+))*(_(Open|LowQ|HighQ))*");
  std::string L1object ="";
  std::string postfix = "";
  std::string muonType = "";
  std::string muonQual = "";
  int muonBx = 0;
  int imuonQual = 2;
  int imuonType = 0;
  float muonpt = -10;

  if (std::regex_match(SeedName, base_match, integer))
  {
    L1object = base_match[1];
    muonpt = std::stoi(base_match[2], nullptr);
    postfix =  base_match[3];
    if (base_match[4] != "")
    {
      muonType = base_match[5];
    }
    if (base_match[6]!="")
    {
      muonBx = std::stoi(base_match[7], nullptr);
    }
    if (base_match[8] != "")
    {
      muonQual = base_match[9];
    }
  } else return false;

  if (!muonQual.empty())
  {
    if (muonQual == "Open") imuonQual = 0;
    if (muonQual == "LowQ") imuonQual = 1;
    if (muonQual == "HighQ") imuonQual = 2;
  }

  if (!muonType.empty())
  {
    if (muonType == "BMTF") imuonType = 1;
    if (muonType == "OMTF") imuonType = 2;
    if (muonType == "EMTF") imuonType = 3;
  }
  assert(L1object == "Mu");
  bool isER = postfix=="er";
  L1SeedFun[SeedName] = std::bind(&L1AlgoFactory::ComplexSingleMu, this, muonpt,  isER, imuonQual, imuonType, muonBx);
  return true;
}       // -----  end of function L1Menu2016::ParseComplexSingleMu  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::FillLumiSection
//  Description:  
// ===========================================================================
bool L1Menu2016::FillLumiSection(int currentLumi)
{
  if (currentLumi == -1) return false;

  for(auto l1 : mL1Seed)
  {
    if(L1LSCount[l1.first].find(currentLumi) == L1LSCount[l1.first].end())
    {
      L1LSCount[l1.first][currentLumi] = 0;
    }
    if (l1.second.eventfire)
    {
      L1LSCount[l1.first][currentLumi]++;
    }
  }

  L1LSCount["Count"][currentLumi]++;
  
  return true;
}       // -----  end of function L1Menu2016::FillLumiSection  -----


// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::FillPileUpSec
//  Description:  
// ===========================================================================
bool L1Menu2016::FillPileUpSec()
{
  float pu = -1;
  //bool eFired = false;
  // Data
  if (event_->run > 1 && DataLSPU.find(event_->run) != DataLSPU.end())
  {
    if (DataLSPU[event_->run].find(event_->lumi) != DataLSPU[event_->run].end())
    {
      pu = DataLSPU[event_->run][event_->lumi];
    }
  }

  // MC
  if (event_->run == 1)
  {
    pu = event_->nPV;
  }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fill Rate per PU ~~~~~
  for(auto col : ColumnMap)
  {
    col.second->FillPileUpSec(pu);
  }

  return true;
}       // -----  end of function L1Menu2016::FillPileUpSec  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::PrintCSV
//  Description:  
// ===========================================================================
bool L1Menu2016::PrintCSV(std::ostream &out)
{
  if (!writecsv) return false;

  // L1Seeds
  std::vector<std::string> csvout;
  std::stringstream ss;
  ss << "L1Bit"         
    << "," << "L1SeedName"<< "," ;
  csvout.push_back(ss.str());
  for(auto sed : vL1Seed)
  {
    auto seed = mL1Seed[sed];
    ss.str("");
    ss << seed.bit
      << "," << sed<< "," ;
    csvout.push_back(ss.str());
  }

  // POG
  ss.str("");
  csvout.push_back(ss.str());

  int idx = 1000;
  for(auto pog : POGMap)
  {
    ss.str("");
    ss << idx++ << "," << pog.first <<",";
    csvout.push_back(ss.str());
  }

  // PAG
  ss.str("");
  csvout.push_back(ss.str());

  idx = 2000;
  for(auto pag : PAGMap)
  {
    ss.str("");
    ss << idx++ << "," << pag.first <<",";
    csvout.push_back(ss.str());
  }

  // Total
  csvout.push_back("");
  csvout.push_back("9999,Total rate,");

  for(auto col : ColumnMap)
  {
    col.second->PrintCSV(csvout, scale);
  }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Print out additional info ~~~~~
  std::vector<std::string>::iterator csvt = csvout.begin();
  ss.str("");
  ss <<"Type"<<","<<"PAG"<<","<<"Comment";
  csvt++->append(ss.str());
  // L1Seeds
  for(auto sed : vL1Seed)
  {
    auto seed = mL1Seed[sed];
    ss.str("");
    std::copy(seed.POG.begin(), seed.POG.end(), std::ostream_iterator<std::string>(ss, "|"));
    ss.seekp(ss.str().length()-1); //remove trailing |
    ss<<",";
    std::copy(seed.PAG.begin(), seed.PAG.end(), std::ostream_iterator<std::string>(ss, "|"));
    ss.seekp(ss.str().length()-1); //remove trailing |
    ss<<","<<seed.comment;
    csvt++->append(ss.str());
  }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Writing out to file ~~~~~
  for(auto i : csvout)
  {
    out << i <<std::endl;
  }
  return true;
}       // -----  end of function L1Menu2016::PrintCSV  -----


// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ReadDataPU
//  Description:  
// ===========================================================================
bool L1Menu2016::ReadDataPU() 
{
  const std::string pucsv = L1ConfigStr["Lumilist"];
  std::ifstream csvfile(pucsv);
  if (!csvfile)
  {
    std::cout << "Data PU CSV File "<<pucsv<<" is not found !"<<std::endl;
    return false;
  }

  std::string line;
  DataLSPU.clear();
  std::getline(csvfile, line); // Skip the first line;
  while (std::getline(csvfile, line))
  {
    std::istringstream iss(line);
    std::string seed;
    char c;
    int Fill, Run, LS;
    float pileup;
    iss >> Fill >> c >> Run >> c >> LS >> c >> pileup;
    DataLSPU[Run][LS] = pileup;
  }

  return true;
}       // -----  end of function L1Menu2016::ReadDataPU  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::CalLocalHT
//  Description:  
// ===========================================================================
void L1Menu2016::CalLocalHT(float &HTTcut, bool withHF)
{
  float sumJetHt= 0;
  for(UInt_t ue=0; ue < upgrade_->nJets; ue++) {
    Int_t bx = upgrade_->jetBx.at(ue);        		
    if(bx != 0) continue;
    Float_t pt = upgrade_->jetEt.at(ue);
    Float_t eta = upgrade_->jetEta.at(ue);
    if (pt >= L1Config["SumJetET"])
    {
      if (withHF)
      {
        sumJetHt += pt;
      }
      else{
        if ( fabs(eta) <= L1Config["SumJetEta"] )
          sumJetHt += pt;
      }
    }
  }
  HTTcut = sumJetHt;
  return;

}       // -----  end of function L1Menu2016::CalLocalHT  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::CalLocalHTM
//  Description:  
// ===========================================================================
void L1Menu2016::CalLocalHTM(float &HTMcut)
{
  
  TLorentzVector temp(0, 0,0,0);
  for(UInt_t ue=0; ue < upgrade_->nJets; ue++) {
    Int_t bx = upgrade_->jetBx.at(ue);        		
    if(bx != 0) continue;
    Float_t pt = upgrade_->jetEt.at(ue);
    Float_t eta = upgrade_->jetEta.at(ue);
    if (pt >= L1Config["SumJetET"] && fabs(eta) <= L1Config["SumJetEta"] )
    {
      TLorentzVector jet(0, 0,0,0);
      jet.SetPtEtaPhiE(upgrade_->jetEt.at(ue),
          upgrade_->jetEta.at(ue),
          upgrade_->jetPhi.at(ue), 
          0);
      temp -= jet;
    }
  }
  HTMcut = temp.Pt();
}       // -----  end of function L1Menu2016::CalLocalHTM  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ParseRanges
//  Description:  
// ===========================================================================
bool L1Menu2016::ParseRanges(const std::string Config,  std::vector<std::pair<unsigned int, unsigned int> >& Container)
{
  if (L1ConfigStr[Config] == "") return false;
  
  std::regex pattern("\\[\\s*([0-9]+,\\s*[0-9]+)\\s*\\]");
  std::regex bounds("([0-9]+),\\s*([0-9]+)");
  std::smatch base_match;

  for (std::sregex_token_iterator i(L1ConfigStr[Config].begin(), L1ConfigStr[Config].end(), pattern, 1); 
      i != std::sregex_token_iterator(); ++i) 
  {
	std::string match_str = i->str(); 
    if (std::regex_match(match_str, base_match, bounds))
    {
      unsigned lead =  std::stoi(base_match[1], nullptr);
      unsigned sec =  std::stoi(base_match[2], nullptr);
      assert( sec >= lead );
      Container.push_back(std::make_pair(lead, sec));
    }
  }

  return true;
}       // -----  end of function L1Menu2016::ParseRanges  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::CheckLS
//  Description:  Return whether to skip this LS
// ===========================================================================
bool L1Menu2016::CheckLS(unsigned int currentLumi) const
{
  for(auto p : pLS)
  {
    if (currentLumi >= p.first && currentLumi <= p.second)
    {
      if (L1Config.at("doScanLS"))
        std::cout << "LS"<<currentLumi << " within "<< L1ConfigStr.at("SelectLS") << " % "<< fChain->GetCurrentFile()->GetName()  << std::endl;
      return false;
    }
  }
  return true;
}       // -----  end of function L1Menu2016::CheckLS  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::CheckBX
//  Description:  Return whether to skip this BX
// ===========================================================================
bool L1Menu2016::CheckBX(unsigned int currentBX) const
{
  for(auto p : pBX)
  {
    if (currentBX >= p.first && currentBX <= p.second)
    {
      return false;
    }
  }
  return true;
}       // -----  end of function L1Menu2016::CheckBX  -----


// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::CalLocalETM
//  Description:  
// ===========================================================================
void L1Menu2016::CalLocalETM(float &ETMcut, bool withHF)
{

  if (!l1CaloTower_) return;

  TVector2 revec(0,0);
  
  // Ignore Tower 28
  int ietamax = 29;
  if (withHF)
    ietamax = 99;
  float metX =0;
  float metY =0;


  for(int jTower=0; jTower< l1CaloTower_ ->nTower; ++jTower){
    Int_t ieta = l1CaloTower_->ieta[jTower];
    Int_t iphi = l1CaloTower_->iphi[jTower];
    Int_t iet  = l1CaloTower_->iet[jTower];
    Double_t phi = (Double_t)iphi * TMath::Pi()/36.;
    Double_t et = 0.5 * (Double_t)iet;
    if (abs(ieta) == 28) continue;
    if(abs(ieta) < ietamax){
      metX -= et * TMath::Cos(phi);
      metY -= et * TMath::Sin(phi);
    }
  }

  revec.Set(metX, metY);
  ETMcut = revec.Mod();
}       // -----  end of function L1Menu2016::CalLocalETM  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::TokenGroups
//  Description:  
// ===========================================================================
std::vector<std::string> L1Menu2016::TokenGroups(std::string instring) const
{
  std::vector<std::string>  temp;
  if (instring.empty()) return temp;

  boost::char_separator<char> sep(",.;|- ");
  tokenizer tokens(instring, sep);
  for(auto &t : tokens)
  {
    temp.push_back(t);
  }
  return temp;
}       // -----  end of function L1Menu2016::TokenGroups  -----

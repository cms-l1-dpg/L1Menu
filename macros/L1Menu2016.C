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
  l1uGT(nullptr)
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
  WriteHistogram();
  outfile->close();
  outcsv->close();
  outrootfile->Close();
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
  L1Config["isData"] = 0;
  L1Config["SumJetET"] = 0;
  L1Config["nBunches"] = 3965.4; //Scaling Run259721 to 1.15E34 Lumi
  //L1Config["nBunches"] = 2736;
  //L1Config["nBunches"] = 0;
  L1Config["AveragePU"] = 0;
  L1Config["Energy"] = 0;
  L1Config["targetlumi"] = 0;
  L1Config["doPlotRate"] = 0;
  L1Config["doPlotEff"] = 0;
  L1Config["doPlotTest"] = 0;
  L1Config["doTnPMuon"] = 0;
  L1Config["doPrintLS"] = 0;
  L1Config["doPrintPU"] = 0;
  L1Config["doCompuGT"] = 0;
  L1Config["maxEvent"] = -1;
  L1Config["SetMuonER"] = -1;
  L1Config["SetNoPrescale"] = 0;
  L1Config["UseUpgradeLyr1"] = -1;
  L1Config["UseL1CaloTower"] = -1;
  L1Config["SelectRun"] = -1;
  L1Config["SelectEvent"] = -1;
  L1Config["UsePFMETNoMuon"] = 0;
  L1Config["UseuGTDecision"] = 0;
  L1Config["UseUnpackTree"] = 0;
  
  L1ObjectMap["Jet"] = &L1Event.JetPt;
  L1ObjectMap["JetC"] = &L1Event.JetCenPt;
  L1ObjectMap["Tau"] = &L1Event.TauPt;
  L1ObjectMap["Tauer"] = &L1Event.TauerPt;
  L1ObjectMap["IsoTau"] = &L1Event.IsoTauPt;
  L1ObjectMap["EG"] = &L1Event.EGPt;
  L1ObjectMap["EGer"] = &L1Event.EGerPt;
  L1ObjectMap["IsoEG"] = &L1Event.IsoEGPt;
  L1ObjectMap["IsoEGer"] = &L1Event.IsoEGerPt;
  L1ObjectMap["Mu"] = &L1Event.MuPt;
  L1ObjectMap["MuOpen"] = &L1Event.MuOpenPt;
  L1ObjectMap["Muer"] = &L1Event.MuerPt;
  L1ObjectMap["HTT"] = &L1Event.HTT;
  L1ObjectMap["HTM"] = &L1Event.HTM;
  L1ObjectMap["ETM"] = &L1Event.ETM;
  L1ObjectMap["ETT"] = &L1Event.ETT;


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Map to old func for now. ~~~~~
  // MutliJets
  L1SeedFun["L1_QuadJetC36_TauJet52"] = std::bind(&L1AlgoFactory::QuadJetCentral_TauJet, this, 36.,52.);
  L1SeedFun["L1_QuadJetC36_Tau52"] = std::bind(&L1AlgoFactory::QuadJetCentral_TauJet, this, 36.,52.);

  // MultiMuon
  L1SeedFun["L1_DoubleMu0er1p6_dEta_Max1p8_OS"] = std::bind(&L1AlgoFactory::Onia2015, this, 0.,0.,true,true,18);
  L1SeedFun["L1_DoubleMu0er1p6_dEta_Max1p8"] = std::bind(&L1AlgoFactory::Onia2015, this, 0.,0.,true,false,18);
  L1SeedFun["L1_DoubleMu0er1p25_dEta_Max1p8_OS"] = std::bind(&L1AlgoFactory::Onia2016, this, 0.,0.,true,true,18);
  L1SeedFun["L1_DoubleMu_10_0_dEta_Max1p8"] = std::bind(&L1AlgoFactory::Onia2015, this, 10.,0.,false,false,18);
  L1SeedFun["L1_DoubleMu0"] = std::bind(&L1AlgoFactory::DoubleMu, this, 0.,0.,true, false);
  L1SeedFun["L1_DoubleMuOpen"] = std::bind(&L1AlgoFactory::DoubleMuOpen, this, 0.);
  L1SeedFun["L1_DoubleMu_10_Open"] = std::bind(&L1AlgoFactory::DoubleMuXOpen, this, 10.);
  L1SeedFun["L1_DoubleMu_10_0"] = std::bind(&L1AlgoFactory::DoubleMu, this, 10.,0.,true, false);
  L1SeedFun["L1_DoubleMu_10_3p5"] = std::bind(&L1AlgoFactory::DoubleMu, this, 10.,3.5,true, false);
  L1SeedFun["L1_DoubleMu_12_5"] = std::bind(&L1AlgoFactory::DoubleMu, this, 12.,5.,true,false);
  L1SeedFun["L1_DoubleMu_11_4"] = std::bind(&L1AlgoFactory::DoubleMu, this, 11.,4.,true,false);
  L1SeedFun["L1_DoubleMu_12_8"] = std::bind(&L1AlgoFactory::DoubleMu, this, 12.,8.,true,false);
  L1SeedFun["L1_DoubleMu_13_6"] = std::bind(&L1AlgoFactory::DoubleMu, this, 13.,6.,true,false);
  L1SeedFun["L1_DoubleMu_15_5"] = std::bind(&L1AlgoFactory::DoubleMu, this, 15.,5.,true,false);
  L1SeedFun["L1_TripleMu0"] = std::bind(&L1AlgoFactory::TripleMu, this, 0.,0.,0.,1);
  L1SeedFun["L1_TripleMuOpen"] = std::bind(&L1AlgoFactory::TripleMu, this, 0.,0.,0.,0);
  L1SeedFun["L1_TripleMu_5_5_3"] = std::bind(&L1AlgoFactory::TripleMu, this, 5.,5.,3.,1);
  L1SeedFun["L1_QuadMu0"] = std::bind(&L1AlgoFactory::QuadMu, this, 0.,0.,0.,0.,1);

  //Cross
  L1SeedFun["L1_IsoEG20er_Tau20er_dEta_Min0p2"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 20.,20.);
  L1SeedFun["L1_IsoEG22er_Tau20er_dEta_Min0p2"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 22.,20.);
  L1SeedFun["L1_IsoEG20er_Tau24er_dEta_Min0p2"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 20.,24.);
  L1SeedFun["L1_IsoEG23er_TauJet20er_NotWdEta0"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 23.,20.); // l1t-tsg-v3:  L1_IsoEG20er_TauJet20er_NotWdEta0
  L1SeedFun["L1_IsoEG22er_TauJet20er_NotWdEta0"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 22.,20.); // l1t-tsg-v3:  L1_IsoEG20er_TauJet20er_NotWdEta0
  L1SeedFun["L1_IsoEG20er_Tau20er_NotWdEta0"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 20.,20.);
  L1SeedFun["L1_IsoEG22er_Tau20er_NotWdEta0"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 22.,20.);
  L1SeedFun["L1_IsoEG20er_Tau24er_NotWdEta0"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 20.,24.);
  L1SeedFun["L1_IsoEG23er_Tau20er_NotWdEta0"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 23.,20.); // l1t-tsg-v3:  L1_IsoEG20er_TauJet20er_NotWdEta0
  L1SeedFun["L1_DoubleMu6_EG6"] = std::bind(&L1AlgoFactory::DoubleMu_EG, this, 6.,6.,true);
  L1SeedFun["L1_DoubleMu6_EG16"] = std::bind(&L1AlgoFactory::DoubleMu_EG, this, 6.,16.,true); // l1t-tsg-v3:  L1_DoubleMu6_EG6
  L1SeedFun["L1_DoubleMu7_EG7"] = std::bind(&L1AlgoFactory::DoubleMu_EG, this, 7,7.,true);
  L1SeedFun["L1_DoubleMu7_EG14"] = std::bind(&L1AlgoFactory::DoubleMu_EG, this, 7,14.,true);  // l1t-tsg-v3:  L1_DoubleMu7_EG7
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
  L1SeedFun["L1_DoubleEG6_HTT150"] = std::bind(&L1AlgoFactory::DoubleEG_HT, this, 6., 150.);
  L1SeedFun["L1_DoubleEG6_HTT255"] = std::bind(&L1AlgoFactory::DoubleEG_HT, this, 6., 255.); // l1t-tsg-v3:  L1_DoubleEG6_HTT150
  L1SeedFun["L1_Jet32_DoubleMuOpen_Mu10_dPhi_Jet_Mu0_Max1p05_dPhi_Mu_Mu_Min1p0"] = std::bind(&L1AlgoFactory::Jet_MuOpen_Mu_dPhiMuMu1, this, 32.,10.);
  L1SeedFun["L1_Jet32_DoubleMuOpen_Mu10_dPhi_Jet_Mu0_Max0p4_dPhi_Mu_Mu_Min1p0"] = std::bind(&L1AlgoFactory::Jet_MuOpen_Mu_dPhiMuMu1, this, 32.,10.);
  L1SeedFun["L1_Jet32_MuOpen_EG10_dPhi_Jet_Mu_Max1p05_dPhi_Mu_EG_Min1p05"] = std::bind(&L1AlgoFactory::Jet_MuOpen_EG_dPhiMuEG1, this, 32.,10.);
  L1SeedFun["L1_Jet32_MuOpen_EG10_dPhi_Jet_Mu_Max0p4_dPhi_Mu_EG_Min1p0"] = std::bind(&L1AlgoFactory::Jet_MuOpen_EG_dPhiMuEG1, this, 32.,10.);
  L1SeedFun["L1_Jet32MuOpen_EG17_dPhiMu_EG1"] = std::bind(&L1AlgoFactory::Jet_MuOpen_EG_dPhiMuEG1, this, 32.,17.); // l1t-tsg-v3:  L1_Jet32MuOpen_EG10_dPhiMu_EG1

  L1SeedFun["L1_DoubleJet8_ForwardBackward"] = std::bind(&L1AlgoFactory::DoubleJet_ForwardBackward, this, 8., 8.); 
  L1SeedFun["L1_DoubleJet12_ForwardBackward"] = std::bind(&L1AlgoFactory::DoubleJet_ForwardBackward, this, 12., 12.); 
  L1SeedFun["L1_DoubleJet16_ForwardBackward"] = std::bind(&L1AlgoFactory::DoubleJet_ForwardBackward, this, 16., 16.); 
  L1SeedFun["L1_ETM60_Jet60_dPhi_Min0p4"] = std::bind(&L1AlgoFactory::ETM_Jet, this, 60., 60., false); 
  L1SeedFun["L1_HTM60_HTT260"] = std::bind(&L1AlgoFactory::HTM_HTT, this, 60., 260.); 
  L1SeedFun["L1_HTM80_HTT220"] = std::bind(&L1AlgoFactory::HTM_HTT, this, 80., 220.); 
  L1SeedFun["L1_Mu3_JetC35"] = std::bind(&L1AlgoFactory::Mu_Jet, this, 3., 35., false, true); 
  L1SeedFun["L1_Mu3_JetC16"] = std::bind(&L1AlgoFactory::Mu_Jet, this, 3., 16., false, true); 
  L1SeedFun["L1_Mu3_JetC60"] = std::bind(&L1AlgoFactory::Mu_Jet, this, 3., 60., false, true); 
  L1SeedFun["L1_Mu3_JetC120"] = std::bind(&L1AlgoFactory::Mu_Jet, this, 3., 120., false, true); 

  L1SeedFun["L1_ZeroBias"] = [](){return true;};
  return true;

}       // -----  end of function L1Menu2016::InitConfig  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::BookHistogram
//  Description:  
// ===========================================================================
bool L1Menu2016::BookHistogram()
{

  Int_t Nbin_max = 50;
  const Int_t NPAGS = 7;
  const Int_t NTRIGPHYS = 6;
  const Int_t N128 = 128;			// could be > 128 for "test seeds"

  HistMap["Cross"]           = new TH1F("h_Cross","h_Cross",Nbin_max,-0.5,(float)Nbin_max-0.5);
  HistMap["MultiCross"]      = new TH1F("h_MultiCross","h_MultiCross",Nbin_max,-0.5,(float)Nbin_max-0.5);
  HistMap["Sums"]            = new TH1F("h_Sums","h_Sums",Nbin_max,-0.5,(float)Nbin_max-0.5);
  HistMap["Jets"]            = new TH1F("h_Jets","h_Jets",Nbin_max,-0.5,(float)Nbin_max-0.5);
  HistMap["MultiJets"]       = new TH1F("h_MultiJets","h_MultiJets",Nbin_max,-0.5,(float)Nbin_max-0.5);
  HistMap["EG"]          = new TH1F("h_Egamma","h_Egamma",Nbin_max,-0.5,(float)Nbin_max-0.5);
  HistMap["MultiEG"]     = new TH1F("h_MultiEgamma","h_MultiEgamma",Nbin_max,-0.5,(float)Nbin_max-0.5);
  HistMap["Muons"]           = new TH1F("h_Muons","h_Muons",Nbin_max,-0.5,(float)Nbin_max-0.5);
  HistMap["MultiMuons"]      = new TH1F("h_MultiMuons","h_MultiMuons",Nbin_max,-0.5,(float)Nbin_max-0.5);
  HistMap["Technical"]       = new TH1F("h_Technical","h_Technical",Nbin_max,-0.5,(float)Nbin_max-0.5);

  HistMap["Block"]           = new TH1F("h_Block","h_Block",11,-0.5,10.5);

  HistMap["PAGS_pure"]       = new TH1F("h_PAGS_pure","h_PAGS_pure",NPAGS,-0.5,(float)NPAGS-0.5);
  HistMap["PAGS_shared"]     = new TH1F("h_PAGS_shared","h_PAGS_shared",NPAGS,-0.5,(float)NPAGS-0.5);

  HistMap["TRIGPHYS_pure"]   = new TH1F("h_TRIGPHYS_pure","h_TRIGPHYS_pure",NTRIGPHYS,-0.5,(float)NTRIGPHYS-0.5);
  HistMap["TRIGPHYS_shared"] = new TH1F("h_TRIGPHYS_shared","h_TRIGPHYS_shared",NTRIGPHYS,-0.5,(float)NTRIGPHYS-0.5);

  HistMap["All"]             = new TH1F("h_All","h_All",N128,-0.5,N128-0.5);
  HistMap["Pure"]           = new TH1F("h_Pure","h_Pure",N128,-0.5,N128-0.5);

  unsigned fbit = BitMap.begin()->first;
  unsigned lbit = BitMap.rbegin()->first;
  Hist2D["cor_Seeds"]       = new TH2F("cor_Seeds","cor_Seeds",lbit-fbit+1, fbit, lbit, lbit-fbit+1, fbit, lbit);
  for(auto i : BitMap)
  {
    int bin = Hist2D["cor_Seeds"]->GetXaxis()->FindBin(i.first);
    Hist2D["cor_Seeds"]->GetXaxis()->SetBinLabel(bin, i.second.c_str());
    Hist2D["cor_Seeds"]->GetYaxis()->SetBinLabel(bin, i.second.c_str());
  }
  Hist2D["cor_Block"]       = new TH2F("cor_Block","cor_Block",10,-0.5,9.5,19,-0.5,9.5);
  Hist2D["cor_PAGS"]        = new TH2F("cor_PAGS","cor_PAGS",NPAGS,-0.5,(float)NPAGS-0.5,NPAGS,-0.5,(float)NPAGS-0.5);
  Hist2D["cor_TRIGPHYS"]    = new TH2F("cor_TRIGPHYS","cor_TRIGPHYS",NTRIGPHYS,-0.5,(float)NTRIGPHYS-0.5,NTRIGPHYS,-0.5,(float)NTRIGPHYS-0.5);
  return true;
}       // -----  end of function L1Menu2016::BookHistogram  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::WriteHistogram
//  Description:  
// ===========================================================================
bool L1Menu2016::WriteHistogram() 
{

  outrootfile->cd();
  for(auto h : HistMap)
  {
    h.second->Write();
  }
  for(auto h : Hist2D)
  {
    h.second->Write();
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

  //Read the prescales table
  boost::char_separator<char> sep(",");
  std::ifstream menufile(menufilename);
  std::string line;
  if (!menufile)
  {
    std::cout << "MenuFile "<<menufilename<<" is not found !"<<std::endl;
    return false;
  }

  if (writefiles)
    *outfile <<  "---------------------------- Input Menu -------------------------" << std::endl;

  while (std::getline(menufile, line))
  {
    if (line.empty()) continue;
    if (writefiles)
      *outfile << line <<std::endl;
    if (line.at(0) == '#')
      continue;
    if (line.at(0) == '%')
    {
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
    temp.bit = bit;
    temp.comment = comline;
    if (prescale < 0)
      temp.prescale = INFTY;
    else
      temp.prescale = prescale;

    if (L1Config["doCompuGT"] || L1Config["SetNoPrescale"] )
      temp.prescale = 1;
      
    if (pog.length() != 0)
    {
      tokenizer tokens(pog, sep);
      for(auto &t : tokens)
      {
        temp.POG.push_back(t);
      }
    }
    if (pag.length() != 0)
    {
      tokenizer tokens2(pag, sep);
      for(auto &t : tokens2)
      {
        temp.PAG.push_back(t);
      }
    }
    if (writefiles)
    {
      assert(outfile != nullptr);
      
    }

    mL1Seed[seed] = temp;

  }
  
  if (writefiles)
    *outfile <<  "---------------------------- Input Menu -------------------------" <<std::endl << std::endl;

  return true;
}       // -----  end of function L1Menu2016::ReadMenu  -----

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
  if (iss.fail())
  {
    std::cout<<"\033[0;31mCan't parse config:\033[0m "<<line<< std::endl; 
  }
  
  if (L1Config.find(key) != L1Config.end())
  {
    L1Config[key] = value;
  }
  else
    std::cout << "Not reconfiguzed config key " << key<< std::endl;
  

  return true;
}       // -----  end of function L1Menu2016::ParseConfig  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::GetRunConfig
//  Description:  Get the running time config from command line
// ===========================================================================
bool L1Menu2016::GetRunConfig(std::map<std::string, float> &config)
{

  for(auto c : config)
  {
    if (L1Config.find(c.first) != L1Config.end())
    {
      L1Config[c.first] = c.second;
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
  std::cout << "Printed configuration ============================ " << std::endl;
  return true;
}       // -----  end of function L1Menu2016::PrintConfig  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::PreLoop
//  Description:  
// ===========================================================================
bool L1Menu2016::PreLoop(std::map<std::string, float> &config)
{
  GetRunConfig(config);

  ReadMenu();
  BuildRelation();
  L1SeedFunc();

  OpenWithList(tuplefilename);
  BookHistogram();

  // Set the Run Config as highest priority
  GetRunConfig(config);
  PrintConfig();

  if (writeplots)
  {
    l1Plot = new L1Plot(outrootfile, event_, upgrade_, recoJet_,
        recoSum_, recoEle_, recoMuon_, recoTau_, recoFilter_, l1CaloTower_, recoVtx_);
    l1Plot->SetTodo(L1Config);
    l1Plot->PreRun(&L1Event, &mL1Seed);
  }

  if (L1Config["doPrintPU"])
  {
    ReadDataPU();
  }
    
  if (L1Config["doTnPMuon"])
  {
    l1TnP = new L1TnP(outrootfile, event_, upgrade_, recoJet_,
        recoSum_, recoEle_, recoMuon_, recoTau_, recoFilter_, l1CaloTower_, recoVtx_);
    if (L1Config["doTnPMuon"])
      l1TnP->DoMuonTnP();
  }

  if (L1Config["doCompuGT"] || L1Config["UseuGTDecision"])
  {
    l1uGT = new L1uGT( outrootfile, event_, l1uGT_, &L1Event, &mL1Seed);
    l1uGT->GetTreeAlias(L1Ntuple::GetuGTAlias());
  }

  if (L1Config["SetMuonER"] != -1) SetMuonER(L1Config["SetMuonER"]);
  if (L1Config["UseUpgradeLyr1"] != -1) SetUseUpgradeLyr1(L1Config["UseUpgradeLyr1"]);
  if (L1Config["UseL1CaloTower"] != -1) SetUseL1CaloTower(L1Config["UseL1CaloTower"]);

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
  if (L1Config["SumJetET"] != 0)
  {
    CalLocalHT(L1Event.HTT);
    CalLocalHTM(L1Event.HTM);
  } else {
    L1AlgoFactory::HTTVal(L1Event.HTT);
    L1AlgoFactory::HTMVal(L1Event.HTM);
  }

  L1AlgoFactory::ETMVal(L1Event.ETM);
  L1AlgoFactory::ETTVal(L1Event.ETT);

  // Mulit
  float dummy = 0;
  L1AlgoFactory::DoubleMuPt(L1Event.doubleMuPt1,L1Event.doubleMuPt2);
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
  nFireevents = 0.;
  nLumi = 0;
  int i = -1;

  //for (Long64_t i=0; i<nevents; i++){     
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
        std::cout << "New Lumi section: " << event_->lumi << std::endl;      
        currentLumi=event_ -> lumi;
        nLumi++;
      } 
    } else if (i % 200000 == 0)
      std::cout << "Processed " << i << " events." << std::endl;

    nZeroBiasevents++;

    GetL1Event();

    if (RunMenu())
      nFireevents++;

    if (L1Config["doPrintLS"])
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

  return true;
}       // -----  end of function L1Menu2016::Loop  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::PostLoop
//  Description:  
// ===========================================================================
bool L1Menu2016::PostLoop()
{
  scale = CalScale(0, 0, true);

  std::cout << "Summary" << std::endl;

  for(auto &seed : mL1Seed)
  {
    seed.second.firerate = seed.second.firecounts *scale;
    seed.second.firerateerror = sqrt(seed.second.firecounts)*scale;
    seed.second.purerate = seed.second.purecounts *scale;
  }
  
  FillDefHist1D();
  FillDefHist2D();
  PrintRates(std::cout);
  if (writefiles)
    PrintRates(*outfile);
  PrintCSV(*outcsv);

  if (l1Plot != NULL)
    l1Plot->PostRun(scale);

  if (l1TnP != NULL)
    l1TnP->PostRun();

  if (L1Config["doPrintPU"])
  {
    PrintPUCSV();
  }
  return true;
}       // -----  end of function L1Menu2016::PostLoop  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::PrintPUCSV
//  Description:  
// ===========================================================================
bool L1Menu2016::PrintPUCSV()
{
  const int nBunches = 2736;
  std::fstream pucsv (outputdir + "/" + outputname+"_PU" +".csv", std::fstream::out );
  pucsv <<"L1Seed,PileUp,Fired,Total,Rate,Error"<<std::endl;
  for(auto l1seed : L1PUCount)
  {
    if (l1seed.first == "Count") continue;
    if (l1seed.first.find("L1A_") !=std::string::npos) continue;

    for(auto pu : l1seed.second)
    {
      pucsv << l1seed.first <<","<< pu.first <<","
        <<pu.second<<","
        <<L1PUCount["Count"][pu.first]<<","
        <<pu.second * CalScale(L1PUCount["Count"][pu.first], nBunches) <<","
        <<sqrt(pu.second) * CalScale(L1PUCount["Count"][pu.first], nBunches)
        <<std::endl;
    }
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
  
  float totalrate = 0.;
  float totalpurerate = 0.;
  bool bybit = true;
  std::size_t L1NameLength = 0;
  for(auto k : mL1Seed)
  {
    L1NameLength = k.first.size() > L1NameLength ? k.first.size() : L1NameLength;
  }

  out << std::left
      << std::setw(10)             << "L1Bit"
      << std::setw(L1NameLength+2) << "L1SeedName"
      << std::setw(10)             << "pre-scale"
      << std::setw(10)             << "rate@13TeV"       << " +/- "
      << std::setw(20)             << "error_rate@13TeV"
      << std::setw(10)             << "pure@13TeV"       
      << "Comments"
      << std::endl;

  if (bybit)
  {
    for(auto i : BitMap)
    {
      auto seed = mL1Seed[i.second];
      out << std::left
          << std::setw(10)             << seed.bit
          << std::setw(L1NameLength+2) << seed.name
          << std::setw(10)             << seed.prescale
          << std::setw(10)             << seed.firerate      << " +/- "
          << std::setw(20)             << seed.firerateerror
          << std::setw(10)             << seed.purerate      
          << seed.comment
          << std::endl;
      totalrate +=seed.firerate;
      totalpurerate +=seed.purerate;
    }
    
  }
  else{
    for(auto seed : mL1Seed)
    {
      out << std::left
          << std::setw(10)             << seed.second.bit
          << std::setw(L1NameLength+2) << seed.first
          << std::setw(10)             << seed.second.prescale
          << std::setw(10)             << seed.second.firerate      << " +/- "
          << std::setw(20)             << seed.second.firerateerror
          << std::setw(10)             << seed.second.purerate      
          << seed.second.comment
          << std::endl;
      totalrate +=seed.second.firerate;
      totalpurerate +=seed.second.purerate;
    }

  }


  out << std::endl << "Total rate  = " << nFireevents / 1000 * scale 
    <<" +/- " << sqrt(nFireevents) * scale / 1000 << " (kHz)" << std::endl;
  out << std::endl << "Total rate (with overlaps) = " << totalrate / 1000 << " (kHz)" << std::endl;
  out << std::endl << "Total pure rate  = " << totalpurerate / 1000 <<" (kHz)" << std::endl;

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

  // Total
  for(auto& pog : POGMap)
  {
    PhyCounts[pog.first] = 0;
    PhyPureCounts[pog.first] = 0;
  }
  for(auto& pag : PAGMap)
  {
    PhyCounts[pag.first] = 0;
    PhyPureCounts[pag.first] = 0;
  }

  return true;
}       // -----  end of function L1Menu2016::BuildRelation  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::L1SeedFunc
//  Description:  
// ===========================================================================
bool L1Menu2016::L1SeedFunc()
{
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


void L1Menu2016::CorrectScale(TH1* h, Float_t scal) {

  Int_t nbins = h -> GetNbinsX();

  for (Int_t i=1; i<= nbins; i++)  {
    Float_t val = h -> GetBinContent(i);
    Float_t er = sqrt(val);
    val = val * scal;
    er = er * scal;
    h -> SetBinContent(i,val);
    h -> SetBinError(i,er);
  }
}

bool L1Menu2016::InsertInMenu(std::string L1name, bool value) {


  bool post_prescale = false;

  if ( mL1Seed.find(L1name) == mL1Seed.end() ) {
    std::cout << "This shouldn't happen!" << std::endl;
    return false;
  }

  mL1Seed[L1name].eventfire = false;
  mL1Seed[L1name].ncounts++;
  if ( mL1Seed[L1name].ncounts % mL1Seed[L1name].prescale == 0) 
    post_prescale = value; 

  mL1Seed[L1name].eventfire = post_prescale;
  if (post_prescale)
  {
    mL1Seed[L1name].firecounts++;
    FireSeed.insert(L1name);
  }

  return true;
}

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::BindAlgo
//  Description:  
// ===========================================================================
bool L1Menu2016::BindAlgo() 
{

  return true;
}       // -----  end of function L1Menu2016::BindAlgo  -----

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
  FireSeed.clear();
  FiredPhy.clear();
  for(auto& seed: mL1Seed)
  {
    if (L1Config["UseuGTDecision"])
    {
      assert(l1uGT != NULL);
      InsertInMenu(seed.first, l1uGT->GetuGTDecision(seed.first));
    }
    else
      InsertInMenu(seed.first, CheckL1Seed(seed.first));
  }

  CheckPhysFire();
  CheckPureFire();

  return FireSeed.size() > 0;
}       // -----  end of function L1Menu2016::RunMenu  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::CheckPureFire
//  Description:  
// ===========================================================================
bool L1Menu2016::CheckPureFire() 
{
  // Trigger path pure rate
  if (FireSeed.size() == 1) 
    mL1Seed[*(FireSeed.begin())].purecounts++;

  // POG pure rate
  std::set<std::string> POGset;
  for(auto fireit : FireSeed)
  {
    for(auto &pog : mL1Seed[fireit].POG)
    {
      POGset.insert(pog);
    }
  }
  if (POGset.size() == 1)
  {
    PhyPureCounts[*(POGset.begin())]++;
  }

  std::set<std::string> PAGset;
  for(auto fireit : FireSeed)
  {
    for(auto &pag : mL1Seed[fireit].PAG)
    {
      PAGset.insert(pag);
    }
  }
  if (PAGset.size() == 1)
  {
    PhyPureCounts[*(PAGset.begin())]++;
  }

  return true;
}       // -----  end of function L1Menu2016::CheckPureFire  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::CheckPhysFire
//  Description:  
// ===========================================================================
bool L1Menu2016::CheckPhysFire()
{
  for(auto fired : FireSeed)
  {
    L1Seed &seed = mL1Seed[fired];
    //if (writefiles)
      //*outfile <<  event_->run <<","<<event_->lumi<<"," <<event_->event<<","<<seed.name << std::endl;
    for(auto pog : seed.POG)
    {
      if (FiredPhy.insert(pog).second) PhyCounts[pog]++;
    }
    for(auto pag : seed.PAG)
    {
      if (FiredPhy.insert(pag).second) PhyCounts[pag]++;
    }
  }

  return true;
}       // -----  end of function L1Menu2016::CheckPhysFire  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::CalScale
//  Description:  
// ===========================================================================
double L1Menu2016::CalScale(int nEvents_, int nBunches_, bool print) 
{
  double scale = 0.0;
  int nEvents = nEvents_ == 0 ? nZeroBiasevents : nEvents_;
  double nBunches = nBunches_ == 0 ?  L1Config["nBunches"] : nBunches_;

  if (L1Config["nBunches"] == -1)
  {
    //scal = (80.*631.)/(1326*23.3);      
    scale = (80.*631.)/(nLumi*23.3);      
    if (print)
      std::cout << "Scale by "   << "(80.*631.)/(nLumi*23.3) with nLumi = " << nLumi      << std::endl;
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
//         Name:  L1Menu2016::FillDefHist1D
//  Description:  Fill in the redantant Histogram, for comparison with the old
//  code
// ===========================================================================
bool L1Menu2016::FillDefHist1D()
{
  for(auto pog : POGMap)
  {
    if (HistMap.find(pog.first) == HistMap.end())
      continue;
    //Sort by bit
    std::sort( pog.second.begin(), pog.second.end());
    int binidx =1;
    for(auto l1bit : pog.second)
    {
      std::string l1name = BitMap[l1bit];
      std::cout << l1bit << " " << l1name<< std::endl;
      HistMap[pog.first]->GetXaxis()->SetBinLabel(binidx, l1name.c_str());
      HistMap[pog.first]->SetBinContent(binidx, mL1Seed[l1name].firecounts);
      binidx++;
    }
    //assert(binidx == pog.second.size());
    CorrectScale(HistMap[pog.first], scale);
  }

  return true;
}       // -----  end of function L1Menu2016::FillDefHist1D  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::FillDefHist2D
//  Description:  
// ===========================================================================
bool L1Menu2016::FillDefHist2D()
{
  for(auto h2d : Hist2D)
  {
    CorrectScale(h2d.second, scale);
  }
  return true;
}       // -----  end of function L1Menu2016::FillDefHist2D  -----
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
  std::regex integerobj("Single([^0-9]+)([0-9]+)([^0-9]*)");
  std::regex integerSum("(ETM|HTT|ETT|HTM)([0-9]+)");
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
    postfix = "";
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
  std::regex integer("L1_DoubleJet([C]*)([0-9]+)");
  if (std::regex_match(SeedName, base_match, integer))
  {
    bool isCentral = base_match.length(1) == 1;
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
  std::regex integer("L1_Double(Iso|)Tau([0-9]+)er");
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
  std::regex integer_sys("L1_QuadJet([C]*)([0-9]+)");
  std::regex integer_asys("L1_QuadJet([C]*)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)");
  if (std::regex_match(SeedName, base_match, integer_sys))
  {
    isCentral = base_match.length(1) == 1;
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
  std::regex integerEGerHTT("L1_(Iso|)EG([0-9]+)er_HTT([0-9]+)");
  if (std::regex_match(SeedName, base_match, integerEGerHTT))
  {
    isIsoEG = base_match.length(1) == 3;
    EGpt =  std::stoi(base_match[2].str(), nullptr);
    Sumpt = std::stoi(base_match[3].str(), nullptr);
    L1SeedFun[SeedName] = std::bind(&L1AlgoFactory::SingleEG_Eta2p1_HTT, this, EGpt, Sumpt, isIsoEG);
    return true;
  }

  return false;
}       // -----  end of function L1Menu2016::ParseEGSum  -----

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
  std::regex integer("L1_Mu([0-9]+)er_(Iso|)Tau(Jet|)([0-9]+)er");
  if (std::regex_match(SeedName, base_match, integer))
  {
    Mupt =  std::stoi(base_match[1].str(), nullptr);
    isIsoTau = base_match.length(2) == 3;
    Taupt = std::stoi(base_match[4].str(), nullptr);
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
  std::regex integerMuHTT("L1_Mu([0-9]+)_HTT([0-9]+)");
  std::regex integerMuerETM("L1_Mu([0-9]+)er_ETM([0-9]+)");
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
    Sumpt = std::stoi(base_match[2].str(), nullptr);
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
  std::regex integer("L1_Single([^0-9_]+)([0-9]+)([^0-9_]*)(_(EMTF|OMTF|BMTF))*(_Bx([-+0-9]+))*(_(Open|LowQ|HighQ))*");
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
  bool eFired = false;
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
    //pu = event_->event
  }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ L1Seed ~~~~~
  for(auto l1 : mL1Seed)
  {
    if(L1PUCount[l1.first].find(pu) == L1PUCount[l1.first].end())
    {
      L1PUCount[l1.first][pu] = 0;
    }
    if (l1.second.eventfire)
    {
      eFired= true;
      L1PUCount[l1.first][pu]++;
    }
  }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ POG ~~~~~
  std::set<std::string> POGset;
  for(auto& pog : POGMap)
  {
    std::string l1pog = "L1T_"+pog.first;
    std::string l1pogpure = "L1T_Pure"+pog.first;
    if(L1PUCount[l1pog].find(pu) == L1PUCount[l1pog].end())
    {
      L1PUCount[l1pog][pu] = 0;
    }
    if(L1PUCount[l1pogpure].find(pu) == L1PUCount[l1pogpure].end())
    {
      L1PUCount[l1pogpure][pu] = 0;
    }

    if (FiredPhy.find(pog.first) != FiredPhy.end())
    {
      L1PUCount[l1pog][pu] ++;
      POGset.insert(pog.first);
    }
  }
  if (POGset.size() == 1)
  {
    std::string l1pogpure = "L1T_Pure"+*(POGset.begin());
    L1PUCount[l1pogpure][pu]++;
  }
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PAG ~~~~~
  std::set<std::string> PAGset;
  for(auto& pag : PAGMap)
  {
    std::string l1pag = "L1A_"+pag.first;
    std::string l1pagpure = "L1A_Pure"+pag.first;
    if(L1PUCount[l1pag].find(pu) == L1PUCount[l1pag].end())
    {
      L1PUCount[l1pag][pu] = 0;
    }
    if(L1PUCount[l1pagpure].find(pu) == L1PUCount[l1pagpure].end())
    {
      L1PUCount[l1pagpure][pu] = 0;
    }
    if (FiredPhy.find(pag.first) != FiredPhy.end())
    {
      L1PUCount[l1pag][pu] ++;
      PAGset.insert(pag.first);
    }
  }
  if (PAGset.size() == 1)
  {
    std::string l1pagpure = "L1A_Pure"+*(PAGset.begin());
    L1PUCount[l1pagpure][pu]++;
  }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Total ~~~~~
  L1PUCount["Count"][pu]++;
  if (eFired)
  {
    L1PUCount["L1APhysics"][pu]++;
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
  out << "L1Bit"
      << "," << "L1SeedName"
      << "," << "pre-scale"
      << "," << "rate"
      << "," << "error_rate"
      << "," << "pure"       
      << "," << "Comments"
      << std::endl;


  bool bybit = true;
  float totalrate = 0.;
  float totalpurerate = 0.;

  if (bybit)
  {
    for(auto i : BitMap)
    {
      auto seed = mL1Seed[i.second];
      out << seed.bit
        << "," << seed.name
        << "," << seed.prescale
        << "," << seed.firerate     
        << "," << seed.firerateerror
        << "," << seed.purerate      
        << ",\""<<seed.comment<<"\""
        << std::endl;
      totalrate +=seed.firerate;
      totalpurerate +=seed.purerate;
    }
  } else{
    for(auto seed : mL1Seed)
    {
      out << seed.second.bit
          << "," << seed.first
          << "," << seed.second.prescale
          << "," << seed.second.firerate
          << "," << seed.second.firerateerror
          << "," << seed.second.purerate      
          << ",\""<<seed.second.comment<<"\""
          << std::endl;
      totalrate +=seed.second.firerate;
      totalpurerate +=seed.second.purerate;
    }
  }
  
  out << std::endl;

  out << "L1Type"
      << "," << "rate(kHz)"
      << "," << "error_rate(kHz)"
      << "," << "pure(kHz)"       
      << "," << "error_rate(kHz)"
      << std::endl;

  for(auto pog : POGMap)
  {
    out << pog.first 
      <<","<< PhyCounts[pog.first]           *scale / 1000.
      <<","<< sqrt(PhyCounts[pog.first])     *scale / 1000.
      <<","<< PhyPureCounts[pog.first]       *scale / 1000.
      <<","<< sqrt(PhyPureCounts[pog.first]) *scale / 1000.
      << std::endl;
  }
  out << std::endl;

  out << std::endl;
  out << std::endl << "Total rate  = " << nFireevents / 1000 * scale 
    <<" +/- " << sqrt(nFireevents) * scale / 1000 << " (kHz)" << std::endl;
  out << std::endl;
  out << std::endl << "Total rate (with overlaps) = " << totalrate / 1000 << " (kHz)" << std::endl;
  out << std::endl;
  out << std::endl << "Total pure rate  = " << totalpurerate / 1000 <<" (kHz)" << std::endl;
  return true;
}       // -----  end of function L1Menu2016::PrintCSV  -----


// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ReadDataPU
//  Description:  
// ===========================================================================
bool L1Menu2016::ReadDataPU() 
{
  const std::string pucsv = "menu/run_lumi.csv";
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
void L1Menu2016::CalLocalHT(float &HTTcut)
{
  float sumJetHt= 0;
  for(UInt_t ue=0; ue < upgrade_->nJets; ue++) {
    Int_t bx = upgrade_->jetBx.at(ue);        		
    if(bx != 0) continue;
    Float_t pt = upgrade_->jetEt.at(ue);
    if (pt >= L1Config["SumJetET"])
      sumJetHt += pt;
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
    if (pt >= L1Config["SumJetET"])
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

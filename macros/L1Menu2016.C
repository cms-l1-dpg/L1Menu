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
  scale(0)
{
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
    outfile = new std::fstream( outputdir + "/" + outputname +".txt",
        std::fstream::in | std::fstream::out );
  if (writecsv)
    outcsv = new std::fstream( outputdir + "/" + outputname +".csv",
        std::fstream::in | std::fstream::out );
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
  L1Config["NumberOfBunches"] = 0;
  L1Config["AveragePU"] = 0;
  L1Config["Energy"] = 0;
  L1Config["targetlumi"] = 0;
  
  L1ObjectMap["Jet"] = &L1Event.JetPt;
  L1ObjectMap["JetC"] = &L1Event.JetCenPt;
  L1ObjectMap["Tau"] = &L1Event.TauPt;
  L1ObjectMap["EG"] = &L1Event.EGPt;
  L1ObjectMap["EGer"] = &L1Event.EGerPt;
  L1ObjectMap["IsoEG"] = &L1Event.IsoEGPt;
  L1ObjectMap["IsoEGer"] = &L1Event.IsoEGerPt;
  L1ObjectMap["Mu"] = &L1Event.MuPt;
  L1ObjectMap["MuOpen"] = &L1Event.MuPt;
  L1ObjectMap["Muer"] = &L1Event.MuerPt;
  L1ObjectMap["HTT"] = &L1Event.HTT;
  L1ObjectMap["ETM"] = &L1Event.ETM;
  L1ObjectMap["ETT"] = &L1Event.ETT;


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Map to old func for now. ~~~~~
  // MutliJets
  //L1SeedFun["L1_DoubleJetC52"] = std::bind(&L1AlgoFactory::DoubleJet, this, 52, 52, true);
  //L1SeedFun["L1_DoubleJetC84"] = std::bind(&L1AlgoFactory::DoubleJet, this, 84, 84, true);
  //L1SeedFun["L1_DoubleJetC100"] = std::bind(&L1AlgoFactory::DoubleJet, this, 100, 100, true);
  //L1SeedFun["L1_DoubleJetC112"] = std::bind(&L1AlgoFactory::DoubleJet, this, 112, 112, true);
  //L1SeedFun["L1_DoubleIsoTau28er"] = std::bind(&L1AlgoFactory::DoubleTauJetEta2p17, this, 28.,28.,true);
  //L1SeedFun["L1_DoubleIsoTau32er"] = std::bind(&L1AlgoFactory::DoubleTauJetEta2p17, this, 32.,32.,true);
  //L1SeedFun["L1_DoubleIsoTau36er"] = std::bind(&L1AlgoFactory::DoubleTauJetEta2p17, this, 36.,36.,true);
  //L1SeedFun["L1_DoubleIsoTau40er"] = std::bind(&L1AlgoFactory::DoubleTauJetEta2p17, this, 40.,40.,true);
  L1SeedFun["L1_DoubleTauJet40er"] = std::bind(&L1AlgoFactory::DoubleTauJetEta2p17, this, 40.,40.,false);
  //L1SeedFun["L1_TripleJet_92_76_64_VBF"] = std::bind(&L1AlgoFactory::TripleJet_VBF, this, 92.,76.,64.,1);
  //L1SeedFun["L1_TripleJet_84_68_48_VBF"] = std::bind(&L1AlgoFactory::TripleJet_VBF, this, 84.,68.,48.,1);
  //L1SeedFun["L1_QuadJetC40"] = std::bind(&L1AlgoFactory::QuadJet, this, 40.,40.,40.,40.,true);
  //L1SeedFun["L1_QuadJetC60"] = std::bind(&L1AlgoFactory::QuadJet, this, 60.,60.,60.,60.,true);
  L1SeedFun["L1_QuadJetC36_TauJet52"] = std::bind(&L1AlgoFactory::QuadJetCentral_TauJet, this, 36.,52.);

  // MultiMuon
  L1SeedFun["L1_DoubleMu0_Eta1p6_WdEta18_OS"] = std::bind(&L1AlgoFactory::Onia2015, this, 0.,0.,true,true,18);
  L1SeedFun["L1_DoubleMu0_Eta1p6_WdEta18"] = std::bind(&L1AlgoFactory::Onia2015, this, 0.,0.,true,false,18);
  L1SeedFun["L1_DoubleMu_10_0_WdEta18"] = std::bind(&L1AlgoFactory::Onia2015, this, 10.,0.,false,false,18);
  L1SeedFun["L1_DoubleMu0"] = std::bind(&L1AlgoFactory::DoubleMu, this, 0.,0.,true, false);
  L1SeedFun["L1_DoubleMuOpen"] = std::bind(&L1AlgoFactory::DoubleMuXOpen, this, 0.);
  L1SeedFun["L1_DoubleMu_10_Open"] = std::bind(&L1AlgoFactory::DoubleMuXOpen, this, 10.);
  L1SeedFun["L1_DoubleMu_10_0"] = std::bind(&L1AlgoFactory::DoubleMu, this, 10.,0.,true, false);
  L1SeedFun["L1_DoubleMu_10_3p5"] = std::bind(&L1AlgoFactory::DoubleMu, this, 10.,3.5,true, false);
  L1SeedFun["L1_DoubleMu_12_5"] = std::bind(&L1AlgoFactory::DoubleMu, this, 12.,5.,true,false);
  L1SeedFun["L1_TripleMu0"] = std::bind(&L1AlgoFactory::TripleMu, this, 0.,0.,0.,4);
  L1SeedFun["L1_TripleMu_5_5_3"] = std::bind(&L1AlgoFactory::TripleMu, this, 5.,5.,3.,4);
  L1SeedFun["L1_QuadMu0"] = std::bind(&L1AlgoFactory::QuadMu, this, 0.,0.,0.,0.,4);

  //Cross
  L1SeedFun["L1_Mu6_HTT100"] = std::bind(&L1AlgoFactory::Mu_HTT, this, 6.,100.);
  L1SeedFun["L1_Mu8_HTT50"] = std::bind(&L1AlgoFactory::Mu_HTT, this, 8.,50.);
  L1SeedFun["L1_Mu0er_ETM40"] = std::bind(&L1AlgoFactory::Muer_ETM, this, 0.,40.);
  L1SeedFun["L1_Mu0er_ETM55"] = std::bind(&L1AlgoFactory::Muer_ETM, this, 0.,55.);
  L1SeedFun["L1_Mu10er_ETM30"] = std::bind(&L1AlgoFactory::Muer_ETM, this, 10.,30.);
  L1SeedFun["L1_Mu10er_ETM50"] = std::bind(&L1AlgoFactory::Muer_ETM, this, 10.,50.);
  L1SeedFun["L1_Mu14er_ETM30"] = std::bind(&L1AlgoFactory::Muer_ETM, this, 14.,30.);
  L1SeedFun["L1_EG25er_HTT100"] = std::bind(&L1AlgoFactory::SingleEG_Eta2p1_HTT, this, 25., 100.,false);
  L1SeedFun["L1_Mu16er_TauJet20er"] = std::bind(&L1AlgoFactory::Muer_TauJetEta2p17, this, 16.,20.,false);
  L1SeedFun["L1_Mu16er_IsoTau28er"] = std::bind(&L1AlgoFactory::Muer_TauJetEta2p17, this, 16.,28.,true);
  L1SeedFun["L1_Mu16er_IsoTau32er"] = std::bind(&L1AlgoFactory::Muer_TauJetEta2p17, this, 16.,32.,true);
  L1SeedFun["L1_Mu12_EG10"] = std::bind(&L1AlgoFactory::Mu_EG, this, 12.,10.,false, 4);
  L1SeedFun["L1_Mu20_EG10"] = std::bind(&L1AlgoFactory::Mu_EG, this, 20.,10.,false, 4);
  L1SeedFun["L1_Mu4_EG18"] = std::bind(&L1AlgoFactory::Mu_EG, this, 4.,18.,false, 4);
  L1SeedFun["L1_Mu5_EG15"] = std::bind(&L1AlgoFactory::Mu_EG, this, 5.,15.,false, 4);
  L1SeedFun["L1_Mu5_EG20"] = std::bind(&L1AlgoFactory::Mu_EG, this, 5.,20.,false, 4);
  L1SeedFun["L1_Mu5_IsoEG18"] = std::bind(&L1AlgoFactory::Mu_EG, this, 5.,18.,true, 4);
  L1SeedFun["L1_IsoEG20er_TauJet20er_NotWdEta0"] = std::bind(&L1AlgoFactory::IsoEGer_TauJetEta2p17, this, 20.,20.);
  L1SeedFun["L1_DoubleMu6_EG6"] = std::bind(&L1AlgoFactory::DoubleMu_EG, this, 6.,6.,true);
  L1SeedFun["L1_DoubleMu7_EG7"] = std::bind(&L1AlgoFactory::DoubleMu_EG, this, 7,7.,true);
  L1SeedFun["L1_Mu5_DoubleEG5"] = std::bind(&L1AlgoFactory::Mu_DoubleEG, this, 5., 5.);
  L1SeedFun["L1_Mu6_DoubleEG10"] = std::bind(&L1AlgoFactory::Mu_DoubleEG, this, 6., 10.);

  //MultiCross
  L1SeedFun["L1_Mu3_JetC16_WdEtaPhi2"] = std::bind(&L1AlgoFactory::Mu_JetCentral_delta, this, 3.,16.);
  L1SeedFun["L1_Mu3_JetC52_WdEtaPhi2"] = std::bind(&L1AlgoFactory::Mu_JetCentral_delta, this, 3.,52.);
  L1SeedFun["L1_DoubleJetC56_ETM60"] = std::bind(&L1AlgoFactory::DoubleJetCentral_ETM, this, 56.,56.,60.);
  L1SeedFun["L1_DoubleEG6_HTT150"] = std::bind(&L1AlgoFactory::DoubleEG_HT, this, 6., 150.);
  L1SeedFun["L1_Jet32MuOpen_Mu10_dPhiMu_Mu1"] = std::bind(&L1AlgoFactory::Jet_MuOpen_Mu_dPhiMuMu1, this, 32.,10.);
  L1SeedFun["L1_Jet32MuOpen_EG10_dPhiMu_EG1"] = std::bind(&L1AlgoFactory::Jet_MuOpen_EG_dPhiMuEG1, this, 32.,10.);

  //MultiEG
  //L1SeedFun["L1_DoubleEG_15_10"] = std::bind(&L1AlgoFactory::DoubleEG, this, 15.,10.,false );
  //L1SeedFun["L1_DoubleEG_22_10"] = std::bind(&L1AlgoFactory::DoubleEG, this, 22.,10.,false );
  //L1SeedFun["L1_TripleEG_14_10_8"] = std::bind(&L1AlgoFactory::TripleEG, this, 14.,10.,8. );

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
    std::cout << "MenuFile "<<menufile<<" is not found !"<<std::endl;
    return false;
  }

  while (std::getline(menufile, line))
  {
    if (line.empty()) continue;
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
    //std::cout << temp.name <<" " << temp.bit <<" "<< temp.prescale << " " 
      //<< temp.POG.size() <<" " << temp.PAG.size() << std::endl;

    mL1Seed[seed] = temp;

  }

  return true;
}       // -----  end of function L1Menu2016::ReadMenu  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ReadFilelist
//  Description:  
// ===========================================================================
bool L1Menu2016::OpenWithList(std::string filelist)
{
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
  
  if (L1Config.find(key) != L1Config.end())
  {
    L1Config[key] = value;
  }
  else
    std::cout << "Not reconfiguzed config key " << key<< std::endl;
  

  return true;
}       // -----  end of function L1Menu2016::ParseConfig  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::PrintConfig
//  Description:  
// ===========================================================================
bool L1Menu2016::PrintConfig() const
{
  std::cout << "Printing configuration ___________________________ " << std::endl;
  for(auto &x: L1Config)
    std::cout << x.first <<" : " << x.second << std::endl;
  std::cout << std::endl;
  return true;
}       // -----  end of function L1Menu2016::PrintConfig  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::PreLoop
//  Description:  
// ===========================================================================
bool L1Menu2016::PreLoop()
{
  InitConfig();
  ReadMenu();
  BuildRelation();
  L1SeedFunc();

  OpenWithList(tuplefilename);
  
  BookHistogram();

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
  L1AlgoFactory::SingleTauPt(L1Event.TauPt, false);

  //Mu
  L1AlgoFactory::SingleMuPt(L1Event.MuPt, false);
  L1AlgoFactory::SingleMuPt(L1Event.MuerPt, true);

  //Sum
  L1AlgoFactory::HTTVal(L1Event.HTT);
  L1AlgoFactory::ETMVal(L1Event.ETM);
  L1AlgoFactory::ETTVal(L1Event.ETT);

  return true;
}       // -----  end of function L1Menu2016::GetL1Event  -----
// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::Loop
//  Description:  
// ===========================================================================
bool L1Menu2016::Loop()
{
  ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Initial ~~~~~
  Int_t nevents = fChain->GetEntriesFast();//GetEntries();
  int currentLumi(-1);
  nZeroBiasevents = 0.;
  nLumi = 0;

  for (Long64_t i=0; i<nevents; i++){     
    Long64_t ientry = LoadTree(i); 
    if (ientry < 0) break;
    GetEntry(i);

    if (event_ != NULL )
    {
      if(event_ -> lumi != currentLumi){
        std::cout << "New Lumi section: " << event_->lumi << std::endl;      
        currentLumi=event_ -> lumi;
        nLumi++;
      } 
    } else if (i % 20000 == 0)
        std::cout << "Processed " << i << " events." << std::endl;

    nZeroBiasevents++;

    GetL1Event();
    RunMenu();
    FillLumiSection(currentLumi);

  }

  return true;
}       // -----  end of function L1Menu2016::Loop  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::PostLoop
//  Description:  
// ===========================================================================
bool L1Menu2016::PostLoop()
{
  CalScale();

  std::cout << "Summary" << std::endl;

  for(auto &seed : mL1Seed)
  {
    seed.second.firerate = seed.second.firecounts *scale;
    seed.second.firerateerror = sqrt(seed.second.firecounts)*scale;
    seed.second.purerate = seed.second.purecounts *scale;
  }
  
  FillDefHist1D();
  PrintRates(std::cout);
  if (writefiles)
    PrintRates(*outfile);
  PrintCSV(*outcsv);

  return true;
}       // -----  end of function L1Menu2016::PostLoop  -----

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


  out << std::endl << "Total rate (without overlaps) = " << totalrate << std::endl;
  out << std::endl << "Total pure rate  = " << totalpurerate << std::endl;

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
  PhyCounts["All"] = 0;
  

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


void L1Menu2016::CorrectScale(TH1F* h, Float_t scal) {

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

Bool_t L1Menu2016::EGamma() {

  //insert_ibin = 0;

  //Float_t SingleEGPtCut = -10.;
  //SingleEGPt(SingleEGPtCut,false, false);

  //Float_t SingleIsoEGerPtCut = -10.;
  //SingleEGPt(SingleIsoEGerPtCut,true,true);


  //InsertInMenu("L1_SingleEG2_BptxAND",SingleEGPtCut >= 2.);
  //InsertInMenu("L1_SingleEG5",SingleEGPtCut >= 5.);
  //InsertInMenu("L1_SingleEG10",SingleEGPtCut >= 10.);
  //InsertInMenu("L1_SingleEG15",SingleEGPtCut >= 15.);
  //InsertInMenu("L1_SingleEG20",SingleEGPtCut >= 20.);
  //InsertInMenu("L1_SingleEG25",SingleEGPtCut >= 25.);
  //InsertInMenu("L1_SingleEG30",SingleEGPtCut >= 30.);
  //InsertInMenu("L1_SingleEG35",SingleEGPtCut >= 35.);
  //InsertInMenu("L1_SingleEG40",SingleEGPtCut >= 40.);
  //InsertInMenu("L1_SingleIsoEG20", SingleEG(20.,true,false) );
  //InsertInMenu("L1_SingleIsoEG25", SingleEG(25.,true,false) );
  //InsertInMenu("L1_SingleIsoEG18er",SingleIsoEGerPtCut >= 18.);
  //InsertInMenu("L1_SingleIsoEG20er",SingleIsoEGerPtCut >= 20.);
  //InsertInMenu("L1_SingleIsoEG22er",SingleIsoEGerPtCut >= 22.);
  //InsertInMenu("L1_SingleIsoEG25er",SingleIsoEGerPtCut >= 25.);
  //InsertInMenu("L1_SingleIsoEG30er",SingleIsoEGerPtCut >= 30.);

  //Int_t NN = insert_ibin;

  //Int_t kOFFSET_old = kOFFSET;
  //for (Int_t k=0; k < NN; k++) {
    //TheTriggerBits[k + kOFFSET_old] = insert_val[k];
  //}
  //kOFFSET += insert_ibin;

  //if (first) {

    //NBITS_EGAMMA = NN;

    //for (Int_t ibin=0; ibin < insert_ibin; ibin++) {
      //TString l1name = (TString)insert_names[ibin];
      //h_Egamma -> GetXaxis() -> SetBinLabel(ibin+1, l1name );
    //}
    //h_Egamma-> GetXaxis() -> SetBinLabel(NN+1,"EGAMMA");

    //for (Int_t k=1; k <= kOFFSET -kOFFSET_old; k++) {
      //h_All -> GetXaxis() -> SetBinLabel(k +kOFFSET_old , h_Egamma -> GetXaxis() -> GetBinLabel(k) );
    //}
  //}                      

  //Bool_t res = false;      
  //for (Int_t i=0; i < NN; i++) {
    //res = res || insert_val[i] ;
    //if (insert_val[i]) h_Egamma -> Fill(i);
  //}      
  //if (res) h_Egamma -> Fill(NN);

  //return res;
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
  assert(L1Seed.find("L1_") != std::string::npos);
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
  for(auto& seed: mL1Seed)
  {
    //std::cout <<  seed.first <<" " << CheckL1Seed(seed.first) << std::endl;
    InsertInMenu(seed.first, CheckL1Seed(seed.first));
    //std::cout <<  seed.first <<" " << ParseSingleObject(seed.first) << std::endl;
  }

  CheckPhysFire();
  CheckPureFire();

  return true;
}       // -----  end of function L1Menu2016::RunMenu  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::CheckPureFire
//  Description:  
// ===========================================================================
bool L1Menu2016::CheckPureFire() 
{
  if (FireSeed.size() != 1) return false;
  std::set<std::string>::const_iterator fireit = FireSeed.begin();
  mL1Seed[*fireit].purecounts++;
  for(auto &pog : mL1Seed[*fireit].POG)
  {
    PhyPureCounts[pog]++;
  }
  for(auto &pag : mL1Seed[*fireit].PAG)
  {
    PhyPureCounts[pag]++;
  }

  return true;
}       // -----  end of function L1Menu2016::CheckPureFire  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::CheckPhysFire
//  Description:  
// ===========================================================================
bool L1Menu2016::CheckPhysFire()
{
  for(auto &phy : PhyCounts)
  {
    if (phy.first == "All") continue;
    bool phyfire = false;
    if (POGMap.find(phy.first) != POGMap.end())
    {
      for(auto &pogseed : POGMap[phy.first])
      {
        phyfire =  phyfire || mL1Seed[BitMap[pogseed]].eventfire;
      }
    }
    if (PAGMap.find(phy.first) != PAGMap.end())
    {
      for(auto &pagseed : PAGMap[phy.first])
        phyfire =  phyfire || mL1Seed[BitMap[pagseed]].eventfire;
    }

    if (phyfire)
    {
      phy.second++;
    }
  }

  if (FireSeed.size() > 0) PhyCounts["All"]++;

  return true;
}       // -----  end of function L1Menu2016::CheckPhysFire  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::CalScale
//  Description:  
// ===========================================================================
bool L1Menu2016::CalScale() 
{
  if (L1Config["NumberOfBunches"] == -1)
  {
    //scal = (80.*631.)/(1326*23.3);      
    scale = (80.*631.)/(nLumi*23.3);      
  } else {
    scale = 11246.; // ZB per bunch in Hz
    //scale /= nZeroBiasevents*1000.; // in kHz
    scale /= nZeroBiasevents; // in Hz
    scale *= L1Config["NumberOfBunches"];
  }
  return true;
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

    //for (int i = 0; i < HistMap[pog.first]->GetNbinsX(); ++i)
    //{
      //std::cout << HistMap[pog.first]->GetXaxis()->GetBinLabel(i) <<" " 
        //<< HistMap[pog.first]->GetBinContent(i) << std::endl;
    //}
  }

  return true;
}       // -----  end of function L1Menu2016::FillDefHist1D  -----

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

  // EGMass
  if (ParseMultiEGMass(SeedName)) return true;

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
  std::regex integerSum("(ETM|HTT|ETT)([0-9]+)");
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
  const int jetclass = 1; 
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
  
  return true;
}       // -----  end of function L1Menu2016::FillLumiSection  -----

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
  }
  return true;
}       // -----  end of function L1Menu2016::PrintCSV  -----

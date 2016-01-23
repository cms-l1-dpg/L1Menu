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
L1Menu2016::L1Menu2016 ():
  writefiles(true), drawplots(true)
{
  InitConfig();
}  // -----  end of method L1Menu2016::L1Menu2016  (constructor)  -----

//----------------------------------------------------------------------------
//       Class:  L1Menu2016
//      Method:  L1Menu2016
// Description:  copy constructor
//----------------------------------------------------------------------------
L1Menu2016::L1Menu2016 ( const L1Menu2016 &other )
{
}  // -----  end of method L1Menu2016::L1Menu2016  (copy constructor)  -----

//----------------------------------------------------------------------------
//       Class:  L1Menu2016
//      Method:  ~L1Menu2016
// Description:  destructor
//----------------------------------------------------------------------------
L1Menu2016::~L1Menu2016 ()
{
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

  L1SeedFun["L1_TripleJet_92_76_64_VBF"] = std::bind(&L1AlgoFactory::TripleJet_VBF, this, 92.,76.,64.,1);
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
  HistMap["Egamma"]          = new TH1F("h_Egamma","h_Egamma",Nbin_max,-0.5,(float)Nbin_max-0.5);
  HistMap["MultiEgamma"]     = new TH1F("h_MultiEgamma","h_MultiEgamma",Nbin_max,-0.5,(float)Nbin_max-0.5);
  HistMap["Muons"]           = new TH1F("h_Muons","h_Muons",Nbin_max,-0.5,(float)Nbin_max-0.5);
  HistMap["MultiMuons"]      = new TH1F("h_MultiMuons","h_MultiMuons",Nbin_max,-0.5,(float)Nbin_max-0.5);
  HistMap["Technical"]       = new TH1F("h_Technical","h_Technical",Nbin_max,-0.5,(float)Nbin_max-0.5);

  HistMap["Block"]           = new TH1F("h_Block","h_Block",11,-0.5,10.5);
  //HistMap["cor_Block"]       = new TH2F("cor_Block","cor_Block",10,-0.5,9.5,19,-0.5,9.5);

  //HistMap["cor_PAGS"]        = new TH2F("cor_PAGS","cor_PAGS",NPAGS,-0.5,(float)NPAGS-0.5,NPAGS,-0.5,(float)NPAGS-0.5);
  HistMap["PAGS_pure"]       = new TH1F("h_PAGS_pure","h_PAGS_pure",NPAGS,-0.5,(float)NPAGS-0.5);
  HistMap["PAGS_shared"]     = new TH1F("h_PAGS_shared","h_PAGS_shared",NPAGS,-0.5,(float)NPAGS-0.5);

  //HistMap["cor_TRIGPHYS"]    = new TH2F("cor_TRIGPHYS","cor_TRIGPHYS",NTRIGPHYS,-0.5,(float)NTRIGPHYS-0.5,NTRIGPHYS,-0.5,(float)NTRIGPHYS-0.5);
  HistMap["TRIGPHYS_pure"]   = new TH1F("h_TRIGPHYS_pure","h_TRIGPHYS_pure",NTRIGPHYS,-0.5,(float)NTRIGPHYS-0.5);
  HistMap["TRIGPHYS_shared"] = new TH1F("h_TRIGPHYS_shared","h_TRIGPHYS_shared",NTRIGPHYS,-0.5,(float)NTRIGPHYS-0.5);

  HistMap["All"]             = new TH1F("h_All","h_All",N128,-0.5,N128-0.5);
  HistMap["Pure"]           = new TH1F("h_Pure","h_Pure",N128,-0.5,N128-0.5);
  return true;
}       // -----  end of function L1Menu2016::BookHistogram  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::ReadMenu
//  Description:  /* cursor */
// ===========================================================================
bool L1Menu2016::ReadMenu(std::string themenufilename)
{
  mL1Seed.clear();

  //Read the prescales table
  boost::char_separator<char> sep(",");
  std::ifstream menufile(themenufilename);
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

    std::istringstream iss(line);

    std::string seed;
    int bit;
    int prescale;
    std::string pog, pag;

    iss >> seed >> bit >> prescale >> pog >> pag;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Form L1Seed ~~~~~
    L1Seed temp;
    temp.name = seed;
    temp.bit = bit;
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
    std::cout << temp.name <<" " << temp.bit <<" "<< temp.prescale << " " << temp.POG.size() <<" " << temp.PAG.size() << std::endl;

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

  std::cout.flush();std::cout<<"Going to init the available trees..."<<std::endl;std::cout.flush();
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
    std::cout << "Not reconfiguzed config: " << key<< std::endl;
  

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
  

  return true;
}       // -----  end of function L1Menu2016::PreLoop  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::GetL1Event
//  Description:  
// ===========================================================================
bool L1Menu2016::GetL1Event()
{
  L1Event = {};

  SingleEGPt(L1Event.EGPt,false, false);
  SingleEGPt(L1Event.EGerPt,false, true);
  SingleEGPt(L1Event.IsoEGPt,true, false);
  SingleEGPt(L1Event.IsoEGerPt,true, true);



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
  //Double_t nZeroBiasevents = 0.;
  int nLumi(0),currentLumi(-1);
  PreLoop();

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


    GetL1Event();
    RunMenu();

    //if (i > 10) break;

  }


    ////FilL1Bits();
    //if(first) MyInit();

    //Bool_t raw = true; //= PhysicsBits[0];  // check for ZeroBias triggered events
    //if(!raw) continue;

    //nZeroBiasevents++;

    //// reset the emulated trigger bits
    //kOFFSET = 0;
    //for (Int_t k=0; k < N128; k++) {
      //TheTriggerBits[k] = false;
    //}

    //Bool_t jets        = Jets() ;
    //Bool_t multijets   = MultiJets() ;
    //Bool_t eg          = EGamma();
    //Bool_t multieg     = MultiEGamma();
    //Bool_t muons       = Muons();
    //Bool_t multimuons  = MultiMuons();
    //Bool_t sums        = Sums();
    //Bool_t technical   = Technical();
    //Bool_t cross       = Cross();
    //Bool_t multicross  = MultiCross();

    //Bool_t pass  = jets || multijets || eg || multieg || sums || muons || multimuons || cross || multicross || technical;

    //if(pass) NPASS ++;

    //if(cross)      NCROSS++;
    //if(multicross) MULTINCROSS++;
    //if(muons)      NMUONS++;
    //if(multimuons) MULTINMUONS++;
    //if(sums)       NSUMS++;
    //if(eg)         NEG++;
    //if(multieg)    MULTINEG++;
    //if(jets)       NJETS++;
    //if(multijets)  MULTINJETS++;
    //if(technical)  TECHNICAL++;

    //if(pass) h_Block->Fill(10.);

    //Bool_t dec[10];
    //dec[0] = eg;
    //dec[1] = multieg;
    //dec[2] = jets;
    //dec[3] = multijets;
    //dec[4] = muons;
    //dec[5] = multimuons;
    //dec[6] = sums;
    //dec[7] = cross;
    //dec[8] = multicross;
    //dec[9] = technical;

    //for (Int_t l=0; l < 9; l++) {
      //if(dec[l]){
	//h_Block -> Fill(l);
	//for (Int_t k=0; k < 5; k++) {
	  //if (dec[k]) cor_Block -> Fill(l,k);
	//}
      //}

    //}

    //first = false;

    //// now the pure rate stuff
    //// kOFFSET now contains the number of triggers we have calculated

    //Bool_t ddd[NPAGS];
    //for (Int_t idd=0; idd < NPAGS; idd++) {
      //ddd[idd] = false; 
    //} 

    //Bool_t eee[NTRIGPHYS];
    //for (Int_t iee=0; iee < NTRIGPHYS; iee++) {
      //eee[iee] = false; 
    //}

    //Float_t weightEventPAGs = 1.;
    //Float_t weightEventTRIGPHYS = 1.;

    //for (Int_t k=0; k < kOFFSET; k++) {
      //if ( ! TheTriggerBits[k] ) continue;
      //h_All -> Fill(k);

      //TString name = h_All -> GetXaxis() -> GetBinLabel(k+1);
      //std::string L1namest = (std::string)name;

      //Bool_t IsTOP   = setTOP.count(L1namest) > 0;
      //Bool_t IsHIGGS = setHIGGS.count(L1namest) > 0;
      //Bool_t IsBPH   = setBPH.count(L1namest) > 0;
      //Bool_t IsEXO   = setEXO.count(L1namest) > 0;
      //Bool_t IsSUSY  = setSUSY.count(L1namest) > 0;
      //Bool_t IsSMP   = setSMP.count(L1namest) > 0;
      //Bool_t IsB2G   = setB2G.count(L1namest) > 0;
      //if(IsHIGGS) ddd[0] = true;
      //if(IsSUSY)  ddd[1] = true;
      //if(IsEXO)   ddd[2] = true;
      //if(IsTOP)   ddd[3] = true;
      //if(IsSMP)   ddd[4] = true;
      //if(IsBPH)   ddd[5] = true;
      //if(IsB2G)   ddd[6] = true;

      //Float_t ww = WeightsPAGs[L1namest];
      //if (ww < weightEventPAGs) weightEventPAGs = ww;

      //Bool_t IsMuon     = setMuon.count(L1namest) > 0;
      //Bool_t IsEG       = setEG.count(L1namest) > 0;
      //Bool_t IsHadronic = setHadronic.count(L1namest) > 0;

      //Bool_t IsMuonEG       = setMuonEG.count(L1namest) > 0;
      //Bool_t IsMuonHadronic = setMuonHadronic.count(L1namest) > 0;
      //Bool_t IsEGHadronic   = setEGHadronic.count(L1namest) > 0;

      //if(IsMuon)     eee[0] = true;
      //if(IsEG)       eee[1] = true;
      //if(IsHadronic) eee[2] = true;

      //if(IsMuonEG)       eee[3] = true;
      //if(IsMuonHadronic) eee[4] = true;
      //if(IsEGHadronic)   eee[5] = true;

      //Float_t www = WeightsTRIGPHYS[L1namest];
      //if(www < weightEventTRIGPHYS) weightEventTRIGPHYS = www;

      ////did the event pass another trigger ?
      //Bool_t pure = true;
      //for (Int_t k2=0; k2 < kOFFSET; k2++) {
	//if (k2 == k) continue;
	//if ( TheTriggerBits[k2] ) pure = false;
      //}
      //if (pure) h_Pure -> Fill(k);
    //}

    //// for the PAG rates
    //Bool_t PAG = false;
    //for (Int_t idd=0; idd < NPAGS; idd++) {
      //if (ddd[idd]) {
	//Bool_t pure = true;
	//PAG = true;
	//for (Int_t jdd=0; jdd < NPAGS; jdd++) {
	  //if (ddd[jdd]) {
		//cor_PAGS -> Fill(idd,jdd);
		//if (jdd != idd) pure = false;
	  //}
	//}   
	//if(pure) h_PAGS_pure -> Fill(idd);
	//h_PAGS_shared -> Fill(idd,weightEventPAGs);

      //}  
    //}
    //if(PAG) nPAG ++;

    ////for the TRIGPHYS rates :
    //Bool_t TRIGPHYS = false;
    //for (Int_t iee=0; iee < NTRIGPHYS; iee++) {
      //if (eee[iee]) {
	//Bool_t pure = true;
	//TRIGPHYS = true;
	//for (Int_t jee=0; jee < NTRIGPHYS; jee++) {
	  //if (eee[jee]) {
		//cor_TRIGPHYS -> Fill(iee,jee);
		//if (jee != iee) pure = false;
	  //}
	//}   
	//if(pure) h_TRIGPHYS_pure -> Fill(iee);
	//h_TRIGPHYS_shared -> Fill(iee,weightEventTRIGPHYS);

      //}  
    //}
    //if(TRIGPHYS) nTRIGPHYS++;

  //}  // end evt loop

  PostLoop();

  return true;
}       // -----  end of function L1Menu2016::Loop  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::PostLoop
//  Description:  
// ===========================================================================
bool L1Menu2016::PostLoop()
{
  

  return true;
}       // -----  end of function L1Menu2016::PostLoop  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::BuildRelation
//  Description:  
// ===========================================================================
bool L1Menu2016::BuildRelation()
{
  

  return true;
}       // -----  end of function L1Menu2016::BuildRelation  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::L1SeedFunc
//  Description:  
// ===========================================================================
bool L1Menu2016::L1SeedFunc() const
{
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

  //bool post_prescale = false;

  if ( mL1Seed.find(L1name) == mL1Seed.end() ) {
    std::cout << "This shouldn't happen!" << std::endl;
    return false;
  }

  //if ( mL1Seed[L1name].ncounts % mL1Seed[L1name].prescale == 0) 
    //post_prescale = value; 

  //void(post_prescale);
  //insert_names[insert_ibin] = L1name;
  //insert_val[insert_ibin] = post_prescale ;
  //insert_ibin++;

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
//         Name:  L1Menu2016::ParseSingleObject
//  Description:  /* cursor */
// ===========================================================================
bool L1Menu2016::ParseSingleObject(const std::string SeedName)
{
  if (SeedName.find("L1_") == std::string::npos)
  {
    return false;
  }
  if (SeedName.find("Single") == std::string::npos)
  {
    return false;
  }

  bool returnbool = true;
  boost::char_separator<char> sep("_");
  tokenizer tokens(SeedName, sep);
  if (std::distance(tokens.begin(), tokens.end()) < 2) return false;
  boost::tokenizer<boost::char_separator<char> >::iterator tokenit = tokens.begin();
  tokenit++;
  std::string Seedtoken(*tokenit);

  std::string L1object ="";
  std::string postfix = "";
  unsigned int pt = -10;

  std::smatch base_match;
  std::regex integer("Single([^0-9]+)([0-9]+)([^0-9]*)");
  if (std::regex_match(Seedtoken, base_match, integer))
  {
	// The first sub_match is the whole string; the next
	// sub_match is the first parenthesized expression.
    if (base_match.size() >= 3) {
      std::ssub_match base_sub_match = base_match[1];
      L1object = base_sub_match.str();
      base_sub_match = base_match[2];
      pt = std::stoi(base_sub_match.str(), nullptr);
      if (base_match.size() > 3)
      {
        base_sub_match = base_match[3];
        postfix = base_sub_match.str();
      }
    }
  }

  L1object += postfix;
  if (L1ObjectMap.find(L1object) != L1ObjectMap.end())
  {
    //std::cout << L1object <<" "<<*(L1ObjectMap[L1object])<<" "<< pt<<" ";
    returnbool = returnbool && *(L1ObjectMap[L1object])> pt;
  }

  std::cout <<  std::distance(tokenit, tokens.end())<< std::endl;
  if (std::distance(tokenit, tokens.end()) > 1) 
  {
    tokenit++;
    Seedtoken = *tokenit;
    if (Seedtoken == "NotBptxOR" || Seedtoken == "BptxAND")
    {
      returnbool = returnbool && true;
    }
  }
  return returnbool;
}       // -----  end of function L1Menu2016::ParseSingleObject  -----

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
  if (L1Seed.find("Single") != std::string::npos)
  {
    return ParseSingleObject(L1Seed);
  }
  std::cout << "No function call for :" << L1Seed << std::endl;
  return false;
}       // -----  end of function L1Menu2016::CheckL1Seed  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Menu2016::RunMenu
//  Description:  
// ===========================================================================
bool L1Menu2016::RunMenu()
{
  for(auto& seed: mL1Seed)
  {
    std::cout <<  seed.first <<" " << CheckL1Seed(seed.first) << std::endl;
    InsertInMenu(seed.first, CheckL1Seed(seed.first));
    //std::cout <<  seed.first <<" " << ParseSingleObject(seed.first) << std::endl;
  }

  return true;
}       // -----  end of function L1Menu2016::RunMenu  -----

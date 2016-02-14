#include "L1Ntuple.h"

Long64_t L1Ntuple::GetEntries()
{
  return nentries_;
}

L1Ntuple::L1Ntuple()
{
  doEvent       = true;
  doreco        = true;
  domuonreco    = true;
  dol1extra     = true;
  dol1emuextra  = true;
  dol1menu      = true;
  doRecoJet     = true;
  doRecoMet     = true;
  doRecoEle     = true;
  doRecoMuon    = true;
  doRecoTau     = true;

  fChain        = nullptr;
  fEvent        = nullptr;
  ftreeEvent    = nullptr;
  ftreemuon     = nullptr;
  ftreereco     = nullptr;
  ftreeExtra    = nullptr;
  ftreeMenu     = nullptr;
  ftreeEmuExtra = nullptr;
  ftreeRecoJet  = nullptr;
  ftreeRecoMet  = nullptr;
  ftreeRecoEle  = nullptr;
  ftreeRecoMuon = nullptr;
  ftreeRecoTau  = nullptr;

  event_        = nullptr;
  upgrade_      = nullptr;
  recoJet_       = nullptr;
  recoSum_       = nullptr;
  recoEle_       = nullptr;
  recoMuon_      = nullptr;
  recoTau_       = nullptr;
  MainTreePath = "l1UpgradeTree/L1UpgradeTree";
}

L1Ntuple::L1Ntuple(const std::string & fname)
{
  L1Ntuple();
  Open(fname);
}

bool L1Ntuple::OpenWithList(const std::string & fname)
{
  if (!OpenNtupleList(fname)) exit(0);
  if (!CheckFirstFile())      exit(0);
  if (!OpenWithoutInit())     exit(0);

  std::cout.flush();std::cout<<"Going to init the available trees..."<<std::endl;std::cout.flush();
  Init();

  return true;
}

bool L1Ntuple::Open(const std::string & fname)
{
  listNtuples.push_back(fname);

  if (!CheckFirstFile())  exit(0);
  if (!OpenWithoutInit()) exit(0);

  std::cout.flush();std::cout<<"Going to init the available trees..."<<std::endl;std::cout.flush();
  Init();

  return true;
}

bool L1Ntuple::OpenNtupleList(const std::string & fname)
{
  std::ifstream flist(fname.c_str());
  if (!flist)
  {
    std::cout << "File "<<fname<<" is not found !"<<std::endl;
    return false;
  }

  while(!flist.eof())
  {
    std::string str;
    getline(flist,str);
    if (!flist.fail())
    {
      if (str!="") listNtuples.push_back(str);
    }
  }

  return true;
}

bool L1Ntuple::CheckFirstFile()
{
  if (listNtuples.size()==0) return false;

  rf = TFile::Open(listNtuples[0].c_str());

  if (rf==0) return false;
  if (rf->IsOpen()==0) return false;

  TTree * myChain        = (TTree*) rf->Get(MainTreePath.c_str());
  if (myChain == NULL) // in case other name
  {
    MainTreePath = "l1UpgradeEmuTree/L1UpgradeTree";
    myChain = (TTree*) rf->Get(MainTreePath.c_str());
  }

  TTree * mytreeEvent    = (TTree*) rf->Get("l1EventTree/L1EventTree");
  TTree * mytreemuon     = (TTree*) rf->Get("l1MuonRecoTreeProducer/MuonRecoTree");
  TTree * mytreejets     = (TTree*) rf->Get("l1RecoTreeProducer/RecoTree");
  TTree * mytreeExtra    = (TTree*) rf->Get("l1ExtraTreeProducer/L1ExtraTree");
  TTree * mytreeEmuExtra = (TTree*) rf->Get("l1EmulatorExtraTree/L1ExtraTree");
  TTree * mytreeMenu     = (TTree*) rf->Get("l1MenuTreeProducer/L1MenuTree");
  TTree * mytreeRecoJet  = (TTree*) rf->Get("l1JetRecoTree/JetRecoTree");
  TTree * mytreeRecoMet  = nullptr;
  //TTree * mytreeRecoMet  = (TTree*) rf->Get("l1MetFilterRecoTree/MetFilterRecoTree");
  TTree * mytreeRecoEle  = (TTree*) rf->Get("l1ElectronRecoTree/ElectronRecoTree");
  TTree * mytreeRecoMuon = (TTree*) rf->Get("l1MuonRecoTree/Muon2RecoTree");
  TTree * mytreeRecoTau  = (TTree*) rf->Get("l1TauRecoTree/TauRecoTree");

  if (!myChain) {
    std::cout<<"L1Tree not found .... "<<std::endl;
    return false;
  } else {
    std::cout<<"Main tree is found .." <<std::endl;
  }

  if (!mytreeEvent) {
    std::cout<<"EventTree not found, it will be skipped..."<<std::endl;
    doEvent=false;
  } else
  {
    std::cout << "EventTree is found ..."<<std::endl;
  }

  if (!mytreeRecoJet) {
    std::cout<<"RecoJet not found, it will be skipped..."<<std::endl;
    doRecoJet=false;
  } else
  {
    std::cout << "RecoJet is found ..."<<std::endl;
  }

  if (!mytreeRecoMet) {
    std::cout<<"RecoMet not found, it will be skipped..."<<std::endl;
    doRecoMet=false;
  } else
  {
    std::cout << "RecoMet is found ..."<<std::endl;
  }
  doRecoMet=false;

  if (!mytreeRecoEle) {
    std::cout<<"RecoEle not found, it will be skipped..."<<std::endl;
    doRecoEle=false;
  } else
  {
    std::cout << "RecoEle is found ..."<<std::endl;
  }

  if (!mytreeRecoMuon) {
    std::cout<<"RecoMuon not found, it will be skipped..."<<std::endl;
    doRecoMuon=false;
  } else
  {
    std::cout << "RecoMuon is found ..."<<std::endl;
  }

  if (!mytreeRecoTau) {
    std::cout<<"RecoTau not found, it will be skipped..."<<std::endl;
    doRecoTau=false;
  } else
  {
    std::cout << "RecoTau is found ..."<<std::endl;
  }

  if (!mytreejets) {
    std::cout<<"RecoTree not found, it will be skipped..."<<std::endl;
    doreco=false;
  } else
  {
    std::cout << "RecoTree is found ..."<<std::endl;
  }


  if (!mytreemuon) {
    std::cout<<"MuonRecoTree not found, it will be skipped..."<<std::endl;
    domuonreco=false;
  }
  else
  {
    std::cout << "MuonRecoTree is found ..."<<std::endl;
  }

  if (!mytreeExtra) {
    std::cout<<"L1ExtraTree not found, it will be skipped..."<<std::endl;
    dol1extra=false;
  }
  else
  {
    std::cout << "L1ExtraTree is found ..."<<std::endl;
  }
  if(!mytreeEmuExtra){
    std::cout<<"L1EmuExtraTree not found, it will be skipped..."<<std::endl;
    dol1emuextra=false;
  }
  else{
    std::cout << "L1EmuExtraTree is found ..."<<std::endl;
  }

  if (!mytreeMenu) {
    std::cout<<"L1MenuTree not found, it will be skipped..."<<std::endl;
    dol1menu=false;
  }
  else
  {
    std::cout << "L1MenuTree is found ..."<<std::endl;
  }

  return true;
}


bool L1Ntuple::OpenWithoutInit()
{
  fChain        = new TChain(MainTreePath.c_str());
  ftreeEvent    = new TChain("l1EventTree/L1EventTree");
  ftreemuon     = new TChain("l1MuonRecoTreeProducer/MuonRecoTree");
  ftreereco     = new TChain("l1RecoTreeProducer/RecoTree");
  ftreeExtra    = new TChain("l1ExtraTreeProducer/L1ExtraTree");
  ftreeEmuExtra = new TChain("l1EmulatorExtraTree/L1ExtraTree");
  ftreeMenu     = new TChain("l1MenuTreeProducer/L1MenuTree");
  ftreeRecoJet  = new TChain("l1JetRecoTree/JetRecoTree");
  //ftreeRecoMet  = new TChain("l1MetFilterRecoTree/MetFilterRecoTree");
  ftreeRecoEle  = new TChain("l1ElectronRecoTree/ElectronRecoTree");
  ftreeRecoMuon = new TChain("l1MuonRecoTree/Muon2RecoTree");
  ftreeRecoTau  = new TChain("l1TauRecoTree/TauRecoTree");


  for (unsigned int i=0;i<listNtuples.size();i++)
  {
    std::cout << " -- Adding " << listNtuples[i] << std::endl;
    fChain->Add(listNtuples[i].c_str());

    if (doEvent)      ftreeEvent    -> Add(listNtuples[i].c_str());
    if (doreco)       ftreereco     -> Add(listNtuples[i].c_str());
    if (domuonreco)   ftreemuon     -> Add(listNtuples[i].c_str());
    if (dol1extra)    ftreeExtra    -> Add(listNtuples[i].c_str());
    if (dol1emuextra) ftreeEmuExtra -> Add(listNtuples[i].c_str());
    if (dol1menu)     ftreeMenu     -> Add(listNtuples[i].c_str());
    if (doRecoJet)    ftreeRecoJet  -> Add(listNtuples[i].c_str());
    if (doRecoMet)    ftreeRecoMet  -> Add(listNtuples[i].c_str());
    if (doRecoEle)    ftreeRecoEle  -> Add(listNtuples[i].c_str());
    if (doRecoMuon)   ftreeRecoMuon -> Add(listNtuples[i].c_str());
    if (doRecoTau)    ftreeRecoTau  -> Add(listNtuples[i].c_str());

  }

  if (doEvent)    fChain->AddFriend(ftreeEvent);
  if (doreco)       fChain->AddFriend(ftreereco);
  if (domuonreco)   fChain->AddFriend(ftreemuon);
  if (dol1extra)    fChain->AddFriend(ftreeExtra);
  if (dol1emuextra) fChain->AddFriend(ftreeEmuExtra);
  if (dol1menu)     fChain->AddFriend(ftreeMenu);
  if (doRecoJet)    fChain->AddFriend(ftreeRecoJet);
  if (doRecoMet)    fChain->AddFriend(ftreeRecoMet);
  if (doRecoEle)    fChain->AddFriend(ftreeRecoEle);
  if (doRecoMuon)   fChain->AddFriend(ftreeRecoMuon);
  if (doRecoTau)    fChain->AddFriend(ftreeRecoTau);

  return true;
}

L1Ntuple::~L1Ntuple()
{
  if (ftreeRecoJet)  delete ftreeRecoJet;
  if (ftreeRecoMet)  delete ftreeRecoMet;
  if (ftreeRecoEle)  delete ftreeRecoEle;
  if (ftreeRecoMuon) delete ftreeRecoMuon;
  if (ftreeRecoTau)  delete ftreeRecoTau;
  if (ftreemuon)     delete ftreemuon;
  if (ftreereco)     delete ftreereco;
  if (ftreeExtra)    delete ftreeExtra;
  if (ftreeEmuExtra) delete ftreeEmuExtra;
  if (ftreeMenu)     delete ftreeMenu;
  if (fChain)        delete fChain;
  if (rf)            delete rf;
}



Int_t L1Ntuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t L1Ntuple::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);

   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
   }
   return centry;
}

void L1Ntuple::Init()
{
   if (!fChain) return;
   fCurrent = -1;
   /*
   fChain->SetMakeClass(1);
   ftreemuon->SetMakeClass(1);
   ftreereco->SetMakeClass(1);
   ftreeExtra->SetMakeClass(1); */

   std::cout << "Estimate the number of entries ... " <<std::endl;

   upgrade_      = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
   std::cout<<"Setting branch addresses for L1Upgrade tree...  "<<std::endl;
   fChain->SetBranchAddress("L1Upgrade", &upgrade_ );

   if (doEvent){
     std::cout<<"Setting branch addresses for Event tree..."<<std::endl;
     event_        = new L1Analysis::L1AnalysisEventDataFormat();
     
     ftreeEvent->SetBranchAddress("Event", &event_ );
     fChain-> AddFriend(ftreeEvent);     
   }

   if (doRecoJet)
   {
     std::cout<<"Setting branch addresses for reco Jets...   "<<std::endl;
     recoJet_ = new L1Analysis::L1AnalysisRecoJetDataFormat();
     recoSum_ = new L1Analysis::L1AnalysisRecoMetDataFormat();
     ftreeRecoJet->SetBranchAddress("Jet",&recoJet_);
     ftreeRecoJet->SetBranchAddress("Sums",&recoSum_);
   }



   if (doRecoMet)
   {
     //std::cout<<"Setting branch addresses for reco Jets...   "<<std::endl;
     //recoJet_ = new L1Analysis::L1AnalysisRecoJetDataFormat();
     //recoSum_ = new L1Analysis::L1AnalysisRecoMetDataFormat();
     //ftreeRecoJet->SetBranchAddress("Jet",&recoJet_);
     //ftreeRecoJet->SetBranchAddress("Sums",&recoSum_);
   }



   if (doRecoEle)
   {
     std::cout<<"Setting branch addresses for reco Eles...   "<<std::endl;
     recoEle_ = new L1Analysis::L1AnalysisRecoElectronDataFormat();
     ftreeRecoEle->SetBranchAddress("Electron",&recoEle_);
   }

   if (doRecoMuon)
   {
     std::cout<<"Setting branch addresses for reco Muons...   "<<std::endl;
     recoMuon_ = new L1Analysis::L1AnalysisRecoMuon2DataFormat();
     ftreeRecoMuon->SetBranchAddress("Muon",&recoMuon_);
   }

   if (doRecoTau)
   {
     std::cout<<"Setting branch addresses for reco Taus...   "<<std::endl;
     recoTau_ = new L1Analysis::L1AnalysisRecoTauDataFormat();
     ftreeRecoTau->SetBranchAddress("Tau",&recoTau_);
   }

   // if (fChain->GetBranch("Simulation"))
   //   fChain->SetBranchAddress("Simulation", &simulation_ );
   // else
   //   std::cout<<"Simulation Branch not added..."<<std::endl;
   // 
   // fChain->SetBranchAddress("GCT",   &gct_   );
   // fChain->SetBranchAddress("GMT",   &gmt_   );
   // fChain->SetBranchAddress("GT",    &gt_    );
   // fChain->SetBranchAddress("RCT",   &rct_   );
   // fChain->SetBranchAddress("CSCTF", &csctf_ );
   // fChain->SetBranchAddress("DTTF",  &dttf_  );



     //if (doreco)
     //  {
     // std::cout<<"Setting branch addresses for reco..."<<std::endl;
     //
     // recoMet_          = new L1Analysis::L1AnalysisRecoMetDataFormat();
     // recoJet_          = new L1Analysis::L1AnalysisRecoJetDataFormat();
     // recoBasicCluster_ = new L1Analysis::L1AnalysisRecoClusterDataFormat();
     // recoSuperCluster_ = new L1Analysis::L1AnalysisRecoClusterDataFormat();
     // recoVertex_       = new L1Analysis::L1AnalysisRecoVertexDataFormat();
     // recoTrack_        = new L1Analysis::L1AnalysisRecoTrackDataFormat();
     //
     // ftreereco->SetBranchAddress("Jet", &recoJet_);
     // ftreereco->SetBranchAddress("BasicClusters", &recoBasicCluster_);
     // ftreereco->SetBranchAddress("SuperClusters", &recoSuperCluster_);
     // ftreereco->SetBranchAddress("Met", &recoMet_);
     // ftreereco->SetBranchAddress("Tracks", &recoTrack_);
     // ftreereco->SetBranchAddress("Vertices", &recoVertex_);
     // fChain->AddFriend(ftreereco);
     // std::cout << "RecoTree "  << ftreereco->GetEntries() << endl;
     //  }
     //
     //if (domuonreco)
     //  {
     //  std::cout<<"Setting branch addresses for muons...   "<<std::endl;
     //  recoMuon_     = new L1Analysis::L1AnalysisRecoMuonDataFormat() ;
     //  recoRpcHit_   = new L1Analysis::L1AnalysisRecoRpcHitDataFormat();
     //  ftreemuon->SetBranchAddress("Muon",&recoMuon_);
     //  ftreemuon->SetBranchAddress("RpcHit",&recoRpcHit_);
     //  }
     //
     //if (dol1extra)
     //  {
     //  std::cout<<"Setting branch addresses for L1Extra... "<<std::endl;
     //  l1extra_ = new L1Analysis::L1AnalysisL1ExtraDataFormat();
     //  ftreeExtra->SetBranchAddress("L1Extra",&l1extra_);
     //  fChain->AddFriend(ftreeExtra);
     //  std::cout  << "L1ExtraTree: "<< ftreeExtra->GetEntries() << std::endl;
     //  }
     //if ( dol1emuextra){
     //      std::cout<<"Setting branch addresses for L1EmuExtra... "<<std::endl;
     //  l1emuextra_ = new L1Analysis::L1AnalysisL1ExtraDataFormat();
     //  ftreeEmuExtra->SetBranchAddress("L1Extra",&l1emuextra_);
     //  fChain->AddFriend(ftreeEmuExtra);
     //  std::cout  << "L1EmuExtraTree: "<< ftreeEmuExtra->GetEntries() << std::endl;
     //}
     //
     //if (dol1menu)
     //  {
     //  std::cout<<"Setting branch addresses for L1Menu... "<<std::endl;
     //  l1menu_ = new L1Analysis::L1AnalysisL1MenuDataFormat();
     //  ftreeMenu->SetBranchAddress("L1Menu",&l1menu_);
     //  }

}

void L1Ntuple::Test()
{ 

  if (fChain == 0)  return;
 
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  unsigned int nevents =0;

  std::cout << nentries << " events to process"<<std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    //fChain->GetEvent(jentry);
  
    nevents++;
    if (nevents<9)  //eight first events
      { 
    std::cout << "--------------------- Event "<<jentry<<" ---------------------"<<std::endl;

    //event_
    std::cout << "L1Tree         : event_->run = "<<event_->run<<std::endl;

    }
  }
   
}

// ===  FUNCTION  ============================================================
//         Name:  L1Ntuple::PrintRecoTree
//  Description:  
// ===========================================================================
bool L1Ntuple::PrintRecoTree() const
{
  if (recoJet_ != NULL)
  {
    for (int i = 0; i < recoJet_->nJets; ++i)
    {
      std::cout << recoJet_->eta.at(i) << " ";
    }
    std::cout << "." << std::endl;
  }
  return true;
}       // -----  end of function L1Ntuple::PrintRecoTree  -----

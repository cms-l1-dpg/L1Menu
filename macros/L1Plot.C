// ===========================================================================
// 
//       Filename:  L1Plot.C
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  01/31/2016 03:54:37 PM
//       Compiler:  g++ -std=c++11
// 
//         Author:  Zhenbin Wu (benwu)
//          Email:  zhenbin.wu@gmail.com
//        Company:  UIC, CMS@LPC, CDF@FNAL
// 
// ===========================================================================

#include "L1Plot.h"

//----------------------------------------------------------------------------
//       Class:  L1Plot
//      Method:  L1Plot
// Description:  constructor
//----------------------------------------------------------------------------
L1Plot::L1Plot (
    TFile* outrootfile_,
    L1Analysis::L1AnalysisEventDataFormat         *event__,
    L1Analysis::L1AnalysisL1UpgradeDataFormat     *upgrade__,
    L1Analysis::L1AnalysisRecoJetDataFormat       *recoJet__,
    L1Analysis::L1AnalysisRecoMetDataFormat       *recoSum__,
    L1Analysis::L1AnalysisRecoElectronDataFormat  *recoEle__,
    L1Analysis::L1AnalysisRecoMuon2DataFormat     *recoMuon__,
    L1Analysis::L1AnalysisRecoTauDataFormat       *recoTau__,
    L1Analysis::L1AnalysisRecoMetFilterDataFormat *recoFilter__,
    L1Analysis::L1AnalysisL1CaloTowerDataFormat   *l1CaloTower__,
    L1Analysis::L1AnalysisRecoVertexDataFormat    *recoVtx__
    ):
  outfile(outrootfile_),
  event_(event__),
  upgrade_(upgrade__),
  recoJet_(recoJet__),
  recoSum_(recoSum__),
  recoEle_(recoEle__),
  recoMuon_(recoMuon__),
  recoTau_(recoTau__),
  recoFilter_(recoFilter__),
  l1CaloTower_(l1CaloTower__),
  recoVtx_(recoVtx__),
  doPlotRate(false),
  doPlotEff(false)
{

}  // -----  end of method L1Plot::L1Plot  (constructor)  -----

//----------------------------------------------------------------------------
//       Class:  L1Plot
//      Method:  ~L1Plot
// Description:  destructor
//----------------------------------------------------------------------------
L1Plot::~L1Plot ()
{
}  // -----  end of method L1Plot::-L1Plot  (destructor)  -----

//----------------------------------------------------------------------------
//       Class:  L1Plot
//      Method:  operator =
// Description:  assignment operator
//----------------------------------------------------------------------------
  L1Plot&
L1Plot::operator = ( const L1Plot &other )
{
  if ( this != &other ) {
  }
  return *this;
}  // -----  end of method L1Plot::operator =  (assignment operator)  ---

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::BookRateHistogram
//  Description:  
// ===========================================================================
bool L1Plot::BookRateHistogram()
{
  if (!doPlotRate) return false;
  //Event Counter
  hRate1F["nEvts"]       = new TH1F("nEvts","Number of Events Processed",1,-0.5,.5);
  //Single stuff
  hRate1F["nJetVsPt"]    = new TH1F("nJetVsPt","SingleJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hRate1F["nJetCenVsPt"] = new TH1F("nJetCenVsPt","SingleJetCentral; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hRate1F["nEGVsPt"]     = new TH1F("nEGVsPt","SingleEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hRate1F["nEGErVsPt"]   = new TH1F("nEGErVsPt","SingleEGer; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hRate1F["nIsoEGVsPt"]  = new TH1F("nIsoEGVsPt","SingleIsoEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hRate1F["nIsoEGerVsPt"]  = new TH1F("nIsoEGerVsPt","SingleIsoEGer; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hRate1F["nMuVsPt"]     = new TH1F("nMuVsPt","SingleMu; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hRate1F["nMuErVsPt"]   = new TH1F("nMuErVsPt","SingleMu |#eta|<2.1; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hRate1F["nMuVsEta"]    = new TH1F("nMuVsEta","nMuVsEta",24,-2.4,2.4);
  hRate1F["nEGVsEta"]    = new TH1F("nEGVsEta","nEGVsEta",50,-3.,3.);
  hRate1F["nIsoEGVsEta"] = new TH1F("nIsoEGVsEta","nIsoEGVsEta",50,-3.,3.);
  hRate1F["nJetVsEta"]   = new TH1F("nJetVsEta","nJetVsEta",50,-5.,5.);
  hRate1F["nJetVsEta_Central"] = new TH1F("nJetVsEta_Central","nJetVsEta_Central",50,-5.,5.);
  hRate1F["nJetVsEta_Fwd"]     = new TH1F("nJetVsEta_Fwd","nJetVsEta_Fwd",50,-5.,5.);
  hRate1F["nTauVsPt"]     = new TH1F("nTauVsPt","SingleTau; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hRate1F["nTauErVsPt"]   = new TH1F("nTauErVsPt","SingleTauer; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hRate1F["nIsoTauVsPt"]  = new TH1F("nIsoTauVsPt","SingleIsoTau; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hRate1F["nIsoTauErVsPt"]  = new TH1F("nIsoTauErVsPt","SingleIsoTauEr; E_{T} cut; rate [Hz]",65,-0.5,64.5);

  //Multistuff
  hRate1F["nDiJetVsPt"]        = new TH1F("nDiJetVsPt","DiJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hRate1F["nDiCenJetVsPt"]     = new TH1F("nDiCenJetVsPt","DiCenJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hRate2F["nAsymDiJetVsPt"]    = new TH2F("nAsymDiJetVsPt","DiJet; E_{T} cut jet 1; E_{T} cut jet 2",256,-0.5,255.5,256,-0.5,255.5);
  hRate2F["nAsymDiCenJetVsPt"] = new TH2F("nAsymDiCenJetVsPt","DiCenJet; E_{T} cut jet 1; E_{T} cut jet 2",256,-0.5,255.5,256,-0.5,255.5);
  hRate1F["nQuadJetVsPt"]      = new TH1F("nQuadJetVsPt","QuadJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hRate1F["nQuadCenJetVsPt"]   = new TH1F("nQuadCenJetVsPt","QuadCenJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hRate1F["nDiTauVsPt"]        = new TH1F("nDiTauVsPt","DiTau; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hRate1F["nDiIsoTauVsPt"]     = new TH1F("nDiIsoTauVsPt","DiIsoTau; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hRate1F["nDiEGVsPt"]         = new TH1F("nDiEGVsPt","DiEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hRate1F["nDiIsoEGVsPt"]      = new TH1F("nDiIsoEGVsPt","DiIsoEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hRate2F["nEGPtVsPt"]         = new TH2F("nEGPtVsPt","DoubleEle; p_{T} cut EG_{1}; p_{T} cut EG_{2}",65,-0.5,64.5,65,-0.5,64.5);
  hRate2F["nIsoEGPtVsPt"]      = new TH2F("nIsoEGPtVsPt","DoubleIsolEle; p_{T} cut EG_{1}; p_{T} cut EG_{2}",65,-0.5,64.5,65,-0.5,64.5);
  hRate2F["nMuPtVsPt"]         = new TH2F("nMuPtVsPt","DoubleMu; p_{T} cut mu_{1}; p_{T} cut mu_{2}",41,-0.25,20.25,41,-0.25,20.25);
  hRate2F["nOniaMuPtVsPt"]     = new TH2F("nOniaMuPtVsPt","DoubleMu_Er_HighQ_WdEta22 (Quarkonia); p_{T} cut mu_{1}; p_{T} cut mu_{2}",41,-0.25,20.25,41,-0.25,20.25);

  //Sums
  hRate1F["nHTTVsHTT"] = new TH1F("nHTTVsHTT","HTT; HTT cut; rate [Hz]",512,-.5,511.5);
  hRate1F["nHTMVsHTM"] = new TH1F("nHTMVsHTM","HTM; HTM cut; rate [Hz]",512,-.5,511.5);
  hRate1F["nETTVsETT"] = new TH1F("nETTVsETT","ETT; ETT cut; rate [Hz]",512,-.5,511.5);
  hRate1F["nETMVsETM"] = new TH1F("nETMVsETM","ETM; ETM cut; rate [Hz]",512,-.5,511.5);

  return true;
}       // -----  end of function L1Plot::BookRateHistogram  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::BookTestHistogram
//  Description:  
// ===========================================================================
bool L1Plot::BookTestHistogram()
{
  if (!doPlotTest) return false;

//**************************************************************************//
//                                MET Studies                               //
//**************************************************************************//
  hTest1F["L1MET"]          = new TH1F("L1MET","L1MET",                                          150, 0,    150);
  hTest1F["PFMET"]          = new TH1F("PFMET","PFMET",                                          150, 0,    150);
  hTest1F["L1METPFMET"]     = new TH1F("L1METPFMET","L1METPFMET;L1MET/PFMET;Events",             100, 0, 5);
  hTest1F["L1METPFMETNoMu"] = new TH1F("L1METPFMETNoMu","L1METPFMETNoMu;L1MET/PFMETNoMu;Events", 100, 0, 5);
  hTest1F["PFMETVNoMu"]     = new TH1F("PFMETVNoMu","PFMETVNoMu;PFMET-PFMETNoMu;Events",         200, -100, 100);
  hTest1F["L1METx"]         = new TH1F("L1METx","L1METx",                                        100, -50,  50);
  hTest1F["L1METy"]         = new TH1F("L1METy","L1METy",                                        100, -50,  50);
  hTest1F["PFMETx"]         = new TH1F("PFMETx","PFMETx",                                        100, -50,  50);
  hTest1F["PFMETy"]         = new TH1F("PFMETy","PFMETy",                                        100, -50,  50);

  hTest2F["L1METResVsRecoAct"] = new TH2F("L1METResVsRecoAct"," METResVsRecoAct; Reco HT activity; L1MET/PFMET", 20,  0, 1,   50, 0, 5);
  hTest2F["L1METResVsPFMET"]   = new TH2F("L1METResVsPFMET"," METResVsPFMET; Reco PFMET; L1MET/PFMET",           100, 0, 200, 50, 0, 5);
  hTest2F["L1METResVsNvtx"]    = new TH2F("L1METResVsNvtx","METResVsNvtx; Nvtx; L1MET/PFMET",                    30,  0, 30,  50, 0, 5);
  hTest2F["L1METxVsNvtx"]   = new TH2F("L1METxVsNvtx","L1METxVsNvtx; Nvtx; L1METx",          30, 0,   30, 100, -50, 50);
  hTest2F["L1METyVsNvtx"]   = new TH2F("L1METyVsNvtx","L1METyVsNvtx; Nvtx; L1METy",          30, 0,   30, 100, -50, 50);
  hTest2F["PFMETxVsNvtx"]   = new TH2F("PFMETxVsNvtx","PFMETxVsNvtx; Nvtx; PFMETx",          30, 0,   30, 100, -50, 50);
  hTest2F["PFMETyVsNvtx"]   = new TH2F("PFMETyVsNvtx","PFMETyVsNvtx; Nvtx; PFMETy",          30, 0,   30, 100, -50, 50);
  hTest2F["L1METxVsRecoAct"]   = new TH2F("L1METxVsRecoAct","L1METxVsRecoAct; RecoAct; L1METx",          20, 0,   1, 100, -50, 50);
  hTest2F["L1METyVsRecoAct"]   = new TH2F("L1METyVsRecoAct","L1METyVsRecoAct; RecoAct; L1METy",          20, 0,   1, 100, -50, 50);
  hTest2F["PFMETxVsRecoAct"]   = new TH2F("PFMETxVsRecoAct","PFMETxVsRecoAct; RecoAct; PFMETx",          20, 0,   1, 100, -50, 50);
  hTest2F["PFMETyVsRecoAct"]   = new TH2F("PFMETyVsRecoAct","PFMETyVsRecoAct; RecoAct; PFMETy",          20, 0,   1, 100, -50, 50);

  hTestPro["ECalVsiEta"]     = new TProfile("ECalVsiEta","ECalVsiEta; iEta; ECalTower",           84, -42, 42);
  hTestPro["HCalVsiEta"]     = new TProfile("HCalVsiEta","HCalVsiEta; iEta; HCalTower",           84, -42, 42);
  hTestPro["CaloVsiEta"]     = new TProfile("CaloVsiEta","CaloVsiEta; iEta; CaloTower",           84, -42, 42);



//**************************************************************************//
//                                Muon Study                                //
//**************************************************************************//
  hTest2F["Mu5BxEta"]  = new TH2F("Mu5BxEta","Mu5BxEta; muonEta (pT>5); muonBX",    50, -2.5, 2.5, 5, -2, 2);
  hTest2F["Mu16BxEta"] = new TH2F("Mu16BxEta","Mu16BxEta; muonEta (pT>16); muonBX", 50, -2.5, 2.5, 5, -2, 2);
  hTest2F["Mu25BxEta"] = new TH2F("Mu25BxEta","Mu25BxEta; muonEta (pT>25); muonBX", 50, -2.5, 2.5, 5, -2, 2);
  return true;
}       // -----  end of function L1Plot::BookTestHistogram  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::WriteTestHistogram()
//  Description:  
// ===========================================================================
bool L1Plot::WriteTestHistogram() const
{
  if (!doPlotTest) return false;
  outfile->mkdir("Test");
  outfile->cd("Test");
  for(auto f : hTest2F)
  {
    f.second->Write();
  }
  for(auto f : hTest1F)
  {
    f.second->Write();
  }
  for(auto f : hTestPro)
  {
    f.second->Write();
  }
  outfile->cd();
  return true;
}       // -----  end of function L1Plot::WriteTestHistogram()  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::WriteRateHistogram
//  Description:  
// ===========================================================================
bool L1Plot::WriteRateHistogram(double scale) const
{
  if (!doPlotRate) return false;
  outfile->mkdir("Rate");
  outfile->cd("Rate");
  for(auto f : hRate1F)
  {
    f.second->Sumw2();
    f.second->Scale(scale);
    f.second->Write();
  }
  
  for(auto f : hRate2F)
  {
    f.second->Sumw2();
    f.second->Scale(scale);
    f.second->Write();
  }

  outfile->cd();
  
  return true;
}       // -----  end of function L1Plot::WriteRateHistogram  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::FillRateHistogram
//  Description:  
// ===========================================================================
bool L1Plot::FillRateHistogram()
{
  if (!doPlotRate) return false;
  for(int ptCut=0; ptCut<256; ++ptCut) {
    if(L1Event->JetPt>=ptCut)	  hRate1F["nJetVsPt"]->Fill(ptCut);
    if(L1Event->JetCenPt>=ptCut) hRate1F["nJetCenVsPt"]->Fill(ptCut);
    if(L1Event->TauPt>=ptCut)	  hRate1F["nTauVsPt"]->Fill(ptCut);
    if(L1Event->TauerPt>=ptCut)	  hRate1F["nTauErVsPt"]->Fill(ptCut);
    if(L1Event->IsoTauPt>=ptCut)	  hRate1F["nIsoTauVsPt"]->Fill(ptCut);
    if(L1Event->IsoTauerPt>=ptCut)	  hRate1F["nIsoTauErVsPt"]->Fill(ptCut);

    if(L1Event->dijetPt2>=ptCut)
    {
      hRate1F["nDiJetVsPt"]->Fill(ptCut);

      for(int ptCut_0=ptCut; ptCut_0<256; ++ptCut_0) 
      {
        if(L1Event->dijetPt1>=ptCut_0) hRate2F["nAsymDiJetVsPt"]->Fill(ptCut_0,ptCut);
      }
    }

    if(L1Event->diCenjetPt2>=ptCut)
    {
      hRate1F["nDiCenJetVsPt"]->Fill(ptCut);
      for(int ptCut_0=ptCut; ptCut_0<256; ++ptCut_0)
      {
        if(L1Event->diCenjetPt1>=ptCut_0) hRate2F["nAsymDiCenJetVsPt"]->Fill(ptCut_0,ptCut);
      }
    }

    if(L1Event->ditauPt>=ptCut)    hRate1F["nDiTauVsPt"]->Fill(ptCut);
    if(L1Event->diIsotauPt>=ptCut) hRate1F["nDiIsoTauVsPt"]->Fill(ptCut);
    if(L1Event->quadjetPt>=ptCut)  hRate1F["nQuadJetVsPt"]->Fill(ptCut);
    if(L1Event->quadjetCPt>=ptCut) hRate1F["nQuadCenJetVsPt"]->Fill(ptCut);

  } //loop on 256




  for(int ptCut=0; ptCut<65; ++ptCut) {
    if(L1Event->EGPt>=ptCut)    hRate1F["nEGVsPt"]->Fill(ptCut);
    if(L1Event->EGerPt>=ptCut)  hRate1F["nEGErVsPt"]->Fill(ptCut);
    if(L1Event->IsoEGPt>=ptCut) hRate1F["nIsoEGVsPt"]->Fill(ptCut);
    if(L1Event->IsoEGerPt>=ptCut) hRate1F["nIsoEGerVsPt"]->Fill(ptCut);
    if(L1Event->diEG2>=ptCut)     hRate1F["nDiEGVsPt"]->Fill(ptCut);
    if(L1Event->diIsolEG2>=ptCut) hRate1F["nDiIsoEGVsPt"]->Fill(ptCut);
     

    for(int ptCut2=0; ptCut2<=65; ++ptCut2) {
      if(L1Event->diEG1>=ptCut && L1Event->diEG2>=ptCut2 && ptCut2 <= ptCut) hRate2F["nEGPtVsPt"]->Fill(ptCut,ptCut2);
      if(L1Event->diIsolEG1>=ptCut && L1Event->diIsolEG2>=ptCut2 && ptCut2<= ptCut) hRate2F["nIsoEGPtVsPt"]->Fill(ptCut,ptCut2);
    }
     
  }//loop on 65
    
   for(int ptCut=0; ptCut<131; ++ptCut) {
     if (L1Event->MuPt>=ptCut)    hRate1F["nMuVsPt"]->Fill(ptCut);
    if (L1Event->MuerPt>=ptCut)  hRate1F["nMuErVsPt"]->Fill(ptCut);
   }
     
   for(int iCut=0; iCut<41; ++iCut) {
     for(int iCut2=0; iCut2<=iCut; ++iCut2) {
       float ptCut = iCut*0.5;
       float ptCut2 = iCut2*0.5;
       if (L1Event->doubleMuPt1>=ptCut && L1Event->doubleMuPt2>=ptCut2) hRate2F["nMuPtVsPt"]->Fill(ptCut,ptCut2);
       if (L1Event->oniaMuPt1>=ptCut && L1Event->oniaMuPt2>=ptCut2)     hRate2F["nOniaMuPtVsPt"]->Fill(ptCut,ptCut2);
     }
   }
   
  for(int httCut=0; httCut<512; ++httCut) {
    if(L1Event->HTT>httCut) hRate1F["nHTTVsHTT"]->Fill(httCut);
    if(L1Event->HTM>httCut) hRate1F["nHTMVsHTM"]->Fill(httCut);
    if(L1Event->ETT>httCut) hRate1F["nETTVsETT"]->Fill(httCut);
    if(L1Event->ETM>httCut) hRate1F["nETMVsETM"]->Fill(httCut);
  }

  //std::cout <<  "sumJetHt " << sumJetHt <<" HTT " << L1Event->HTT << std::endl;
  return true;
}       // -----  end of function L1Plot::FillRateHistogram  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::RunPlot
//  Description:  
// ===========================================================================
bool L1Plot::RunPlot()
{
  if (doPlotRate)
  {
    FillRateHistogram();
  }

  if (doPlotEff) 
  {
    FillEffHistogram();
  }

  if (doPlotTest) 
  {
    TestMETActivity();
    TestMuon();
  }

  return true;
}       // -----  end of function L1Plot::RunPlot  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::PreRun
//  Description:  
// ===========================================================================
bool L1Plot::PreRun( StructL1Event *L1Event_, std::map<std::string, L1Seed> *mL1Seed_)
{
  L1Event  = L1Event_;
  mL1Seed = mL1Seed_;
  BookRateHistogram();
  BookEffHistogram();
  BookTestHistogram();
  return true;
}       // -----  end of function L1Plot::PreRun  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::PostRun
//  Description:  
// ===========================================================================
bool L1Plot::PostRun(double scale)
{
  WriteRateHistogram(scale);
  WriteEffHistogram();
  WriteTestHistogram();
  return true;
}       // -----  end of function L1Plot::PostRun  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::SortVTLVs
//  Description:  /* cursor */ 
// ===========================================================================
bool L1Plot::SortVTLVs(std::vector<TLorentzVector> &reTLVs) const
{
  std::sort(reTLVs.begin(), reTLVs.end(), 
      [](TLorentzVector a, TLorentzVector b){ return a.Pt() > b.Pt();});
  return true;
}       // -----  end of function L1Plot::SortVTLVs  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::GetRecoJet
//  Description:  /* cursor */
// ===========================================================================
std::vector<TLorentzVector> L1Plot::GetRecoJet(bool isCent) const
{
  std::vector<TLorentzVector> reTLVs;
  const float jetERcut = 3.0;
  if (recoJet_ == NULL) return reTLVs;
 
  for (int i = 0; i < recoJet_->nJets; ++i)
  {
    if (isCent && fabs(recoJet_->eta.at(i)) > jetERcut )
      continue;

    if (!recoJet_->isPF.at(i)) continue;
    
    TLorentzVector temp(0, 0, 0, 0);


    temp.SetPtEtaPhiE( recoJet_->etCorr.at(i), recoJet_->eta.at(i), 
        recoJet_->phi.at(i), recoJet_->e.at(i));
    reTLVs.push_back(temp);
  }

  SortVTLVs(reTLVs);

  return reTLVs;
}       // -----  end of function L1Plot::GetRecoJet  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::GoodRecoJet
//  Description:  
// ===========================================================================
bool L1Plot::GoodRecoJet(int ijet) const
{
  bool tmp = true;
  tmp &= recoJet_ ->nhef[ijet] < 0.9 ;
  tmp &= recoJet_ ->nemef[ijet] < 0.9 ;
  tmp &= (recoJet_ ->cMult[ijet] + recoJet_ ->nMult[ijet]) > 1 ;
  tmp &= recoJet_ ->mef[ijet] < 0.8 ;
  if (fabs(recoJet_ ->eta[ijet]) < 2.4) {
    tmp &= recoJet_ ->chef[ijet] > 0.0 ;
    tmp &= recoJet_ ->cMult[ijet] > 0 ;
    tmp &= recoJet_ ->cemef[ijet] < 0.9 ;
  }
  if (fabs(recoJet_ ->eta[ijet]) > 3.0) {
    tmp &= recoJet_ ->nemef[ijet] < 0.9 ;
    tmp &= recoJet_ ->nMult[ijet] > 10 ;
  }

  // our custom selection
  tmp &= recoJet_ ->muMult[ijet] == 0;
  tmp &= recoJet_ ->elMult[ijet] == 0;

  return tmp;
}       // -----  end of function L1Plot::GoodRecoJet  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::GetRecoHTLocal
//  Description:  From Jim Brooke
//  offline HT is computed by hand from PF jets with V6 JEC, eta < 3 and pT > 30.
// ===========================================================================
std::vector<TLorentzVector> L1Plot::GetRecoHTLocal() const
{
  double ht = 0;

  for (int i = 0; i < recoJet_->nJets; ++i)
  {
    if (!recoJet_->isPF.at(i)) continue;
    if (fabs(recoJet_->eta.at(i)) > 3) continue;
    if (recoJet_->etCorr.at(i) < 30) continue;
    ht += recoJet_->etCorr.at(i);
  }
  std::vector<TLorentzVector> reTLVs;
  TLorentzVector temp(0, 0, 0, 0);
  temp.SetPtEtaPhiE(ht, 0, 0, 0);
  reTLVs.push_back(temp);

  return reTLVs;
}       // -----  end of function L1Plot::GetRecoHTLocal  -----
// ===  FUNCTION  ============================================================
//         Name:  L1Plot::GetRecoSum
//  Description:  
// ===========================================================================
std::vector<TLorentzVector> L1Plot::GetRecoSum(std::string type ) const
{
  std::vector<TLorentzVector> reTLVs;
  if (!recoSum_) return reTLVs;

  TLorentzVector temp(0, 0, 0, 0);

  if (type == "HTT")
  {
    temp.SetPtEtaPhiE(recoSum_->Ht, 0, 0, 0);
  }

  if (type == "ETM")
  {
    temp.SetPtEtaPhiE(recoSum_->met, 0, recoSum_->metPhi, 0);
  }

  if (type == "HTM")
  {
    temp.SetPtEtaPhiE(recoSum_->mHt, 0, recoSum_->mHtPhi, 0);
  }
  if (type == "ETT")
  {
    temp.SetPtEtaPhiE(recoSum_->sumEt, 0, 0, 0);
  }

  reTLVs.push_back(temp);
  return reTLVs;
}       // -----  end of function L1Plot::GetRecoSum  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::GetRecoEle
//  Description:  /* cursor */
// ===========================================================================
std::vector<TLorentzVector> L1Plot::GetRecoEle(bool isER, float IsoCut, int qual) const
{
  std::vector<TLorentzVector> reTLVs;
  const float EleERcut = 2.1;
  if (recoEle_ == NULL) return reTLVs;

  for (unsigned int i = 0; i < recoEle_->nElectrons; ++i)
  {
    if (qual == -1 && ! recoEle_ -> isVetoElectron.at(i)) continue;
    if (qual == 1 && ! recoEle_ -> isLooseElectron.at(i)) continue;
    if (qual == 2 && ! recoEle_ -> isMediumElectron.at(i)) continue;
    if (qual == 3 && ! recoEle_ -> isTightElectron.at(i)) continue;

    if (isER && fabs(recoEle_->eta.at(i)) > EleERcut )
      continue;
    if (recoEle_->iso.at(i) < IsoCut) continue;

    TLorentzVector temp(0, 0, 0, 0);
    temp.SetPtEtaPhiE( recoEle_->pt.at(i),
        recoEle_->eta.at(i),
        recoEle_->phi.at(i),
        recoEle_->e.at(i));
    reTLVs.push_back(temp);
  }

  SortVTLVs(reTLVs);

  return reTLVs;
}       // -----  knd kf kunction k1Plot::ketRecoEle  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::GetRecoMuon
//  Description:  
// ===========================================================================
std::vector<TLorentzVector> L1Plot::GetRecoMuon(float MuERcut, float IsoCut, int qual) const
{
  std::vector<TLorentzVector> reTLVs;
  if (recoMuon_ == NULL) return reTLVs;

  for (int i = 0; i < recoMuon_->nMuons; ++i)
  {
    if (qual == 1 && ! recoMuon_ -> isLooseMuon.at(i)) continue;
    if (qual == 2 && ! recoMuon_ -> isMediumMuon.at(i)) continue;

    if (fabs(recoMuon_->eta.at(i)) > MuERcut )
      continue;

    if (recoMuon_->iso.at(i) >= IsoCut) continue;

    TLorentzVector temp(0, 0, 0, 0);
    temp.SetPtEtaPhiE( recoMuon_->pt.at(i),
        recoMuon_->eta.at(i),
        recoMuon_->phi.at(i),
        recoMuon_->e.at(i));
    reTLVs.push_back(temp);
  }

  SortVTLVs(reTLVs);

  return reTLVs;
}       // -----  end of function L1Plot::GetRecoMuon  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::GetRecoTau
//  Description:  
// ===========================================================================
std::vector<TLorentzVector> L1Plot::GetRecoTau(bool isER, int Iso) const
{
  std::vector<TLorentzVector> reTLVs;
  const float TauERcut = 3.0;
  if (recoTau_ == NULL) return reTLVs;

  for (unsigned int i = 0; i < recoTau_->nTaus; ++i)
  {
    if (Iso == 1 && ! recoTau_ -> LooseIsoFlag.at(i)) continue;
    if (Iso == 2 && ! recoTau_ -> TightIsoFlag.at(i)) continue;

    if (isER && fabs(recoTau_->eta.at(i)) > TauERcut )
      continue;


    TLorentzVector temp(0, 0, 0, 0);
    temp.SetPtEtaPhiE( recoTau_->pt.at(i),
        recoTau_->eta.at(i),
        recoTau_->phi.at(i),
        recoTau_->e.at(i));
    reTLVs.push_back(temp);
  }

  SortVTLVs(reTLVs);

  return reTLVs;
}       // -----  end of function L1Plot::GetRecoTau  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::GetRecoEvent
//  Description:  
// ===========================================================================
bool L1Plot::GetRecoEvent()
{
  recoEvent.clear();
  recoEvent["Jet"]     = GetRecoJet();
  recoEvent["JetC"]    = GetRecoJet(true);
  recoEvent["Tau"]     = GetRecoTau();
  recoEvent["IsoTau"]     = GetRecoTau(false, 1);
  recoEvent["Tauer"]   = GetRecoTau(true);
  recoEvent["EG"]      = GetRecoEle();
  recoEvent["EGer"]    = GetRecoEle(true,   false, 0 );
  recoEvent["IsoEG"]   = GetRecoEle(false,  true,  0 );
  recoEvent["IsoEGer"] = GetRecoEle(true,   true,  0 );
  recoEvent["Mu"]      = GetRecoMuon(99, 0.15, 2);
  recoEvent["MuOpen"]  = GetRecoMuon(99, 0.15, 1);
  recoEvent["Muer"]    = GetRecoMuon(2.1);
  recoEvent["HTT"]     = GetRecoSum("HTT");
  recoEvent["ETM"]     = GetRecoSum("ETM");
  recoEvent["ETT"]     = GetRecoSum("ETT");
  //recoEvent["ETM"]     = GetRecoHTMLocal();
  recoEvent["HTM"]     = GetRecoSum("HTM");
  //recoEvent["HTM"]     = GetRecoHTMLocal();
  //recoEvent["HTM"]     = GetRecoSum("ETM");
  return true;
}       // -----  end of function L1Plot::GetRecoEvent  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::BookEffHistogram
//  Description:  
// ===========================================================================
bool L1Plot::BookEffHistogram() 
{
  if (!doPlotEff) return false;
  // map L1Seed to hEff, func of X, object? 
  hEff.clear();

  for(const auto& l1 : *mL1Seed)
  {
    if (l1.second.singleObj != "")
    {
      std::string objstr = l1.second.singleObj;
      if (objstr.find("Jet") != std::string::npos || objstr.find("Tau") != std::string::npos )
      {
        std::string hname = l1.first +"_Pt";
        hEff[hname] = new TEfficiency(hname.c_str(), l1.first.c_str(), 60, 0, 300);
        hEffFun[hname] = std::bind(&L1Plot::FunLeadingPt, this, l1.second.singleObj);
      }

      if (objstr.find("EG") != std::string::npos)
      {
        std::string hname = l1.first +"_Pt";
        hEff[hname] = new TEfficiency(hname.c_str(), l1.first.c_str(), 50, 0, 100);
        hEffFun[hname] = std::bind(&L1Plot::FunLeadingPt, this, l1.second.singleObj);
      }

      if (objstr.find("Mu") != std::string::npos)
      {
        std::string hname = l1.first +"_Pt";
        hEff[hname] = new TEfficiency(hname.c_str(), l1.first.c_str(), 50, 0, 150);
        hEffFun[hname] = std::bind(&L1Plot::FunLeadingPt, this, l1.second.singleObj);
      }

      if (objstr.find("HTT") != std::string::npos ||
          objstr.find("ETT") != std::string::npos ||
          objstr.find("ETM") != std::string::npos ||
          objstr.find("HTM") != std::string::npos )
      {
        std::string hname = l1.first +"_Pt";
        if (objstr.find("ETM") != std::string::npos ||
            objstr.find("HTM") != std::string::npos)
          hEff[hname] = new TEfficiency(hname.c_str(), l1.first.c_str(), 100,0,500);
        else
          hEff[hname] = new TEfficiency(hname.c_str(), l1.first.c_str(), 30,0,600);
        hEffFun[hname] = std::bind(&L1Plot::FunLeadingPt, this, l1.second.singleObj);
      }
    } // End of SingleObj
  }
  return true;
}       // -----  end of function L1Plot::BookEffHistogram  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::FunLeadingPt
//  Description:  /* cursor */
// ===========================================================================
double L1Plot::FunLeadingPt(std::string obj)
{
  const std::vector<TLorentzVector> &vs = recoEvent[obj];
  if (vs.size() == 0) return -1;
  else
    return vs.front().Pt();
}       // -----  end of function L1Plot::FunLeadingPt  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::FillEffHistogram
//  Description:  
// ===========================================================================
bool L1Plot::FillEffHistogram()
{
  if (!doPlotEff) return false;
  if (!GetRecoFilter()) return false;

  GetRecoEvent();
  for(auto h : hEff)
  {
    std::string l1seed = h.second->GetTitle();
    std::string objname = (*mL1Seed)[l1seed].singleObj;
    if (recoEvent.find(objname) == recoEvent.end())
      continue;
    //std::cout << l1seed <<" " << (*mL1Seed)[l1seed].eventfire <<"  " << hEffFun[h.first]() << std::endl;
    h.second->Fill((*mL1Seed)[l1seed].eventfire, hEffFun[h.first]());
  }

  return true;
}       // -----  end of function L1Plot::FillEffHistogram  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::WriteEffHistogram
//  Description:  
// ===========================================================================
bool L1Plot::WriteEffHistogram()
{
  if (!doPlotEff) return false;
  outfile->mkdir("Eff");
  outfile->cd("Eff");
  TCanvas c1("Eff", "Eff", 600, 500);

  for(auto f : hEff)
  {
    f.second->GetPassedHistogram()->Write();
    f.second->GetTotalHistogram()->Write();
    c1.cd();
    f.second->Paint("AP");
    f.second->GetPaintedGraph()->Write( f.first.c_str());
  }

  outfile->cd();
  return true;
}       // -----  end of function L1Plot::WriteEffHistogram  -----


// ===  FUNCTION  ============================================================
//         Name:  L1Plot::SetTodo
//  Description:  
// ===========================================================================
void L1Plot::SetTodo (bool doPlotRate_, bool doPlotEff_, bool doPlotTest_)
{
  doPlotRate = doPlotRate_;
  doPlotEff = doPlotEff_;
  doPlotTest = doPlotTest_;
}       // -----  end of function L1Plot::SetTodo  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::GetRecoFilter
//  Description:  
// ===========================================================================
bool L1Plot::GetRecoFilter() const
{
  if (!recoFilter_) return true;

  bool pass = true && recoFilter_->goodVerticesFilter
                   && recoFilter_->cscTightHalo2015Filter
                   && recoFilter_->eeBadScFilter
                   && recoFilter_->ecalDeadCellTPFilter
                   && recoFilter_->hbheNoiseIsoFilter
                   && recoFilter_->hbheNoiseFilter
                   && recoFilter_->chHadTrackResFilter
                   && recoFilter_->muonBadTrackFilter;

  if (!pass) return pass; // If event is skip aleady

  for (int i = 0; i < recoJet_->nJets; ++i)
  {
    if (recoJet_->etCorr.at(i) < 30) continue;
    pass = pass && GoodRecoJet(i);
    if (!pass) return pass; // If event is skip aleady
  }

  return pass;
}       // -----  end of function L1Plot::GetRecoFilter  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::TestMETActivity
//  Description:  
// ===========================================================================
bool L1Plot::TestMETActivity()
{
  TVector2 L1MET(0, 0);
  L1MET = GetL1METCalo();
  double TheETM=L1MET.Mod();

  if (recoSum_)
  {
    double PFMET = recoSum_->met;
    double PFMETphi = recoSum_->metPhi;
    TVector2 PFMETnoMu = GetRecoMETNoMu();
    double RecoAct =TestRecoAct(1.6);
    hTest2F["L1METResVsRecoAct"]->Fill(RecoAct, TheETM/PFMET);
    hTest2F["PFMETxVsRecoAct"] -> Fill(RecoAct, PFMET * TMath::Cos(PFMETphi));
    hTest2F["PFMETyVsRecoAct"] -> Fill(RecoAct, PFMET * TMath::Sin(PFMETphi));
    hTest2F["L1METxVsRecoAct"] -> Fill(RecoAct, L1MET.Px());
    hTest2F["L1METyVsRecoAct"] -> Fill(RecoAct, L1MET.Py());

    hTest2F["L1METResVsPFMET"]->Fill(PFMET, TheETM/PFMET);
    hTest1F["PFMET"]->Fill(PFMET);
    hTest1F["L1METPFMET"]->Fill(TheETM/PFMET);
    hTest1F["PFMETVNoMu"]->Fill(PFMET-PFMETnoMu.Mod());
    hTest1F["L1METPFMETNoMu"]->Fill(TheETM/PFMETnoMu.Mod());
    hTest1F["PFMETx"]->Fill(PFMET * TMath::Cos(PFMETphi));
    hTest1F["PFMETy"]->Fill(PFMET * TMath::Sin(PFMETphi));
    if (recoVtx_)
    {
      hTest2F["PFMETxVsNvtx"] -> Fill(recoVtx_->nVtx, PFMET * TMath::Cos(PFMETphi));
      hTest2F["PFMETyVsNvtx"] -> Fill(recoVtx_->nVtx, PFMET * TMath::Sin(PFMETphi));
      hTest2F["L1METxVsNvtx"] -> Fill(recoVtx_->nVtx, L1MET.Px());
      hTest2F["L1METyVsNvtx"] -> Fill(recoVtx_->nVtx, L1MET.Py());
      hTest2F["L1METResVsNvtx"] -> Fill(recoVtx_->nVtx, TheETM / PFMET );
    }
  }

  hTest1F["L1METx"]->Fill(L1MET.Px());
  hTest1F["L1METy"]->Fill(L1MET.Py());
  hTest1F["L1MET"]->Fill(TheETM);
  return true;
}       // -----  end of function L1Plot::TestMETActivity  -----


// ===  FUNCTION  ============================================================
//         Name:  L1Plot::GetRecoMETNoMu
//  Description:  
// ===========================================================================
TVector2 L1Plot::GetRecoMETNoMu() const
{
  TVector2 revec(0,0);
  if (!recoSum_ || !recoMuon_) return revec;

  double pfmetx =  recoSum_->met * TMath::Cos(recoSum_->metPhi);
  double pfmety =  recoSum_->met * TMath::Sin(recoSum_->metPhi);
  std::vector<TLorentzVector> muons = GetRecoMuon(999, 0.15, 2);
  for(auto m : muons)
  {
    pfmetx += m.Px();
    pfmety += m.Py();
  }
  revec.Set(pfmetx, pfmety);

  return revec;
}       // -----  end of function L1Plot::GetRecoMETNoMu  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::GetL1METCalo
//  Description:  
// ===========================================================================
TVector2 L1Plot::GetL1METCalo()
{
  TVector2 revec(0,0);
  if (!l1CaloTower_) return revec;
  
  const int ietamax = 29;
  float metX =0;
  float metY =0;


  for(int jTower=0; jTower< l1CaloTower_ ->nTower; ++jTower){
    Int_t ieta = l1CaloTower_->ieta[jTower];
    Int_t iphi = l1CaloTower_->iphi[jTower];
    Int_t iet  = l1CaloTower_->iet[jTower];
    Int_t iem  = l1CaloTower_->iem[jTower];
    Int_t ihad = l1CaloTower_->ihad[jTower];
    hTestPro["ECalVsiEta"] -> Fill(ieta, (float)iem / 2.0);
    hTestPro["HCalVsiEta"] -> Fill(ieta, (float)ihad / 2.0);
    hTestPro["CaloVsiEta"] -> Fill(ieta, (float)iet / 2.0);
    Double_t phi = (Double_t)iphi * TMath::Pi()/36.;
    Double_t et = 0.5 * (Double_t)iet;
    if( abs(ieta) < ietamax){
      metX += et * TMath::Cos(phi);
      metY += et * TMath::Sin(phi);
    }
  }

  revec.Set(metX, metY);
  return revec;
}       // -----  end of function L1Plot::GetL1METCalo  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::TestRecoAct
//  Description:  
// ===========================================================================
float L1Plot::TestRecoAct(float eta) const
{
  if (recoJet_ == NULL) return -1 ;
  float selHT = 0.0;
  float totHT = 0.0;

  for (int i = 0; i < recoJet_->nJets; ++i)
  {
    if (!recoJet_->isPF.at(i)) continue;
    totHT += recoJet_->etCorr.at(i);

    if (fabs(recoJet_->eta.at(i)) >= eta && fabs(recoJet_->eta.at(i)) <= 3.0 ) 
      selHT+= recoJet_->etCorr.at(i);
  }

  return selHT/totHT;
}       // -----  end of function L1Plot::TestRecoAct  -----

// ===  FUNCTION  ============================================================
//         Name:  L1Plot::TestMuon
//  Description:  
// ===========================================================================
bool L1Plot::TestMuon()
{
  
  for(UInt_t imu=0; imu < upgrade_->nMuons; imu++) {
    if (upgrade_->muonQual.at(imu) < 12) continue;
    if(upgrade_->muonEt.at(imu) > 5)  
      hTest2F["Mu5BxEta"] -> Fill(upgrade_->muonEta.at(imu), upgrade_->muonBx.at(imu));
    if(upgrade_->muonEt.at(imu) > 16)  
      hTest2F["Mu16BxEta"] -> Fill(upgrade_->muonEta.at(imu), upgrade_->muonBx.at(imu));
    if(upgrade_->muonEt.at(imu) > 25)  
      hTest2F["Mu25BxEta"] -> Fill(upgrade_->muonEta.at(imu), upgrade_->muonBx.at(imu));
  }

  return true;
}       // -----  end of function L1Plot::TestMuon  -----

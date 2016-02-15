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
    L1Analysis::L1AnalysisEventDataFormat        *event__,
    L1Analysis::L1AnalysisL1UpgradeDataFormat    *upgrade__,
    L1Analysis::L1AnalysisRecoJetDataFormat      *recoJet__,
    L1Analysis::L1AnalysisRecoMetDataFormat      *recoSum__,
    L1Analysis::L1AnalysisRecoElectronDataFormat *recoEle__,
    L1Analysis::L1AnalysisRecoMuon2DataFormat    *recoMuon__,
    L1Analysis::L1AnalysisRecoTauDataFormat      *recoTau__
    ):
  outfile(outrootfile_),
  event_(event__),
  upgrade_(upgrade__),
  recoJet_(recoJet__),
  recoSum_(recoSum__),
  recoEle_(recoEle__),
  recoMuon_(recoMuon__),
  recoTau_(recoTau__),
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
  hRate1F["nTauVsPt"]    = new TH1F("nTauVsPt","SingleTau; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hRate1F["nJetCenVsPt"] = new TH1F("nJetCenVsPt","SingleJetCentral; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hRate1F["nEGVsPt"]     = new TH1F("nEGVsPt","SingleEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hRate1F["nEGErVsPt"]   = new TH1F("nEGErVsPt","SingleEGer; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hRate1F["nIsoEGVsPt"]  = new TH1F("nIsoEGVsPt","SingleIsoEGer; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hRate1F["nMuVsPt"]     = new TH1F("nMuVsPt","SingleMu; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hRate1F["nMuErVsPt"]   = new TH1F("nMuErVsPt","SingleMu |#eta|<2.1; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hRate1F["nMuVsEta"]    = new TH1F("nMuVsEta","nMuVsEta",24,-2.4,2.4);
  hRate1F["nEGVsEta"]    = new TH1F("nEGVsEta","nEGVsEta",50,-3.,3.);
  hRate1F["nIsoEGVsEta"] = new TH1F("nIsoEGVsEta","nIsoEGVsEta",50,-3.,3.);
  hRate1F["nJetVsEta"]   = new TH1F("nJetVsEta","nJetVsEta",50,-5.,5.);
  hRate1F["nJetVsEta_Central"] = new TH1F("nJetVsEta_Central","nJetVsEta_Central",50,-5.,5.);
  hRate1F["nJetVsEta_Fwd"]     = new TH1F("nJetVsEta_Fwd","nJetVsEta_Fwd",50,-5.,5.);

  //Multistuff
  hRate1F["nDiJetVsPt"]        = new TH1F("nDiJetVsPt","DiJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hRate1F["nDiCenJetVsPt"]     = new TH1F("nDiCenJetVsPt","DiCenJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hRate2F["nAsymDiJetVsPt"]    = new TH2F("nAsymDiJetVsPt","DiJet; E_{T} cut jet 1; E_{T} cut jet 2",256,-0.5,255.5,256,-0.5,255.5);
  hRate2F["nAsymDiCenJetVsPt"] = new TH2F("nAsymDiCenJetVsPt","DiCenJet; E_{T} cut jet 1; E_{T} cut jet 2",256,-0.5,255.5,256,-0.5,255.5);
  hRate1F["nQuadJetVsPt"]      = new TH1F("nQuadJetVsPt","QuadJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hRate1F["nQuadCenJetVsPt"]   = new TH1F("nQuadCenJetVsPt","QuadCenJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hRate1F["nDiTauVsPt"]        = new TH1F("nDiTauVsPt","DiTau; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hRate1F["nDiEGVsPt"]         = new TH1F("nDiEGVsPt","DiEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hRate1F["nDiIsoEGVsPt"]      = new TH1F("nDiIsoEGVsPt","DiIsoEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hRate2F["nEGPtVsPt"]         = new TH2F("nEGPtVsPt","DoubleEle; p_{T} cut EG_{1}; p_{T} cut EG_{2}",65,-0.5,64.5,65,-0.5,64.5);
  hRate2F["nIsoEGPtVsPt"]      = new TH2F("nIsoEGPtVsPt","DoubleIsolEle; p_{T} cut EG_{1}; p_{T} cut EG_{2}",65,-0.5,64.5,65,-0.5,64.5);
  hRate2F["nMuPtVsPt"]         = new TH2F("nMuPtVsPt","DoubleMu; p_{T} cut mu_{1}; p_{T} cut mu_{2}",41,-0.25,20.25,41,-0.25,20.25);
  hRate2F["nOniaMuPtVsPt"]     = new TH2F("nOniaMuPtVsPt","DoubleMu_Er_HighQ_WdEta22 (Quarkonia); p_{T} cut mu_{1}; p_{T} cut mu_{2}",41,-0.25,20.25,41,-0.25,20.25);

  //Sums
  hRate1F["nHTTVsHTT"] = new TH1F("nHTTVsHTT","HTT; HTT cut; rate [Hz]",512,-.5,511.5);
  hRate1F["nETTVsETT"] = new TH1F("nETTVsETT","ETT; ETT cut; rate [Hz]",512,-.5,511.5);
  hRate1F["nETMVsETM"] = new TH1F("nETMVsETM","ETM; ETM cut; rate [Hz]",512,-.5,511.5);

  return true;
}       // -----  end of function L1Plot::BookRateHistogram  -----

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
    f.second->Scale(scale);
    f.second->Write();
  }
  
  for(auto f : hRate2F)
  {
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

    //if(dijetPt2>=ptCut){
    //hTH1F["nDiJetVsPt"]->Fill(ptCut);

    //for(int ptCut_0=ptCut; ptCut_0<256; ++ptCut_0) {
    //if(dijetPt1>=ptCut_0) hTH2F["nAsymDiJetVsPt"]->Fill(ptCut_0,ptCut);
    //}
    //}

    //if(diCenjetPt2>=ptCut){
    //hTH1F["nDiCenJetVsPt"]->Fill(ptCut);

    //for(int ptCut_0=ptCut; ptCut_0<256; ++ptCut_0) {
    //if(diCenjetPt1>=ptCut_0) hTH2F["nAsymDiCenJetVsPt"]->Fill(ptCut_0,ptCut);
    //}
    //}

    //if(ditauPt>=ptCut)    hTH1F["nDiTauVsPt"]->Fill(ptCut);
    //if(quadjetPt>=ptCut)  hTH1F["nQuadJetVsPt"]->Fill(ptCut);
    //if(quadjetCPt>=ptCut) hTH1F["nQuadCenJetVsPt"]->Fill(ptCut);

  } //loop on 256




  for(int ptCut=0; ptCut<65; ++ptCut) {
    if(L1Event->EGPt>=ptCut)    hRate1F["nEGVsPt"]->Fill(ptCut);
    if(L1Event->EGerPt>=ptCut)  hRate1F["nEGErVsPt"]->Fill(ptCut);
    if(L1Event->IsoEGerPt>=ptCut) hRate1F["nIsoEGVsPt"]->Fill(ptCut);
    // 
    //   if(diEG2>=ptCut)     hTH1F["nDiEGVsPt"]->Fill(ptCut);
    //   if(diIsolEG2>=ptCut) hTH1F["nDiIsoEGVsPt"]->Fill(ptCut);
    // 
    // 
    //   for(int ptCut2=0; ptCut2<=65; ++ptCut2) {
    // 	if(diEG1>=ptCut && diEG2>=ptCut2 && ptCut2 <= ptCut) hTH2F["nEGPtVsPt"]->Fill(ptCut,ptCut2);
    // 	if(diIsolEG1>=ptCut && diIsolEG2>=ptCut2 && ptCut2<= ptCut) hTH2F["nIsoEGPtVsPt"]->Fill(ptCut,ptCut2);
    //   }
    // 
  }//loop on 65
  //  
  // for(int ptCut=0; ptCut<131; ++ptCut) {
  //   if (muPt>=ptCut)    hTH1F["nMuVsPt"]->Fill(ptCut);
  //  if (muErPt>=ptCut)  hTH1F["nMuErVsPt"]->Fill(ptCut);
  // }
  //   
  // 
  // for(int iCut=0; iCut<41; ++iCut) {
  //   for(int iCut2=0; iCut2<=iCut; ++iCut2) {
  // 	float ptCut = iCut*0.5;
  // 	float ptCut2 = iCut2*0.5;
  // 	if (doubleMuPt1>=ptCut && doubleMuPt2>=ptCut2) hTH2F["nMuPtVsPt"]->Fill(ptCut,ptCut2);
  // 	if (oniaMuPt1>=ptCut && oniaMuPt2>=ptCut2)     hTH2F["nOniaMuPtVsPt"]->Fill(ptCut,ptCut2);
  //   }
  // }
  // 
  for(int httCut=0; httCut<512; ++httCut) {
    if(L1Event->HTT>httCut) hRate1F["nHTTVsHTT"]->Fill(httCut);
    if(L1Event->ETT>httCut) hRate1F["nETTVsETT"]->Fill(httCut);
    if(L1Event->ETM>httCut) hRate1F["nETMVsETM"]->Fill(httCut);
  }

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
    GetRecoEvent();
    FillEffHistogram();
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
//         Name:  L1Plot::GetRecoSum
//  Description:  
// ===========================================================================
std::vector<TLorentzVector> L1Plot::GetRecoSum(std::string type ) const
{
  std::vector<TLorentzVector> reTLVs;
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
std::vector<TLorentzVector> L1Plot::GetRecoMuon(bool isER, float IsoCut, int qual) const
{
  std::vector<TLorentzVector> reTLVs;
  const float MuERcut = 2.1;
  if (recoMuon_ == NULL) return reTLVs;

  for (int i = 0; i < recoMuon_->nMuons; ++i)
  {
    if (qual == 1 && ! recoMuon_ -> isLooseMuon.at(i)) continue;
    if (qual == 2 && ! recoMuon_ -> isMediumMuon.at(i)) continue;

    if (isER && fabs(recoMuon_->eta.at(i)) > MuERcut )
      continue;

    if (recoMuon_->iso.at(i) < IsoCut) continue;

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
  recoEvent["Tauer"]   = GetRecoTau(true);
  recoEvent["EG"]      = GetRecoEle();
  recoEvent["EGer"]    = GetRecoEle(true,   false, 0 );
  recoEvent["IsoEG"]   = GetRecoEle(false,  true,  0 );
  recoEvent["IsoEGer"] = GetRecoEle(true,   true,  0 );
  recoEvent["Mu"]      = GetRecoMuon();
  recoEvent["MuOpen"]  = GetRecoMuon();
  recoEvent["Muer"]    = GetRecoMuon(true);
  recoEvent["HTT"]     = GetRecoSum("HTT");
  recoEvent["ETM"]     = GetRecoSum("ETM");
  recoEvent["ETT"]     = GetRecoSum("ETT");
  recoEvent["HTM"]     = GetRecoSum("HTM");
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
        if (objstr.find("ETM") != std::string::npos )
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
void L1Plot::SetTodo (bool doPlotRate_, bool doPlotEff_)
{
  doPlotRate = doPlotRate_;
  doPlotEff = doPlotEff_;
}       // -----  end of function L1Plot::SetTodo  -----

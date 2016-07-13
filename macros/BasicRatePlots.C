#include "L1Ntuple.h"
#include "L1AlgoFactory.h"
#include <algorithm>
#include<map>
#include<iostream>

#include "TH1F.h"
#include "TH2F.h"

class BasicRatePlots : public L1AlgoFactory
{
public :
  
  //constructor    
  BasicRatePlots(std::string filename){
    if (filename.find(".root") != std::string::npos) {
      std::cout << "Reading RootFile: " << filename << std::endl;
      L1Ntuple::Open(filename);
    }else{
      std::cout << "Reading Filelist: " << filename << std::endl;
      if (! L1Ntuple::OpenWithList(filename)) exit(0);
    }
  }
  ~BasicRatePlots() {}
  
  void run(bool runOnData, std::string resultTag, float crossSec, int nBunches, int isCrossSec, int nEvents = 0, bool isReweight = false);
  
private :
  
  float ScaleFactor(float nZeroBias, float nBunches);
  void  SingleJetPt(Float_t& ptcut, Bool_t isCentral = false);
  float SingleMuEta(float ptcut);
  void  SingleMuPt(Float_t& ptcut, Bool_t isER, float ErMu);
  float SecondMuEta(float ptcut);
  void  SingleMuEtaRate(Float_t& etaCut, float ptCut, bool doIso);
  void  SingleEGPt(Float_t& ptcut, Bool_t isIsolated, Bool_t isER, float Er);
  float SingleEGEta(float ptCut, bool doIso);
  float SingleEGIEta(float ptCut, bool doIso);
  void  SingleEGEtaRate(Float_t& etaCut, float ptCut, bool doIso);
  float SingleJetEta(float pt, Int_t accept_flag = 0);

  void setRateError(TH1F* histo);
  
  void ETMVal(Float_t& ETMcut);
  void HTTVal(Float_t& HTTcut);
  void HTMVal(Float_t& HTMcut);
  void ETTVal(Float_t& ETTcut);
  
  std::map<std::string,TH1F*> hTH1F;
  std::map<std::string,TH2F*> hTH2F;
  std::map<std::string,TH1F*> hCountTH1F;
};

// ------------------------------------------------------------------
// BasicRatePlots::BasicRatePlots(TTree *tree) : L1AlgoFactory(tree) {}

// scale factor computed w.r.t. ZeroBias rate fratcion and # bunches 
float BasicRatePlots::ScaleFactor(float nZeroBias, float nBunches) {

  float scal = 11246.; // ZB per bunch in Hz
  scal /= nZeroBias;
  scal *= nBunches;

  return scal;
}

float BasicRatePlots::SingleMuEta(float ptCut ) {

  float maxPt = -10;
  float iMuMaxPt = -10;

  UInt_t nMuons = upgrade_ -> nMuons;
  for (UInt_t imu=0; imu < nMuons; imu++) 
    {
    Int_t bx = upgrade_ -> muonBx.at(imu);
    if(bx != 0) continue;
    if((upgrade_->muonQual.at(imu)==4) or (upgrade_->muonQual.at(imu)==8)) continue;
    Float_t pt = upgrade_ -> muonEt.at(imu);                       
      if ( pt > maxPt) 
	{
	  maxPt = pt;
	  iMuMaxPt = imu;
	}
    }
  
  float eta = maxPt>ptCut ? upgrade_ -> muonEta.at(iMuMaxPt) : -999.; 
  return eta;
}

float BasicRatePlots::SecondMuEta(float ptCut) {

  Float_t mu1ptmax = -10.;
  Float_t mu2ptmax = -10.;
  float eta1 = -10.;
  float eta2 = -10.;
  float aux = 0.;

  if(upgrade_->nMuons < 2) return -10.;

  for(UInt_t um=0; um < upgrade_ -> nMuons; um++) {
    Int_t bx = upgrade_ -> muonBx.at(um);  
    if(bx != 0) continue;
    if(upgrade_->muonQual.at(um)==4) continue;

    Float_t pt = upgrade_ -> muonEt.at(um);    // the rank of the electron
    if(pt >= mu1ptmax) {
      mu2ptmax = mu1ptmax;
      mu1ptmax = pt;
    }
    else if(pt >= mu2ptmax) {
      mu2ptmax = pt;
    }
  }

  if(mu2ptmax < ptCut) return -10.;

  for(UInt_t um=0; um < upgrade_ -> nMuons; um++) {
    if (upgrade_->muonEt.at(um)==mu1ptmax) {
      eta1=upgrade_ -> muonEta.at(um);
      aux = um;
      break;
    }
  }

  //std::cout<<"aux: "<<std::to_string(aux)<<std::endl;

  for(UInt_t um=0; um < upgrade_ -> nMuons; um++) {
    if (upgrade_->muonEt.at(um)==mu2ptmax && um != aux) {
      eta2=upgrade_ -> muonEta.at(um);
      break;
    }
  }

  //std::cout<<"pt1: "<<std::to_string(mu1ptmax)<<std::endl;
  //std::cout<<"pt2: "<<std::to_string(mu2ptmax)<<std::endl;
  //std::cout<<"eta1: "<<std::to_string(eta1)<<std::endl;
  //std::cout<<"eta2: "<<std::to_string(eta2)<<std::endl;
  
  return eta2;
}

void BasicRatePlots::SingleMuPt(Float_t& cut, Bool_t isER, float ErMu) {

  Float_t ptmax = -10.;

  for(UInt_t um=0; um < upgrade_ -> nMuons; um++) {
    Int_t bx = upgrade_ -> muonBx.at(um);  
    if(bx != 0) continue;
    Float_t eta = upgrade_ -> muonEta.at(um);
    if(fabs(eta) > ErMu && isER) continue;
    if((upgrade_->muonQual.at(um)==4) or (upgrade_->muonQual.at(um)==8)) continue;

    Float_t pt = upgrade_ -> muonEt.at(um);    // the rank of the electron
    if(pt >= ptmax) ptmax = pt;
  }

  cut = ptmax;
  return;
}

void BasicRatePlots::SingleMuEtaRate(Float_t& etaCut, float ptCut, bool doIso) {

  Float_t etamax = 10.;

  for(UInt_t um=0; um < upgrade_ -> nMuons; um++) {
    Int_t bx = upgrade_ -> muonBx.at(um);  
    if(bx != 0) continue;
    //if(upgrade_->muonQual.at(um)==0) continue;
    //if(doIso && !(upgrade_ -> egIso.at(um))) continue;
    if((upgrade_->muonQual.at(um)==4) or (upgrade_->muonQual.at(um)==8)) continue;
    Float_t pt = upgrade_ -> muonEt.at(um);
    if(pt < ptCut) continue;
	
    Float_t eta = fabs(upgrade_ -> muonEta.at(um));
    if(eta <= etamax) etamax = eta;
  }

  etaCut = etamax;
  return;
}

float BasicRatePlots::SingleEGEta(float ptCut, bool doIso) {

  float maxPt = -10;
  float iEGMaxPt = -10;

  for (UInt_t ue=0; ue < upgrade_ -> nEGs; ue++) {
    Int_t bx = upgrade_ -> egBx.at(ue);        		
    if (bx != 0) continue;
    Bool_t iso = upgrade_ -> egIso.at(ue);
    if (!iso && doIso) continue;
    Float_t pt = upgrade_ -> egEt.at(ue);
    if ( pt >= maxPt) 
      {
	maxPt = pt;
	iEGMaxPt = ue;
      }
  }

  return iEGMaxPt>=0 && maxPt>ptCut ? upgrade_ -> egEta.at(iEGMaxPt) : -10.; 
}

float BasicRatePlots::SingleEGIEta(float ptCut, bool doIso) {

  float maxPt = -10;
  float iEGMaxPt = -999;

  for (UInt_t ue=0; ue < upgrade_ -> nEGs; ue++) {
    Int_t bx = upgrade_ -> egBx.at(ue);        		
    if (bx != 0) continue;
    Bool_t iso = upgrade_ -> egIso.at(ue);
    if (!iso && doIso) continue;
    Float_t pt = upgrade_ -> egEt.at(ue);
    if ( pt >= maxPt) 
      {
	maxPt = pt;
	iEGMaxPt = ue;
      }
  }

  return iEGMaxPt>=0 && maxPt>ptCut ? upgrade_ -> egIEta.at(iEGMaxPt) : -999.; 
}

void BasicRatePlots::SingleEGPt(Float_t& cut, Bool_t isIsolated , Bool_t isER, float Er) {

  Float_t ptmax = -10.;

  for(UInt_t ue=0; ue < upgrade_ -> nEGs; ue++) {
    Int_t bx = upgrade_ -> egBx.at(ue);  
    if(bx != 0) continue;
    if(isIsolated && !(upgrade_ -> egIso.at(ue))) continue;
    Float_t eta = upgrade_ -> egEta.at(ue);
    if(fabs(eta) > Er && isER) continue;  // eta = 5 - 16

    Float_t pt = upgrade_ -> egEt.at(ue);    // the rank of the electron
    if(pt >= ptmax) ptmax = pt;
  }

  cut = ptmax;
  return;
}

void BasicRatePlots::SingleEGEtaRate(Float_t& etaCut, float ptCut, bool doIso) {

  Float_t etamax = 10.;
  
  for(UInt_t ue=0; ue<upgrade_ -> nEGs; ue++) {
    Int_t bx = upgrade_ -> egBx.at(ue);
    if (bx !=0) continue;
    if(doIso && !(upgrade_ -> egIso.at(ue))) continue;
    Float_t pt = upgrade_ -> egEt.at(ue);
    if(pt < ptCut) continue;
	
    Float_t eta = fabs(upgrade_ -> egEta.at(ue));
    if(eta <= etamax) etamax = eta;
  }
  
  etaCut = etamax;
  return;
}

void BasicRatePlots::SingleJetPt(Float_t& cut, Bool_t isCentral) {

  Float_t ptmax = -10.;
  Int_t Nj = upgrade_ -> nJets ;
  for(Int_t ue=0; ue < Nj; ue++) {
    Int_t bx = upgrade_ -> jetBx.at(ue);
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(upgrade_ -> jetEta.at(ue)) > 3. ? true : false;
    if(isCentral && isFwdJet) continue;
    //if(NOTauInJets && upgrade_->Taujet[ue]) continue;
    //if(isCentral && noHF && (upgrade_->jetEta.at(ue) < 5 || upgrade_->jetEta.at(ue) > 17)) continue;

    Float_t pt = upgrade_ -> jetEt.at(ue);
    if(pt >= ptmax) ptmax = pt;
  }

  cut = ptmax;
  return;
}

float BasicRatePlots::SingleJetEta(float ptCut, Int_t accept_flag) {

  float maxPt = -10;
  float iJetMaxPt = -10;
  
  for(UInt_t ue=0; ue < upgrade_ -> nJets; ue++) {
    Int_t bx = upgrade_ -> jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(upgrade_ -> jetEta.at(ue)) > 3. ? true : false;

    if(accept_flag == 1 && isFwdJet) continue;
    if(accept_flag == 2 && !isFwdJet) continue;

    Float_t pt = upgrade_ -> jetEt.at(ue);
    if(pt >= maxPt){
      maxPt = pt;
      iJetMaxPt = ue;
    }
  }

  return iJetMaxPt>=0 && maxPt>ptCut ? upgrade_ -> jetEta.at(iJetMaxPt) : -10.;
}

void BasicRatePlots::ETMVal(Float_t& ETMcut ) {
  
  Float_t TheETM = -10;
  int idx = GetSumEtIdx(EtSumType::ETM);
  assert(upgrade_->sumType.at(idx) == EtSumType::ETM);
  //std::cout<<"ETM: "<<std::to_string(idx)<<std::endl;
  if(upgrade_->sumBx.at(idx)==0) TheETM =upgrade_->sumEt.at(idx);
  ETMcut = TheETM;
  return;
}

void BasicRatePlots::HTTVal(Float_t& HTTcut) {

  Float_t TheHTT = -10;
  int idx= GetSumEtIdx(EtSumType::HTT);
  assert(upgrade_->sumType.at(idx) == EtSumType::HTT);
  if(upgrade_->sumBx.at(idx)==0) TheHTT =upgrade_->sumEt.at(idx);
  HTTcut = TheHTT;
  return;
}

void BasicRatePlots::HTMVal(Float_t& HTMcut) {

  Float_t TheHTM = -10;
  int idx= GetSumEtIdx(EtSumType::HTM);
  assert(upgrade_->sumType.at(idx) == EtSumType::HTM);
  if(upgrade_->sumBx.at(idx)==0) TheHTM =upgrade_->sumEt.at(idx);
  HTMcut = TheHTM;
  return;
}

void BasicRatePlots::ETTVal(Float_t& ETTcut) {

  Float_t TheETT = -10;
  int idx= GetSumEtIdx(EtSumType::ETT);
  assert(upgrade_->sumType.at(idx) == EtSumType::ETT);
  if(upgrade_->sumBx.at(idx)==0) TheETT =upgrade_->sumEt.at(idx);
  ETTcut = TheETT;
  return;
}

void BasicRatePlots::setRateError(TH1F* histo) {

  int nBins = histo->GetNbinsX();

  for (int iBin=1; iBin<=nBins; ++iBin) {
    float value = histo->GetBinContent(iBin);
    float error = sqrt(value);

    histo->SetBinError(iBin,error);
  }  
}

// --------------------------------------------------------------------
//                             run function 
// --------------------------------------------------------------------


void BasicRatePlots::run(bool runOnData, std::string resultTag, float crossSec, int nBunches, int isCrossSec, int nEvents, bool isReweight) {

  system("mkdir -p results");
  std::string resultName = "/afs/cern.ch/user/z/zhangj/test/L1Menu/80x_stage2/v62p3_MC2016/CMSSW_8_0_9/src/L1TriggerDPG/L1Menu/macros/results/results_" + resultTag + (isCrossSec ? "_XSEC" : "_RATE") + ".root";
  TFile *outFile = new TFile((resultName).c_str(),"recreate");
  outFile->cd();

  float triggerTowerBinsWithHF[65]={-5.,-4.5,-4.,-3.5,-3.,-2.65,-2.5,-2.322,-2.172,-2.043,-1.93,-1.83,-1.74,-1.653,-1.566,-1.479,-1.392,-1.305,-1.218,-1.131,-1.044,-0.957,-0.87,-0.783,-0.696,-0.609,-0.522,-0.435,-0.348,-0.261,-0.174,-0.087,0.,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696,0.783,0.87,0.957,1.044,1.131,1.218,1.305,1.392,1.479,1.566,1.653,1.74,1.83,1.93,2.043,2.172,2.322,2.5,2.65,3.,3.5,4.,4.5,5.};
  float triggerTowerBins[57]={-3.,-2.65,-2.5,-2.322,-2.172,-2.043,-1.93,-1.83,-1.74,-1.653,-1.566,-1.479,-1.392,-1.305,-1.218,-1.131,-1.044,-0.957,-0.87,-0.783,-0.696,-0.609,-0.522,-0.435,-0.348,-0.261,-0.174,-0.087,0.,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696,0.783,0.87,0.957,1.044,1.131,1.218,1.305,1.392,1.479,1.566,1.653,1.74,1.83,1.93,2.043,2.172,2.322,2.5,2.65,3.};
  float absTriggerTowerBins[29]={0.,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696,0.783,0.87,0.957,1.044,1.131,1.218,1.305,1.392,1.479,1.566,1.653,1.74,1.83,1.93,2.043,2.172,2.322,2.5,2.65,3.};

  //Event Counter
  hTH1F["nEvts"]       = new TH1F("nEvts","Number of Events Processed",1,-0.5,.5);
  //Single stuff
  hTH1F["nJetVsPt"]    = new TH1F("nJetVsPt","SingleJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hTH1F["nJetCenVsPt"] = new TH1F("nJetCenVsPt","SingleJetCentral; E_{T} cut; rate [Hz]",256,-0.5,255.5);

  hTH1F["nTauVsPt"]     = new TH1F("nTauVsPt","SingleTau; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH1F["nTauErVsPt"]   = new TH1F("nTauErVsPt","SingleTauer; E_{T} cut; rate [Hz]",129,-0.5,128.5);
  hTH1F["nIsoTauVsPt"]  = new TH1F("nIsoTauVsPt","SingleIsoTau; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH1F["nIsoTauErVsPt"]  = new TH1F("nIsoTauErVsPt","SingleIsoTauEr; E_{T} cut; rate [Hz]",65,-0.5,64.5);

  hTH1F["nEGVsPt"]     = new TH1F("nEGVsPt", "SingleEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH1F["nEGErVsPt"]   = new TH1F("nEGErVsPt", "SingleEGEr; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH1F["nIsoEGVsPt"]  = new TH1F("nIsoEGVsPt","SingleIsoEGer; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH1F["nMuVsPt"]     = new TH1F("nMuVsPt","SingleMu; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hTH1F["nMuErVsPt"]   = new TH1F("nMuErVsPt","SingleMu |#eta|<2.1; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hTH1F["nMuVsEta"]    = new TH1F("nMuVsEta","nMuVsEta",24,-2.4,2.4);
  hTH1F["nEGVsEta"]    = new TH1F("nEGVsEta","nEGVsEta",56,triggerTowerBins);
  hTH1F["nIsoEGVsEta"] = new TH1F("nIsoEGVsEta","nIsoEGVsEta",56,triggerTowerBins);
  
  hCountTH1F["nMuCountVsEta"]    = new TH1F("nMuCountVsEta","nMuVsEta",24,-2.4,2.4);
  hCountTH1F["nEGCountVsEta"]    = new TH1F("nEGCountVsEta","nEGVsEta",56,triggerTowerBins);
  hCountTH1F["nIsoEGCountVsEta"] = new TH1F("nIsoEGCountVsEta","nIsoEGVsEta",56,triggerTowerBins);

  //--------Rate w.r.t. different eta--------//
  float ER[5] = {1.5,1.8,2.1,2.5,3.};
  float ERMu[3] = {0.8,1.2,2.5};
  for(int i=0; i<5; i++){
    hTH1F["nEGEr"+std::to_string(ER[i])+"VsPt"] = new TH1F(("nEGEr"+std::to_string(ER[i])+"VsPt").c_str(),("SingleEGEr"+std::to_string(ER[i])+"; E_{T} cut; rate [Hz]").c_str(),65,-0.5,64.5);
    hTH1F["nIsoEGEr"+std::to_string(ER[i])+"VsPt"] = new TH1F(("nIsoEGEr"+std::to_string(ER[i])+"VsPt").c_str(),("SingleIsoEGEr"+std::to_string(ER[i])+"; E_{T} cut; rate [Hz]").c_str(),65,-0.5,64.5);
  }
  for(int j=0; j<3; j++){
    hTH1F["nMuEr"+std::to_string(ERMu[j])+"VsPt"] = new TH1F(("nMuEr"+std::to_string(ERMu[j])+"VsPt").c_str(),("SingleMuonEr"+std::to_string(ERMu[j])+"; E_{T} cut; rate [Hz]").c_str(),131,-0.5,130.5);
  }

  //--------Eta Rate w.r.t. different Pt--------//
  float PT[6] = {10.,15.,20.,25.,30.,35.};
  float PTMu[4] = {14.,16.,18.,20.,};
  for(unsigned i=0; i<6; i++) {
    hTH1F["nEGEr"+std::to_string(PT[i])+"VsEta"] = new TH1F(("nEGEr"+std::to_string(PT[i])+"VsEta").c_str(),("SingleEGEr"+std::to_string(PT[i])+"; #eta restriction cut; rate [Hz]").c_str(),28,absTriggerTowerBins);
    hTH1F["nIsoEGEr"+std::to_string(PT[i])+"VsEta"] = new TH1F(("nIsoEGEr"+std::to_string(PT[i])+"VsEta").c_str(),("SingleIsoEGEr"+std::to_string(PT[i])+"; #eta restriction cut; rate [Hz]").c_str(),28,absTriggerTowerBins);
  }
  for(unsigned j=0; j<4; j++) {
    hTH1F["nMuEr"+std::to_string(PTMu[j])+"VsEta"] = new TH1F(("nMuEr"+std::to_string(PTMu[j])+"VsEta").c_str(),("SingleMuEr"+std::to_string(PTMu[j])+"; #eta restriction cut; rate [Hz]").c_str(),28,-0.05,2.55);
  }
  
  hTH1F["nJetVsEta"]   = new TH1F("nJetVsEta","nJetVsEta",64,triggerTowerBinsWithHF);
  hTH1F["nJetVsEta_Central"] = new TH1F("nJetVsEta_Central","nJetVsEta_Central",50,-5.,5.);
  hTH1F["nJetVsEta_Fwd"]     = new TH1F("nJetVsEta_Fwd","nJetVsEta_Fwd",50,-5.,5.);

  //Multistuff
  hTH1F["nDiJetVsPt"]        = new TH1F("nDiJetVsPt","DiJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hTH1F["nDiCenJetVsPt"]     = new TH1F("nDiCenJetVsPt","DiCenJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hTH2F["nAsymDiJetVsPt"]    = new TH2F("nAsymDiJetVsPt","DiJet; E_{T} cut jet 1; E_{T} cut jet 2",256,-0.5,255.5,256,-0.5,255.5);
  hTH2F["nAsymDiCenJetVsPt"] = new TH2F("nAsymDiCenJetVsPt","DiCenJet; E_{T} cut jet 1; E_{T} cut jet 2",256,-0.5,255.5,256,-0.5,255.5);
  hTH1F["nQuadJetVsPt"]      = new TH1F("nQuadJetVsPt","QuadJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hTH1F["nQuadCenJetVsPt"]   = new TH1F("nQuadCenJetVsPt","QuadCenJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hTH1F["nDiTauVsPt"]        = new TH1F("nDiTauVsPt","DiTau; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hTH1F["nDiIsoTauErVsPt"]     = new TH1F("nDiIsoTauErVsPt","DiIsoTauEr; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  
  hTH1F["nDiEGVsPt"]         = new TH1F("nDiEGVsPt","DiEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH1F["nDiIsoEGVsPt"]      = new TH1F("nDiIsoEGVsPt","DiIsoEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH2F["nEGPtVsPt"]         = new TH2F("nEGPtVsPt","DoubleEle; p_{T} cut EG_{1}; p_{T} cut EG_{2}",65,-0.5,64.5,65,-0.5,64.5);
  hTH2F["nIsoEGPtVsPt"]      = new TH2F("nIsoEGPtVsPt","DoubleIsolEle; p_{T} cut EG_{1}; p_{T} cut EG_{2}",65,-0.5,64.5,65,-0.5,64.5);
  hTH2F["nMuPtVsPt"]         = new TH2F("nMuPtVsPt","DoubleMu; p_{T} cut mu_{1}; p_{T} cut mu_{2}",41,-0.25,20.25,41,-0.25,20.25);
  hTH2F["nOniaMuPtVsPt"]     = new TH2F("nOniaMuPtVsPt","DoubleMu_Er_HighQ_WdEta22 (Quarkonia); p_{T} cut mu_{1}; p_{T} cut mu_{2}",41,-0.25,20.25,41,-0.25,20.25);
  
  hTH1F["nDiMuVsPt"]         = new TH1F("nDiMuVsPt","DoubleMu; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hTH1F["nSecondMuVsEta"]    = new TH1F("nSecondMuVsEta","nSecondMuVsEta",24,-2.4,2.4);

  //Sums
  hTH1F["nHTTVsHTT"] = new TH1F("nHTTVsHTT","HTT; HTT cut; rate [Hz]",512,-.5,511.5);
  hTH1F["nETTVsETT"] = new TH1F("nETTVsETT","ETT; ETT cut; rate [Hz]",512,-.5,511.5);
  hTH1F["nETMVsETM"] = new TH1F("nETMVsETM","ETM; ETM cut; rate [Hz]",512,-.5,511.5);

  TH1F *egIEta = new TH1F("nEGVsIEta","nEGVsIEta",229,-114.5,114.5);
  TH1F *isoEgIEta = new TH1F("nIsoEGVsIEta","nIsoEGVsIEta",229,-114.5,114.5);

  egIEta->Sumw2();
  isoEgIEta->Sumw2();

  if (isCrossSec) {
    std::map<std::string,TH1F*>::iterator hTH1FIt  = hTH1F.begin();
    std::map<std::string,TH1F*>::iterator hTH1FEnd = hTH1F.end();

    for(; hTH1FIt!=hTH1FEnd; ++hTH1FIt)
      {
	hTH1FIt->second->GetYaxis()->SetTitle("cross section [#mubarn]");
      }
  }

  //outFile->Close();
  Double_t nZeroBias = 0.;

  int nLumi(0);

  if (nEvents <= 0){
    nEvents=fChain->GetEntriesFast();
  }

  std::cout << "Tree contains " << fChain->GetEntriesFast() << " events." << std::endl;
  std::cout << "Running on " << nEvents << " events." << std::endl;

  for (Long64_t event=0; event<nEvents; ++event) {
    Long64_t eventEntry = LoadTree(event);
    if (eventEntry < 0) break;
    GetEntry(event);
    if (event%200000 == 0) {
      if (event_!=NULL)
        std::cout << "Processed " << event << " events. Current run number: " << event_ -> run << " lumi: " << event_ -> lumi << std::endl;
      else
        std::cout << "Processed " << event << " events." << std::endl;
    }

    //if (event_->nPV_True!=19 and event_->nPV_True!=18 and event_->nPV_True!=17 and event_->nPV_True!=19 and event_->nPV_True!=18) continue;
    if (event_->nPV_True!=29 and event_->nPV_True!=28 and event_->nPV_True!=27) continue;
    
    hTH1F["nEvts"]->Fill(0.);  // count number of events processed
    
    float jetPt     = 0.; SingleJetPt(jetPt);
    float jetCenPt  = 0.; SingleJetPt(jetCenPt,true);
    float jetEta    = SingleJetEta(60.);
    float jetEta_Central = SingleJetEta(60.,1);
    float jetEta_Fwd    = SingleJetEta(60.,2);

    float TauPt      = 0.;
    float TauErPt    = 0.;
    float isoTauPt   = 0.;
    float isoTauErPt = 0.;
    SingleTauPt(TauPt      , false, false);
    SingleTauPt(TauErPt    , true,  false);
    SingleTauPt(isoTauPt   , false, true);
    SingleTauPt(isoTauErPt , true,  true);


    float htt       = 0.; HTTVal(htt);
    float ett       = 0.; ETTVal(ett);
    float etm       = 0.; ETMVal(etm);
    float egPt      = 0.; SingleEGPt(egPt,false,false, 3.);
    float egErPt    = 0.; SingleEGPt(egErPt,false,true, 2.1);
    //---fixed this---//
    float isoEgPt   = 0.; SingleEGPt(isoEgPt,true,false, 3.);
    
    //--------Rate w.r.t. different eta--------//
    float egErVPt[5] = {0., 0., 0., 0., 0.,};
    float isoEgErVPt[5] = {0., 0., 0., 0., 0.,};
    for(int i=0; i<5; i++){
      SingleEGPt(egErVPt[i],false,true, ER[i]);
      SingleEGPt(isoEgErVPt[i],true,true, ER[i]);
    }
    
    //--------Eta Rate w.r.t. different Pt--------//
    float egErEta[6] = {0.,0.,0.,0.,0.,0.};
    float isoEgErEta[6] = {0.,0.,0.,0.,0.,0.};
    for(unsigned i=0; i<6; i++){
      SingleEGEtaRate(egErEta[i],PT[i],false);
    }
    for(unsigned i=0; i<6; i++){
      SingleEGEtaRate(isoEgErEta[i],PT[i],true);
    }
    
    float egEta     = SingleEGEta(20.,false);
    float isoegEta  = SingleEGEta(20.,true);

    float egIeta     = SingleEGIEta(20.,false);
    float isoegIeta  = SingleEGIEta(20.,true);

    float muPt     = -10.; SingleMuPt(muPt,false,2.5);
    float muErPt   = -10.; SingleMuPt(muErPt,true,2.1);

    //---Muon Eta Tricks---//
    float muErVPt[3] = {-10., -10., -10.};
    for(int j=0; j<3; j++){
      SingleMuPt(muErVPt[j], true, ERMu[j]);
    }
    float muErVEta[4] = {0.,0.,0.,0.};
    for(unsigned j=0; j<4; j++){
      SingleMuEtaRate(muErVEta[j],PTMu[j],false);
    }
    
    float muEta    = SingleMuEta(16.);
    float secondMuEta = SecondMuEta(0.);
    float doubleMuPt1 = -10.;
    float doubleMuPt2 = -10.;
    DoubleMuPt(doubleMuPt1,doubleMuPt2,true,false);

    float diMu1    = -10.;
    float diMu2    = -10.;
    DoubleMuPt(diMu1,diMu2,true,false);    
     
    float oniaMuPt1 = 0.;
    float oniaMuPt2 = 0.;
    Onia2015Pt(oniaMuPt1,oniaMuPt2,true, false, 22);
     
    float dijetPt1    = -10.;
    float dijetPt2    = -10.;
    float diCenjetPt1 = -10.;
    float diCenjetPt2 = -10.;
    DoubleJetPt(dijetPt1,dijetPt2);
    DoubleJetPt(diCenjetPt1,diCenjetPt2,true);
    Float_t dummy = -1;
    float diTauErPt    = -10.; DoubleTauJetEta2p17Pt(dummy,diTauErPt);
    float diIsoTauEr = -10.; DoubleTauJetEta2p17Pt(dummy,diIsoTauEr,true);
    dummy = -1.;
    float quadjetPt  = -10.; QuadJetPt(dummy,dummy,dummy,quadjetPt);
    dummy = -1.;
    float quadjetCPt = -10.; QuadJetPt(dummy,dummy,dummy,quadjetCPt,true);
    dummy = -1.;
    // 
    float diEG1     = -10.;
    float diEG2     = -10.;
    float diIsolEG1 = -10.;
    float diIsolEG2 = -10.;
    DoubleEGPt(diEG1,diEG2,false);
    DoubleEGPt(diIsolEG1,diIsolEG2,true);

    //outFile->Open(resultName.c_str(),"update");

    // pile up value from 0 to 50,
    float puWeight[51] = {1638.538183856806, 4.693674415610433, 1.5300909160502814, 0.7182406823718569, 0.4461270187978222, 0.39812703022422885, 0.4432122086751145, 0.5742601670647391, 0.7739373760577931, 1.0025956642167804, 1.3122680368593147, 1.605655051145846, 1.8929024225833195, 2.1284381102257535, 2.2556165182260717, 2.244482200364411, 2.154571520396462, 1.9528537819111809, 1.695897273333598, 1.4282923583611258, 1.1369113035284417, 0.8862206448254017, 0.6536240971561117, 0.46856364438557374, 0.33179561861669243, 0.22382832523838006, 0.14730363368878097, 0.09882046047249428, 0.060769117256921634, 0.04258674568273182, 0.024471799847250696, 0.017208722825574143, 0.009074708975294003, 0.006641952074695454, 0.0033690053530415084, 0.0020709022278971793, 0.0018436379632363492, 0.0009373317263744043, 0.0, 0.0007986637595886011, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0};

    float weight = 1;
    
    if (isReweight){
      int PU = recoVtx_->nVtx;
      weight = puWeight[PU];
      //weight = 1;
      //std::cout<<std::to_string(weight)<<std::endl;
    } else {
      weight = 1;
    }

    //std::cout<<std::to_string(event_->nPV_True)<<std::endl;
    
    nZeroBias+=weight;
    //std::cout<<std::to_string(weight)<<std::endl;
    //std::cout<<std::to_string(nZeroBias)<<std::endl;

    hTH1F["nMuVsEta"]->Fill(muEta,weight);
    hTH1F["nSecondMuVsEta"]->Fill(secondMuEta,weight);
    hTH1F["nEGVsEta"]->Fill(egEta,weight);
    hTH1F["nIsoEGVsEta"]->Fill(isoegEta,weight);
    hCountTH1F["nMuCountVsEta"]->Fill(muEta,weight);
    hCountTH1F["nEGCountVsEta"]->Fill(egEta,weight);
    hCountTH1F["nIsoEGCountVsEta"]->Fill(isoegEta,weight);
    
    egIEta->Fill(egIeta,weight);
    isoEgIEta->Fill(isoegIeta,weight);
    
    hTH1F["nJetVsEta"]->Fill(jetEta,weight);
    hTH1F["nJetVsEta_Central"]->Fill(jetEta_Central,weight);
    hTH1F["nJetVsEta_Fwd"]->Fill(jetEta_Fwd,weight);

    for(int ptCut=0; ptCut<256; ++ptCut) {
      if(jetPt>=ptCut)	  hTH1F["nJetVsPt"]->Fill(ptCut,weight);
      if(jetCenPt>=ptCut) hTH1F["nJetCenVsPt"]->Fill(ptCut,weight);
      if(TauPt>=ptCut)	  hTH1F["nTauVsPt"]->Fill(ptCut,weight);
      if(TauErPt>=ptCut)	  hTH1F["nTauErVsPt"]->Fill(ptCut,weight);
      if(isoTauPt>=ptCut)	  hTH1F["nIsoTauVsPt"]->Fill(ptCut,weight);
      if(isoTauErPt>=ptCut)	  hTH1F["nIsoTauErVsPt"]->Fill(ptCut,weight);
     
       if(dijetPt2>=ptCut){
         hTH1F["nDiJetVsPt"]->Fill(ptCut,weight);
     
         for(int ptCut_0=ptCut; ptCut_0<256; ++ptCut_0) {
           if(dijetPt1>=ptCut_0) hTH2F["nAsymDiJetVsPt"]->Fill(ptCut_0,ptCut,weight);
         }
       }
     
       if(diCenjetPt2>=ptCut){
         hTH1F["nDiCenJetVsPt"]->Fill(ptCut,weight);
     
         for(int ptCut_0=ptCut; ptCut_0<256; ++ptCut_0) {
           if(diCenjetPt1>=ptCut_0) hTH2F["nAsymDiCenJetVsPt"]->Fill(ptCut_0,ptCut,weight);
         }
       }
     
       if(diTauErPt>=ptCut)    hTH1F["nDiTauVsPt"]->Fill(ptCut,weight);
       if(diIsoTauEr>=ptCut)   hTH1F["nDiIsoTauErVsPt"]->Fill(ptCut,weight);
       if(quadjetPt>=ptCut)  hTH1F["nQuadJetVsPt"]->Fill(ptCut,weight);
       if(quadjetCPt>=ptCut) hTH1F["nQuadCenJetVsPt"]->Fill(ptCut,weight);
     
    }//loop on 256
    //
    for(int ptCut=0; ptCut<65; ++ptCut) {
      if(egPt>=ptCut)    hTH1F["nEGVsPt"]->Fill(ptCut,weight);
      if(egErPt>=ptCut)  hTH1F["nEGErVsPt"]->Fill(ptCut,weight);
      if(isoEgPt>=ptCut) hTH1F["nIsoEGVsPt"]->Fill(ptCut,weight);
      
      for(int i=0; i<5; ++i){	
	if(egErVPt[i]>=ptCut)  hTH1F["nEGEr"+std::to_string(ER[i])+"VsPt"]->Fill(ptCut,weight);
	if(isoEgErVPt[i]>=ptCut)  hTH1F["nIsoEGEr"+std::to_string(ER[i])+"VsPt"]->Fill(ptCut,weight);
      }

      if(diEG2>=ptCut)     hTH1F["nDiEGVsPt"]->Fill(ptCut,weight);
      if(diIsolEG2>=ptCut) hTH1F["nDiIsoEGVsPt"]->Fill(ptCut,weight);           
      
      for(int ptCut2=0; ptCut2<=65; ++ptCut2) {
	if(diEG1>=ptCut && diEG2>=ptCut2 && ptCut2 <= ptCut) hTH2F["nEGPtVsPt"]->Fill(ptCut,ptCut2,weight);
	if(diIsolEG1>=ptCut && diIsolEG2>=ptCut2 && ptCut2<= ptCut) hTH2F["nIsoEGPtVsPt"]->Fill(ptCut,ptCut2,weight);
      }
    }//loop on 65
    
    //---fill Eta Rate Histogram---//
    //for(int etaCut=0; etaCut<31; ++etaCut) {
    for(unsigned j=1; j<29; j++) {
      for(unsigned i=0; i<6; ++i) {
	//for(unsigned j=1; j<29; j++){
	if (egErEta[i]<=absTriggerTowerBins[j]) hTH1F["nEGEr"+std::to_string(PT[i])+"VsEta"]->Fill(absTriggerTowerBins[j]-0.001,weight);
	if (isoEgErEta[i]<=absTriggerTowerBins[j]) hTH1F["nIsoEGEr"+std::to_string(PT[i])+"VsEta"]->Fill(absTriggerTowerBins[j]-0.001,weight);
	  //}
      }
    }
    for(int etaCut=0; etaCut<26; ++etaCut) {
      for(unsigned j=0; j<4; ++j) {
	if (muErVEta[j]<=0.1*etaCut) hTH1F["nMuEr"+std::to_string(PTMu[j])+"VsEta"]->Fill(0.1*etaCut,weight);
      }
    }
    //std::cout<<std::to_string(weight)<<std::endl;
    //---End fill Eta Rate histogram---//
    
    for(int ptCut=0; ptCut<131; ++ptCut) {
      if (muPt>=ptCut)    hTH1F["nMuVsPt"]->Fill(ptCut,weight);
      if (muErPt>=ptCut)  hTH1F["nMuErVsPt"]->Fill(ptCut,weight);
      if (diMu2>=ptCut)   hTH1F["nDiMuVsPt"]->Fill(ptCut,weight);

      for(int j=0; j<3; j++) {
	if (muErVPt[j]>=ptCut)   hTH1F["nMuEr"+std::to_string(ERMu[j])+"VsPt"]->Fill(ptCut,weight);
      }
    }
    for(int iCut=0; iCut<41; ++iCut) {
      for(int iCut2=0; iCut2<=iCut; ++iCut2) {
	float ptCut = iCut*0.5;
	float ptCut2 = iCut2*0.5;
	if (doubleMuPt1>=ptCut && doubleMuPt2>=ptCut2) hTH2F["nMuPtVsPt"]->Fill(ptCut,ptCut2,weight);
	if (oniaMuPt1>=ptCut && oniaMuPt2>=ptCut2)     hTH2F["nOniaMuPtVsPt"]->Fill(ptCut,ptCut2,weight);
      }
    }
    for(int httCut=0; httCut<512; ++httCut) {
      if(htt>httCut) hTH1F["nHTTVsHTT"]->Fill(httCut,weight);
      if(ett>httCut) hTH1F["nETTVsETT"]->Fill(httCut,weight);
      if(etm>httCut) hTH1F["nETMVsETM"]->Fill(httCut,weight);
    }
  } // end event loop
  
  float scaleFactor(1.);
  if (runOnData){
    std::cout << "# of lumis sections used for rate computation : " << nLumi << std::endl;
    //scaleFactor = (80.*631.)/(1326*23.3);      
    scaleFactor = (80.*631.)/(nLumi*23.3);      
  }else{
    std::cout << "# of zero bias events (weighted) used for rate computation : " << nZeroBias << std::endl;
    scaleFactor = ScaleFactor(nZeroBias,nBunches);    
  }
  float scaleFactorTest(1.);
  scaleFactorTest = ScaleFactor(nZeroBias,nBunches);
  std::cout << "Scale factor applied to histograms = " << scaleFactor << std::endl;
  std::cout << "Scale factor derived from nBunches applied to histograms = " << scaleFactorTest << std::endl;

  std::map<std::string,TH1F*>::iterator hTH1FIt  = hTH1F.begin();
  std::map<std::string,TH1F*>::iterator hTH1FEnd = hTH1F.end();
  
  for(; hTH1FIt!=hTH1FEnd; ++hTH1FIt) {
    TH1F* histo = hTH1FIt->second;
    //setRateError(histo);
    histo->Sumw2();
    histo->Scale(scaleFactor);
  }

  std::map<std::string,TH2F*>::iterator hTH2FIt  = hTH2F.begin();
  std::map<std::string,TH2F*>::iterator hTH2FEnd = hTH2F.end();

  for(; hTH2FIt!=hTH2FEnd; ++hTH2FIt) {
    TH2F* histo = hTH2FIt->second;
    histo->Sumw2();
    histo->Scale(scaleFactor);
  }
  
  std::map<std::string,TH1F*>::iterator hCountTH1FIt  = hCountTH1F.begin();
  std::map<std::string,TH1F*>::iterator hCountTH1FEnd = hCountTH1F.end();

  for(; hCountTH1FIt!=hCountTH1FEnd; ++hCountTH1FIt) {
    TH1F* histo = hCountTH1FIt->second;
    //setRateError(histo);
    histo->Sumw2();
    //histo->Scale(scaleFactor);
  }
  
  for (int i=0; i<57; i++){
    float eg =  hTH1F["nEGVsEta"]->GetBinContent(i)/hTH1F["nEGVsEta"]->GetXaxis()->GetBinWidth(i);
    float isoeg = hTH1F["nIsoEGVsEta"]->GetBinContent(i)/hTH1F["nIsoEGVsEta"]->GetXaxis()->GetBinWidth(i);
    hTH1F["nEGVsEta"]->SetBinContent(i,eg);
    hTH1F["nIsoEGVsEta"]->SetBinContent(i,isoeg);
  }
  for (int i=0; i<65; i++){
    float jet = hTH1F["nJetVsEta"]->GetBinContent(i)/hTH1F["nJetVsEta"]->GetXaxis()->GetBinWidth(i);
    hTH1F["nJetVsEta"]->SetBinContent(i,jet);
  } 
  
  outFile->Write();
  outFile->Close();
  delete outFile;
}

// --------------------------------------------------------------------

void goRatePlots(std::string fileType, int isCrossSec = false, int nEvents = 0) 
{
  //int nBunches50ns = 1368;
  int nBunches25ns = 2508; //2508 is what agreed with TSG for # bunches
  //int nBunches25ns_run256843 = 1021;
  int nBunches = -1;

  float xSec13TeV = isCrossSec ? 78.26 : 80.; // Using McM for cross section comparison and 80 (agreed with TSG) for rates
  //float xSec8TeV  = 72.7; 

  std::string filename;
  bool isData(true);
  bool isReweight(false);
  
  if (fileType == "13TEV_40PU_2016_RE-EMUL")
    {
      isData = false;
      nBunches = nBunches25ns;
      filename = "/afs/cern.ch/user/p/pellicci/data2/L1DPG/root/2016/v2/40PU_25ns_Stage2/40PU_25ns_Stage2_1.root";
    }
  //*************************************************************************//
  else if (fileType == "run80xMC_zerobias_l1tint_v63p1_muonFix_PU28")
    {
      isData = false;
      nBunches = 2736;
      isReweight = false;
      filename = "/afs/cern.ch/user/z/zhangj/test/L1Menu/80x_stage2/v62p3_MC2016/CMSSW_8_0_9/src/L1TriggerDPG/L1Menu/macros/r80xMC_809_zerobias_l1tint_v62p3_muonFix.txt";
    }
  else if (fileType == "run274241_zerobias_l1tint_v61p1_test_v2")
    {
      isData = false;      
      filename = "/afs/cern.ch/user/z/zhangj/test/L1Menu/80x_stage2/v62p3_MC2016/CMSSW_8_0_9/src/L1TriggerDPG/L1Menu/macros/r274241_809_zerobias_l1tint_v61p1.txt";
      nBunches = 1752;
      isReweight = false;
    }
  else if (fileType == "run273725_stage2_l1tint_v58p1")
    {
      isData = false;
      filename = "/afs/cern.ch/user/z/zhangj/test/L1Menu/80x_stage2/v58p1_MC2015/CMSSW_8_0_8/src/L1TriggerDPG/L1Menu/macros/r273725_808_zerobias_l1t-int-v58p1.txt";
      nBunches = 2736;
      isReweight = false;
    }  
  //**************************************************************************//
  else if (fileType == "RUN256843_Stage2_8X")
    {
      isData = true;      
      //isData = false;      
      //nBunches = 1021;
      // filename = "/data/user/gennai/L1Ntuple/l1t_debug-stage-2_256843.root";
      // filename = "root://cmseos.fnal.gov//store/user/lpctrig/apana/Stage2/ZeroBias1/crab_ZeroBias1_Run2015D-v1/151230_012024/0000/l1t_stage2_2.root";
      filename = "ntuple/Run256843_stage2_8X.list";
      //filename = "ntuples_256843_stage2_full.list";
    }
  else if (fileType == "RUN256843_Stage2_76")
    {
      isData = true;      
      // filename = "/data/user/gennai/L1Ntuple/l1t_debug-stage-2_256843.root";
      // filename = "root://cmseos.fnal.gov//store/user/lpctrig/apana/Stage2/ZeroBias1/crab_ZeroBias1_Run2015D-v1/151230_012024/0000/l1t_stage2_2.root";
      filename = "ntuple/Run256843_stage2_Len.list";
      //filename = "ntuples_256843_stage2_full.list";
    }
  else if (fileType == "Stage2_Simone")
    {
      isData = false;      
      nBunches = 1021;
      // filename = "/data/user/gennai/L1Ntuple/l1t_debug-stage-2_256843.root";
      filename = "ntuple/Run256843_stage2_Simone.list";
      //filename = "ntuples_256843_stage2_Simone.list";
    }
  else if (fileType == "RUN259721_Stage2")
  {
    isData = false;      
    nBunches = 517; 
    // filename = "/data/user/gennai/L1Ntuple/l1t_debug-stage-2_256843.root";
    filename = "ntuple/Run259721_stage2_Simone.list";
    //filename = "ntuples_256843_stage2_Simone.list";
    }
  else if (fileType == "RUN260627_Aaron")
    {
      isData = false;      
      nBunches = 1021;
      // filename = "/data/user/gennai/L1Ntuple/l1t_debug-stage-2_256843.root";
      filename = "ntuples_260627_Aaron.list";
    }
  else 
    {
      std::cout << "Config param " << fileType << " invalid! \n"
		<< "Valid fileType values are : DATA, 8TEV_TF_DATA, 8TEV_TF_2012_RE-EMUL, "
		<< "8TEV_25PU_ORIG_RE-EMUL, 8TEV_25PU_2012_RE-EMUL, 8TEV_25PU_2012GCT10GEV_RE-EMUL, 8TEV_25PU_2015_RE-EMUL, "
		<< "13TEV_25PU_ORIG_RE-EMUL, 13TEV_25PU_2012_RE-EMUL, 13TEV_25PU_2012GCT10GEV_RE-EMUL, 13TEV_25PU_2015_RE-EMUL\n";
    }

  // TTree *tree;
  // TFile f(filename.c_str());
  // TDirectory * dir = (TDirectory*)f.Get("l1UpgradeTree");
  // dir->GetObject("L1UpgradeTree",tree);

  BasicRatePlots basicRatePlots(filename); 
  basicRatePlots.run(isData,fileType,xSec13TeV,nBunches,isCrossSec,nEvents,isReweight);    
}

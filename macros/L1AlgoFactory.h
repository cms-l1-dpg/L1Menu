#include "L1Ntuple.h"
#include<iostream>

class L1AlgoFactory: public L1Ntuple{
 public:
  L1AlgoFactory(TTree *tree);

  void SingleMuPt(Float_t& ptcut, Bool_t isER, Int_t qualmin=4);
  void DoubleMuPt(Float_t& mu1pt, Float_t& mu2pt, Bool_t isHighQual = false, Bool_t isER = false);
  void OniaPt(Float_t& ptcut1, Float_t& ptcut2, Int_t delta);
  void Onia2015Pt(Float_t& ptcut1, Float_t& ptcut2, Bool_t isER, Bool_t isOS, Int_t delta);
  void TripleMuPt(Float_t& mu1pt, Float_t& mu2pt, Float_t& mu3pt, Int_t qualmin = 4);
  void QuadMuPt(Float_t& mu1pt, Float_t& mu2pt, Float_t& mu3pt, Float_t& mu4pt, Int_t qualmin = 4);

  void SingleEGPt(Float_t& ptcut, Bool_t isIsolated, Bool_t isER);
  void DoubleEGPt(Float_t& ele1pt, Float_t& ele2pt, Bool_t isIsolated = false, Bool_t isER = false);
  void TripleEGPt(Float_t& ele1pt, Float_t& ele2pt, Float_t& ele3pt);

  void SingleJetPt(Float_t& ptcut, Bool_t isCentral = false);
  void DoubleJetPt(Float_t& cut1, Float_t& cut2, Bool_t isCentral = false);
  void DoubleJet_Eta1p7_deltaEta4Pt(Float_t& cut1, Float_t& cut2 );
  void DoubleTauJetEta2p17Pt(Float_t& cut1, Float_t& cut2, Bool_t isIsolated = false);
  void TripleJetPt(Float_t& cut1, Float_t& cut2, Float_t& cut3, Bool_t isCentral = false);
  Bool_t TripleJet_VBF(Float_t jet1, Float_t jet2, Float_t jet3, Int_t jetclass = 0);
  void QuadJetPt(Float_t& cut1, Float_t& cut2, Float_t& cut3, Float_t& cut4, Bool_t isCentral = false);

  void Mu_EGPt(Float_t& mucut, Float_t& EGcut, Bool_t isIsolated = false, Int_t qualmin=4);
  void DoubleMu_EGPt(Float_t& mucut, Float_t& EGcut, Bool_t isMuHighQual = false );
  void Mu_DoubleEGPt(Float_t& mucut, Float_t& EGcut );

  void Muer_JetCentralPt(Float_t& mucut, Float_t& jetcut);
  void Mu_JetCentral_deltaPt(Float_t& mucut, Float_t& jetcut);
  void Mu_DoubleJetCentralPt(Float_t& mucut, Float_t& jetcut);
  void Muer_TauJetEta2p17Pt(Float_t& mucut, Float_t& taucut, Bool_t isIsolated = false);

  void EG_FwdJetPt(Float_t& EGcut, Float_t& FWcut);
  void EG_DoubleJetCentralPt(Float_t& EGcut, Float_t& jetcut);
  void EGer_TripleJetCentralPt(Float_t& EGcut, Float_t& jetcut);
  void IsoEGer_TauJetEta2p17Pt(Float_t& egcut, Float_t& taucut);

  void QuadJetCentral_TauJetPt(Float_t& jetcut, Float_t& taucut);

  void ETMVal(Float_t& ETMcut);
  void HTTVal(Float_t& HTTcut);
  void HTMVal(Float_t& HTMcut);
  void ETTVal(Float_t& ETTcut);


  Bool_t SingleMu(Float_t ptcut, Bool_t isER, Int_t qualmin=4);
  Bool_t DoubleMu(Float_t mu1pt, Float_t mu2pt, Bool_t isHighQual = false, Bool_t isER = false);
  Bool_t Onia(Float_t mu1pt, Float_t mu2pt, Int_t delta);
  Bool_t Onia2015(Float_t mu1pt, Float_t mu2pt, Bool_t isER, Bool_t isOS, Int_t delta);
  Bool_t TripleMu(Float_t mu1pt, Float_t mu2pt, Float_t mu3pt, Int_t qualmin);
  Bool_t QuadMu(Float_t mu1pt, Float_t mu2pt, Float_t mu3pt, Float_t mu4pt, Int_t qualmin);

  Bool_t SingleEG(Float_t ptcut, Bool_t isIsolated, Bool_t isER);
  Bool_t DoubleEG(Float_t ptcut1, Float_t ptcut2, Bool_t isIsolated = false);
  Bool_t TripleEG(Float_t ptcut1, Float_t ptcut2, Float_t ptcut3);

  Bool_t SingleJet(Float_t ptcut, Bool_t isCentral = false);
  Bool_t DoubleJet(Float_t cut1, Float_t cut2, Bool_t isCentral = false);
  Bool_t DoubleJet_Eta1p7_deltaEta4(Float_t cut1, Float_t cut2 );
  Bool_t DoubleTauJetEta2p17(Float_t cut1, Float_t cut2, Bool_t isIsolated = false);
  Bool_t TripleJet(Float_t cut1, Float_t cut2, Float_t cut3, Bool_t isCentral = false);
  Bool_t QuadJet(Float_t cut1, Float_t cut2, Float_t cut3, Float_t cut4, Bool_t isCentral);

  Bool_t Mu_EG(Float_t mucut, Float_t EGcut, Bool_t isIsolated = false, Int_t qualmin=4);
  Bool_t DoubleMu_EG(Float_t mucut, Float_t EGcut, Bool_t isMuHighQual = false);
  Bool_t Mu_DoubleEG(Float_t mucut, Float_t EGcut);

  Bool_t Muer_JetCentral(Float_t mucut, Float_t jetcut);
  Bool_t Mu_JetCentral_delta(Float_t mucut, Float_t jetcut);
  Bool_t Mu_DoubleJetCentral(Float_t mucut, Float_t jetcut);
  Bool_t Muer_TauJetEta2p17(Float_t mucut, Float_t taucut, Bool_t isIsolated = false);

  Bool_t EG_FwdJet(Float_t EGcut, Float_t FWcut);
  Bool_t EG_DoubleJetCentral(Float_t EGcut, Float_t jetcut);
  Bool_t EGer_TripleJetCentral(Float_t EGcut, Float_t jetcut);
  Bool_t IsoEGer_TauJetEta2p17(Float_t egcut, Float_t taucut);

  Bool_t QuadJetCentral_TauJet(Float_t jetcut, Float_t taucut);


  Bool_t ETM(Float_t ETMcut);
  Bool_t HTT(Float_t HTTcut);
  Bool_t HTM(Float_t HTMcut);
  Bool_t ETT(Float_t ETTcut);

  inline Bool_t correlateInPhi(Int_t jetphi, Int_t muphi, Int_t delta=1);
  inline Bool_t correlateInEta(Int_t mueta, Int_t jeteta, Int_t delta=1);
  Int_t etaMuIdx(Double_t eta);
  Int_t phiINjetCoord(Double_t phi);
  Int_t etaINjetCoord(Double_t eta);
  inline Double_t degree(Double_t radian);

 private:
};

const size_t PHIBINS = 18;
const Double_t PHIBIN[] = {10,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350};
const size_t ETABINS = 23;
const Double_t ETABIN[] = {-5.,-4.5,-4.,-3.5,-3.,-2.172,-1.74,-1.392,-1.044,-0.696,-0.348,0.,0.348,0.696,1.044,1.392,1.74,2.172,3.,3.5,4.,4.5,5.};

const size_t ETAMUBINS = 65;
const Double_t ETAMU[] = { -2.45,-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45 };

L1AlgoFactory::L1AlgoFactory(TTree *tree) : L1Ntuple(tree) {}

Bool_t L1AlgoFactory::SingleMu(Float_t ptcut, Bool_t isER, Int_t qualmin) {
  Float_t tmp_cut = -10.;
  SingleMuPt(tmp_cut,isER, qualmin);
  if(tmp_cut >= ptcut) return true;
  return false;
}

Bool_t L1AlgoFactory::DoubleMu(Float_t mu1pt, Float_t mu2pt, Bool_t isHighQual, Bool_t isER) {
  Float_t tmp_cut1 = -10.;
  Float_t tmp_cut2 = -10.;
  DoubleMuPt(tmp_cut1,tmp_cut2,isHighQual,isER);
  if(tmp_cut1 >= mu1pt && tmp_cut2 >= mu2pt) return true;
  return false;
}

Bool_t L1AlgoFactory::Onia(Float_t mu1pt, Float_t mu2pt, Int_t delta) {
  Float_t tmp_cut1 = -10.;
  Float_t tmp_cut2 = -10.;
  OniaPt(tmp_cut1,tmp_cut2,delta);
  if(tmp_cut1 >= mu1pt && tmp_cut2 >= mu2pt) return true;
  return false;
}

Bool_t L1AlgoFactory::Onia2015(Float_t mu1pt, Float_t mu2pt, Bool_t isER, Bool_t isOS, Int_t delta) {
  Float_t tmp_cut1 = -10.;
  Float_t tmp_cut2 = -10.;
  Onia2015Pt(tmp_cut1,tmp_cut2,isER,isOS,delta);
  if(tmp_cut1 >= mu1pt && tmp_cut2 >= mu2pt) return true;
  return false;
}

Bool_t L1AlgoFactory::TripleMu(Float_t mu1pt, Float_t mu2pt, Float_t mu3pt, Int_t qualmin) {
  Float_t tmp_cut1 = -10.;
  Float_t tmp_cut2 = -10.;
  Float_t tmp_cut3 = -10.;
  TripleMuPt(tmp_cut1,tmp_cut2,tmp_cut3,qualmin);
  if(tmp_cut1 >= mu1pt && tmp_cut2 >= mu2pt &&  tmp_cut3 >= mu3pt) return true;
  return false;
}

Bool_t L1AlgoFactory::QuadMu(Float_t mu1pt, Float_t mu2pt, Float_t mu3pt, Float_t mu4pt, Int_t qualmin) {
  Float_t tmp_cut1 = -10.;
  Float_t tmp_cut2 = -10.;
  Float_t tmp_cut3 = -10.;
  Float_t tmp_cut4 = -10.;
  QuadMuPt(tmp_cut1,tmp_cut2,tmp_cut3,tmp_cut4,qualmin);
  if(tmp_cut1 >= mu1pt && tmp_cut2 >= mu2pt &&  tmp_cut3 >= mu3pt && tmp_cut4 >= mu4pt) return true;
  return false;
}

Bool_t L1AlgoFactory::SingleEG(Float_t ptcut, Bool_t isIsolated, Bool_t isER) {
  Float_t tmp_cut1 = -10.;
  SingleEGPt(tmp_cut1,isIsolated,isER);
  if(tmp_cut1 >= ptcut) return true;
  return false;
}

Bool_t L1AlgoFactory::DoubleEG(Float_t ptcut1, Float_t ptcut2, Bool_t isIsolated) {
  Float_t tmp_cut1 = -10.;
  Float_t tmp_cut2 = -10.;
  DoubleEGPt(tmp_cut1,tmp_cut2,isIsolated);
  if(tmp_cut1 >= ptcut1 && tmp_cut2 >= ptcut2) return true;
  return false;
}

Bool_t L1AlgoFactory::TripleEG(Float_t ptcut1, Float_t ptcut2, Float_t ptcut3) {
  Float_t tmp_cut1 = -10.;
  Float_t tmp_cut2 = -10.;
  Float_t tmp_cut3 = -10.;
  TripleEGPt(tmp_cut1,tmp_cut2,tmp_cut3);
  if(tmp_cut1 >= ptcut1 && tmp_cut2 >= ptcut2 && tmp_cut3 >= ptcut3) return true;
  return false;
}

Bool_t L1AlgoFactory::SingleJet(Float_t ptcut, Bool_t isCentral) {
  Float_t tmp_cut1 = -10.;
  SingleJetPt(tmp_cut1,isCentral);
  if(tmp_cut1 >= ptcut) return true;
  return false;
}

Bool_t L1AlgoFactory::DoubleJet(Float_t ptcut1, Float_t ptcut2, Bool_t isCentral) {
  Float_t tmp_cut1 = -10.;
  Float_t tmp_cut2 = -10.;
  DoubleJetPt(tmp_cut1,tmp_cut2,isCentral);
  if(tmp_cut1 >= ptcut1 && tmp_cut2 >= ptcut2) return true;
  return false;
}

Bool_t L1AlgoFactory::DoubleJet_Eta1p7_deltaEta4(Float_t ptcut1, Float_t ptcut2) {
  Float_t tmp_cut1 = -10.;
  Float_t tmp_cut2 = -10.;
  DoubleJet_Eta1p7_deltaEta4Pt(tmp_cut1,tmp_cut2);
  if(tmp_cut1 >= ptcut1 && tmp_cut2 >= ptcut2) return true;
  return false;
}

Bool_t L1AlgoFactory::DoubleTauJetEta2p17(Float_t ptcut1, Float_t ptcut2, Bool_t isIsolated) {
  Float_t tmp_cut1 = -10.;
  Float_t tmp_cut2 = -10.;
  DoubleTauJetEta2p17Pt(tmp_cut1,tmp_cut2,isIsolated);
  if(tmp_cut1 >= ptcut1 && tmp_cut2 >= ptcut2) return true;
  return false;
}

Bool_t L1AlgoFactory::TripleJet(Float_t ptcut1, Float_t ptcut2, Float_t ptcut3, Bool_t isCentral) {
  Float_t tmp_cut1 = -10.;
  Float_t tmp_cut2 = -10.;
  Float_t tmp_cut3 = -10.;
  TripleJetPt(tmp_cut1,tmp_cut2,tmp_cut3,isCentral);
  if(tmp_cut1 >= ptcut1 && tmp_cut2 >= ptcut2 && tmp_cut3 >= ptcut3) return true;
  return false;
}

Bool_t L1AlgoFactory::QuadJet(Float_t ptcut1, Float_t ptcut2, Float_t ptcut3, Float_t ptcut4, Bool_t isCentral) {
  Float_t tmp_cut1 = -10.;
  Float_t tmp_cut2 = -10.;
  Float_t tmp_cut3 = -10.;
  Float_t tmp_cut4 = -10.;
  QuadJetPt(tmp_cut1,tmp_cut2,tmp_cut3,tmp_cut4,isCentral);
  if(tmp_cut1 >= ptcut1 && tmp_cut2 >= ptcut2 && tmp_cut3 >= ptcut3 && tmp_cut4 >= ptcut4) return true;
  return false;
}

Bool_t L1AlgoFactory::Mu_EG(Float_t mucut, Float_t EGcut, Bool_t isIsolated, Int_t qualmin) {
  Float_t tmp_mucut = -10.;
  Float_t tmp_elecut = -10.;
  Mu_EGPt(tmp_mucut,tmp_elecut,isIsolated,qualmin);
  if(tmp_mucut >= mucut && tmp_elecut >= EGcut) return true;
  return false;
}

Bool_t L1AlgoFactory::DoubleMu_EG(Float_t mucut, Float_t EGcut, Bool_t isMuHighQual) {
  Float_t tmp_mucut = -10.;
  Float_t tmp_elecut = -10.;
  DoubleMu_EGPt(tmp_mucut,tmp_elecut,isMuHighQual);
  if(tmp_mucut >= mucut && tmp_elecut >= EGcut) return true;
  return false;
}

Bool_t L1AlgoFactory::Mu_DoubleEG(Float_t mucut, Float_t EGcut) {
  Float_t tmp_mucut = -10.;
  Float_t tmp_elecut = -10.;
  Mu_DoubleEGPt(tmp_mucut,tmp_elecut);
  if(tmp_mucut >= mucut && tmp_elecut >= EGcut) return true;
  return false;
}

Bool_t L1AlgoFactory::Muer_JetCentral(Float_t mucut, Float_t jetcut) {
  Float_t tmp_mucut = -10.;
  Float_t tmp_jetcut = -10.;
  Muer_JetCentralPt(tmp_mucut,tmp_jetcut);
  if(tmp_mucut >= mucut && tmp_jetcut >= jetcut) return true;
  return false;
}

Bool_t L1AlgoFactory::Mu_JetCentral_delta(Float_t mucut, Float_t jetcut) {
  Float_t tmp_mucut = -10.;
  Float_t tmp_jetcut = -10.;
  Mu_JetCentral_deltaPt(tmp_mucut,tmp_jetcut);
  if(tmp_mucut >= mucut && tmp_jetcut >= jetcut) return true;
  return false;
}

Bool_t L1AlgoFactory::Mu_DoubleJetCentral(Float_t mucut, Float_t jetcut) {
  Float_t tmp_mucut = -10.;
  Float_t tmp_jetcut = -10.;
  Mu_DoubleJetCentralPt(tmp_mucut,tmp_jetcut);
  if(tmp_mucut >= mucut && tmp_jetcut >= jetcut) return true;
  return false;
}

Bool_t L1AlgoFactory::EG_FwdJet(Float_t egcut, Float_t FWcut) {
  Float_t tmp_egcut = -10.;
  Float_t tmp_FWcut = -10.;
  EG_FwdJetPt(tmp_egcut,tmp_FWcut);
  if(tmp_egcut >= egcut && tmp_FWcut >= FWcut) return true;
  return false;
}

Bool_t L1AlgoFactory::EG_DoubleJetCentral(Float_t egcut, Float_t jetcut) {
  Float_t tmp_egcut = -10.;
  Float_t tmp_jetcut = -10.;
  EG_DoubleJetCentralPt(tmp_egcut,tmp_jetcut);
  if(tmp_egcut >= egcut && tmp_jetcut >= jetcut) return true;
  return false;
}

Bool_t L1AlgoFactory::EGer_TripleJetCentral(Float_t egcut, Float_t jetcut) {
  Float_t tmp_egcut = -10.;
  Float_t tmp_jetcut = -10.;
  EGer_TripleJetCentralPt(tmp_egcut,tmp_jetcut);
  if(tmp_egcut >= egcut && tmp_jetcut >= jetcut) return true;
  return false;
}

Bool_t L1AlgoFactory::Muer_TauJetEta2p17(Float_t mucut, Float_t taucut, Bool_t isIsolated){
  Float_t tmp_mucut  = -10.;
  Float_t tmp_taucut = -10.;
  Muer_TauJetEta2p17Pt(tmp_mucut, tmp_taucut,isIsolated);
  if(tmp_mucut >= mucut && tmp_taucut >= taucut) return true;
  return false;
}

Bool_t L1AlgoFactory::IsoEGer_TauJetEta2p17(Float_t egcut, Float_t taucut){
  Float_t tmp_egcut  = -10.;
  Float_t tmp_taucut = -10.;
  IsoEGer_TauJetEta2p17Pt(tmp_egcut, tmp_taucut);
  if(tmp_egcut >= egcut && tmp_taucut >= taucut) return true;
  return false;
}

Bool_t L1AlgoFactory::QuadJetCentral_TauJet(Float_t jetcut, Float_t taucut){
  Float_t tmp_jetcut = -10.;
  Float_t tmp_taucut = -10.;
  QuadJetCentral_TauJetPt(tmp_jetcut,tmp_taucut);
  if(tmp_jetcut >= jetcut && tmp_taucut >= taucut) return true;
  return false;
}

inline Bool_t L1AlgoFactory::correlateInPhi(Int_t jetphi, Int_t muphi, Int_t delta){
  return fabs(muphi-jetphi) < fabs(1 + delta) || fabs(muphi-jetphi) > fabs(PHIBINS - 1 - delta) ;
}

inline Bool_t L1AlgoFactory::correlateInEta(Int_t mueta, Int_t jeteta, Int_t delta) {
  return fabs(mueta-jeteta) < 1 + delta;
}

Int_t L1AlgoFactory::etaMuIdx(Double_t eta) {
  size_t etaIdx = 0.;
  for (size_t idx=0; idx<ETAMUBINS; idx++) {
    if (eta>=ETAMU[idx] and eta<ETAMU[idx+1])
      etaIdx = idx;
  }

  return int(etaIdx);
}

Int_t L1AlgoFactory::phiINjetCoord(Double_t phi) {
  size_t phiIdx = 0;
  Double_t phidegree = degree(phi);
  for (size_t idx=0; idx<PHIBINS; idx++) {
    if (phidegree>=PHIBIN[idx] and phidegree<PHIBIN[idx+1])
      phiIdx = idx;
    else if (phidegree>=PHIBIN[PHIBINS-1] || phidegree<=PHIBIN[0])
      phiIdx = idx;
  }
  phiIdx = phiIdx + 1;
  if (phiIdx == 18)  phiIdx = 0;

  return int(phiIdx);
}

Int_t L1AlgoFactory::etaINjetCoord(Double_t eta) {
  size_t etaIdx = 0.;
  for (size_t idx=0; idx<ETABINS; idx++) {
    if (eta>=ETABIN[idx] and eta<ETABIN[idx+1])
      etaIdx = idx;
  }

  return int(etaIdx);
}

inline Double_t L1AlgoFactory::degree(Double_t radian) {
  if (radian<0.)
    return 360.+(radian/acos(-1.)*180.);
  else
    return radian/acos(-1.)*180.;
}


Bool_t L1AlgoFactory::ETM(Float_t ETMcut) {
  Float_t tmp_cut = -10.;
  ETMVal(tmp_cut);
  if(tmp_cut >= ETMcut) return true;
  return false;
}

Bool_t L1AlgoFactory::HTT(Float_t HTTcut) {
  Float_t tmp_cut = -10.;
  HTTVal(tmp_cut);
  if(tmp_cut >= HTTcut) return true;
  return false;
}

Bool_t L1AlgoFactory::HTM(Float_t HTMcut) {
  Float_t tmp_cut = -10.;
  HTMVal(tmp_cut);
  if(tmp_cut >= HTMcut) return true;
  return false;
}

Bool_t L1AlgoFactory::ETT(Float_t ETTcut) {
  Float_t tmp_cut = -10.;
  ETTVal(tmp_cut);
  if(tmp_cut >= ETTcut) return true;
  return false;
}





void L1AlgoFactory::SingleMuPt(Float_t& ptcut, Bool_t isER, Int_t qualmin) {

  if(nMuons < 1) return;

  Float_t ptmax = -10.;

  for(UInt_t imu=0; imu < nMuons; imu++) {
    Int_t bx = muonBx.at(imu);
    if(bx != 0) continue;
    Float_t eta = muonEta.at(imu);        
    if(fabs(eta) > 2.1 && isER) continue;
    Float_t pt = muonEt.at(imu);                       

    //if(isER == false) std::cout << "pt = " << pt << std::endl;
    if(pt >= ptmax) ptmax = pt;
 }

  ptcut = ptmax;
  return;
}

void L1AlgoFactory::DoubleMuPt(Float_t& cut1, Float_t& cut2, Bool_t isHighQual, Bool_t isER) {

  Float_t mu1ptmax = -10.;
  Float_t mu2ptmax = -10.;

  if(nMuons < 2) return;

  for (UInt_t imu=0; imu < nMuons; imu++) {
    Int_t bx = muonBx.at(imu);
    if(bx != 0) continue;
    Float_t pt = muonEt.at(imu);			
    Float_t eta = muonEta.at(imu);
    if(isER && fabs(eta) > 2.1) continue;

    if(pt >= mu1ptmax)
      {
	mu2ptmax = mu1ptmax;
	mu1ptmax = pt;
      }
    else if(pt >= mu2ptmax) mu2ptmax = pt;
  }

  if(mu2ptmax >= 0.){
    cut1 = mu1ptmax;
    cut2 = mu2ptmax;
  }

  return;
}

void L1AlgoFactory::OniaPt(Float_t& ptcut1, Float_t& ptcut2, Int_t delta) {

  if(nMuons < 2) return;

  Float_t maxpt1 = -10.;
  Float_t maxpt2 = -10.;
  Float_t corr = false;

  std::vector<std::pair<Float_t,Float_t> > muonPairs;

  for (UInt_t imu=0; imu < nMuons; imu++) {
    Int_t bx = muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = muonEt.at(imu);			
    Float_t eta = muonEta.at(imu);
    if (fabs(eta) > 2.1) continue;
    Int_t ieta1 = etaMuIdx(eta);

    for (UInt_t imu2=0; imu2 < nMuons; imu2++) {
      if (imu2 == imu) continue;
      Int_t bx2 = muonBx.at(imu2);		
      if(bx2 != 0) continue;
      Float_t pt2 = muonEt.at(imu2);			
      Float_t eta2 = muonEta.at(imu2);
      if (fabs(eta2) > 2.1) continue;
      Int_t ieta2 = etaMuIdx(eta2);

      Float_t deta = ieta1 - ieta2; 
      if ( fabs(deta) <= delta){
	corr = true;
	muonPairs.push_back(std::pair<Float_t,Float_t>(pt,pt2));
      }

    }
  }

  if(corr){
    std::vector<std::pair<Float_t,Float_t> >::const_iterator muonPairIt  = muonPairs.begin();
    std::vector<std::pair<Float_t,Float_t> >::const_iterator muonPairEnd = muonPairs.end();
    for (; muonPairIt != muonPairEnd; ++muonPairIt) {
      Float_t pt1 = muonPairIt->first;
      Float_t pt2 = muonPairIt->second;
      
      if ( pt1 > maxpt1 || (fabs(maxpt1-pt1)<10E-2 && pt2>maxpt2) ) 
	{
	  maxpt1 = pt1;
	  maxpt2 = pt2;
	}
    }

  }

  if(corr && maxpt2 >= 0.){
    ptcut1 = maxpt1;
    ptcut2 = maxpt2;
  }

  return;
}

void L1AlgoFactory::Onia2015Pt(Float_t& ptcut1, Float_t& ptcut2, Bool_t isER, Bool_t isOS, Int_t delta) {

  if(nMuons < 2) return;

  Float_t maxpt1 = -10.;
  Float_t maxpt2 = -10.;
  Float_t corr = false;

  std::vector<std::pair<Float_t,Float_t> > muonPairs;

  for (UInt_t imu=0; imu < nMuons; imu++) {
    Int_t bx = muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = muonEt.at(imu);			
    Float_t eta = muonEta.at(imu);
    if(isER && fabs(eta) > 1.6) continue;
    Int_t ieta1 = etaMuIdx(eta);
    Int_t charge1 = muonChg.at(imu);

    for (UInt_t imu2=0; imu2 < nMuons; imu2++) {
      if (imu2 == imu) continue;
      Int_t bx2 = muonBx.at(imu2);		
      if(bx2 != 0) continue;
      Float_t pt2 = muonEt.at(imu2);			
      Float_t eta2 = muonEta.at(imu2);
      if(isER && fabs(eta2) > 1.6) continue;
      Int_t ieta2 = etaMuIdx(eta2);
      Int_t charge2 = muonChg.at(imu2);

      if(isOS && charge1*charge2 > 0) continue;

      Float_t deta = ieta1 - ieta2; 
      if(fabs(deta) <= delta){
	corr = true;
	muonPairs.push_back(std::pair<Float_t,Float_t>(pt,pt2));
      }

    }
  }

  if(corr){
    std::vector<std::pair<Float_t,Float_t> >::const_iterator muonPairIt  = muonPairs.begin();
    std::vector<std::pair<Float_t,Float_t> >::const_iterator muonPairEnd = muonPairs.end();
    for(; muonPairIt != muonPairEnd; ++muonPairIt){
      Float_t pt1 = muonPairIt->first;
      Float_t pt2 = muonPairIt->second;
      
      if(pt1 > maxpt1 || (fabs(maxpt1-pt1)<10E-2 && pt2>maxpt2) ) 
	{
	  maxpt1 = pt1;
	  maxpt2 = pt2;
	}
    }

  }

  if(corr && maxpt2 >= 0.){
    ptcut1 = maxpt1;
    ptcut2 = maxpt2;
  }

  return;
}

void L1AlgoFactory::TripleMuPt(Float_t& cut1, Float_t& cut2, Float_t& cut3, Int_t qualmin) {

  Float_t mu1ptmax = -10.;
  Float_t mu2ptmax = -10.;
  Float_t mu3ptmax = -10.;

  if(nMuons < 3) return;

  for (UInt_t imu=0; imu < nMuons; imu++) {
    Int_t bx = muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = muonEt.at(imu);			

    if(pt >= mu1ptmax)
      {
	mu3ptmax = mu2ptmax;
	mu2ptmax = mu1ptmax;
	mu1ptmax = pt;
      }
    else if(pt >= mu2ptmax){
      mu3ptmax = mu2ptmax;
      mu2ptmax = pt;
    }
    else if(pt >= mu3ptmax) mu3ptmax = pt;
  }

  if(mu3ptmax >= 0.){
    cut1 = mu1ptmax;
    cut2 = mu2ptmax;
    cut3 = mu3ptmax;
  }

  return;
}

void L1AlgoFactory::QuadMuPt(Float_t& cut1, Float_t& cut2, Float_t& cut3, Float_t& cut4, Int_t qualmin) {

  Float_t mu1ptmax = -10.;
  Float_t mu2ptmax = -10.;
  Float_t mu3ptmax = -10.;
  Float_t mu4ptmax = -10.;

  if(nMuons < 4) return;

  for (UInt_t imu=0; imu < nMuons; imu++) {
    Int_t bx = muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = muonEt.at(imu);			

    if(pt >= mu1ptmax)
      {
	mu4ptmax = mu3ptmax;
	mu3ptmax = mu2ptmax;
	mu2ptmax = mu1ptmax;
	mu1ptmax = pt;
      }
    else if(pt >= mu2ptmax){
      mu4ptmax = mu3ptmax;
      mu3ptmax = mu2ptmax;
      mu2ptmax = pt;
    }
    else if(pt >= mu3ptmax){
      mu4ptmax = mu3ptmax;
      mu3ptmax = pt;
    }
    else if(pt >= mu4ptmax) mu4ptmax = pt;
  }

  if(mu4ptmax >= 0.){
    cut1 = mu1ptmax;
    cut2 = mu2ptmax;
    cut3 = mu3ptmax;
    cut4 = mu4ptmax;
  }

  return;
}

void L1AlgoFactory::SingleEGPt(Float_t& cut, Bool_t isIsolated , Bool_t isER) {

  if(nEGs < 1) return;

  Float_t ptmax = -10.;

  for(UInt_t ue=0; ue < nEGs; ue++) {
    Int_t bx = egBx.at(ue);  
    if(bx != 0) continue;
    if(isIsolated && !egIso.at(ue)) continue;
    Float_t eta = egEta.at(ue);
    if(fabs(eta) > 2.1 && isER) continue;  // eta = 5 - 16

    Float_t pt = egEt.at(ue);    // the rank of the electron
    if(pt >= ptmax) ptmax = pt;
  }

  cut = ptmax;

  return;
}

void L1AlgoFactory::DoubleEGPt(Float_t& cut1, Float_t& cut2, Bool_t isIsolated, Bool_t isER ) {

  if(nEGs < 2) return;

  Float_t ele1ptmax = -10.;
  Float_t ele2ptmax = -10.;

  Float_t ele1Phimax = -1000.;
  Float_t ele1Etamax = -1000.;

  Bool_t EG1_ER = false;
  Bool_t EG2_ER = false;

  Bool_t EG1_isol = false;
  Bool_t EG2_isol = false;

  for(UInt_t ue=0; ue < nEGs; ue++) {               
    Int_t bx = egBx.at(ue);  
    if(bx != 0) continue;
    Float_t pt = egEt.at(ue);
    Float_t phi = egPhi.at(ue);
    Float_t eta = egEta.at(ue);

    if(fabs(pt-ele1ptmax) < 0.001 && fabs(phi-ele1Phimax) < 0.001 && fabs(eta-ele1Etamax) < 0.001) continue; //to avoid double counting in noniso/relaxiso lists

    if(pt >= ele1ptmax)
      {
	ele2ptmax = ele1ptmax;
	EG2_ER = EG1_ER;
	EG2_isol = EG1_isol;
	ele1ptmax = pt;
	ele1Phimax = phi;
	ele1Etamax = eta;
	if(fabs(eta) < 2.1) EG1_ER = true;
	else EG1_ER = false;
	EG1_isol = egIso.at(ue);
      }
    else if(pt >= ele2ptmax){
      ele2ptmax = pt;
      if(fabs(eta) < 2.1) EG2_ER = true;
      else EG2_ER = false;
      EG1_isol = egIso.at(ue);
    }
  }

  if(isER && (!EG1_ER || !EG2_ER)) return;
  if(isIsolated && (!EG1_isol && !EG2_isol)) return;

  if(ele2ptmax >= 0.){
    cut1 = ele1ptmax;
    cut2 = ele2ptmax;
  }

  return;
}

void L1AlgoFactory::TripleEGPt(Float_t& cut1, Float_t& cut2, Float_t& cut3 ) {

  if(nEGs < 2) return;

  Float_t ele1ptmax = -10.;
  Float_t ele2ptmax = -10.;
  Float_t ele3ptmax = -10.;

  Float_t ele1Phimax = -1000.;
  Float_t ele1Etamax = -1000.;

  Float_t ele2Phimax = -1000.;
  Float_t ele2Etamax = -1000.;

  for (UInt_t ue=0; ue < nEGs; ue++) {
    Int_t bx = egBx.at(ue);  
    if(bx != 0) continue;
    Float_t pt = egEt.at(ue);
    Float_t phi = egPhi.at(ue);
    Float_t eta = egEta.at(ue);

    if(fabs(pt-ele1ptmax) < 0.001 && fabs(phi-ele1Phimax) < 0.001 && fabs(eta-ele1Etamax) < 0.001) continue; //to avoid double counting in noniso/relaxiso lists
    if(fabs(pt-ele2ptmax) < 0.001 && fabs(phi-ele2Phimax) < 0.001 && fabs(eta-ele2Etamax) < 0.001) continue; //to avoid double counting in noniso/relaxiso lists

    if(pt >= ele1ptmax)
      {
	ele3ptmax = ele2ptmax;
	ele2ptmax = ele1ptmax;
	ele1ptmax = pt;
 	ele1Phimax = phi;
	ele1Etamax = eta;
      }
    else if(pt >= ele2ptmax){
      ele3ptmax = ele2ptmax;
      ele2ptmax = pt;
      ele2Phimax = phi;
      ele2Etamax = eta;
    }
    else if(pt >= ele3ptmax) ele3ptmax = pt;
  }

  if(ele3ptmax >= 0.){
    cut1 = ele1ptmax;
    cut2 = ele2ptmax;
    cut3 = ele3ptmax;
  }

  return;
}

void L1AlgoFactory::SingleJetPt(Float_t& cut, Bool_t isCentral) {

  Float_t ptmax = -10.;
  for(UInt_t ue=0; ue < nJets; ue++) {
    Int_t bx = jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(jetEta.at(ue)) > 3. ? true : false;
    if(isCentral && isFwdJet) continue;

    Float_t pt = jetEt.at(ue);
    if(pt >= ptmax) ptmax = pt;
  }

  cut = ptmax;
  return;
}

void L1AlgoFactory::DoubleJetPt(Float_t& cut1, Float_t& cut2, Bool_t isCentral ) {

  Float_t maxpt1 = -10.;
  Float_t maxpt2 = -10.;

  if(nJets < 2) return;

  for(UInt_t ue=0; ue < nJets; ue++) {
    Int_t bx = jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(jetEta.at(ue)) > 3. ? true : false;
    if(isCentral && isFwdJet) continue;

    Float_t pt = jetEt.at(ue);

    if(pt >= maxpt1)
      {
	maxpt2 = maxpt1;
	maxpt1 = pt;
      }
    else if(pt >= maxpt2) maxpt2 = pt;
  }

  if(maxpt2 >= 0.){
    cut1 = maxpt1;
    cut2 = maxpt2;
  }

  return;
}

void L1AlgoFactory::DoubleJet_Eta1p7_deltaEta4Pt(Float_t& cut1, Float_t& cut2 ) {

  if(nJets < 2) return;

  Float_t maxpt1 = -10.;
  Float_t maxpt2 = -10.;
  Bool_t corr = false;
  std::vector<std::pair<Float_t,Float_t> > jetPairs;

  for (UInt_t ue=0; ue < nJets; ue++) {
    Int_t bx = jetBx.at(ue);        		
    if(bx != 0) continue;
    Float_t pt = jetEt.at(ue);
    Float_t eta1 = jetEta.at(ue);
    if (fabs(eta1) > 1.7) continue;  // eta = 6 - 16

    for(UInt_t ve=0; ve < nJets; ve++) {
      if(ve == ue) continue;
      Int_t bx2 = jetBx.at(ve);        		
      if(bx2 != 0) continue;
      Float_t pt2 = jetEt.at(ve);
      Float_t eta2 = jetEta.at(ve);
      if (fabs(eta2) > 1.7) continue;  // eta = 6 - 16

      if(correlateInEta((int)eta1, (int)eta2, 4)){
	corr = true;
	jetPairs.push_back(std::pair<Float_t,Float_t>(pt,pt2));
      }

    }
  }

  if(corr){
    std::vector<std::pair<Float_t,Float_t> >::const_iterator jetPairIt  = jetPairs.begin();
    std::vector<std::pair<Float_t,Float_t> >::const_iterator jetPairEnd = jetPairs.end();
    for (; jetPairIt != jetPairEnd; ++jetPairIt) {
      Float_t pt1 = jetPairIt->first;
      Float_t pt2 = jetPairIt->second;
      
      if ( pt1 > maxpt1 || (fabs(maxpt1-pt1)<10E-2 && pt2>maxpt2) ) 
	{
	  maxpt1 = pt1;
	  maxpt2 = pt2;
	}
    }
  }

  cut1 = maxpt1;
  cut2 = maxpt2;

  return;
}

void L1AlgoFactory::DoubleTauJetEta2p17Pt(Float_t& cut1, Float_t& cut2, Bool_t isIsolated) {

  Float_t maxpt1 = -10.;
  Float_t maxpt2 = -10.;

  if(nTaus < 2) return;

  for(UInt_t ue=0; ue < nTaus; ue++) {
    Int_t bx = tauBx.at(ue);
    if(bx != 0) continue; 
    if(!isIsolated && !tauIso.at(ue)) continue;
    Float_t pt = tauEt.at(ue);
    Float_t eta = tauEta.at(ue);
    if(fabs(eta) > 2.17) continue;  // eta = 5 - 16

    if(pt >= maxpt1)
      {
	maxpt2 = maxpt1;
	maxpt1 = pt;
      }
    else if(pt >= maxpt2) maxpt2 = pt;
  }

  if(maxpt2 >= 0.){
    cut1 = maxpt1;
    cut2 = maxpt2;
  }

  return;
}

void L1AlgoFactory::TripleJetPt(Float_t& cut1, Float_t& cut2, Float_t& cut3, Bool_t isCentral) {

  Float_t jet1ptmax = -10.;
  Float_t jet2ptmax = -10.;
  Float_t jet3ptmax = -10.;

  if(nJets < 3) return;

  for (UInt_t ue=0; ue < nJets; ue++) {
    Int_t bx = jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(jetEta.at(ue)) > 3. ? true : false;
    if(isCentral && isFwdJet) continue;

    Float_t pt = jetEt.at(ue);

    if(pt >= jet1ptmax)
      {
	jet3ptmax = jet2ptmax;
	jet2ptmax = jet1ptmax;
	jet1ptmax = pt;
      }
    else if(pt >= jet2ptmax){
      jet3ptmax = jet2ptmax;
      jet2ptmax = pt;
    }
    else if(pt >= jet3ptmax) jet3ptmax = pt;
  }

  if(jet3ptmax >= 0.){
    cut1 = jet1ptmax;
    cut2 = jet2ptmax;
    cut3 = jet3ptmax;
  }

  return;
}

//For now, only usable in Menu mode
Bool_t L1AlgoFactory::TripleJet_VBF(Float_t jet1, Float_t jet2, Float_t jet3, Int_t jetclass ) {

  Bool_t jet=false;
  Bool_t jetf=false;

  Bool_t jetc1=false;
  Bool_t jetc2=false;
  Bool_t jetc3=false;

  Bool_t jetf1=false;           
  Bool_t jetf2=false;
  Bool_t jetf3=false;

  Int_t n1=0;
  Int_t n2=0;
  Int_t n3=0;

  Int_t f1=0;
  Int_t f2=0;
  Int_t f3=0;

  if(nJets < 3) return false;

  for (UInt_t ue=0; ue < nJets; ue++) {
    Int_t bx = jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(jetEta.at(ue)) > 3. ? true : false;
    Float_t pt = jetEt.at(ue);

    if (isFwdJet) {
      if(pt >= jet1) f1++;
      if(pt >= jet2) f2++;
      if(pt >= jet3) f3++;              
    } 
    else {
      if(pt >= jet1) n1++;
      if(pt >= jet2) n2++;
      if(pt >= jet3) n3++;
    }    
  }

  jet   = ( n1 >= 1 && n2 >= 2 && n3 >= 3 );
  jetf  = ( f1 >= 1 && f2 >= 2 && f3 >= 3 );

  jetc1 = ( n1 >= 1 && f2 >= 1 && f3 >= 2 );
  jetc2 = ( f1 >= 1 && n2 >= 1 && f3 >= 2 );
  jetc3 = ( f1 >= 1 && f2 >= 1 && n3 >= 2 );

  jetf1 = ( f1 >= 1 && n2 >= 1 && n3 >= 2 );
  jetf2 = ( n1 >= 1 && f2 >= 1 && n3 >= 2 );  
  jetf3 = ( n1 >= 1 && n2 >= 1 && f3 >= 2 );  

  if(jetclass == 1) return jet;
  else if(jetclass == 2) return jetf1;
  else if(jetclass == 3) return jetf2;
  else if(jetclass == 4) return jetf3;
  else if(jetclass == 5) return jetc1;
  else if(jetclass == 6) return jetc2;
  else if(jetclass == 7) return jetc3;
  else if(jetclass == 8) return jetf;

  return ( jet || jetf1 || jetf2 || jetf3 );
  //return ( jet || jetf1 || jetf2 || jetf3 || jetc1 || jetc2 || jetc3);
}

void L1AlgoFactory::QuadJetPt(Float_t& cut1, Float_t& cut2, Float_t& cut3, Float_t& cut4, Bool_t isCentral){

  Float_t jet1ptmax = -10.;
  Float_t jet2ptmax = -10.;
  Float_t jet3ptmax = -10.;
  Float_t jet4ptmax = -10.;

  if(nJets < 4) return;

  for (UInt_t ue=0; ue < nJets; ue++) {
    Int_t bx = jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(jetEta.at(ue)) > 3. ? true : false;
    if(isCentral && isFwdJet) continue;

    Float_t pt = jetEt.at(ue);

    if(pt >= jet1ptmax)
      {
	jet4ptmax = jet3ptmax;
	jet3ptmax = jet2ptmax;
	jet2ptmax = jet1ptmax;
	jet1ptmax = pt;
      }
    else if(pt >= jet2ptmax){
      jet4ptmax = jet3ptmax;
      jet3ptmax = jet2ptmax;
      jet2ptmax = pt;
    }
    else if(pt >= jet3ptmax){
      jet4ptmax = jet3ptmax;
      jet3ptmax = pt;
    }
    else if(pt >= jet4ptmax) jet4ptmax = pt;
  }

  if(jet4ptmax >= 0.){
    cut1 = jet1ptmax;
    cut2 = jet2ptmax;
    cut3 = jet3ptmax;
    cut4 = jet4ptmax;
  }

  return;
}

void L1AlgoFactory::Mu_EGPt(Float_t& mucut, Float_t& EGcut, Bool_t isIsolated, Int_t qualmin) {

  Float_t muptmax = -10.;

  for(UInt_t imu=0; imu < nMuons; imu++) {   
    Int_t bx = muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = muonEt.at(imu);			
    if(pt >= muptmax) muptmax = pt;
  }

  Float_t eleptmax = -10.;

  for(UInt_t ue=0; ue < nEGs; ue++) {
    Int_t bx = egBx.at(ue);        		
    if(bx != 0) continue;
    if(isIsolated && egIso.at(ue)) continue;
    Float_t pt = egEt.at(ue);    // the rank of the electron
    if(pt >= eleptmax) eleptmax = pt;
  }

  if(muptmax >= 0. && eleptmax >= 0.){
    mucut = muptmax;
    EGcut = eleptmax;
  }

  return;
}


void L1AlgoFactory::DoubleMu_EGPt(Float_t& mucut, Float_t& EGcut, Bool_t isMuHighQual ) {

  Float_t muptmax = -10.;
  Float_t second_muptmax = -10.;
  Float_t EGmax = -10.;

  if(nMuons < 2) return;

  for (UInt_t imu=0; imu < nMuons; imu++) {
    Int_t bx = muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = muonEt.at(imu);			
    if(pt >= muptmax){
      second_muptmax = muptmax;
      muptmax = pt;
    }
    else if(pt >= second_muptmax) second_muptmax = pt;
  }

  for (UInt_t ue=0; ue < nEGs; ue++) {
    Int_t bx = egBx.at(ue);        		
    if(bx != 0) continue;
    Float_t pt = egEt.at(ue);
    if (pt >= EGmax) EGmax = pt;
  }  // end loop over EM objects

  if(second_muptmax >= 0.){
    mucut = second_muptmax;
    EGcut = EGmax;
  }

  return;
}

void L1AlgoFactory::Mu_DoubleEGPt(Float_t& mucut, Float_t& EGcut ) {

  Float_t muptmax    = -10.;
  Float_t eleptmax1  = -10.;
  Float_t eleptmax2  = -10.;
  Float_t ele1Phimax = -1000.;
  Float_t ele1Etamax = -1000.;

  for (UInt_t imu=0; imu < nMuons; imu++) {
    Int_t bx = muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = muonEt.at(imu);			
    if(pt >= muptmax) muptmax = pt; 
  }

  if(nEGs < 2) return;

  for (UInt_t ue=0; ue < nEGs; ue++) {
    Int_t bx = egBx.at(ue);        		
    if(bx != 0) continue;
    Float_t pt  = egEt.at(ue);
    Float_t phi = egPhi.at(ue);
    Float_t eta = egEta.at(ue);

    if(fabs(pt-eleptmax1) < 0.001 && fabs(phi-ele1Phimax) < 0.001 && fabs(eta-ele1Etamax) < 0.001) continue; //to avoid double counting in noniso/relaxiso lists

    if(pt >= eleptmax1){
      eleptmax2 = eleptmax1;
      eleptmax1 = pt;
      ele1Phimax = phi;
      ele1Etamax = eta;
    }
    else if(pt >= eleptmax2) eleptmax2 = pt;
  }

  if(muptmax >= 0. && eleptmax2 >= 0.){
    mucut = muptmax;
    EGcut = eleptmax2;
  }

  return;
}

void L1AlgoFactory::Mu_DoubleJetCentralPt(Float_t& mucut, Float_t& jetcut) {

  Float_t muptmax = -10.;
  Float_t jetptmax1 = -10.;
  Float_t jetptmax2 = -10.;

  for (UInt_t imu=0; imu < nMuons; imu++) {
    Int_t bx = muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = muonEt.at(imu);			
    if (pt >= muptmax) muptmax = pt;
  }

  if(nJets < 2) return;

  for (UInt_t ue=0; ue < nJets; ue++) {
    Int_t bx = jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(jetEta.at(ue)) > 3. ? true : false;
    if(isFwdJet) continue;

    Float_t pt = jetEt.at(ue);
    if(pt >= jetptmax1){
      jetptmax2 = jetptmax1;
      jetptmax1 = pt;
    }
    else if(pt >= jetptmax2) jetptmax2 = pt;
  }

  if(muptmax >= 0. && jetptmax2 >= 0.){
    mucut = muptmax;
    jetcut = jetptmax2;
  }

  return;
}

void L1AlgoFactory::Muer_JetCentralPt(Float_t& mucut, Float_t& jetcut) {

  Float_t muptmax = -10.;
  Float_t jetptmax = -10.;

  for (UInt_t imu=0; imu < nMuons; imu++) {
    Int_t bx = muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = muonEt.at(imu);
    Float_t eta = muonEta.at(imu); 
    if(fabs(eta) > 2.1) continue;

    if (pt >= muptmax) muptmax = pt;
  }

  for (UInt_t ue=0; ue < nJets; ue++) {
    Int_t bx = jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(jetEta.at(ue)) > 3. ? true : false;
    if(isFwdJet) continue;

    Float_t pt = jetEt.at(ue);
    if (pt >= jetptmax) jetptmax = pt;
  }

  if(muptmax >= 0. && jetptmax >= 0.){
    mucut = muptmax;
    jetcut = jetptmax;
  }

  return;
}

void L1AlgoFactory::Mu_JetCentral_deltaPt(Float_t& mucut, Float_t& jetcut) {

  Float_t muptmax = -10.;
  Float_t jetptmax = -10.;
  Bool_t correlate = false;

  for (UInt_t imu=0; imu < nMuons; imu++) {
    Int_t bx = muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = muonEt.at(imu);
    Float_t eta = muonEta.at(imu); 
    Float_t phi = muonPhi.at(imu); 

    for (UInt_t ue=0; ue < nJets; ue++) {
    Int_t bxj = jetBx.at(ue);        		
    if(bxj != 0) continue;
    Bool_t isFwdJet = fabs(jetEta.at(ue)) > 3. ? true : false;
    if(isFwdJet) continue;

    Float_t ptj = jetEt.at(ue);
    Float_t phijet = jetPhi.at(ue);
    Float_t etajet = jetEta.at(ue);

      if(pt < mucut || ptj < jetcut) continue;

      if(fabs(phi-phijet) < 0.4 && fabs(eta-etajet) < 0.4){
	correlate = true;
	if(pt >= muptmax) muptmax = pt;
	if(ptj >= jetptmax) jetptmax = ptj;
      }
    }
  }

  if(correlate){
    mucut = muptmax;
    jetcut = jetptmax;
  }

  return;
}

void L1AlgoFactory::EG_FwdJetPt(Float_t& EGcut, Float_t& FWcut) {

  Float_t eleptmax = -10.;
  Float_t jetptmax = -10.;

  for (UInt_t ue=0; ue < nEGs; ue++) {
    Int_t bx = egBx.at(ue);        		
    if(bx != 0) continue;
    Float_t pt  = egEt.at(ue);
    if(pt >= eleptmax) eleptmax = pt;
  }

  for (UInt_t ue=0; ue < nJets; ue++) {        
    Int_t bx = jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(jetEta.at(ue)) > 3. ? true : false;
    if(!isFwdJet) continue;
    Float_t pt = jetEt.at(ue);
    if (pt >= jetptmax) jetptmax = pt;
  }

  if(eleptmax >= 0. && jetptmax >= 0.){
    EGcut = eleptmax;
    FWcut = jetptmax;
  }

  return;
}

void L1AlgoFactory::EG_DoubleJetCentralPt(Float_t& EGcut, Float_t& jetcut) {

  Float_t eleptmax = -10.;
  Float_t jetptmax1 = -10.;
  Float_t jetptmax2 = -10.;

  for (UInt_t ue=0; ue < nEGs; ue++) {
    Int_t bx = egBx.at(ue);        		
    if(bx != 0) continue;
    Float_t pt  = egEt.at(ue);
    if(pt >= eleptmax) eleptmax = pt; 
  }

  if(nJets < 2) return;

  for (UInt_t ue=0; ue < nJets; ue++) {
    Int_t bx = jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(jetEta.at(ue)) > 3. ? true : false;
    if(isFwdJet) continue;
    Float_t pt = jetEt.at(ue);

    if(pt >= jetptmax1){
      jetptmax2 = jetptmax1;
      jetptmax1 = pt;
    }
    else if(pt >= jetptmax2) jetptmax2 = pt;
  }

  if(eleptmax >= 0. && jetptmax2 >= 0.){
    EGcut = eleptmax;
    jetcut = jetptmax2;
  }

  return;
}

void L1AlgoFactory::EGer_TripleJetCentralPt(Float_t& EGcut, Float_t& jetcut) {

  Float_t eleptmax = -10.;
  Float_t elemaxeta = -10.;
  Float_t jetptmax1 = -10.;
  Float_t jetptmax2 = -10.;
  Float_t jetptmax3 = -10.;

  for (UInt_t ue=0; ue < nEGs; ue++){
    Int_t bx = egBx.at(ue);        		
    if(bx != 0) continue;
    Float_t pt  = egEt.at(ue);
    Float_t eta = egEta.at(ue);
    if(fabs(eta) > 2.1) continue;  // eta = 5 - 16
    if(pt >= eleptmax){
      eleptmax = pt; 
      elemaxeta = egEta.at(ue);
    }
  }  // end loop over EM objects

  if(nJets < 3) return;

  for (UInt_t ue=0; ue < nJets; ue++) {
    Int_t bx = jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(jetEta.at(ue)) > 3. ? true : false;
    if(isFwdJet) continue;
    Float_t pt = jetEt.at(ue);
    Float_t jeteta = jetEta.at(ue);

    if(jeteta == elemaxeta) continue;   //both are binned with the same binning

    if(pt >= jetptmax1){
      jetptmax3 = jetptmax2;
      jetptmax2 = jetptmax1;
      jetptmax1 = pt;
    }
    else if(pt >= jetptmax2){
      jetptmax3 = jetptmax2;
      jetptmax2 = pt;
    }
    else if(pt >= jetptmax3) jetptmax3 = pt;

  }

  if(eleptmax >= 0. && jetptmax3 >= 0.){
    EGcut = eleptmax;
    jetcut = jetptmax3;
  }

  return;
}

void L1AlgoFactory::Muer_TauJetEta2p17Pt(Float_t& mucut, Float_t& taucut, Bool_t isIsolated) {

  Float_t maxptmu  = -10.;
  Float_t maxpttau = -10.;

  if(nMuons < 1) return;
  for (UInt_t imu=0; imu < nMuons; imu++) {
    Int_t bx = muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = muonEt.at(imu);
    Float_t eta = muonEta.at(imu); 
    if(fabs(eta) > 2.1) continue;
    if(pt >= maxptmu) maxptmu = pt;
  }

  for(UInt_t ue=0; ue < nTaus; ue++) {
    Int_t bx = tauBx.at(ue);        		
    if(bx != 0) continue; 
    if(!isIsolated && !tauIso.at(ue)) continue;

    Float_t pt = tauEt.at(ue);
    Float_t eta = tauEta.at(ue);
    if(fabs(eta) > 2.17) continue;

    if(pt >= maxpttau) maxpttau = pt;
  }

  if(maxptmu >= 0.){
    mucut  = maxptmu;
    taucut = maxpttau;
  }

  return;

}

void L1AlgoFactory::IsoEGer_TauJetEta2p17Pt(Float_t& egcut, Float_t& taucut) {

  Float_t eleptmax  = -10.;
  Float_t eleetamax = -999.;
  Float_t maxpttau  = -10.;

  if(nEGs < 1) return;
  for (UInt_t ue=0; ue < nEGs; ue++) {
    Int_t bx = egBx.at(ue);        		
    if(bx != 0) continue;
    if(!egIso.at(ue)) continue;
    Float_t pt  = egEt.at(ue);
    Float_t eta = egEta.at(ue);
    if(fabs(eta) > 2.1) continue;  // eta = 5 - 16
    if(pt >= eleptmax){
      eleptmax = pt;
      eleetamax = eta;
    }
  }

  for(UInt_t ue=0; ue < nTaus; ue++) {
    Int_t bx = tauBx.at(ue);        		
    if(bx != 0) continue; 
    Float_t pt = tauEt.at(ue);
    Float_t eta = tauEta.at(ue);
    if(fabs(eta) > 2.17) continue;  // eta = 5 - 16

    if(fabs(eta-eleetamax) < 0.2) continue;

    if(pt >= maxpttau) maxpttau = pt;
  }

  if(eleptmax >= 0.){
    egcut  = eleptmax;
    taucut = maxpttau;
  }

  return;
}

void L1AlgoFactory::QuadJetCentral_TauJetPt(Float_t& jetcut, Float_t& taucut){

  Float_t jet1ptmax = -10.;
  Float_t jet2ptmax = -10.;
  Float_t jet3ptmax = -10.;
  Float_t jet4ptmax = -10.;
  Float_t maxpttau  = -10.;

  if(nJets < 4) return;
  if(nTaus < 1) return;

  for (UInt_t ue=0; ue < nJets; ue++) {
    Int_t bx = jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(jetEta.at(ue)) > 3. ? true : false;
    if(isFwdJet) continue;
    Float_t pt = jetEt.at(ue);

    if(pt >= jet1ptmax)
      {
	jet4ptmax = jet3ptmax;
	jet3ptmax = jet2ptmax;
	jet2ptmax = jet1ptmax;
	jet1ptmax = pt;
      }
    else if(pt >= jet2ptmax){
      jet4ptmax = jet3ptmax;
      jet3ptmax = jet2ptmax;
      jet2ptmax = pt;
    }
    else if(pt >= jet3ptmax){
      jet4ptmax = jet3ptmax;
      jet3ptmax = pt;
    }
    else if(pt >= jet4ptmax) jet4ptmax = pt;
  }

  for(UInt_t ue=0; ue < nTaus; ue++) {
    Int_t bx = tauBx.at(ue);        		
    if(bx != 0) continue; 
    Float_t pt = tauEt.at(ue);
    if(pt >= maxpttau) maxpttau = pt;
  }

  if(jet4ptmax >= 0. && maxpttau >= 0.){
    jetcut = jet4ptmax;
    taucut = maxpttau;
  }

  return;
}


void L1AlgoFactory::ETMVal(Float_t& ETMcut ) {

  Float_t TheETM = -10;
  if(sumBx[1]==0) TheETM =sumEt[1];
  ETMcut = TheETM;
  return;
}

void L1AlgoFactory::HTTVal(Float_t& HTTcut) {

  Float_t TheHTT = -10;
  if(sumBx[2]==0) TheHTT =sumEt[2];
  HTTcut = TheHTT;
  return;
}

void L1AlgoFactory::HTMVal(Float_t& HTMcut) {

  Float_t TheHTM = -10;
  if (sumBx[3]==0) TheHTM = sumEt[3];
  HTMcut = TheHTM;
  return;
}

void L1AlgoFactory::ETTVal(Float_t& ETTcut) {

  Float_t TheETT = -10;
  if(sumBx[0]==0) TheETT = sumEt[0];
  ETTcut = TheETT;
  return;
}

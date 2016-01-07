#ifndef  __L1ALGOFACTORY_INC__
#define  __L1ALGOFACTORY_INC__

#include "L1Ntuple.h"
#include<iostream>

enum EtSumType { ETT, HTT, ETM, HTM }; // Base on "DataFormats/L1Trigger/interface/EtSum.h"
const size_t PHIBINS = 18;
const Double_t PHIBIN[] = {10,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350};
const size_t ETABINS = 23;
const Double_t ETABIN[] = {-5.,-4.5,-4.,-3.5,-3.,-2.172,-1.74,-1.392,-1.044,-0.696,-0.348,0.,0.348,0.696,1.044,1.392,1.74,2.172,3.,3.5,4.,4.5,5.};

const size_t ETAMUBINS = 65;
const Double_t ETAMU[] = { -2.45,-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45 };

class L1AlgoFactory: public L1Ntuple{
 public:
  L1AlgoFactory():jetCentFwd(3.0)
  {};
  //L1AlgoFactory(TTree *tree);

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
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Hard coded cut value ~~~~~
 float jetCentFwd;
 

};

#endif   // ----- #ifndef __L1ALGOFACTORY_INC__  -----

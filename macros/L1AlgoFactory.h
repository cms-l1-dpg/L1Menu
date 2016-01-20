#ifndef  __L1ALGOFACTORY_INC__
#define  __L1ALGOFACTORY_INC__

#include "L1Ntuple.h"
#include <iostream>
#include <cassert>

enum EtSumType { ETT, HTT, ETM, HTM }; // Base on "DataFormats/L1Trigger/interface/EtSum.h"

class L1AlgoFactory: public L1Ntuple{
 public:
  L1AlgoFactory():jetCentFwd(3.0),muonER(2.1),eleER(2.1),tauER(2.17),
  MuJetCordPhi(0.4), MuJetCordEta(0.4), MuOpenJetCordPhi(3.0), MuMudPhi(1.0), Onia2015ER(1.6)
  {};
  //L1AlgoFactory(TTree *tree);

  void SingleMuPt(Float_t& ptcut, Bool_t isER, Int_t qualmin=4);
  void DoubleMuPt(Float_t& mu1pt, Float_t& mu2pt, Bool_t isHighQual = false, Bool_t isER = false);
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
  void SingleTauPt(Float_t& cut, Bool_t isCentral);

  void Mu_EGPt(Float_t& mucut, Float_t& EGcut, Bool_t isIsolated = false, Int_t qualmin=4);
  void DoubleMu_EGPt(Float_t& mucut, Float_t& EGcut, Bool_t isMuHighQual = false );
  void Mu_DoubleEGPt(Float_t& mucut, Float_t& EGcut );

  void Muer_JetCentralPt(Float_t& mucut, Float_t& jetcut);
  void Mu_JetCentral_deltaPt(Float_t& mucut, Float_t& jetcut);
  void Mu_DoubleJetCentralPt(Float_t& mucut, Float_t& jetcut);

  void EG_FwdJetPt(Float_t& EGcut, Float_t& FWcut);
  void EG_DoubleJetCentralPt(Float_t& EGcut, Float_t& jetcut);
  void EGer_TripleJetCentralPt(Float_t& EGcut, Float_t& jetcut);
  void IsoEGer_TauJetEta2p17Pt(Float_t& egcut, Float_t& taucut);

  void QuadJetCentral_TauJetPt(Float_t& jetcut, Float_t& taucut);

  int GetSumEtIdx(EtSumType type);
  void ETMVal(Float_t& ETMcut);
  void HTTVal(Float_t& HTTcut);
  void HTMVal(Float_t& HTMcut);
  void ETTVal(Float_t& ETTcut);


  Bool_t SingleMu(Float_t ptcut, Bool_t isER, Int_t qualmin=4);
  Bool_t DoubleMu(Float_t mu1pt, Float_t mu2pt, Bool_t isHighQual = false, Bool_t isER = false);
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

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MultiMuon ~~~~~
  Bool_t DoubleMuXOpen(Float_t mu1pt);
  void DoubleMuXOpenPt(Float_t& cut);
  Bool_t Onia2015(Float_t mu1pt, Float_t mu2pt, Bool_t isER, Bool_t isOS, Int_t delta);
  void Onia2015Pt(Float_t& ptcut1, Float_t& ptcut2, Bool_t isER, Bool_t isOS, Int_t delta);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Cross ~~~~~
  Bool_t Mu_HTT(Float_t mucut, Float_t HTcut);
  void Mu_HTTPt(Float_t& mucut, Float_t& HTcut );
  Bool_t Muer_ETM(Float_t mucut, Float_t ETMcut);
  void Muer_ETMPt(Float_t& mucut, Float_t& ETMcut );
  Bool_t SingleEG_Eta2p1_HTT(Float_t egcut, Float_t HTTcut, Bool_t isIsolated = false);
  void SingleEG_Eta2p1_HTTPt(Float_t& egcut, Float_t& HTTcut, Bool_t isIsolated = false);
  Bool_t Muer_TauJetEta2p17(Float_t mucut, Float_t taucut, Bool_t isIsolated = false);
  void Muer_TauJetEta2p17Pt(Float_t& mucut, Float_t& taucut, Bool_t isIsolated = false);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MultiCross ~~~~~
  Bool_t DoubleJetCentral_ETM(Float_t jetcut1, Float_t jetcut2, Float_t ETMcut);
  void DoubleJetCentral_ETMPt(Float_t& jetcut1, Float_t& jetcut2, Float_t& ETMcut);
  Bool_t DoubleEG_HT(Float_t EGcut, Float_t HTcut);
  void DoubleEG_HTPt(Float_t& EGcut, Float_t& HTcut);
  Bool_t Jet_MuOpen_Mu_dPhiMuMu1(Float_t jetcut, Float_t mucut);
  void Jet_MuOpen_Mu_dPhiMuMu1Pt(Float_t& jetcut, Float_t& mucut);
  Bool_t Jet_MuOpen_EG_dPhiMuEG1(Float_t jetcut, Float_t egcut);
  void Jet_MuOpen_EG_dPhiMuEG1Pt(Float_t& jetcut, Float_t& egcut);

 private:
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Hard coded cut value ~~~~~
 float jetCentFwd;
 float muonER;
 float eleER;
 float tauER;
 float MuJetCordPhi;
 float MuJetCordEta;
 float MuOpenJetCordPhi;
 float MuMudPhi;
 float Onia2015ER;
 
 bool PassMuonQual(int imu, bool isMuHighQual=true) const;

};

#endif   // ----- #ifndef __L1ALGOFACTORY_INC__  -----

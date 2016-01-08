// ===========================================================================
// 
//       Filename:  L1AlgoFactory.C
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  01/07/2016 01:41:46 PM
//       Compiler:  g++ -std=c++11
// 
//         Author:  Zhenbin Wu (benwu)
//          Email:  zhenbin.wu@gmail.com
//        Company:  UIC, CMS@LPC, CDF@FNAL
// 
// ===========================================================================

#include "L1AlgoFactory.h"

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

inline Bool_t L1AlgoFactory::correlateInEta(Int_t mueta, Int_t jeteta, Int_t delta) {
  return fabs(mueta-jeteta) < 1 + delta;
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

  if(upgrade_->nMuons < 1) return;

  Float_t ptmax = -10.;

  for(UInt_t imu=0; imu < upgrade_->nMuons; imu++) {
    Int_t bx = upgrade_->muonBx.at(imu);
    if(bx != 0) continue;
    Float_t eta = upgrade_->muonEta.at(imu);        
    if(fabs(eta) > muonER && isER) continue;
    Float_t pt = upgrade_->muonEt.at(imu);                       

    if(pt >= ptmax) ptmax = pt;
  }

  ptcut = ptmax;
  return;
}

void L1AlgoFactory::DoubleMuPt(Float_t& cut1, Float_t& cut2, Bool_t isMuHighQual, Bool_t isER) {

  Float_t mu1ptmax = -10.;
  Float_t mu2ptmax = -10.;

  if(upgrade_->nMuons < 2) return;

  for (UInt_t imu=0; imu < upgrade_->nMuons; imu++) {
    Int_t bx = upgrade_->muonBx.at(imu);
    if(bx != 0) continue;
    Float_t pt = upgrade_->muonEt.at(imu);			
    Float_t eta = upgrade_->muonEta.at(imu);
    if(isER && fabs(eta) > muonER) continue;
    if(!PassMuonQual(imu, isMuHighQual)) continue;

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

  if(upgrade_->nMuons < 2) return;

  Float_t maxpt1 = -10.;
  Float_t maxpt2 = -10.;
  Float_t corr = false;

  std::vector<std::pair<Float_t,Float_t> > muonPairs;

  for (UInt_t imu=0; imu < upgrade_->nMuons; imu++) {
    Int_t bx = upgrade_->muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = upgrade_->muonEt.at(imu);			
    Float_t eta = upgrade_->muonEta.at(imu);
    if (fabs(eta) > 2.1) continue;
    Int_t ieta1 = etaMuIdx(eta);

    for (UInt_t imu2=0; imu2 < upgrade_->nMuons; imu2++) {
      if (imu2 == imu) continue;
      Int_t bx2 = upgrade_->muonBx.at(imu2);		
      if(bx2 != 0) continue;
      Float_t pt2 = upgrade_->muonEt.at(imu2);			
      Float_t eta2 = upgrade_->muonEta.at(imu2);
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

  if(upgrade_->nMuons < 2) return;

  Float_t maxpt1 = -10.;
  Float_t maxpt2 = -10.;
  Float_t corr = false;

  std::vector<std::pair<Float_t,Float_t> > muonPairs;

  for (UInt_t imu=0; imu < upgrade_->nMuons; imu++) {
    Int_t bx = upgrade_->muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = upgrade_->muonEt.at(imu);			
    Float_t eta = upgrade_->muonEta.at(imu);
    if(isER && fabs(eta) > 1.6) continue;
    Int_t ieta1 = etaMuIdx(eta);
    Int_t charge1 = upgrade_->muonChg.at(imu);

    for (UInt_t imu2=0; imu2 < upgrade_->nMuons; imu2++) {
      if (imu2 == imu) continue;
      Int_t bx2 = upgrade_->muonBx.at(imu2);		
      if(bx2 != 0) continue;
      Float_t pt2 = upgrade_->muonEt.at(imu2);			
      Float_t eta2 = upgrade_->muonEta.at(imu2);
      if(isER && fabs(eta2) > 1.6) continue;
      Int_t ieta2 = etaMuIdx(eta2);
      Int_t charge2 = upgrade_->muonChg.at(imu2);

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

  if(upgrade_->nMuons < 3) return;

  for (UInt_t imu=0; imu < upgrade_->nMuons; imu++) {
    Int_t bx = upgrade_->muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = upgrade_->muonEt.at(imu);			

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

  if(upgrade_->nMuons < 4) return;

  for (UInt_t imu=0; imu < upgrade_->nMuons; imu++) {
    Int_t bx = upgrade_->muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = upgrade_->muonEt.at(imu);			

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

  if(upgrade_->nEGs < 1) return;

  Float_t ptmax = -10.;

  for(UInt_t ue=0; ue < upgrade_->nEGs; ue++) {
    Int_t bx = upgrade_->egBx.at(ue);  
    if(bx != 0) continue;
    if(isIsolated && !upgrade_->egIso.at(ue)) continue;
    Float_t eta = upgrade_->egEta.at(ue);
    if(fabs(eta) > eleER && isER) continue;  // eta = 5 - 16

    Float_t pt = upgrade_->egEt.at(ue);    // the rank of the electron
    if(pt >= ptmax) ptmax = pt;
  }

  cut = ptmax;

  return;
}

void L1AlgoFactory::DoubleEGPt(Float_t& cut1, Float_t& cut2, Bool_t isIsolated, Bool_t isER ) {

  if(upgrade_->nEGs < 2) return;

  Float_t ele1ptmax = -10.;
  Float_t ele2ptmax = -10.;

  Float_t ele1Phimax = -1000.;
  Float_t ele1Etamax = -1000.;

  Bool_t EG1_ER = false;
  Bool_t EG2_ER = false;

  Bool_t EG1_isol = false;
  Bool_t EG2_isol = false;

  for(UInt_t ue=0; ue < upgrade_->nEGs; ue++) {               
    Int_t bx = upgrade_->egBx.at(ue);  
    if(bx != 0) continue;
    Float_t pt = upgrade_->egEt.at(ue);
    Float_t phi = upgrade_->egPhi.at(ue);
    Float_t eta = upgrade_->egEta.at(ue);

    if(fabs(pt-ele1ptmax) < 0.001 && fabs(phi-ele1Phimax) < 0.001 && fabs(eta-ele1Etamax) < 0.001) continue; //to avoid double counting in noniso/relaxiso lists

    if(pt >= ele1ptmax)
    {
      ele2ptmax = ele1ptmax;
      EG2_ER = EG1_ER;
      EG2_isol = EG1_isol;
      ele1ptmax = pt;
      ele1Phimax = phi;
      ele1Etamax = eta;
      if(fabs(eta) < eleER) EG1_ER = true;
      else EG1_ER = false;
      EG1_isol = upgrade_->egIso.at(ue);
    }
    else if(pt >= ele2ptmax){
      ele2ptmax = pt;
      if(fabs(eta) < eleER) EG2_ER = true;
      else EG2_ER = false;
      EG1_isol = upgrade_->egIso.at(ue);
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

  if(upgrade_->nEGs < 2) return;

  Float_t ele1ptmax = -10.;
  Float_t ele2ptmax = -10.;
  Float_t ele3ptmax = -10.;

  Float_t ele1Phimax = -1000.;
  Float_t ele1Etamax = -1000.;

  Float_t ele2Phimax = -1000.;
  Float_t ele2Etamax = -1000.;

  for (UInt_t ue=0; ue < upgrade_->nEGs; ue++) {
    Int_t bx = upgrade_->egBx.at(ue);  
    if(bx != 0) continue;
    Float_t pt = upgrade_->egEt.at(ue);
    Float_t phi = upgrade_->egPhi.at(ue);
    Float_t eta = upgrade_->egEta.at(ue);

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
  for(UInt_t ue=0; ue < upgrade_->nJets; ue++) {
    Int_t bx = upgrade_->jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(upgrade_->jetEta.at(ue)) > jetCentFwd ? true : false;
    if(isCentral && isFwdJet) continue;

    Float_t pt = upgrade_->jetEt.at(ue);
    if(pt >= ptmax) ptmax = pt;
  }

  cut = ptmax;
  return;
}

void L1AlgoFactory::DoubleJetPt(Float_t& cut1, Float_t& cut2, Bool_t isCentral ) {

  Float_t maxpt1 = -10.;
  Float_t maxpt2 = -10.;

  if(upgrade_->nJets < 2) return;

  for(UInt_t ue=0; ue < upgrade_->nJets; ue++) {
    Int_t bx = upgrade_->jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(upgrade_->jetEta.at(ue)) > jetCentFwd ? true : false;
    if(isCentral && isFwdJet) continue;

    Float_t pt = upgrade_->jetEt.at(ue);

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

  if(upgrade_->nJets < 2) return;

  Float_t maxpt1 = -10.;
  Float_t maxpt2 = -10.;
  Bool_t corr = false;
  std::vector<std::pair<Float_t,Float_t> > jetPairs;

  for (UInt_t ue=0; ue < upgrade_->nJets; ue++) {
    Int_t bx = upgrade_->jetBx.at(ue);        		
    if(bx != 0) continue;
    Float_t pt = upgrade_->jetEt.at(ue);
    Float_t eta1 = upgrade_->jetEta.at(ue);
    if (fabs(eta1) > 1.7) continue;  // eta = 6 - 16

    for(UInt_t ve=0; ve < upgrade_->nJets; ve++) {
      if(ve == ue) continue;
      Int_t bx2 = upgrade_->jetBx.at(ve);        		
      if(bx2 != 0) continue;
      Float_t pt2 = upgrade_->jetEt.at(ve);
      Float_t eta2 = upgrade_->jetEta.at(ve);
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

  if(upgrade_->nTaus < 2) return;

  for(UInt_t ue=0; ue < upgrade_->nTaus; ue++) {
    Int_t bx = upgrade_->tauBx.at(ue);
    if(bx != 0) continue; 
    if(!isIsolated && !upgrade_->tauIso.at(ue)) continue;
    Float_t pt = upgrade_->tauEt.at(ue);
    Float_t eta = upgrade_->tauEta.at(ue);
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

  if(upgrade_->nJets < 3) return;

  for (UInt_t ue=0; ue < upgrade_->nJets; ue++) {
    Int_t bx = upgrade_->jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(upgrade_->jetEta.at(ue)) > 3. ? true : false;
    if(isCentral && isFwdJet) continue;

    Float_t pt = upgrade_->jetEt.at(ue);

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

  if(upgrade_->nJets < 3) return false;

  for (UInt_t ue=0; ue < upgrade_->nJets; ue++) {
    Int_t bx = upgrade_->jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(upgrade_->jetEta.at(ue)) > jetCentFwd ? true : false;
    Float_t pt = upgrade_->jetEt.at(ue);

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

  if(upgrade_->nJets < 4) return;

  for (UInt_t ue=0; ue < upgrade_->nJets; ue++) {
    Int_t bx = upgrade_->jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(upgrade_->jetEta.at(ue)) > jetCentFwd ? true : false;
    if(isCentral && isFwdJet) continue;

    Float_t pt = upgrade_->jetEt.at(ue);

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

  for(UInt_t imu=0; imu < upgrade_->nMuons; imu++) {   
    Int_t bx = upgrade_->muonBx.at(imu);		
    if(bx != 0) continue;
    if(!PassMuonQual(imu)) continue;
    Float_t pt = upgrade_->muonEt.at(imu);			
    if(pt >= muptmax) muptmax = pt;
  }

  Float_t eleptmax = -10.;

  for(UInt_t ue=0; ue < upgrade_->nEGs; ue++) {
    Int_t bx = upgrade_->egBx.at(ue);        		
    if(bx != 0) continue;
    if(isIsolated && upgrade_->egIso.at(ue)) continue;
    Float_t pt = upgrade_->egEt.at(ue);    // the rank of the electron
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

  if(upgrade_->nMuons < 2) return;

  for (UInt_t imu=0; imu < upgrade_->nMuons; imu++) {
    Int_t bx = upgrade_->muonBx.at(imu);		
    if(bx != 0) continue;
    if(!PassMuonQual(imu, isMuHighQual)) continue;
    Float_t pt = upgrade_->muonEt.at(imu);			
    if(pt >= muptmax){
      second_muptmax = muptmax;
      muptmax = pt;
    }
    else if(pt >= second_muptmax) second_muptmax = pt;
  }

  for (UInt_t ue=0; ue < upgrade_->nEGs; ue++) {
    Int_t bx = upgrade_->egBx.at(ue);        		
    if(bx != 0) continue;
    Float_t pt = upgrade_->egEt.at(ue);
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

  for (UInt_t imu=0; imu < upgrade_->nMuons; imu++) {
    Int_t bx = upgrade_->muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = upgrade_->muonEt.at(imu);			
    if(!PassMuonQual(imu)) continue;
    if(pt >= muptmax) muptmax = pt; 
  }

  if(upgrade_->nEGs < 2) return;

  for (UInt_t ue=0; ue < upgrade_->nEGs; ue++) {
    Int_t bx = upgrade_->egBx.at(ue);        		
    if(bx != 0) continue;
    Float_t pt  = upgrade_->egEt.at(ue);
    Float_t phi = upgrade_->egPhi.at(ue);
    Float_t eta = upgrade_->egEta.at(ue);

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

  for (UInt_t imu=0; imu < upgrade_->nMuons; imu++) {
    Int_t bx = upgrade_->muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = upgrade_->muonEt.at(imu);			
    if (pt >= muptmax) muptmax = pt;
  }

  if(upgrade_->nJets < 2) return;

  for (UInt_t ue=0; ue < upgrade_->nJets; ue++) {
    Int_t bx = upgrade_->jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(upgrade_->jetEta.at(ue)) > 3. ? true : false;
    if(isFwdJet) continue;

    Float_t pt = upgrade_->jetEt.at(ue);
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

  for (UInt_t imu=0; imu < upgrade_->nMuons; imu++) {
    Int_t bx = upgrade_->muonBx.at(imu);		
    if(bx != 0) continue;
    Float_t pt = upgrade_->muonEt.at(imu);
    Float_t eta = upgrade_->muonEta.at(imu); 
    if(fabs(eta) > 2.1) continue;

    if (pt >= muptmax) muptmax = pt;
  }

  for (UInt_t ue=0; ue < upgrade_->nJets; ue++) {
    Int_t bx = upgrade_->jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(upgrade_->jetEta.at(ue)) > 3. ? true : false;
    if(isFwdJet) continue;

    Float_t pt = upgrade_->jetEt.at(ue);
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

  for (UInt_t imu=0; imu < upgrade_->nMuons; imu++) {
    Int_t bx = upgrade_->muonBx.at(imu);		
    if(bx != 0) continue;
    if(!PassMuonQual(imu)) continue;
    Float_t pt = upgrade_->muonEt.at(imu);
    Float_t eta = upgrade_->muonEta.at(imu); 
    Float_t phi = upgrade_->muonPhi.at(imu); 

    for (UInt_t ue=0; ue < upgrade_->nJets; ue++) {
      Int_t bxj = upgrade_->jetBx.at(ue);        		
      if(bxj != 0) continue;
      Bool_t isFwdJet = fabs(upgrade_->jetEta.at(ue)) > jetCentFwd ? true : false;
      if(isFwdJet) continue;

      Float_t ptj = upgrade_->jetEt.at(ue);
      Float_t phijet = upgrade_->jetPhi.at(ue);
      Float_t etajet = upgrade_->jetEta.at(ue);

      if(pt < mucut || ptj < jetcut) continue;

      //MuJetCordPhi,  MuJetCordEta
      if(fabs(phi-phijet) < MuJetCordPhi && fabs(eta-etajet) < MuJetCordEta){
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

  for (UInt_t ue=0; ue < upgrade_->nEGs; ue++) {
    Int_t bx = upgrade_->egBx.at(ue);        		
    if(bx != 0) continue;
    Float_t pt  = upgrade_->egEt.at(ue);
    if(pt >= eleptmax) eleptmax = pt;
  }

  for (UInt_t ue=0; ue < upgrade_->nJets; ue++) {        
    Int_t bx = upgrade_->jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(upgrade_->jetEta.at(ue)) > 3. ? true : false;
    if(!isFwdJet) continue;
    Float_t pt = upgrade_->jetEt.at(ue);
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

  for (UInt_t ue=0; ue < upgrade_->nEGs; ue++) {
    Int_t bx = upgrade_->egBx.at(ue);        		
    if(bx != 0) continue;
    Float_t pt  = upgrade_->egEt.at(ue);
    if(pt >= eleptmax) eleptmax = pt; 
  }

  if(upgrade_->nJets < 2) return;

  for (UInt_t ue=0; ue < upgrade_->nJets; ue++) {
    Int_t bx = upgrade_->jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(upgrade_->jetEta.at(ue)) > 3. ? true : false;
    if(isFwdJet) continue;
    Float_t pt = upgrade_->jetEt.at(ue);

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

  for (UInt_t ue=0; ue < upgrade_->nEGs; ue++){
    Int_t bx = upgrade_->egBx.at(ue);        		
    if(bx != 0) continue;
    Float_t pt  = upgrade_->egEt.at(ue);
    Float_t eta = upgrade_->egEta.at(ue);
    if(fabs(eta) > 2.1) continue;  // eta = 5 - 16
    if(pt >= eleptmax){
      eleptmax = pt; 
      elemaxeta = upgrade_->egEta.at(ue);
    }
  }  // end loop over EM objects

  if(upgrade_->nJets < 3) return;

  for (UInt_t ue=0; ue < upgrade_->nJets; ue++) {
    Int_t bx = upgrade_->jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(upgrade_->jetEta.at(ue)) > 3. ? true : false;
    if(isFwdJet) continue;
    Float_t pt = upgrade_->jetEt.at(ue);
    Float_t jeteta = upgrade_->jetEta.at(ue);

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

void L1AlgoFactory::IsoEGer_TauJetEta2p17Pt(Float_t& egcut, Float_t& taucut) {

  Float_t eleptmax  = -10.;
  Float_t eleetamax = -999.;
  Float_t maxpttau  = -10.;

  if(upgrade_->nEGs < 1) return;
  for (UInt_t ue=0; ue < upgrade_->nEGs; ue++) {
    Int_t bx = upgrade_->egBx.at(ue);        		
    if(bx != 0) continue;
    if(!upgrade_->egIso.at(ue)) continue;
    Float_t pt  = upgrade_->egEt.at(ue);
    Float_t eta = upgrade_->egEta.at(ue);
    if(fabs(eta) > eleER) continue;  // eta = 5 - 16
    if(pt >= eleptmax){
      eleptmax = pt;
      eleetamax = eta;
    }
  }

  for(UInt_t ue=0; ue < upgrade_->nTaus; ue++) {
    Int_t bx = upgrade_->tauBx.at(ue);        		
    if(bx != 0) continue; 
    Float_t pt = upgrade_->tauEt.at(ue);
    Float_t eta = upgrade_->tauEta.at(ue);
    if(fabs(eta) > tauER) continue;  // eta = 5 - 16

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

  if(upgrade_->nJets < 4) return;
  if(upgrade_->nTaus < 1) return;

  for (UInt_t ue=0; ue < upgrade_->nJets; ue++) {
    Int_t bx = upgrade_->jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(upgrade_->jetEta.at(ue)) > jetCentFwd ? true : false;
    if(isFwdJet) continue;
    Float_t pt = upgrade_->jetEt.at(ue);

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

  for(UInt_t ue=0; ue < upgrade_->nTaus; ue++) {
    Int_t bx = upgrade_->tauBx.at(ue);        		
    if(bx != 0) continue; 
    Float_t pt = upgrade_->tauEt.at(ue);
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
  if(upgrade_->sumBx[EtSumType::ETM]==0) TheETM =upgrade_->sumEt[EtSumType::ETM];
  ETMcut = TheETM;
  return;
}

void L1AlgoFactory::HTTVal(Float_t& HTTcut) {

  Float_t TheHTT = -10;
  if(upgrade_->sumBx[EtSumType::HTT]==0) TheHTT =upgrade_->sumEt[EtSumType::HTT];
  HTTcut = TheHTT;
  return;
}

void L1AlgoFactory::HTMVal(Float_t& HTMcut) {

  Float_t TheHTM = -10;
  if (upgrade_->sumBx[EtSumType::HTM]==0) TheHTM = upgrade_->sumEt[EtSumType::HTM];
  HTMcut = TheHTM;
  return;
}

void L1AlgoFactory::ETTVal(Float_t& ETTcut) {

  Float_t TheETT = -10;
  if(upgrade_->sumBx[EtSumType::ETT]==0) TheETT = upgrade_->sumEt[EtSumType::ETT];
  ETTcut = TheETT;
  return;
}

//**************************************************************************//
//                                   Cross                                  //
//**************************************************************************//
// ===  FUNCTION  ============================================================
//         Name:  L1AlgoFactory::PassMuonQual
//  Description:  
// ===========================================================================
bool L1AlgoFactory::PassMuonQual(int imu, bool isMuHighQual) const
{
  (void)isMuHighQual; // Not used for now as the L1UpgradeTree only store 0 or 1
  return upgrade_->muonQual.at(imu) == 1;
}       // -----  end of function L1AlgoFactory::PassMuonQual  -----

Bool_t L1AlgoFactory::Mu_HTT(Float_t mucut, Float_t HTcut) {
  Float_t tmp_mucut = -10.;
  Float_t tmp_HTcut = -10.;
  Mu_HTTPt(tmp_mucut,tmp_HTcut);
  if(tmp_mucut >= mucut && tmp_HTcut >= HTcut) return true;
  return false;
}

void L1AlgoFactory::Mu_HTTPt(Float_t& mucut, Float_t& HTcut ) {

  Float_t muptmax = -10.;

  Int_t Nmu = upgrade_->nMuons;
  for (Int_t imu=0; imu < Nmu; imu++) {
    Int_t bx = upgrade_->muonBx.at(imu);
    if(bx != 0) continue;
    Float_t pt = upgrade_->muonEt.at(imu);                       
    if(!PassMuonQual(imu)) continue;
    if(pt >= muptmax) muptmax = pt;
  }

  Float_t TheHTT = -10;
  if(upgrade_->sumBx[EtSumType::HTT]==0) 
    TheHTT =upgrade_->sumEt[EtSumType::HTT];

  if(muptmax >= 0.){
    mucut = muptmax;
    HTcut = TheHTT;
  }

  return;
}

Bool_t L1AlgoFactory::Muer_ETM(Float_t mucut, Float_t ETMcut) {
  Float_t tmp_mucut = -10.;
  Float_t tmp_ETMcut = -10.;
  Muer_ETMPt(tmp_mucut,tmp_ETMcut);
  if(tmp_mucut >= mucut && tmp_ETMcut >= ETMcut) return true;
  return false;
}

void L1AlgoFactory::Muer_ETMPt(Float_t& mucut, Float_t& ETMcut ) {

  Float_t muptmax = -10.;

  Int_t Nmu = upgrade_->nMuons;
  for (Int_t imu=0; imu < Nmu; imu++) {
    Int_t bx = upgrade_->muonBx.at(imu);
    if(bx != 0) continue;
    Float_t pt = upgrade_->muonEt.at(imu);                       
    if(!PassMuonQual(imu)) continue;
    Float_t eta = upgrade_->muonEta.at(imu); 
    if(fabs(eta) > muonER) continue;
    if(pt >= muptmax) muptmax = pt;
  }

  Float_t TheETM = -10;
  if(upgrade_->sumBx[EtSumType::ETM]==0) 
    TheETM =upgrade_->sumEt[EtSumType::ETM];

  if(muptmax >= 0.){
    mucut = muptmax;
    ETMcut = TheETM;
  }

  return;
}

Bool_t L1AlgoFactory::SingleEG_Eta2p1_HTT(Float_t egcut, Float_t HTTcut, Bool_t isIsolated) {
  Float_t tmp_egcut = -10.;
  Float_t tmp_HTTcut = -10.;
  SingleEG_Eta2p1_HTTPt(tmp_egcut,tmp_HTTcut,isIsolated);
  if(tmp_egcut >= egcut && tmp_HTTcut >= HTTcut) return true;
  return false;
}

void L1AlgoFactory::SingleEG_Eta2p1_HTTPt(Float_t& egcut, Float_t& HTTcut, Bool_t isIsolated) {

  Float_t eleptmax = -10.;

  Int_t Nele = upgrade_->nEGs;
  for (Int_t ue=0; ue < Nele; ue++) {
    Int_t bx = upgrade_->egBx.at(ue);  
    if(bx != 0) continue;
    if(isIsolated && !upgrade_->egIso.at(ue)) continue;
    Float_t eta = upgrade_->egEta.at(ue);
    if(fabs(eta) > eleER ) continue;  
    Float_t pt = upgrade_->egEt.at(ue);   
    if(pt >= eleptmax) eleptmax = pt;
  }

  Float_t TheHTT = -10;
  if(upgrade_->sumBx[EtSumType::HTT]==0) 
    TheHTT =upgrade_->sumEt[EtSumType::HTT];

  if(eleptmax >= 0.){
    egcut = eleptmax;
    HTTcut = TheHTT;
  }

  return;
}

Bool_t L1AlgoFactory::Muer_TauJetEta2p17(Float_t mucut, Float_t taucut, Bool_t isIsolated){
  Float_t tmp_mucut  = -10.;
  Float_t tmp_taucut = -10.;
  Muer_TauJetEta2p17Pt(tmp_mucut, tmp_taucut,isIsolated);
  if(tmp_mucut >= mucut && tmp_taucut >= taucut) return true;
  return false;
}

void L1AlgoFactory::Muer_TauJetEta2p17Pt(Float_t& mucut, Float_t& taucut, Bool_t isIsolated) {

  Float_t maxptmu  = -10.;
  Float_t maxpttau = -10.;

  Int_t Nmu = upgrade_->nMuons;
  if(Nmu < 1) return;
  for (Int_t imu=0; imu < Nmu; imu++) {
    Int_t bx = upgrade_->muonBx.at(imu);
    if(bx != 0) continue;
    Float_t pt = upgrade_->muonEt.at(imu);                       
    if(!PassMuonQual(imu)) continue;
    Float_t eta = upgrade_->muonEta.at(imu);        
    if(fabs(eta) > muonER ) continue;
    if(pt >= maxptmu) maxptmu = pt;
  }

  for(Int_t ue=0; ue < upgrade_->nTaus; ue++) {
    Int_t bx = upgrade_->tauBx.at(ue);        		
    if(bx != 0) continue; 
    if(!isIsolated && !upgrade_->tauIso.at(ue)) continue;

    Float_t pt = upgrade_->tauEt.at(ue);
    Float_t eta = upgrade_->tauEta.at(ue);
    if(fabs(eta) > tauER) continue;  // eta = 5 - 16

    if(pt >= maxpttau) maxpttau = pt;
  }

  if(maxptmu >= 0.){
    mucut  = maxptmu;
    taucut = maxpttau;
  }

  return;

}

Bool_t L1AlgoFactory::DoubleJetCentral_ETM(Float_t jetcut1, Float_t jetcut2, Float_t ETMcut) {
  Float_t tmp_jetcut1 = -10.;
  Float_t tmp_jetcut2 = -10.;
  Float_t tmp_ETMcut  = -10.;
  DoubleJetCentral_ETMPt(tmp_jetcut1,tmp_jetcut2,tmp_ETMcut);
  if(tmp_jetcut1 >= jetcut1 && tmp_jetcut2 >= jetcut2 && tmp_ETMcut >= ETMcut) return true;
  return false;
}

void L1AlgoFactory::DoubleJetCentral_ETMPt(Float_t& jetcut1, Float_t& jetcut2, Float_t& ETMcut){

  Float_t jetptmax1 = -10.;
  Float_t jetptmax2 = -10.;

  Float_t TheETM = -10;
  if(upgrade_->sumBx[EtSumType::ETM]==0) 
    TheETM =upgrade_->sumEt[EtSumType::ETM];

  if(upgrade_->nJets < 2) return;

  for (Int_t ue=0; ue < upgrade_->nJets; ue++) {
    Int_t bx = upgrade_->jetBx.at(ue);        		
    if(bx != 0) continue;
    Bool_t isFwdJet = fabs(upgrade_->jetEta.at(ue)) > jetCentFwd ? true : false;
    if(isFwdJet) continue;
    Float_t pt = upgrade_->jetEt.at(ue);

    if(pt >= jetptmax1){
      jetptmax2 = jetptmax1;
      jetptmax1 = pt;
    }
    else if(pt >= jetptmax2) jetptmax2 = pt;
  }       

  if(jetptmax2 >= 0.){
    jetcut1 = jetptmax1;
    jetcut2 = jetptmax2;
    ETMcut = TheETM;
  }

  return;
}

Bool_t L1AlgoFactory::DoubleEG_HT(Float_t egcut, Float_t HTcut) {
  Float_t tmp_egcut = -10.;
  Float_t tmp_HTcut = -10.;
  DoubleEG_HTPt(tmp_egcut,tmp_HTcut);
  if(tmp_egcut >= egcut && tmp_HTcut >= HTcut) return true;
  return false;
}

void L1AlgoFactory::DoubleEG_HTPt(Float_t& EGcut, Float_t& HTcut) {

  Float_t eleptmax1  = -10.;
  Float_t eleptmax2  = -10.;
  Float_t ele1Phimax = -1000.;
  Float_t ele1Etamax = -1000.;

  if(upgrade_->nEGs < 2) return;

  for(Int_t ue=0; ue < upgrade_->nEGs ; ue++) {
    Int_t bx = upgrade_->egBx.at(ue);  
    if(bx != 0) continue;
    Float_t pt = upgrade_->egEt.at(ue);
    Float_t phi = upgrade_->egPhi.at(ue);
    Float_t eta = upgrade_->egEta.at(ue);

    if(fabs(pt-eleptmax1) < 0.001 && fabs(phi-ele1Phimax) < 0.001 && fabs(eta-ele1Etamax) < 0.001) continue; //to avoid double counting in noniso/relaxiso lists

    if(pt >= eleptmax1){
      eleptmax2 = eleptmax1;
      eleptmax1 = pt;
      ele1Phimax = phi;
      ele1Etamax = eta;
    }
    else if(pt >= eleptmax2) eleptmax2 = pt;
  }

  Float_t TheHTT = -10;
  if(upgrade_->sumBx[EtSumType::HTT]==0) 
    TheHTT =upgrade_->sumEt[EtSumType::HTT];

  if(eleptmax2 >= 0.){
    EGcut = eleptmax2;
    HTcut = TheHTT;
  }

  return;
}
Bool_t L1AlgoFactory::Jet_MuOpen_Mu_dPhiMuMu1(Float_t jetcut, Float_t mucut) {
  Float_t tmp_jetcut = -10.;
  Float_t tmp_mucut = -10.;
  Jet_MuOpen_Mu_dPhiMuMu1Pt(tmp_jetcut,tmp_mucut);
  if(tmp_jetcut >= jetcut && tmp_mucut >= mucut) return true;
  return false;
}

void L1AlgoFactory::Jet_MuOpen_Mu_dPhiMuMu1Pt(Float_t& jetcut, Float_t& mucut) {

  //Find the highest pt jet with deltaphi condition
  Float_t jetptmax = -10.;
  Int_t Nj = upgrade_->nJets;
  Int_t Nmu = upgrade_->nMuons;

  for(Int_t ue=0; ue < Nj; ue++){
    Int_t bxjet = upgrade_->jetBx.at(ue);        		
    if(bxjet != 0) continue;
    Float_t pt = upgrade_->jetEt.at(ue);
    if(pt < jetptmax) continue;
    Float_t phijet = upgrade_->jetPhi.at(ue);

    Bool_t corr = false;

    for(Int_t imu=0; imu < Nmu; imu++){
      Int_t bx = upgrade_->muonBx.at(imu);
      if (bx != 0) continue;
      if(!PassMuonQual(imu)) continue;
      Float_t muonpt = upgrade_->muonEt.at(imu);			
      if(muonpt < 0.) continue;
      Float_t muphi = upgrade_->muonPhi.at(imu);			
      if(fabs(muphi-phijet) < MuOpenJetCordPhi) corr = true;
    }

    if(corr) jetptmax = pt;
  }

  //Loop over the muons list twice and save all pairs that satisfy the deltaphi veto
  Bool_t corr = false;
  Float_t maxpt1 = -10.;
  Float_t maxpt2 = -10.;

  std::vector<std::pair<Float_t,Float_t> > muonPairs;
  for (Int_t imu=0; imu < Nmu; imu++) {
    Int_t bx = upgrade_->muonBx.at(imu);
    if (bx != 0) continue;
    Float_t pt = upgrade_->muonEt.at(imu);			
    if(!PassMuonQual(imu)) continue;
    Float_t phi1 = upgrade_->muonPhi.at(imu);			

    for (Int_t imu2=0; imu2 < Nmu; imu2++) {
      if (imu2 == imu) continue;
      Int_t bx2 = upgrade_->muonBx.at(imu2);
      if (bx2 != 0) continue;
      Float_t pt2 = upgrade_->muonEt.at(imu2);			
      if(!PassMuonQual(imu)) continue;
      Float_t phi2 = upgrade_->muonPhi.at(imu2);			

      Float_t dphi = phi1 - phi2; //Should get the binning, but for GMT is quite fine
      if(fabs(dphi) > MuMudPhi){
        corr = true;
        muonPairs.push_back(std::pair<Float_t,Float_t>(pt,pt2));
      }
    }
  }

  //Select the muon pair in which one of the two muons is the highest pt muon satisfying the deltaphi veto
  if(corr){
    std::vector<std::pair<Float_t,Float_t> >::const_iterator muonPairIt  = muonPairs.begin();
    std::vector<std::pair<Float_t,Float_t> >::const_iterator muonPairEnd = muonPairs.end();
    for(; muonPairIt != muonPairEnd; ++muonPairIt){
      Float_t pt1 = muonPairIt->first;
      Float_t pt2 = muonPairIt->second;

      if(pt1 > maxpt1 || (fabs(maxpt1-pt1)<10E-2 && pt2>maxpt2) ) {
        maxpt1 = pt1;
        maxpt2 = pt2;
      }
    }
  }

  Float_t maxptmu = maxpt1 > maxpt2 ? maxpt1 : maxpt2; //only the highest pt muon counts for the correlation, the second muon is Open

  if(jetptmax > 0. && maxptmu > 0.){
    jetcut = jetptmax;
    mucut = maxptmu;
  }

  return;
}

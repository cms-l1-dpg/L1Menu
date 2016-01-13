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
  
  void run(bool runOnData, std::string resultTag, float crossSec, int nBunches, int isCrossSec, int nEvents = 0);

private :
  
  float ScaleFactor(float nZeroBias, float nBunches);
  //void  SingleJetPt(Float_t& ptcut, Bool_t isCentral = false);
  float SingleTauPt();
  float SingleMuEta(float eta);
  //void  SingleEGPt(Float_t& ptcut, Bool_t isIsolated, Bool_t isER);
  float SingleEGEta(float ptCut, bool doIso);
  float SingleJetEta(float pt, Int_t accept_flag = 0);

  void setRateError(TH1F* histo);
  
  //void ETMVal(Float_t& ETMcut);
  //void HTTVal(Float_t& HTTcut);
  //void HTMVal(Float_t& HTMcut);
  //void ETTVal(Float_t& ETTcut);
  
  std::map<std::string,TH1F*> hTH1F;
  std::map<std::string,TH2F*> hTH2F;
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
    Float_t pt = upgrade_ -> muonEt.at(imu);                       
      if ( pt > maxPt) 
	{
	  maxPt = pt;
	  iMuMaxPt = imu;
	}
    }
  
  float eta = maxPt>ptCut ? upgrade_ -> muonEta.at(iMuMaxPt) : -10.; 

  //cout << "max mu pt = " << maxPt << "  max eta = " << eta << endl;

  return eta;
}

float BasicRatePlots::SingleTauPt() {
  
  float maxPt = -10;

  for (UInt_t ue=0; ue < upgrade_ -> nTaus; ue++) {
    Int_t bx = upgrade_ -> tauBx.at(ue);        		
    if (bx != 0) continue; 
    Float_t pt = upgrade_ -> tauEt.at(ue);
    if (pt >= maxPt) maxPt = pt;
  } 
  
  return maxPt;
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

//void BasicRatePlots::SingleEGPt(Float_t& cut, Bool_t isIsolated , Bool_t isER) {

  
  ////if(nEGs < 1) return;

  //Float_t ptmax = -10.;

  //for(UInt_t ue=0; ue < upgrade_ -> nEGs; ue++) {
    //Int_t bx = upgrade_ -> egBx.at(ue);  
    //if(bx != 0) continue;
    //if(isIsolated && !(upgrade_ -> egIso.at(ue))) continue;
    //Float_t eta = upgrade_ -> egEta.at(ue);
    //if(fabs(eta) > 2.1 && isER) continue;  // eta = 5 - 16

    //Float_t pt = upgrade_ -> egEt.at(ue);    // the rank of the electron
    //if(pt >= ptmax) ptmax = pt;
  //}

  //cut = ptmax;

  //return;
//}

//void BasicRatePlots::SingleJetPt(Float_t& cut, Bool_t isCentral) {

  //Float_t ptmax = -10.;
  //Int_t Nj = upgrade_ -> nJets ;
  //for(Int_t ue=0; ue < Nj; ue++) {
    //Int_t bx = upgrade_ -> jetBx.at(ue);
    //if(bx != 0) continue;
    //Bool_t isFwdJet = fabs(upgrade_ -> jetEta.at(ue)) > 3. ? true : false;
    //if(isCentral && isFwdJet) continue;
    ////if(NOTauInJets && upgrade_->Taujet[ue]) continue;
    ////if(isCentral && noHF && (upgrade_->jetEta.at(ue) < 5 || upgrade_->jetEta.at(ue) > 17)) continue;

    //Float_t pt = upgrade_ -> jetEt.at(ue);
    //if(pt >= ptmax) ptmax = pt;
  //}

  //cut = ptmax;
  //return;
//}

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

//void BasicRatePlots::ETMVal(Float_t& ETMcut ) {

  //Float_t TheETM = -10;
  //if(upgrade_ ->sumBx[2]==0) TheETM =upgrade_ ->sumEt[2];
  //ETMcut = TheETM;
  //return;
//}

//void BasicRatePlots::HTTVal(Float_t& HTTcut) {

  //Float_t TheHTT = -10;
  //if(upgrade_ ->sumBx[1]==0) TheHTT =upgrade_ ->sumEt[1];
  //HTTcut = TheHTT;
  //return;
//}

//void BasicRatePlots::HTMVal(Float_t& HTMcut) {

  //Float_t TheHTM = -10;
  //if (upgrade_ ->sumBx[3]==0) TheHTM = upgrade_ ->sumEt[3];
  //HTMcut = TheHTM;
  //return;
//}

//void BasicRatePlots::ETTVal(Float_t& ETTcut) {

  //Float_t TheETT = -10;
  //if(upgrade_ ->sumBx[0]==0) TheETT = upgrade_ ->sumEt[0];
  //ETTcut = TheETT;
  //return;
//}

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


void BasicRatePlots::run(bool runOnData, std::string resultTag, float crossSec, int nBunches, int isCrossSec, int nEvents) {

  system("mkdir -p results");
  std::string resultName = "results/results_" + resultTag + (isCrossSec ? "_XSEC" : "_RATE") + ".root";
  TFile *outFile = new TFile((resultName).c_str(),"recreate");
  outFile->cd();

  //Event Counter
  hTH1F["nEvts"]       = new TH1F("nEvts","Number of Events Processed",1,-0.5,.5);
  //Single stuff
  hTH1F["nJetVsPt"]    = new TH1F("nJetVsPt","SingleJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hTH1F["nTauVsPt"]    = new TH1F("nTauVsPt","SingleTau; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hTH1F["nJetCenVsPt"] = new TH1F("nJetCenVsPt","SingleJetCentral; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hTH1F["nEGVsPt"]     = new TH1F("nEGVsPt","SingleEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH1F["nEGErVsPt"]   = new TH1F("nEGErVsPt","SingleEGer; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH1F["nIsoEGVsPt"]  = new TH1F("nIsoEGVsPt","SingleIsoEGer; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH1F["nMuVsPt"]     = new TH1F("nMuVsPt","SingleMu; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hTH1F["nMuErVsPt"]   = new TH1F("nMuErVsPt","SingleMu |#eta|<2.1; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hTH1F["nMuVsEta"]    = new TH1F("nMuVsEta","nMuVsEta",24,-2.4,2.4);
  hTH1F["nEGVsEta"]    = new TH1F("nEGVsEta","nEGVsEta",50,-3.,3.);
  hTH1F["nIsoEGVsEta"] = new TH1F("nIsoEGVsEta","nIsoEGVsEta",50,-3.,3.);
  hTH1F["nJetVsEta"]   = new TH1F("nJetVsEta","nJetVsEta",50,-5.,5.);
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
  hTH1F["nDiEGVsPt"]         = new TH1F("nDiEGVsPt","DiEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH1F["nDiIsoEGVsPt"]      = new TH1F("nDiIsoEGVsPt","DiIsoEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH2F["nEGPtVsPt"]         = new TH2F("nEGPtVsPt","DoubleEle; p_{T} cut EG_{1}; p_{T} cut EG_{2}",65,-0.5,64.5,65,-0.5,64.5);
  hTH2F["nIsoEGPtVsPt"]      = new TH2F("nIsoEGPtVsPt","DoubleIsolEle; p_{T} cut EG_{1}; p_{T} cut EG_{2}",65,-0.5,64.5,65,-0.5,64.5);
  hTH2F["nMuPtVsPt"]         = new TH2F("nMuPtVsPt","DoubleMu; p_{T} cut mu_{1}; p_{T} cut mu_{2}",41,-0.25,20.25,41,-0.25,20.25);
  hTH2F["nOniaMuPtVsPt"]     = new TH2F("nOniaMuPtVsPt","DoubleMu_Er_HighQ_WdEta22 (Quarkonia); p_{T} cut mu_{1}; p_{T} cut mu_{2}",41,-0.25,20.25,41,-0.25,20.25);

  //Sums
  hTH1F["nHTTVsHTT"] = new TH1F("nHTTVsHTT","HTT; HTT cut; rate [Hz]",512,-.5,511.5);
  hTH1F["nETTVsETT"] = new TH1F("nETTVsETT","ETT; ETT cut; rate [Hz]",512,-.5,511.5);
  hTH1F["nETMVsETM"] = new TH1F("nETMVsETM","ETM; ETM cut; rate [Hz]",512,-.5,511.5);

  if (isCrossSec) {
    std::map<std::string,TH1F*>::iterator hTH1FIt  = hTH1F.begin();
    std::map<std::string,TH1F*>::iterator hTH1FEnd = hTH1F.end();

    for(; hTH1FIt!=hTH1FEnd; ++hTH1FIt)
      {
	hTH1FIt->second->GetYaxis()->SetTitle("cross section [#mubarn]");
      }
  }

  Double_t nZeroBias = 0.;

  int nLumi(0),currentLumi(-1);

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

    if (event_!=NULL && event_ -> lumi != currentLumi){
      std::cout << "New Lumi section: " << event_->lumi << std::endl;      
      currentLumi=event_ -> lumi;
      nLumi++;
    }
    
    hTH1F["nEvts"]->Fill(0.);  // count number of events processed

    nZeroBias += 1.;

    float jetPt     = 0.; SingleJetPt(jetPt);
    float jetCenPt  = 0.; SingleJetPt(jetCenPt,true);
    float jetEta    = SingleJetEta(36.);
    float jetEta_Central = SingleJetEta(36.,1);
    float jetEta_Fwd    = SingleJetEta(36.,2);

    float tauPt     = SingleTauPt();

    float htt       = 0.; HTTVal(htt);
    float ett       = 0.; ETTVal(ett);
    float etm       = 0.; ETMVal(etm);

    float egPt      = 0.; SingleEGPt(egPt,false,false);
    float egErPt    = 0.; SingleEGPt(egErPt,false,true);
    float isoEgPt   = 0.; SingleEGPt(isoEgPt,true,false);
    float egEta     = SingleEGEta(16.,false);
    float isoegEta  = SingleEGEta(16.,true);

    //cout << "Event number = " << nZeroBias << endl;
    //float muPt     = -10.; SingleMuPt(muPt,false);
    //cout << "muPt = " << muPt << endl; 
    //float muErPt   = -10.; SingleMuPt(muErPt,true);

    //cout << "Event number = " << nZeroBias << endl;
    float muEta    = SingleMuEta(16.);
    //cout << "muEta = " << muEta << endl;

    // float doubleMuPt1 = -10.; 
    // float doubleMuPt2 = -10.;
    // DoubleMuPt(doubleMuPt1,doubleMuPt2);
    // 
    // float oniaMuPt1 = 0.;
    // float oniaMuPt2 = 0.;
    // OniaPt(oniaMuPt1,oniaMuPt2,22);
    // 
    // 
    // float dijetPt1    = -10.;
    // float dijetPt2    = -10.;
    // float diCenjetPt1 = -10.;
    // float diCenjetPt2 = -10.;
    // DoubleJetPt(dijetPt1,dijetPt2);
    // DoubleJetPt(diCenjetPt1,diCenjetPt2,true);
    // Float_t dummy = -1;
    // float ditauPt    = -10.; DoubleTauJetEta2p17Pt(dummy,ditauPt);
    // dummy = -1.;
    // float quadjetPt  = -10.; QuadJetPt(dummy,dummy,dummy,quadjetPt);
    // dummy = -1.;
    // float quadjetCPt = -10.; QuadJetPt(dummy,dummy,dummy,quadjetCPt,true);
    // dummy = -1.;
    // 
    // float diEG1     = -10.;
    // float diEG2     = -10.;
    // float diIsolEG1 = -10.;
    // float diIsolEG2 = -10.;
    // DoubleEGPt(diEG1,diEG2,false);
    // DoubleEGPt(diIsolEG1,diIsolEG2,true);

    if(muEta > -9.) hTH1F["nMuVsEta"]->Fill(muEta);
    hTH1F["nEGVsEta"]->Fill(egEta);
    hTH1F["nIsoEGVsEta"]->Fill(isoegEta);
    hTH1F["nJetVsEta"]->Fill(jetEta);
    hTH1F["nJetVsEta_Central"]->Fill(jetEta_Central);
    hTH1F["nJetVsEta_Fwd"]->Fill(jetEta_Fwd);

    for(int ptCut=0; ptCut<256; ++ptCut) {
      if(jetPt>=ptCut)	  hTH1F["nJetVsPt"]->Fill(ptCut);
      if(jetCenPt>=ptCut) hTH1F["nJetCenVsPt"]->Fill(ptCut);
      if(tauPt>=ptCut)	  hTH1F["nTauVsPt"]->Fill(ptCut);
    // 
    //   if(dijetPt2>=ptCut){
    // 	hTH1F["nDiJetVsPt"]->Fill(ptCut);
    // 
    // 	for(int ptCut_0=ptCut; ptCut_0<256; ++ptCut_0) {
    // 	  if(dijetPt1>=ptCut_0) hTH2F["nAsymDiJetVsPt"]->Fill(ptCut_0,ptCut);
    // 	}
    //   }
    // 
    //   if(diCenjetPt2>=ptCut){
    // 	hTH1F["nDiCenJetVsPt"]->Fill(ptCut);
    // 
    // 	for(int ptCut_0=ptCut; ptCut_0<256; ++ptCut_0) {
    // 	  if(diCenjetPt1>=ptCut_0) hTH2F["nAsymDiCenJetVsPt"]->Fill(ptCut_0,ptCut);
    // 	}
    //   }
    // 
    //   if(ditauPt>=ptCut)    hTH1F["nDiTauVsPt"]->Fill(ptCut);
    //   if(quadjetPt>=ptCut)  hTH1F["nQuadJetVsPt"]->Fill(ptCut);
    //   if(quadjetCPt>=ptCut) hTH1F["nQuadCenJetVsPt"]->Fill(ptCut);
    // 
    }//loop on 256
    //   
    for(int ptCut=0; ptCut<65; ++ptCut) {
       if(egPt>=ptCut)    hTH1F["nEGVsPt"]->Fill(ptCut);
       if(egErPt>=ptCut)  hTH1F["nEGErVsPt"]->Fill(ptCut);
       if(isoEgPt>=ptCut) hTH1F["nIsoEGVsPt"]->Fill(ptCut);
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
       if(htt>httCut) hTH1F["nHTTVsHTT"]->Fill(httCut);
       if(ett>httCut) hTH1F["nETTVsETT"]->Fill(httCut);
       if(etm>httCut) hTH1F["nETMVsETM"]->Fill(httCut);
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
  std::cout << "Scale factor applied to histograms = " << scaleFactor << std::endl;

  std::map<std::string,TH1F*>::iterator hTH1FIt  = hTH1F.begin();
  std::map<std::string,TH1F*>::iterator hTH1FEnd = hTH1F.end();

  for(; hTH1FIt!=hTH1FEnd; ++hTH1FIt) {
    TH1F* histo = hTH1FIt->second;
    setRateError(histo);
    histo->Scale(scaleFactor);
  }

  std::map<std::string,TH2F*>::iterator hTH2FIt  = hTH2F.begin();
  std::map<std::string,TH2F*>::iterator hTH2FEnd = hTH2F.end();

  for(; hTH2FIt!=hTH2FEnd; ++hTH2FIt) {
    TH2F* histo = hTH2FIt->second;
    histo->Scale(scaleFactor);
  }

  outFile->Write();
  outFile->Close();
  delete outFile;
}

// --------------------------------------------------------------------

void goRatePlots(std::string fileType, int isCrossSec = false, int nEvents = 0) 
{
  int nBunches50ns = 1368;
  int nBunches25ns = 2508; //2508 is what agreed with TSG for # bunches
  int nBunches25ns_run256843 = 1021;
  int nBunches = -1;

  float xSec13TeV = isCrossSec ? 78.26 : 80.; // Using McM for cross section comparison and 80 (agreed with TSG) for rates
  float xSec8TeV  = 72.7; 

  std::string filename;
  bool isData(true);
  
  if (fileType == "13TEV_40PU_2016_RE-EMUL")
    {
      isData = false;
      nBunches = nBunches25ns;
      filename = "/afs/cern.ch/user/p/pellicci/data2/L1DPG/root/2016/v2/40PU_25ns_Stage2/40PU_25ns_Stage2_1.root";
    }
  else if (fileType == "13TEV_20PU_2016_RE-EMUL")
    {
      isData = false;
      nBunches = nBunches25ns;
      filename = "/afs/cern.ch/user/p/pellicci/data2/L1DPG/root/2016/v2/20PU_25ns_Stage2/20PU_25ns_Stage2_1.root";
    }
  else if (fileType == "RUN256843_Stage2")
    {
      isData = true;      
      // filename = "/data/user/gennai/L1Ntuple/l1t_debug-stage-2_256843.root";
      // filename = "root://cmseos.fnal.gov//store/user/lpctrig/apana/Stage2/ZeroBias1/crab_ZeroBias1_Run2015D-v1/151230_012024/0000/l1t_stage2_2.root";
      filename = "ntuples_256843_stage2_full.list";
    }
  else if (fileType == "RUN256843_Stage1")
    {
      isData = true;      
      // filename = "/data/user/gennai/L1Ntuple/l1t_debug-stage-2_256843.root";
      filename = "ntuples_256843_stage1.list";
    }
  else if (fileType == "Stage2_Simone")
    {
      isData = true;      
      // filename = "/data/user/gennai/L1Ntuple/l1t_debug-stage-2_256843.root";
      filename = "ntuples_256843_stage2_Simone.list";
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
  basicRatePlots.run(isData,fileType,xSec13TeV,nBunches,isCrossSec,nEvents);
    
}


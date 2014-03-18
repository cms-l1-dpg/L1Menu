#include "../../L1Ntuples/macros/L1Ntuple.h"
#include <algorithm>
#include<map>

Int_t etaMuIdx(Double_t eta) {

  static const size_t ETAMUBINS = 65;
  static const Double_t ETAMU[] = { -2.45,-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45 };
  size_t etaIdx = 0.;
  for (size_t idx=0; idx<ETAMUBINS; idx++) {
    if (eta>=ETAMU[idx] and eta<ETAMU[idx+1])
      etaIdx = idx;
  }
  return int(etaIdx);
}

class BasicRatePlots : public L1Ntuple
{
public :
  
  //constructor    
  BasicRatePlots(std::string filename) : L1Ntuple(filename) {}
  BasicRatePlots()  {}
  ~BasicRatePlots() {}
  
  void run(bool runOnData, std::string resultTag, int minLs, int maxLs, 
	   float crossSec, float avPU, int nBunches, int isCrossSec, int nEvents = 0);
  
private :
  
  void FillBits();

  float ScaleFactor(float nZeroBias, float nBunches);

  //Single Stuff
  float SingleJetCentralPt();
  float SingleJetCentralEta(Float_t cut );
  float SingleJetPt();
  float SingleTauPt();
  float SingleEGPt();
  float SingleIsoEGPt();
  float SingleMuPt();
  float SingleMuErPt();
  float SingleMuEta(float eta);

  //MultiStuff
  float QuadJetCentral();
  void DoubleMu(Float_t & muPt1, Float_t & muPt2);
  void Onia(Float_t & muPt1, Float_t & muPt2, Float_t etacut, Int_t delta);
  void EGIsoEGPt(Float_t &myIsoEGPt, Float_t &myEGPt);

  //Sums
  float HTT();
  float ETT();
  float ETM();
  
  //Cross
  void Mu_HTT(Float_t & muPt, Float_t & myHTT);

  //GMT stuff
  float DttfPt();
  float RpcbPt();
  float RpcfPt();
  float CsctfPt();

  void setRateError(TH1F* histo);

  vector<float> SortJets(bool doCentral, bool doTau);
  vector<float> SortEGs(bool doIsol);


  float computeAvgLumi(float xSec, float avPU, int nBunches) { return 11246. * avPU * nBunches / (1E7 * xSec); };

  bool PhysicsBits[128];

  std::map<std::string,TH1F*> hTH1F;
  std::map<std::string,TH2F*> hTH2F;

  vector<float> sortedCenJets,sortedJets,sortedTaus,sortedEGs,sortedIsolEGs;

};

void BasicRatePlots::FillBits() {

  //Fill the physics bits:
  for (int ibit=0; ibit < 128; ibit++) 
    {
      PhysicsBits[ibit] = 0;
      if (ibit<64) PhysicsBits[ibit] = (gt_->tw1[2]>>ibit)&1;
      else PhysicsBits[ibit] = (gt_->tw2[2]>>(ibit-64))&1;
    }

}       

// ------------------------------------------------------------------

// scale factor computed w.r.t. ZeroBias rate fratcion and # bunches 
float BasicRatePlots::ScaleFactor(float nZeroBias, float nBunches) {

  float scal = 11246.; // ZB per bunch in Hz
  scal /= nZeroBias;
  scal *= nBunches;

  return scal;
}

float BasicRatePlots::SingleMuPt(){

  float maxPt = -10;
  
  int Nmu = gmt_ -> N;
  for (int imu=0; imu < Nmu; imu++) 
    {
      int bx = gmt_ -> CandBx[imu];       // BX = 0, +/- 1 or +/- 2
      if(bx != 0) continue;
      int qual = gmt_ -> Qual[imu];
      if(qual < 4) continue;
      float pt = gmt_ -> Pt[imu];         // get pt to get highest pt one
      if( pt > maxPt) maxPt = pt;
    }
  
  return maxPt;
}

float BasicRatePlots::SingleMuErPt()  {

  float maxPt = -10;
  
  int Nmu = gmt_ -> N;
  for (int imu=0; imu < Nmu; imu++) 
    {
      int bx = gmt_ -> CandBx[imu];       // BX = 0, +/- 1 or +/- 2
      if(bx != 0) continue;
      int qual = gmt_ -> Qual[imu];
      if(qual < 4) continue;
      float eta = gmt_ -> Eta[imu];
      if(fabs(eta) > 2.1) continue;
      float pt = gmt_ -> Pt[imu];         // get pt to get highest pt one
      if( pt > maxPt) maxPt = pt;
    }
  
  return maxPt;
}

float BasicRatePlots::SingleMuEta(float ptCut ) {

  float maxPt = -10;
  float iMuMaxPt = -10;
  
  int Nmu = gmt_ -> N;
  for (int imu=0; imu < Nmu; imu++) 
    {
      int bx = gmt_ -> CandBx[imu];       // BX = 0, +/- 1 or +/- 2
      if (bx != 0) continue;
      int qual = gmt_ -> Qual[imu];
      if (qual < 4) continue;
      float pt = gmt_ -> Pt[imu];         // get pt to get highest pt one
      if ( pt > maxPt) 
	{
	  maxPt = pt;
	  iMuMaxPt = imu;
	}
    }
  
  float eta = iMuMaxPt>=0 && maxPt>ptCut ? gmt_ -> Eta[iMuMaxPt] : -10; 
  return eta;
}

float BasicRatePlots::DttfPt()  {

  float maxPt = -10;
  
  int Nmu = gmt_ -> Ndt;
  for (int imu=0; imu < Nmu; imu++) 
    {
      int bx = gmt_ -> Bxdt[imu];           // BX = 0, +/- 1 or +/- 2
      if (bx != 0) continue;
      float pt = gmt_ -> Ptdt[imu];         // get pt to get highest pt one
      if ( pt > maxPt) maxPt = pt;
    }
  
  return maxPt;
}

float BasicRatePlots::RpcbPt(){

  float maxPt = -10;
  
  int Nmu = gmt_ -> Nrpcb;
  for (int imu=0; imu < Nmu; imu++) 
    {
      int bx = gmt_ -> Bxrpcb[imu];           // BX = 0, +/- 1 or +/- 2
      if (bx != 0) continue;
      float pt = gmt_ -> Ptrpcb[imu];         // get pt to get highest pt one
      if ( pt > maxPt)  maxPt = pt;
    }
  
  return maxPt;
}

float BasicRatePlots::CsctfPt(){

  float maxPt = -10;
  
  int Nmu = gmt_ -> Ncsc;
  for (int imu=0; imu < Nmu; imu++) 
    {
      int bx = gmt_ -> Bxcsc[imu];           // BX = 0, +/- 1 or +/- 2
      if (bx != 0) continue;
      float pt = gmt_ -> Ptcsc[imu];         // get pt to get highest pt one
      if ( pt > maxPt) maxPt = pt;
    }
  
  return maxPt;
}

float BasicRatePlots::RpcfPt()  {

  float maxPt = -10;
  
  int Nmu = gmt_ -> Nrpcf;
  for (int imu=0; imu < Nmu; imu++) 
    {
      int bx = gmt_ -> Bxrpcf[imu];           // BX = 0, +/- 1 or +/- 2
      if (bx != 0) continue;
      float pt = gmt_ -> Ptrpcf[imu];         // get pt to get highest pt one
      if ( pt > maxPt) maxPt = pt;
    }
  
  return maxPt;
}

void BasicRatePlots::DoubleMu(Float_t & muPt1, Float_t & muPt2) {

  muPt1 = -10.;
  muPt2 = -10.;
  Int_t Nmu = gmt_ -> N;
  for (Int_t imu=0; imu < Nmu; imu++) {
    Int_t bx = gmt_ -> CandBx[imu];		
    if (bx != 0) continue;
    Float_t pt = gmt_ -> Pt[imu];			
    Int_t qual = gmt_ -> Qual[imu];        
    if (qual < 4  && qual != 3 ) continue;
    if (pt >= muPt1)
      {
	muPt2 = muPt1;
	muPt1 = pt;
      }
    else if (pt >= muPt2) muPt2 = pt;
  }
  
}


void BasicRatePlots::Onia(Float_t & muPt1, Float_t & muPt2, Float_t etacut, Int_t delta) {

  std::vector<std::pair<Float_t,Float_t> > muonPairs;

  Int_t Nmu = gmt_ -> N;
  for (Int_t imu=0; imu < Nmu; imu++) {
    Int_t bx = gmt_ -> CandBx[imu];		
    if (bx != 0) continue;
    Float_t pt = gmt_ -> Pt[imu];			
    Int_t qual = gmt_ -> Qual[imu];        
    if ( qual < 4) continue;
    Float_t eta = gmt_  -> Eta[imu];        
    if (fabs(eta) > etacut) continue;
    Int_t ieta1 = etaMuIdx(eta);
    
    for (Int_t imu2=0; imu2 < Nmu; imu2++) {
      if (imu2 == imu) continue;
      Int_t bx2 = gmt_ -> CandBx[imu2];		
      if (bx2 != 0) continue;
      Float_t pt2 = gmt_ -> Pt[imu2];			
      Int_t qual2 = gmt_ -> Qual[imu2];        
      if ( qual2 < 4) continue;
      Float_t eta2 = gmt_  -> Eta[imu2];        
      if (fabs(eta2) > etacut) continue;
      Int_t ieta2 = etaMuIdx(eta2);
      
      bool hasCorrCond = ( fabs(ieta1 - ieta2) <= delta);

      if (hasCorrCond) muonPairs.push_back(std::pair<Float_t,Float_t>(pt,pt2));
      
    }
    
  }
  
  // CB loop on muon pairs and get highest pT pair
  // (in order look for first and second muon, can be arbitrary, it's a choice)
  
  muPt1 = -10;
  muPt2 = -10;
  
  std::vector<std::pair<Float_t,Float_t> >::const_iterator muonPairIt  = muonPairs.begin();
  std::vector<std::pair<Float_t,Float_t> >::const_iterator muonPairEnd = muonPairs.end();
  for (; muonPairIt != muonPairEnd; ++muonPairIt)
    {
      Float_t pt1 = muonPairIt->first;
      Float_t pt2 = muonPairIt->second;
      
      if ( (pt1 > muPt1) ||
	   (fabs(muPt1-pt1)<10E-2 && pt2>muPt2) ) 
	{
	  muPt1 = pt1;
	  muPt2 = pt2;
	}
    }

}


float BasicRatePlots::SingleJetCentralPt() {
  
  float maxPt = -10;

  Int_t Nj = gt_ -> Njet ;
  for (Int_t ue=0; ue < Nj; ue++) {
    Int_t bx = gt_ -> Bxjet[ue];        		
    if (bx != 0) continue; 
    Bool_t isFwdJet = gt_ -> Fwdjet[ue];
    if (isFwdJet) continue;
    //Bool_t isTauJet = gt_ -> Taujet[ue];
    //if (isTauJet) continue;
    Float_t rank = gt_ -> Rankjet[ue];
    Float_t pt = rank*4.;
    if (pt >= maxPt) maxPt = pt;
  } 
  
  return maxPt;  
}

float BasicRatePlots::SingleTauPt() {
  
  float maxPt = -10;

  Int_t Nj = gt_ -> Njet ;
  for (Int_t ue=0; ue < Nj; ue++) {
    Int_t bx = gt_ -> Bxjet[ue];        		
    if (bx != 0) continue; 
    Bool_t isTauJet = gt_ -> Taujet[ue];
    if (!isTauJet) continue;
    Float_t rank = gt_ -> Rankjet[ue];
    Float_t pt = rank*4.;
    if (pt >= maxPt) maxPt = pt;
  } 
  
  return maxPt;
}

float BasicRatePlots::SingleJetPt() {

  float maxPt = -10;
  Int_t Nj = gt_ -> Njet ;
  for (Int_t ue=0; ue < Nj; ue++) {
    Int_t bx = gt_ -> Bxjet[ue];        		
    if (bx != 0) continue;
    Float_t rank = gt_ -> Rankjet[ue];
    Float_t pt = rank*4.;
    if (pt >= maxPt) maxPt = pt;
  }
  
  return maxPt;  
}

inline float BasicRatePlots::HTT() {
  
  Float_t adc = gt_ -> RankHTT ;
  Float_t theHTT = adc / 2. ;
  
  return theHTT;
}

inline Float_t BasicRatePlots::ETT() {

  Float_t adc = gt_ -> RankETT ;
  Float_t theETT = adc / 2. ;
  
  return theETT;  
}

inline float BasicRatePlots::ETM() {

  Float_t adc = gt_ -> RankETM ;
  Float_t TheETM = adc / 2. ;

  return TheETM;
}

float BasicRatePlots::SingleEGPt() {

  float maxPt = -10; 

  Int_t Nele = gt_ -> Nele;
  for (Int_t ue=0; ue < Nele; ue++) {               
    Int_t bx = gt_ -> Bxel[ue];        		
    if (bx != 0) continue;
    Float_t pt = gt_ -> Rankel[ue];    // the rank of the electron
    if (pt >= maxPt) maxPt = pt;
  }  // end loop over EM objects
  
  return maxPt;   
}


float BasicRatePlots::SingleIsoEGPt() {

  float maxPt = -10; 

  Int_t Nele = gt_ -> Nele;
  for (Int_t ue=0; ue < Nele; ue++) {               
    Int_t bx = gt_ -> Bxel[ue];        		
    if (bx != 0) continue;
    Bool_t iso = gt_ -> Isoel[ue];
    if (! iso) continue;
    Float_t pt = gt_ -> Rankel[ue];    // the rank of the electron
    if (pt >= maxPt) maxPt = pt;
  }  // end loop over EM objects
  
  return maxPt;
}


void BasicRatePlots::EGIsoEGPt(Float_t &myIsoEGPt, Float_t &myEGPt) {

  float maxPtEGIso = -10.; 
  float maxPtEG    = -10.; 

  Int_t iEGIso = -1;

  Int_t Nele = gt_ -> Nele;
  for (Int_t ue=0; ue < Nele; ue++) {
    Int_t bx = gt_ -> Bxel[ue];        		
    if (bx != 0) continue;
    Bool_t iso = gt_ -> Isoel[ue];
    if (! iso) continue;
    Float_t ptIso = gt_ -> Rankel[ue];    // the rank of the first (isolated) electron
    if (ptIso >= maxPtEGIso){
      maxPtEGIso = ptIso;
      iEGIso = ue;
    }
  }

  for (Int_t ue2=0; ue2 < Nele; ue2++) {
    if(ue2 == iEGIso) continue;
    Int_t bx2 = gt_ -> Bxel[ue2];
    if (bx2 != 0) continue;

    Float_t pt2 = gt_ -> Rankel[ue2];    // the rank of the second electron
    if (pt2 >= maxPtEG) maxPtEG = pt2;
  }
  
  myIsoEGPt = maxPtEGIso;
  myEGPt = maxPtEG;

  return;
}


void BasicRatePlots::Mu_HTT(Float_t & mymuPt, Float_t & myHTT) {

  Float_t adc = gt_ -> RankHTT ;
  myHTT = adc / 2. ;

  Float_t maxPt = -1.;
  int Nmu = gmt_ -> N;
  for (int imu=0; imu < Nmu; imu++) 
    {
      int bx = gmt_ -> CandBx[imu];       // BX = 0, +/- 1 or +/- 2
      if (bx != 0) continue;
      int qual = gmt_ -> Qual[imu];
      if (qual < 4) continue;
      float pt = gmt_ -> Pt[imu];         // get pt to get highest pt one
      if ( pt > maxPt) maxPt = pt;
    }

  mymuPt = maxPt;

  return;
}


vector<float> BasicRatePlots::SortJets(bool doCentral, bool doTau){

  if (doCentral && doTau){
    cout << "Incompatible arguements to SortJets -- Exiting macro" << endl;
    exit(1);
  }
  vector<float> theJets, theSortedJets;

  Int_t Nj = gt_ -> Njet ;
  for (Int_t ue=0; ue < Nj; ue++) {
    Int_t bx = gt_ -> Bxjet[ue];        		
    if (bx != 0) continue; 
    if (doCentral){
      Bool_t isFwdJet = gt_ -> Fwdjet[ue];
      if (isFwdJet) continue;
    }else if (doTau){
      Bool_t isTauJet = gt_ -> Taujet[ue];
      if (! isTauJet) continue;
    }
    Float_t rank = gt_ -> Rankjet[ue];
    Float_t pt = rank*4.;
    theJets.push_back(float(pt));
  } 
  //cout << "Size of jets: " << Nj << " central: " << theJets.size() << endl;
  
  theSortedJets=theJets;
  for (unsigned i=0; i < theJets.size(); i++) {
    if (theJets[i] != theSortedJets[i] ){
      cout << "\t ERROR XXXX: " << i << "\t " << theSortedJets[i] << "\t " << theJets[i] << "\t " << theJets.size() << "\t " << theSortedJets.size() << endl;
      //}else{
      //cout << "\t " << theSortedJets[i] << "\t " << theJets[i] << endl;
    }
  }
  std::stable_sort(theSortedJets.begin(),theSortedJets.end(),std::greater<float>());
  return theSortedJets;

}

vector<float> BasicRatePlots::SortEGs(bool doIsol){

  vector<float> theEGs, theSortedEGs;


  Int_t Nele = gt_ -> Nele;
  for (Int_t ue=0; ue < Nele; ue++) {               
    Int_t bx = gt_ -> Bxel[ue];        		
    if (bx != 0) continue;

    if (doIsol){
      Bool_t iso = gt_ -> Isoel[ue];
      if (!iso) continue;
    }

    Float_t pt = gt_ -> Rankel[ue];    // the rank of the electron
    theEGs.push_back(float(pt));
  }  // end loop over EM objects

  
  theSortedEGs=theEGs;
  for (unsigned i=0; i < theEGs.size(); i++) {
    if (theEGs[i] != theSortedEGs[i] ){
      cout << "\t ERROR XXXX: " << i << "\t " << theSortedEGs[i] << "\t " << theEGs[i] << "\t " << theEGs.size() << "\t " << theSortedEGs.size() << endl;
      //}else{
      //cout << "\t " << theSortedEGs[i] << "\t " << theEGs[i] << endl;
    }
  }
  std::stable_sort(theSortedEGs.begin(),theSortedEGs.end(),std::greater<float>());
  return theSortedEGs;

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


void BasicRatePlots::run(bool runOnData, std::string resultTag, int minLs, int maxLs, float crossSec, float avPU, int nBunches, int isCrossSec, int nEvents) {

  system("mkdir -p results");
  std::string resultName = "results_" + resultTag + (isCrossSec ? "_XSEC" : "_RATE") + ".root";
  TFile *outFile = new TFile(("results/" + resultName).c_str(),"recreate");
  outFile->cd();

  //Single stuff
  hTH1F["nJetVsPt"]    = new TH1F("nJetVsPt","SingleJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hTH1F["nTauVsPt"]    = new TH1F("nTauVsPt","SingleTau; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hTH1F["nJetCenVsPt"] = new TH1F("nJetCenVsPt","SingleJetCentral; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hTH1F["nEGVsPt"]     = new TH1F("nEGVsPt","SingleEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH1F["nIsoEGVsPt"]  = new TH1F("nIsoEGVsPt","SingleIsoEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH1F["nMuVsPt"]     = new TH1F("nMuVsPt","SingleMu; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hTH1F["nMuErVsPt"]   = new TH1F("nMuErVsPt","SingleMu |#eta|<2.1; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hTH1F["nMuVsEta"]    = new TH1F("nMuVsEta","nMuVsEta",24,-2.4,2.4);

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
  hTH2F["nEGIsoEGVsPt"]      = new TH2F("nEGIsoEGVsPt","IsoEle_Ele; p_{T} cut iso EG_{1}; p_{T} cut EG_{2}",65,-0.5,64.5,65,-0.5,64.5);

  //Sums
  hTH1F["nHTTVsHTT"] = new TH1F("nHTTVsHTT","HTT; HTT cut; rate [Hz]",512,-.5,511.5);
  hTH1F["nETTVsETT"] = new TH1F("nETTVsETT","ETT; ETT cut; rate [Hz]",512,-.5,511.5);
  hTH1F["nETMVsETM"] = new TH1F("nETMVsETM","ETM; ETM cut; rate [Hz]",512,-.5,511.5);

  //Cross  
  hTH2F["nMuVsHTT"]  = new TH2F("nMuVsHTT","Mu_HTT; p_{T} cut mu_{1}; HTT",61,-0.25,30.25,512,-.5,511.5);
  hTH2F["nMuVsEG"]   = new TH2F("nMuVsEG","Mu_EG; p_{T} cut mu; p_{T} cut EG",61,-0.25,30.25,65,-0.5,64.5);

  hTH1F["nDttfVsPt"]  = new TH1F("nDttfVsPt","DTTF; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hTH1F["nRpcbVsPt"]  = new TH1F("nRpcbVsPt","RPCb; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hTH1F["nRpcfVsPt"]  = new TH1F("nRpcfVsPt","RPCf; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hTH1F["nCsctfVsPt"] = new TH1F("nCsctfVsPt","CSCTF; p_{T} cut; rate [Hz]",131,-0.5,130.5);

  if (isCrossSec) {
    std::map<std::string,TH1F*>::iterator hTH1FIt  = hTH1F.begin();
    std::map<std::string,TH1F*>::iterator hTH1FEnd = hTH1F.end();

    for(; hTH1FIt!=hTH1FEnd; ++hTH1FIt)
      {
	hTH1FIt->second->GetYaxis()->SetTitle("cross section [#mubarn]");
      }
  }

  hTH1F["nPUvsPU"]  = new TH1F("nPUvsPU","Num. of PU iteractions; Num of iteractions; Reweighted couns [arb units]",101,-0.5,100.5);

  float nZeroBias = 0;

  int nevents = nEvents == 0 ? GetEntries() : nEvents;
    
  std::cout << "Running on " << nevents << " events." << std::endl;
  for (Long64_t event=0; event<nevents; ++event) { 
    Long64_t eventEntry = LoadTree(event); 
    if (eventEntry < 0) break;
    GetEntry(event);
      
    //if ( event_->lumi < FIRST_VALID_LS || event_->lumi > LAST_VALID_LS ) continue;
      
    if (event%200000 == 0) {
      std::cout << "Processed " << event << " events." << std::endl;
    }
      
    if ( event_->lumi < minLs || event_->lumi > maxLs ) continue;

    double weight = event_->puWeight > -0.001 ? event_->puWeight : 1; 
      
    FillBits();

    if(runOnData && !PhysicsBits[0]) continue;
      
    nZeroBias += weight;

    float jetPt     = SingleJetPt();
    float jetCenPt  = SingleJetCentralPt();

    float tauPt     = SingleTauPt();

    float htt       = HTT();
    float ett       = ETT();
    float etm       = ETM();

    float egPt      = SingleEGPt();
    float isoEgPt   = SingleIsoEGPt();

    float muPt     = SingleMuPt();
    float muErPt   = SingleMuErPt();

    float muEta    = SingleMuEta(16.);
      
    float doubleMuPt1 = -10;
    float doubleMuPt2 = -10;
    DoubleMu(doubleMuPt1,doubleMuPt2);

    float oniaMuPt1 = -10;
    float oniaMuPt2 = -10;
    Onia(oniaMuPt1,oniaMuPt2,2.1,22);

    float EGIsoPt1 = -10;
    float EGPt2    = -10;
    EGIsoEGPt(EGIsoPt1,EGPt2);

    float dttfPt   = DttfPt();
    float rpcbPt   = RpcbPt();
    float rpcfPt   = RpcfPt();
    float csctfPt  = CsctfPt();

    float muPt_forHTT = -10.;
    float HTT_forMu = -10.;
    Mu_HTT(muPt_forHTT,HTT_forMu);

    // creates sorted jet vectors for all jet cands and central jets
    sortedJets = SortJets(false,false);  // no cuts --> doCentral=false, doTau=false
    sortedCenJets = SortJets(true,false); // for central jet rate  --> doCentral=true, doTau=false
    sortedTaus = SortJets(false,true); // for tau jet rate  --> doCentral=false, doTau=true

    sortedEGs     = SortEGs(false);
    sortedIsolEGs = SortEGs(true);

    bool dijets    = sortedJets.size()>1;
    bool dijetsC   = sortedCenJets.size()>1;
    bool ditaus    = sortedTaus.size()>1;
    bool quadjets  = sortedJets.size()>3;
    bool quadjetsC = sortedCenJets.size()>3;

    bool diEG     = sortedEGs.size()>1;
    bool diIsolEG = sortedIsolEGs.size()>1;

    float dijetPt    = -10.;
    float diCenjetPt = -10.;
    float ditauPt    = -10.;
    float quadjetPt  = -10.;
    float quadjetCPt = -10.;
    if(dijets)    dijetPt     = sortedJets.at(1);
    if(dijetsC)   diCenjetPt  = sortedCenJets.at(1);
    if(ditaus)    ditauPt     = sortedTaus.at(1);
    if(quadjets)  quadjetPt   = sortedJets.at(3);
    if(quadjetsC) quadjetCPt  = sortedCenJets.at(3);

    hTH1F["nPUvsPU"]->Fill(simulation_->actualInt,weight);
    hTH1F["nMuVsEta"]->Fill(muEta,weight);

    for(int ptCut=0; ptCut<256; ++ptCut) {
      if(jetPt>ptCut)	   hTH1F["nJetVsPt"]->Fill(ptCut,weight);
      if(jetCenPt>ptCut) hTH1F["nJetCenVsPt"]->Fill(ptCut,weight);
      if(tauPt>ptCut)	   hTH1F["nTauVsPt"]->Fill(ptCut,weight);

      if(dijetPt>ptCut){
	hTH1F["nDiJetVsPt"]->Fill(ptCut,weight);

	for(int ptCut_0=ptCut; ptCut_0<256; ++ptCut_0) {
	  if(sortedJets.at(0)>ptCut_0) hTH2F["nAsymDiJetVsPt"]->Fill(ptCut_0,ptCut,weight);
	}
      }

      if(diCenjetPt>ptCut){
	hTH1F["nDiCenJetVsPt"]->Fill(ptCut,weight);

	for(int ptCut_0=ptCut; ptCut_0<256; ++ptCut_0) {
	  if(sortedCenJets.at(0)>ptCut_0) hTH2F["nAsymDiCenJetVsPt"]->Fill(ptCut_0,ptCut,weight);
	}
      }

      if(ditauPt>ptCut)    hTH1F["nDiTauVsPt"]->Fill(ptCut,weight);
      if(quadjetPt>ptCut)  hTH1F["nQuadJetVsPt"]->Fill(ptCut,weight);
      if(quadjetCPt>ptCut) hTH1F["nQuadCenJetVsPt"]->Fill(ptCut,weight);

    }//loop on 256
      
    for(int ptCut=0; ptCut<65; ++ptCut) {
      if(egPt>ptCut)    hTH1F["nEGVsPt"]->Fill(ptCut,weight);
      if(isoEgPt>ptCut) hTH1F["nIsoEGVsPt"]->Fill(ptCut,weight);

      if(diEG && sortedEGs[1]>ptCut)         hTH1F["nDiEGVsPt"]->Fill(ptCut,weight);
      if(diIsolEG && sortedIsolEGs[1]>ptCut) hTH1F["nDiIsoEGVsPt"]->Fill(ptCut,weight);


      for(int ptCut2=0; ptCut2<=65; ++ptCut2) {
	if(diEG && sortedEGs[0]>ptCut && sortedEGs[1]>ptCut2 && ptCut2 <= ptCut) hTH2F["nEGPtVsPt"]->Fill(ptCut,ptCut2,weight);
	if(diIsolEG && sortedIsolEGs[0]>ptCut && sortedIsolEGs[1]>ptCut2 && ptCut2<= ptCut) hTH2F["nIsoEGPtVsPt"]->Fill(ptCut,ptCut2,weight);
	if(EGIsoPt1>ptCut && EGPt2>ptCut2) hTH2F["nEGIsoEGVsPt"]->Fill(ptCut,ptCut2,weight);
      }

    }//loop on 65
      
    for(int ptCut=0; ptCut<131; ++ptCut) {
      if (muPt>=ptCut)    hTH1F["nMuVsPt"]->Fill(ptCut,weight);
      if (muErPt>=ptCut)  hTH1F["nMuErVsPt"]->Fill(ptCut,weight);
      if (dttfPt>=ptCut)  hTH1F["nDttfVsPt"]->Fill(ptCut,weight);
      if (rpcbPt>=ptCut)  hTH1F["nRpcbVsPt"]->Fill(ptCut,weight);
      if (rpcfPt>=ptCut)  hTH1F["nRpcfVsPt"]->Fill(ptCut,weight);
      if (csctfPt>=ptCut) hTH1F["nCsctfVsPt"]->Fill(ptCut,weight);
    }
      

    for(int iCut=0; iCut<41; ++iCut) {
      for(int iCut2=0; iCut2<=iCut; ++iCut2) {
	float ptCut = iCut*0.5;
	float ptCut2 = iCut2*0.5;
	if (doubleMuPt1>ptCut && doubleMuPt2>ptCut2) hTH2F["nMuPtVsPt"]->Fill(ptCut,ptCut2,weight);
	if (oniaMuPt1>ptCut && oniaMuPt2>ptCut2)     hTH2F["nOniaMuPtVsPt"]->Fill(ptCut,ptCut2,weight);
      }
    }

    for(int iCut=0; iCut<61; ++iCut) {
      float ptCutMu = iCut*0.5;

      for(int ptCutEG=0; ptCutEG<65; ++ptCutEG){
	if(muPt>ptCutMu && egPt>ptCutEG) hTH2F["nMuVsEG"]->Fill(ptCutMu,ptCutEG,weight);
      }
    }

      
    for(int httCut=0; httCut<512; ++httCut) {
      if(htt>httCut) hTH1F["nHTTVsHTT"]->Fill(httCut,weight);
      if(ett>httCut) hTH1F["nETTVsETT"]->Fill(httCut,weight);
      if(etm>httCut) hTH1F["nETMVsETM"]->Fill(httCut,weight);

      for(int iCut=0; iCut<61; ++iCut) {
      float ptCutMu = iCut*0.5;
	if (muPt_forHTT>ptCutMu && HTT_forMu>httCut) hTH2F["nMuVsHTT"]->Fill(ptCutMu,httCut,weight);
      }
    }
      

  } // end event loop

  cout << "# of zero bias events (weighted) used for rate computation : " << nZeroBias << std::endl;

  float scaleFactor = ScaleFactor(nZeroBias,nBunches);

  if (isCrossSec) 
    scaleFactor /= (computeAvgLumi(crossSec,avPU,nBunches)*10000) ; // CB lumi is in 1E34 units
  
  map<string,TH1F*>::iterator hTH1FIt  = hTH1F.begin();
  map<string,TH1F*>::iterator hTH1FEnd = hTH1F.end();

  for(; hTH1FIt!=hTH1FEnd; ++hTH1FIt) {
    TH1F* histo = hTH1FIt->second;
    if (hTH1FIt->first == "nPUvsPU") histo->Scale(1./nZeroBias);      
    setRateError(histo);
    histo->Scale(scaleFactor);
  }

  map<string,TH2F*>::iterator hTH2FIt  = hTH2F.begin();
  map<string,TH2F*>::iterator hTH2FEnd = hTH2F.end();

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

  float xSec13TeV = isCrossSec ? 78.26 : 80.; // Using McM for cross section comparison and 80 (agreed with TSG) for rates
  float xSec8TeV  = 72.7; 

  if (fileType == "DATA") {
      BasicRatePlots basicRatePlots("/afs/cern.ch/user/h/heistera/scratch1/L1Ntuples/L1TreeL1Accept_207477_LS_57_133.root");
      basicRatePlots.run(true,"DATA_207477",57,133,xSec8TeV,999.,nBunches50ns,isCrossSec,nEvents); // 999 is dummy do not use for cross-section
    }
  else if (fileType == "13TEV_25PU_ORIG_RE-EMUL")
    {
      BasicRatePlots basicRatePlots("/data2/battilan/L1Trigger/L1T2015Menu/L1Tree_13TeV_25PU_53X_OrigEmul_v4.root");
      basicRatePlots.run(false,"13TEV_25PU_ORIG_RE-EMUL",0,500000000,xSec13TeV,25,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "13TEV_25PU_2012_RE-EMUL" )
    {
      BasicRatePlots basicRatePlots("/data2/battilan/L1Trigger/L1T2015Menu/L1Tree_13TeV_25PU_53X_ReEmul2012_v4.root");
      basicRatePlots.run(false,"13TEV_25PU_2012_RE-EMUL",0,500000000,xSec13TeV,25,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "13TEV_25PU_2012GCT10GEV_RE-EMUL" )
    {
      BasicRatePlots basicRatePlots("/data2/battilan/L1Trigger/L1T2015Menu/L1Tree_13TeV_25PU_53X_ReEmul2012Gct10GeV_v4.root");
      basicRatePlots.run(false,"13TEV_25PU_2012GCT10GEV_RE-EMUL",0,500000000,xSec13TeV,25,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "13TEV_25PU_2015_RE-EMUL")
    {
      BasicRatePlots basicRatePlots("/data2/battilan/L1Trigger/L1T2015Menu/L1Tree_13TeV_25PU_53X_ReEmul2015_v4.root"); 
      basicRatePlots.run(false,"13TEV_25PU_2015_RE-EMUL",0,500000000,xSec13TeV,25,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "8TEV_TF_2012_RE-EMUL")
    {
      BasicRatePlots basicRatePlots("/data2/battilan/L1Trigger/L1T2015Menu/L1Ntuple_8TeV_53X_ReEmul2012_v2.root"); 
      basicRatePlots.run(false,"8TEV_TF_2012_RE-EMUL",0,500000000,xSec8TeV,999.,nBunches50ns,isCrossSec,nEvents);  // 999 is dummy do not use for cross-section
    }
  else if (fileType == "8TEV_TF_DATA")
    {
      BasicRatePlots basicRatePlots("/data2/battilan/L1Trigger/L1T2015Menu/L1Ntuple_8TeV_53X_202299_v2.root"); 
      basicRatePlots.run(true,"TF_DATA_202299",70,550,xSec8TeV,999.,nBunches50ns,isCrossSec,nEvents); // CB need to find PU for lumi calculation
    }
  else if (fileType == "8TEV_25PU_ORIG_RE-EMUL")
    {
      BasicRatePlots basicRatePlots("/data2/battilan/L1Trigger/L1T2015Menu/L1Tree_8TeV_25PU_53X_OrigEmul_v4.root"); 
      basicRatePlots.run(false,"8TEV_25PU_ORIG_RE-EMUL",0,500000000,xSec8TeV,25,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "8TEV_25PU_2012_RE-EMUL")
    {
      BasicRatePlots basicRatePlots("/data2/battilan/L1Trigger/L1T2015Menu/L1Tree_8TeV_25PU_53X_ReEmul2012_v4.root"); 
      basicRatePlots.run(false,"8TEV_25PU_2012_RE-EMUL",0,500000000,xSec8TeV,25,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "8TEV_25PU_2012GCT10GEV_RE-EMUL")
    {
      BasicRatePlots basicRatePlots("/data2/battilan/L1Trigger/L1T2015Menu/L1Tree_8TeV_25PU_53X_ReEmul2012Gct10GeV_v4.root"); 
      basicRatePlots.run(false,"8TEV_25PU_2012GCT10GEV_RE-EMUL",0,500000000,xSec8TeV,25,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "8TEV_25PU_2015_RE-EMUL")
    {
      BasicRatePlots basicRatePlots("/data2/battilan/L1Trigger/L1T2015Menu/L1Tree_8TeV_25PU_53X_ReEmul2015_v4.root"); 
      basicRatePlots.run(false,"8TEV_25PU_2015_RE-EMUL",0,500000000,xSec8TeV,25,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "13TEV_40PU_2012_RE-EMUL")
    {
      BasicRatePlots basicRatePlots("/data2/battilan/L1Trigger/L1T2015Menu/L1Tree_v5_62X_13TeV_40PU_25bx_ReEmul2012.root"); 
      basicRatePlots.run(false,"13TEV_40PU_2012_RE-EMUL",0,500000000,xSec13TeV,40,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "13TEV_40PU_2012GCT10GEV_RE-EMUL")
    {
      BasicRatePlots basicRatePlots("/data2/p/pellicci/L1DPG/root/v4_62X_40PU_25bx_ReEmul2012Gct10GeV/L1Tree.root"); 
      basicRatePlots.run(false,"13TEV_40PU_2012_RE-EMUL",0,500000000,xSec13TeV,40,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "13TEV_40PU_2015_RE-EMUL")
    {
      BasicRatePlots basicRatePlots("/data2/p/pellicci/L1DPG/root/v4_62X_40PU_25bx_ReEmul2015/L1Tree_v3.root"); 
      basicRatePlots.run(false,"13TEV_40PU_2015_RE-EMUL",0,500000000,xSec13TeV,40,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "13TEV_45p4PU_2012_RE-EMUL")
    {
      BasicRatePlots basicRatePlots("/data2/battilan/L1Trigger/L1T2015Menu/L1Tree_v5_62X_13TeV_45p4PU_25bx_ReEmul2012.root"); 
      basicRatePlots.run(false,"13TEV_45p4PU_2012_RE-EMUL",0,500000000,xSec13TeV,45,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "TEST")
    {
      BasicRatePlots basicRatePlots("../test/L1Tree.root"); 
      basicRatePlots.run(false,"TEST",0,500000000,xSec13TeV,25,nBunches25ns,isCrossSec,nEvents);
    }
  else 
    {
      std::cout << "Config param " << fileType << " invalid! \n"
		<< "Valid fileType values are : DATA, 8TEV_TF_DATA, 8TEV_TF_2012_RE-EMUL, "
		<< "8TEV_25PU_ORIG_RE-EMUL, 8TEV_25PU_2012_RE-EMUL, 8TEV_25PU_2012GCT10GEV_RE-EMUL, 8TEV_25PU_2015_RE-EMUL, "
		<< "13TEV_25PU_ORIG_RE-EMUL, 13TEV_25PU_2012_RE-EMUL, 13TEV_25PU_2012GCT10GEV_RE-EMUL, 13TEV_25PU_2015_RE-EMUL\n";
    }
    
}


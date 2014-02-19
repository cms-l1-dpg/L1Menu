#include "../../L1Ntuples/macros/L1Ntuple.h"

size_t ETAMUBINS = 65;
Double_t ETAMU[] = { -2.45,-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45 };

Int_t etaMuIdx(Double_t eta) {
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
  float SingleJetCentralPt();
  float SingleJetCentralEta(Float_t cut );
  float SingleJetPt();

  float HTT();
  float ETT();
  
  float SingleEGPt();
  float SingleIsoEGPt();

  float SingleMuPt();
  float SingleMuErPt();

  void DoubleMu(Float_t & muPt1, Float_t & muPt2);
  void Onia(Float_t & muPt1, Float_t & muPt2, Float_t etacut, Int_t delta);

  float SingleMuEta(float eta);

  float DttfPt();
  float RpcbPt();
  float RpcfPt();
  float CsctfPt();

  void setRateError(TH1F* histo);

  float computeAvgLumi(float xSec, float avPU, int nBunches) { return 11246. * avPU * nBunches / (1E7 * xSec); };

  bool PhysicsBits[128];

  std::map<std::string,TH1F*> hTH1F;
  std::map<std::string,TH2F*> hTH2F;

};

void BasicRatePlots::FillBits() {

    //      Fill the physics bits:

    for (int ibit=0; ibit < 128; ibit++) 
    {
        PhysicsBits[ibit] = 0;
        if (ibit<64) 
        {
            PhysicsBits[ibit] = (gt_->tw1[2]>>ibit)&1;
        }
        else 
        {
            PhysicsBits[ibit] = (gt_->tw2[2]>>(ibit-64))&1;
        }
    }

}       

// ------------------------------------------------------------------

// scale factor computed w.r.t. ZeroBias rate fratcion and # bunches 
float BasicRatePlots::ScaleFactor(float nZeroBias, float nBunches) {

  float scal = 11246.; // ZB per bunch in kHz
  scal /= nZeroBias;
  scal *= nBunches;

  return scal;

}

float BasicRatePlots::SingleMuPt()  {

  float maxPt = -10;
  
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
	}
    }
  
  return maxPt;

}

float BasicRatePlots::SingleMuErPt()  {

  float maxPt = -10;
  
  int Nmu = gmt_ -> N;
  for (int imu=0; imu < Nmu; imu++) 
    {
      int bx = gmt_ -> CandBx[imu];       // BX = 0, +/- 1 or +/- 2
      if (bx != 0) continue;
      int qual = gmt_ -> Qual[imu];
      if (qual < 4) continue;
      float eta = gmt_ -> Eta[imu];
      if (fabs(eta) > 2.1) continue;
      float pt = gmt_ -> Pt[imu];         // get pt to get highest pt one
      if ( pt > maxPt) 
	{
	  maxPt = pt;
	}
    }
  
  return maxPt;

}

float BasicRatePlots::SingleMuEta(float ptCut )  {

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
      if ( pt > maxPt) 
	{
	  maxPt = pt;
	}
    }
  
  return maxPt;

}

float BasicRatePlots::RpcbPt()  {

  float maxPt = -10;
  
  int Nmu = gmt_ -> Nrpcb;
  for (int imu=0; imu < Nmu; imu++) 
    {
      int bx = gmt_ -> Bxrpcb[imu];           // BX = 0, +/- 1 or +/- 2
      if (bx != 0) continue;
      float pt = gmt_ -> Ptrpcb[imu];         // get pt to get highest pt one
      if ( pt > maxPt) 
	{
	  maxPt = pt;
	}
    }
  
  return maxPt;

}

float BasicRatePlots::CsctfPt()  {

  float maxPt = -10;
  
  int Nmu = gmt_ -> Ncsc;
  for (int imu=0; imu < Nmu; imu++) 
    {
      int bx = gmt_ -> Bxcsc[imu];           // BX = 0, +/- 1 or +/- 2
      if (bx != 0) continue;
      float pt = gmt_ -> Ptcsc[imu];         // get pt to get highest pt one
      if ( pt > maxPt) 
	{
	  maxPt = pt;
	}
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
      if ( pt > maxPt) 
	{
	  maxPt = pt;
	}
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

      if (hasCorrCond) muonPairs.push_back(std::make_pair<Float_t& ,Float_t&>(pt,pt2));
      
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
    Bool_t isTauJet = gt_ -> Taujet[ue];
    if (isTauJet) continue;
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

float BasicRatePlots::HTT() {
  
  Float_t adc = gt_ -> RankHTT ;
  Float_t theHTT = adc / 2. ;
  
  return theHTT;
  
}

Float_t BasicRatePlots::ETT() {

  Float_t adc = gt_ -> RankETT ;
  Float_t theETT = adc / 2. ;
  
  return theETT;
  
}


float BasicRatePlots::SingleEGPt() {

  float maxPt = -10; 

  Int_t Nele = gt_ -> Nele;
  for (Int_t ue=0; ue < Nele; ue++) {               
    Int_t bx = gt_ -> Bxel[ue];        		
    if (bx != 0) continue;
    Float_t rank = gt_ -> Rankel[ue];    // the rank of the electron
    Float_t pt = rank ; 
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
    Float_t rank = gt_ -> Rankel[ue];    // the rank of the electron
    Float_t pt = rank ; 
    if (pt >= maxPt) maxPt = pt;
  }  // end loop over EM objects
  
  return maxPt; 
  
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

  hTH1F["nJetVsPt"]    = new TH1F("nJetVsPt","SingleJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hTH1F["nJetCenVsPt"] = new TH1F("nJetCenVsPt","SingleJetCentral; E_{T} cut; rate [Hz]",256,-0.5,255.5);

  hTH1F["nHTTVsHTT"] = new TH1F("nHTTVsHTT","HTT; HTT cut; rate [Hz]",512,-.5,511.5);
  hTH1F["nETTVsETT"] = new TH1F("nETTVsETT","ETT; ETT cut; rate [Hz]",512,-.5,511.5);

  hTH1F["nEGVsPt"]    = new TH1F("nEGVsPt","SingleEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH1F["nIsoEGVsPt"] = new TH1F("nIsoEGVsPt","SingleIsoEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);

  hTH1F["nMuVsPt"]   = new TH1F("nMuVsPt","SingleMu; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hTH1F["nMuErVsPt"] = new TH1F("nMuErVsPt","SingleMu |#eta|<2.1; p_{T} cut; rate [Hz]",131,-0.5,130.5);

  hTH1F["nMuVsEta"] = new TH1F("nMuVsEta","nMuVsEta",24,-2.4,2.4);

  hTH1F["nDttfVsPt"]  = new TH1F("nDttfVsPt","DTTF; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hTH1F["nRpcbVsPt"]  = new TH1F("nRpcbVsPt","RPCb; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hTH1F["nRpcfVsPt"]  = new TH1F("nRpcfVsPt","RPCf; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hTH1F["nCsctfVsPt"] = new TH1F("nCsctfVsPt","CSCTF; p_{T} cut; rate [Hz]",131,-0.5,130.5);

  hTH2F["nMuPtVsPt"]    = new TH2F("nMuPtVsPt","DoubleMu; p_{T} cut mu_{1}; p_{T} cut mu_{2}",41,-0.25,20.25,41,-0.25,20.25);
  hTH2F["nOniaMuPtVsPt"]= new TH2F("nOniaMuPtVsPt","DoubleMu_Er_HighQ_WdEta22 (Quarkonia); p_{T} cut mu_{1}; p_{T} cut mu_{2}",41,-0.25,20.25,41,-0.25,20.25);
  
  if (isCrossSec)
    {
      std::map<std::string,TH1F*>::iterator hTH1FIt  = hTH1F.begin();
      std::map<std::string,TH1F*>::iterator hTH1FEnd = hTH1F.end();

      for(; hTH1FIt!=hTH1FEnd; ++hTH1FIt)
	{
	  hTH1FIt->second->GetYaxis()->SetTitle("cross section [#mubarn]");
	}
    }
      

  float nZeroBias = 0;

  int nevents = nEvents == 0 ? GetEntries() : nEvents;
    
  std::cout << "Running on " << nevents << " events." << std::endl;
  for (Long64_t event=0; event<nevents; ++event)
    { 
      Long64_t eventEntry = LoadTree(event); 
      if (eventEntry < 0) break;
      GetEntry(event);
      
      //if ( event_->lumi < FIRST_VALID_LS || event_->lumi > LAST_VALID_LS ) continue;
      
      if (event%200000 == 0) {
	std::cout << "Processed " << event << " events." << std::endl;
      }
      
      if ( event_->lumi < minLs || event_->lumi > maxLs ) continue;
      
      FillBits();

      if(runOnData && !PhysicsBits[0]) continue;
      
      nZeroBias += 1 ;

      float jetPt     = SingleJetPt();
      float jetCenPt  = SingleJetCentralPt();

      float htt       = HTT();
      float ett       = ETT();

      float egPt    = SingleEGPt();
      float isoEgPt = SingleIsoEGPt();

      float muPt     = SingleMuPt();
      float muErPt   = SingleMuErPt();

      float muEta    = SingleMuEta(16.);
      
      float doubleMuPt1 = -10;
      float doubleMuPt2 = -10;

      DoubleMu(doubleMuPt1,doubleMuPt2);

      float oniaMuPt1 = -10;
      float oniaMuPt2 = -10;

      Onia(oniaMuPt1,oniaMuPt2,2.1,22);

      float dttfPt   = DttfPt();
      float rpcbPt   = RpcbPt();
      float rpcfPt   = RpcfPt();
      float csctfPt  = CsctfPt();

      hTH1F["nMuVsEta"]->Fill(muEta);

      for(int ptCut=0; ptCut<256; ++ptCut) {
	if(jetPt>ptCut)
	  hTH1F["nJetVsPt"]->Fill(ptCut);
	if(jetCenPt>ptCut)
	  hTH1F["nJetCenVsPt"]->Fill(ptCut);
      }

      for(int ptCut=0; ptCut<65; ++ptCut) {
	if(egPt>ptCut)
	  hTH1F["nEGVsPt"]->Fill(ptCut);
	if(isoEgPt>ptCut)
	  hTH1F["nIsoEGVsPt"]->Fill(ptCut);
      }
      
      for(int ptCut=0; ptCut<131; ++ptCut) {
	if (muPt>ptCut)
	  hTH1F["nMuVsPt"]->Fill(ptCut);
	if (muErPt>ptCut)
	  hTH1F["nMuErVsPt"]->Fill(ptCut);
	if (dttfPt>ptCut)
	  hTH1F["nDttfVsPt"]->Fill(ptCut);
	if (rpcbPt>ptCut)
	  hTH1F["nRpcbVsPt"]->Fill(ptCut);
	if (rpcfPt>ptCut)
	  hTH1F["nRpcfVsPt"]->Fill(ptCut);
	if (csctfPt>ptCut)
	  hTH1F["nCsctfVsPt"]->Fill(ptCut);
      }

      for(int iCut=0; iCut<41; ++iCut) {
	for(int iCut2=0; iCut2<=iCut; ++iCut2) {
	  float ptCut = iCut*0.5;
	  float ptCut2 = iCut2*0.5;
	  if (doubleMuPt1>ptCut && doubleMuPt2>ptCut2)
	    hTH2F["nMuPtVsPt"]->Fill(ptCut,ptCut2);
	  if (oniaMuPt1>ptCut && oniaMuPt2>ptCut2)
	    hTH2F["nOniaMuPtVsPt"]->Fill(ptCut,ptCut2);
	}
      }
      
      for(int httCut=0; httCut<512; ++httCut) {
	if(htt>httCut)
	  hTH1F["nHTTVsHTT"]->Fill(httCut);
      }
      
      for(int ettCut=0; ettCut<512; ++ettCut) {
	if(ett>ettCut)
	  hTH1F["nETTVsETT"]->Fill(ettCut);
      }  

    } // end event loop

  cout << "# of zero bias events used for rate computation : " << nZeroBias << std::endl;

  float scaleFactor = ScaleFactor(nZeroBias,nBunches);

  if (isCrossSec) 
    scaleFactor /= (computeAvgLumi(crossSec,avPU,nBunches)*10000) ; // CB lumi is in 1E34 units

  map<string,TH1F*>::iterator hTH1FIt  = hTH1F.begin();
  map<string,TH1F*>::iterator hTH1FEnd = hTH1F.end();

  for(; hTH1FIt!=hTH1FEnd; ++hTH1FIt) 
    {

      TH1F* histo = hTH1FIt->second;
      setRateError(histo);
      histo->Scale(scaleFactor);

    }

  hTH2F["nMuPtVsPt"]->Scale(scaleFactor);
  hTH2F["nOniaMuPtVsPt"]->Scale(scaleFactor);
  
  outFile->Write();
  outFile->Close();
  delete outFile;
  
}

// --------------------------------------------------------------------

void goRatePlots(std::string fileType, int isCrossSec = false, int nEvents = 0) 
{

  int nBunches50ns = 1368;
  int nBunches25ns = 2508; //2508 is what agreed with TSG for # bunches

  float xSec13TeV = isCrossSec ? 78.26 : 80.; // Using McM for cross section comparison and 80 (agreed with TSG) for rates)
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
      BasicRatePlots basicRatePlots("/data2/p/pellicci/L1DPG/root/v4_62X_40PU_25bx_ReEmul2012/L1Tree.root"); 
      basicRatePlots.run(false,"13TEV_25PU_2012_RE-EMUL",0,500000000,xSec13TeV,25,nBunches25ns,isCrossSec,nEvents);
    }
  else 
    {
      std::cout << "Config param " << fileType << " invalid! \n"
		<< "Valid fileType values are : DATA, 8TEV_TF_DATA, 8TEV_TF_2012_RE-EMUL, "
		<< "8TEV_25PU_ORIG_RE-EMUL, 8TEV_25PU_2012_RE-EMUL, 8TEV_25PU_2012GCT10GEV_RE-EMUL, 8TEV_25PU_2015_RE-EMUL, "
		<< "13TEV_25PU_ORIG_RE-EMUL, 13TEV_25PU_2012_RE-EMUL, 13TEV_25PU_2012GCT10GEV_RE-EMUL, 13TEV_25PU_2015_RE-EMUL\n";
    }
    
}


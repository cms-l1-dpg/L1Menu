#include "L1Ntuple.h"
#include "L1AlgoFactory.h"
#include <algorithm>
#include<map>

Double_t convertRegionEta(int iEta) {
  static const double rgnEtaValues[11] = {
     0.174, // HB and inner HE bins are 0.348 wide
     0.522,
     0.870,
     1.218,
     1.566,
     1.956, // Last two HE bins are 0.432 and 0.828 wide
     2.586,
     3.250, // HF bins are 0.5 wide
     3.750,
     4.250,
     4.750
  };
  if(iEta < 11) {
    return -rgnEtaValues[-(iEta - 10)]; // 0-10 are negative eta values
  }
  else if (iEta < 22) {
    return rgnEtaValues[iEta - 11];     // 11-21 are positive eta values
  }
  return -9;
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

  float SingleTauPt();
  float SingleMuEta(float eta);
  float SingleEGEta(float ptCut, bool doIso);

  //GMT stuff
  float DttfPt();
  float RpcbPt();
  float RpcfPt();
  float CsctfPt();

  void setRateError(TH1F* histo);

  float computeAvgLumi(float xSec, float avPU, int nBunches) { return 11246. * avPU * nBunches / (1E7 * xSec); };

  bool PhysicsBits[128];

  L1AlgoFactory *algoFactory;

  std::map<std::string,TH1F*> hTH1F;
  std::map<std::string,TH2F*> hTH2F;

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

float BasicRatePlots::SingleEGEta(float ptCut, bool doIso) {

  float maxPt = -10;
  float iEGMaxPt = -10;

  Int_t Nele = gt_ -> Nele;
  for (Int_t ue=0; ue < Nele; ue++) {
    Int_t bx = gt_ -> Bxel[ue];        		
    if (bx != 0) continue;
    Bool_t iso = gt_ -> Isoel[ue];
    if (!iso && doIso) continue;
    Float_t pt = gt_ -> Rankel[ue];    // the rank of the electron
    if ( pt >= maxPt) 
      {
	maxPt = pt;
	iEGMaxPt = ue;
      }
  }

  float eta = -10.;
  int ieta = iEGMaxPt>=0 && maxPt>ptCut ? gt_ -> Etael[iEGMaxPt] : -10; 
  if(ieta > 0) eta = convertRegionEta(ieta);
  return eta;
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

  algoFactory = new L1AlgoFactory(gt_,gmt_);

  //Single stuff
  hTH1F["nJetVsPt"]    = new TH1F("nJetVsPt","SingleJet; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hTH1F["nTauVsPt"]    = new TH1F("nTauVsPt","SingleTau; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hTH1F["nJetCenVsPt"] = new TH1F("nJetCenVsPt","SingleJetCentral; E_{T} cut; rate [Hz]",256,-0.5,255.5);
  hTH1F["nEGVsPt"]     = new TH1F("nEGVsPt","SingleEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH1F["nIsoEGVsPt"]  = new TH1F("nIsoEGVsPt","SingleIsoEG; E_{T} cut; rate [Hz]",65,-0.5,64.5);
  hTH1F["nMuVsPt"]     = new TH1F("nMuVsPt","SingleMu; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hTH1F["nMuErVsPt"]   = new TH1F("nMuErVsPt","SingleMu |#eta|<2.1; p_{T} cut; rate [Hz]",131,-0.5,130.5);
  hTH1F["nMuVsEta"]    = new TH1F("nMuVsEta","nMuVsEta",24,-2.4,2.4);
  hTH1F["nEGVsEta"]    = new TH1F("nEGVsEta","nEGVsEta",50,-3.,3.);
  hTH1F["nIsoEGVsEta"] = new TH1F("nIsoEGVsEta","nIsoEGVsEta",50,-3.,3.);

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

  //GMT stuff
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
      
    if (event%200000 == 0) {
      std::cout << "Processed " << event << " events." << std::endl;
    }
      
    if ( event_->lumi < minLs || event_->lumi > maxLs ) continue;

    double weight = event_->puWeight > -0.001 ? event_->puWeight : 1; 
      
    FillBits();

    if(runOnData && !PhysicsBits[0]) continue;
      
    nZeroBias += weight;

    float jetPt     = 0.; algoFactory->SingleJetPt(jetPt);
    float jetCenPt  = 0.; algoFactory->SingleJetPt(jetCenPt,true);

    float tauPt     = SingleTauPt();

    float htt       = 0.; algoFactory->HTTVal(htt);
    float ett       = 0.; algoFactory->ETTVal(ett);
    float etm       = 0.; algoFactory->ETMVal(etm);

    float egPt      = 0.; algoFactory->SingleEGPt(egPt);
    float isoEgPt   = 0.; algoFactory->SingleEGPt(isoEgPt,true);
    float egEta     = SingleEGEta(16.,false);
    float isoegEta  = SingleEGEta(16.,true);

    float muPt     = -10.; algoFactory->SingleMuPt(muPt);
    float muErPt   = -10.; algoFactory->SingleMuEta2p1Pt(muErPt);
    float muEta    = SingleMuEta(16.);

    float doubleMuPt1 = -10.; 
    float doubleMuPt2 = -10.;
    algoFactory->DoubleMuPt(doubleMuPt1,doubleMuPt2);

    float oniaMuPt1 = 0.;
    float oniaMuPt2 = 0.;
    algoFactory->OniaPt(oniaMuPt1,oniaMuPt2,22);

    float EGIsoPt1 = -10;
    float EGPt2    = -10;
    algoFactory->DoubleIsoEGEGPt(EGIsoPt1,EGPt2);

    float dttfPt   = DttfPt();
    float rpcbPt   = RpcbPt();
    float rpcfPt   = RpcfPt();
    float csctfPt  = CsctfPt();

    float muPt_forHTT = -10.;
    float HTT_forMu = -10.;
    algoFactory->Mu_HTTPt(muPt_forHTT,HTT_forMu);


    float dijetPt1    = -10.;
    float dijetPt2    = -10.;
    float diCenjetPt1 = -10.;
    float diCenjetPt2 = -10.;
    algoFactory->DoubleJetPt(dijetPt1,dijetPt2);
    algoFactory->DoubleJetPt(diCenjetPt1,diCenjetPt2,true);

    Float_t dummy = -1;
    float ditauPt    = -10.; algoFactory->DoubleTauJetEta2p17Pt(dummy,ditauPt);
    dummy = -1.;
    float quadjetPt  = -10.; algoFactory->QuadJetPt(dummy,dummy,dummy,quadjetPt);
    dummy = -1.;
    float quadjetCPt = -10.; algoFactory->QuadJetPt(dummy,dummy,dummy,quadjetCPt,true);
    dummy = -1.;

    float diEG1     = -10.;
    float diEG2     = -10.;
    float diIsolEG1 = -10.;
    float diIsolEG2 = -10.;
    algoFactory->DoubleEGPt(diEG1,diEG2);
    algoFactory->DoubleEGPt(diIsolEG1,diIsolEG2,true);

    hTH1F["nPUvsPU"]->Fill(simulation_->actualInt,weight);
    hTH1F["nMuVsEta"]->Fill(muEta,weight);
    hTH1F["nEGVsEta"]->Fill(egEta,weight);
    hTH1F["nIsoEGVsEta"]->Fill(isoegEta,weight);

    for(int ptCut=0; ptCut<256; ++ptCut) {
      if(jetPt>=ptCut)	  hTH1F["nJetVsPt"]->Fill(ptCut,weight);
      if(jetCenPt>=ptCut) hTH1F["nJetCenVsPt"]->Fill(ptCut,weight);
      if(tauPt>=ptCut)	  hTH1F["nTauVsPt"]->Fill(ptCut,weight);

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

      if(ditauPt>=ptCut)    hTH1F["nDiTauVsPt"]->Fill(ptCut,weight);
      if(quadjetPt>=ptCut)  hTH1F["nQuadJetVsPt"]->Fill(ptCut,weight);
      if(quadjetCPt>=ptCut) hTH1F["nQuadCenJetVsPt"]->Fill(ptCut,weight);

    }//loop on 256
      
    for(int ptCut=0; ptCut<65; ++ptCut) {
      if(egPt>=ptCut)    hTH1F["nEGVsPt"]->Fill(ptCut,weight);
      if(isoEgPt>=ptCut) hTH1F["nIsoEGVsPt"]->Fill(ptCut,weight);

      if(diEG2>=ptCut)     hTH1F["nDiEGVsPt"]->Fill(ptCut,weight);
      if(diIsolEG2>=ptCut) hTH1F["nDiIsoEGVsPt"]->Fill(ptCut,weight);


      for(int ptCut2=0; ptCut2<=65; ++ptCut2) {
	if(diEG1>=ptCut && diEG2>=ptCut2 && ptCut2 <= ptCut) hTH2F["nEGPtVsPt"]->Fill(ptCut,ptCut2,weight);
	if(diIsolEG1>=ptCut && diIsolEG2>=ptCut2 && ptCut2<= ptCut) hTH2F["nIsoEGPtVsPt"]->Fill(ptCut,ptCut2,weight);
	if(EGIsoPt1>=ptCut && EGPt2>=ptCut2) hTH2F["nEGIsoEGVsPt"]->Fill(ptCut,ptCut2,weight);
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
	if (doubleMuPt1>=ptCut && doubleMuPt2>=ptCut2) hTH2F["nMuPtVsPt"]->Fill(ptCut,ptCut2,weight);
	if (oniaMuPt1>=ptCut && oniaMuPt2>=ptCut2)     hTH2F["nOniaMuPtVsPt"]->Fill(ptCut,ptCut2,weight);
      }
    }

    for(int iCut=0; iCut<61; ++iCut) {
      float ptCutMu = iCut*0.5;

      for(int ptCutEG=0; ptCutEG<65; ++ptCutEG){
	if(muPt>=ptCutMu && egPt>=ptCutEG) hTH2F["nMuVsEG"]->Fill(ptCutMu,ptCutEG,weight);
      }
    }

      
    for(int httCut=0; httCut<512; ++httCut) {
      if(htt>httCut) hTH1F["nHTTVsHTT"]->Fill(httCut,weight);
      if(ett>httCut) hTH1F["nETTVsETT"]->Fill(httCut,weight);
      if(etm>httCut) hTH1F["nETMVsETM"]->Fill(httCut,weight);

      for(int iCut=0; iCut<61; ++iCut) {
      float ptCutMu = iCut*0.5;
	if (muPt_forHTT>=ptCutMu && HTT_forMu>=httCut) hTH2F["nMuVsHTT"]->Fill(ptCutMu,httCut,weight);
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
      BasicRatePlots basicRatePlots("/data2/p/pellicci/L1DPG/root/v4_62X_40PU_25bx_ReEmul2015/L1Tree.root"); 
      basicRatePlots.run(false,"13TEV_40PU_2015_RE-EMUL",0,500000000,xSec13TeV,40,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "13TEV_45p4PU_2012_RE-EMUL")
    {
      BasicRatePlots basicRatePlots("/data2/battilan/L1Trigger/L1T2015Menu/L1Tree_v5_62X_13TeV_45p4PU_25bx_ReEmul2012.root"); 
      basicRatePlots.run(false,"13TEV_45p4PU_2012_RE-EMUL",0,500000000,xSec13TeV,45,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "13TEV_45p4PU_2015_RE-EMUL")
    {
      BasicRatePlots basicRatePlots("/afs/cern.ch/user/p/pellicci/data2/L1DPG/root/v4_62X_45PU_25bx_ReEmul2015/L1Tree.root");
      basicRatePlots.run(false,"13TEV_45p4PU_2015_RE-EMUL",0,500000000,xSec13TeV,45,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "13TEV_20PU_2012_RE-EMUL")
    {
      BasicRatePlots basicRatePlots("/afs/cern.ch/user/p/pellicci/data2/L1DPG/root/v4_62X_20PU_25bx_ReEmul2012/L1Tree.root"); 
      basicRatePlots.run(false,"13TEV_20PU_2012_RE-EMUL",0,500000000,xSec13TeV,20,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "13TEV_20PU_2015_RE-EMUL")
    {
      BasicRatePlots basicRatePlots("/afs/cern.ch/user/p/pellicci/data2/L1DPG/root/v4_62X_20PU_25bx_ReEmul2015/L1Tree.root"); 
      basicRatePlots.run(false,"13TEV_20PU_2015_RE-EMUL",0,500000000,xSec13TeV,20,nBunches25ns,isCrossSec,nEvents);
    }
  else if (fileType == "TEST")
    {
      BasicRatePlots basicRatePlots("/afs/cern.ch/user/p/pellicci/data2/L1DPG/root/v4_62X_20PU_25bx_ReEmul2015/L1Tree.root"); 
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


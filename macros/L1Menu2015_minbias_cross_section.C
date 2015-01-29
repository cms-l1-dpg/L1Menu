#include "L1AlgoFactory.h"
#include "L1Ntuple.h"

#include "TH2.h"
#include "TString.h"
#include "Style.C"

#include <iostream>
#include <sstream>
#include <map>
#include <set>

// huge prescale value for seeds changed on-the-fly
#define INFTY 262139

TH1F *h_Cross;
TH1F *h_MultiCross;
TH1F *h_Jets;
TH1F *h_MultiJets;
TH1F *h_Sums;
TH1F *h_Egamma;
TH1F *h_MultiEgamma;
TH1F *h_Muons;
TH1F *h_MultiMuons;
TH1F *h_Technical;

TH1F *h_Block;
TH2F *cor_Block;

const Int_t NPAGS = 7;
TH2F *cor_PAGS;
TH1F *h_PAGS_pure;
TH1F *h_PAGS_shared;

const Int_t NTRIGPHYS = 6;
TH2F *cor_TRIGPHYS;
TH1F *h_TRIGPHYS_pure;
TH1F *h_TRIGPHYS_shared;

const Int_t N128 = 128;			// could be > 128 for "test seeds"
Int_t kOFFSET = 0;
Bool_t TheTriggerBits[N128];	// contains the emulated triggers for each event
TH1F *h_All;		// one bin for each trigger. Fill bin i if event fires trigger i.
TH1F *h_Pure;		// one bin for each trigger. Fill bin i if event fires trigger i and NO OTHER TRIGGER.

// set the errors properly
void CorrectScale(TH1F* h, Float_t scal) {

  Int_t nbins = h -> GetNbinsX();

  for (Int_t i=1; i<= nbins; i++)  {
    Float_t val = h -> GetBinContent(i);
    Float_t er = sqrt(val);
    val = val * scal;
    er = er * scal;
    h -> SetBinContent(i,val);
    h -> SetBinError(i,er);
  }
}

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

class L1Menu2012 : public L1Ntuple {

 public :

  L1Menu2012(string menufile, Float_t aNumberOfBunches, Bool_t aL1JetCorrection, Bool_t anoHF, Bool_t aisTauInJet, Int_t AveragePU) : 
    themenufilename(menufile), 
    theNumberOfBunches(aNumberOfBunches),
    theL1JetCorrection(aL1JetCorrection),
    noHF(anoHF),
    noTauInJet(aisTauInJet),
    theAveragePU(AveragePU)
  {
    isOneSeedNotDefined = false;
  }

  ~L1Menu2012() {}

  // the setting below are/will be specific for each L1Ntuple file used
  string themenufilename;
  Float_t theNumberOfBunches;
  Bool_t theL1JetCorrection;
  Bool_t noHF;
  Bool_t noTauInJet;
  Int_t theAveragePU;

  std::stringstream output;
  TString GetPrintout() { return output.str(); };

  void MyInit();
  void FilL1Bits();

  L1AlgoFactory *algoFactory;

  std::map<std::string, int> Counts;
  std::map<std::string, int> Prescales;
  std::map<std::string, bool> Biased;

  std::map<std::string, int> BitMapping;

  std::map<std::string, float> WeightsPAGs;
  std::map<std::string, float> WeightsTRIGPHYS;

  void InsertInMenu(std::string L1name, Bool_t value);

  Int_t L1BitNumber(std::string l1name);

  Bool_t Cross();
  Bool_t MultiCross();
  Bool_t Jets();
  Bool_t MultiJets();
  Bool_t EGamma();
  Bool_t MultiEGamma();
  Bool_t Muons();
  Bool_t MultiMuons();
  Bool_t Sums();
  Bool_t Technical();

  void Loop();

private :

  Bool_t PhysicsBits[N128];
  Bool_t first;

  Int_t insert_ibin;
  Bool_t insert_val[100];
  std::string insert_names[100];

  Int_t NBITS_MUONS;
  Int_t NBITS_MULTIMUONS;
  Int_t NBITS_EGAMMA;
  Int_t NBITS_MULTIEGAMMA;
  Int_t NBITS_JETS;
  Int_t NBITS_MULTIJETS;
  Int_t NBITS_SUMS;
  Int_t NBITS_CROSS;
  Int_t NBITS_MULTICROSS;

  std::set<std::string> setTOP;
  std::set<std::string> setHIGGS;
  std::set<std::string> setEXO;
  std::set<std::string> setSMP;
  std::set<std::string> setBPH;
  std::set<std::string> setSUSY;
  std::set<std::string> setB2G;

  std::set<std::string> setMuon;
  std::set<std::string> setEG;
  std::set<std::string> setHadronic;
  std::set<std::string> setMuonEG;
  std::set<std::string> setMuonHadronic;
  std::set<std::string> setEGHadronic;

  Bool_t isOneSeedNotDefined;
  Int_t MissingBits[N128];

};

void L1Menu2012::InsertInMenu(std::string L1name, Bool_t value) {

  Bool_t post_prescale = false;

  Int_t prescale = 1;

  std::map<std::string, int>::const_iterator it = Prescales.find(L1name);
  if (it == Prescales.end() ) {
    isOneSeedNotDefined = true;
    MissingBits[L1BitNumber(L1name)] = 1;
    return;
  }
  else {
    prescale = Prescales[L1name];
  }

  if(prescale > 0){
    Counts[L1name] ++;
    Int_t n = Counts[L1name];
    if ( n % prescale == 0) post_prescale = value; 
  }

  insert_names[insert_ibin] = L1name;
  insert_val[insert_ibin] = post_prescale ;

  insert_ibin++;

  return;
}

Int_t L1Menu2012::L1BitNumber(std::string l1name) {

  std::map<std::string, int>::const_iterator it = BitMapping.find(l1name);
  if (it == BitMapping.end() ) {
    std::cout << " Wrong L1 name, not in BitMapping " << l1name << std::endl;
    return -1;
  }

  return BitMapping[l1name];
}

void L1Menu2012::FilL1Bits() {
  for (Int_t ibit=0; ibit < N128; ibit++) {
    PhysicsBits[ibit] = 0;
    if (ibit<64) {
      PhysicsBits[ibit] = (gt_->tw1[2]>>ibit)&1;
    }
    else {
      PhysicsBits[ibit] = (gt_->tw2[2]>>(ibit-64))&1;
    }
  }

  return;
}       

void L1Menu2012::MyInit() {

  algoFactory = new L1AlgoFactory(gt_,gmt_,gct_);
  algoFactory->setL1JetCorrection(theL1JetCorrection);
  algoFactory->setHF(noHF);
  algoFactory->setTau(noTauInJet);

  // the seeds per PAG
  setHIGGS.insert("L1_SingleMu5");
  setHIGGS.insert("L1_SingleMu12");
  setHIGGS.insert("L1_SingleMu16");
  setHIGGS.insert("L1_SingleMu20");
  setHIGGS.insert("L1_SingleMu25");
  setHIGGS.insert("L1_SingleMu30");
  setHIGGS.insert("L1_SingleMu14er");
  setHIGGS.insert("L1_SingleMu16er");
  setHIGGS.insert("L1_SingleMu18er");
  setHIGGS.insert("L1_SingleMu20er");
  setHIGGS.insert("L1_SingleMu25er");
  setHIGGS.insert("L1_SingleMu30er");
  setHIGGS.insert("L1_DoubleMu_10_3p5");
  setHIGGS.insert("L1_DoubleMu_12_5");
  setHIGGS.insert("L1_TripleMu0");
  setHIGGS.insert("L1_TripleMu_5_5_3");
  setHIGGS.insert("L1_SingleEG5");
  setHIGGS.insert("L1_SingleEG10");
  setHIGGS.insert("L1_SingleEG15");
  setHIGGS.insert("L1_SingleEG20");
  setHIGGS.insert("L1_SingleEG25");
  setHIGGS.insert("L1_SingleEG35er");
  setHIGGS.insert("L1_SingleEG40");
  setHIGGS.insert("L1_SingleIsoEG18");
  setHIGGS.insert("L1_SingleIsoEG22er");
  setHIGGS.insert("L1_SingleIsoEG25er");
  setHIGGS.insert("L1_SingleIsoEG30er");
  setHIGGS.insert("L1_DoubleEG_15_10");
  setHIGGS.insert("L1_DoubleEG_22_10");
  setHIGGS.insert("L1_DoubleEG_20_10_1LegIso");
  setHIGGS.insert("L1_TripleEG_14_10_8");
  setHIGGS.insert("L1_ETM50");
  setHIGGS.insert("L1_ETM60");
  setHIGGS.insert("L1_ETM70");
  setHIGGS.insert("L1_ETM100");
  setHIGGS.insert("L1_ETM60_NoJet52WdPhi2");
  setHIGGS.insert("L1_ETM70_NoJet52WdPhi2");
  setHIGGS.insert("L1_Mu12_EG10");
  setHIGGS.insert("L1_Mu20_EG8");
  setHIGGS.insert("L1_Mu20_EG10");
  setHIGGS.insert("L1_Mu4_EG18");
  setHIGGS.insert("L1_Mu5_EG15");
  setHIGGS.insert("L1_Mu5_EG20");
  setHIGGS.insert("L1_Mu5_IsoEG18");
  setHIGGS.insert("L1_Mu5_DoubleEG5");
  setHIGGS.insert("L1_Mu6_DoubleEG10");
  setHIGGS.insert("L1_DoubleMu6_EG6");
  setHIGGS.insert("L1_DoubleMu7_EG7");
  setHIGGS.insert("L1_Mu16er_TauJet20er");
  setHIGGS.insert("L1_IsoEG20er_TauJet20er");
  setHIGGS.insert("L1_DoubleIsoTau36er");
  setHIGGS.insert("L1_DoubleIsoTau40er");
  setHIGGS.insert("L1_DoubleIsoTau44er");
  setHIGGS.insert("L1_QuadJetC36_Tau52");
  setHIGGS.insert("L1_DoubleJetC56_ETM60");
  setHIGGS.insert("L1_DoubleJetC60_ETM60");
  setHIGGS.insert("L1_TripleJet_92_76_64");
  setHIGGS.insert("L1_EG25er_HTT125");

  setEXO.insert("L1_SingleMuOpen");
  setEXO.insert("L1_SingleMu5");
  setEXO.insert("L1_SingleMu16");
  setEXO.insert("L1_SingleMu20");
  setEXO.insert("L1_SingleMu25");
  setEXO.insert("L1_SingleMu30");
  setEXO.insert("L1_SingleMu14er");
  setEXO.insert("L1_SingleMu16er");
  setEXO.insert("L1_SingleMu20er");
  setEXO.insert("L1_SingleMu25er");
  setEXO.insert("L1_SingleMu30er");
  setEXO.insert("L1_DoubleMu_10_3p5");
  setEXO.insert("L1_DoubleMu_12_5");
  setEXO.insert("L1_SingleEG25");
  setEXO.insert("L1_SingleEG30");
  setEXO.insert("L1_SingleEG35");
  setEXO.insert("L1_SingleEG35er");
  setEXO.insert("L1_SingleIsoEG22er");
  setEXO.insert("L1_SingleIsoEG30er");
  setEXO.insert("L1_DoubleEG_15_10");
  setEXO.insert("L1_DoubleEG_22_10");
  setEXO.insert("L1_DoubleEG_20_10_1LegIso");
  setEXO.insert("L1_SingleJet52");
  setEXO.insert("L1_SingleJet128");
  setEXO.insert("L1_SingleJet176");
  setEXO.insert("L1_SingleJet200");
  setEXO.insert("L1_SingleJet240");
  setEXO.insert("L1_DoubleIsoTau36er");
  setEXO.insert("L1_DoubleIsoTau40er");
  setEXO.insert("L1_DoubleIsoTau44er");
  setEXO.insert("L1_DoubleJetC84");
  setEXO.insert("L1_DoubleJetC100");
  setEXO.insert("L1_DoubleJetC112");
  setEXO.insert("L1_DoubleJetC120");
  setEXO.insert("L1_QuadJetC40");
  setEXO.insert("L1_QuadJetC60");
  setEXO.insert("L1_QuadJetC84");
  setEXO.insert("L1_ETM30");
  setEXO.insert("L1_ETM50");
  setEXO.insert("L1_ETM60");
  setEXO.insert("L1_ETM70");
  setEXO.insert("L1_ETM100");
  setEXO.insert("L1_HTT125");
  setEXO.insert("L1_HTT150");
  setEXO.insert("L1_HTT175");
  setEXO.insert("L1_HTT200");
  setEXO.insert("L1_HTT250");
  setEXO.insert("L1_Mu20_EG8");
  setEXO.insert("L1_Mu20_EG10");
  setEXO.insert("L1_Mu4_EG18");
  setEXO.insert("L1_Mu5_EG20");
  setEXO.insert("L1_Mu5_IsoEG18");
  setEXO.insert("L1_Mu6_HTT150");
  setEXO.insert("L1_Mu14er_ETM30");
  setEXO.insert("L1_Mu10er_ETM50");
  setEXO.insert("L1_Mu3_JetC52_WdEtaPhi2");

  setSMP.insert("L1_SingleMu16");
  setSMP.insert("L1_SingleMu25");
  setSMP.insert("L1_SingleMu30");
  setSMP.insert("L1_SingleMu14er");
  setSMP.insert("L1_SingleMu16er");
  setSMP.insert("L1_SingleMu20er");
  setSMP.insert("L1_SingleMu25er");
  setSMP.insert("L1_SingleMu30er");
  setSMP.insert("L1_DoubleMu_10_3p5");
  setSMP.insert("L1_DoubleMu_12_5");
  setSMP.insert("L1_SingleEG25");
  setSMP.insert("L1_SingleEG30");
  setSMP.insert("L1_SingleEG35er");
  setSMP.insert("L1_SingleEG40");
  setSMP.insert("L1_SingleIsoEG22er");
  setSMP.insert("L1_SingleIsoEG25er");
  setSMP.insert("L1_SingleIsoEG30er");
  setSMP.insert("L1_DoubleEG_15_10");
  setSMP.insert("L1_DoubleEG_22_10");
  setSMP.insert("L1_DoubleEG_20_10_1LegIso");
  setSMP.insert("L1_SingleJet52");
  setSMP.insert("L1_SingleJet68");
  setSMP.insert("L1_SingleJet92");
  setSMP.insert("L1_SingleJet128");
  setSMP.insert("L1_SingleJet176");
  setSMP.insert("L1_SingleJet200");
  setSMP.insert("L1_SingleJet240");
  setSMP.insert("L1_Mu20_EG8");
  setSMP.insert("L1_Mu20_EG10");
  setSMP.insert("L1_Mu4_EG18");
  setSMP.insert("L1_Mu5_EG20");
  setSMP.insert("L1_Mu5_IsoEG18");

  setBPH.insert("L1_SingleMu16");
  setBPH.insert("L1_SingleMu25");
  setBPH.insert("L1_SingleMu30");
  setBPH.insert("L1_DoubleMu0");
  setBPH.insert("L1_TripleMu0");
  setBPH.insert("L1_DoubleMu0er16_WdEta18_OS");
  setBPH.insert("L1_DoubleMu0er16_WdEta18");
  setBPH.insert("L1_DoubleMu_10_0_WdEta18");
  setBPH.insert("L1_QuadMu0");

  setTOP.insert("L1_SingleMu16");
  setTOP.insert("L1_SingleMu25");
  setTOP.insert("L1_SingleMu30");
  setTOP.insert("L1_DoubleMu_10_3p5");
  setTOP.insert("L1_DoubleMu_12_5");
  setTOP.insert("L1_SingleIsoEG20er");
  setTOP.insert("L1_SingleIsoEG22er");
  setTOP.insert("L1_SingleIsoEG30er");
  setTOP.insert("L1_DoubleEG_15_10");
  setTOP.insert("L1_DoubleEG_22_10");
  setTOP.insert("L1_DoubleEG_20_10_1LegIso");
  setTOP.insert("L1_Mu20_EG8");
  setTOP.insert("L1_Mu20_EG10");
  setTOP.insert("L1_Mu4_EG18");
  setTOP.insert("L1_Mu5_EG20");
  setTOP.insert("L1_Mu5_IsoEG18");

  setB2G.insert("L1_SingleMu16");
  setB2G.insert("L1_SingleMu25");
  setB2G.insert("L1_SingleMu30");
  setB2G.insert("L1_SingleMu14er");
  setB2G.insert("L1_SingleMu20er");
  setB2G.insert("L1_SingleMu25er");
  setB2G.insert("L1_SingleMu30er");
  setB2G.insert("L1_DoubleMu_10_3p5");
  setB2G.insert("L1_DoubleMu_12_5");
  setB2G.insert("L1_SingleEG25");
  setB2G.insert("L1_SingleEG35er");
  setB2G.insert("L1_SingleIsoEG22er");
  setB2G.insert("L1_SingleIsoEG30er");
  setB2G.insert("L1_DoubleEG_15_10");
  setB2G.insert("L1_DoubleEG_22_10");
  setB2G.insert("L1_SingleJet128");
  setB2G.insert("L1_SingleJet176");
  setB2G.insert("L1_SingleJet200");
  setB2G.insert("L1_SingleJet240");
  setB2G.insert("L1_QuadJetC40");
  setB2G.insert("L1_QuadJetC60");
  setB2G.insert("L1_QuadJetC84");
  setB2G.insert("L1_ETM50");
  setB2G.insert("L1_ETM70");
  setB2G.insert("L1_ETM100");
  setB2G.insert("L1_HTT125");
  setB2G.insert("L1_HTT150");
  setB2G.insert("L1_HTT175");
  setB2G.insert("L1_HTT200");
  setB2G.insert("L1_HTT250");
  setB2G.insert("L1_Mu4_EG18");
  setB2G.insert("L1_Mu5_EG20");
  setB2G.insert("L1_Mu5_IsoEG18");
  setB2G.insert("L1_Mu20_EG8");
  setB2G.insert("L1_Mu20_EG10");

  setSUSY.insert("L1_SingleMu5");
  setSUSY.insert("L1_SingleMu12");
  setSUSY.insert("L1_SingleMu16");
  setSUSY.insert("L1_SingleMu14er");
  setSUSY.insert("L1_SingleMu20er");
  setSUSY.insert("L1_SingleMu25er");
  setSUSY.insert("L1_SingleMu30er");
  setSUSY.insert("L1_DoubleMu_10_3p5");
  setSUSY.insert("L1_DoubleMu_12_5");
  setSUSY.insert("L1_SingleIsoEG22er");
  setSUSY.insert("L1_SingleIsoEG30er");
  setSUSY.insert("L1_DoubleEG_15_10");
  setSUSY.insert("L1_DoubleEG_22_10");
  setSUSY.insert("L1_DoubleEG_20_10_1LegIso");
  setSUSY.insert("L1_SingleJet128");
  setSUSY.insert("L1_SingleJet176");
  setSUSY.insert("L1_SingleJet200");
  setSUSY.insert("L1_SingleJet240");
  setSUSY.insert("L1_DoubleIsoTau36er");
  setSUSY.insert("L1_DoubleIsoTau40er");
  setSUSY.insert("L1_DoubleIsoTau44er");
  setSUSY.insert("L1_DoubleJetC84");
  setSUSY.insert("L1_DoubleJetC100");
  setSUSY.insert("L1_DoubleJetC112");
  setSUSY.insert("L1_DoubleJetC120");
  setSUSY.insert("L1_QuadJetC40");
  setSUSY.insert("L1_QuadJetC60");
  setSUSY.insert("L1_QuadJetC84");
  setSUSY.insert("L1_ETM50");
  setSUSY.insert("L1_ETM70");
  setSUSY.insert("L1_ETM100");
  setSUSY.insert("L1_HTT125");
  setSUSY.insert("L1_HTT150");
  setSUSY.insert("L1_HTT175");
  setSUSY.insert("L1_HTT200");
  setSUSY.insert("L1_HTT250");
  setSUSY.insert("L1_Mu20_EG8");
  setSUSY.insert("L1_Mu20_EG10");
  setSUSY.insert("L1_Mu4_EG18");
  setSUSY.insert("L1_Mu5_EG20");
  setSUSY.insert("L1_Mu5_IsoEG18");
  setSUSY.insert("L1_DoubleJetC56_ETM60");
  setSUSY.insert("L1_DoubleJetC60_ETM60");
  setSUSY.insert("L1_DoubJetC32_WdPhi7_HTT125");
  setSUSY.insert("L1_Mu0er_ETM55");
  setSUSY.insert("L1_Mu8_HTT125");
  setSUSY.insert("L1_DoubleEG6_HTT150");
  setSUSY.insert("L1_Jet32MuOpen_Mu10_dPhiMu_Mu1");
  setSUSY.insert("L1_Jet32MuOpen_EG10_dPhiMu_EG1");

  // the seeds per physics triggers (TRIGPHYS));
  setMuon.insert("L1_SingleMuOpen");
  setMuon.insert("L1_SingleMu5");
  setMuon.insert("L1_SingleMu12");
  setMuon.insert("L1_SingleMu16");
  setMuon.insert("L1_SingleMu20");
  setMuon.insert("L1_SingleMu25");
  setMuon.insert("L1_SingleMu30");
  setMuon.insert("L1_SingleMu14er");
  setMuon.insert("L1_SingleMu16er");
  setMuon.insert("L1_SingleMu18er");
  setMuon.insert("L1_SingleMu20er");
  setMuon.insert("L1_SingleMu25er");
  setMuon.insert("L1_SingleMu30er");
  setMuon.insert("L1_DoubleMu0");
  setMuon.insert("L1_DoubleMu0er16_WdEta18_OS");
  setMuon.insert("L1_DoubleMu0er16_WdEta18");
  setMuon.insert("L1_DoubleMu_10_0_WdEta18");
  setMuon.insert("L1_DoubleMu_10_Open");
  setMuon.insert("L1_DoubleMu_10_3p5");
  setMuon.insert("L1_DoubleMu_12_5");
  setMuon.insert("L1_TripleMu0");
  setMuon.insert("L1_TripleMu_5_5_3");
  setMuon.insert("L1_SingleMu6_NotBptxOR");
  setMuon.insert("L1_QuadMu0");
  setMuonHadronic.insert("L1_Mu3_JetC16_WdEtaPhi2");
  setMuonHadronic.insert("L1_Mu3_JetC52_WdEtaPhi2");
  setMuonHadronic.insert("L1_Mu6_HTT150");
  setMuonHadronic.insert("L1_Mu8_HTT125");
  setMuonHadronic.insert("L1_Mu0er_ETM55");
  setMuonHadronic.insert("L1_Mu14er_ETM30");
  setMuonHadronic.insert("L1_Mu10er_ETM50");
  setMuonHadronic.insert("L1_Mu16er_TauJet20er");
  setMuonEG.insert("L1_Mu4_EG18");
  setMuonEG.insert("L1_Mu5_EG15");
  setMuonEG.insert("L1_Mu5_EG20");
  setMuonEG.insert("L1_Mu5_IsoEG18");
  setMuonEG.insert("L1_Mu12_EG10");
  setMuonEG.insert("L1_Mu20_EG8");
  setMuonEG.insert("L1_Mu20_EG10");
  setMuonEG.insert("L1_Mu5_DoubleEG5");
  setMuonEG.insert("L1_Mu6_DoubleEG10");
  setMuonEG.insert("L1_DoubleMu6_EG6");
  setMuonEG.insert("L1_DoubleMu7_EG7");
  setEG.insert("L1_SingleEG5");
  setEG.insert("L1_SingleEG10");
  setEG.insert("L1_SingleEG15");
  setEG.insert("L1_SingleEG20");
  setEG.insert("L1_SingleEG25");
  setEG.insert("L1_SingleEG30");
  setEG.insert("L1_SingleEG35");
  setEG.insert("L1_SingleEG35er");
  setEG.insert("L1_SingleEG40");
  setEG.insert("L1_SingleIsoEG18");
  setEG.insert("L1_SingleIsoEG25");
  setEG.insert("L1_SingleIsoEG20er");
  setEG.insert("L1_SingleIsoEG22er");
  setEG.insert("L1_SingleIsoEG25er");
  setEG.insert("L1_SingleIsoEG28er");
  setEG.insert("L1_SingleIsoEG30er");
  setEG.insert("L1_DoubleEG_15_10");
  setEG.insert("L1_DoubleEG_22_10");
  setEG.insert("L1_DoubleEG_20_10_1LegIso");
  setEG.insert("L1_TripleEG_14_10_8");
  setHadronic.insert("L1_SingleJet36");
  setHadronic.insert("L1_SingleJet52");
  setHadronic.insert("L1_SingleJet68");
  setHadronic.insert("L1_SingleJet92");
  setHadronic.insert("L1_SingleJet128");
  setHadronic.insert("L1_SingleJet176"); 
  setHadronic.insert("L1_SingleJet200"); 
  setHadronic.insert("L1_SingleJet240");
  setHadronic.insert("L1_DoubleJetC52");
  setHadronic.insert("L1_DoubleJetC72");
  setHadronic.insert("L1_DoubleJetC84");
  setHadronic.insert("L1_DoubleJetC100");
  setHadronic.insert("L1_DoubleJetC112");
  setHadronic.insert("L1_DoubleJetC120");
  setHadronic.insert("L1_DoubleIsoTau36er");
  setHadronic.insert("L1_DoubleIsoTau40er");
  setHadronic.insert("L1_DoubleIsoTau44er");
  setHadronic.insert("L1_TripleJet_92_76_64");
  setHadronic.insert("L1_QuadJetC40");
  setHadronic.insert("L1_QuadJetC60");
  setHadronic.insert("L1_QuadJetC84");
  setHadronic.insert("L1_HTT125");
  setHadronic.insert("L1_HTT150");
  setHadronic.insert("L1_HTT175");
  setHadronic.insert("L1_HTT200");
  setHadronic.insert("L1_HTT250");
  setHadronic.insert("L1_ETM30");
  setHadronic.insert("L1_ETM40");
  setHadronic.insert("L1_ETM50");
  setHadronic.insert("L1_ETM60");
  setHadronic.insert("L1_ETM70");
  setHadronic.insert("L1_ETM100");
  setHadronic.insert("L1_ETM60_NoJet52WdPhi2");
  setHadronic.insert("L1_ETM70_NoJet52WdPhi2");
  setHadronic.insert("L1_QuadJetC36_Tau52");
  setHadronic.insert("L1_DoubleJetC56_ETM60");
  setHadronic.insert("L1_DoubleJetC60_ETM60");
  setHadronic.insert("L1_DoubJetC32_WdPhi7_HTT125");
  setEGHadronic.insert("L1_IsoEG20er_TauJet20er");
  setEGHadronic.insert("L1_DoubleEG6_HTT150");
  setEGHadronic.insert("L1_EG25er_HTT125");
  setEGHadronic.insert("L1_Jet32MuOpen_Mu10_dPhiMu_Mu1");
  setEGHadronic.insert("L1_Jet32MuOpen_EG10_dPhiMu_EG1");

  BitMapping["L1_ZeroBias"] = 0;
  BitMapping["L1_AlwaysTrue"] = 1;
  BitMapping["FREE2"] = 2;
  BitMapping["FREE3"] = 3;
  BitMapping["L1_DoubleEG6_HTT150"] = 4;
  BitMapping["L1_Jet32MuOpen_Mu10_dPhiMu_Mu1"] = 5;
  BitMapping["L1_Jet32MuOpen_EG10_dPhiMu_EG1"] = 6;
  BitMapping["L1_EG25er_HTT125"] = 7;
  BitMapping["FREE8"] = 8;
  BitMapping["FREE9"] = 9;
  BitMapping["FREE10"] = 10;
  BitMapping["L1_DoubleJetC72"] = 11;
  BitMapping["L1_IsoEG20er_TauJet20er"] = 12;
  BitMapping["L1_Mu16er_TauJet20er"] = 13;
  BitMapping["L1_QuadJetC84"] = 14;
  BitMapping["L1_SingleJet36"] = 15;
  BitMapping["L1_DoubleJetC120"] = 16;
  BitMapping["L1_SingleJet52"] = 17;
  BitMapping["L1_SingleJet68"] = 18;
  BitMapping["L1_SingleJet92"] = 19;
  BitMapping["L1_SingleJet128"] = 20;
  BitMapping["L1_SingleJet176"] = 21;
  BitMapping["L1_SingleJet200"] = 22;
  BitMapping["L1_DoubleIsoTau36er"] = 23;
  BitMapping["L1_DoubleIsoTau40er"] = 24;
  BitMapping["L1_DoubleIsoTau44er"] = 25;
  BitMapping["L1_DoubleMu0"] = 26;
  BitMapping["L1_SingleMu30er"] = 27;
  BitMapping["L1_Mu3_JetC16_WdEtaPhi2"] = 28;
  BitMapping["L1_Mu3_JetC52_WdEtaPhi2"] = 29;
  BitMapping["L1_SingleEG15"] = 30;
  BitMapping["L1_SingleEG35er"] = 31;
  BitMapping["L1_SingleEG40"] = 32;
  BitMapping["L1_SingleIsoEG25er"] = 33;
  BitMapping["L1_SingleIsoEG25"] = 34;
  BitMapping["L1_SingleIsoEG28er"] = 35;
  BitMapping["L1_SingleIsoEG30er"] = 36;
  BitMapping["L1_SingleEG10"] = 37;
  BitMapping["L1_SingleIsoEG20er"] = 38;
  BitMapping["L1_SingleJet240"] = 39;
  BitMapping["L1_DoubleJetC52"] = 40;
  BitMapping["L1_DoubleJetC84"] = 41;
  BitMapping["L1_SingleMu14er"] = 42 ;
  BitMapping["L1_DoubleJetC112"] = 43;
  BitMapping["L1_DoubleMu_10_Open"] = 44 ;
  BitMapping["L1_DoubleMu_10_3p5"] = 45 ;
  BitMapping["L1_QuadJetC40"] = 46;
  BitMapping["L1_SingleEG5"] = 47;
  BitMapping["L1_SingleEG25"] = 48;
  BitMapping["L1_SingleIsoEG22er"] = 49;
  BitMapping["L1_SingleIsoEG18"] = 50;
  BitMapping["L1_Mu8_HTT125"] = 51;
  BitMapping["L1_SingleEG20"] = 52;
  BitMapping["L1_SingleEG30"] = 53;
  BitMapping["L1_SingleEG35"] = 54;
  BitMapping["L1_SingleMuOpen"] = 55;
  BitMapping["L1_SingleMu16"] = 56;
  BitMapping["FREE57"] = 57;
  BitMapping["L1_SingleMu5"] = 58;
  BitMapping["FREE59"] = 59;
  BitMapping["L1_SingleMu20er"] = 60;
  BitMapping["L1_SingleMu12"] = 61;
  BitMapping["L1_SingleMu20"] = 62;
  BitMapping["L1_SingleMu25er"] = 63;
  BitMapping["L1_SingleMu25"] = 64;
  BitMapping["L1_SingleMu30"] = 65;
  BitMapping["L1_ETM30"] = 66;
  BitMapping["L1_ETM50"] = 67;
  BitMapping["L1_ETM70"] = 68;
  BitMapping["L1_ETM100"] = 69;
  BitMapping["L1_HTT125"] = 70;
  BitMapping["L1_HTT150"] = 71;
  BitMapping["L1_HTT175"] = 72;
  BitMapping["L1_HTT200"] = 73;
  BitMapping["L1_Mu20_EG10"] = 74;
  BitMapping["L1_Mu5_EG20"] = 75;
  BitMapping["L1_Mu5_IsoEG18"] = 76;
  BitMapping["L1_Mu6_DoubleEG10"] = 77;
  BitMapping["L1_SingleJetC32_NotBptxOr"] = 78;
  BitMapping["L1_ETM40"] = 79;
  BitMapping["L1_HTT250"] = 80;
  BitMapping["L1_Mu20_EG8"] = 81;
  BitMapping["L1_Mu6_HTT150"] = 82;
  BitMapping["L1_Mu10er_ETM50"] = 83;
  BitMapping["L1_Mu14er_ETM30"] = 84;
  BitMapping["L1_DoubleMu7_EG7"] = 85;
  BitMapping["L1_SingleMu16er"] = 86;
  BitMapping["L1_Mu5_EG15"] = 87;
  BitMapping["L1_Mu4_EG18"] = 88;
  BitMapping["L1_SingleMu6_NotBptxOR"] = 89;
  BitMapping["L1_DoubleMu6_EG6"] = 90;
  BitMapping["L1_Mu5_DoubleEG5"] = 91;
  BitMapping["L1_Mu12_EG10"] = 92;
  BitMapping["L1_QuadJetC36_Tau52"] = 93;
  BitMapping["FREE94"] = 94;
  BitMapping["FREE95"] = 95;
  BitMapping["FREE96"] = 96;
  BitMapping["L1_TripleMu0"] = 97;
  BitMapping["L1_TripleMu_5_5_3"] = 98;
  BitMapping["FREE99"] = 99;
  BitMapping["L1_TripleEG_14_10_8"] = 100;
  BitMapping["L1_DoubleEG_15_10"] = 101;
  BitMapping["L1_DoubleEG_22_10"] = 102;
  BitMapping["L1_DoubleEG_20_10_1LegIso"] = 103;
  BitMapping["L1_Mu0er_ETM55"] = 104;
  BitMapping["L1_DoubleJetC60_ETM60"] = 105;
  BitMapping["L1_DoubJetC32_WdPhi7_HTT125"] = 106;
  BitMapping["FREE107"] = 107;
  BitMapping["FREE108"] = 108;
  BitMapping["L1_ETM60"] = 109;
  BitMapping["L1_DoubleJetC100"] = 110;
  BitMapping["L1_QuadJetC60"] = 111;
  BitMapping["L1_SingleJetC20_NotBptxOr"] = 112;
  BitMapping["L1_DoubleJetC56_ETM60"] = 113;
  BitMapping["L1_DoubleMu0er16_WdEta18"] = 114;
  BitMapping["L1_ETM60_NoJet52WdPhi2"] = 115;
  BitMapping["L1_ETM70_NoJet52WdPhi2"] = 116;
  BitMapping["FREE117"] = 117;
  BitMapping["FREE118"] = 118;
  BitMapping["L1_TripleJet_92_76_64"] = 119;
  BitMapping["FREE120"] = 120;
  BitMapping["L1_QuadMu0"] = 121;
  BitMapping["L1_SingleMu18er"] = 122;
  BitMapping["L1_DoubleMu0er16_WdEta18_OS"] = 123;
  BitMapping["L1_DoubleMu_12_5"] = 124;
  BitMapping["FREE125"] = 125;
  BitMapping["L1_DoubleMu_10_0_WdEta18"] = 126;
  BitMapping["L1_SingleMuBeamHalo"] = 127;

  //Read the prescales table
  ifstream menufile;
  menufile.open(themenufilename);
  char path_name[50];
  Int_t prescale_value = -1;

  while(menufile >> path_name >> prescale_value){
    if(prescale_value < 0) Prescales[path_name] = INFTY;
    else Prescales[path_name] = prescale_value;
  }

  for (std::map<std::string, int>::iterator it=Prescales.begin(); it != Prescales.end(); it++) {
    std::string name = it -> first;

    // each seed gets a "weight" according to how many PAGS are using it
    Int_t UsedPernPAG = 0;
    if(setTOP.count(name)   > 0) UsedPernPAG ++;
    if(setHIGGS.count(name) > 0) UsedPernPAG ++;
    if(setSUSY.count(name)  > 0) UsedPernPAG ++;
    if(setEXO.count(name)   > 0) UsedPernPAG ++;
    if(setSMP.count(name)   > 0) UsedPernPAG ++;
    if(setBPH.count(name)   > 0) UsedPernPAG ++;
    if(setB2G.count(name)   > 0) UsedPernPAG ++;
    WeightsPAGs[name] = 1./(float)UsedPernPAG;

    // each seed gets a "weight" according to how many trigger groups are using it
    Int_t UsedPernTrigPhyGroup = 0;
    if(setMuon.count(name) > 0)         UsedPernTrigPhyGroup ++;
    if(setEG.count(name) > 0)           UsedPernTrigPhyGroup ++;
    if(setHadronic.count(name) > 0)     UsedPernTrigPhyGroup ++;
    if(setMuonEG.count(name) > 0)       UsedPernTrigPhyGroup ++;
    if(setMuonHadronic.count(name) > 0) UsedPernTrigPhyGroup ++;
    if(setEGHadronic.count(name) > 0)   UsedPernTrigPhyGroup ++;
    WeightsTRIGPHYS[name] = 1./(float)UsedPernTrigPhyGroup;
  }

  // The "Biased" table is only used for the final print-out set true for seeds for which the rate estimation is biased by the sample (because of the seeds enabled in the high PU run)
  for (std::map<std::string, int>::iterator it=Prescales.begin(); it != Prescales.end(); it++) {
    std::string name = it -> first;
    Counts[name] = 0;
    Biased[name] = false; 
  }


}

Bool_t L1Menu2012::Muons() {

  insert_ibin = 0;

  Float_t SingleMuPt = -10.;
  algoFactory->SingleMuPt(SingleMuPt);

  InsertInMenu("L1_SingleMuOpen",algoFactory->SingleMu(0.,0));
  InsertInMenu("L1_SingleMu5",SingleMuPt >= 5.);
  InsertInMenu("L1_SingleMu12",SingleMuPt >= 12.);
  InsertInMenu("L1_SingleMu16",SingleMuPt >= 16.);
  InsertInMenu("L1_SingleMu20",SingleMuPt >= 20.);
  InsertInMenu("L1_SingleMu25",SingleMuPt >= 25.);
  InsertInMenu("L1_SingleMu30",SingleMuPt >= 30.);

  Float_t SingleMuEta2p1Pt = -10.;
  algoFactory->SingleMuEta2p1Pt(SingleMuEta2p1Pt);

  InsertInMenu("L1_SingleMu14er",SingleMuEta2p1Pt >= 14.);
  InsertInMenu("L1_SingleMu16er",SingleMuEta2p1Pt >= 16.);
  InsertInMenu("L1_SingleMu18er",SingleMuEta2p1Pt >= 18.);
  InsertInMenu("L1_SingleMu20er",SingleMuEta2p1Pt >= 20.);
  InsertInMenu("L1_SingleMu25er",SingleMuEta2p1Pt >= 25.);
  InsertInMenu("L1_SingleMu30er",SingleMuEta2p1Pt >= 30.);

  Int_t NN = insert_ibin;

  Int_t kOFFSET_old = kOFFSET;
  for (Int_t k=0; k < NN; k++) {
    TheTriggerBits[k + kOFFSET_old] = insert_val[k];
  }
  kOFFSET += insert_ibin;

  if(first){

    NBITS_MUONS = NN;

    for (Int_t ibin=0; ibin < insert_ibin; ibin++) {
      TString l1name = (TString)insert_names[ibin];
      h_Muons -> GetXaxis() -> SetBinLabel(ibin+1, l1name );
    }
    h_Muons -> GetXaxis() -> SetBinLabel(NN+1, "MUONS") ;

    for (Int_t k=1; k <= kOFFSET - kOFFSET_old ; k++) {
      h_All -> GetXaxis() -> SetBinLabel(k +kOFFSET_old , h_Muons -> GetXaxis() -> GetBinLabel(k) );
    }
  }

  Bool_t res = false;
  for (Int_t i=0; i < NN; i++) {
    res = res || insert_val[i] ;
    if (insert_val[i]) h_Muons -> Fill(i);
  }
  if (res) h_Muons -> Fill(NN);

  return res;
}

Bool_t L1Menu2012::MultiMuons() {

  insert_ibin = 0;
  InsertInMenu("L1_DoubleMu0er16_WdEta18_OS",algoFactory->Onia2015(0.,0.,true,true,18));
  InsertInMenu("L1_DoubleMu0er16_WdEta18",algoFactory->Onia2015(0.,0.,true,false,18));
  InsertInMenu("L1_DoubleMu_10_0_WdEta18",algoFactory->Onia2015(10.,0.,false,false,18));

  InsertInMenu("L1_DoubleMu0",algoFactory->DoubleMu(0.,0.,true));
  InsertInMenu("L1_DoubleMu_10_Open",algoFactory->DoubleMuXOpen(10.));
  InsertInMenu("L1_DoubleMu_10_3p5",algoFactory->DoubleMu(10.,3.5,true));
  InsertInMenu("L1_DoubleMu_12_5",algoFactory->DoubleMu(12.,5.,true));

  InsertInMenu("L1_TripleMu0",algoFactory->TripleMu(0.,0.,0.,4));
  InsertInMenu("L1_TripleMu_5_5_3",algoFactory->TripleMu(5.,5.,3.,4));

  InsertInMenu("L1_QuadMu0",algoFactory->QuadMu(0.,0.,0.,0.,4));

  Int_t NN = insert_ibin;

  Int_t kOFFSET_old = kOFFSET;
  for (Int_t k=0; k < NN; k++) {
    TheTriggerBits[k + kOFFSET_old] = insert_val[k];
  }
  kOFFSET += insert_ibin;

  if(first){

    NBITS_MULTIMUONS = NN;

    for (Int_t ibin=0; ibin < insert_ibin; ibin++) {
      TString l1name = (TString)insert_names[ibin];
      h_MultiMuons -> GetXaxis() -> SetBinLabel(ibin+1, l1name );
    }
    h_MultiMuons -> GetXaxis() -> SetBinLabel(NN+1, "MULTIMUONS") ;

    for (Int_t k=1; k <= kOFFSET - kOFFSET_old ; k++) {
      h_All -> GetXaxis() -> SetBinLabel(k +kOFFSET_old , h_MultiMuons -> GetXaxis() -> GetBinLabel(k) );
    }
  }

  Bool_t res = false;
  for (Int_t i=0; i < NN; i++) {
    res = res || insert_val[i] ;
    if (insert_val[i]) h_MultiMuons -> Fill(i);
  }
  if (res) h_MultiMuons -> Fill(NN);

  return res;
}

Bool_t L1Menu2012::Cross() {

  insert_ibin = 0;

  InsertInMenu("L1_Mu6_HTT150", algoFactory->Mu_HTT(6.,150.));
  InsertInMenu("L1_Mu8_HTT125", algoFactory->Mu_HTT(8.,125.));
  InsertInMenu("L1_Mu0er_ETM55", algoFactory->Muer_ETM(0.,55.));
  InsertInMenu("L1_Mu14er_ETM30", algoFactory->Muer_ETM(14.,30.));
  InsertInMenu("L1_Mu10er_ETM50", algoFactory->Muer_ETM(10.,50.));
  InsertInMenu("L1_EG25er_HTT125", algoFactory->SingleEG_Eta2p1_HTT(25., 125.,false));
  InsertInMenu("L1_Mu16er_TauJet20er", algoFactory->Muer_TauJetEta2p17(16.,20.));
  InsertInMenu("L1_IsoEG20er_TauJet20er", algoFactory->IsoEGer_TauJetEta2p17(20.,20.));
  InsertInMenu("L1_Mu12_EG10", algoFactory->Mu_EG(12.,10.));
  InsertInMenu("L1_Mu20_EG8", algoFactory->Mu_EG(20.,8.));
  InsertInMenu("L1_Mu20_EG10", algoFactory->Mu_EG(20.,10.));
  InsertInMenu("L1_Mu4_EG18", algoFactory->Mu_EG(4.,18.));
  InsertInMenu("L1_Mu5_EG15", algoFactory->Mu_EG(5.,15.));
  InsertInMenu("L1_Mu5_EG20", algoFactory->Mu_EG(5.,20.));
  InsertInMenu("L1_Mu5_IsoEG18", algoFactory->Mu_EG(5.,18.,true));
  InsertInMenu("L1_DoubleMu6_EG6", algoFactory->DoubleMu_EG(6.,6.,true));
  InsertInMenu("L1_DoubleMu7_EG7", algoFactory->DoubleMu_EG(7,7.,true));
  InsertInMenu("L1_Mu5_DoubleEG5", algoFactory->Mu_DoubleEG(5., 5.));
  InsertInMenu("L1_Mu6_DoubleEG10", algoFactory->Mu_DoubleEG(6., 10.));

  Int_t NN = insert_ibin;
  Int_t kOFFSET_old = kOFFSET;
  for (Int_t k=0; k < NN; k++) {
    TheTriggerBits[k + kOFFSET_old] = insert_val[k];
  }
  kOFFSET += NN;

  if (first) {

    NBITS_CROSS = NN;

    for (Int_t ibin=0; ibin < insert_ibin; ibin++) {
      TString l1name = (TString)insert_names[ibin];
      h_Cross -> GetXaxis() -> SetBinLabel(ibin+1, l1name );
    }
    h_Cross-> GetXaxis() -> SetBinLabel(NN+1,"CROSS");

    for (Int_t k=1; k <= kOFFSET - kOFFSET_old; k++) {
      h_All -> GetXaxis() -> SetBinLabel(k +kOFFSET_old , h_Cross -> GetXaxis() -> GetBinLabel(k) );
    }

  }

  Bool_t res = false;
  for (Int_t i=0; i < NN; i++) {
    res = res || insert_val[i] ;
    if (insert_val[i]) h_Cross -> Fill(i);
  }
  if (res) h_Cross -> Fill(NN);

  return res;
}

Bool_t L1Menu2012::MultiCross() {

  insert_ibin = 0;

  InsertInMenu("L1_Mu3_JetC16_WdEtaPhi2", algoFactory->Mu_JetCentral_delta(3.,16.));
  InsertInMenu("L1_Mu3_JetC52_WdEtaPhi2", algoFactory->Mu_JetCentral_delta(3.,52.));

  InsertInMenu("L1_DoubleJetC56_ETM60", algoFactory->DoubleJetCentral_ETM(56.,56.,60.));
  InsertInMenu("L1_DoubleJetC60_ETM60", algoFactory->DoubleJetCentral_ETM(60.,60.,60.));
  InsertInMenu("L1_DoubJetC32_WdPhi7_HTT125", algoFactory->DoubleJetC_deltaPhi7_HTT(32.,125.));

  InsertInMenu("L1_DoubleEG6_HTT150", algoFactory->DoubleEG_HT(6., 150.));

  InsertInMenu("L1_Jet32MuOpen_Mu10_dPhiMu_Mu1", algoFactory->Jet_MuOpen_Mu_dPhiMuMu1(32.,10.));
  InsertInMenu("L1_Jet32MuOpen_EG10_dPhiMu_EG1", algoFactory->Jet_MuOpen_EG_dPhiMuEG1(32.,10.));

  Int_t NN = insert_ibin;
  Int_t kOFFSET_old = kOFFSET;
  for (Int_t k=0; k < NN; k++) {
    TheTriggerBits[k + kOFFSET_old] = insert_val[k];
  }
  kOFFSET += NN;

  if (first) {

    NBITS_MULTICROSS = NN;

    for (Int_t ibin=0; ibin < insert_ibin; ibin++) {
      TString l1name = (TString)insert_names[ibin];
      h_MultiCross -> GetXaxis() -> SetBinLabel(ibin+1, l1name );
    }
    h_MultiCross-> GetXaxis() -> SetBinLabel(NN+1,"MULTICROSS");

    for (Int_t k=1; k <= kOFFSET - kOFFSET_old; k++) {
      h_All -> GetXaxis() -> SetBinLabel(k +kOFFSET_old , h_MultiCross -> GetXaxis() -> GetBinLabel(k) );
    }

  }

  Bool_t res = false;
  for (Int_t i=0; i < NN; i++) {
    res = res || insert_val[i] ;
    if (insert_val[i]) h_MultiCross -> Fill(i);
  }
  if (res) h_MultiCross -> Fill(NN);

  return res;
}

Bool_t L1Menu2012::Jets() {

  insert_ibin = 0;

  Float_t SingleJetPt = -10.;
  algoFactory->SingleJetPt(SingleJetPt);

  InsertInMenu("L1_SingleJet36",SingleJetPt >= 36.);
  InsertInMenu("L1_SingleJet52",SingleJetPt >= 52.);
  InsertInMenu("L1_SingleJet68",SingleJetPt >= 68.);
  InsertInMenu("L1_SingleJet92",SingleJetPt >= 92.);
  InsertInMenu("L1_SingleJet128",SingleJetPt >= 128.);
  InsertInMenu("L1_SingleJet176",SingleJetPt >= 176.);
  InsertInMenu("L1_SingleJet200",SingleJetPt >= 200.);
  InsertInMenu("L1_SingleJet240",SingleJetPt >= 240.);

  Int_t NN = insert_ibin;

  Int_t kOFFSET_old = kOFFSET;
  for (Int_t k=0; k < NN; k++) {
    TheTriggerBits[k + kOFFSET_old] = insert_val[k];
  }
  kOFFSET += insert_ibin;


  if (first) {

    NBITS_JETS = NN;

    for (Int_t ibin=0; ibin < insert_ibin; ibin++) {
      TString l1name = (TString)insert_names[ibin];
      h_Jets -> GetXaxis() -> SetBinLabel(ibin+1, l1name );
    }

    h_Jets-> GetXaxis() -> SetBinLabel(NN+1,"JETS");

    for (Int_t k=1; k <= kOFFSET -kOFFSET_old; k++) {
      h_All -> GetXaxis() -> SetBinLabel(k +kOFFSET_old , h_Jets -> GetXaxis() -> GetBinLabel(k) );
    }

  }

  Bool_t res = false;
  for (Int_t i=0; i < NN; i++) {
    res = res || insert_val[i] ;
    if (insert_val[i]) h_Jets -> Fill(i);
  }
  if (res) h_Jets -> Fill(NN);

  return res;
}

Bool_t L1Menu2012::MultiJets() {

  insert_ibin = 0;

  Float_t DoubleJet1 = -10.;
  Float_t DoubleJet2 = -10.;
  algoFactory->DoubleJetPt(DoubleJet1,DoubleJet2,true);

  InsertInMenu("L1_DoubleJetC52",DoubleJet1 >= 52. && DoubleJet2 >= 52.);
  InsertInMenu("L1_DoubleJetC72",DoubleJet1 >= 72. && DoubleJet2 >= 72.);
  InsertInMenu("L1_DoubleJetC84",DoubleJet1 >= 84. && DoubleJet2 >= 84.);
  InsertInMenu("L1_DoubleJetC100",DoubleJet1 >= 100. && DoubleJet2 >= 100.);
  InsertInMenu("L1_DoubleJetC112",DoubleJet1 >= 112. && DoubleJet2 >= 112.);
  InsertInMenu("L1_DoubleJetC120",DoubleJet1 >= 120. && DoubleJet2 >= 120.);

  InsertInMenu("L1_DoubleIsoTau36er", algoFactory->DoubleTauJetEta2p17(36.,36.,true));
  InsertInMenu("L1_DoubleIsoTau40er", algoFactory->DoubleTauJetEta2p17(40.,40.,true));
  InsertInMenu("L1_DoubleIsoTau44er", algoFactory->DoubleTauJetEta2p17(44.,44.,true));

  InsertInMenu("L1_TripleJet_92_76_64", algoFactory->TripleJet(92.,76.,64.,false));

  InsertInMenu("L1_QuadJetC40", algoFactory->QuadJet(40.,40.,40.,40.,true));
  InsertInMenu("L1_QuadJetC60", algoFactory->QuadJet(60.,60.,60.,60.,true));
  InsertInMenu("L1_QuadJetC84", algoFactory->QuadJet(84.,84.,84.,84.,true));

  InsertInMenu("L1_QuadJetC36_Tau52", algoFactory->QuadJetCentral_TauJet(36.,52.));

  Int_t NN = insert_ibin;

  Int_t kOFFSET_old = kOFFSET;
  for (Int_t k=0; k < NN; k++) {
    TheTriggerBits[k + kOFFSET_old] = insert_val[k];
  }
  kOFFSET += insert_ibin;


  if (first) {

    NBITS_MULTIJETS = NN;

    for (Int_t ibin=0; ibin < insert_ibin; ibin++) {
      TString l1name = (TString)insert_names[ibin];
      h_MultiJets -> GetXaxis() -> SetBinLabel(ibin+1, l1name );
    }

    h_MultiJets-> GetXaxis() -> SetBinLabel(NN+1,"MULTIJETS");

    for (Int_t k=1; k <= kOFFSET -kOFFSET_old; k++) {
      h_All -> GetXaxis() -> SetBinLabel(k +kOFFSET_old , h_MultiJets -> GetXaxis() -> GetBinLabel(k) );
    }

  }

  Bool_t res = false;
  for (Int_t i=0; i < NN; i++) {
    res = res || insert_val[i] ;
    if (insert_val[i]) h_MultiJets -> Fill(i);
  }
  if (res) h_MultiJets -> Fill(NN);

  return res;
}

Bool_t L1Menu2012::Sums() {

  insert_ibin = 0;

  InsertInMenu("L1_ETM30",  algoFactory->ETM(30.));
  InsertInMenu("L1_ETM40",  algoFactory->ETM(40.));
  InsertInMenu("L1_ETM50",  algoFactory->ETM(50.));
  InsertInMenu("L1_ETM60",  algoFactory->ETM(60.));
  InsertInMenu("L1_ETM70",  algoFactory->ETM(70.));
  InsertInMenu("L1_ETM100", algoFactory->ETM(100.));

  InsertInMenu("L1_ETM60_NoJet52WdPhi2",  algoFactory->ETM_NoQCD(60.));
  InsertInMenu("L1_ETM70_NoJet52WdPhi2",  algoFactory->ETM_NoQCD(70.));

  InsertInMenu("L1_HTT125", algoFactory->HTT(125.));
  InsertInMenu("L1_HTT150", algoFactory->HTT(150.));
  InsertInMenu("L1_HTT175", algoFactory->HTT(175.));
  InsertInMenu("L1_HTT200", algoFactory->HTT(200.));
  InsertInMenu("L1_HTT250", algoFactory->HTT(250.));

  Int_t NN = insert_ibin;

  Int_t kOFFSET_old = kOFFSET;
  for (Int_t k=0; k < NN; k++) {
    TheTriggerBits[k + kOFFSET_old] = insert_val[k];
  }
  kOFFSET += insert_ibin;


  if (first) {

    NBITS_SUMS = NN;

    for (Int_t ibin=0; ibin < insert_ibin; ibin++) {
      TString l1name = (TString)insert_names[ibin]; 
      h_Sums -> GetXaxis() -> SetBinLabel(ibin+1, l1name );
    }

    h_Sums -> GetXaxis() -> SetBinLabel(NN+1,"SUMS");

    for (Int_t k=1; k <= kOFFSET -kOFFSET_old; k++) {
      h_All -> GetXaxis() -> SetBinLabel(k +kOFFSET_old , h_Sums -> GetXaxis() -> GetBinLabel(k) );
    }
  }

  Bool_t res = false; 
  for (Int_t i=0; i < NN; i++) {
    res = res || insert_val[i] ; 
    if (insert_val[i]) h_Sums -> Fill(i);
  }
  if (res) h_Sums -> Fill(NN);

  return res;
}

Bool_t L1Menu2012::EGamma() {

  insert_ibin = 0;

  Float_t SingleEGPt = -10.;
  algoFactory->SingleEGPt(SingleEGPt);

  Float_t SingleIsoEGerPt = -10.;
  algoFactory->SingleEGEta2p1Pt(SingleIsoEGerPt,true);

  InsertInMenu("L1_SingleEG5",SingleEGPt >= 5.);
  InsertInMenu("L1_SingleEG10",SingleEGPt >= 10.);
  InsertInMenu("L1_SingleEG15",SingleEGPt >= 15.);
  InsertInMenu("L1_SingleEG20",SingleEGPt >= 20.);
  InsertInMenu("L1_SingleEG25",SingleEGPt >= 25.);
  InsertInMenu("L1_SingleEG30",SingleEGPt >= 30.);
  InsertInMenu("L1_SingleEG35",SingleEGPt >= 35.);
  InsertInMenu("L1_SingleEG40",SingleEGPt >= 40.);
  InsertInMenu("L1_SingleEG35er", algoFactory->SingleEGEta2p1(35.) );
  InsertInMenu("L1_SingleIsoEG18", algoFactory->SingleEG(18.,true) );
  InsertInMenu("L1_SingleIsoEG25", algoFactory->SingleEG(25.,true) );
  InsertInMenu("L1_SingleIsoEG20er",SingleIsoEGerPt >= 20.);
  InsertInMenu("L1_SingleIsoEG22er",SingleIsoEGerPt >= 22.);
  InsertInMenu("L1_SingleIsoEG25er",SingleIsoEGerPt >= 25.);
  InsertInMenu("L1_SingleIsoEG28er",SingleIsoEGerPt >= 28.);
  InsertInMenu("L1_SingleIsoEG30er",SingleIsoEGerPt >= 30.);

  Int_t NN = insert_ibin;

  Int_t kOFFSET_old = kOFFSET;
  for (Int_t k=0; k < NN; k++) {
    TheTriggerBits[k + kOFFSET_old] = insert_val[k];
  }
  kOFFSET += insert_ibin;

  if (first) {

    NBITS_EGAMMA = NN;

    for (Int_t ibin=0; ibin < insert_ibin; ibin++) {
      TString l1name = (TString)insert_names[ibin];
      h_Egamma -> GetXaxis() -> SetBinLabel(ibin+1, l1name );
    }
    h_Egamma-> GetXaxis() -> SetBinLabel(NN+1,"EGAMMA");

    for (Int_t k=1; k <= kOFFSET -kOFFSET_old; k++) {
      h_All -> GetXaxis() -> SetBinLabel(k +kOFFSET_old , h_Egamma -> GetXaxis() -> GetBinLabel(k) );
    }
  }                      

  Bool_t res = false;      
  for (Int_t i=0; i < NN; i++) {
    res = res || insert_val[i] ;
    if (insert_val[i]) h_Egamma -> Fill(i);
  }      
  if (res) h_Egamma -> Fill(NN);

  return res;
}       

Bool_t L1Menu2012::MultiEGamma() {

  insert_ibin = 0;

  InsertInMenu("L1_DoubleEG_15_10", algoFactory->DoubleEG(15.,10.) );
  InsertInMenu("L1_DoubleEG_22_10", algoFactory->DoubleEG(22.,10.) );
  InsertInMenu("L1_DoubleEG_20_10_1LegIso", algoFactory->DoubleEG(20.,10.,true) );
  InsertInMenu("L1_TripleEG_14_10_8", algoFactory->TripleEG(14.,10.,8.) );

  Int_t NN = insert_ibin;

  Int_t kOFFSET_old = kOFFSET;
  for (Int_t k=0; k < NN; k++) {
    TheTriggerBits[k + kOFFSET_old] = insert_val[k];
  }
  kOFFSET += insert_ibin;


  if (first) {

    NBITS_MULTIEGAMMA = NN;

    for (Int_t ibin=0; ibin < insert_ibin; ibin++) {
      TString l1name = (TString)insert_names[ibin];
      h_MultiEgamma -> GetXaxis() -> SetBinLabel(ibin+1, l1name );
    }
    h_MultiEgamma-> GetXaxis() -> SetBinLabel(NN+1,"MULTIEGAMMA");

    for (Int_t k=1; k <= kOFFSET -kOFFSET_old; k++) {
      h_All -> GetXaxis() -> SetBinLabel(k +kOFFSET_old , h_MultiEgamma -> GetXaxis() -> GetBinLabel(k) );
    }
  }                      

  Bool_t res = false;      
  for (Int_t i=0; i < NN; i++) {
    res = res || insert_val[i] ;
    if (insert_val[i]) h_MultiEgamma -> Fill(i);
  }      
  if (res) h_MultiEgamma -> Fill(NN);

  return res;
}       

Bool_t L1Menu2012::Technical() {

  insert_ibin = 0;

  InsertInMenu("L1_ZeroBias", 1 );

  Int_t NN = insert_ibin;

  Int_t kOFFSET_old = kOFFSET;
  for (Int_t k=0; k < NN; k++) {
    TheTriggerBits[k + kOFFSET_old] = insert_val[k];
  }
  kOFFSET += insert_ibin;

  if (first) {
    for (Int_t ibin=0; ibin < insert_ibin; ibin++) {
      TString l1name = (TString)insert_names[ibin];
      h_Technical->GetXaxis()->SetBinLabel(ibin+1, l1name );
    }
    h_Technical-> GetXaxis() -> SetBinLabel(NN+1,"Technical");

    for (Int_t k=1; k <= kOFFSET -kOFFSET_old; k++) {
      h_All->GetXaxis()->SetBinLabel(k +kOFFSET_old , h_Technical->GetXaxis()->GetBinLabel(k) );
    }
  }                      

  Bool_t res = false;      
  for (Int_t i=0; i < NN; i++) {
    res = res || insert_val[i] ;
    if (insert_val[i]) h_Technical -> Fill(i);
  }      
  if(res)h_Technical -> Fill(NN);

  return res;
}       


void L1Menu2012::Loop() {

  Int_t nevents = GetEntries();
  Double_t nZeroBiasevents = 0.;
  if(nevents > 6000000) nevents = 6000000;

  Int_t NPASS = 0; 

  Int_t NJETS = 0;
  Int_t MULTINJETS = 0;
  Int_t NEG = 0;
  Int_t MULTINEG = 0;
  Int_t NSUMS =0;
  Int_t NMUONS = 0;
  Int_t MULTINMUONS = 0;
  Int_t NCROSS = 0;
  Int_t MULTINCROSS = 0;
  Int_t TECHNICAL = 0;

  Int_t nPAG      = 0;
  Int_t nTRIGPHYS = 0;
	
  first = true;

  for (Long64_t i=0; i<nevents; i++){     
    Long64_t ientry = LoadTree(i); if (ientry < 0) break;
    GetEntry(i);

    FilL1Bits();
    if(first) MyInit();

    Bool_t raw = PhysicsBits[0];  // check for ZeroBias triggered events
    if(!raw) continue;

    nZeroBiasevents++;

    // reset the emulated trigger bits
    kOFFSET = 0;
    for (Int_t k=0; k < N128; k++) {
      TheTriggerBits[k] = false;
    }

    Bool_t cross       = Cross();
    Bool_t multicross  = MultiCross();
    Bool_t eg          = EGamma();
    Bool_t multieg     = MultiEGamma();
    Bool_t muons       = Muons();
    Bool_t multimuons  = MultiMuons();
    Bool_t jets        = Jets() ;
    Bool_t multijets   = MultiJets() ;
    Bool_t sums        = Sums();
    Bool_t technical   = Technical();

    Bool_t pass  = jets || multijets || eg || multieg || sums || muons || multimuons || cross || multicross || technical;

    if(pass) NPASS ++;

    if(cross)      NCROSS++;
    if(multicross) MULTINCROSS++;
    if(muons)      NMUONS++;
    if(multimuons) MULTINMUONS++;
    if(sums)       NSUMS++;
    if(eg)         NEG++;
    if(multieg)    MULTINEG++;
    if(jets)       NJETS++;
    if(multijets)  MULTINJETS++;
    if(technical)  TECHNICAL++;

    if(pass) h_Block->Fill(10.);

    Bool_t dec[10];
    dec[0] = eg;
    dec[1] = multieg;
    dec[2] = jets;
    dec[3] = multijets;
    dec[4] = muons;
    dec[5] = multimuons;
    dec[6] = sums;
    dec[7] = cross;
    dec[8] = multicross;
    dec[9] = technical;

    for (Int_t l=0; l < 9; l++) {
      if(dec[l]){
	h_Block -> Fill(l);
	for (Int_t k=0; k < 5; k++) {
	  if (dec[k]) cor_Block -> Fill(l,k);
	}
      }

    }

    first = false;

    // now the pure rate stuff
    // kOFFSET now contains the number of triggers we have calculated

    Bool_t ddd[NPAGS];
    for (Int_t idd=0; idd < NPAGS; idd++) {
      ddd[idd] = false; 
    } 

    Bool_t eee[NTRIGPHYS];
    for (Int_t iee=0; iee < NTRIGPHYS; iee++) {
      eee[iee] = false; 
    }

    Float_t weightEventPAGs = 1.;
    Float_t weightEventTRIGPHYS = 1.;

    for (Int_t k=0; k < kOFFSET; k++) {
      if ( ! TheTriggerBits[k] ) continue;
      h_All -> Fill(k);

      TString name = h_All -> GetXaxis() -> GetBinLabel(k+1);
      std::string L1namest = (std::string)name;

      Bool_t IsTOP   = setTOP.count(L1namest) > 0;
      Bool_t IsHIGGS = setHIGGS.count(L1namest) > 0;
      Bool_t IsBPH   = setBPH.count(L1namest) > 0;
      Bool_t IsEXO   = setEXO.count(L1namest) > 0;
      Bool_t IsSUSY  = setSUSY.count(L1namest) > 0;
      Bool_t IsSMP   = setSMP.count(L1namest) > 0;
      Bool_t IsB2G   = setB2G.count(L1namest) > 0;
      if(IsHIGGS) ddd[0] = true;
      if(IsSUSY)  ddd[1] = true;
      if(IsEXO)   ddd[2] = true;
      if(IsTOP)   ddd[3] = true;
      if(IsSMP)   ddd[4] = true;
      if(IsBPH)   ddd[5] = true;
      if(IsB2G)   ddd[6] = true;

      Float_t ww = WeightsPAGs[L1namest];
      if (ww < weightEventPAGs) weightEventPAGs = ww;

      Bool_t IsMuon     = setMuon.count(L1namest) > 0;
      Bool_t IsEG       = setEG.count(L1namest) > 0;
      Bool_t IsHadronic = setHadronic.count(L1namest) > 0;

      Bool_t IsMuonEG       = setMuonEG.count(L1namest) > 0;
      Bool_t IsMuonHadronic = setMuonHadronic.count(L1namest) > 0;
      Bool_t IsEGHadronic   = setEGHadronic.count(L1namest) > 0;

      if(IsMuon)     eee[0] = true;
      if(IsEG)       eee[1] = true;
      if(IsHadronic) eee[2] = true;

      if(IsMuonEG)       eee[3] = true;
      if(IsMuonHadronic) eee[4] = true;
      if(IsEGHadronic)   eee[5] = true;

      Float_t www = WeightsTRIGPHYS[L1namest];
      if(www < weightEventTRIGPHYS) weightEventTRIGPHYS = www;

      //did the event pass another trigger ?
      Bool_t pure = true;
      for (Int_t k2=0; k2 < kOFFSET; k2++) {
	if (k2 == k) continue;
	if ( TheTriggerBits[k2] ) pure = false;
      }
      if (pure) h_Pure -> Fill(k);
    }

    // for the PAG rates
    Bool_t PAG = false;
    for (Int_t idd=0; idd < NPAGS; idd++) {
      if (ddd[idd]) {
	Bool_t pure = true;
	PAG = true;
	for (Int_t jdd=0; jdd < NPAGS; jdd++) {
	  if (ddd[jdd]) {
	    cor_PAGS -> Fill(idd,jdd);
	    if (jdd != idd) pure = false;
	  }
	}   
	if(pure) h_PAGS_pure -> Fill(idd);
	h_PAGS_shared -> Fill(idd,weightEventPAGs);

      }  
    }
    if(PAG) nPAG ++;

    //for the TRIGPHYS rates :
    Bool_t TRIGPHYS = false;
    for (Int_t iee=0; iee < NTRIGPHYS; iee++) {
      if (eee[iee]) {
	Bool_t pure = true;
	TRIGPHYS = true;
	for (Int_t jee=0; jee < NTRIGPHYS; jee++) {
	  if (eee[jee]) {
	    cor_TRIGPHYS -> Fill(iee,jee);
	    if (jee != iee) pure = false;
	  }
	}   
	if(pure) h_TRIGPHYS_pure -> Fill(iee);
	h_TRIGPHYS_shared -> Fill(iee,weightEventTRIGPHYS);

      }  
    }
    if(TRIGPHYS) nTRIGPHYS++;

  }  // end evt loop

  std::cout << std::endl << std::endl << " ... number of usable events: " << nevents << std::endl;
  std::cout << std::endl << "... number of zero bias triggered events: " << nZeroBiasevents << std::endl << std::endl;

  Float_t scal = 11246.; // ZB per bunch in kHz
  scal /= nZeroBiasevents*1000.;
  scal *= theNumberOfBunches;

  Float_t extrarate = 10.;

  h_Cross       -> Scale(scal);
  h_MultiCross  -> Scale(scal);
  h_Jets        -> Scale(scal);
  h_MultiJets   -> Scale(scal);
  h_Egamma      -> Scale(scal);
  h_MultiEgamma -> Scale(scal);
  h_Sums        -> Scale(scal);
  h_Muons       -> Scale(scal);
  h_MultiMuons  -> Scale(scal);
  h_Technical   -> Scale(scal);
  CorrectScale(h_All, scal);
  h_Pure  -> Scale(scal);

  std::cout << " Prescales for: " << themenufilename << std::endl;
  std::cout << std::endl << " --------------------------------------------------------- " << std::endl << std::endl;
  std::cout << " Rate that pass L1 " << NPASS * scal << " kHz  ( claimed by a PAG " << nPAG * scal << " kHz  i.e. " << 100.*(float)nPAG/(float)NPASS << "%. ) " << ") adding " << extrarate << " kHz = " << NPASS * scal + extrarate << " kHz " << std::endl << std::endl;
  std::cout << " --------------------------------------------------------- " << std::endl << std::endl;
  std::cout << " Rate that pass L1 jets: " << NJETS * scal << " kHz" << std::endl;
  std::cout << " Rate that pass L1 multi-jets: " << MULTINJETS * scal << " kHz" << std::endl;
  std::cout << " Rate that pass L1 EG: " << NEG * scal << " kHz" << std::endl;
  std::cout << " Rate that pass L1 multi-EG: " << MULTINEG * scal << " kHz" << std::endl;
  std::cout << " Rate that pass L1 Sums: " << NSUMS * scal << " kHz" << std::endl;
  std::cout << " Rate that pass L1 Muons: " << NMUONS * scal << " kHz" << std::endl;
  std::cout << " Rate that pass L1 multi-Muons: " << MULTINMUONS * scal << " kHz" << std::endl;
  std::cout << " Rate that pass L1 Cross: " << NCROSS * scal << " kHz" << std::endl;
  std::cout << " Rate that pass L1 multi-Cross: " << MULTINCROSS * scal << " kHz" << std::endl;
  std::cout << " Rate that pass L1 Technical: " << TECHNICAL * scal << " kHz" << std::endl;

  output << " Prescales for: L1NtupleFileName = " << themenufilename << std::endl;
  output << std::endl << " --------------------------------------------------------- " << std::endl << std::endl;
  output << " Rate that pass L1 " << NPASS * scal << " kHz  ( claimed by a PAG " << nPAG * scal << " kHz  i.e. " << 100.*(float)nPAG/(float)NPASS << "%. ) and adding " << extrarate << " kHz = " << NPASS * scal + extrarate << " kHz " << std::endl << std::endl;
  output << " --------------------------------------------------------- " << std::endl << std::endl;
  output << " Rate that pass L1 jets: " << NJETS * scal << " kHz" << std::endl;
  output << " Rate that pass L1 multi-jets: " << MULTINJETS * scal << " kHz" << std::endl;
  output << " Rate that pass L1 EG: " << NEG * scal << " kHz" << std::endl;
  output << " Rate that pass L1 multi-EG: " << MULTINEG * scal << " kHz" << std::endl;
  output << " Rate that pass L1 Sums: " << NSUMS * scal << " kHz" << std::endl;
  output << " Rate that pass L1 Muons: " << NMUONS * scal << " kHz" << std::endl;
  output << " Rate that pass L1 multi-Muons: " << MULTINMUONS * scal << " kHz" << std::endl;
  output << " Rate that pass L1 Cross: " << NCROSS * scal << " kHz" << std::endl;
  output << " Rate that pass L1 multi-Cross: " << MULTINCROSS * scal << " kHz" << std::endl;
  output << " Rate that pass L1 Technical: " << TECHNICAL * scal << " kHz" << std::endl;

  for (Int_t i=1; i<= 5; i++) {
    Float_t nev = h_Block -> GetBinContent(i);
    for (Int_t j=1; j<= 5; j++) {
      Int_t ibin = cor_Block -> FindBin(i-1,j-1);
      Float_t val = cor_Block -> GetBinContent(ibin);
      val = val / nev;
      cor_Block -> SetBinContent(ibin,val);
    }
  }   

  h_Block->Scale(scal);

  cor_PAGS->Scale(scal);
  h_PAGS_pure->Scale(scal);
  h_PAGS_shared->Scale(scal);

  cor_TRIGPHYS->Scale(scal);
  h_TRIGPHYS_pure->Scale(scal);
  h_TRIGPHYS_shared->Scale(scal);

  Int_t NBITS_ALL = NBITS_MUONS + NBITS_MULTIMUONS + NBITS_EGAMMA + NBITS_MULTIEGAMMA + NBITS_JETS + NBITS_MULTIJETS + NBITS_SUMS + NBITS_CROSS + NBITS_MULTICROSS;

  std::cout << std::endl << " --- TOTAL NUMBER OF BITS: " << std::endl;
  std::cout << "  USED : " << NBITS_ALL << std::endl;
  std::cout << "  MUONS : " << NBITS_MUONS << std::endl;
  std::cout << "  MULTIMUONS : " << NBITS_MULTIMUONS << std::endl;
  std::cout << "  EGAMMA : " << NBITS_EGAMMA << std::endl;
  std::cout << "  MULTIEGAMMA : " << NBITS_MULTIEGAMMA << std::endl;
  std::cout << "  JETS : " << NBITS_JETS << std::endl;
  std::cout << "  MULTIJETS : " << NBITS_MULTIJETS << std::endl;
  std::cout << "  SUMS : " << NBITS_SUMS << std::endl;
  std::cout << "  CROSS : " << NBITS_CROSS << std::endl;
  std::cout << "  MULTICROSS : " << NBITS_MULTICROSS << std::endl << std::endl;

  output << std::endl << " --- TOTAL NUMBER OF BITS: " << std::endl;
  output << "  USED : " << NBITS_ALL << std::endl;
  output << "  MUONS : " << NBITS_MUONS << std::endl;
  output << "  MULTIMUONS : " << NBITS_MULTIMUONS << std::endl;
  output << "  EGAMMA : " << NBITS_EGAMMA << std::endl;
  output << "  MULTIEGAMMA : " << NBITS_MULTIEGAMMA << std::endl;
  output << "  JETS : " << NBITS_JETS << std::endl;
  output << "  MULTIJETS : " << NBITS_MULTIJETS << std::endl;
  output << "  SUMS : " << NBITS_SUMS << std::endl;
  output << "  CROSS : " << NBITS_CROSS << std::endl;
  output << "  MULTICROSS : " << NBITS_MULTICROSS << std::endl << std::endl;

  cout << "##########################" << endl;
  if(isOneSeedNotDefined) cout << "Warning, undefined seeds! Perhaps is normal, but check!" << endl;
  cout << "Missing seeds are in bit numbers ";
  for(int i=0;i<N128;i++){
    if(MissingBits[i] == 1) cout << i << ", ";
  }
  cout << endl << "##########################" << endl;
}

void RunL1(Bool_t drawplots=true, Bool_t writefiles=true, Int_t whichFileAndLumiToUse=1){

  Int_t Nbin_max = 50;
  h_Cross       = new TH1F("h_Cross","h_Cross",Nbin_max,-0.5,(float)Nbin_max-0.5);
  h_MultiCross  = new TH1F("h_MultiCross","h_MultiCross",Nbin_max,-0.5,(float)Nbin_max-0.5);
  h_Sums        = new TH1F("h_Sums","h_Sums",Nbin_max,-0.5,(float)Nbin_max-0.5);
  h_Jets        = new TH1F("h_Jets","h_Jets",Nbin_max,-0.5,(float)Nbin_max-0.5);
  h_MultiJets   = new TH1F("h_MultiJets","h_MultiJets",Nbin_max,-0.5,(float)Nbin_max-0.5);
  h_Egamma      = new TH1F("h_Egamma","h_Egamma",Nbin_max,-0.5,(float)Nbin_max-0.5);
  h_MultiEgamma = new TH1F("h_MultiEgamma","h_MultiEgamma",Nbin_max,-0.5,(float)Nbin_max-0.5);
  h_Muons       = new TH1F("h_Muons","h_Muons",Nbin_max,-0.5,(float)Nbin_max-0.5);
  h_MultiMuons  = new TH1F("h_MultiMuons","h_MultiMuons",Nbin_max,-0.5,(float)Nbin_max-0.5);
  h_Technical   = new TH1F("h_Technical","h_Technical",Nbin_max,-0.5,(float)Nbin_max-0.5);

  h_Block   = new TH1F("h_Block","h_Block",11,-0.5,10.5);
  cor_Block = new TH2F("cor_Block","cor_Block",10,-0.5,9.5,19,-0.5,9.5);

  cor_PAGS      = new TH2F("cor_PAGS","cor_PAGS",NPAGS,-0.5,(float)NPAGS-0.5,NPAGS,-0.5,(float)NPAGS-0.5);
  h_PAGS_pure   = new TH1F("h_PAGS_pure","h_PAGS_pure",NPAGS,-0.5,(float)NPAGS-0.5);
  h_PAGS_shared = new TH1F("h_PAGS_shared","h_PAGS_shared",NPAGS,-0.5,(float)NPAGS-0.5);

  cor_TRIGPHYS      = new TH2F("cor_TRIGPHYS","cor_TRIGPHYS",NTRIGPHYS,-0.5,(float)NTRIGPHYS-0.5,NTRIGPHYS,-0.5,(float)NTRIGPHYS-0.5);
  h_TRIGPHYS_pure   = new TH1F("h_TRIGPHYS_pure","h_TRIGPHYS_pure",NTRIGPHYS,-0.5,(float)NTRIGPHYS-0.5);
  h_TRIGPHYS_shared = new TH1F("h_TRIGPHYS_shared","h_TRIGPHYS_shared",NTRIGPHYS,-0.5,(float)NTRIGPHYS-0.5);

  h_All  = new TH1F("h_All","h_All",N128,-0.5,N128-0.5);
  h_Pure = new TH1F("h_Pure","h_Pure",N128,-0.5,N128-0.5);

  Float_t NumberOfBunches = 0.; 
  std::string L1NtupleFileName = "";
  Float_t AveragePU = 0.;
  Bool_t L1JetCorrection = false;
  Bool_t noHF = false;
  Bool_t noTauInJet = false;
  Float_t Energy = 0.;
  Float_t targetlumi = 1.;

  string themenufilename;

  if(whichFileAndLumiToUse==1){
    // 13 TeV ZeroBias 72X sample 30 PU 50 ns, 2012 re-emulation with 10 GeV cut on jet seed
    NumberOfBunches = 1368; 
    L1NtupleFileName = "root://lxcms02//data2/p/pellicci/L1DPG/root/v11/50ns_30PU_ReEmul2012Gct10GeV/L1Tree.root";
    themenufilename = "Menu_30PU_50bx.txt";
    //themenufilename = "Menu_Noprescales.txt";
    AveragePU = 30;
    L1JetCorrection=false;
    Energy = 13;
    noHF = false;
    targetlumi= 70.;
  }
  else if(whichFileAndLumiToUse==2){
    // 13 TeV ZeroBias 62X sample 40 PU 50 ns, 2012 re-emulation with 10 GeV cut on jet seed
    NumberOfBunches = 1368; 
    L1NtupleFileName = "root://lxcms02//data2/p/pellicci/L1DPG/root/v9/50ns_40PU_ReEmul2012Gct10GeV/L1Tree.root";
    themenufilename = "Menu_40PU_50bx.txt";
    //themenufilename = "Menu_Noprescales.txt";
    AveragePU = 40;
    L1JetCorrection=false;
    Energy = 13;
    noHF = false;
    targetlumi= 70.;
  }
  else if(whichFileAndLumiToUse==3){
    // 13 TeV ZeroBias 62X sample 20PU 25 ns, 2015 re-emulation
    NumberOfBunches = 2508; 
    L1NtupleFileName = "root://lxcms02//data2/p/pellicci/L1DPG/root/v11/25ns_20PU_ReEmul2015/L1Tree.root";
    themenufilename = "Menu_20PU_25bx.txt";
    //themenufilename = "Menu_Noprescales.txt";
    AveragePU = 20;
    L1JetCorrection=false;
    Energy = 13;
    noTauInJet = true;
    targetlumi= 70.;
  }
  else if(whichFileAndLumiToUse==4){
    // 13 TeV ZeroBias 62X sample 40 PU 25 ns, 2015 re-emulation
    NumberOfBunches = 2508; 
    L1NtupleFileName = "root://lxcms02//data2/p/pellicci/L1DPG/root/v11/25ns_40PU_ReEmul2015/L1Tree.root";
    themenufilename = "Menu_40PU_25bx.txt";
    //themenufilename = "Menu_Noprescales.txt";
    AveragePU = 40;
    L1JetCorrection=false;
    Energy = 13;
    noTauInJet = true;
    targetlumi= 140.;
  }
  else{
    std::cout << "##########################" << std::endl;
    std::cout << "ERROR: Please define a ntuple file which is in the allowed range! You did use: whichFileAndLumiToUse = " << whichFileAndLumiToUse << " This is not in the allowed range" << std::endl;
    std::cout << "##########################" << std::endl;
  }

  std::stringstream txtos;
  txtos << "L1Menu_" << Energy << "_" << AveragePU << "_rates.txt";
  TString TXTOutPutFileName = txtos.str();
  std::ofstream TXTOutfile(TXTOutPutFileName);

  std::stringstream csvos;
  csvos << "L1Menu_" << Energy << "_" << AveragePU << "_rates.csv";
  TString CSVOutPutFileName = csvos.str();
  std::ofstream CSVOutfile(CSVOutPutFileName);

  std::cout << std::endl << "L1 menu used = " << themenufilename << std::endl;
  std::cout << std::endl << "Using: whichFileAndLumiToUse = " << whichFileAndLumiToUse << std::endl;
  std::cout << "  L1NtupleFileName                 = " << L1NtupleFileName << std::endl;
  std::cout << "  L1MenuFileName                   = " << themenufilename << std::endl;
  std::cout << "  AveragePU                        = " << AveragePU << std::endl;

  if (writefiles) { 
    std::cout << std::endl << "Writing CSV and txt files as well ..." << std::endl;
    std::cout << "  TXTOutPutFileName: " << TXTOutPutFileName << std::endl;
    std::cout << "  CVSOutPutFileName: " << CSVOutPutFileName << std::endl << std::endl;

    TXTOutfile << std::endl << "L1 menu used = " << themenufilename << std::endl;

    TXTOutfile << std::endl << "Target Luminosity = " << targetlumi << std::endl << std::endl;

    TXTOutfile << std::endl << "Using: whichFileAndLumiToUse = " << whichFileAndLumiToUse << std::endl;
    TXTOutfile << "  L1NtupleFileName                 = " << L1NtupleFileName << std::endl;
  }

  L1Menu2012 a(themenufilename,NumberOfBunches,L1JetCorrection,noHF,noTauInJet,AveragePU);
  a.Open(L1NtupleFileName);
  a.Loop();

  if(drawplots){

    TFile fOut("plots_menu.root","RECREATE");
    fOut.cd();

    TString YaxisName;
    if(targetlumi == 1.)   YaxisName = "Rate (kHz) at 1e32";
    if(targetlumi == 2.)   YaxisName = "Rate (kHz) at 2e32";
    if(targetlumi == 50)   YaxisName = "Rate (kHz) at 5e33";
    if(targetlumi == 60.)  YaxisName = "Rate (kHz) at 6e33";
    if(targetlumi == 70.)  YaxisName = "Rate (kHz) at 7e33";
    if(targetlumi == 140.) YaxisName = "Rate (kHz) at 1.4e34";

    TCanvas* c1 = new TCanvas("c1","c1");
    c1 -> cd();
    gStyle -> SetOptStat(0);
    h_Cross -> SetLineColor(4);
    h_Cross -> GetXaxis() -> SetLabelSize(0.035);
    h_Cross -> SetYTitle(YaxisName);
    h_Cross -> Draw();
    c1->Write();

    TCanvas* c2 = new TCanvas("c2","c2");
    c2 -> cd();
    h_MultiCross -> SetLineColor(4);
    h_MultiCross -> GetXaxis() -> SetLabelSize(0.035);
    h_MultiCross -> SetYTitle(YaxisName);
    h_MultiCross -> Draw();
    c2->Write();

    TCanvas* c3 = new TCanvas("c3","c3");
    c3 -> cd();
    h_Sums -> SetLineColor(4);
    h_Sums -> GetXaxis() -> SetLabelSize(0.035);
    h_Sums -> SetYTitle(YaxisName);
    h_Sums -> Draw();
    c3->Write();

    TCanvas* c4 = new TCanvas("c4","c4");
    c4 -> cd();
    h_Egamma -> SetLineColor(4);
    h_Egamma -> GetXaxis() -> SetLabelSize(0.035);
    h_Egamma -> SetYTitle(YaxisName);
    h_Egamma -> Draw();
    c4->Write();

    TCanvas* c5 = new TCanvas("c5","c5");
    c5 -> cd();
    h_MultiEgamma -> SetLineColor(4);
    h_MultiEgamma -> GetXaxis() -> SetLabelSize(0.035);
    h_MultiEgamma -> SetYTitle(YaxisName);
    h_MultiEgamma -> Draw();
    c5->Write();

    TCanvas* c6 = new TCanvas("c6","c6");
    c6 -> cd();
    h_Jets -> SetLineColor(4);
    h_Jets -> GetXaxis() -> SetLabelSize(0.035);
    h_Jets -> SetYTitle(YaxisName);
    h_Jets -> Draw();
    c6->Write();

    TCanvas* c7 = new TCanvas("c7","c7");
    c7 -> cd();
    h_MultiJets -> SetLineColor(4);
    h_MultiJets -> GetXaxis() -> SetLabelSize(0.035);
    h_MultiJets -> SetYTitle(YaxisName);
    h_MultiJets -> Draw();
    c7->Write();

    TCanvas* c8 = new TCanvas("c8","c8");
    c8 -> cd();
    h_Muons -> SetLineColor(4);
    h_Muons -> GetXaxis() -> SetLabelSize(0.035);
    h_Muons -> SetYTitle(YaxisName);
    h_Muons -> Draw();
    c8->Write();

    TCanvas* c9 = new TCanvas("c9","c9");
    c9 -> cd();
    h_MultiMuons -> SetLineColor(4);
    h_MultiMuons -> GetXaxis() -> SetLabelSize(0.035);
    h_MultiMuons -> SetYTitle(YaxisName);
    h_MultiMuons -> Draw();
    c9->Write();

    TCanvas* c10 = new TCanvas("c10","c10");
    c10 -> cd();
    cor_Block -> GetXaxis() -> SetBinLabel(1,"EG");
    cor_Block -> GetXaxis() -> SetBinLabel(2,"MultiEG");
    cor_Block -> GetXaxis() -> SetBinLabel(3,"Jets");
    cor_Block -> GetXaxis() -> SetBinLabel(4,"MultiJets");
    cor_Block -> GetXaxis() -> SetBinLabel(5,"Muons");
    cor_Block -> GetXaxis() -> SetBinLabel(6,"MultiMuons");
    cor_Block -> GetXaxis() -> SetBinLabel(7,"Sums");
    cor_Block -> GetXaxis() -> SetBinLabel(8,"Cross");
    cor_Block -> GetXaxis() -> SetBinLabel(9,"MultiCross");
    cor_Block -> GetXaxis() -> SetBinLabel(10,"Technical");

    cor_Block -> GetYaxis() -> SetBinLabel(1,"EG");
    cor_Block -> GetYaxis() -> SetBinLabel(2,"MultiEG");
    cor_Block -> GetYaxis() -> SetBinLabel(3,"Jets");
    cor_Block -> GetYaxis() -> SetBinLabel(4,"MultiJets");
    cor_Block -> GetYaxis() -> SetBinLabel(5,"Muons");
    cor_Block -> GetYaxis() -> SetBinLabel(6,"MultiMuons");
    cor_Block -> GetYaxis() -> SetBinLabel(7,"Sums");
    cor_Block -> GetYaxis() -> SetBinLabel(8,"Cross");
    cor_Block -> GetYaxis() -> SetBinLabel(9,"MultiCross");
    cor_Block -> GetYaxis() -> SetBinLabel(10,"Technical");

    cor_Block -> Draw("colz");
    cor_Block -> Draw("same,text");
    c10->Write();

    TCanvas* c11 = new TCanvas("c11","c11");
    c11 -> cd();
    cor_PAGS -> GetXaxis() -> SetBinLabel(1,"HIGGS");
    cor_PAGS -> GetXaxis() -> SetBinLabel(2,"SUSY");
    cor_PAGS -> GetXaxis() -> SetBinLabel(3,"EXO");
    cor_PAGS -> GetXaxis() -> SetBinLabel(4,"TOP");
    cor_PAGS -> GetXaxis() -> SetBinLabel(5,"SMP");
    cor_PAGS -> GetXaxis() -> SetBinLabel(6,"BPH");
    cor_PAGS -> GetXaxis() -> SetBinLabel(7,"B2G");

    cor_PAGS -> GetYaxis() -> SetBinLabel(1,"HIGGS");
    cor_PAGS -> GetYaxis() -> SetBinLabel(2,"SUSY");
    cor_PAGS -> GetYaxis() -> SetBinLabel(3,"EXO");
    cor_PAGS -> GetYaxis() -> SetBinLabel(4,"TOP");
    cor_PAGS -> GetYaxis() -> SetBinLabel(5,"SMP");
    cor_PAGS -> GetYaxis() -> SetBinLabel(6,"BPH");
    cor_PAGS -> GetYaxis() -> SetBinLabel(7,"B2G");

    cor_PAGS -> Draw("colz");
    cor_PAGS -> Draw("same,text"); 
    c11->Write();

    TCanvas* c12 = new TCanvas("c12","c12");
    c12 -> cd();
    h_PAGS_pure -> GetXaxis() -> SetBinLabel(1,"HIGGS");
    h_PAGS_pure -> GetXaxis() -> SetBinLabel(2,"SUSY");
    h_PAGS_pure -> GetXaxis() -> SetBinLabel(3,"EXO");
    h_PAGS_pure -> GetXaxis() -> SetBinLabel(4,"TOP");
    h_PAGS_pure -> GetXaxis() -> SetBinLabel(5,"SMP");
    h_PAGS_pure -> GetXaxis() -> SetBinLabel(6,"BPH");
    h_PAGS_pure -> GetXaxis() -> SetBinLabel(7,"B2G");
    h_PAGS_pure -> SetYTitle("Pure rate (kHz)");
    h_PAGS_pure -> Draw();
    c12->Write();

    TCanvas* c13 = new TCanvas("c13","c13");
    c13 -> cd();
    h_PAGS_shared -> GetXaxis() -> SetBinLabel(1,"HIGGS");
    h_PAGS_shared -> GetXaxis() -> SetBinLabel(2,"SUSY");
    h_PAGS_shared -> GetXaxis() -> SetBinLabel(3,"EXO");
    h_PAGS_shared -> GetXaxis() -> SetBinLabel(4,"TOP");
    h_PAGS_shared -> GetXaxis() -> SetBinLabel(5,"SMP");
    h_PAGS_shared -> GetXaxis() -> SetBinLabel(6,"BPH");
    h_PAGS_shared -> GetXaxis() -> SetBinLabel(7,"B2G");
    h_PAGS_shared -> SetYTitle("Shared rate (kHz)");
    h_PAGS_shared -> Draw();
    c13->Write();
		
    TCanvas* c14 = new TCanvas("c14","c14");
    c14 -> cd();
    cor_TRIGPHYS -> GetXaxis() -> SetBinLabel(1,"Muon");
    cor_TRIGPHYS -> GetXaxis() -> SetBinLabel(2,"EG");
    cor_TRIGPHYS -> GetXaxis() -> SetBinLabel(3,"Hadronic");
    cor_TRIGPHYS -> GetXaxis() -> SetBinLabel(4,"Muon+EG");
    cor_TRIGPHYS -> GetXaxis() -> SetBinLabel(5,"Muon+Hadronic");
    cor_TRIGPHYS -> GetXaxis() -> SetBinLabel(6,"EG+Hadronic");

    cor_TRIGPHYS -> GetYaxis() -> SetBinLabel(1,"Muon");
    cor_TRIGPHYS -> GetYaxis() -> SetBinLabel(2,"EG");
    cor_TRIGPHYS -> GetYaxis() -> SetBinLabel(3,"Hadronic");
    cor_TRIGPHYS -> GetYaxis() -> SetBinLabel(4,"Muon+EG");
    cor_TRIGPHYS -> GetYaxis() -> SetBinLabel(5,"Muon+Hadronic");
    cor_TRIGPHYS -> GetYaxis() -> SetBinLabel(6,"EG+Hadronic");

    cor_TRIGPHYS -> Draw("colz");
    cor_TRIGPHYS -> Draw("same,text"); 
    c14->Write();

    TCanvas* c15 = new TCanvas("c15","c15");
    c15 -> cd();
    h_TRIGPHYS_pure -> GetXaxis() -> SetBinLabel(1,"Muon");
    h_TRIGPHYS_pure -> GetXaxis() -> SetBinLabel(2,"EG");
    h_TRIGPHYS_pure -> GetXaxis() -> SetBinLabel(3,"Hadronic");
    h_TRIGPHYS_pure -> GetXaxis() -> SetBinLabel(4,"Muon+EG");
    h_TRIGPHYS_pure -> GetXaxis() -> SetBinLabel(5,"Muon+Hadronic");
    h_TRIGPHYS_pure -> GetXaxis() -> SetBinLabel(6,"EG+Hadronic");
		
    h_TRIGPHYS_pure -> SetYTitle("Pure rate (kHz)");
    h_TRIGPHYS_pure -> Draw();
    c15->Write();

    TCanvas* c16 = new TCanvas("c16","c16");
    c16 -> cd();
    h_TRIGPHYS_shared -> GetXaxis() -> SetBinLabel(1,"Muon");
    h_TRIGPHYS_shared -> GetXaxis() -> SetBinLabel(2,"EG");
    h_TRIGPHYS_shared -> GetXaxis() -> SetBinLabel(3,"Hadronic");
    h_TRIGPHYS_shared -> GetXaxis() -> SetBinLabel(4,"Muon+EG");
    h_TRIGPHYS_shared -> GetXaxis() -> SetBinLabel(5,"Muon+Hadronic");
    h_TRIGPHYS_shared -> GetXaxis() -> SetBinLabel(6,"EG+Hadronic");
	
    h_TRIGPHYS_shared -> SetYTitle("Shared rate (kHz)");
    h_TRIGPHYS_shared -> Draw();
    c16->Write();

    fOut.Close();
  }

  std::cout << "L1Bit" << "\t" << "L1SeedName" << "\t" << "pre-scale" << "\t" << "rate@13TeV" << "\t +/- \t" << "error_rate@13TeV" << "\t " << "pure@13TeV" << std::endl;
  TXTOutfile << "L1Bit" << "\t" << "L1SeedName" << "\t" << "pre-scale" << "\t" << "rate@13TeV" << "\t +/- \t" << "error_rate@13TeV" << "\t " << "pure@13TeV" << std::endl;

  if(writefiles) CSVOutfile << "L1Bit" << "\t" << "TriggerName" << "\t" << "Prescale" << std::endl; 

  Float_t totalrate = 0.;

  for(Int_t k=1; k < kOFFSET+1; k++) {  // -- kOFFSET now contains the number of triggers we have calculated
    TString name = h_All->GetXaxis()->GetBinLabel(k);

    Float_t rate = h_All->GetBinContent(k);
    Float_t err_rate  = h_All->GetBinError(k);
    Float_t pure = h_Pure->GetBinContent(k);

    std::string L1namest = (std::string)name;
    std::map<std::string, int>::const_iterator it = a.Prescales.find(L1namest);
    Float_t pre;
    if (it == a.Prescales.end() ) {
      std::cout << " --- SET P = 1 FOR SEED :  " << L1namest << std::endl;
      if (writefiles) { TXTOutfile << " --- SET P = 1 FOR SEED :  " << L1namest << std::endl; }
      pre = 1;
    }
    else {
      pre = it -> second;
    }
    Bool_t bias = a.Biased[L1namest];

    if (bias) std::cout << a.L1BitNumber(L1namest) << "\t" << name << "\t" << pre << "\t" << rate << "\t +/- \t" << err_rate << "\t " << pure << "\t" << " ***  BIAS  *** " << std::endl;
    else
      std::cout << a.L1BitNumber(L1namest) << "\t" << name << "\t" << pre << "\t" << rate << "\t +/- \t" << err_rate << "\t " << pure << std::endl;

    if (writefiles) { 
      if (bias) { TXTOutfile << a.L1BitNumber(L1namest) << "\t" << name << "\t" << pre << "\t" << rate << "\t +/- \t" << err_rate << "\t " << pure << "\t" << " ***  BIAS  *** " << std::endl; }
      else { TXTOutfile << a.L1BitNumber(L1namest) << "\t" << name << "\t" << pre << "\t" << rate << "\t +/- \t" << err_rate << "\t " << pure << std::endl; }

      CSVOutfile << a.L1BitNumber(L1namest) << ":\t" << name << ":\t" << pre << std::endl; 
    }

    totalrate +=rate;
  }

  std::cout << std::endl << "Total rate (without overlaps) = " << totalrate << std::endl;

  if (writefiles) {
    TXTOutfile << std::endl << a.GetPrintout() << std::endl;
    TXTOutfile << std::endl << "Total rate (without overlaps) = " << totalrate << std::endl;
  }

  TXTOutfile.close();
  CSVOutfile.close();
}

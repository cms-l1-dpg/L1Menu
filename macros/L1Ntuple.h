//////////////////////////////////////
// Example root macro for l1 ntuples
//////////////////////////////////////

#ifndef L1Ntuple_h
#define L1Ntuple_h

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TFriendElement.h>
#include <TList.h>
#include <TMatrix.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TCanvas.h>

#include "L1Trigger/L1TNtuples/interface/L1AnalysisEventDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"

class L1Ntuple {
public:
  TChain          *fChain;   //!pointer to the analyzed TTree or TChain
  TChain          *fEvent;
  TChain          *ftreeEvent;
  TChain          *ftreemuon;
  TChain          *ftreereco;
  TChain          *ftreeExtra;
  TChain          *ftreeMenu;
  TChain          *ftreeEmuExtra;
  Int_t            fCurrent; //!current Tree number in a TChain

  bool doreco;
  bool doEvent;
  bool domuonreco;
  bool dol1extra;
  bool dol1emuextra;
  bool dol1menu;  

  L1Analysis::L1AnalysisEventDataFormat        *event_;
  L1Analysis::L1AnalysisL1UpgradeDataFormat    *upgrade_;

  L1Ntuple();
  L1Ntuple(const std::string & fname);

  virtual ~L1Ntuple();

  bool Open(const std::string & fname);
  bool OpenWithList(const std::string & fname);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init();
  //virtual void     Loop();
  void     Test();
  void     Test2();
  Long64_t GetEntries();

 private :
  bool CheckFirstFile();
  bool OpenWithoutInit();
  bool OpenNtupleList(const std::string & fname);

  std::vector<std::string> listNtuples;
  Long64_t nentries_;
  TFile* rf;
};

#endif

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec  2 12:38:27 2015 by ROOT version 6.02/13
// from TTree L1UpgradeTree/L1UpgradeTree
// found on file: l1t_ntuple_test.root
//////////////////////////////////////////////////////////
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"

class L1Ntuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
 //L1Analysis::L1AnalysisL1UpgradeDataFormat *L1Upgrade;
    unsigned short int nEGs;
    std::vector<float> egEt;
    std::vector<float> egEta;
    std::vector<float> egPhi;
    std::vector<short int> egIEt;
    std::vector<short int> egIEta;
    std::vector<short int> egIPhi;
    std::vector<short int>    egIso;
    std::vector<short int>    egBx;
 
    unsigned short int nTaus;
    std::vector<float> tauEt;
    std::vector<float> tauEta;
    std::vector<float> tauPhi;
    std::vector<short int> tauIEt;
    std::vector<short int> tauIEta;
    std::vector<short int> tauIPhi;
    std::vector<short int>    tauIso;
    std::vector<short int>    tauBx;

    unsigned short int nJets;
    std::vector<float> jetEt;
    std::vector<float> jetEta;
    std::vector<float> jetPhi;
    std::vector<short int> jetIEt;
    std::vector<short int> jetIEta;
    std::vector<short int> jetIPhi;
    std::vector<short int>    jetBx;

    unsigned short int nMuons;
    std::vector<float>   muonEt;
    std::vector<float>   muonEta;
    std::vector<float>   muonPhi;
    std::vector<short int>   muonIEt;
    std::vector<short int>   muonIEta;
    std::vector<short int>   muonIPhi;
    std::vector<short int>      muonChg;
    std::vector<unsigned short int> muonIso;
    std::vector<unsigned short int> muonQual;
    std::vector<short int>      muonBx;

    
    unsigned short int nSums;
    std::vector<short int> sumType;
    std::vector<float> sumEt;
    std::vector<float> sumPhi;
    std::vector<short int> sumIEt;
    std::vector<short int> sumIPhi;
    std::vector<float> sumBx;
   
   // List of branches
   TBranch        *b_L1Upgrade_nEGs;   //!
   TBranch        *b_L1Upgrade_egEt;   //!
   TBranch        *b_L1Upgrade_egEta;   //!
   TBranch        *b_L1Upgrade_egPhi;   //!
   TBranch        *b_L1Upgrade_egIEt;   //!
   TBranch        *b_L1Upgrade_egIEta;   //!
   TBranch        *b_L1Upgrade_egIPhi;   //!
   TBranch        *b_L1Upgrade_egIso;   //!
   TBranch        *b_L1Upgrade_egBx;   //!
   TBranch        *b_L1Upgrade_nTaus;   //!
   TBranch        *b_L1Upgrade_tauEt;   //!
   TBranch        *b_L1Upgrade_tauEta;   //!
   TBranch        *b_L1Upgrade_tauPhi;   //!
   TBranch        *b_L1Upgrade_tauIEt;   //!
   TBranch        *b_L1Upgrade_tauIEta;   //!
   TBranch        *b_L1Upgrade_tauIPhi;   //!
   TBranch        *b_L1Upgrade_tauIso;   //!
   TBranch        *b_L1Upgrade_tauBx;   //!
   TBranch        *b_L1Upgrade_nJets;   //!
   TBranch        *b_L1Upgrade_jetEt;   //!
   TBranch        *b_L1Upgrade_jetEta;   //!
   TBranch        *b_L1Upgrade_jetPhi;   //!
   TBranch        *b_L1Upgrade_jetIEt;   //!
   TBranch        *b_L1Upgrade_jetIEta;   //!
   TBranch        *b_L1Upgrade_jetIPhi;   //!
   TBranch        *b_L1Upgrade_jetBx;   //!
   TBranch        *b_L1Upgrade_nMuons;   //!
   TBranch        *b_L1Upgrade_muonEt;   //!
   TBranch        *b_L1Upgrade_muonEta;   //!
   TBranch        *b_L1Upgrade_muonPhi;   //!
   TBranch        *b_L1Upgrade_muonIEt;   //!
   TBranch        *b_L1Upgrade_muonIEta;   //!
   TBranch        *b_L1Upgrade_muonIPhi;   //!
   TBranch        *b_L1Upgrade_muonChg;   //!
   TBranch        *b_L1Upgrade_muonIso;   //!
   TBranch        *b_L1Upgrade_muonBx;   //!
   TBranch        *b_L1Upgrade_nSums;   //!
   TBranch        *b_L1Upgrade_sumEt;   //!
   TBranch        *b_L1Upgrade_sumPhi;   //!
   TBranch        *b_L1Upgrade_sumIEt;   //!
   TBranch        *b_L1Upgrade_sumIPhi;   //!
   TBranch        *b_L1Upgrade_sumBx;   //!

   L1Ntuple(TTree *tree=0);
   virtual ~L1Ntuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

L1Ntuple::L1Ntuple(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("l1t_ntuple_test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("l1t_ntuple_test.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("l1t_ntuple_test.root:/l1UpgradeTree");
      dir->GetObject("L1UpgradeTree",tree);

   }
   Init(tree);
}

L1Ntuple::~L1Ntuple()
{
   if (!fChain) return;
   //delete fChain->GetCurrentFile();
}

Int_t L1Ntuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t L1Ntuple::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void L1Ntuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   
   fChain->SetBranchAddress("nEGs", &nEGs, &b_L1Upgrade_nEGs);
   fChain->SetBranchAddress("egEt", &egEt, &b_L1Upgrade_egEt);
   fChain->SetBranchAddress("egEta", &egEta, &b_L1Upgrade_egEta);
   fChain->SetBranchAddress("egPhi", &egPhi, &b_L1Upgrade_egPhi);
   fChain->SetBranchAddress("egIEt", &egIEt, &b_L1Upgrade_egIEt);
   fChain->SetBranchAddress("egIEta", &egIEta, &b_L1Upgrade_egIEta);
   fChain->SetBranchAddress("egIPhi", &egIPhi, &b_L1Upgrade_egIPhi);
   fChain->SetBranchAddress("egIso", &egIso, &b_L1Upgrade_egIso);
   fChain->SetBranchAddress("egBx", &egBx, &b_L1Upgrade_egBx);
   fChain->SetBranchAddress("nTaus", &nTaus, &b_L1Upgrade_nTaus);
   fChain->SetBranchAddress("tauEt", &tauEt, &b_L1Upgrade_tauEt);
   fChain->SetBranchAddress("tauEta", &tauEta, &b_L1Upgrade_tauEta);
   fChain->SetBranchAddress("tauPhi", &tauPhi, &b_L1Upgrade_tauPhi);
   fChain->SetBranchAddress("tauIEt", &tauIEt, &b_L1Upgrade_tauIEt);
   fChain->SetBranchAddress("tauIEta", &tauIEta, &b_L1Upgrade_tauIEta);
   fChain->SetBranchAddress("tauIPhi", &tauIPhi, &b_L1Upgrade_tauIPhi);
   fChain->SetBranchAddress("tauIso", &tauIso, &b_L1Upgrade_tauIso);
   fChain->SetBranchAddress("tauBx", &tauBx, &b_L1Upgrade_tauBx);
   fChain->SetBranchAddress("nJets", &nJets, &b_L1Upgrade_nJets);
   fChain->SetBranchAddress("jetEt", &jetEt, &b_L1Upgrade_jetEt);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_L1Upgrade_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_L1Upgrade_jetPhi);
   fChain->SetBranchAddress("jetIEt", &jetIEt, &b_L1Upgrade_jetIEt);
   fChain->SetBranchAddress("jetIEta", &jetIEta, &b_L1Upgrade_jetIEta);
   fChain->SetBranchAddress("jetIPhi", &jetIPhi, &b_L1Upgrade_jetIPhi);
   fChain->SetBranchAddress("jetBx", &jetBx, &b_L1Upgrade_jetBx);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_L1Upgrade_nMuons);
   fChain->SetBranchAddress("muonEt", &muonEt, &b_L1Upgrade_muonEt);
   fChain->SetBranchAddress("muonEta", &muonEta, &b_L1Upgrade_muonEta);
   fChain->SetBranchAddress("muonPhi", &muonPhi, &b_L1Upgrade_muonPhi);
   fChain->SetBranchAddress("muonIEt", &muonIEt, &b_L1Upgrade_muonIEt);
   fChain->SetBranchAddress("muonIEta", &muonIEta, &b_L1Upgrade_muonIEta);
   fChain->SetBranchAddress("muonIPhi", &muonIPhi, &b_L1Upgrade_muonIPhi);
   fChain->SetBranchAddress("muonChg", &muonChg, &b_L1Upgrade_muonChg);
   fChain->SetBranchAddress("muonIso", &muonIso, &b_L1Upgrade_muonIso);
   fChain->SetBranchAddress("muonBx", &muonBx, &b_L1Upgrade_muonBx);
   fChain->SetBranchAddress("nSums", &nSums, &b_L1Upgrade_nSums);
   fChain->SetBranchAddress("sumEt", &sumEt, &b_L1Upgrade_sumEt);
   fChain->SetBranchAddress("sumPhi", &sumPhi, &b_L1Upgrade_sumPhi);
   fChain->SetBranchAddress("sumIEt", &sumIEt, &b_L1Upgrade_sumIEt);
   fChain->SetBranchAddress("sumIPhi", &sumIPhi, &b_L1Upgrade_sumIPhi);
   fChain->SetBranchAddress("sumBx", &sumBx, &b_L1Upgrade_sumBx);
   Notify();
}

Bool_t L1Ntuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void L1Ntuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t L1Ntuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

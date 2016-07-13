#include "L1Ntuple.h"
#include "L1AlgoFactory.h"
#include <algorithm>
#include<map>
#include<iostream>

#include "TH1F.h"
#include "TH2F.h"

class getDistributions : public L1AlgoFactory
{
public :
  
  //constructor    
  getDistributions(std::string filename){
    if (filename.find(".root") != std::string::npos) {
      std::cout << "Reading RootFile: " << filename << std::endl;
      L1Ntuple::Open(filename);
    }else{
      std::cout << "Reading Filelist: " << filename << std::endl;
      if (! L1Ntuple::OpenWithList(filename)) exit(0);
    }
  }
  ~getDistributions() {}
  
  void run(std::string outfile, int nEvents);

private:
  std::map<std::string,TH1F*> hTH1F;
};

void getDistributions::run(std::string outfilename, int nEvents)
{

  TFile *outFile = new TFile((outfilename).c_str(),"recreate");
  outFile->cd();
  
  //L1Ntuple::OpenWithList(filename);

  hTH1F["nVtex"]   = new TH1F("nVtex","Number of Vertex",51,0,50);

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

    hTH1F["nVtex"]->Fill(recoVtx_->nVtx);

  }

  outFile->Write();
}
    
void goGetDistributions(int nEvents)
{
  std::string filename;
  filename="/afs/cern.ch/user/z/zhangj/test/L1Menu/80x_stage2/v62p3_MC2016/CMSSW_8_0_9/src/L1TriggerDPG/L1Menu/macros/r274199_809_zerobias_l1tint_v61p1.txt";

  std::string outfilename;
  outfilename="/afs/cern.ch/user/z/zhangj/test/L1Menu/80x_stage2/v62p3_MC2016/CMSSW_8_0_9/src/L1TriggerDPG/L1Menu/macros/myDataDistributions_274199.root";

  getDistributions getDistributions(filename);
  getDistributions.run(outfilename, nEvents);
}

    
  

    

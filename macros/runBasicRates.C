{
  std::cout << "L1Ntuple library loading ..." <<std::endl;
  gROOT->ProcessLine(".L L1Ntuple.C+");
  gROOT->ProcessLine(".L BasicRatePlots.C+");
  // gROOT->ProcessLine(".L L1AlgoFactory.h+");

  goRatePlots("RUN256843_Stage2",0,0)
}

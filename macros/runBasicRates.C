{
  std::cout << "L1Ntuple library loading ..." <<std::endl;
  gROOT->ProcessLine(".L L1Ntuple.C+");
  gROOT->ProcessLine(".L BasicRatePlots.C+");
  // gROOT->ProcessLine(".L L1AlgoFactory.h+");

  goRatePlots("RUN256843_Stage2",0,0)

  //("RUN256843_Stage2" is the sample identifier, 0=plot rates (1=plot cross sections), 0=number of events to process (0=all))
}

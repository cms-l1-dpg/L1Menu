{

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  //TFile fIn("results/results_13TEV_40PU_50bx_2012GCT10GEV_RE-EMUL_RATE.root");
  TFile fIn("results/results_13TEV_40PU_2015_RE-EMUL_RATE.root");
  fIn.cd();

  TH2F *nMuPtVsPt         = (TH2F*)fIn.Get("nMuPtVsPt");
  TH2F *nIsoElePtVsPt     = (TH2F*)fIn.Get("nIsoEGPtVsPt");
  TH2F *nEGPtVsPt         = (TH2F*)fIn.Get("nEGPtVsPt");
  TH2F *nOniaMuPtVsPt     = (TH2F*)fIn.Get("nOniaMuPtVsPt");
  TH2F *nMuVsHTT          = (TH2F*)fIn.Get("nMuVsHTT");
  TH2F *nAsymDiJetVsPt    = (TH2F*)fIn.Get("nAsymDiJetVsPt");
  TH2F *nAsymDiCenJetVsPt = (TH2F*)fIn.Get("nAsymDiCenJetVsPt");
  TH2F *nMuVsEG           = (TH2F*)fIn.Get("nMuVsEG");
  TH2F *nEGIsoEGVsPt      = (TH2F*)fIn.Get("nEGIsoEGVsPt");

  nMuPtVsPt->GetXaxis()->SetRangeUser(6.5,20.);

  nEGPtVsPt->GetXaxis()->SetRangeUser(12.,30.);
  nEGPtVsPt->GetYaxis()->SetRangeUser(6.,30.);

  nOniaMuPtVsPt->GetXaxis()->SetRangeUser(3.,15.);

  nMuVsHTT->GetXaxis()->SetRangeUser(5.,15.);
  nMuVsHTT->GetYaxis()->SetRangeUser(60.,200.);

  nAsymDiJetVsPt->GetXaxis()->SetRangeUser(70.,120.);
  nAsymDiJetVsPt->GetYaxis()->SetRangeUser(40.,120.);

  nAsymDiCenJetVsPt->GetXaxis()->SetRangeUser(60.,120.);
  nAsymDiCenJetVsPt->GetYaxis()->SetRangeUser(30.,120.);

  nMuVsEG->GetXaxis()->SetRangeUser(3.,15.);
  nMuVsEG->GetYaxis()->SetRangeUser(16.,25.);

  TCanvas c1("c1","c",1200,600); c1.cd();
  nMuPtVsPt->Draw("COLZ");
  c1.SaveAs("results/comparePlots/nMuPtVsPt.gif");

  TCanvas c2; c2.cd();
  nIsoElePtVsPt->Draw("COLZ");
  c2.SaveAs("results/comparePlots/nIsoEGPtVsPt.gif");

  TCanvas c3; c3.cd();
  nEGPtVsPt->Draw("COLZ");
  c3.SaveAs("results/comparePlots/nEGPtVsPt.gif");

  TCanvas c4("c4","c",1200,600); c4.cd();
  nOniaMuPtVsPt->Draw("COLZ");
  c4.SaveAs("results/comparePlots/nOniaMuPtVsPt.gif");

  TCanvas c5; c5.cd();
  nMuVsHTT->Draw("COLZ");
  c5.SaveAs("results/comparePlots/nMuVsHTT.gif");

  TCanvas c6("c6","c",1200,600); c6.cd();
  nAsymDiJetVsPt->Draw("COLZ");
  c6.SaveAs("results/comparePlots/nAsymDiJetVsPt.gif");

  TCanvas c7("c7","c",1200,600); c7.cd();
  nAsymDiCenJetVsPt->Draw("COLZ");
  c7.SaveAs("results/comparePlots/nAsymDiCenJetVsPt.gif");

  TCanvas c8; c8.cd();
  nMuVsEG->Draw("COLZ");
  c8.SaveAs("results/comparePlots/nMuVsEG.gif");

  TCanvas c9("c9","c",1200,600); c9.cd();
  nEGIsoEGVsPt->Draw("COLZ");
  c9.SaveAs("results/comparePlots/nEGIsoEGVsPt.gif");


}

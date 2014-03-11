{

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  TFile fIn("results/results_13TEV_40PU_2015_RE-EMUL_RATE.root");
  fIn.cd();

  TH2F *nMuPtVsPt         = (TH2F*)fIn.Get("nMuPtVsPt");
  TH2F *nIsolElePtVsPt    = (TH2F*)fIn.Get("nIsolElePtVsPt");
  TH2F *nElePtVsPt        = (TH2F*)fIn.Get("nElePtVsPt");
  TH2F *nOniaMuPtVsPt     = (TH2F*)fIn.Get("nOniaMuPtVsPt");
  TH2F *nMuVsHTT          = (TH2F*)fIn.Get("nMuVsHTT");
  TH2F *nAsymDiJetVsPt    = (TH2F*)fIn.Get("nAsymDiJetVsPt");
  TH2F *nAsymDiCenJetVsPt = (TH2F*)fIn.Get("nAsymDiCenJetVsPt");

  nMuPtVsPt->GetXaxis()->SetRangeUser(5.,30.);

  nElePtVsPt->GetXaxis()->SetRangeUser(13.,30.);
  nElePtVsPt->GetYaxis()->SetRangeUser(3.,30.);

  nMuVsHTT->GetXaxis()->SetRangeUser(2.,15.);
  nMuVsHTT->GetYaxis()->SetRangeUser(60.,200.);

  nAsymDiJetVsPt->GetXaxis()->SetRangeUser(70.,100.);
  nAsymDiJetVsPt->GetYaxis()->SetRangeUser(50.,100.);

  nAsymDiCenJetVsPt->GetXaxis()->SetRangeUser(70.,100.);
  nAsymDiCenJetVsPt->GetYaxis()->SetRangeUser(50.,100.);

  TCanvas c1; c1.cd();
  nMuPtVsPt->Draw("COLZ");
  c1.SaveAs("results/comparePlots/nMuPtVsPt.gif");

  TCanvas c2; c2.cd();
  nIsolElePtVsPt->Draw("COLZ");
  c2.SaveAs("results/comparePlots/nIsolElePtVsPt.gif");

  TCanvas c3; c3.cd();
  nElePtVsPt->Draw("COLZ");
  c3.SaveAs("results/comparePlots/nElePtVsPt.gif");

  TCanvas c4; c4.cd();
  nOniaMuPtVsPt->Draw("COLZ");
  c4.SaveAs("results/comparePlots/nOniaMuPtVsPt.gif");

  TCanvas c5; c5.cd();
  nMuVsHTT->Draw("COLZ");
  c5.SaveAs("results/comparePlots/nMuVsHTT.gif");

  TCanvas c6; c6.cd();
  nAsymDiJetVsPt->Draw("COLZ");
  c6.SaveAs("results/comparePlots/nAsymDiJetVsPt.gif");

  TCanvas c7; c7.cd();
  nAsymDiCenJetVsPt->Draw("COLZ");
  c7.SaveAs("results/comparePlots/nAsymDiCenJetVsPt.gif");


}

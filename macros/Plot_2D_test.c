int Plot_2D_test()
{
TCanvas* mycanvas = new TCanvas("mycanvas", "mycanvas", 900, 700);
//TCanvas* mycanvas = new TCanvas();

//TString sp = "T2bwC_mStop500_mLSP490";
//TString sp = "T2bwC_mStop600_mLSP520";
//TString sp = "T2fbd_mStop600_mLSP520";
TString sp = "Prescale_2018_v1_0_0_Col_2.0_opt2_Pure_Rate_nanodst";

//TString folder = "Baseline_Only/";
TString folder = "Col0/";
//TString var = "mtb_mt2_h";
TString var = "cor_Seeds";
gStyle->SetOptStat(kFALSE);
gStyle->SetPaintTextFormat("4.2f");
TFile *f1 = new TFile("results/" + sp + ".root");
TH2D *h1 = (TH2D*)f1->Get(folder + var);

h1->SetTitle("v1_0_0_Col_2.0_opt2 seeds correlation");
//h1->GetXaxis()->SetRangeUser(240,800);
//h1->GetYaxis()->SetRangeUser(240,800);
h1->GetYaxis()->SetLabelSize(0.002);
h1->GetXaxis()->SetLabelSize(0.002);
//h1->GetYaxis()->SetTitle("ISR pt");
//h1->GetXaxis()->SetTitle("MET");
//h1->SetMarkerSize(2);
//h1->Draw("COLZTEXT");
h1->Draw("COLZ");

//TLegend* leg = new TLegend(0.7,0.9,0.89,0.95);
//leg->SetBorderSize(0);
//leg->SetTextFont(42);
//leg->SetFillColor(0);
//leg->AddEntry((TObject*)0,"36 fb^{-1} (13TeV)","");
//leg->Draw("SAME");

mycanvas->SetBottomMargin(0.1);
mycanvas->SetLeftMargin(0.1);
mycanvas->SaveAs(sp + "_" + var + ".pdf");

return 0;
}


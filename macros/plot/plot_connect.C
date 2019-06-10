int plot_connect()
{
	TH1F* h_8b4e = new TH1F("","",120, 0, 120);

	h_8b4e->SetBinContent( 31, 43.131283);
	h_8b4e->SetBinError( 31, 2.312079);
	h_8b4e->SetBinContent( 32, 45.428141);
	h_8b4e->SetBinError( 32, 1.938826);
	h_8b4e->SetBinContent( 33, 46.575148);
	h_8b4e->SetBinError( 33, 1.998731);
	h_8b4e->SetBinContent( 34, 49.815156);
	h_8b4e->SetBinError( 34, 2.167926);
	h_8b4e->SetBinContent( 35, 47.180764);
	h_8b4e->SetBinError( 35, 2.151258);
	h_8b4e->SetBinContent( 36, 51.249068);
	h_8b4e->SetBinError( 36, 2.247420);
	h_8b4e->SetBinContent( 37, 59.739853);
	h_8b4e->SetBinError( 37, 2.422771);
	h_8b4e->SetBinContent( 38, 55.781678);
	h_8b4e->SetBinError( 38, 2.411652);
	h_8b4e->SetBinContent( 39, 63.153634);
	h_8b4e->SetBinError( 39, 2.576090);
	h_8b4e->SetBinContent( 40, 54.079462);
	h_8b4e->SetBinError( 40, 2.440572);
	h_8b4e->SetBinContent( 41, 59.553251);
	h_8b4e->SetBinError( 41, 2.671328);
	h_8b4e->SetBinContent( 42, 64.944150);
	h_8b4e->SetBinError( 42, 2.715451);
	h_8b4e->SetBinContent( 43, 62.976726);
	h_8b4e->SetBinError( 43, 3.526020);
	h_8b4e->SetBinContent( 44, 63.812938);
	h_8b4e->SetBinError( 44, 5.029164);
	h_8b4e->SetBinContent( 45, 73.924206);
	h_8b4e->SetBinError( 45, 3.202013);
	h_8b4e->SetBinContent( 46, 70.039733);
	h_8b4e->SetBinError( 46, 3.071445);
	h_8b4e->SetBinContent( 47, 76.802345);
	h_8b4e->SetBinError( 47, 3.208464);
	h_8b4e->SetBinContent( 48, 83.564478);
	h_8b4e->SetBinError( 48, 3.408666);
	h_8b4e->SetBinContent( 49, 79.625499);
	h_8b4e->SetBinError( 49, 3.280912);
	h_8b4e->SetBinContent( 50, 80.645637);
	h_8b4e->SetBinError( 50, 3.278711);
	h_8b4e->SetBinContent( 51, 87.473108);
	h_8b4e->SetBinError( 51, 3.676769);
	h_8b4e->SetBinContent( 52, 89.371536);
	h_8b4e->SetBinError( 52, 3.508140);
	h_8b4e->SetBinContent( 53, 86.423603);
	h_8b4e->SetBinError( 53, 3.192129);
	h_8b4e->SetBinContent( 54, 90.583953);
	h_8b4e->SetBinError( 54, 4.022975);
	h_8b4e->SetBinContent( 55, 97.173881);
	h_8b4e->SetBinError( 55, 6.003422);


	h_8b4e->SetBinContent( 94, 713.084815);
	h_8b4e->SetBinError( 94, 14.644529);
	h_8b4e->SetBinContent( 95, 805.004373);
	h_8b4e->SetBinError( 95, 6.712800);
	h_8b4e->SetBinContent( 96, 822.245425);
	h_8b4e->SetBinError( 96, 6.050831);
	h_8b4e->SetBinContent( 97, 859.797896);
	h_8b4e->SetBinError( 97, 5.251268);
	h_8b4e->SetBinContent( 98, 905.713356);
	h_8b4e->SetBinError( 98, 4.766714);
	h_8b4e->SetBinContent( 99, 844.763789);
	h_8b4e->SetBinError( 99, 6.482284);
	h_8b4e->SetBinContent( 100, 1056.457283);
	h_8b4e->SetBinError( 100, 7.456128);
	h_8b4e->SetBinContent( 101, 1168.631958);
	h_8b4e->SetBinError( 101, 6.744964);
	h_8b4e->SetBinContent( 102, 1201.002198);
	h_8b4e->SetBinError( 102, 4.748046);
	h_8b4e->SetBinContent( 103, 1292.425598);
	h_8b4e->SetBinError( 103, 9.010225);
	h_8b4e->SetBinContent( 104, 1377.492603);
	h_8b4e->SetBinError( 104, 14.742030);
	h_8b4e->SetBinContent( 105, 1477.195725);
	h_8b4e->SetBinError( 105, 14.325511);
	h_8b4e->SetBinContent( 106, 1537.969760);
	h_8b4e->SetBinError( 106, 14.653313);
	h_8b4e->SetBinContent( 107, 1606.824437);
	h_8b4e->SetBinError( 107, 14.582209);
	h_8b4e->SetBinContent( 108, 1683.141110);
	h_8b4e->SetBinError( 108, 62.168235);

	TH1F* h_12b_all = new TH1F("","",120, 0, 120);

	h_12b_all->SetBinContent( 94, 532.087215);
	h_12b_all->SetBinError( 94, 9.875505);
	h_12b_all->SetBinContent( 95, 578.433330);
	h_12b_all->SetBinError( 95, 4.450409);
	h_12b_all->SetBinContent( 96, 591.326697);
	h_12b_all->SetBinError( 96, 4.009202);
	h_12b_all->SetBinContent( 97, 618.353910);
	h_12b_all->SetBinError( 97, 3.472091);
	h_12b_all->SetBinContent( 98, 645.300977);
	h_12b_all->SetBinError( 98, 3.138711);
	h_12b_all->SetBinContent( 99, 607.546697);
	h_12b_all->SetBinError( 99, 4.298799);
	h_12b_all->SetBinContent( 100, 749.024374);
	h_12b_all->SetBinError( 100, 4.924915);
	h_12b_all->SetBinContent( 101, 822.294294);
	h_12b_all->SetBinError( 101, 4.437321);
	h_12b_all->SetBinContent( 102, 840.385492);
	h_12b_all->SetBinError( 102, 3.108427);
	h_12b_all->SetBinContent( 103, 891.799773);
	h_12b_all->SetBinError( 103, 5.812145);
	h_12b_all->SetBinContent( 104, 940.013546);
	h_12b_all->SetBinError( 104, 9.449878);
	h_12b_all->SetBinContent( 105, 1002.184227);
	h_12b_all->SetBinError( 105, 9.148267);
	h_12b_all->SetBinContent( 106, 1041.640134);
	h_12b_all->SetBinError( 106, 9.344413);
	h_12b_all->SetBinContent( 107, 1087.057200);
	h_12b_all->SetBinError( 107, 9.287356);
	h_12b_all->SetBinContent( 108, 1106.372311);
	h_12b_all->SetBinError( 108, 39.018743);


	TH1F* h_12b_5to10 = new TH1F("","",120, 0, 120);

	h_12b_5to10->SetBinContent( 94, 344.527146);
	h_12b_5to10->SetBinError( 94, 11.673857);
	h_12b_5to10->SetBinContent( 95, 327.448410);
	h_12b_5to10->SetBinError( 95, 4.923062);
	h_12b_5to10->SetBinContent( 96, 335.441923);
	h_12b_5to10->SetBinError( 96, 4.436814);
	h_12b_5to10->SetBinContent( 97, 351.017739);
	h_12b_5to10->SetBinError( 97, 3.838839);
	h_12b_5to10->SetBinContent( 98, 354.673923);
	h_12b_5to10->SetBinError( 98, 3.415382);
	h_12b_5to10->SetBinContent( 99, 349.420536);
	h_12b_5to10->SetBinError( 99, 4.790630);
	h_12b_5to10->SetBinContent( 100, 394.174540);
	h_12b_5to10->SetBinError( 100, 5.258934);
	h_12b_5to10->SetBinContent( 101, 418.052992);
	h_12b_5to10->SetBinError( 101, 4.656256);
	h_12b_5to10->SetBinContent( 102, 420.491597);
	h_12b_5to10->SetBinError( 102, 3.232066);
	h_12b_5to10->SetBinContent( 103, 443.988445);
	h_12b_5to10->SetBinError( 103, 5.998748);
	h_12b_5to10->SetBinContent( 104, 453.202786);
	h_12b_5to10->SetBinError( 104, 9.590654);
	h_12b_5to10->SetBinContent( 105, 474.250965);
	h_12b_5to10->SetBinError( 105, 9.191877);
	h_12b_5to10->SetBinContent( 106, 487.885771);
	h_12b_5to10->SetBinError( 106, 9.339344);
	h_12b_5to10->SetBinContent( 107, 505.511484);
	h_12b_5to10->SetBinError( 107, 9.243210);
	h_12b_5to10->SetBinContent( 108, 503.758671);
	h_12b_5to10->SetBinError( 108, 38.411239);

	TH1F* h_48b = new TH1F("","",120, 0, 120);

	h_48b->SetBinContent( 31, 42.917117);
	h_48b->SetBinError( 31, 0.941926);
	h_48b->SetBinContent( 32, 44.464358);
	h_48b->SetBinError( 32, 0.783217);
	h_48b->SetBinContent( 33, 44.630384);
	h_48b->SetBinError( 33, 0.799268);
	h_48b->SetBinContent( 34, 48.485340);
	h_48b->SetBinError( 34, 0.873645);
	h_48b->SetBinContent( 35, 50.027837);
	h_48b->SetBinError( 35, 0.904676);
	h_48b->SetBinContent( 36, 50.568133);
	h_48b->SetBinError( 36, 0.911767);
	h_48b->SetBinContent( 37, 54.731730);
	h_48b->SetBinError( 37, 0.947318);
	h_48b->SetBinContent( 38, 55.240470);
	h_48b->SetBinError( 38, 0.980205);
	h_48b->SetBinContent( 39, 58.210273);
	h_48b->SetBinError( 39, 1.010254);
	h_48b->SetBinContent( 40, 56.701056);
	h_48b->SetBinError( 40, 1.020688);
	h_48b->SetBinContent( 41, 60.453473);
	h_48b->SetBinError( 41, 1.099336);
	h_48b->SetBinContent( 42, 60.675852);
	h_48b->SetBinError( 42, 1.072105);
	h_48b->SetBinContent( 43, 64.520475);
	h_48b->SetBinError( 43, 1.457369);
	h_48b->SetBinContent( 44, 68.484662);
	h_48b->SetBinError( 44, 2.127715);
	h_48b->SetBinContent( 45, 70.816454);
	h_48b->SetBinError( 45, 1.279978);
	h_48b->SetBinContent( 46, 69.019308);
	h_48b->SetBinError( 46, 1.245463);
	h_48b->SetBinContent( 47, 72.627415);
	h_48b->SetBinError( 47, 1.274361);
	h_48b->SetBinContent( 48, 75.381889);
	h_48b->SetBinError( 48, 1.322082);
	h_48b->SetBinContent( 49, 76.325446);
	h_48b->SetBinError( 49, 1.311868);
	h_48b->SetBinContent( 50, 78.606493);
	h_48b->SetBinError( 50, 1.321911);
	h_48b->SetBinContent( 51, 81.085180);
	h_48b->SetBinError( 51, 1.446106);
	h_48b->SetBinContent( 52, 83.037298);
	h_48b->SetBinError( 52, 1.380890);
	h_48b->SetBinContent( 53, 84.329279);
	h_48b->SetBinError( 53, 1.287808);
	h_48b->SetBinContent( 54, 89.186110);
	h_48b->SetBinError( 54, 1.630211);
	h_48b->SetBinContent( 55, 92.090624);
	h_48b->SetBinError( 55, 2.387338);


	gStyle->SetOptStat(0);
	gStyle->SetPadLeftMargin(0.15);


	TCanvas* myCanvas = new TCanvas("myCanvas","myCanvas", 600, 600);
	myCanvas->SetGrid();
	myCanvas->SetLogy();

	TLegend* leg = new TLegend(0.16,0.6,0.6,0.89);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);


	//h_8b4e->GetYaxis()->SetLimits(0,2000);
	h_8b4e->GetXaxis()->SetTitle("PileUp");
	h_8b4e->GetYaxis()->SetRangeUser(1,1000000);
	h_8b4e->GetYaxis()->SetTitle("Rate (nBunches = 2544) [kHz]");
	h_8b4e->GetYaxis()->SetTitleOffset(1.8);

	h_8b4e->SetLineColor(kRed);
	h_8b4e->SetMarkerStyle(20);
	h_8b4e->SetMarkerColor(kRed);

	h_8b4e->Draw("pe");

	//TString fName = "expo";
	TF1 *my_func = new TF1("my_func","exp([0]*x + [1]) + [2]*x + [3]",0,120);

	//h_8b4e->Fit(fName, "0", "");
	h_8b4e->Fit(my_func, "0", "");
	//TF1* fit_8b4e = (TF1*)h_8b4e->GetFunction(fName)->Clone();
	TF1* fit_8b4e = (TF1*)h_8b4e->GetFunction("my_func")->Clone();
	fit_8b4e->SetLineColor(kRed);
	fit_8b4e->SetLineStyle(5);
	fit_8b4e->Draw("same");

	char func_8b4e[100];
	//sprintf (func_8b4e, "y = exp(%.2f + %.3fx)", fit_8b4e->GetParameter(0), fit_8b4e->GetParameter(1));
	sprintf (func_8b4e, "y = exp(%.3fx + %.2f) + %.2fx + %.1f", fit_8b4e->GetParameter(0), fit_8b4e->GetParameter(1), fit_8b4e->GetParameter(2), fit_8b4e->GetParameter(3));
	char chi2_8b4e[100];
	sprintf (chi2_8b4e, "8b4e (chi2/NDF = %.2f)", fit_8b4e->GetChisquare() / fit_8b4e->GetNDF());

	leg->AddEntry(h_8b4e,chi2_8b4e,"lep");
	leg->AddEntry((TObject*)0, func_8b4e, "");

	h_48b->Sumw2();
	h_48b->Add(h_12b_all, 12./48.);
	h_48b->Add(h_12b_5to10, 36./48.);

	h_48b->SetLineColor(kBlue);
	h_48b->SetMarkerStyle(22);
	h_48b->SetMarkerColor(kBlue);

	h_48b->Draw("pesame");

	//h_48b->Fit(fName, "0", "");
	h_48b->Fit(my_func, "0", "");
	//TF1* fit_48b = (TF1*)h_48b->GetFunction(fName)->Clone();
	TF1* fit_48b = (TF1*)h_48b->GetFunction("my_func")->Clone();
	fit_48b->SetLineColor(kBlue);
	fit_48b->SetLineStyle(5);
	fit_48b->Draw("same");

	char func_48b[100];
	//sprintf (func_48b, "y = exp(%.2f + %.3fx)", fit_48b->GetParameter(0), fit_48b->GetParameter(1));
	sprintf (func_48b, "y = exp(%.3fx + %.2f) + %.2fx + %.1f", fit_48b->GetParameter(0), fit_48b->GetParameter(1), fit_48b->GetParameter(2), fit_48b->GetParameter(3));
	char chi2_48b[100];
	sprintf (chi2_48b, "48b (chi2/NDF = %.2f)", fit_48b->GetChisquare() / fit_48b->GetNDF());

	leg->AddEntry(h_48b,chi2_48b,"lep");
	leg->AddEntry((TObject*)0, func_48b, "");

	leg->Draw("same");

	TLatex latex;
	latex.SetTextSize(0.04);
	latex.SetNDC();
	//latex.DrawLatex(0.5,ymax+0.4,"#bf{CMS} Preliminary, 2017 data");
	latex.DrawLatex(0.15,0.91,"CMS #bf{Preliminary}");
	//TString lumi_and_energy = "#bf{" + std::to_string(lumi) + " fb^{-1} (13TeV)}";
	latex.SetTextAlign(31);  //align at right bottom
	latex.DrawLatex(0.9,0.91,"#bf{2018 data (13TeV)}");

	//gPad->SetLogy(1);
	//myCanvas->SetLogy();

	myCanvas->SaveAs("8b4e_vs_48b_logY.png");

	return 0;
}

// dPhi / dEta rectangle definition 
// for GMT - RECO matching
#define MAX_MU_GMT_DPHI .2
#define MAX_MU_GMT_DETA .5

const float MuPtBins[26] = {0., 3., 4., 5., 6., 7., 8., 10., 12., 14., 16., 18., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 90., 100., 120., 140.};

#include "L1Ntuple.h"
#include "TriggeredMuon.C"

#include "TROOT.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"

class analyze_L1 : public L1Ntuple {

public :

  analyze_L1(Bool_t isReco)
  {
    _isReco = isReco;
  }
  ~analyze_L1() {}

  void Loop();
  float binRecoMuPt(float theMuPt);

private :

  Bool_t _isReco;

};

float analyze_L1::binRecoMuPt(float theMuPt){
  float binnedpt = -999.;

  if(theMuPt < 8.) binnedpt = trunc(theMuPt);
  else if(theMuPt < 10.) binnedpt = 8.;
  else if(theMuPt < 12.) binnedpt = 10.;
  else if(theMuPt < 14.) binnedpt = 12.;
  else if(theMuPt < 16.) binnedpt = 14.;
  else if(theMuPt < 18.) binnedpt = 16.;
  else if(theMuPt < 20.) binnedpt = 18.;
  else if(theMuPt < 25.) binnedpt = 20.;
  else if(theMuPt < 30.) binnedpt = 25.;
  else if(theMuPt < 35.) binnedpt = 30.;
  else if(theMuPt < 40.) binnedpt = 35.;
  else if(theMuPt < 45.) binnedpt = 40.;
  else if(theMuPt < 50.) binnedpt = 45.;
  else if(theMuPt < 60.) binnedpt = 50.;
  else if(theMuPt < 70.) binnedpt = 60.;
  else if(theMuPt < 80.) binnedpt = 70.;
  else if(theMuPt < 90.) binnedpt = 80.;
  else if(theMuPt < 100.) binnedpt = 90.;
  else if(theMuPt < 120.) binnedpt = 100.;
  else if(theMuPt < 140.) binnedpt = 120.;
  else binnedpt = 140.;

  return binnedpt;
}


void analyze_L1::Loop(){

  Int_t nevents = GetEntries();
  std::cout << "Running on " << nevents << " events." << std::endl;
  if(nevents > 2000000) nevents = 2000000;

  TH1F hTT11bx("hTT11bx","hTT11bx",7,-3.5,3.5);

  TH1F hDTTFbx("hDTTFbx","hDTTFbx",7,-3.5,3.5);
  TH1F hCSCTFbx("hCSCTFbx","hCSCTFbx",7,-3.5,3.5);
  TH1F hRPCbbx("hRPCbbx","hRPCbbx",7,-3.5,3.5);
  TH1F hRPCfbx("hRPCfbx","hRPCfbx",7,-3.5,3.5);
  TH1F hGMTbx("hGMTbx","hGMTbx",7,-3.5,3.5);

  TH1F hGMTpt("hGMTpt","hGMTpt",100,0.,100.);
  TH1F hGMTeta("hGMTeta","hGMTeta",30,-2.5,2.5);
  TH1F hGMTphi("hGMTphi","hGMTphi",30,0.,6.3);

  TH1F hGMT2mubx("hGMT2mubx","hGMT2mubx",7,-3.5,3.5);
  TH1F hGMT2mupt("hGMT2mupt","hGMT2mupt",100,0.,100.);
  TH1F hGMT2muphi("hGMT2muphi","hGMT2muphi",30,0.,6.3);
  TH1F hGMT2mueta("hGMT2mueta","hGMT2mueta",30,-2.5,2.5);
  TH1F hGMT2mu1BX("hGMT2mu1BX","hGMT2mu1BX",7,-3.5,3.5);

  TH2F hDTRPCBX("hDTRPCBX","hDTRPCBX",7,-3.5,3.5,7,-3.5,3.5);
  TH2F hCSCRPCBX("hCSCRPCBX","hCSCRPCBX",7,-3.5,3.5,7,-3.5,3.5);
  TH2F hDTCSCBX("hDTCSCBX","hDTCSCBX",7,-3.5,3.5,7,-3.5,3.5);

  TH1F hpTresidual("hpTresidual","hpTresidual",100.,-10.,10.);
  TProfile hpTresidual_r("hpTresidual_r","hpTresidual_r",30.,0.,50.);
  TProfile hpTresidual_z("hpTresidual_z","hpTresidual_z",30.,-50.,50.);
  TProfile hpTresidual_eta("hpTresidual_eta","hpTresidual_eta",30.,-2.5,2.5);

  TH1F hRecoMupT("hRecoMupT","hRecoMupT",100.,0.,100.);
  TH1F hRecoMuBinpT("hRecoMuBinpT","hRecoMuBinpT",25,MuPtBins);
  TH1F hGMTMuBinpT("hGMTMuBinpT","hGMTMuBinpT",25,MuPtBins);
  TH2F hRecoMuGMTpT("hRecoMuGMTpT","hRecoMuGMTpT",25,MuPtBins,25,MuPtBins);

  for (Long64_t event=0; event<nevents; event++)
    {     
      Long64_t ientry = LoadTree(event); if (ientry < 0) break;
      GetEntry(event);

      if (event%200000 == 0) {
	std::cout << "Processed " << event << " events." << std::endl;
      }

      //select events by technical trigger bit 11 - HO
      //Int_t Bit11 = (gt_->tt[2]>>11)&1;
      //if(Bit11 != 1) continue;

      Int_t Ndt = gmt_->Ndt;
      Int_t Ncsc = gmt_->Ncsc;
      Int_t Nrpcb = gmt_->Nrpcb;
      Int_t Nrpcf = gmt_->Nrpcf;
      Int_t Nmu = gmt_->N;

      if(_isReco){
	for(int iMu=0;iMu<recoMuon_->nMuons;iMu++){
	  if(recoMuon_->type.at(iMu) != 0) continue;
	  if(recoMuon_->pt[iMu] < 3.) continue;
	  bool isDisplaced = false;
	  if( fabs(recoMuon_->tr_imp_point_z[iMu]) > 20. || sqrt( pow(recoMuon_->tr_imp_point_x[iMu],2.) + pow(recoMuon_->tr_imp_point_y[iMu],2) ) > 10.) isDisplaced = true;

	  for (int iGMTmu=0; iGMTmu < Nmu; iGMTmu++){
	    TriggeredMuon trigMuons(recoMuon_,iMu,gmt_,iGMTmu);
	    if(!trigMuons.hasTriggerMatch() || !trigMuons.hasGmtBX0()) continue;

	    float iMu_binnedpt = binRecoMuPt(recoMuon_->pt[iMu]);
	    float theresidual = (iMu_binnedpt - gmt_->Pt[iGMTmu]) / iMu_binnedpt;

	    if(!isDisplaced) hpTresidual.Fill(theresidual);
	    hpTresidual_r.Fill(sqrt( pow(recoMuon_->tr_imp_point_x[iMu],2.) + pow(recoMuon_->tr_imp_point_y[iMu],2) ), theresidual);
	    hpTresidual_z.Fill(recoMuon_->tr_imp_point_z[iMu] , theresidual);
	    if(!isDisplaced) hpTresidual_eta.Fill(recoMuon_->eta[iMu] , theresidual);

	    hRecoMupT.Fill(recoMuon_->pt[iMu]);
	    if(iMu_binnedpt < 140. && gmt_->Pt[iGMTmu] < 140. && !isDisplaced){
	      hRecoMuBinpT.Fill(iMu_binnedpt);
	      hGMTMuBinpT.Fill(gmt_->Pt[iGMTmu]);
	      hRecoMuGMTpT.Fill(iMu_binnedpt,gmt_->Pt[iGMTmu]);
	    }

	    if( fabs((iMu_binnedpt - gmt_->Pt[iGMTmu]) / iMu_binnedpt) > 3. && !isDisplaced){
	      cout << "###################" << endl;
	      cout << "GMT pt = " << gmt_->Pt[iGMTmu] << endl;
	      cout << "GMT quality = " << gmt_->Qual[iGMTmu] << endl;
	      cout << "mu pt = " << recoMuon_->pt[iMu] << endl;
	      cout << "binned mu pt = " << iMu_binnedpt << endl;
	      cout << "Mu phi = " << recoMuon_->phi[iMu] << endl;
	      cout << "Mu eta = " << recoMuon_->eta[iMu] << endl;
	      cout << "Mu charge = " << recoMuon_->ch[iMu] << endl;
	      cout << "Mu x = " << recoMuon_->tr_imp_point_x[iMu] << endl;
	      cout << "Mu y = " << recoMuon_->tr_imp_point_y[iMu] << endl;
	      cout << "Mu z = " << recoMuon_->tr_imp_point_z[iMu] << endl;
	      cout << "Mu ch2 = " << recoMuon_->tr_normchi2[iMu] << endl;
	      cout << "Mu ntrhits = " << recoMuon_->tr_validhits[iMu] << endl;
	      cout << "Mu npixhits = " << recoMuon_->tr_validpixhits[iMu] << endl;
	      cout << "N reco Mu = " << recoMuon_->nMuons << endl;
	      cout << "N GMT Mu = " << Nmu << endl;
	      cout << "Run number = " << event_->run << endl;
	      cout << "Event number = " << event_->event << endl;
	      cout << "###################" << endl;
	    }


	  }
	}
      }

      for(int i11=0;i11<5;i11++){
	if( (gt_->tt[i11-2]>>11)&1 ) hTT11bx.Fill(i11-2.);
      }

      for (int idt=0; idt < Ndt; idt++){
	int bx = gmt_ -> Bxdt[idt];
	hDTTFbx.Fill(bx);
      }
 
      for (int icsc=0; icsc < Ncsc; icsc++){
	int bx = gmt_ -> Bxcsc[icsc];
	hCSCTFbx.Fill(bx);
      }

      for (int irpcb=0; irpcb < Nrpcb; irpcb++){
	int bx = gmt_ -> Bxrpcb[irpcb];
	hRPCbbx.Fill(bx);
      }

      for (int irpcf=0; irpcf < Nrpcf; irpcf++){
	int bx = gmt_ -> Bxrpcf[irpcf];
	hRPCfbx.Fill(bx);
      }
 
      for (int imu=0; imu < Nmu; imu++){
	int bx      = gmt_->CandBx[imu];
	float ptmu  = gmt_->Pt[imu];
	float etamu = gmt_->Eta[imu];
	float phimu = gmt_->Phi[imu];

	hGMTbx.Fill(bx);
	hGMTpt.Fill(ptmu);
	hGMTeta.Fill(etamu);
	hGMTphi.Fill(phimu);
      }

      if(Nmu > 1 && fabs((gmt_->Pt[0]-gmt_->Pt[1])/gmt_->Pt[0]) < 0.05 && fabs((gmt_->Eta[0]-gmt_->Eta[1])/gmt_->Eta[0]) < 0.1 ){
	hGMT2mubx.Fill(gmt_->CandBx[0] - gmt_->CandBx[1]);
	hGMT2mupt.Fill(gmt_->Pt[0]);
	hGMT2muphi.Fill(gmt_->Phi[0]);
	hGMT2mueta.Fill(gmt_->Eta[0]);
	if(gmt_->CandBx[0] <= gmt_->CandBx[1]) hGMT2mu1BX.Fill(gmt_->CandBx[0]);
	else if(gmt_->CandBx[0] > gmt_->CandBx[1]) hGMT2mu1BX.Fill(gmt_->CandBx[1]);
      }

      if(Ndt > 0 && Nrpcb > 0) hDTRPCBX.Fill(gmt_->Bxdt[0],gmt_->Bxrpcb[0]);
      if(Ncsc > 0 && Nrpcf > 0) hCSCRPCBX.Fill(gmt_->Bxcsc[0],gmt_->Bxrpcf[0]);
      if(Ndt > 0 && Ncsc > 0) hDTCSCBX.Fill(gmt_->Bxdt[0],gmt_->Bxcsc[0]);
    }

  //Aestethics
  hTT11bx.SetTitle("");

  hDTTFbx.SetTitle("");
  hCSCTFbx.SetTitle("");
  hRPCbbx.SetTitle("");
  hRPCfbx.SetTitle("");
  hGMTbx.SetTitle("");

  hGMTpt.SetTitle("");
  hGMTeta.SetTitle("");
  hGMTphi.SetTitle("");

  hGMT2mubx.SetTitle("");
  hGMT2mupt.SetTitle("");
  hGMT2muphi.SetTitle("");
  hGMT2mueta.SetTitle("");
  hGMT2mu1BX.SetTitle("");

  hDTRPCBX.SetTitle("");
  hCSCRPCBX.SetTitle("");
  hDTCSCBX.SetTitle("");

  hpTresidual.SetTitle("");
  hpTresidual_r.SetTitle("");
  hpTresidual_z.SetTitle("");
  hpTresidual_eta.SetTitle("");
  hRecoMupT.SetTitle("");
  hRecoMuBinpT.SetTitle("");
  hGMTMuBinpT.SetTitle("");
  hRecoMuGMTpT.SetTitle("");

  hTT11bx.GetXaxis()->SetTitle("TT11 BX");

  hDTTFbx.GetXaxis()->SetTitle("DTTF BX");
  hCSCTFbx.GetXaxis()->SetTitle("CSCTF BX");
  hRPCbbx.GetXaxis()->SetTitle("RPCb BX");
  hRPCfbx.GetXaxis()->SetTitle("RPCf BX");
  hGMTbx.GetXaxis()->SetTitle("GMT BX");

  hGMTpt.GetXaxis()->SetTitle("L1 Muon #p_{T} [GeV]");
  hGMTeta.GetXaxis()->SetTitle("L1 Muon #eta");
  hGMTphi.GetXaxis()->SetTitle("L1 Muon #phi");

  hGMT2mubx.GetXaxis()->SetTitle("#Delta_{BX} for identical muons");
  hGMT2mupt.GetXaxis()->SetTitle("p_{T} for identical muons");
  hGMT2muphi.GetXaxis()->SetTitle("#phi for identical muons");
  hGMT2mueta.GetXaxis()->SetTitle("#eta for identical muons");
  hGMT2mu1BX.GetXaxis()->SetTitle("BX for first muon for identical muons");

  hDTRPCBX.GetXaxis()->SetTitle("DT BX");
  hCSCRPCBX.GetXaxis()->SetTitle("CSC BX");
  hDTCSCBX.GetXaxis()->SetTitle("DT BX");

  hpTresidual.GetXaxis()->SetTitle("#Delta p_{T} / p_{T}");
  hpTresidual_r.GetXaxis()->SetTitle("Radius on transverse plane (cm)");
  hpTresidual_z.GetXaxis()->SetTitle("Distance on z (cm)");
  hpTresidual_eta.GetXaxis()->SetTitle("#eta");
  hRecoMupT.GetXaxis()->SetTitle("Reco Muon p_{T}");
  hRecoMuBinpT.GetXaxis()->SetTitle("Muon p_{T}");
  hGMTMuBinpT.GetXaxis()->SetTitle("Muon p_{T}");
  hRecoMuGMTpT.GetXaxis()->SetTitle("Reco Muon p_{T}");

  hDTRPCBX.GetYaxis()->SetTitle("RPC BX");
  hCSCRPCBX.GetYaxis()->SetTitle("RPC BX");
  hDTCSCBX.GetYaxis()->SetTitle("CSC BX");

  hpTresidual_r.GetYaxis()->SetTitle("Average #Delta p_{T} / p_{T}");
  hpTresidual_z.GetYaxis()->SetTitle("Average #Delta p_{T} / p_{T}");
  hpTresidual_eta.GetYaxis()->SetTitle("Average #Delta p_{T} / p_{T}");
  hRecoMuGMTpT.GetYaxis()->SetTitle("GMT Muon p_{T}");

  hGMTMuBinpT.SetLineColor(kRed);
  hRecoMuBinpT.SetLineWidth(2);
  hGMTMuBinpT.SetLineWidth(2);

  TCanvas c0;
  c0.cd(); hTT11bx.Draw();
  c0.SaveAs("results/c0.gif");

  TCanvas c1;
  c1.Divide(2,2);
  c1.cd(1); hDTTFbx.Draw();
  c1.cd(2); hCSCTFbx.Draw();
  c1.cd(3); hRPCbbx.Draw();
  c1.cd(4); hRPCfbx.Draw();
  c1.SaveAs("results/c1.gif");

  TCanvas c2;
  c2.Divide(2,2);
  c2.cd(1); hGMTbx.Draw();
  c2.cd(2); hGMTpt.Draw();
  c2.cd(3); hGMTeta.Draw();
  c2.cd(4); hGMTphi.Draw();
  c2.SaveAs("results/c2.gif");

  TCanvas c3;
  c3.Divide(2,2);
  c3.cd(1); hGMT2mubx.Draw();
  c3.cd(2); hGMT2mupt.Draw();
  c3.cd(3); hGMT2mueta.Draw();
  c3.cd(4); hGMT2muphi.Draw();
  c3.SaveAs("results/c3.gif");

  TCanvas c4;
  c4.cd(); hGMT2mu1BX.Draw();
  c4.SaveAs("results/c4.gif");

  TCanvas c5;
  c5.Divide(2,2);
  c5.cd(1); hDTRPCBX.Draw("COLZ");
  c5.cd(2); hCSCRPCBX.Draw("COLZ");
  c5.cd(3); hDTCSCBX.Draw("COLZ");
  c5.SaveAs("results/c5.gif");

  if(_isReco){

    TCanvas c6;
    c6.cd();
    hpTresidual.Draw();
    c6.SaveAs("results/c6.gif");

    gStyle->SetOptStat(0);

    TCanvas c7;
    c7.cd();
    hRecoMupT.Draw();
    c7.SaveAs("results/c7.gif");

    TCanvas c8("c8","c8",1200,1200);
    c8.cd();
    c8.SetLogx(true);
    c8.SetLogy(true);
    hRecoMuGMTpT.Draw("COLZ");
    c8.SaveAs("results/c8.gif");

    TCanvas c9;
    c9.cd();
    hGMTMuBinpT.Draw("E");
    hRecoMuBinpT.Draw("ESAME");
    c9.SaveAs("results/c9.gif");

    TCanvas c10;
    c10.Divide(2,1);
    c10.cd(1);
    hpTresidual_r.Draw();
    c10.cd(2);
    hpTresidual_z.Draw();
    c10.SaveAs("results/c10.gif");

    TCanvas c11;
    c11.cd();
    hpTresidual_eta.Draw();
    c11.SaveAs("results/c11.gif");
  }

  return;
}

void RunL1(Int_t whichFileAndLumiToUse=1){

  std::string L1NtupleFileName = "";

  Bool_t isReco = false;

  if(whichFileAndLumiToUse==1){
    L1NtupleFileName = "root://lxcms02//data2/p/pellicci/L1DPG/root/Data/Cosmics/229713_MinimumBias/L1Tree.root";
  }
  else if(whichFileAndLumiToUse==2){
    L1NtupleFileName = "root://lxcms02//data2/p/pellicci/L1DPG/root/Data/Cosmics/229713_Cosmics/L1Tree.root";
  }
  else if(whichFileAndLumiToUse==3){
    L1NtupleFileName = "root://lxcms02//data2/p/pellicci/L1DPG/root/Data/Cosmics/232956_Cosmics/L1Tree.root";
  }
  else if(whichFileAndLumiToUse==4){
    L1NtupleFileName = "root://lxcms02//data2/p/pellicci/L1DPG/root/Data/Cosmics/233238_Cosmics/L1Tree.root";
  }
  else if(whichFileAndLumiToUse==5){
    L1NtupleFileName = "file:///data2/p/pellicci/L1DPG/root/Data/Cosmics/238492_CosmicsSP/L1Tree.root";
    isReco = true;
  }
  else if(whichFileAndLumiToUse==6){
    L1NtupleFileName = "file:///data2/p/pellicci/L1DPG/root/Data/Cosmics/Full_SP_310315/L1Tree.root";
    isReco = true;
  }
  else if(whichFileAndLumiToUse==7){
    L1NtupleFileName = "root://lxcms02//data2/p/pellicci/L1DPG/root/Data/Collisions/256843_ZeroBias.root";
  }
  else{
    std::cout << std::endl << "ERROR: Please define a ntuple file which is in the allowed range! You did use: whichFileAndLumiToUse = " << whichFileAndLumiToUse << " This is not in the allowed range" << std::endl << std::endl;
  }

  analyze_L1 a(isReco);
  a.Open(L1NtupleFileName);
  a.Loop();

 return;
}

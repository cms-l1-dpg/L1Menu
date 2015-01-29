#include "L1Ntuple.h"

#include "TH1F.h"

class analyze_L1 : public L1Ntuple {

 public :

  analyze_L1() {}
  ~analyze_L1() {}

  void Loop();

private :

};

void analyze_L1::Loop(){

  Int_t nevents = GetEntries();
  std::cout << "Running on " << nevents << " events." << std::endl;
  if(nevents > 1000000) nevents = 1000000;

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

      for(int i11=0;i11<5;i11++){
	if( (gt_->tt[i11-2]>>11)&1 ) hTT11bx->Fill(i11-2.);
      }

      Int_t Ndt = gmt_->Ndt;
      for (int idt=0; idt < Ndt; idt++){
	int bx = gmt_ -> Bxdt[idt];
	hDTTFbx.Fill(bx);
      }
 
      Int_t Ncsc = gmt_->Ncsc;
      for (int icsc=0; icsc < Ncsc; icsc++){
	int bx = gmt_ -> Bxcsc[icsc];
	hCSCTFbx.Fill(bx);
      }

      Int_t Nrpcb = gmt_->Nrpcb;
      for (int irpcb=0; irpcb < Nrpcb; irpcb++){
	int bx = gmt_ -> Bxrpcb[irpcb];
	hRPCbbx.Fill(bx);
      }

      Int_t Nrpcf = gmt_->Nrpcf;
      for (int irpcf=0; irpcf < Nrpcf; irpcf++){
	int bx = gmt_ -> Bxrpcf[irpcf];
	hRPCfbx.Fill(bx);
      }
 
      Int_t Nmu = gmt_->N;
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

      if(Nmu > 1 && fabs((gmt_->Pt[0]-gmt_->Pt[1])/gmt_->Pt[0]) < 0.1){
	hGMT2mubx->Fill(gmt_CandBx[0]);
	hGMT2mubx->Fill(gmt_CandBx[1]);
      }

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

  hTT11bx.GetXaxis()->SetTitle("TT11 BX");

  hDTTFbx.GetXaxis()->SetTitle("DTTF BX");
  hCSCTFbx.GetXaxis()->SetTitle("CSCTF BX");
  hRPCbbx.GetXaxis()->SetTitle("RPCb BX");
  hRPCfbx.GetXaxis()->SetTitle("RPCf BX");
  hGMTbx.GetXaxis()->SetTitle("GMT BX");

  hGMTpt.GetXaxis()->SetTitle("L1 Muon #p_{T} [GeV]");
  hGMTeta.GetXaxis()->SetTitle("L1 Muon #eta");
  hGMTph.GetXaxis()->SetTitle("L1 Muon #phi");

  hGMT2mubx.GetXaxis()->SetTitle("BX for identical muons");

  TCanvas c0;
  c0.cd(); hTT11bx->Draw();
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
  c3.cd(); hGMT2mubx->Draw();
  c3.SaveAs("results/c3.gif");

  return;
}

void RunL1(Int_t whichFileAndLumiToUse=1){

  std::string L1NtupleFileName = "";

  if(whichFileAndLumiToUse==1){
    L1NtupleFileName = "root://lxcms02//data2/p/pellicci/L1DPG/root/Data/Cosmics/229713_MinimumBias/L1Tree.root";
  }
  else if(whichFileAndLumiToUse==2){
    L1NtupleFileName = "root://lxcms02//data2/p/pellicci/L1DPG/root/Data/Cosmics/229713_Cosmics/L1Tree.root";
  }
  else{
    std::cout << std::endl << "ERROR: Please define a ntuple file which is in the allowed range! You did use: whichFileAndLumiToUse = " << whichFileAndLumiToUse << " This is not in the allowed range" << std::endl << std::endl;
  }

  analyze_L1 a;
  a.Open(L1NtupleFileName);
  a.Loop();

 return;
}

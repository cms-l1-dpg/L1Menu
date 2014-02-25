#include<string>
#include<iostream>

#include"TF1.h"
#include"TH1F.h"
#include"TMath.h"
#include"TFile.h"

TH1F * getPoissonTH1F(double avPU, std::string histoName)
{

  TH1F * histo = new TH1F(histoName.c_str(),histoName.c_str(),101,-0.5,100.5);
  
  int nBins = histo->GetNbinsX();
  
  for (int iBin=1; iBin<nBins; ++iBin) 
    {
      float binCenter = histo->GetBinCenter(iBin);
      histo->SetBinContent(iBin,TMath::PoissonI(binCenter,avPU));
    }
  
  return histo;

}

void printHisto(TH1F *histo)
{

  std::string hName = histo->GetName();
  float mean = histo->GetMean();
  float rms  = histo->GetRMS();
 
  std::cout << "Histo : "  << hName
	    << "\tmean : " << mean
	    << "\tRMS^2 : "  << rms*rms
	    << std::endl;

}    

void genPUReweightinFile(double origAvPU, double targetAvPU, 
			 std::string fileName)
{
  
  TFile * outputFile = new TFile(fileName.c_str(),"RECREATE");
  
  outputFile->cd();
  
  TH1F *orig   = getPoissonTH1F(origAvPU,"productionPileUpHisto"); 
  TH1F *target = getPoissonTH1F(targetAvPU,"targetPileUpHisto");

  printHisto(orig);
  printHisto(target);

  TH1F *ratio = static_cast<TH1F*>(orig->Clone("ratio"));
  ratio->Divide(target,orig);
  
  orig->Write();
  target->Write(); 

  ratio->Write(); 

  outputFile->Close();
  
  delete outputFile;

}

    

    
  

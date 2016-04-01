#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>

#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TList.h"
#include "THStack.h"
#include "TIterator.h"
#include "TObject.h"
#include "TClass.h"
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>

#include "tdrstyle.C"

std::string BuildLeg(std::string filename);
void getPlotList( const std::string & fileName,
		  std::map<std::string, std::vector<std::string> > & plotList )
{
  TH1::AddDirectory(kFALSE);

  // build list of histogram names
  TFile rootfile( fileName.c_str() );

  std::vector< std::string > dirList;
	  
  TList *dirs = rootfile.GetListOfKeys();
  TIterator *itdir = dirs->MakeIterator();
  TObject *nextdir;

  while ( (nextdir = itdir->Next()) ) {
    
    if( nextdir->IsFolder())
      dirList.push_back( nextdir->GetName() );
    else 
      plotList[""].push_back( ( nextdir->GetName() ) );

    
  }
  
  std::vector<std::string>::const_iterator dirIt  = dirList.begin();
  std::vector<std::string>::const_iterator dirEnd = dirList.end();
  
  for (;dirIt!=dirEnd;++dirIt){
    
    TDirectory * thisdir = (TDirectory*)( rootfile.Get( dirIt->c_str() ) );
    
    TList * dircontent = thisdir->GetListOfKeys();
    TIterator * thisplot = dircontent->MakeIterator();
    TObject * nextplot;
 
    const std::string & dirName = (*dirIt); 
    
    while ( (nextplot = thisplot->Next()) ) {
      plotList[dirName].push_back( (  dirName + "/" + nextplot->GetName() ) );
    }

  }

  rootfile.Close();

}


void getRange( TH1* plot, float & minY, float & maxY )
{

  minY = 1E10;
  maxY = 0.;

  int nBins = plot->GetNbinsX();
  
  for ( int iBin=0; iBin<=nBins; ++iBin )
    {

      float val = plot->GetBinContent( iBin );
      minY = ( val > 0.001 && val < minY ) ? val : minY;
      maxY = ( val > 0.001 && val > maxY ) ? val : maxY;

    }

  minY = ( fabs(minY - 1E10) > 0.01 ) ? minY*0.5 : 0.;
  maxY = ( fabs(maxY - 0.  ) > 0.01 ) ? maxY*2. : 1.;

  return;

}


void plot( std::vector<TH1*> plots,
	   std::string &baseDir, std::string outputDir,
       std::vector<std::string> legs) 
{
  
  if (plots.at(0))
    {
      std::cout << "Plotting : " << plots.at(0)->GetName() << std::endl;
  
      // plot everything
      TCanvas *c = new TCanvas();

      c->cd();

      if (dynamic_cast<TH1F*>(plots.at(0))) 
	{

         float minY = 0;
         float maxY = 0;

	  for (size_t iPlot=0; iPlot<plots.size(); ++iPlot){
	    float tmp_minY = 999.;
	    float tmp_maxY = -999.;
            getRange( plots.at(iPlot), tmp_minY, tmp_maxY );
	    if(tmp_minY < minY) minY = tmp_minY;
	    if(tmp_maxY > maxY) maxY = tmp_maxY;
	  }

	  std::cout << "min " << minY << std::endl;
	  std::cout << "max " << maxY << std::endl;
	  
	  TPad *pPlot = ( plots.size()>1 ) ? new TPad("pPlot","",0.05,0.26,0.99,0.99) :
	    new TPad("pPlot","",0.05,0.01,0.99,0.99) ;
	  
	  pPlot->Draw();
	  pPlot->SetGrid();
      TLegend *leg = new TLegend(0.4055405,0.7411604,0.9767243,0.896555,NULL,"brNDC");
      leg->SetFillColor(0);
      leg->SetBorderSize(0);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(62);
	  
	  c->cd();
	  
	  TPad *pRatio = ( plots.size()>1 ) ? new TPad("pRatio","",0.05,0.01,0.99,0.25) : 0;
	  if(pRatio)
	    {
	      pRatio->Draw();
	      pRatio->SetGrid();
	    }
	  
	  for (size_t iPlot=0; iPlot<plots.size(); ++iPlot) 
	    {
	      pPlot->cd();
	      plots.at(iPlot)->SetLineColor( iPlot+1 );
	      plots.at(iPlot)->SetFillColor( iPlot+1 );
	      plots.at(iPlot)->SetMarkerColor( iPlot+1 );
	      plots.at(iPlot)->SetMarkerStyle( 21 + iPlot );
	       
          getRange( plots.at(iPlot), minY, maxY );
	      plots.at(iPlot)->GetYaxis()->SetRangeUser( minY, maxY );
	      plots.at(iPlot)->GetYaxis()->SetTitleSize(0.06);
	      plots.at(iPlot)->GetYaxis()->SetTitleOffset(0.7);
	      plots.at(iPlot)->Draw( iPlot ? "samePE1" : "PE1" );
          leg->AddEntry(plots.at(iPlot), BuildLeg(legs.at(iPlot)).c_str(), "P");
	      

	      if ( iPlot>0 ) 
		{
		  pRatio->cd();
		  TH1* den = plots.at(0);
		  std::string name = std::string(plots.at(iPlot)->GetName()) + "_Eff";
		  TH1* eff = static_cast<TH1*>( plots.at(iPlot)->Clone( name.c_str() ) );
		  eff->Divide(den);
		  eff->SetTitle(";;Ratio");
		  
		  eff->SetLineColor( iPlot+1 );
		  eff->SetFillColor( iPlot+1 );
		  eff->SetMarkerColor( iPlot+1 );
		  eff->SetMarkerStyle( 21 + iPlot );
	      eff->GetYaxis()->SetTitleSize(0.14);
          eff->GetYaxis()->SetTitleOffset(0.4);
		  eff->GetYaxis()->SetRangeUser( .1, 5.);
		  eff->GetYaxis()->SetLabelSize( .11);
		  
		  eff->Draw( iPlot>1 ? "samePE1" : "PE1" );
		}
	    }
	  
      pPlot->cd();
	  pPlot->SetLogy();
      leg->Draw();
	}
      else if(dynamic_cast<TH2F*>(plots.at(0)) && plots.size() == 2) 
	{ 
	  c->SetGrid();
	  
	  TH1* den = plots.at(0);
	  std::string name = std::string(plots.at(1)->GetName()) + "_Eff";
	  TH1* eff = static_cast<TH1*>( plots.at(1)->Clone( name.c_str() ) );
	  eff->Divide(den);
	  eff->SetTitle(";;Ratio");
	  eff->SetMinimum(1);
	  eff->SetMaximum(2.5);
	  eff->Draw("colz");
	}
      
      std::string path = baseDir + "/" + outputDir;
      system( (std::string("mkdir -p ") + path).c_str() );
      
      c->Update();
      std::string printname = path + "/" + plots.at(0)->GetName();
      c->Print ( ( printname + ".gif" ).c_str() ); 
      c->Print ( ( printname + ".root" ).c_str() ); 
      //c->Print ( ( printname + ".C" ).c_str() ); 
    }
  
}

void plotAll(std::vector<std::string> &files,
	     std::string &baseDir) 
{

  size_t nFiles = 0;
  std::vector<TFile*> filesRoot;

  for (size_t iFile=0; iFile<files.size(); ++iFile) {
    filesRoot.push_back(new TFile(files.at(iFile).c_str(),"READONLY"));
  }

  system( (std::string("mkdir -p ") + baseDir).c_str() );
  
  std::map<std::string, std::vector<std::string> > plotNames;

  getPlotList(files.at(0),plotNames);

  std::map<std::string, std::vector<std::string> >::const_iterator plotDirIt  = plotNames.begin();
  std::map<std::string, std::vector<std::string> >::const_iterator plotDirEnd = plotNames.end();

  for(;plotDirIt!=plotDirEnd;++plotDirIt) {

    std::vector<std::string>::const_iterator plotIt  = plotDirIt->second.begin();
    std::vector<std::string>::const_iterator plotEnd = plotDirIt->second.end();

    for(;plotIt!=plotEnd;++plotIt) {

      std::vector<TH1*> plots;
      std::vector<std::string> legs;

      for (size_t iFile=0; iFile<filesRoot.size(); ++iFile) {
        plots.push_back(static_cast<TH1*>( filesRoot.at(iFile)->Get( plotIt->c_str() )  )); 
        legs.push_back(files.at(iFile));
      }
      
      plot(plots,baseDir,plotDirIt->first, legs);
      
    }
  }
  
}

// ===  FUNCTION  ============================================================
//         Name:  BuildLeg
//  Description:  
// ===========================================================================
std::string BuildLeg(std::string filename)
{
  std::string temp=filename;
  size_t found = filename.find_last_of("/");
  if (found != std::string::npos)
  {
    temp = filename.substr(found+1);
  }

  found = temp.find_last_of(".");
  if (found != std::string::npos)
  {
    temp = temp.substr(0, found);
  }

  return temp;
}       // -----  end of function BuildLeg  -----

int main(int argc, char* argv[]) 
{  

  setTDRStyle();
  
  if ( argc<2 ) {
    std::cout << "Error in number of arguments: " << argc << std::endl;
    std::cout << "Passed args: " << argc << std::endl;
    for ( int i = 1; i < argc; ++i ) {
      std::cout << "\t" << argv[i] << std::endl;
    }
    std::cout << "Usage: \n\t\t " <<  argv[0] << " <first inputfile> <second inputfile> ... "
	      << std::endl << std::endl;
    return -1;
  }

  
  std::vector<std::string> files;
  for (int iArg=1; iArg<argc; ++iArg) {
    files.push_back(argv[iArg]);
  }

  std::string baseDir = "results/comparePlots";

  plotAll(files,baseDir);

}

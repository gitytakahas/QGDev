#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>
#include "TChain.h"
#include "TFile.h"
#include "TH3D.h"
#include "TProfile2D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "binClass.h"
#include "binningConfigurations.h"
#include "treeLooper.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc"


void setCorrectMinimum(TProfile* h){
  bool first = true;
  float minimum = 0;
  for(int i = 0; i < h->GetSize(); ++i){
    if(h->GetBinContent(i) == 0) continue;
    if(first || h->GetBinContent(i) < minimum){
      first = false;
      minimum = h->GetBinContent(i);
    }
  }
  h->SetMinimum(minimum);
}


int main(int argc, char**argv){
 for(bool useBins : {true}){ 
 for(TString xVar : {"rho"}){

  std::vector<TString> files	= {"QCD_AllPtBins"};
  //std::vector<TString> files	= {"QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14"};
//std::vector<TString> files	= {"TTJets"};
  //std::vector<TString> jetTypes = {"AK4","AK4chs"};
  std::vector<TString> jetTypes = {"AK4chs"};

  binClass bins = getSmallEtaBinning();
  bins.printBinRanges();

  // Define x and y-axis variable and bins
  std::vector<float> xAxisBins = {0,1}; TString xTitle = ""; bool xLog = false;
  if(xVar == "eta"){ 	      		xAxisBins = bins.getBins(40, 0, 4, false, {4.2,4.4,4.7});		xTitle = "#eta";}
  if(xVar == "additionalJets"){		xAxisBins = bins.getBins(6, -0.5, 5.5, false, {});		 	xTitle = "N_{jets>20GeV} within #DeltaR < 0.8";}
  if(xVar == "rho"){  			xAxisBins = bins.getBins(9, 5, 50, false, {});		 		xTitle = "#rho";}
  if(xVar == "nPileUp"){  		xAxisBins = bins.getBins(12, 9.5, 69.5, false, {});		 	xTitle = "nPileUp";}
  if(xVar == "nPriVtxs"){  		xAxisBins = bins.getBins(9, 9.5, 54.5, false, {});		 	xTitle = "nPriVtxs";}

  double xAxisBinsArray[1000]; std::copy(xAxisBins.begin(), xAxisBins.end(), xAxisBinsArray);

  // For different samples and jet types
  for(TString file : files){
    for(TString jetType : jetTypes){
      std::cout << "Making plots for " << jetType << " in file " << file << "..." << std::endl;
      system("rm -rf plots/profile/" + file + "/" + jetType + "/" + xVar + (useBins ? "" : "_nobins"));

      treeLooper t(file, jetType);							// Init tree
      bins.setReference("pt",  &t.pt);
      bins.setReference("eta", &t.eta);
      bins.setReference("rho", &t.rho);

      // Set x-axis variable
      float& xAxisVar = (xVar == "additionalJets" ? 		t.additionalJets : 
                        (xVar == "rho" ? 			t.rho : t.eta));
//    int& xAxisVar =   (xVar == "nPileUp" ? 			t.nPileUp : t.nPriVtxs);

      // Init local QGLikelihoodCalculator
      QGLikelihoodCalculator localQG("../data/pdfQG_" + jetType + "_13TeV_76X.root");

      // Creation of histos
      std::map<TString, TProfile*> plots;
      for(TString type : {"quark","gluon"}){
        for(TString binName : bins.getAllBinNames()){
          TString histName = "_" + type + "_" + binName;
          plots["axis2" + histName] = new TProfile("axis2" + histName, "mean -log(#sigma_{2})",           xAxisBins.size()-1, xAxisBinsArray);
          plots["ptD"   + histName] = new TProfile("ptD"   + histName, "mean p_{T}D",                     xAxisBins.size()-1, xAxisBinsArray);
          plots["mult"  + histName] = new TProfile("mult"  + histName, "mean multiplicity",               xAxisBins.size()-1, xAxisBinsArray);
          plots["qg"    + histName] = new TProfile("qg"    + histName, "mean quark-gluon likelihood",     xAxisBins.size()-1, xAxisBinsArray);
        }
      }

      // Fill histos
      TString binName;
      while(t.next()){
        if(!bins.getBinName(binName)) 	continue;									// Find bin and return false if outside ranges 
        if(t.jetIdLevel < 3) 		continue;									// Select tight jets
        if(!t.matchedJet) 		continue; 									// Only matched jets
        if(t.bTag) 			continue;									// Anti-b tagging
	if(!t.balanced) 		continue;									// Take only two leading jets with pt3 < 0.15*(pt1+pt2)
        if(t.nGenJetsInCone != 1 || t.nJetsForGenParticle != 1 || t.nGenJetsForGenParticle != 1) continue;		// Use only jets matched to exactly one gen jet and gen particle, and no other jet candidates

        TString type;
        if(t.partonId == 21) 	 	type = "gluon";
        else if(fabs(t.partonId) < 4) 	type = "quark";
        else continue;

        if(t.mult < 3 || t.axis2 > 99) continue; 									// Take jets with at least 3 particles

        float qg = localQG.computeQGLikelihood(          t.pt, t.eta, t.rho, {(float)t.mult, t.ptD, t.axis2});

        TString histName = "_" + type + "_" + binName;
        plots["axis2" + histName]->Fill(xAxisVar, t.axis2, 	t.weight);
        plots["ptD"   + histName]->Fill(xAxisVar, t.ptD, 	t.weight);
        plots["mult"  + histName]->Fill(xAxisVar, t.mult, 	t.weight);
        plots["qg"    + histName]->Fill(xAxisVar, qg, 		t.weight);
      }

      // Stacking, cosmetics and saving
      for(auto& plot : plots){
        if(!plot.first.Contains("gluon")) continue;
        TCanvas c;
        if(xLog) c.SetLogx();

        TLegend l(0.3,0.91,0.7,0.98);
        l.SetNColumns(2);
        l.SetFillColor(kWhite);
        l.SetBorderSize(0);
       
        float maximum = 0, minimum = 99999; 
        for(TString type : {"quark", "gluon"}){
          TString histName = plot.first;
          histName.ReplaceAll("gluon", type);
          setCorrectMinimum(plots[histName]);
          if(maximum < plots[histName]->GetMaximum()) maximum = plots[histName]->GetMaximum();
          if(minimum > plots[histName]->GetMinimum()) minimum = plots[histName]->GetMinimum();
        }

        plot.second->GetXaxis()->SetTitle(xTitle);
        plot.second->GetYaxis()->SetTitle(plot.second->GetTitle());
        plot.second->SetTitle("");
        plot.second->SetStats(0);
        plot.second->SetMaximum(maximum*1.01);
        plot.second->SetMinimum(minimum*0.99);
        plot.second->Draw("HIST LP");
        std::vector<int> color = {2,46,3,30,4,39};
        int i = 0;
        for(TString type : {"quark", "gluon"}){
          TString histName = plot.first;
          histName.ReplaceAll("gluon", type);
          plots[histName]->SetMarkerStyle(20);
          plots[histName]->SetMarkerSize(0.5);
          plots[histName]->SetMarkerColor(color[i]);
          plots[histName]->SetLineColor(color[i++]);
          plots[histName]->Draw("HIST LP same");
          l.AddEntry(plots[histName], type + "s", "p");
        }
        l.Draw();

        bins.printInfoOnPlot(plot.first, jetType);

        TString variable = plot.first(0, plot.first.First("_"));
        TString pdfDir = "./plots/profile/" + file + "/" + jetType + "/" + xVar + (useBins ? "" : "_nobins") + "/" + variable + "/";
        TString pdfName = pdfDir + plot.first + ".pdf";
        pdfName.ReplaceAll("_gluon","");
        system("mkdir -p " + pdfDir);
        c.SaveAs(pdfName);
      }

      for(auto& plot : plots) delete plot.second;
    }
  }
 }
 }
 return 0;
}

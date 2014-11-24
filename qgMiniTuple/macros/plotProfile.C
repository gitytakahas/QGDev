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
#include "treeLooper2.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc"


int main(int argc, char**argv){
 bool useBins = false;
 for(TString xVar : {"deltaRmin", "closebyJetsInCone"}){

//  std::vector<TString> files	= {"QCD_AllPtBins"};
  std::vector<TString> files	= {"QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14"};
  std::vector<TString> jetTypes = {"AK4chs"};

  // Define binning for plots (i.e. separate plots for each bin)
  binClass bins;
  if(useBins) bins.setBinRange("eta", 	"#eta",		{0,1.3,1.5,2,2.5,3,4.7});
  if(useBins) bins.setBinRange("pt" , 	"p_{T}",	bins.getBins(20, 20, 2000, true, {6500}));				// i.e. 20 bins from 20 to 2000 with log=true and with an additional bin up to 6500
  bins.printBinRanges();

  // Link some bins to be merged because of low statistics (for example higher pT bins at large eta)
  if(useBins){
    for(int i=10; i < 21; ++i) bins.setLinks("eta5_pt9", {TString::Format("eta5_pt%d",i)});
    for(int i=14; i < 21; ++i) bins.setLinks("eta4_pt13", {TString::Format("eta4_pt%d",i)});
    for(int i=16; i < 21; ++i) bins.setLinks("eta3_pt15", {TString::Format("eta3_pt%d",i)});
    for(int i=18; i < 21; ++i) bins.setLinks("eta2_pt17", {TString::Format("eta2_pt%d",i)});
    for(int i=19; i < 21; ++i) bins.setLinks("eta1_pt18", {TString::Format("eta1_pt%d",i)});
  }

  // Define x and y-axis variable and bins
  std::vector<float> xAxisBins = {0,1}; TString xTitle = ""; bool xLog = false;
  if(xVar == "eta"){ 	      		xAxisBins = bins.getBins(40, 0, 4, false, {4.2,4.4,4.7});		xTitle = "#eta";}
  if(xVar == "deltaRmin"){ 		xAxisBins = bins.getBins(16, .4, .8, false, {0.85,0.9,1,1.2}); 		xTitle = "#DeltaR_{min}";}
  if(xVar == "closebyJetsInCone"){  	xAxisBins = bins.getBins(9, -0.5, 8.5, false, {});		 	xTitle = "N_{jets} within #DeltaR < 0.8";}
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

      // Set x-axis variable
      float& xAxisVar = (xVar == "deltaRmin" ? 		t.closestJetdR : 
                        (xVar == "closebyJetsInCone" ? 	t.closebyJetsInCone : 
                        (xVar == "rho" ? 		t.rho : 
			t.eta)));
//      int& xAxisVar =   (xVar == "nPileUp" ? 		t.nPileUp : t.nPriVtxs);

      // Init local QGLikelihoodCalculator
      QGLikelihoodCalculator localQG("../data/pdfQG_" + jetType + "_13TeV.root");
      QGLikelihoodCalculator localQG_cdf("../data/pdfQG_" + jetType + "_fineBinning_13TeV.root");

      // Creation of histos
      std::map<TString, TProfile*> plots;
      for(TString type : {"quark","gluon"}){
        for(TString binName : bins.getAllBinNames()){
          TString histName = "_" + type + "_" + binName;
          for(TString dRcut : {"", "_dR2", "_dR3"}){
            plots["axis2" + dRcut + histName] = new TProfile("axis2" + dRcut + histName, "mean -log(#sigma_{2})",           xAxisBins.size()-1, xAxisBinsArray);
            plots["ptD"   + dRcut + histName] = new TProfile("ptD"   + dRcut + histName, "mean p_{T}D",                     xAxisBins.size()-1, xAxisBinsArray);
            plots["nChg"  + dRcut + histName] = new TProfile("nChg"  + dRcut + histName, "mean multiplicity",               xAxisBins.size()-1, xAxisBinsArray);
            plots["mult"  + dRcut + histName] = new TProfile("mult"  + dRcut + histName, "mean multiplicity",               xAxisBins.size()-1, xAxisBinsArray);
            plots["qg"    + dRcut + histName] = new TProfile("qg"    + dRcut + histName, "mean quark-gluon likelihood",     xAxisBins.size()-1, xAxisBinsArray);
            plots["cdf"   + dRcut + histName] = new TProfile("cdf"   + dRcut + histName, "mean quark-gluon CDF-likelihood", xAxisBins.size()-1, xAxisBinsArray);
          }
        }
      }

      // Fill histos
      while(t.next()){
        if(!bins.update()) 	continue;										// Find bin and return false if outside ranges 
        if(t.jetIdLevel < 3) 	continue;										// Select tight jets
        if(!t.matchedJet) 	continue; 										// Only matched jets
        if(t.bTag) 		continue;										// Anti-b tagging
//	if(!t.balanced) 	continue;										// Take only two leading jets with pt3 < 0.15*(pt1+pt2)
//      if(t.nGenJetsInCone != 1 || t.nJetsForGenParticle != 1 || t.nGenJetsForGenParticle != 1) continue;		// Use only jets matched to exactly one gen jet and gen particle, and no other jet candidates

        TString type;
        if(t.partonId == 21) 	 	type = "gluon";
        else if(fabs(t.partonId) < 4) 	type = "quark";
        else continue;

        for(TString dRcut : {"","_dR2","_dR3"}){
          int& mult    = (dRcut == ""? t.mult  : (dRcut == "_dR2"? t.mult_dR2  : t.mult_dR3));
          int& nChg    = (dRcut == ""? t.nChg  : (dRcut == "_dR2"? t.nChg_dR2  : t.nChg_dR3));
          float& ptD   = (dRcut == ""? t.ptD   : (dRcut == "_dR2"? t.ptD_dR2   : t.ptD_dR3));
          float& axis2 = (dRcut == ""? t.axis2 : (dRcut == "_dR2"? t.axis2_dR2 : t.axis2_dR3));
          if(mult < 3 || axis2 > 99) continue; 										// Take jets with at least 3 particles

          float qg = localQG.computeQGLikelihood(          t.pt, t.eta, t.rho, {(float)mult, ptD, axis2});
          float qgcdf = localQG_cdf.computeQGLikelihoodCDF(t.pt, t.eta, t.rho, {(float)mult, ptD, axis2});

          TString histName = "_" + type + "_" + bins.name;
          plots["axis2" + dRcut + histName]->Fill(xAxisVar, axis2, 	t.weight);
          plots["ptD"   + dRcut + histName]->Fill(xAxisVar, ptD, 	t.weight);
          plots["nChg"  + dRcut + histName]->Fill(xAxisVar, nChg, 	t.weight);
          plots["mult"  + dRcut + histName]->Fill(xAxisVar, mult, 	t.weight);
          plots["qg"    + dRcut + histName]->Fill(xAxisVar, qg, 	t.weight);
          plots["cdf"   + dRcut + histName]->Fill(xAxisVar, qgcdf, 	t.weight);
        }
      }

      // Stacking, cosmetics and saving
      for(auto& plot : plots){
        if(!plot.first.Contains("gluon") || !plot.first.Contains("dR3")) continue;
        TCanvas c;
        if(xLog) c.SetLogx();

        TLegend l(0.3,0.91,0.7,0.98);
        l.SetNColumns(2);
        l.SetFillColor(kWhite);
        l.SetBorderSize(0);
       
        float maximum = 0, minimum = 99999; 
        for(TString type : {"quark", "gluon"}){
          for(TString dRcut : {"", "_dR2", "_dR3"}){
            TString histName = plot.first;
            histName.ReplaceAll("gluon", type);
            histName.ReplaceAll("_dR3", dRcut);
            if(maximum < plots[histName]->GetMaximum()) maximum = plots[histName]->GetMaximum();
            if(plots[histName]->GetMinimum() > 0.01 && minimum > plots[histName]->GetMinimum()) minimum = plots[histName]->GetMinimum();
          }
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
        for(TString dRcut : {"", "_dR2", "_dR3"}){
          for(TString type : {"quark", "gluon"}){
            TString histName = plot.first;
            histName.ReplaceAll("gluon", type);
            histName.ReplaceAll("_dR3", dRcut);
            plots[histName]->SetMarkerStyle(20);
            plots[histName]->SetMarkerSize(0.5);
            plots[histName]->SetMarkerColor(color[i]);
            plots[histName]->SetLineColor(color[i++]);
            plots[histName]->Draw("HIST LP same");
            l.AddEntry(plots[histName], type + "s" + (dRcut == "_dR3" ? " (#DeltaR < 0.3)" : (dRcut == "_dR2" ? " (#DeltaR < 0.2)" : "")), "p");
          }
        }
        l.Draw();

        bins.printInfoOnPlot(plot.first, jetType);

        TString variable = plot.first(0, plot.first.First("_"));
        TString pdfDir = "./plots/profile/" + file + "/" + jetType + "/" + xVar + (useBins ? "" : "_nobins") + "/" + variable + "/";
        TString pdfName = pdfDir + plot.first + ".pdf";
        pdfName.ReplaceAll("_gluon","").ReplaceAll("_dR3","");
        system("mkdir -p " + pdfDir);
        c.SaveAs(pdfName);
      }

      for(auto& plot : plots) delete plot.second;
    }
  }
 }
 return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>
#include "TH3D.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "binClass.h"
#include "binningConfigurations.h"
#include "treeLooper.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc"

void setCorrectMinimum(TH2D* h){
  bool first = true;
  float minimum = 0;
  for(int i = 0; i < h->GetSize(); ++i){
    if(h->GetBinContent(i) == 0) continue;
    if(first || h->GetBinContent(i) < minimum){
      first = false;
      minimum = h->GetBinContent(i);
    }
  }
  h->SetMinimum(minimum*0.99);
}


void makeProfile2D(TString xAndyVar, bool useBins = false, bool extraBinInformation = false){
//std::vector<TString> files	= {"QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_PU20"};
  std::vector<TString> files	= {"TTJets"};
  std::vector<TString> jetTypes = {"AK4chs"};

  // Define binning for plots (i.e. separate plots for each bin)
  binClass bins;
  if(useBins) bins = getSmallEtaBinning();

  // Define x and y-axis variable and bins
  std::vector<float> xAxisBins, yAxisBins; TString xVar, yVar; bool xLog, yLog;
  if(xAndyVar.Contains("eta_")){ 			xAxisBins = bins.getBins(40, 0, 4, false, {4.2,4.4,4.7}); 	xVar = "#eta";					xLog = false;}
  if(xAndyVar.Contains("rho_")){ 			xAxisBins = bins.getBins(8, 0, 40, false, {50});		xVar = "#rho";					xLog = false;}
  if(xAndyVar.Contains("additionalJets_")){		xAxisBins = bins.getBins(6, -0.5, 5.5, false, {});		xVar = "N_{jets>20GeV} within #DeltaR < 0.8";	xLog = false;}
  if(xAndyVar.Contains("_eta")){ 			yAxisBins = bins.getBins(40, 0, 4, false, {4.2,4.4,4.7}); 	yVar = "#eta";					yLog = false;}
  if(xAndyVar.Contains("_pt")){				yAxisBins = bins.getBins(20, 20, 2000, true, {6500});		yVar = "p_{T}"; 				yLog = true;}
  if(xAndyVar.Contains("_nPileUp")){			yAxisBins = bins.getBins(10, 10, 70, false);			yVar = "nPileUp"; 				yLog = false;}
  if(xAndyVar.Contains("_nPriVtxs")){			yAxisBins = bins.getBins(10, 5, 55, false);			yVar = "nPriVtxs"; 				yLog = false;}

  // Needed for 3D plots (which will be needed to get kurtosis and skewness)
  std::vector<float> axis2Bins	= bins.getBins(100, 0, 8, false);
  std::vector<float> multBins	= bins.getBins(100, 0.5, 100.5, false);
  std::vector<float> otherBins	= bins.getBins(100, -0.001, 1.001, false);

  double xAxisBinsArray[1000]; std::copy(xAxisBins.begin(), xAxisBins.end(), xAxisBinsArray);
  double yAxisBinsArray[1000]; std::copy(yAxisBins.begin(), yAxisBins.end(), yAxisBinsArray);
  double otherBinsArray[1000]; std::copy(otherBins.begin(), otherBins.end(), otherBinsArray);
  double axis2BinsArray[1000]; std::copy(axis2Bins.begin(), axis2Bins.end(), axis2BinsArray);
  double multBinsArray[1000];  std::copy(multBins.begin(),  multBins.end(),  multBinsArray);

  // For different samples and jet types
  for(TString file : files){
    for(TString jetType : jetTypes){
      std::cout << "Making plots for " << jetType << " in file " << file << "..." << std::endl;
      system("rm -rf plots/profile2D/" + file + "/" + jetType + "/" + xAndyVar + (useBins ? "" : "_nobins"));

      treeLooper t(file, jetType);											// Init tree
      bins.setReference("pt",  &t.pt);
      bins.setReference("eta", &t.eta);
      bins.setReference("rho", &t.rho);
      bins.setReference("aj",  &t.additionalJets);

      // Set x and y-axis variables
      float& xAxisVar 	= (xAndyVar.Contains("eta_")?			t.eta :
			  (xAndyVar.Contains("rho_")?			t.rho :
			  (xAndyVar.Contains("additionalJets_")?	t.additionalJets : t.pt)));
      float& yAxisVar 	= (xAndyVar.Contains("_eta")? 			t.eta : t.pt);
//    int& yAxisVar 	= (xAndyVar.Contains("_nPileUp"?		t.nPileUp : t.nPriVtxs));

      // Init local QGLikelihoodCalculator
      QGLikelihoodCalculator localQG("../data/pdfQG_" + jetType + "_13TeV_76X.root");

      // Creation of histos
      std::map<TString, TProfile2D*> plots;
      std::map<TString, TH3D*> plots3D;
      for(TString type : {"quark","gluon"}){
        for(TString binName : bins.getAllBinNames()){
          TString histName = "_" + type + "_" + binName;
          plots["axis2"    + histName] = new TProfile2D("axis2"    + histName, "mean -log(#sigma_{2})",           xAxisBins.size()-1, xAxisBinsArray, yAxisBins.size()-1, yAxisBinsArray);
          plots["ptD"      + histName] = new TProfile2D("ptD"      + histName, "mean p_{T}D",                     xAxisBins.size()-1, xAxisBinsArray, yAxisBins.size()-1, yAxisBinsArray);
          plots["mult"     + histName] = new TProfile2D("mult"     + histName, "mean multiplicity",               xAxisBins.size()-1, xAxisBinsArray, yAxisBins.size()-1, yAxisBinsArray);
          plots["qg"       + histName] = new TProfile2D("qg"       + histName, "mean quark-gluon likelihood",     xAxisBins.size()-1, xAxisBinsArray, yAxisBins.size()-1, yAxisBinsArray);
          if(extraBinInformation){
            plots3D["axis2"    + histName] = new TH3D("axis2"    + histName + "3D", "-log(#sigma_{2})",           xAxisBins.size()-1, xAxisBinsArray, yAxisBins.size()-1, yAxisBinsArray, axis2Bins.size()-1, axis2BinsArray);
            plots3D["ptD"      + histName] = new TH3D("ptD"      + histName + "3D", "p_{T}D",                     xAxisBins.size()-1, xAxisBinsArray, yAxisBins.size()-1, yAxisBinsArray, otherBins.size()-1, otherBinsArray);
            plots3D["mult"     + histName] = new TH3D("mult"     + histName + "3D", "multiplicity",               xAxisBins.size()-1, xAxisBinsArray, yAxisBins.size()-1, yAxisBinsArray, multBins.size()-1, multBinsArray);
            plots3D["qg"       + histName] = new TH3D("qg"       + histName + "3D", "quark-gluon likelihood",     xAxisBins.size()-1, xAxisBinsArray, yAxisBins.size()-1, yAxisBinsArray, otherBins.size()-1, otherBinsArray);
          }
        }
      }

      // Fill histos
      TString binName;
      while(t.next()){
        if(!bins.getBinName(binName)) 	continue;									// Find bin and return false if outside ranges 
        if(t.jetIdLevel < 3) 		continue;									// Select tight jets
        if(!t.matchedJet) 		continue; 									// Only matched jets
        if(t.bTag) 			continue;									// Anti-b tagging
        if(t.mult < 3 || t.axis2 > 99)	continue; 									// Take only jets with at least 3 particles
	if(!t.balanced) 		continue;									// Take only two leading jets with pt3 < 0.15*(pt1+pt2)
        if(t.nGenJetsInCone != 1 || t.nJetsForGenParticle != 1 || t.nGenJetsForGenParticle != 1) continue;		// Use only jets matched to exactly one gen jet and gen particle, and no other jet candidates

        TString type;
        if(t.partonId == 21) 	 	type = "gluon";
        else if(fabs(t.partonId) < 4) 	type = "quark";
        else continue;

        float qg = localQG.computeQGLikelihood(          t.pt, t.eta, t.rho, {(float) t.mult, t.ptD, t.axis2});

        TString histName = "_" + type + "_" + binName;
        plots["axis2"    + histName]->Fill(xAxisVar, yAxisVar, t.axis2,    t.weight);
        plots["ptD"      + histName]->Fill(xAxisVar, yAxisVar, t.ptD,      t.weight);
        plots["mult"     + histName]->Fill(xAxisVar, yAxisVar, t.mult, 	   t.weight);
        plots["qg"       + histName]->Fill(xAxisVar, yAxisVar, qg,         t.weight);
        if(extraBinInformation){
          plots3D["axis2"    + histName]->Fill(xAxisVar, yAxisVar, t.axis2,    t.weight);
          plots3D["ptD"      + histName]->Fill(xAxisVar, yAxisVar, t.ptD,      t.weight);
          plots3D["mult"     + histName]->Fill(xAxisVar, yAxisVar, t.mult,     t.weight);
          plots3D["qg"       + histName]->Fill(xAxisVar, yAxisVar, qg, 	       t.weight);
        }
      }

      // Stacking, cosmetics and saving
      for(auto& plot : plots){
        TCanvas c(xAndyVar);
        if(xLog) c.SetLogx();
        if(yLog) c.SetLogy();
        
        TString type = (plot.first.Contains("gluon") ? "gluon" : "quark");

        plot.second->SetTitle(TString(plot.second->GetTitle()) + " for " + type + "s;" + xVar + ";" + yVar);
        plot.second->SetStats(0);
        plot.second->SetContour(50);
        setCorrectMinimum((TH2D*) plot.second);
        if(plot.first.Contains("qg") || plot.first.Contains("ptD")){
          plot.second->SetMinimum(0.);
          plot.second->SetMaximum(1.);
        } else if(plot.first.Contains("mult")){
          plot.second->SetMinimum(0.);
          plot.second->SetMaximum(50.);
        } else if(plot.first.Contains("axis2")){
          plot.second->SetMinimum(2.);
          plot.second->SetMaximum(7.);
        }
        plot.second->Draw("COLZ");
        c.Modified();
        bins.printInfoOnPlot(plot.first, jetType);

        TString variable = plot.first(0, plot.first.First("_"));
        TString pdfDir = "./plots/profile2D/" + file + "/" + jetType + "/"  + xAndyVar + (useBins ? "" : "_nobins") + "/" + variable + "/";
        TString pdfName = pdfDir + plot.first + ".pdf";
        system("mkdir -p " + pdfDir);
        c.SaveAs(pdfName);

        // RMS, kurtosis and skewness
        if(extraBinInformation){
          TH2D* rmsInBin = plot.second->ProjectionXY("RMS " + plot.first, "C=E");
          rmsInBin->SetTitle(TString(plot.second->GetTitle()).ReplaceAll("mean","RMS") + ";" + xVar + ";" + yVar);
          rmsInBin->SetStats(0);
          setCorrectMinimum(rmsInBin);
          rmsInBin->SetContour(50);
          rmsInBin->Draw("COLZ");
          c.Modified();
          bins.printInfoOnPlot(plot.first, jetType);

          pdfName = pdfDir + plot.first + "_RMS.pdf";
          c.SaveAs(pdfName);

          TH2D* kurtosisInBin = (TH2D*) rmsInBin->Clone("kurtosis" + plot.first);
          TH2D* skewnessInBin = (TH2D*) rmsInBin->Clone("skewness" + plot.first);
          kurtosisInBin->Reset();
          skewnessInBin->Reset();
          for(int i = 1; i < xAxisBins.size(); ++i){
            for(int j = 1; j < yAxisBins.size(); ++j){
              TH1D *temp = plots3D[plot.first]->ProjectionZ("", i, i, j, j);
              if(temp->GetEntries() != 0 && temp->GetKurtosis() < 10) kurtosisInBin->SetBinContent(i, j, temp->GetKurtosis());
              if(temp->GetEntries() != 0 && temp->GetSkewness() < 10) skewnessInBin->SetBinContent(i, j, temp->GetSkewness());
              delete temp;
            }
          }
          kurtosisInBin->SetTitle(TString(plot.second->GetTitle()).ReplaceAll("mean","kurtosis"));
          skewnessInBin->SetTitle(TString(plot.second->GetTitle()).ReplaceAll("mean","skewness"));

          setCorrectMinimum(kurtosisInBin);
          kurtosisInBin->Draw("COLZ");
          c.Modified();
          bins.printInfoOnPlot(plot.first, jetType);
          pdfName = pdfDir + plot.first + "_kurtosis.pdf";
          c.SaveAs(pdfName);

          
          setCorrectMinimum(skewnessInBin);
          skewnessInBin->Draw("COLZ");
          c.Modified();
          bins.printInfoOnPlot(plot.first, jetType);
          pdfName = pdfDir + plot.first + "_skewness.pdf";
          c.SaveAs(pdfName);

          delete rmsInBin;
          delete kurtosisInBin;
          delete skewnessInBin;
        }

        // Entries/bin
        TH2D* entriesInBin = plot.second->ProjectionXY("jets/bin " + plot.first, "B");
        entriesInBin->SetTitle("jets / bin for " + type + "s;" + xVar + ";" + yVar);
        entriesInBin->SetStats(0);
        entriesInBin->SetContour(50);
        setCorrectMinimum(entriesInBin);
        entriesInBin->Draw("COLZ");
        c.Modified();

        pdfDir = "./plots/profile2D/" + file + "/" + jetType + "/" + xAndyVar + (useBins ? "" : "_nobins") + "/jetsPerBin/";
        pdfName = pdfDir + "jetsPerBin_" + type + ".pdf";
        system("mkdir -p " + pdfDir);
        c.SaveAs(pdfName);
        delete entriesInBin;
      }

      for(auto& plot : plots) delete plot.second;
      for(auto& plot : plots3D) delete plot.second;
    }
  }
}


int main(int argc, char**argv){
  if(argc == 1) return 1;
  makeProfile2D(TString(argv[1]), false, false);
  return 0;
}

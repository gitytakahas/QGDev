#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "binClass.h"
#include "binningConfigurations.h"
#include "treeLooper.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc"


// Function to switch between identification string of quark and gluon pdf
TString switchQG(TString inputBin){
  if(inputBin.Contains("gluon")) inputBin.ReplaceAll("gluon","quark");
  else                           inputBin.ReplaceAll("quark","gluon");
  return inputBin;
}

// Function to replace in string
TString replace(TString input, TString a, TString b){
  input.ReplaceAll(a,b);
  return input;
}


int main(int argc, char**argv){
  std::vector<TString> files	= {"QCD_AllPtBins"};
  //std::vector<TString> files	= {"QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14"};
  std::vector<TString> jetTypes = {"AK4chs"};

  binClass bins = getCentralPtSlices();											// This is the binning/selection of the ROC plots

  // Loop over different samples and jet types
  for(TString file : files){
    for(TString jetType : jetTypes){
      std::cout << "Making ROC curves for " << jetType << " in file " << file << "..." << std::endl;
      system("rm -rf plots/ROC/" + file + "/" + jetType);

      treeLooper t(file, jetType);											// Init tree
      bins.setReference("pt",  &t.pt);
      bins.setReference("eta", &t.eta);
      bins.setReference("rho", &t.rho);

      // Init local QGLikelihoodCalculators to compare, now we are comparing the 76X pdf's with the 80X pdf's
      std::map<TString, QGLikelihoodCalculator*> localQG;
      localQG["76X"] = new QGLikelihoodCalculator("../data/pdfQG_" + jetType + "_13TeV_76X.root");
      localQG["80X"] = new QGLikelihoodCalculator("../data/pdfQG_" + jetType + "_13TeV_80X.root");

      // Keep track of different ROC types
      std::vector<TString> rocTypes;
      rocTypes.push_back("");							// Empty string means these are ROC curves of the variables themselves
      for(auto& l : localQG) rocTypes.push_back("_" + l.first);			// These are ROC curves for each of the defined localQG

      // Creation of histos: for each type, plotting bin and pdf bin
      std::map<TString, TH1D*> plots;
      for(TString binName : bins.getAllBinNames()){
          for(TString type : {"quark","gluon"}){
            TString histName = "_" + type + "_" + binName;
            plots["axis2"  + histName] 	= new TH1D("axis2"  + histName, "axis2"     + histName, 200, 0, 8);
            plots["ptD"    + histName]	= new TH1D("ptD"    + histName, "ptD"       + histName, 200, 0, 1);
            plots["mult"   + histName]	= new TH1D("mult"   + histName, "mult"      + histName, 140, 2.5, 142.5);
            for(TString var : {"qg","axis2","ptD","mult"}){
              for(auto& l : localQG){
                TString id = var + "_" + l.first + histName;
                plots[id] = new TH1D(id, id, 100, -0.0001, 1.0001);
              }
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
        if(t.mult < 3) 			continue; 									// Need at least three particles in the jet
	if(!t.balanced) 		continue;									// Take only two leading jets with pt3 < 0.15*(pt1+pt2)
        if(t.nGenJetsInCone != 1 || t.nJetsForGenParticle != 1 || t.nGenJetsForGenParticle != 1) continue;		// Use only jets matched to exactly one gen jet and gen particle, and no other jet candidates

        TString type;
        if(t.partonId == 21) 	 	type = "gluon";
        else if(fabs(t.partonId) < 4) 	type = "quark";
        else continue;

        TString histName = "_" + type + "_" + binName;
        plots["axis2"   + histName]->Fill(t.axis2, t.weight);
        plots["ptD"     + histName]->Fill(t.ptD,   t.weight);
        plots["mult"    + histName]->Fill(t.mult,  t.weight);
        for(auto& l : localQG){
          plots["qg_"      + l.first + histName]->Fill(l.second->computeQGLikelihood(t.pt, t.eta, t.rho, {(float) t.mult, t.ptD, t.axis2}),	t.weight);
          plots["axis2_"   + l.first + histName]->Fill(l.second->computeQGLikelihood(t.pt, t.eta, t.rho, {-1, -1, t.axis2}), 			t.weight);
          plots["ptD_"     + l.first + histName]->Fill(l.second->computeQGLikelihood(t.pt, t.eta, t.rho, {-1, t.ptD}), 				t.weight);
          plots["mult_"    + l.first + histName]->Fill(l.second->computeQGLikelihood(t.pt, t.eta, t.rho, {(float) t.mult}), 			t.weight);
        }
      }
      for(auto& plot : plots) plot.second->Scale(1./plot.second->Integral(0, plot.second->GetNbinsX() + 1));

      // Stacking, cosmetics and saving
      for(auto& plot : plots){
        if(!plot.first.Contains("gluon") || !plot.first.Contains("qg" + rocTypes[1]) || plot.second->GetEntries() == 0) continue;
        TCanvas c;

        TLegend l(0.12,0.2,0.4,0.5);
        l.SetFillColor(kWhite);
        l.SetBorderSize(0);

        std::map<TString, TGraph*> roc;
        for(TString var : {"qg","axis2","ptD","mult"}){
          for(TString type : rocTypes){
            if(var.Contains("qg")  && type == "") continue;
            if(!var.Contains("qg") && type != "") continue; // skip single-variable likelihoods
            TH1D *pdfGluon = plots[replace(plot.first,"qg" + rocTypes[1], var+type)]; 
            TH1D *pdfQuark = plots[switchQG(replace(plot.first,"qg" + rocTypes[1], var+type))];
            if(!pdfQuark || !pdfGluon) continue;

            roc[var+type] = new TGraph(plot.second->GetNbinsX() + 2);
            for(int bin = 0; bin <= plot.second->GetNbinsX() + 1; ++bin){
              double gluonRej = pdfGluon->Integral(0, bin);
              double quarkEff = 1.-pdfQuark->Integral(0, bin);
              if(var+type == "mult"){ gluonRej = 1.-gluonRej; quarkEff = 1.-quarkEff;}
              roc[var+type]->SetPoint(bin, gluonRej, quarkEff);
            }

            if(var == "axis2")	roc[var+type]->SetLineColor(type != rocTypes[0] ? kGreen+4   : kYellow);
            if(var == "ptD")	roc[var+type]->SetLineColor(type != rocTypes[0] ? kMagenta+4 : kAzure+10);
            if(var == "mult")	roc[var+type]->SetLineColor(type != rocTypes[0] ? kOrange    : kRed);
            if(var == "qg")	roc[var+type]->SetLineColor(type != rocTypes[1] ? kBlack     : kGray);

            TString entryName = "quark-gluon";
            if(var == "axis2") 	entryName = "-log(#sigma_{2})";
            if(var == "ptD") 	entryName = "p_{T}D";
            if(var == "mult") 	entryName = "multiplicity";
            if(type == "_76X")	entryName += " likelihood 76X";
            if(type == "_80X")	entryName += " likelihood 80X";
            l.AddEntry(roc[var+type], entryName, "l");

            if(roc.size() == 1){
              roc[var+type]->GetYaxis()->SetTitle("quark-jet efficiency");
              roc[var+type]->GetXaxis()->SetTitle("gluon-jet rejection");
              roc[var+type]->GetXaxis()->SetRangeUser(0,1);
              roc[var+type]->GetYaxis()->SetRangeUser(0,1);
              roc[var+type]->SetTitle("ROC");
              roc[var+type]->Draw("AL");
            } else roc[var+type]->Draw("l");
          }
        }
        l.Draw();
        c.Modified();

        bins.printInfoOnPlot(plot.first, jetType);

        TString pdfDir = "./plots/ROC/" + file + "/" + jetType + "/";
        TString pdfName = pdfDir + plot.first + ".pdf";
        pdfName.ReplaceAll("_gluon","");
        pdfName.ReplaceAll("qg" + rocTypes[1],"ROC");
        system("mkdir -p " + pdfDir);
        c.SaveAs(pdfName);
        for(auto& r : roc) delete r.second;
      }

      for(auto& plot : plots) delete plot.second;
      for(auto& l : localQG) delete l.second;
    }
  }
  return 0;
}

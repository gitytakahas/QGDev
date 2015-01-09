#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TVector.h"
#include "binClass.h"
#include "binningConfigurations.h"
#include "treeLooper.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc"


int main(int argc, char**argv){
  std::vector<TString> files	= {"QCD_AllPtBins"};
//std::vector<TString> files	= {"QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14"};
  std::vector<TString> jetTypes = {"AK4","AK4chs"};

  // Define binning for pdfs
  binClass bins = getV1Binning();

  // Loop over different samples and jet types
  for(TString file : files){
    for(TString jetType : jetTypes){
      std::cout << "Making plots for " << jetType << " in file " << file << "..." << std::endl;
      system("rm -rf plots/ROC/" + file + "/" + jetType + "_");

      treeLooper t(file, jetType);						// Init tree
      bins.setReference("pt",  &t.pt);
      bins.setReference("eta", &t.eta);
      bins.setReference("rho", &t.rho);

      // Init local QGLikelihoodCalculator
      QGLikelihoodCalculator localQG("../data/pdfQG_" + jetType + "_13TeV_v1.root");

      // Creation of histos
      std::map<TString, TH1D*> plots;
      for(TString binName : bins.getAllBinNames()){
        for(TString type : {"quark","gluon"}){
          TString histName = "_" + type + "_" + binName;
          plots["axis2"  + histName] 	= new TH1D("axis2"  + histName, "axis2"     + histName, 200, 0, 8);
          plots["ptD"    + histName]	= new TH1D("ptD"    + histName, "ptD"       + histName, 200, 0, 1);
          plots["mult"   + histName]	= new TH1D("mult"   + histName, "mult"      + histName, 140, 2.5, 142.5);
          for(TString var : {"qg","axis2","ptD","mult"}){
            plots[var + "_l" + histName] = new TH1D(var + "_l" + histName, var + "_l" + histName, 100, -0.0001, 1.0001);
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
        plots["axis2"   + histName]->Fill(t.axis2, 										t.weight);
        plots["ptD"     + histName]->Fill(t.ptD, 										t.weight);
        plots["mult"    + histName]->Fill(t.mult, 										t.weight);
        plots["qg_l"    + histName]->Fill(localQGa.computeQGLikelihood(t.pt, t.eta, t.rho, {(float) t.mult, t.ptD, t.axis2}),	t.weight);
        plots["axis2_l" + histName]->Fill(localQGa.computeQGLikelihood(t.pt, t.eta, t.rho, {-1, -1, t.axis2}), 			t.weight);
        plots["ptD_l"   + histName]->Fill(localQGa.computeQGLikelihood(t.pt, t.eta, t.rho, {-1, t.ptD}), 			t.weight);
        plots["mult_l"  + histName]->Fill(localQGa.computeQGLikelihood(t.pt, t.eta, t.rho, {(float) t.mult}), 			t.weight);
      }
      for(auto& plot : plots) plot.second->Scale(1./plot.second->Integral(0, plot.second->GetNbinsX() + 1));

      // Stacking, cosmetics and saving
      for(auto& plot : plots){
        if(!plot.first.Contains("gluon") || !plot.first.Contains("qg_l") || plot.second->GetEntries() == 0) continue;
        TCanvas c;

        TLegend l(0.12,0.2,0.4,0.5);
        l.SetFillColor(kWhite);
        l.SetBorderSize(0);

        std::map<TString, TGraph*> roc;
        for(TString var : {"qg","axis2","ptD","mult"}){
          for(TString type : {"_l",""}){
            if(var.Contains("qg") && type == "") continue;

            roc[var+type] = new TGraph(plot.second->GetNbinsX() + 2);
            for(int bin = 0; bin <= plot.second->GetNbinsX() + 1; ++bin){
              TString histName = plot.first;
              histName.ReplaceAll("qga_l",var+type);
              double gluonRej = plots[histName]->Integral(0, bin);
              histName.ReplaceAll("gluon","quark");
              double quarkEff = 1.-plots[histName]->Integral(0, bin);
              if(var+type == "mult"){ gluonRej = 1.-gluonRej; quarkEff = 1.-quarkEff;}
              roc[var+type]->SetPoint(bin, gluonRej, quarkEff);
            }

            if(var == "axis2")		roc[var+type]->SetLineColor(type == "_l" ? kGreen+4   : kYellow);
            if(var == "ptD")		roc[var+type]->SetLineColor(type == "_l" ? kMagenta+4 : kAzure+10);
            if(var == "mult")		roc[var+type]->SetLineColor(type == "_l" ? kRed :       kOrange);
            if(var == "qg")		roc[var+type]->SetLineColor(kBlack);
            roc[var+type]->SetLineWidth(type == "_l" ? 3. : 1.);
            roc[var+type]->SetLineStyle(type == "" ? 3 : 1);


            TString entryName = "quark-gluon";
            if(var == "axis2") 	entryName = "-log(#sigma_{2})";
            if(var == "ptD") 	entryName = "p_{T}D";
            if(var == "mult") 	entryName = "multiplicity";
            if(type == "_l")	entryName += " likelihood";
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
        pdfName.ReplaceAll("qg_l","ROC");
        system("mkdir -p " + pdfDir);
        c.SaveAs(pdfName);
        for(auto& r : roc) delete r.second;
      }

      for(auto& plot : plots) delete plot.second;
    }
  }
  return 0;
}

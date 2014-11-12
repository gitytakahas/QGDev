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
#include "treeLooper.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc"


int main(int argc, char**argv){
  std::vector<TString> files	= {"QCD_AllPtBins"};
//std::vector<TString> files	= {"QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14"};
  std::vector<TString> jetTypes = {"AK4chs"};

  // Define binning for pdfs
  binClass bins;
  bins.setBinRange("eta", {0,1.3,1.5,2,2.5,3,4.7});
  bins.setBinRange("pt" , bins.getBins(20, 20, 2000, true, {6500}));				// i.e. 20 bins from 20 to 2000 with log=true and with an additional bin up to 6500
  bins.setBinRange("rho", {0, 9999});
  bins.printBinRanges();

  // Link some bins to be merged because of low statistics (for example higher pT bins at large eta)
  for(int i=10; i < 21; ++i) bins.setLinks("eta5_pt9_rho0", {TString::Format("eta5_pt%d_rho0",i)});
  for(int i=14; i < 21; ++i) bins.setLinks("eta4_pt13_rho0", {TString::Format("eta4_pt%d_rho0",i)});
  for(int i=16; i < 21; ++i) bins.setLinks("eta3_pt15_rho0", {TString::Format("eta3_pt%d_rho0",i)});
  for(int i=18; i < 21; ++i) bins.setLinks("eta2_pt17_rho0", {TString::Format("eta2_pt%d_rho0",i)});
  for(int i=19; i < 21; ++i) bins.setLinks("eta1_pt18_rho0", {TString::Format("eta1_pt%d_rho0",i)});

  // Loop over different samples and jet types
  for(TString file : files){
    for(TString jetType : jetTypes){
      std::cout << "Making plots for " << jetType << " in file " << file << "..." << std::endl;
      system("rm -rf plots/ROC/" + file + "/" + jetType);

      treeLooper t(file, jetType);						// Init tree
      bins.setReference("rho", &t.rho);
      bins.setReference("pt",  &t.pt);
      bins.setReference("eta", &t.eta);

      // Init local QGLikelihoodCalculator
      QGLikelihoodCalculator localQG("../data/pdfQG_" + jetType + "_13TeV.root");
      QGLikelihoodCalculator localQG_cdf("../data/pdfQG_" + jetType + "_fineBinning_13TeV.root");
      QGLikelihoodCalculator localQG2("../data/pdfQG_" + jetType + "_13TeV_testVarTransform.root", true);
      QGLikelihoodCalculator localQG2_cdf("../data/pdfQG_" + jetType + "_fineBinning_13TeV_testVarTransform.root", true);

      // Creation of histos
      std::map<TString, TH1D*> plots;
      for(TString binName : bins.getAllBinNames()){
        for(TString type : {"quark","gluon"}){
          TString histName = "_" + type + "_" + binName;
          plots["axis2"  + histName] 	= new TH1D("axis2"  + histName, "axis2"     + histName, 1000, 0, 8);
          plots["ptD"    + histName]	= new TH1D("ptD"    + histName, "ptD"       + histName, 1000, 0, 1);
          plots["mult"   + histName]	= new TH1D("mult"   + histName, "mult"      + histName, 100, 0.5, 100.5);
          for(TString var : {"qg","qg2","axis2","ptD","mult"}){
            for(TString type : {"_l","_c"}) plots[var + type + histName] = new TH1D(var + type + histName, var + type + histName, 1000, -0.0001, 1.0001);
          }
        }
      }

      // Fill histos
      while(t.next()){
        if(!bins.update()) 	continue;										// Find bin and return false if outside ranges 
        if(t.jetIdLevel < 3) 	continue;										// Select tight jets
        if(!t.matchedJet) 	continue; 										// Only matched jets
        if(t.bTag) 		continue;										// Anti-b tagging
        if(t.mult < 3) 		continue; 										// Need at least three particles in the jet
	if(!t.balanced) 	continue;										// Take only two leading jets with pt3 < 0.15*(pt1+pt2)
        if(t.nGenJetsInCone != 1 || t.nJetsForGenParticle != 1 || t.nGenJetsForGenParticle != 1) continue;		// Use only jets matched to exactly one gen jet and gen particle, and no other jet candidates

        TString type;
        if(t.partonId == 21) 	 	type = "gluon";
        else if(fabs(t.partonId) < 4) 	type = "quark";
        else continue;

        TString histName = "_" + type + "_" + bins.name;
        plots["axis2"   + histName]->Fill(t.axis2, 											t.weight);
        plots["ptD"     + histName]->Fill(t.ptD, 											t.weight);
        plots["mult"    + histName]->Fill(t.mult, 											t.weight);
        plots["qg_l"    + histName]->Fill(localQG.computeQGLikelihood(       t.pt, t.eta, t.rho, {(float) t.mult, t.ptD, t.axis2}),	t.weight);
        plots["qg2_l"    + histName]->Fill(localQG2.computeQGLikelihood(     t.pt, t.eta, t.rho, {t.axis2, (float) t.mult, t.ptD}),	t.weight);
        plots["axis2_l" + histName]->Fill(localQG.computeQGLikelihood(       t.pt, t.eta, t.rho, {-1, -1, t.axis2}), 			t.weight);
        plots["ptD_l"   + histName]->Fill(localQG.computeQGLikelihood(       t.pt, t.eta, t.rho, {-1, t.ptD}), 				t.weight);
        plots["mult_l"  + histName]->Fill(localQG.computeQGLikelihood(       t.pt, t.eta, t.rho, {(float) t.mult}), 			t.weight);
        plots["qg_c"    + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(t.pt, t.eta, t.rho, {(float) t.mult, t.ptD, t.axis2}), 	t.weight);
        plots["qg2_c"   + histName]->Fill(localQG2_cdf.computeQGLikelihoodCDF(t.pt, t.eta, t.rho, {t.axis2, (float) t.mult, t.ptD}), 	t.weight);
        plots["axis2_c" + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(t.pt, t.eta, t.rho, {-1, -1, t.axis2}), 			t.weight);
        plots["ptD_c"   + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(t.pt, t.eta, t.rho, {-1, t.ptD}), 				t.weight);
        plots["mult_c"  + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(t.pt, t.eta, t.rho, {(float) t.mult}), 			t.weight);
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
        for(TString var : {"qg","qg2","axis2","ptD","mult"}){
          for(TString type : {"_l","_c",""}){
            if(var.Contains("qg") && type == "") continue;

            roc[var+type] = new TGraph(plot.second->GetNbinsX() + 2);
            for(int bin = 0; bin <= plot.second->GetNbinsX() + 1; ++bin){
              TString histName = plot.first;
              histName.ReplaceAll("qg_l",var+type);
              double gluonRej = plots[histName]->Integral(0, bin);
              histName.ReplaceAll("gluon","quark");
              double quarkEff = 1.-plots[histName]->Integral(0, bin);
              if(var+type == "mult"){ gluonRej = 1.-gluonRej; quarkEff = 1.-quarkEff;}
              roc[var+type]->SetPoint(bin, gluonRej, quarkEff);
            }

            if(type == "_l")  	 	roc[var+type]->SetLineColor(var == "qg2" ? 30 : (var == "qg"? kGray+1 : (var == "axis2"? kGreen+4 : (var == "ptD"? kMagenta+4 : kRed))));
            else if(type == "_c")  	roc[var+type]->SetLineColor(var == "qg2" ? 3 : (var == "qg"? kBlack : (var == "axis2"? kYellow : (var == "ptD"? kAzure+10 : kOrange))));
            else		 	roc[var+type]->SetLineColor(var == "axis2"? kGreen+4 : (var == "ptD"? kMagenta+4 : kRed));
	    roc[var+type]->SetLineWidth(type == "_l" ? 3 : 1);
            roc[var+type]->SetLineStyle(type == "" ? 3 : 1);

            TString entryName = "quark-gluon";
            if(var == "axis2") 	entryName = "-log(#sigma_{2})";
            if(var == "ptD") 	entryName = "p_{T}D";
            if(var == "mult") 	entryName = "multiplicity";
            if(type == "_l")	entryName += " likelihood";
            if(type == "_c")	entryName += " CDF-likelihood";
            if(var == "qg2") 	entryName += " (uncorrelated)";
            l.AddEntry(roc[var+type], entryName, "l");

            if(var+type == "qg_l"){
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

        TLatex t;
        t.SetNDC(kTRUE);
        t.SetTextAlign(33);
        t.SetTextSize(0.02);
        t.DrawLatex(0.9,0.98,  TString::Format("%.1f < #rho < %.1f", bins.getLowerEdge(plot.first, "rho"), bins.getUpperEdge(plot.first, "rho")));
        t.DrawLatex(0.9,0.955, TString::Format("%.1f < #eta < %.1f", bins.getLowerEdge(plot.first, "eta"), bins.getUpperEdge(plot.first, "eta")));
        t.DrawLatex(0.9,0.93,  TString::Format("%.1f < p_{T} < %.1f", bins.getLowerEdge(plot.first, "pt"), bins.getUpperEdge(plot.first, "pt")));
        t.SetTextAlign(13);
        t.DrawLatex(0.1,0.93,  jetType);

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

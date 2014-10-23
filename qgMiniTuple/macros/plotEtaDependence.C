#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "binFunctions.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc"


int main(int argc, char**argv){
  bool norm = true;

  // Define binning for plots
  std::vector<float> etaBins = {0,1,1.2,1.4,1.6};
  std::vector<float> ptBins; getBins(ptBins, 20, 20, 2000, true); ptBins.push_back(6500);
  std::vector<float> rhoBins = {0,9999};
  printBins("eta", etaBins);
  printBins("pt", ptBins);
  printBins("rho", rhoBins); std::cout << std::endl;

  // For different samples and jet types
  for(TString file : {"QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14"}){
    for(TString jetType : {"AK4chs"}){
      std::cout << "Making plots for " << jetType << " in file " << file << "..." << std::endl;
      system("rm -rf plots/etaDependence/" + file + "/" + jetType);

      // Init local QGLikelihoodCalculator
      QGLikelihoodCalculator localQG("../data/pdfQG_" + jetType + "_13TeV.root");
      QGLikelihoodCalculator localQG_cdf("../data/pdfQG_" + jetType + "_fineBinning_13TeV.root");

      // Init qgMiniTuple
      TFile *qgMiniTupleFile = new TFile("~/public/merged/QGMiniTuple/qgMiniTuple_" + file + ".root");
      TTree *qgMiniTuple; qgMiniTupleFile->GetObject("qgMiniTuple"+jetType+"/qgMiniTuple",qgMiniTuple);
      float rho, pt, eta, axis2, ptD, bTag; 
      int event, mult, partonId, jetIdLevel, nGenJetsInCone, nJetsForGenParticle, nGenJetsForGenParticle;
      bool balanced, matchedJet;
      qgMiniTuple->SetBranchAddress("rho", 			&rho);
      qgMiniTuple->SetBranchAddress("nEvent", 			&event);
      qgMiniTuple->SetBranchAddress("pt", 			&pt);
      qgMiniTuple->SetBranchAddress("eta",	 		&eta);
      qgMiniTuple->SetBranchAddress("axis2", 			&axis2);
      qgMiniTuple->SetBranchAddress("ptD",	 		&ptD);
      qgMiniTuple->SetBranchAddress("mult",	 		&mult);
      qgMiniTuple->SetBranchAddress("bTag", 			&bTag);
      qgMiniTuple->SetBranchAddress("partonId", 		&partonId);
      qgMiniTuple->SetBranchAddress("jetIdLevel", 		&jetIdLevel);
      qgMiniTuple->SetBranchAddress("balanced",			&balanced);
      qgMiniTuple->SetBranchAddress("matchedJet",		&matchedJet);
      qgMiniTuple->SetBranchAddress("nGenJetsInCone",		&nGenJetsInCone);
      qgMiniTuple->SetBranchAddress("nGenJetsForGenParticle",	&nGenJetsForGenParticle);
      qgMiniTuple->SetBranchAddress("nJetsForGenParticle",	&nJetsForGenParticle);

      // Creation of histos
      std::map<TString, TH1F*> plots;
      for(int etaBin = 0; etaBin < getNBins(etaBins); ++etaBin){
        for(int ptBin = 0; ptBin < getNBins(ptBins); ++ptBin){
          for(int rhoBin = 0; rhoBin < getNBins(rhoBins); ++rhoBin){
            for(TString type : {"quark","gluon"}){
              TString histName = "_" + type + TString::Format("_eta-%d_pt-%d_rho-%d", etaBin, ptBin, rhoBin);
              plots["axis2" + histName] = new TH1F("axis2" + histName, "axis2" + histName, 100, 0, 8);
              plots["ptD"   + histName]	= new TH1F("ptD"   + histName, "ptD"   + histName, 100, 0, 1);
              plots["mult"  + histName]	= new TH1F("mult"  + histName, "mult"  + histName, 100, 0.5, 100.5);
              plots["qg"    + histName]	= new TH1F("qg"    + histName, "qg"    + histName, 100, -0.001, 1.001);
              plots["cdf"   + histName]	= new TH1F("cdf"   + histName, "cdf"   + histName, 100, -0.001, 1.001);
            }
          }
        }
      }

      // Fill histos
      for(int i = 0; i < qgMiniTuple->GetEntries(); ++i){
        qgMiniTuple->GetEntry(i);
        int rhoBin, etaBin, ptBin;
        if(!getBinNumber(rhoBins, rho, rhoBin)) 				continue;
        if(!getBinNumber(etaBins, fabs(eta), etaBin)) 				continue;
        if(!getBinNumber(ptBins, pt, ptBin)) 					continue;

        if(jetIdLevel < 3) continue;											// Select tight jets
        if(mult < 3) continue; 												//Need at least three particles in the jet
        if(!matchedJet || nGenJetsInCone != 1 || nJetsForGenParticle != 1 || nGenJetsForGenParticle != 1) continue;	// Use only jets matched to exactly one gen jet and gen particle, and no other jet candidates
        if(!balanced) continue;												// Take only two leading jets with pt3 < 0.15*(pt1+pt2)
        TString type;
        if(partonId == 21) 	 	type = "gluon";
        else if(fabs(partonId) < 4) 	type = "quark";
        else continue;

        TString histName = "_" + type + TString::Format("_eta-%d_pt-%d_rho-%d", etaBin, ptBin, rhoBin);
        plots["axis2" + histName]->Fill(axis2);
        plots["ptD"   + histName]->Fill(ptD);
        plots["mult"  + histName]->Fill(mult);
        plots["qg"    + histName]->Fill(localQG.computeQGLikelihood(       pt, eta, rho, {(float) mult, ptD, axis2}));
        plots["cdf"   + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {(float) mult, ptD, axis2}));
      }
      if(norm) for(auto& plot : plots) plot.second->Scale(1./plot.second->Integral());

      // Stacking, cosmetics and saving
      for(auto& plot : plots){
        if(!plot.first.Contains("eta-0") || !plot.first.Contains("gluon")) continue;
        TCanvas c;
        TLegend l(0.25,0.91,0.75,0.98);
        l.SetNColumns(2);
        l.SetFillColor(kWhite);
        l.SetBorderSize(0);

        TString axisTitle = "";
        if(plot.first.Contains("axis2")) axisTitle = "-log(#sigma_{2})";
        if(plot.first.Contains("ptD")) 	 axisTitle = "p_{T}D";
        if(plot.first.Contains("mult"))  axisTitle = "multiplicity";
        if(plot.first.Contains("qg"))    axisTitle = "quark-gluon likelihood";
        if(plot.first.Contains("cdf"))   axisTitle = "quark-gluon CDF-likelihood";

        bool emptyStack = true;
        THStack stack(plot.first,";"+axisTitle+";"+(norm?"fraction of ":"")+"jets/bin");
        for(TString type : {"gluon","quark"}){
          for(int etaBin = 0; etaBin < etaBins.size() - 1; ++etaBin){
            TString histName = plot.first; 
            histName.ReplaceAll("eta-0", TString::Format("eta-%d", etaBin)).ReplaceAll("gluon",type);
            if(plots[histName]->GetEntries() == 0) continue;
            else emptyStack = false;
            int color = (type == "gluon"? 30 + etaBin : 49 - etaBin);
            plots[histName]->SetLineColor(color);
            plots[histName]->SetLineWidth(2);
            l.AddEntry(plots[histName], TString::Format(type + ", %.1f < #eta < %.1f", etaBins[etaBin], etaBins[etaBin+1]), "l");
            stack.Add(plots[histName]);
          }
        }
        if(emptyStack) continue;
        stack.Draw("nostack L");
        stack.GetXaxis()->CenterTitle();
        stack.GetYaxis()->CenterTitle();
        stack.GetYaxis()->SetTitleOffset(1.3);
        stack.GetXaxis()->SetTitleSize(0.7);
        stack.GetYaxis()->SetTitleSize(0.7);
        c.Modified();
        l.Draw();

        int ptBin  = getBinFromString(plot.first, "pt");
        int rhoBin = getBinFromString(plot.first, "rho");

        TLatex t;
        t.SetNDC(kTRUE);
        t.SetTextAlign(33);
        t.SetTextSize(0.02);
        t.DrawLatex(0.9,0.955, TString::Format("%.1f < #rho < %.1f", rhoBins[rhoBin], rhoBins[rhoBin+1]));
        t.DrawLatex(0.9,0.93,  TString::Format("%.1f < p_{T} < %.1f", ptBins[ptBin], ptBins[ptBin+1]));
        t.SetTextAlign(13);
        t.DrawLatex(0.1,0.93,  jetType);

        TString variable = plot.first(0, plot.first.First("_"));
        TString pdfDir = "./plots/etaDependence/" + file + "/" + jetType + "/" + variable + "/";
        TString pdfName = pdfDir + plot.first + ".pdf";
        pdfName.ReplaceAll("_eta-0","").ReplaceAll("_gluon","");
        system("mkdir -p " + pdfDir);
        c.SaveAs(pdfName);
      }

      for(auto& plot : plots) delete plot.second;
      delete qgMiniTupleFile;
    }
  }
  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TH3D.h"
#include "TProfile2D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "binFunctions.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc"


int main(int argc, char**argv){
  bool norm = true;
  bool extraBinInformation = true;

  // Define binning for plots
  std::vector<float> etaBins; getBins(etaBins, 40, 0, 4, false); etaBins.push_back(4.2); etaBins.push_back(4.4); etaBins.push_back(4.7);
  std::vector<float> ptBins; getBins(ptBins, 20, 20, 2000, true); ptBins.push_back(6500);
  std::vector<float> rhoBins = {0,9999};
  printBins("eta", etaBins);
  printBins("pt", ptBins);
  printBins("rho", rhoBins); std::cout << std::endl;

  // Needed for 3D plots (which will be needed to get kurtosis and skewness)
  std::vector<float> axis2Bins; getBins(axis2Bins, 100, 0, 8, false);
  std::vector<float> multBins;  getBins(multBins, 100, 0.5, 100.5, false);
  std::vector<float> otherBins; getBins(otherBins, 100, -0.001, 1.001, false);

  double etaBinsArray[1000];   std::copy(etaBins.begin(), etaBins.end(), etaBinsArray);
  double ptBinsArray[1000];    std::copy(ptBins.begin(),  ptBins.end(),  ptBinsArray);
  double otherBinsArray[1000]; std::copy(otherBins.begin(),  otherBins.end(),  otherBinsArray);
  double axis2BinsArray[1000]; std::copy(axis2Bins.begin(),  axis2Bins.end(),  axis2BinsArray);
  double multBinsArray[1000];  std::copy(multBins.begin(),  multBins.end(),  multBinsArray);

  // For different samples and jet types
  for(TString file : {"QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14"}){
    for(TString jetType : {"AK4chs"}){
      std::cout << "Making plots for " << jetType << " in file " << file << "..." << std::endl;
      system("rm -rf plots/profile2D/" + file + "/" + jetType);

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
      std::map<TString, TProfile2D*> plots;
      std::map<TString, TH3D*> plots3D;
      for(TString type : {"quark","gluon"}){
        for(int rhoBin = 0; rhoBin < getNBins(rhoBins); ++rhoBin){
          TString histName = "_" + type + TString::Format("_rho-%d", rhoBin);
          plots["axis2" + histName] = new TProfile2D("axis2" + histName, "mean -log(#sigma_{2})",           getNBins(etaBins), etaBinsArray, getNBins(ptBins), ptBinsArray);
          plots["ptD"   + histName] = new TProfile2D("ptD"   + histName, "mean p_{T}D",                     getNBins(etaBins), etaBinsArray, getNBins(ptBins), ptBinsArray);
          plots["mult"  + histName] = new TProfile2D("mult"  + histName, "mean multiplicity",               getNBins(etaBins), etaBinsArray, getNBins(ptBins), ptBinsArray);
          plots["qg"    + histName] = new TProfile2D("qg"    + histName, "mean quark-gluon likelihood",     getNBins(etaBins), etaBinsArray, getNBins(ptBins), ptBinsArray);
          plots["cdf"   + histName] = new TProfile2D("cdf"   + histName, "mean quark-gluon CDF-likelihood", getNBins(etaBins), etaBinsArray, getNBins(ptBins), ptBinsArray);
          if(extraBinInformation){
            plots3D["axis2" + histName] = new TH3D("axis2" + histName + "3D", "-log(#sigma_{2})",           getNBins(etaBins), etaBinsArray, getNBins(ptBins), ptBinsArray, getNBins(axis2Bins), axis2BinsArray);
            plots3D["ptD"   + histName] = new TH3D("ptD"   + histName + "3D", "p_{T}D",                     getNBins(etaBins), etaBinsArray, getNBins(ptBins), ptBinsArray, getNBins(otherBins), otherBinsArray);
            plots3D["mult"  + histName] = new TH3D("mult"  + histName + "3D", "multiplicity",               getNBins(etaBins), etaBinsArray, getNBins(ptBins), ptBinsArray, getNBins(multBins), multBinsArray);
            plots3D["qg"    + histName] = new TH3D("qg"    + histName + "3D", "quark-gluon likelihood",     getNBins(etaBins), etaBinsArray, getNBins(ptBins), ptBinsArray, getNBins(otherBins), otherBinsArray);
            plots3D["cdf"   + histName] = new TH3D("cdf"   + histName + "3D", "quark-gluon CDF-likelihood", getNBins(etaBins), etaBinsArray, getNBins(ptBins), ptBinsArray, getNBins(otherBins), otherBinsArray);
          }
        }
      }

      // ill histos
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

        TString histName = "_" + type + TString::Format("_rho-%d", rhoBin);
        plots["axis2" + histName]->Fill(fabs(eta), pt, axis2);
        plots["ptD"   + histName]->Fill(fabs(eta), pt, ptD);
        plots["mult"  + histName]->Fill(fabs(eta), pt, mult);
        plots["qg"    + histName]->Fill(fabs(eta), pt, localQG.computeQGLikelihood(       pt, eta, rho, {(float) mult, ptD, axis2}));
        plots["cdf"   + histName]->Fill(fabs(eta), pt, localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {(float) mult, ptD, axis2}));
        if(extraBinInformation){
          plots3D["axis2" + histName]->Fill(fabs(eta), pt, axis2);
          plots3D["ptD"   + histName]->Fill(fabs(eta), pt, ptD);
          plots3D["mult"  + histName]->Fill(fabs(eta), pt, mult);
          plots3D["qg"    + histName]->Fill(fabs(eta), pt, localQG.computeQGLikelihood(       pt, eta, rho, {(float) mult, ptD, axis2}));
          plots3D["cdf"   + histName]->Fill(fabs(eta), pt, localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {(float) mult, ptD, axis2}));
        }
      }

      // Stacking, cosmetics and saving
      for(auto& plot : plots){
        TCanvas c;
        c.SetLogy();
        
        TString type = (plot.first.Contains("gluon") ? "gluon" : "quark");

        plot.second->SetTitle(TString(plot.second->GetTitle()) + " for " + type + "s;#eta;p_{T}");
        plot.second->SetStats(0);
        plot.second->SetContour(50);
        plot.second->Draw("COLZ");
        c.Modified();

        TString variable = plot.first(0, plot.first.First("_"));
        TString pdfDir = "./plots/profile2D/" + file + "/" + jetType + "/" + variable + "/";
        TString pdfName = pdfDir + plot.first + ".pdf";
        system("mkdir -p " + pdfDir);
        c.SaveAs(pdfName);

        // RMS, kurtosis and skewness
        if(extraBinInformation){
          TH2D* rmsInBin = plot.second->ProjectionXY("RMS " + plot.first, "C=E");
          rmsInBin->SetTitle(TString(plot.second->GetTitle()).ReplaceAll("mean","RMS") + ";#eta;p_{T}");
          rmsInBin->SetStats(0);
          rmsInBin->SetContour(50);
          rmsInBin->Draw("COLZ");
          c.Modified();

          pdfName = pdfDir + plot.first + "_RMS.pdf";
          c.SaveAs(pdfName);

          TH2D* kurtosisInBin = (TH2D*) rmsInBin->Clone("kurtosis" + plot.first);
          TH2D* skewnessInBin = (TH2D*) rmsInBin->Clone("skewness" + plot.first);
          for(int i = 1; i < getNBins(etaBins) + 1; ++i){
            for(int j = 1; j < getNBins(ptBins) + 1; ++j){
              TH1D *temp = plots3D[plot.first]->ProjectionZ("", i, i, j, j);
              if(temp->GetEntries() != 0){
                kurtosisInBin->SetBinContent(i, j, temp->GetKurtosis());
                skewnessInBin->SetBinContent(i, j, temp->GetSkewness());
              }
              delete temp;
            }
          }
          kurtosisInBin->SetTitle(TString(plot.second->GetTitle()).ReplaceAll("mean","kurtosis"));
          kurtosisInBin->SetTitle(TString(plot.second->GetTitle()).ReplaceAll("mean","skewness"));

          kurtosisInBin->Draw("COLZ");
          c.Modified();
          pdfName = pdfDir + plot.first + "_kurtosis.pdf";
          c.SaveAs(pdfName);

          skewnessInBin->Draw("COLZ");
          c.Modified();
          pdfName = pdfDir + plot.first + "_skewness.pdf";
          c.SaveAs(pdfName);

          delete rmsInBin;
          delete kurtosisInBin;
          delete skewnessInBin;
        }

        // Entries/bin
        TH2D* entriesInBin = plot.second->ProjectionXY("jets/bin " + plot.first, "B");
        entriesInBin->SetTitle("jets / bin for " + type + "s;#eta;p_{T}");
        entriesInBin->SetStats(0);
        entriesInBin->SetContour(50);
        entriesInBin->Draw("COLZ");
        c.Modified();

        pdfDir = "./plots/profile2D/" + file + "/" + jetType + "/jetsPerBin/";
        pdfName = pdfDir + "jetsPerBin_" + type + ".pdf";
        system("mkdir -p " + pdfDir);
        c.SaveAs(pdfName);
        delete entriesInBin;
      }

      for(auto& plot : plots) delete plot.second;
      delete qgMiniTupleFile;
    }
  }
  return 0;
}

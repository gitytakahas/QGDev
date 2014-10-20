#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TVector.h"
#include "binFunctions.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc"


int main(int argc, char**argv){
  // Define binning for plots
  std::vector<float> etaBins = {0,2.5,4.7};
  std::vector<float> ptBinsC; getBins(ptBinsC, 20, 20, 2000, true); ptBinsC.push_back(4000);
  std::vector<float> ptBinsF; getBins(ptBinsF, 20, 20, 2000, true); ptBinsF.erase(ptBinsF.end() - 12, ptBinsF.end()); ptBinsF.push_back(4000);
  std::vector<float> rhoBins = {0,50};

  printBins("eta", etaBins);
  printBins("pt (central)", ptBinsC);
  printBins("pt (forward)", ptBinsF);
  printBins("rho", rhoBins); std::cout << std::endl;

  // For different samples and jet types
//  for(TString file : {"VBF_HToBB_M-125_13TeV-powheg-pythia6","EWKZjj_mqq120_mll50_13TeV_madgraph-pythia8","QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8"}){
//  for(TString file : {"QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14"}){
  for(TString file : {"QCD_AllPtBins"}){
//    for(TString jetType : {"AK4","AK4chs","AK5","AK5chs"}){
    for(TString jetType : {"AK4chs"}){
      std::cout << "Making plots for " << jetType << " in file " << file << "..." << std::endl;
      system("rm -rf plots/ROC/" + file + "/" + jetType);
      system("mkdir -p plots/ROC/" + file + "/" + jetType);

      // Init local QGLikelihoodCalculator
      QGLikelihoodCalculator localQG("../data/pdfQG_" + jetType + "_13TeV.root");
      QGLikelihoodCalculator localQG_cdf("../data/pdfQG_" + jetType + "_fineBinning_13TeV.root");

      // Init qgMiniTuple
      TFile *qgMiniTupleFile = new TFile(file == "test" ? "../test/qgMiniTuple.root" : "~/public/merged/QGMiniTuple/qgMiniTuple_" + file + ".root");
      TTree *qgMiniTuple; qgMiniTupleFile->GetObject("qgMiniTuple"+jetType+"/qgMiniTuple",qgMiniTuple);
      float rho, pt, eta, axis2, ptD, bTag;
      int event, mult, partonId, jetIdLevel, nGenJetsInCone, nJetsForGenParticle, nGenJetsForGenParticle;
      bool balanced, matchedJet;
      qgMiniTuple->SetBranchAddress("nEvent", 			&event);
      qgMiniTuple->SetBranchAddress("rho", 			&rho);
      qgMiniTuple->SetBranchAddress("pt", 			&pt);
      qgMiniTuple->SetBranchAddress("eta", 			&eta);
      qgMiniTuple->SetBranchAddress("axis2", 			&axis2);
      qgMiniTuple->SetBranchAddress("ptD", 			&ptD);
      qgMiniTuple->SetBranchAddress("mult", 			&mult);
      qgMiniTuple->SetBranchAddress("partonId", 		&partonId);
      qgMiniTuple->SetBranchAddress("jetIdLevel",		&jetIdLevel);
      qgMiniTuple->SetBranchAddress("balanced",			&balanced);
      qgMiniTuple->SetBranchAddress("matchedJet",		&matchedJet);
      qgMiniTuple->SetBranchAddress("nGenJetsInCone",		&nGenJetsInCone);
      qgMiniTuple->SetBranchAddress("nGenJetsForGenParticle",	&nGenJetsForGenParticle);
      qgMiniTuple->SetBranchAddress("nJetsForGenParticle",	&nJetsForGenParticle);

      // Creation of histos
      std::map<TString, TH1D*> plots;
      for(int etaBin = 0; etaBin < getNBins(etaBins); ++etaBin){
        for(int ptBin = 0; ptBin < getNBins(etaBin == 0? ptBinsC : ptBinsF); ++ptBin){
          for(int rhoBin = 0; rhoBin < getNBins(rhoBins); ++rhoBin){
            for(TString type : {"quark","gluon","bquark","cquark"}){
              TString histName = "_" + type + TString::Format("_eta-%d_pt-%d_rho-%d", etaBin, ptBin, rhoBin);
              plots["axis2"  + histName] 	= new TH1D("axis2"  + histName, "axis2"     + histName, 1000, 0, 8);
              plots["ptD"    + histName]	= new TH1D("ptD"    + histName, "ptD"       + histName, 1000, 0, 1);
              plots["mult"   + histName]	= new TH1D("mult"   + histName, "mult"      + histName, 100, 0.5, 100.5);
              for(TString var : {"qg","qg2","axis2","ptD","mult"}){
                for(TString type : {"_l","_c"}) plots[var + type + histName] = new TH1D(var + type + histName, var + type + histName, 1000, -0.0001, 1.0001);
              }
            }
          }
        }
      }

      // Fill histos
      for(int i = 0; i < qgMiniTuple->GetEntries(); ++i){
        qgMiniTuple->GetEntry(i);
        int rhoBin, etaBin, ptBin;
        if(!getBinNumber(rhoBins, rho, rhoBin)) 			continue;
        if(!getBinNumber(etaBins, fabs(eta), etaBin)) 			continue;
        if(!getBinNumber(etaBin == 0? ptBinsC : ptBinsF, pt, ptBin)) 	continue;

        if(jetIdLevel < 3) 		continue;
        if(mult < 3) 			continue; 										//Need at least three particles in the jet
        TString type;
        if(partonId == 21) 	 			type = "gluon";
        else if(fabs(partonId) < 4 && partonId != 0) 	type = "quark";
        else continue;

        TString histName = "_" + type + TString::Format("_eta-%d_pt-%d_rho-%d", etaBin, ptBin, rhoBin);
        plots["axis2"   + histName]->Fill(axis2);
        plots["ptD"     + histName]->Fill(ptD);
        plots["mult"    + histName]->Fill(mult);
        plots["qg_l"    + histName]->Fill(localQG.computeQGLikelihood(       pt, eta, rho, {(float) mult, ptD, axis2}));
        plots["qg2_l"   + histName]->Fill(localQG.computeQGLikelihood(       pt, eta, rho, {(float) mult, ptD}));
        plots["axis2_l" + histName]->Fill(localQG.computeQGLikelihood(       pt, eta, rho, {-1, -1, axis2}));
        plots["ptD_l"   + histName]->Fill(localQG.computeQGLikelihood(       pt, eta, rho, {-1, ptD}));
        plots["mult_l"  + histName]->Fill(localQG.computeQGLikelihood(       pt, eta, rho, {(float) mult}));
        plots["qg_c"    + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {(float) mult, ptD, axis2}));
        plots["qg2_c"   + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {(float) mult, ptD}));
        plots["axis2_c" + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {-1, -1, axis2}));
        plots["ptD_c"   + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {-1, ptD}));
        plots["mult_c"  + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {(float) mult}));
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
            if(var == "qg2") 	entryName = "(mult, p_{T}D)";
            if(type == "_l")	entryName += " likelihood";
            if(type == "_c")	entryName += " CDF-likelihood";
            l.AddEntry(roc[var+type], entryName, "l");

            if(var+type == "qg_l"){
              roc[var+type]->GetYaxis()->SetTitle("quark-jet efficiency");
              roc[var+type]->GetXaxis()->SetTitle("gluon-jet rejection");
              roc[var+type]->GetXaxis()->SetRangeUser(0,1);
              roc[var+type]->GetYaxis()->SetRangeUser(0,1);
              roc[var+type]->SetTitle("RoC");
              roc[var+type]->Draw("AL");
            } else roc[var+type]->Draw("l");
          }
        }
        l.Draw();
        c.Modified();

        int rhoBin = getBinFromString(plot.first, "rho");
        int ptBin  = getBinFromString(plot.first, "pt");
        int etaBin = getBinFromString(plot.first, "eta");

        TLatex t;
        t.SetNDC(kTRUE);
        t.SetTextAlign(33);
        t.SetTextSize(0.02);
        t.DrawLatex(0.9,0.98,  TString::Format("%.1f < #rho < %.1f", rhoBins[rhoBin], rhoBins[rhoBin+1]));
        t.DrawLatex(0.9,0.955, TString::Format("%.1f < #eta < %.1f", etaBins[etaBin], etaBins[etaBin+1]));
        t.DrawLatex(0.9,0.93,  TString::Format("%.1f < p_{T} < %.1f", etaBin == 0? ptBinsC[ptBin] : ptBinsF[ptBin], etaBin == 0? ptBinsC[ptBin+1] : ptBinsF[ptBin+1]));
        t.SetTextAlign(13);
        t.DrawLatex(0.1,0.93,  jetType);

        TString pdfName = "./plots/ROC/" + file + "/" + jetType + "/" + plot.first + ".pdf";
        pdfName.ReplaceAll("_gluon","");
        pdfName.ReplaceAll("qg_l","ROC");
        c.SaveAs(pdfName);
        for(auto& r : roc) delete r.second;
      }

      for(auto& plot : plots) delete plot.second;
      delete qgMiniTupleFile;
    }
  }
  return 0;
}

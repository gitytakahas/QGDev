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
#include "binFunctions.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc"


int main(int argc, char**argv){
  TString qgMiniTuplesDir 	= "~tomc/public/merged/QGMiniTuple/"; // On T2B
  std::vector<TString> files	= {"QCD_AllPtBins"};
  std::vector<TString> jetTypes = {"AK4chs"};

  // To be used in case of file == QCD_AllPtBins:
  std::vector<TString> ptHatBins = {"15to30","30to50","50to80","80to120","120to170","170to300","300to470","470to600","600to800","800to1000","1000to1400","1400to1800","1800to2400", "2400to3200","3200"};
  std::vector<int> ptHatMin      = { 15,      30,      50,      80,       120,       170,       300,       470,       600,       800,        1000,        1400,        1800,         2400,        3200};
  std::vector<float>   nEvents   = { 2498841, 2449363, 2500315, 2500098,  2491398,   1490834,   1498032,   1498417,   1465278,   1500369,    1500642,     1500040,     2953210 ,     2958105,     2953431};
  std::vector<float>   xsec      = { 2237e6,  1615e5,  2211e4,  3000114,  493200,    120300,    7475,      587.1,     167,       28.25,      8.195,       0.7346,      0.102,        0.00644,     0.000163};

  // Define binning for plots
  std::vector<float> etaBins = {0,2.5,4.7};
  std::vector<float> ptBinsC; getBins(ptBinsC, 50, 20, 2500, true); ptBinsC.push_back(6500);
  std::vector<float> ptBinsF; getBins(ptBinsF, 5, 20, 200, true); ptBinsF.push_back(1000);
  std::vector<float> rhoBins = {0,9999};

  printBins("eta", etaBins);
  printBins("pt (central)", ptBinsC);
  printBins("pt (forward)", ptBinsF);
  printBins("rho", rhoBins); std::cout << std::endl;
  makeTexLoop(ptBinsC, "ptBinsC.tex");
  makeTexLoop(ptBinsF, "ptBinsF.tex");

  // Loop over different samples and jet types
  for(TString file : files){
    for(TString jetType : jetTypes){
      std::cout << "Making plots for " << jetType << " in file " << file << "..." << std::endl;
      system("rm -rf plots/ROC/" + file + "/" + jetType);

      // Init local QGLikelihoodCalculator
      QGLikelihoodCalculator localQG("../data/pdfQG_" + jetType + "_13TeV.root");
      QGLikelihoodCalculator localQG_cdf("../data/pdfQG_" + jetType + "_fineBinning_13TeV.root");

      // Init qgMiniTuple
      TChain *qgMiniTuple = new TChain("qgMiniTuple"+jetType+"/qgMiniTuple");
      if(file == "QCD_AllPtBins"){
        for(TString ptHatBin : ptHatBins) 	qgMiniTuple->Add(qgMiniTuplesDir + "/qgMiniTuple_QCD_Pt-"+ ptHatBin +"_Tune4C_13TeV_pythia8.root", -1);
      } else if(file == "test"){ 		qgMiniTuple->Add("../test/qgMiniTuple.root", -1);
      } else 					qgMiniTuple->Add(qgMiniTuplesDir + "qgMiniTuple_" + file + ".root", -1);

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
	if(!balanced) 			continue;										// Take only two leading jets with pt3 < 0.15*(pt1+pt2)
        if(!matchedJet || nGenJetsInCone != 1 || nJetsForGenParticle != 1 || nGenJetsForGenParticle != 1) continue;		// Use only jets matched to exactly one gen jet and gen particle, and no other jet candidates
        if(fabs(eta) > 2 && fabs(eta) < 3) continue;									// Don't use 2.1 < |eta| < 2.9

        TString type;
        if(partonId == 21) 	 			type = "gluon";
        else if(fabs(partonId) < 4 && partonId != 0) 	type = "quark";
        else continue;

        double weight;
        if(file == "QCD_AllPtBins"){												// Try to avoid high weights from jets with  pT >>> ptHat
          int ptIndex = 0;
          while(pt > ptHatMin[ptIndex]) ++ptIndex;
          int treeIndex = qgMiniTuple->GetTreeNumber();
          if(pt < 10*ptHatMin[treeIndex]) weight = xsec[treeIndex]/nEvents[treeIndex];
          else	                          weight = xsec[ptIndex]/nEvents[ptIndex];
        } else 		                  weight = 1.;

        TString histName = "_" + type + TString::Format("_eta-%d_pt-%d_rho-%d", etaBin, ptBin, rhoBin);
        plots["axis2"   + histName]->Fill(axis2, weight);
        plots["ptD"     + histName]->Fill(ptD, weight);
        plots["mult"    + histName]->Fill(mult, weight);
        plots["qg_l"    + histName]->Fill(localQG.computeQGLikelihood(       pt, eta, rho, {(float) mult, ptD, axis2}), weight);
        plots["qg2_l"   + histName]->Fill(localQG.computeQGLikelihood(       pt, eta, rho, {(float) mult, ptD}), 	weight);
        plots["axis2_l" + histName]->Fill(localQG.computeQGLikelihood(       pt, eta, rho, {-1, -1, axis2}), 		weight);
        plots["ptD_l"   + histName]->Fill(localQG.computeQGLikelihood(       pt, eta, rho, {-1, ptD}), 			weight);
        plots["mult_l"  + histName]->Fill(localQG.computeQGLikelihood(       pt, eta, rho, {(float) mult}), 		weight);
        plots["qg_c"    + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {(float) mult, ptD, axis2}), weight);
        plots["qg2_c"   + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {(float) mult, ptD}), 	weight);
        plots["axis2_c" + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {-1, -1, axis2}), 		weight);
        plots["ptD_c"   + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {-1, ptD}), 			weight);
        plots["mult_c"  + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {(float) mult}), 		weight);
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
              roc[var+type]->SetTitle("ROC");
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

        TString pdfDir = "./plots/ROC/" + file + "/" + jetType + "/";
        TString pdfName = pdfDir + plot.first + ".pdf";
        pdfName.ReplaceAll("_gluon","");
        pdfName.ReplaceAll("qg_l","ROC");
        system("mkdir -p " + pdfDir);
        c.SaveAs(pdfName);
        for(auto& r : roc) delete r.second;
      }

      for(auto& plot : plots) delete plot.second;
      delete qgMiniTuple;
    }
  }
  return 0;
}

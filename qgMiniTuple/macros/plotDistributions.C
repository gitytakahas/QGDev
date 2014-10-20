#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <iomanip>
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TVector.h"
#include "binFunctions.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc"

int main(int argc, char**argv){
  bool overlay = true;
  bool norm = true;

  // Define binning for plots
  std::vector<float> etaBins = {0,2.5,4.7};
  std::vector<float> ptBinsC; getBins(ptBinsC, 20, 20, 2000, true); ptBinsC.push_back(4000);
  std::vector<float> ptBinsF; getBins(ptBinsF, 20, 20, 2000, true); ptBinsF.erase(ptBinsF.end() - 12, ptBinsF.end()); ptBinsF.push_back(4000);
  std::vector<float> rhoBins = {0,9999};

  printBins("eta", etaBins);
  printBins("pt (central)", ptBinsC);
  printBins("pt (forward)", ptBinsF);
  printBins("rho", rhoBins); std::cout << std::endl;

  // For different samples and jet types
//  for(TString file : {"VBF_HToBB_M-125_13TeV-powheg-pythia6","EWKZjj_mqq120_mll50_13TeV_madgraph-pythia8","QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14"}){
  for(TString file : {"QCD_AllPtBins"}){
//    for(TString jetType : {"AK4","AK4chs","AK5","AK5chs"}){
    for(TString jetType : {"AK4chs"}){
      std::cout << "Making plots for " << jetType << " in file " << file << "..." << std::endl;
      system("rm -rf plots/distributions/" + file + "/" + jetType);

      // Init local QGLikelihoodCalculator
      QGLikelihoodCalculator localQG("../data/pdfQG_" + jetType + "_13TeV.root");
      QGLikelihoodCalculator localQG_cdf("../data/pdfQG_" + jetType + "_fineBinning_rhoBinning_13TeV.root");

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
      qgMiniTuple->SetBranchAddress("mult",	 		&mult);
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
            for(TString type : {"quark","gluon","bquark","cquark","pu","undefined"}){
              TString histName = "_" + type + TString::Format("_eta-%d_pt-%d_rho-%d", etaBin, ptBin, rhoBin);
              plots["axis2"  + histName] 	= new TH1D("axis2" + histName, "axis2" + histName, 100, 1, 9);
              plots["ptD"    + histName]	= new TH1D("ptD"   + histName, "ptD"   + histName, 100, 0, 1);
              plots["mult"   + histName]	= new TH1D("mult"  + histName, "mult"  + histName, 100, 0.5, 100.5);
              plots["qg"     + histName]	= new TH1D("qg"    + histName, "qg"    + histName, 250, -0.00001, 1.0001);
              plots["axis2L" + histName]	= new TH1D("axis2L" + histName, "axis2 (L)" + histName, 250, -0.0001, 1.0001);
              plots["ptDL"   + histName]	= new TH1D("ptDL"   + histName, "ptD (L)"   + histName, 250, -0.0001, 1.0001);
              plots["multL"  + histName]	= new TH1D("multL"  + histName, "mult (L)"  + histName, 250, -0.0001, 1.0001);
              plots["cdf"    + histName]	= new TH1D("qgC"    + histName, "qg (C)"    + histName, 250, -0.0001, 1.0001);
              plots["axis2C" + histName]	= new TH1D("axis2C" + histName, "axis2 (C)" + histName, 250, -0.0001, 1.0001);
              plots["ptDC"   + histName]	= new TH1D("ptDC"   + histName, "ptD (C)"   + histName, 250, -0.0001, 1.0001);
              plots["multC"  + histName]	= new TH1D("multC"  + histName, "mult (C)"  + histName, 250, -0.0001, 1.0001);
            }
          }
        }
      }

      // Fill histos
      for(int i = 0; i < qgMiniTuple->GetEntries(); ++i){
        qgMiniTuple->GetEntry(i);
        int rhoBin,etaBin,ptBin;
        if(!getBinNumber(rhoBins, rho, rhoBin)) continue;
        if(!getBinNumber(etaBins, fabs(eta), etaBin)) 			continue;
        if(!getBinNumber(etaBin == 0? ptBinsC : ptBinsF, pt, ptBin)) 	continue;

        if(jetIdLevel < 3) 		continue;
        if(mult < 3) 			continue; 										//Need at least three particles in the jet
        TString type;
        if(nGenJetsInCone < 1)		type = "pu";
        else if(!matchedJet)		type = "undefined";
        else if(partonId == 21) 	type = "gluon";
        else if(fabs(partonId) < 4) 	type = "quark";
        else if(fabs(partonId) == 4) 	type = "cquark";
        else if(fabs(partonId) == 5) 	type = "bquark";
        else continue;

        float qgfly 		= localQG.computeQGLikelihood(pt, eta, rho, {(float) mult, ptD, axis2});
        float qgflyMult 	= localQG.computeQGLikelihood(pt, eta, rho, {(float) mult});
        float qgflyPtD 		= localQG.computeQGLikelihood(pt, eta, rho, {-1, ptD});
        float qgflyAxis2 	= localQG.computeQGLikelihood(pt, eta, rho, {-1, -1, axis2});
        float qgflyC 		= localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {(float) mult, ptD, axis2});
        float qgflyCMult 	= localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {(float) mult});
        float qgflyCPtD 	= localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {-1, ptD});
        float qgflyCAxis2 	= localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {-1, -1, axis2});

        TString histName = "_" + type + TString::Format("_eta-%d_pt-%d_rho-%d", etaBin, ptBin, rhoBin);
        plots["axis2"  + histName]->Fill(axis2);
        plots["ptD"    + histName]->Fill(ptD);
        plots["mult"   + histName]->Fill(mult);
        plots["qg"     + histName]->Fill(qgfly);
        plots["axis2L" + histName]->Fill(qgflyAxis2);
        plots["ptDL"   + histName]->Fill(qgflyPtD);
        plots["multL"  + histName]->Fill(qgflyMult);
        plots["axis2C" + histName]->Fill(qgflyCAxis2);
        plots["ptDC"   + histName]->Fill(qgflyCPtD);
        plots["multC"  + histName]->Fill(qgflyCMult);
        plots["cdf"    + histName]->Fill(qgflyC);
      }
      if(norm) for(auto& plot : plots) plot.second->Scale(1./plot.second->Integral(0, plot.second->GetNbinsX() + 1));

      // Stacking, cosmetics and saving
      for(auto& plot : plots){
        if(!plot.first.Contains("gluon")) continue;
        TCanvas c;
        TLegend l(0.3,0.91,0.7,0.98);
        l.SetNColumns(2);
        l.SetFillColor(kWhite);
        l.SetBorderSize(0);

        TString axisTitle = "";
        if(plot.first.Contains("axis2")) axisTitle = "-log(#sigma_{2})";
        if(plot.first.Contains("ptD")) 	 axisTitle = "p_{T}D";
        if(plot.first.Contains("mult"))  axisTitle = "multiplicity";
        if(plot.first.Contains("qg"))    axisTitle = "quark-gluon likelihood";
        if(plot.first.Contains("cdf"))   axisTitle = "quark-gluon likelihood (CDF)";

        THStack stack(plot.first,";"+axisTitle+";");
        for(TString type : {"gluon","quark","bquark","cquark","pu","undefined"}){
          TString histName = plot.first; histName.ReplaceAll("gluon", type);
          if(plots[histName]->GetEntries() == 0) continue;
          int color = (type == "gluon"? 46 : (type == "quark"? 38 : (type == "bquark"? 32 : (type == "cquark" ? 42 : (type == "pu" ? 39 : 49)))));
          if(!overlay || type == "gluon" || type == "quark") plots[histName]->SetFillColor(color);
          plots[histName]->SetLineColor(color);
          plots[histName]->SetLineWidth(type == "quark" || type == "gluon" ? 3 : 2);
          if(overlay) plots[histName]->SetFillStyle(type == "quark"? 3004 : 3005);
          l.AddEntry(plots[histName], (type == "gluon"? "gluon" : (type == "quark"? "uds" : (type == "bquark"? "b" : (type == "cquark"? "c": type)))), "f");
          stack.Add(plots[histName]);
        }
        stack.Draw(overlay ? "nostack" : "");
        stack.GetXaxis()->CenterTitle();
        stack.GetYaxis()->CenterTitle();
        c.Modified();
        l.Draw();

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

        TString variable = "";
        for(TString var: {"axis2","ptD","mult","qg","cdf"}) if(plot.first.Contains(var)) variable = var;
        if(plot.first.Contains("L")) variable += "_likelihood";
        if(plot.first.Contains("C")) variable += "_cdf";
        if(plot.first.Contains("_NR")) variable += "_NR";
        if(plot.first.Contains("_R")) variable += "_R";
        if(plot.first.Contains("_b")) variable += "_balanced";
        if(plot.first.Contains("_nb")) variable += "_notbalanced";
        TString pdfName = "./plots/distributions/" + file + "/" + jetType + "/" + variable + "/" + plot.first + ".pdf";
        pdfName.ReplaceAll("_gluon","");
        system("mkdir -p plots/distributions/" + file + "/" + jetType + "/" + variable);
        c.SaveAs(pdfName);
      }

      for(auto& plot : plots) delete plot.second;
      delete qgMiniTupleFile;
    }
  }
  return 0;
}

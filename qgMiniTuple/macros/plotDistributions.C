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
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TVector.h"
#include "binFunctions.h"
#include "localQGLikelihoodCalculator.h"

int main(int argc, char**argv){
  bool overlay = true;
  bool norm = true;

  // Define binning for plots
  std::vector<float> etaBins 	= {0,2.5,4.7};
  std::vector<float> ptBinsC; getBins(ptBinsC, 10, 20, 2000, true); ptBinsC.push_back(4000);
  std::vector<float> ptBinsF; getBins(ptBinsF, 10, 20, 2000, true); ptBinsF.erase(ptBinsF.end() - 5, ptBinsF.end()); ptBinsF.push_back(4000);
  std::vector<float> rhoBins; getBins(rhoBins, 1, 4, 46, false);
  printBins("eta", etaBins);
  printBins("pt (central)", ptBinsC);
  printBins("pt (forward)", ptBinsF);
  printBins("rho", rhoBins); std::cout << std::endl;

  // For different samples and jet types
  for(TString file : {"VBF_HToBB_M-125_13TeV-powheg-pythia6","EWKZjj_mqq120_mll50_13TeV_madgraph-pythia8","QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8"}){
    for(TString jetType : {"AK4","AK4chs","AK5","AK5chs"}){
      std::cout << "Making plots for " << jetType << " in file " << file << "..." << std::endl;
      system("rm -rf plots/distributions/" + file + "/" + jetType);
      for(TString var: {"axis2","ptD","mult","qg"}) system("mkdir -p plots/distributions/" + file + "/" + jetType + "/" + var);

      // Init local QGLikelihoodCalculator
      QGLikelihoodCalculator localQG;
      localQG.init("../data/pdfQG_" + jetType + "_13TeV.root");

      // Init qgMiniTuple
      TFile *qgMiniTupleFile = new TFile(file == "test" ? "../test/qgMiniTuple.root" : "~/public/merged/QGMiniTuple/qgMiniTuple_" + file + ".root");
      TTree *qgMiniTuple; qgMiniTupleFile->GetObject("qgMiniTuple"+jetType+"/qgMiniTuple",qgMiniTuple);
      float rho = 0;
      std::vector<float> *qg 		= nullptr;
      std::vector<float> *pt 		= nullptr;
      std::vector<float> *eta 		= nullptr;
      std::vector<float> *axis2 	= nullptr;
      std::vector<float> *ptD 		= nullptr;
      std::vector<int> *mult 		= nullptr;
      std::vector<int> *partonId 	= nullptr;
      std::vector<bool> *jetIdLoose 	= nullptr;
      qgMiniTuple->SetBranchAddress("rho", 	&rho);
      qgMiniTuple->SetBranchAddress("qg", 	&qg);
      qgMiniTuple->SetBranchAddress("pt", 	&pt);
      qgMiniTuple->SetBranchAddress("eta", 	&eta);
      qgMiniTuple->SetBranchAddress("axis2", 	&axis2);
      qgMiniTuple->SetBranchAddress("ptD", 	&ptD);
      qgMiniTuple->SetBranchAddress("mult", 	&mult);
      qgMiniTuple->SetBranchAddress("partonId", &partonId);
      qgMiniTuple->SetBranchAddress("jetIdLoose", &jetIdLoose);

      // Creation of histos
      std::map<TString, TH1F*> plots;
      for(int etaBin = 0; etaBin < getNBins(etaBins); ++etaBin){
        for(int ptBin = 0; ptBin < getNBins(etaBin == 0? ptBinsC : ptBinsF); ++ptBin){
          for(int rhoBin = 0; rhoBin < getNBins(rhoBins); ++rhoBin){
            for(TString type : {"quark","gluon","bquark","cquark"}){
              TString histName = "_" + type + TString::Format("_eta-%d_pt-%d_rho-%d", etaBin, ptBin, rhoBin);
              plots["axis2" + histName] = new TH1F("axis2" + histName, "axis2" + histName, 50, 0, 8);
              plots["ptD"   + histName]	= new TH1F("ptD"   + histName, "ptD"   + histName, 50, 0, 1);
              plots["mult"  + histName]	= new TH1F("mult"  + histName, "mult"  + histName, 50, 0.5, 100.5);
              plots["qg"    + histName]	= new TH1F("qg"    + histName, "qg"    + histName, 50, 0, 1);
            }
          }
        }
      }

      // Fill histos
      for(int i = 0; i < qgMiniTuple->GetEntries(); ++i){
        qgMiniTuple->GetEntry(i);
        for(int j = 0; j < pt->size(); ++j){
          if(!jetIdLoose->at(j)) continue;
          TString type;
          if(partonId->at(j) == 21) 	 	type = "gluon";
          else if(fabs(partonId->at(j)) < 4) 	type = "quark";
          else if(fabs(partonId->at(j)) == 4) 	type = "cquark";
          else if(fabs(partonId->at(j)) == 5) 	type = "bquark";
          else continue;

          int etaBin, ptBin, rhoBin;
          if(!getBinNumber(etaBins, fabs(eta->at(j)), etaBin)) 			continue;
          if(!getBinNumber(etaBin == 0? ptBinsC : ptBinsF, pt->at(j), ptBin)) 	continue;
          if(!getBinNumber(rhoBins, rho, rhoBin)) 				continue;

          float qgcmssw = qg->at(j);
          float qgfly = localQG.computeQGLikelihood(pt->at(j), eta->at(j), rho, {(float) mult->at(j), ptD->at(j), -std::log(axis2->at(j))});

          TString histName = "_" + type + TString::Format("_eta-%d_pt-%d_rho-%d", etaBin, ptBin, rhoBin);
          plots["axis2" + histName]->Fill(-std::log(axis2->at(j)));
          plots["ptD"   + histName]->Fill(ptD->at(j));
          plots["mult"  + histName]->Fill(mult->at(j));
          plots["qg"    + histName]->Fill(qgcmssw);
        }
      }
      if(norm) for(auto& plot : plots) plot.second->Scale(1./plot.second->Integral(0, pdf.second->GetNbinsX() + 1));

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

        THStack stack(plot.first,";"+axisTitle+";jets/bin");
        for(TString type : {"bquark","cquark","quark","gluon"}){
          TString histName = plot.first; histName.ReplaceAll("gluon", type);
          int color = (type == "gluon"? 46 : (type == "quark"? 38 : (type == "bquark"? 32 : 42)));
          if(!overlay || type == "gluon" || type == "quark") plots[histName]->SetFillColor(color);
          plots[histName]->SetLineColor(color);
          plots[histName]->SetLineWidth(type == "quark" || type == "gluon" ? 3 : 1);
          if(overlay) plots[histName]->SetFillStyle(type == "quark"? 3004 : 3005);
          if(plots[histName]->GetEntries()>0) l.AddEntry(plots[histName], (type == "gluon"? "gluon" : (type == "quark"? "uds" : (type == "bquark"? "b" : "c"))), "f");
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
        for(TString var: {"axis2","ptD","mult","qg"}) if(plot.first.Contains(var)) variable = var;
        TString pdfName = "./plots/distributions/" + file + "/" + jetType + "/" + variable + "/" + plot.first + ".pdf";
        pdfName.ReplaceAll("_gluon","");
        c.SaveAs(pdfName);
      }

      for(auto& plot : plots) delete plot.second;
      delete qgMiniTupleFile;
    }
  }
  return 0;
}

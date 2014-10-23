#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TVector.h"
#include "binFunctions.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc"

int main(int argc, char**argv){
  bool overlay 			= true;
  bool norm 			= true;
  TString qgMiniTuplesDir 	= "~tomc/public/merged/QGMiniTuple/"; // On T2B
//  std::vector<TString> files	= {"QCD_Pt-1800_Tune4C_13TeV_pythia8"};
  std::vector<TString> files	= {"QCD_AllPtBins"};
//  std::vector<TString> files	= {"QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14"};
  std::vector<TString> jetTypes = {"AK4chs"};

  // To be used in case of file == QCD_AllPtBins:
  std::vector<TString> ptHatBins = {"15to30","30to50","50to80","80to120","120to170","170to300","300to470","470to600","600to800","800to1000","1000to1400","1400to1800","1800to2400", "2400to3200","3200"};
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

  // Loop over different samples and jet types
  for(TString file : files){
    for(TString jetType : jetTypes){
      std::cout << "Making plots for " << jetType << " in file " << file << "..." << std::endl;
      system("rm -rf plots/distributions/" + file + "/" + jetType);

      // Init local localQGikelihoodCalculator
      QGLikelihoodCalculator localQG(    "../data/pdfQG_" + jetType + "_13TeV.root");
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
      std::map<TString, TH2D*> plots2D;
      for(int etaBin = 0; etaBin < getNBins(etaBins); ++etaBin){
        for(int ptBin = 0; ptBin < getNBins(etaBin == 0? ptBinsC : ptBinsF); ++ptBin){
          for(int rhoBin = 0; rhoBin < getNBins(rhoBins); ++rhoBin){
            for(TString type : {"quark","gluon","bquark","cquark","pu","undefined"}){
              TString histName = "_" + type + TString::Format("_eta-%d_pt-%d_rho-%d", etaBin, ptBin, rhoBin);
              plots["axis2"  + histName] 	= new TH1D("axis2" + histName, "axis2" + histName, 100, 1, 9);
              plots["ptD"    + histName]	= new TH1D("ptD"   + histName, "ptD"   + histName, 100, 0, 1);
              plots["mult"   + histName]	= new TH1D("mult"  + histName, "mult"  + histName, 100, 0.5, 100.5);
              for(TString var : {"qg","axis2","ptD","mult"}){
                for(TString type : {"_l","_c"}) plots[var + type + histName] = new TH1D(var + type + histName, var + type + histName, 100, -0.0001, 1.0001);
              }
              for(TString var2D : {"qg-axis2","qg-ptD","qg-mult","axis2-ptD","axis2-mult","ptD-mult"}){
                plots2D[var2D + histName] =new TH2D(var2D + histName, var2D + histName, 100, -0.0001, 1.0001, 100, -0.0001,1.001);
              }
            }
          }
        }
      }

      // Fill histos
      for(int i = 0; i < qgMiniTuple->GetEntries(); ++i){
        qgMiniTuple->GetEntry(i);

        int rhoBin,etaBin,ptBin;
        if(!getBinNumber(rhoBins, rho, rhoBin)) 			continue;
        if(!getBinNumber(etaBins, fabs(eta), etaBin)) 			continue;
        if(!getBinNumber(etaBin == 0? ptBinsC : ptBinsF, pt, ptBin)) 	continue;
        if(fabs(eta) > 2 && fabs(eta) < 3) continue;										// Don't use 2 < |eta| < 3

        if(jetIdLevel < 3) 		continue;										// Tight jets
        if(mult < 3) 			continue; 										// Need at least three particles in the jet

        TString type;
        if(nGenJetsInCone < 1)			type = "pu";
        else if(!matchedJet)			type = "undefined";
        else if(!balanced)			type = "undefined";
        else if(nGenJetsInCone > 1)		type = "undefined";
        else if(nJetsForGenParticle != 1)	type = "undefined";
        else if(nGenJetsForGenParticle != 1)	type = "undefined";
        else if(partonId == 21) 		type = "gluon";
        else if(fabs(partonId) < 4) 		type = "quark";
        else if(fabs(partonId) == 4) 		type = "cquark";
        else if(fabs(partonId) == 5) 		type = "bquark";
        else 					type = "undefined";

        double weight;
        if(file == "QCD_AllPtBins"){
          int treeIndex = qgMiniTuple->GetTreeNumber();
          weight = xsec[treeIndex]/nEvents[treeIndex];
        } else weight = 1.;

        TString histName = "_" + type + TString::Format("_eta-%d_pt-%d_rho-%d", etaBin, ptBin, rhoBin);
        plots["axis2"   + histName]->Fill(axis2, 									weight);
        plots["ptD"     + histName]->Fill(ptD, 										weight);
        plots["mult"    + histName]->Fill(mult, 									weight);
        plots["qg_l"    + histName]->Fill(localQG.computeQGLikelihood(       pt, eta, rho, {(float) mult, ptD, axis2}), weight);
        plots["axis2_l" + histName]->Fill(localQG.computeQGLikelihood(       pt, eta, rho, {-1, -1, axis2}), 		weight);
        plots["ptD_l"   + histName]->Fill(localQG.computeQGLikelihood(       pt, eta, rho, {-1, ptD}), 			weight);
        plots["mult_l"  + histName]->Fill(localQG.computeQGLikelihood(       pt, eta, rho, {(float) mult}), 		weight);
        plots["qg_c"    + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {(float) mult, ptD, axis2}), weight);
        plots["axis2_c" + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {-1, -1, axis2}), 		weight);
        plots["ptD_c"   + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {-1, ptD}), 			weight);
        plots["mult_c"  + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(pt, eta, rho, {(float) mult}), 		weight);

        float qg_l 	= localQG.computeQGLikelihood(       pt, eta, rho, {(float) mult, ptD, axis2});
        float axis2_l 	= localQG.computeQGLikelihood(       pt, eta, rho, {-1, -1, axis2});
        float ptD_l 	= localQG.computeQGLikelihood(       pt, eta, rho, {-1, ptD});
        float mult_l 	= localQG.computeQGLikelihood(       pt, eta, rho, {(float) mult});
        plots2D["qg-axis2"   + histName]->Fill(qg_l, 	axis2_l, weight);
        plots2D["qg-ptD"     + histName]->Fill(qg_l, 	ptD_l, 	 weight);
        plots2D["qg-mult"    + histName]->Fill(qg_l, 	mult_l,  weight);
        plots2D["axis2-ptD"  + histName]->Fill(axis2_l, ptD_l, 	 weight);
        plots2D["axis2-mult" + histName]->Fill(axis2_l, mult_l,  weight);
        plots2D["ptD-mult"   + histName]->Fill(ptD_l, 	mult_l,  weight);
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
        if(plot.first.Contains("qg"))    axisTitle = "quark-gluon";
        if(plot.first.Contains("_l"))	 axisTitle += " likelihood";
        if(plot.first.Contains("_c"))	 axisTitle += " CDF-likelihood";

        THStack stack(plot.first,";"+axisTitle+";"+(norm?"fraction of ":"")+"jets/bin");
        double maximum = 0;
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
          if(overlay && (type == "gluon" || type == "quark") && plots[histName]->GetMaximum() > maximum) maximum = plots[histName]->GetMaximum();
        }
        if(!maximum) continue;
        stack.Draw(overlay ? "nostack" : "");
        if(overlay) stack.SetMaximum(maximum*1.1);
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
        if(plot.first.Contains("_l")) variable += "_likelihood";
        if(plot.first.Contains("_c")) variable += "_cdf";
        TString pdfDir = "./plots/distributions/" + file + "/" + jetType + "/" + variable + "/";
        TString pdfName = pdfDir + plot.first + ".pdf";
        pdfName.ReplaceAll("_gluon","");
        system("mkdir -p " + pdfDir);
        c.SaveAs(pdfName);
      }


      for(auto& plot : plots2D){
        if(!plot.first.Contains("gluon")) continue;

        TCanvas c;
        TString histName = plot.first; histName.ReplaceAll("gluon", "quark");
        plot.second->SetMarkerColor(kRed);
        plot.second->SetMarkerSize(0.04);
        plot.second->SetMarkerStyle(20);
        plots2D[histName]->SetMarkerColor(kBlue);
        plots2D[histName]->SetMarkerSize(0.04);
        plots2D[histName]->SetMarkerStyle(20);
        if(plot.first.Contains("qg")) 		plot.second->GetXaxis()->SetTitle("quark-gluon likelihood");
        else if(plot.first.Contains("axis2-")) 	plot.second->GetXaxis()->SetTitle("-log(#sigma_{2}) likelihood");
        else if(plot.first.Contains("ptD-")) 	plot.second->GetXaxis()->SetTitle("p_{T}D likelihood");
        if(plot.first.Contains("mult")) 	plot.second->GetYaxis()->SetTitle("multiplicity likelihood");
        else if(plot.first.Contains("-axis2")) 	plot.second->GetYaxis()->SetTitle("-log(#sigma_{2}) likelihood");
        else if(plot.first.Contains("-ptD")) 	plot.second->GetYaxis()->SetTitle("p_{T}D likelihood");
        plot.second->GetXaxis()->CenterTitle();
        plot.second->GetYaxis()->CenterTitle();
        plot.second->SetStats(0);
        plot.second->SetTitle("");
        plot.second->Draw("SCAT");
        plots2D[histName]->Draw("SCAT same");

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
        t.DrawLatex(0.35,0.98, "#color[2]{gluons}");
        t.DrawLatex(0.35,0.955,TString::Format("#color[2]{mean: %f}", plot.second->GetMean()));
        t.DrawLatex(0.35,0.93, TString::Format("#color[2]{RMS: %f}", plot.second->GetRMS()));
        t.DrawLatex(0.55,0.98, "#color[4]{quarks}");
        t.DrawLatex(0.55,0.955,TString::Format("#color[4]{mean: %f}", plots2D[histName]->GetMean()));
        t.DrawLatex(0.55,0.93, TString::Format("#color[4]{RMS: %f}", plots2D[histName]->GetRMS()));

        TString variables = histName(0, histName.First("_"));
        TString pdfDir = "./plots/distributions_2D/" + file + "/" + jetType + "/" + variables + "/";
        TString pdfName = pdfDir + plot.first + ".jpg";
        pdfName.ReplaceAll("_gluon","");
        system("mkdir -p " + pdfDir);
        c.SaveAs(pdfName);
      }

      for(auto& plot : plots) delete plot.second;
      for(auto& plot : plots2D) delete plot.second;
      delete qgMiniTuple;
    }
  }
  return 0;
}

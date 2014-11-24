#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TVector.h"
#include "binClass.h"
#include "treeLooper2.h"
#include "getTransform.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc"

int main(int argc, char**argv){
  bool overlay 			= true;
  bool norm 			= true;
  bool plot2D			= false;
  bool useDecorrelation		= false;

  // Define binning for plots
  binClass bins;
  bins.setBinRange("eta",	"#eta",		{0,1.3,1.5,2,2.5,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",	bins.getBins(20, 20, 2000, true, {6500}));				// i.e. 20 bins from 20 to 2000 with log=true and with an additional bin up to 6500
//  bins.setBinRange("cbjdR", {0, 0.55, 0.65,0.75,1,2,9999});
  bins.printBinRanges();

  // Link some bins to be merged because of low statistics (for example higher pT bins at large eta)
  for(int i=10; i < 21; ++i) bins.setLinks("eta5_pt9", {TString::Format("eta5_pt%d",i)});
  for(int i=14; i < 21; ++i) bins.setLinks("eta4_pt13", {TString::Format("eta4_pt%d",i)});
  for(int i=16; i < 21; ++i) bins.setLinks("eta3_pt15", {TString::Format("eta3_pt%d",i)});
  for(int i=18; i < 21; ++i) bins.setLinks("eta2_pt17", {TString::Format("eta2_pt%d",i)});
  for(int i=19; i < 21; ++i) bins.setLinks("eta1_pt18", {TString::Format("eta1_pt%d",i)});

  // For different jet types (if _antib is added bTag is applied)
  for(TString file : {"QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14"}){
    for(TString jetType : {"AK4chs"}){//,"AK5","AK5chs","AK7","AK7chs"}){
      std::cout << "Making plots for " << jetType << " in file " << file << "..." << std::endl;
      system("rm -rf plots/distributions/" + file + "/" + jetType);

      treeLooper t(file, jetType);						// Init tree
      bins.setReference("pt",  &t.pt);
      bins.setReference("eta", &t.eta);
//      bins.setReference("cbjdR", &t.closestJetdR);

      // Init local localQGikelihoodCalculator
      QGLikelihoodCalculator localQG(    "../data/pdfQG_" + jetType + "_13TeV.root");
      QGLikelihoodCalculator localQG_cdf("../data/pdfQG_" + jetType + "_fineBinning_13TeV.root");

      std::map<TString, std::vector<std::vector<double>>> decorrelationMatrices;
      if(useDecorrelation && !getMatricesFromFile(decorrelationMatrices, "../data/pdfQG_" + jetType + "_13TeV.root")) exit(1);

      // Creation of histos
      std::map<TString, TH1D*> plots;
      std::map<TString, TH2D*> plots2D;
      for(TString binName : bins.getAllBinNames()){
        for(TString type : {"quark","gluon","bquark","cquark","pu","undefined"}){
          TString histName = "_" + type + "_" + binName;
          plots["axis2"  + histName] 	= new TH1D("axis2" + histName, "axis2" + histName, 100, 1, 9);
          plots["ptD"    + histName]	= new TH1D("ptD"   + histName, "ptD"   + histName, 100, 0, 1);
          plots["mult"   + histName]	= new TH1D("mult"  + histName, "mult"  + histName, 100, 0.5, 100.5);
          plots["axis2_dR2"  + histName] 	= new TH1D("axis2_dR2" + histName, "axis2" + histName, 100, 1, 9);
          plots["ptD_dR2"    + histName]	= new TH1D("ptD_dR2"   + histName, "ptD"   + histName, 100, 0, 1);
          plots["mult_dR2"   + histName]	= new TH1D("mult_dR2"  + histName, "mult"  + histName, 100, 0.5, 100.5);
          plots["axis2_dR3"  + histName] 	= new TH1D("axis2_dR3" + histName, "axis2" + histName, 100, 1, 9);
          plots["ptD_dR3"    + histName]	= new TH1D("ptD_dR3"   + histName, "ptD"   + histName, 100, 0, 1);
          plots["mult_dR3"   + histName]	= new TH1D("mult_dR3"  + histName, "mult"  + histName, 100, 0.5, 100.5);
          for(TString var : {"qg","axis2","ptD","mult"}){
            for(TString type : {"_l","_cdf"}) plots[var + type + histName] = new TH1D(var + type + histName, var + type + histName, 100, -0.0001, 1.0001);
          }
          if(useDecorrelation){
            auto varRanges = calcRangeTransformation(decorrelationMatrices[binName], {0, 0.5, 0}, {8, 100.5, 1});
            plots["var1" + histName] = new TH1D("var1" + histName, "var1" + histName, 100, varRanges[0][0], varRanges[1][0]);
            plots["var2" + histName] = new TH1D("var2" + histName, "var2" + histName, 100, varRanges[0][1], varRanges[1][1]);
            plots["var3" + histName] = new TH1D("var3" + histName, "var3" + histName, 100, varRanges[0][2], varRanges[1][2]);
          }
          if(plot2D){
            for(TString var2D : {"qg-axis2","qg-ptD","qg-mult","axis2-ptD","axis2-mult","ptD-mult","ptD-dR2"}){
              plots2D[var2D + histName] = new TH2D(var2D + histName, var2D + histName, 100, -0.0001, 1.0001, 100, -0.0001,1.001);
            }
            plots2D["axis2-dR2" + histName] = new TH2D("axis2-dR2" + histName, "axis2-dR2" + histName, 100, 1, 9, 100, 1, 9);
            plots2D["mult-dR2" + histName] = new TH2D("mult-dR2" + histName, "mult-dR2" + histName, 100, 0.5, 100.5, 100, 0.5, 100.5);
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

        TString type;
        if(t.nGenJetsInCone < 1)		type = "pu";
        else if(!t.matchedJet)			type = "undefined";
//        else if(!t.balanced)			type = "undefined";
        else if(t.nGenJetsInCone > 1)		type = "undefined";
        else if(t.nJetsForGenParticle != 1)	type = "undefined";
        else if(t.nGenJetsForGenParticle != 1)	type = "undefined";
        else if(t.partonId == 21) 		type = "gluon";
        else if(fabs(t.partonId) < 4) 		type = "quark";
        else if(fabs(t.partonId) == 4) 		type = "cquark";
        else if(fabs(t.partonId) == 5) 		type = "bquark";
        else 					type = "undefined";

        TString histName = "_" + type + "_" + bins.name; 
        plots["axis2"   + histName]->Fill(t.axis2, 											t.weight);
        plots["axis2_dR2"   + histName]->Fill(t.axis2_dR2, 										t.weight);
        plots["axis2_dR3"   + histName]->Fill(t.axis2_dR3, 										t.weight);
        plots["ptD"     + histName]->Fill(t.ptD, 											t.weight);
        plots["ptD_dR2"     + histName]->Fill(t.ptD_dR2, 										t.weight);
        plots["ptD_dR3"     + histName]->Fill(t.ptD_dR3, 										t.weight);
        plots["mult"    + histName]->Fill(t.mult, 											t.weight);
        plots["mult_dR2"    + histName]->Fill(t.mult_dR2, 										t.weight);
        plots["mult_dR3"    + histName]->Fill(t.mult_dR3, 										t.weight);
        plots["qg_l"    + histName]->Fill(localQG.computeQGLikelihood(       t.pt, t.eta, t.rho, {(float) t.mult, t.ptD, t.axis2}),	t.weight);
        plots["axis2_l" + histName]->Fill(localQG.computeQGLikelihood(       t.pt, t.eta, t.rho, {-1, -1, t.axis2}), 			t.weight);
        plots["ptD_l"   + histName]->Fill(localQG.computeQGLikelihood(       t.pt, t.eta, t.rho, {-1, t.ptD}), 				t.weight);
        plots["mult_l"  + histName]->Fill(localQG.computeQGLikelihood(       t.pt, t.eta, t.rho, {(float) t.mult}), 			t.weight);
        plots["qg_cdf"    + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(t.pt, t.eta, t.rho, {(float) t.mult, t.ptD, t.axis2}), 	t.weight);
        plots["axis2_cdf" + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(t.pt, t.eta, t.rho, {-1, -1, t.axis2}), 			t.weight);
        plots["ptD_cdf"   + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(t.pt, t.eta, t.rho, {-1, t.ptD}), 			t.weight);
        plots["mult_cdf"  + histName]->Fill(localQG_cdf.computeQGLikelihoodCDF(t.pt, t.eta, t.rho, {(float) t.mult}), 			t.weight);

        if(plot2D){
          float qg_l 	= localQG.computeQGLikelihood(       t.pt, t.eta, t.rho, {(float) t.mult, t.ptD, t.axis2});
          float axis2_l = localQG.computeQGLikelihood(       t.pt, t.eta, t.rho, {-1, -1, t.axis2});
          float ptD_l 	= localQG.computeQGLikelihood(       t.pt, t.eta, t.rho, {-1, t.ptD});
          float mult_l 	= localQG.computeQGLikelihood(       t.pt, t.eta, t.rho, {(float) t.mult});
          plots2D["qg-axis2"   + histName]->Fill(qg_l, 		axis2_l, t.weight);
          plots2D["qg-ptD"     + histName]->Fill(qg_l, 		ptD_l, 	 t.weight);
          plots2D["qg-mult"    + histName]->Fill(qg_l, 		mult_l,  t.weight);
          plots2D["axis2-ptD"  + histName]->Fill(axis2_l, 	ptD_l, 	 t.weight);
          plots2D["axis2-mult" + histName]->Fill(axis2_l, 	mult_l,  t.weight);
          plots2D["ptD-mult"   + histName]->Fill(ptD_l, 	mult_l,  t.weight);
          plots2D["axis2-dR2"  + histName]->Fill(t.axis2, 	t.axis2_dR2,  	t.weight);
          plots2D["mult-dR2"   + histName]->Fill(t.mult, 	t.mult_dR2,  	t.weight);
          plots2D["ptD-dR2"    + histName]->Fill(t.ptD, 	t.ptD_dR2,  	t.weight);
        }

        if(useDecorrelation){
          std::vector<double> vars = {t.axis2, (double) t.mult, t.ptD};
          std::vector<double> uncorrVars = decorrelate(decorrelationMatrices[bins.name], vars); 
          plots["var1" + histName]->Fill(uncorrVars[0]);
          plots["var2" + histName]->Fill(uncorrVars[1]);
          plots["var3" + histName]->Fill(uncorrVars[2]);
        }
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
        if(plot.first.Contains("_cdf"))	 axisTitle += " CDF-likelihood";

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

        bins.printInfoOnPlot(plot.first, jetType);

        TString variable = "";
        for(TString var: {"axis2","ptD","mult","qg"}) if(plot.first.Contains(var)) variable = var;
        if(plot.first.Contains("_l")) variable += "_likelihood";
        if(plot.first.Contains("_cdf")) variable += "_cdf";
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

        bins.printInfoOnPlot(plot.first, jetType);

        TString variables = histName(0, histName.First("_"));
        TString pdfDir = "./plots/distributions_2D/" + file + "/" + jetType + "/" + variables + "/";
        TString pdfName = pdfDir + plot.first + ".pdf";
        pdfName.ReplaceAll("_gluon","");
        system("mkdir -p " + pdfDir);
        c.SaveAs(pdfName);
      }

      for(auto& plot : plots) delete plot.second;
      for(auto& plot : plots2D) delete plot.second;
    }
  }
  return 0;
}

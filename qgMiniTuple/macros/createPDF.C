#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "binClass.h"
#include "treeLooper.h"
#include "getTransform.h"


int main(int argc, char**argv){
  bool fineBinning = true;
  bool useDecorrelation = true;

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

  // For different jet types (if _antib is added bTag is applied)
  for(TString jetType : {"AK4chs","AK4chs_antib"}){//,"AK5","AK5chs","AK7","AK7chs"}){
    std::cout << "Building pdf's for " << jetType << "..." << std::endl;

    treeLooper t("QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14", jetType);						// Init tree
    bins.setReference("rho", &t.rho);
    bins.setReference("pt",  &t.pt);
    bins.setReference("eta", &t.eta);

    std::map<TString, std::vector<std::vector<double>>> decorrelationMatrices;
    if(useDecorrelation) decorrelationMatrices = getTransform(t, bins);

    // Creation of the pdfs
    std::map<TString, TH1D*> pdfs;
    for(TString binName : bins.getAllBinNames()){
      for(TString type : {"quark","gluon"}){
        TString histName = "_" + type + "_" + binName;
        pdfs["axis2" + histName] = new TH1D("axis2" + histName, "axis2" + histName, fineBinning ? 1000 : 200, 0, 8);
        pdfs["mult"  + histName] = new TH1D("mult"  + histName, "mult"  + histName, 100, 0.5, 100.5);
        pdfs["ptD"   + histName] = new TH1D("ptD"   + histName, "ptD"   + histName, fineBinning ? 1000 : 200, 0, 1);
        if(useDecorrelation){
          auto varRanges = calcRangeTransformation(decorrelationMatrices[binName], {0, 0.5, 0}, {8, 100.5, 1});
          pdfs["var1" + histName] = new TH1D("var1" + histName, "var1" + histName, fineBinning ? 1000 : 200, varRanges[0][0], varRanges[1][0]);
          pdfs["var2" + histName] = new TH1D("var2" + histName, "var2" + histName, fineBinning ? 1000 : 200, varRanges[0][1], varRanges[1][1]);
          pdfs["var3" + histName] = new TH1D("var3" + histName, "var3" + histName, fineBinning ? 1000 : 200, varRanges[0][2], varRanges[1][2]);
        }
      }
    }

    // Fill pdfs
    while(t.next()){
      if(!bins.update()) 	continue;										// Find bin and return false if outside ranges 
      if(t.jetIdLevel < 3) 	continue;										// Select tight jets
      if(!t.matchedJet) 	continue; 										// Only matched jets
      if(t.nGenJetsInCone != 1 || t.nJetsForGenParticle != 1 || t.nGenJetsForGenParticle != 1) continue;		// Use only jets matched to exactly one gen jet and gen particle, and no other jet candidates
      if((fabs(t.partonId) > 3 && t.partonId != 21)) continue; 								// Keep only udsg
      if(t.bTag) continue;												// Anti-b tagging
      if(!t.balanced) continue;												// Take only two leading jets with pt3 < 0.15*(pt1+pt2)
      TString type = (t.partonId == 21? "gluon" : "quark");								// Define q/g
      TString histName = "_" + type + "_" + bins.name;

      pdfs["axis2" + histName]->Fill(t.axis2);										// "axis2" already contains the log
      pdfs["mult"  + histName]->Fill(t.mult);
      pdfs["ptD"   + histName]->Fill(t.ptD);

      if(useDecorrelation){
        std::vector<double> vars = {t.axis2, (double) t.mult, t.ptD};
        std::vector<double> uncorrVars = decorrelate(decorrelationMatrices[bins.name], vars); 
        pdfs["var1" + histName]->Fill(uncorrVars[0]);
        pdfs["var2" + histName]->Fill(uncorrVars[1]);
        pdfs["var3" + histName]->Fill(uncorrVars[2]);
      }
    }

    // Make file and write binnings
    TFile *pdfFile = new TFile("../data/pdfQG_"+jetType + (fineBinning ? "_fineBinning":"") + "_13TeV_testWithVarTransform.root","RECREATE");
    pdfFile->cd();
    bins.writeBinsToFile();

    if(useDecorrelation){
      pdfFile->mkdir("decorrelationMatrices");
      pdfFile->cd("decorrelationMatrices");
      writeMatricesToFile(decorrelationMatrices, bins);
      pdfFile->cd();
    }

    // Write pdf's
    for(TString var : {"axis2","ptD","mult"}) pdfFile->mkdir(var);
    if(useDecorrelation) for(TString var : {"var1","var2","var3"}) pdfFile->mkdir(var);
    for(auto& pdf : pdfs){
      for(TString var: {"axis2","ptD","mult","var1","var2","var3"}) if(pdf.first.Contains(var)) pdfFile->cd(var);
      if(pdf.second->GetEntries() == 0) std::cout << "Warning: no entries in " << pdf.first << std::endl;		// Give warning for empty pdfs

      if(!fineBinning){													// Try to make pdf more stable:
        TString theOther = pdf.first;											// A = maximum consecutive empty bins between low_mean - RMS <--> high_mean + RMS
        if(theOther.Contains("gluon")) theOther.ReplaceAll("gluon","quark");						// B = maximum consecutive bins with < 0.6*peak between first and last bin with > 0.8*peak
        else theOther.ReplaceAll("quark","gluon");									// Rebin with factor > max(A,B)
        float thisMean 	= pdf.second->GetMean();
        float otherMean = pdfs[theOther]->GetMean();
        float thisRMS 	= pdf.second->GetRMS();
        float otherRMS 	= pdfs[theOther]->GetRMS();
        int leftBin 	= pdf.second->FindBin(thisMean < otherMean ? thisMean - thisRMS : otherMean - otherRMS) - 1;
        int rightBin 	= pdf.second->FindBin(thisMean > otherMean ? thisMean + thisRMS : otherMean + otherRMS) + 1;
        int peakBinLeft = 0;
        int peakBinRight = 0;
        for(int bin = 1; bin < pdf.second->GetNbinsX(); ++bin){
          if(pdf.second->GetBinContent(bin) > pdf.second->GetMaximum()*0.8){
            if(!peakBinLeft) peakBinLeft = bin;
            else peakBinRight = bin;
          }
        }
        int emptyBins = 0;
        int maxEmptyBins = 0;
        for(int bin = 1; bin < pdf.second->GetNbinsX(); ++bin){
          if(bin >= leftBin && bin <= rightBin && pdf.second->GetBinContent(bin) <= 0) emptyBins++;
          else if(bin > peakBinLeft && bin < peakBinRight && pdf.second->GetBinContent(bin) <= pdf.second->GetMaximum()*0.6) emptyBins++;
          else {
            if(emptyBins > maxEmptyBins) maxEmptyBins = emptyBins;
            emptyBins = 0;
          }
        }
        if(emptyBins > maxEmptyBins) maxEmptyBins = emptyBins;
        if(maxEmptyBins > 24) std::cout << "This bin is really empty: " << pdf.first << std::endl;
        else if(maxEmptyBins > 19){ pdf.second->Rebin(25); std::cout << "Rebinned (25): " << pdf.first << std::endl;}
        else if(maxEmptyBins > 9){  pdf.second->Rebin(20); std::cout << "Rebinned (20): " << pdf.first << std::endl;}
        else if(maxEmptyBins > 4){  pdf.second->Rebin(10); std::cout << "Rebinned (10): " << pdf.first << std::endl;}
        else if(maxEmptyBins > 3){  pdf.second->Rebin(5);  std::cout << "Rebinned (5): " << pdf.first << std::endl;}
        else if(maxEmptyBins > 1){  pdf.second->Rebin(4);  std::cout << "Rebinned (4): " << pdf.first << std::endl;}
        else if(maxEmptyBins > 0){  pdf.second->Rebin(2);  std::cout << "Rebinned (2): " << pdf.first << std::endl;}
      }

      pdf.second->Scale(1./pdf.second->Integral(0, pdf.second->GetNbinsX() + 1));					// Scale to integral=1 (also include underflow/overflow)
      pdf.second->Write();

      TString thisBin = pdf.first;
      for(TString i : {"gluon","quark","axis2","mult","ptD"}) thisBin.ReplaceAll(i + "_","");
      for(auto i : bins.getLinkedBins(thisBin)){									// Store copies for merged bins
        TString copyBin = pdf.first;
        copyBin.ReplaceAll(thisBin, i);
        pdf.second->Write(copyBin);
      }
    }

    for(auto& pdf : pdfs) delete pdf.second;
    for(auto& file : {pdfFile}){ file->Close(); delete file;}
  }
  return 0;
}

#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "binClass.h"
#include "binningConfigurations.h"
#include "treeLooper.h"

void rebin(TH1* hist, int rebinFactor){
  hist->Rebin(rebinFactor);
  std::cout << std::left << std::setw(20) << TString::Format("Rebinned (%d):", rebinFactor);
  std::cout << std::left << std::setw(40) << hist->GetTitle() << "(entries: " << hist->GetEntries() << ")" << std::endl;
}


int main(int argc, char**argv){
 for(bool fineBinning : {false}){

  TString binning = "v1";

  // Define binning for pdfs
  binClass bins;
  if(binning == "defaultBinning") 		bins = getDefaultBinning();
  if(binning == "8TeVBinning") 			bins = get8TeVBinning();
  if(binning == "smallEtaBinning") 		bins = getSmallEtaBinning();
  if(binning == "v1") 				bins = getV1Binning();

  // For different jet types (if _antib is added bTag is applied)
  for(TString jetType : {"AK4","AK4_antib","AK4chs","AK4chs_antib"}){
    std::cout << "Building pdf's for " << jetType << "..." << std::endl;

    treeLooper t("QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14", jetType);						// Init tree
    bins.setReference("pt",  &t.pt);
    bins.setReference("eta", &t.eta);
    bins.setReference("rho", &t.rho);

    // Creation of the pdfs
    std::map<TString, TH1D*> pdfs;
    for(TString binName : bins.getAllBinNames()){
      for(TString type : {"quark","gluon"}){
        TString histName = "_" + type + "_" + binName;
        pdfs["axis2" + histName] = new TH1D("axis2" + histName, "axis2" + histName, fineBinning ? 1000 : 200, 0, 8);
        pdfs["mult"  + histName] = new TH1D("mult"  + histName, "mult"  + histName, 140, 0.5, 140.5);
        pdfs["ptD"   + histName] = new TH1D("ptD"   + histName, "ptD"   + histName, fineBinning ? 1000 : 200, 0, 1);
      }
    }


    // Fill pdfs
    TString binName;
    while(t.next()){
      if(!bins.getBinName(binName)) 	continue;									// Find bin and return false if outside ranges
      if(t.jetIdLevel < 3) 		continue;									// Select tight jets
      if(!t.matchedJet) 	 	continue; 									// Only matched jets
      if(t.nGenJetsInCone != 1 || t.nJetsForGenParticle != 1 || t.nGenJetsForGenParticle != 1) continue;		// Use only jets matched to exactly one gen jet and gen particle, and no other jet candidates
      if((fabs(t.partonId) > 3 && t.partonId != 21)) continue; 								// Keep only udsg
      if(t.bTag) continue;												// Anti-b tagging (onluy if jetType.Contains("antib")
      if(!t.balanced) continue;												// Take only two leading jets with pt3 < 0.15*(pt1+pt2)
      if(binning == "8TeVBinning" && fabs(t.eta) > 2 && fabs(t.eta) < 3) continue;					// 8 TeV binning didn't use the intermediate
      TString type = (t.partonId == 21? "gluon" : "quark");								// Define q/g
      TString histName = "_" + type + "_" + binName;

      pdfs["axis2" + histName]->Fill(t.axis2);										// "axis2" already contains the log
      pdfs["mult"  + histName]->Fill(t.mult);
      pdfs["ptD"   + histName]->Fill(t.ptD);
    }

    // Make file and write binnings
    TFile *pdfFile = new TFile("../data/pdfQG_"+jetType + (fineBinning ? "_fineBinning":"") + "_13TeV_" + binning + "_newTest.root","RECREATE");
    pdfFile->cd();
    bins.writeBinsToFile();
    bins.writeWeightsToFile(pdfFile);

    // Write pdf's
    for(TString var : {"axis2","ptD","mult"}) pdfFile->mkdir(var);
    for(auto& pdf : pdfs){
      for(TString var: {"axis2","ptD","mult"}) if(pdf.first.Contains(var)) pdfFile->cd(var);
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
        if(maxEmptyBins > 19) std::cout << "This bin is really empty:" << std::endl;
        if(maxEmptyBins > 9)      rebin(pdf.second, 20);
        else if(maxEmptyBins > 4) rebin(pdf.second, 10);
        else if(maxEmptyBins > 3) rebin(pdf.second, 5);
        else if(maxEmptyBins > 1) rebin(pdf.second, 4);
        else if(maxEmptyBins > 0) rebin(pdf.second, 2);
      }

      pdf.second->Scale(1./pdf.second->Integral(0, pdf.second->GetNbinsX() + 1));					// Scale to integral=1 (also include underflow/overflow)
      pdf.second->SetTitle(pdf.first);
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
 }
 return 0;
}

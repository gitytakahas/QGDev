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

// Function to rebin a histogram, and print information
void rebin(TH1* hist, int rebinFactor){
  hist->Rebin(rebinFactor);
  std::cout << std::left << std::setw(20) << TString::Format("Rebinned (%d):", rebinFactor);
  std::cout << std::left << std::setw(40) << hist->GetTitle() << "(entries: " << hist->GetEntries() << ")" << std::endl;
}

// Function to switch between identification string of quark and gluon pdf
TString switchQG(TString inputBin){
  if(inputBin.Contains("gluon")) inputBin.ReplaceAll("gluon","quark");
  else                           inputBin.ReplaceAll("quark","gluon");
  return inputBin;
}

// Main function to create the pdf's
int main(int argc, char**argv){
  TString version = "v1_PU40bx50";

  // Define binning for pdfs (details and more options in binningConfigurations.h)
  binClass bins;
  if(version.Contains("v1")) 	bins = getV1Binning();
  else return;

  // For different jet types (if _antib is added bTag is applied)
  for(TString jetType : {"AK4","AK4_antib","AK4chs","AK4chs_antib"}){
    std::cout << "Building pdf's for " << jetType << "..." << std::endl;

    treeLooper t("QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14", jetType);						// Init tree (third argument is the directory path, if other than default in treeLooper.h)
    bins.setReference("pt",  &t.pt);											// Give the binning class a pointer to the variables used to bin in
    bins.setReference("eta", &t.eta);
    bins.setReference("rho", &t.rho);

    // Creation of the pdfs
    std::map<TString, TH1D*> pdfs;
    for(TString binName : bins.getAllBinNames()){
      for(TString type : {"quark","gluon"}){
        TString histName = "_" + type + "_" + binName;
        pdfs["axis2" + histName] = new TH1D("axis2" + histName, "axis2" + histName, 200, 0, 8);
        pdfs["mult"  + histName] = new TH1D("mult"  + histName, "mult"  + histName, 140, 2.5, 142.5);
        pdfs["ptD"   + histName] = new TH1D("ptD"   + histName, "ptD"   + histName, 200, 0, 1);
      }
    }


    // Fill pdfs
    TString binName;
    while(t.next()){
      if(!bins.getBinName(binName)) 								continue;		// Find bin and return false if outside ranges
      if(t.jetIdLevel < 3) 									continue;		// Select tight jets
      if(!t.matchedJet) 								 	continue; 		// Only matched jets
      if(t.nGenJetsInCone != 1 || t.nJetsForGenParticle != 1 || t.nGenJetsForGenParticle != 1) 	continue;		// Use only jets matched to exactly one gen jet and gen particle, and no other jet candidates
      if((fabs(t.partonId) > 3 && t.partonId != 21)) 						continue; 		// Keep only udsg
      if(t.bTag) 										continue;		// Anti-b tagging (always false if jetType does not contain "antib")
      if(!t.balanced) 										continue;		// Take only two leading jets with pt3 < 0.15*(pt1+pt2)  (surpresses small radiated jets with pt <<< pthat)
      if(t.mult < 3)										continue;		// Avoid jets with less than 3 particles (otherwise axis2=0)
      TString type = (t.partonId == 21? "gluon" : "quark");								// Define q/g
      TString histName = "_" + type + "_" + binName;

      pdfs["axis2" + histName]->Fill(t.axis2);										// "axis2" already contains the log
      pdfs["mult"  + histName]->Fill(t.mult);
      pdfs["ptD"   + histName]->Fill(t.ptD);
    }


    // Store the mean and RMS of the original histogram (because they could be changed by rebinning operations)
    std::map<TString, float> mean;
    std::map<TString, float> rms;
    for(auto& pdf : pdfs){
      if(pdf.second->GetEntries() == 0){ std::cout << "Error: no entries in " << pdf.first << std::endl; exit(1);}	// Force to exit when no entries in pdfs: the binning configuration should be altered to avoid this
      mean[pdf.first] = pdf.second->GetMean();
      rms[pdf.first]  = pdf.second->GetRMS();
    }


    // Check "smoothness" of the pdf, and rebin if needed
    for(auto& pdf : pdfs){
      bool isBelow      = (mean[pdf.first] < mean[switchQG(pdf.first)]);						// Define region between low(meanQ, meanG) - RMS <--> high(meanQ, meanG) + RMS
      TString low	= isBelow? pdf.first : switchQG(pdf.first);							// Most events will be within those borders, so we should not allow empty bins here
      TString high	= isBelow? switchQG(pdf.first) : pdf.first;							// (an empty bin for 1 of the three variables already results in L = 0 or 1)
      int leftBin	= pdf.second->FindBin(mean[low] - rms[low]) - 1;
      int rightBin	= pdf.second->FindBin(mean[high] + rms[high]) + 1;
     
      int leftBin2      = 0;												// Define region of the peak: moste extreme bins which exceed 80% of the maximum
      int rightBin2     = 0;												// We will check for bins within this region which go below 60%
      for(int bin = 1; bin < pdf.second->GetNbinsX(); ++bin){								// In such cases the fluctuations are really high and a larger bin width is preferred
        if(pdf.second->GetBinContent(bin) > pdf.second->GetMaximum()*0.8){						// (Maybe the thresholds could still be optimized a bit more though)
          if(!leftBin2) leftBin2  = bin;
          else          rightBin2 = bin;
        }
      }

      int emptyBins     = 0;
      int maxEmptyBins  = 0;
      for(int bin = 1; bin < pdf.second->GetNbinsX(); ++bin){
        if(     bin >= leftBin  && bin <= rightBin  && pdf.second->GetBinContent(bin) <= 0)                            ++emptyBins;
        else if(bin >= leftBin2 && bin <= rightBin2 && pdf.second->GetBinContent(bin) <= pdf.second->GetMaximum()*0.6) ++emptyBins;
        else {
          if(emptyBins > maxEmptyBins) maxEmptyBins = emptyBins;
          emptyBins = 0;
        }
      }
      if(emptyBins > maxEmptyBins) maxEmptyBins = emptyBins;

      if(maxEmptyBins > 19)     std::cout << "This pdf has a lot of emtpy bins and fluctuations:" << std::endl;
      if(maxEmptyBins > 9)      rebin(pdf.second, 20);
      else if(maxEmptyBins > 4) rebin(pdf.second, 10);
      else if(maxEmptyBins > 3) rebin(pdf.second, 5);
      else if(maxEmptyBins > 1) rebin(pdf.second, 4);
      else if(maxEmptyBins > 0) rebin(pdf.second, 2);
    }


    // Use right normalization, and try to average out leftover empty bins (originally applied in QGLikelihoodCalculator.cc but more efficient to them already here):
    for(auto& pdf : pdfs){
      pdf.second->Scale(1./pdf.second->Integral(0, pdf.second->GetNbinsX() + 1));					// Scale to integral=1 (also include underflow/overflow)

      TH1* temp = (TH1*) pdf.second->Clone();
      for(int i = 1; i < pdf.second->GetNbinsX() + 1; ++i){								// Do not consider underflow/overflow
        float content = temp->GetBinContent(i);
        float width   = temp->GetBinWidth(i);
        int j = 1;
        if(temp->Integral(0, i) > 0 && temp->Integral(i, temp->GetNbinsX()+1) > 0){					// Don't average on the edges of the pdf
          while(content <= 0 && i-j > 0 && i+j <= pdf.second->GetNbinsX()){
            content += temp->GetBinContent(i-j) + temp->GetBinContent(i+j);
            width   += temp->GetBinWidth(i-j)   + temp->GetBinWidth(i+j);
            ++j;
          }
        }
        pdf.second->SetBinContent(i, content/width);
      }
      delete temp;
    }


    // Even after rebinning and averaging the empty bins out, we still could have some empty bins at the extreme edges
    // If both quark and gluon pdf have an empty bin, we assign extreme values
    // (relative though, so it does not dominate the pdf's when we want to inspect them in the ROOT file; 0.000001/0.000999 has the same effect as 0.001/0.999)
    for(auto& pdf : pdfs){
      for(int i = 0; i <= pdf.second->GetNbinsX() + 1; ++i){
        if(pdf.second->GetBinContent(i) <= 0 && pdfs[switchQG(pdf.first)]->GetBinContent(i) <= 0){
          bool isBelow = (mean[pdf.first] < mean[switchQG(pdf.first)]);
          if(isBelow) pdf.second->SetBinContent(i, 0.000001);
          else        pdf.second->SetBinContent(i, 0.000999);
        }
      }
    }


    // Apply the likelihood weight
    for(auto& pdf : pdfs){
      TString thisBin = pdf.first(pdf.first.Index(TRegexp("eta")), pdf.first.Length());
      TString thisVar = pdf.first(0, pdf.first.Index(TRegexp("_")));
      for(int i = 0; i < pdf.second->GetNbinsX() + 1; ++i){
        pdf.second->SetBinContent(i, std::pow(pdf.second->GetBinContent(i), bins.getWeight(thisBin, (thisVar == "mult"? 0 : (thisVar == "ptD" ? 1 : 2 )))));
      }
    }


    // Make file and write binnings
    TFile *pdfFile = new TFile("../data/pdfQG_"+jetType + "_13TeV_" + version + ".root","RECREATE");
    pdfFile->cd();
    bins.writeBinsToFile();

    // Write to file
    for(TString var : {"axis2","ptD","mult"}) pdfFile->mkdir(var);
    for(auto& pdf : pdfs){
      for(TString var: {"axis2","ptD","mult"}) if(pdf.first.Contains(var)) pdfFile->cd(var);
      pdf.second->SetTitle(pdf.first);
      pdf.second->Write();

      TString thisBin = pdf.first(pdf.first.Index(TRegexp("eta")), pdf.first.Length());
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

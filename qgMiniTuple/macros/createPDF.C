#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
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
  TString version = "80X";

  // Define binning for pdfs (details and more options in binningConfigurations.h)
  binClass bins;
  if(version.Contains("v2")) 		bins = getV2Binning();
  if(version.Contains("80X")) 		bins = get76XBinning();
  else return 1;

  // For different jet types (if _antib is added bTag is applied)
  for(TString jetType : {"AK4chs","AK4chs_antib"}){ //,"AK4","AK4_antib"}){
    std::cout << "Building pdf's for " << jetType << "..." << std::endl;

    treeLooper t("QCD_AllPtBins", jetType);										// Init tree (third argument is the directory path, if other than default in treeLooper.h)
    bins.setReference("pt",  &t.pt);											// Give the binning class a pointer to the variables used to bin in
    bins.setReference("eta", &t.eta);
    bins.setReference("rho", &t.rho);

    // Creation of the pdfs
    std::map<TString, TH1D*> pdfs;
    for(TString binName : bins.getAllBinNames()){
      for(TString type : {"quark","gluon"}){
        TString histName = "_" + type + "_" + binName;
        pdfs["axis2" + histName] = new TH1D("axis2" + histName, "axis2" + histName, 100, 0, 8);				// Has been 200 bins before, but seemed to have a bit too much fluctuations still
        pdfs["mult"  + histName] = new TH1D("mult"  + histName, "mult"  + histName, 140, 2.5, 142.5);
        pdfs["ptD"   + histName] = new TH1D("ptD"   + histName, "ptD"   + histName, 100, 0, 1);				// Also 200 before
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

      pdfs["axis2" + histName]->Fill(t.axis2, t.weight);								// "axis2" already contains the log
      pdfs["mult"  + histName]->Fill(t.mult,  t.weight);
      pdfs["ptD"   + histName]->Fill(t.ptD,   t.weight);
    }

    // Try to add statistics from neighbours (first make copy, so you don't get an iterative effect)
    std::map<TString, TH1D*> pdfsCopy;
    for(auto& pdf : pdfs) pdfsCopy[pdf.first] = (TH1D*) pdf.second->Clone(pdf.first + "clone");
    for(TString binName : bins.getAllBinNames()){
      for(TString var : {"axis2","mult","ptD"}){
        for(TString neighbour : bins.getNeighbourBins(binName, var)){							// If neighbours are defined: add their statistics
          for(TString type : {"quark","gluon"}){
            pdfs[var + "_" + type + "_" + binName]->Add(pdfsCopy[var + "_" + type + "_" + neighbour]);
          }
        }
      }
    }
    for(auto& copy : pdfsCopy) delete copy.second;

    // Store the mean and RMS of the original histogram (because they could be changed by rebinning operations)
    std::map<TString, float> mean;
    std::map<TString, float> rms;
    for(auto& pdf : pdfs){
      if(pdf.second->GetEntries() == 0){ std::cout << "Error: no entries in " << pdf.first << std::endl; exit(1);}	// Force to exit when no entries in pdfs: the binning configuration should be altered to avoid this
      mean[pdf.first] = pdf.second->GetMean();
      rms[pdf.first]  = pdf.second->GetRMS();
    }


    // Check "smoothness" of the pdf: if fluctuations seem really big, we do a rebinning
    std::map<TString, int> rebinFactor;
    for(auto& pdf : pdfs){
      bool isBelow      = (mean[pdf.first] < mean[switchQG(pdf.first)]);						// Define region between low(meanQ, meanG) - RMS <--> high(meanQ, meanG) + RMS
      TString low	= isBelow? pdf.first : switchQG(pdf.first);							// Most events will be within those borders, so we should not allow empty bins here
      TString high	= isBelow? switchQG(pdf.first) : pdf.first;							// (an empty bin for 1 of the three variables already results in L = 0 or 1)
      int leftBin	= pdf.second->FindBin(mean[low] - rms[low]) - 1;
      int rightBin	= pdf.second->FindBin(mean[high] + rms[high]) + 1;
     
      int leftBin2      = 0;												// Define region of the peak: most extreme bins which exceed 80% of the maximum
      int rightBin2     = 0;												// We will check for bins within this region which go below 70%
      for(int bin = 1; bin < pdf.second->GetNbinsX(); ++bin){								// In such cases the fluctuations are really high and a larger bin width is preferred
        if(pdf.second->GetBinContent(bin) > pdf.second->GetMaximum()*0.8){						// (Maybe the thresholds could still be optimized a bit more though)
          if(!leftBin2) leftBin2  = bin;
          else          rightBin2 = bin;
        }
      }


      int leftBin3      = 0;												// Define region of the peak: most extreme bins which exceed 90% of the maximum
      int rightBin3     = 0;												// We will check for bins within this region which go below 80%
      for(int bin = 1; bin < pdf.second->GetNbinsX(); ++bin){								// In such cases the fluctuations are really high and a larger bin width is preferred
        if(pdf.second->GetBinContent(bin) > pdf.second->GetMaximum()*0.9){						// (Maybe the thresholds could still be optimized a bit more though)
          if(!leftBin3) leftBin3  = bin;
          else          rightBin3 = bin;
        }
      }

      int emptyBins     = 0;
      int maxEmptyBins  = 0;
      for(int bin = 1; bin < pdf.second->GetNbinsX(); ++bin){
        if(     bin >= leftBin  && bin <= rightBin  && pdf.second->GetBinContent(bin) <= 0)                            ++emptyBins;
        else if(bin >= leftBin2 && bin <= rightBin2 && pdf.second->GetBinContent(bin) <= pdf.second->GetMaximum()*0.7) ++emptyBins;
        else if(bin >= leftBin3 && bin <= rightBin3 && pdf.second->GetBinContent(bin) <= pdf.second->GetMaximum()*0.8) ++emptyBins;
        else {
          if(emptyBins > maxEmptyBins) maxEmptyBins = emptyBins;
          emptyBins = 0;
        }
      }
      rebinFactor[pdf.first] = std::max(maxEmptyBins, emptyBins);
    }

    for(auto& pdf : pdfs) rebinFactor[pdf.first] = std::max(rebinFactor[pdf.first], rebinFactor[switchQG(pdf.first)]);	// Use same rebin factor in quark and gluon pdf (otherwise bias if second derivate of pdf is non-zero)
    
    for(auto& pdf : pdfs){
      if(rebinFactor[pdf.first] > 19)     std::cout << "This pdf has a lot of emtpy bins and fluctuations:" << std::endl;
      if(rebinFactor[pdf.first] > 9)      rebin(pdf.second, 20);
      else if(rebinFactor[pdf.first] > 4) rebin(pdf.second, 10);
      else if(rebinFactor[pdf.first] > 3) rebin(pdf.second, 5);
      else if(rebinFactor[pdf.first] > 1) rebin(pdf.second, 4);
      else if(rebinFactor[pdf.first] > 0) rebin(pdf.second, 2);
    }


    // Normalization of the pdf's
    for(auto& pdf : pdfs) pdf.second->Scale(1./pdf.second->Integral(0, pdf.second->GetNbinsX() + 1));			// Scale to integral=1 (also include underflow/overflow)


    // Try to average out leftover fluctuations and empty bins
    for(auto& pdf : pdfs){
      if(pdf.first.Contains("gluon")) continue;

      TH1* tempQ = (TH1*) pdf.second->Clone();
      TH1* tempG = (TH1*) pdfs[switchQG(pdf.first)]->Clone();
      for(int i = 1; i < pdf.second->GetNbinsX() + 1; ++i){								// Do not consider underflow/overflow
        float contentQ = tempQ->GetBinContent(i);
        float contentG = tempG->GetBinContent(i);
        float width    = tempQ->GetBinWidth(i);

        if((1.5*contentQ < tempQ->GetBinContent(i-1) && 1.5*contentQ < tempQ->GetBinContent(i+1)) ||			// Try to average out some extreme fluctuations (i.e. only allow a difference of max 50% between two neighbouring bins)
           (1.5*contentG < tempG->GetBinContent(i-1) && 1.5*contentG < tempG->GetBinContent(i+1)) ||
           (contentQ < 1.5*tempQ->GetBinContent(i-1) && contentQ < 1.5*tempQ->GetBinContent(i+1)) ||
           (contentG < 1.5*tempG->GetBinContent(i-1) && contentG < 1.5*tempG->GetBinContent(i+1))){
          if(i-1 > 0 && i+1 < pdf.second->GetNbinsX() + 1){
            contentQ += tempQ->GetBinContent(i-1) + tempQ->GetBinContent(i+1);
            contentG += tempG->GetBinContent(i-1) + tempG->GetBinContent(i+1);
            width    += tempQ->GetBinWidth(i-1)   + tempQ->GetBinWidth(i+1);
          }
        }

        int j = 1;
        while(contentQ <= 0 || contentG <= 0){										// Average empty bins
          if(tempQ->Integral(0, i-j) <= 0) break;									// but not when surpassing the extreme edges of the pdf (see next part)
          if(tempG->Integral(0, i-j) <= 0) break;
          if(tempQ->Integral(i+j, tempQ->GetNbinsX()+1) <= 0) break;
          if(tempG->Integral(i+j, tempG->GetNbinsX()+1) <= 0) break;
          if(i-j == 0) break;
          if(i+j == pdf.second->GetNbinsX()) break;
          contentQ += tempQ->GetBinContent(i-j) + tempQ->GetBinContent(i+j);
          contentG += tempG->GetBinContent(i-j) + tempG->GetBinContent(i+j);
          width    += tempQ->GetBinWidth(i-j)   + tempQ->GetBinWidth(i+j);
          ++j;
        }

        pdf.second->SetBinContent(i, contentQ/width);
        pdfs[switchQG(pdf.first)]->SetBinContent(i, contentG/width);
      }
      delete tempQ;
      delete tempG;
    }


    // Now there are still empty bins left on the edges of the pdf, for which we assign extreme values to avoid a 0
    // (relative though, so it does not dominate the pdf's when we want to inspect them in the ROOT file; 0.000001/0.000999 has the same effect as 0.001/0.999)
    for(auto& pdf : pdfs){
      if(pdf.first.Contains("gluon")) continue;
      for(int i = 0; i <= pdf.second->GetNbinsX() + 1; ++i){
        if(pdf.second->GetBinContent(i) <= 0 || pdfs[switchQG(pdf.first)]->GetBinContent(i) <= 0){
          bool isBelow = (mean[pdf.first] < mean[switchQG(pdf.first)]);
          if(isBelow == pdf.second->GetBinCenter(i) < mean[pdf.first]){
            pdf.second->SetBinContent(i, 0.000999);
            pdfs[switchQG(pdf.first)]->SetBinContent(i, 0.000001);
          } else {
            pdf.second->SetBinContent(i, 0.000001);
            pdfs[switchQG(pdf.first)]->SetBinContent(i, 0.000999);
          }
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
    TFile *pdfFile = new TFile("pdfQG_"+jetType + "_13TeV_" + version + ".root","RECREATE");
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

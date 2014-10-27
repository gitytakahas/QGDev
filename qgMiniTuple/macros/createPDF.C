#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "binFunctions.h"

// Function to merge bins
void addTo(std::multimap<TString, std::vector<std::vector<int>>>& mymultimap, TString variableAndType, std::vector<std::vector<int>> binList){
  for(auto& i : mymultimap){
    if(i.first == variableAndType && i.second[0] == binList[0]){
      i.second.insert(i.second.end(), binList.begin() + 1, binList.end());
      return;
    }
  }
  mymultimap.insert(std::pair<TString, std::vector<std::vector<int>>>(variableAndType, binList));
}


// Main program
int main(int argc, char**argv){
  bool fineBinning = true;

  // Define binning for pdfs
  std::vector<float> etaBins = {0,1.3,1.5,2,2.5,3,4.7};
  std::vector<float> ptBins; getBins(ptBins, 20, 20, 2000, true); ptBins.push_back(6500);
  std::vector<float> rhoBins = {0,9999};

  // Check binning
  printBins("rho", rhoBins);
  printBins("eta", etaBins);
  printBins("pt", ptBins); std::cout << std::endl;

  // Link some bins to be merged (order eta - pt - rho) because of low statistics (for example higher pT bins at large eta)
  // Should be redefined when changes are made to eta-pt-rho grid
  std::multimap<TString, std::vector<std::vector<int>>> associatedBins;
  for(TString type : {"gluon","quark"}){
    for(TString var : {"axis2","mult","ptD"}){
      for(int i=10; i <= 21; ++i) addTo(associatedBins, var + "_" + type, {{5,9,0}, {5,i,0}});
      for(int i=14; i <= 21; ++i) addTo(associatedBins, var + "_" + type, {{4,13,0}, {4,i,0}});
      for(int i=16; i <= 21; ++i) addTo(associatedBins, var + "_" + type, {{3,15,0}, {3,i,0}});
      for(int i=18; i <= 21; ++i) addTo(associatedBins, var + "_" + type, {{2,17,0}, {2,i,0}});
      for(int i=19; i <= 21; ++i) addTo(associatedBins, var + "_" + type, {{1,18,0}, {1,i,0}});
    }
  }

  // For different jet types
  system("rm -r ./plots/mergedBins/");
  for(TString jetType : {"AK4chs","AK4chs_antib"}){//,"AK5","AK5chs","AK7","AK7chs"}){
    std::cout << "Building pdf's for " << jetType << "..." << std::endl;
    TString treePath = "qgMiniTuple"+jetType+"/qgMiniTuple";
    treePath.ReplaceAll("_antib","");

    // Init qgMiniTuple
    TFile *qgMiniTupleFile = new TFile("~tomc/public/merged/QGMiniTuple/qgMiniTuple_QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14.root");
    TTree *qgMiniTuple; qgMiniTupleFile->GetObject(treePath, qgMiniTuple);
    float rho, pt, eta, axis2, ptD, bTag; 
    int event, mult, partonId, jetIdLevel, nGenJetsInCone, nJetsForGenParticle, nGenJetsForGenParticle;
    bool balanced, matchedJet;
    qgMiniTuple->SetBranchAddress("rho", 			&rho);
    qgMiniTuple->SetBranchAddress("nEvent", 			&event);
    qgMiniTuple->SetBranchAddress("pt", 			&pt);
    qgMiniTuple->SetBranchAddress("eta",	 		&eta);
    qgMiniTuple->SetBranchAddress("axis2", 			&axis2);
    qgMiniTuple->SetBranchAddress("ptD",	 		&ptD);
    qgMiniTuple->SetBranchAddress("mult",	 		&mult);
    qgMiniTuple->SetBranchAddress("bTag", 			&bTag);
    qgMiniTuple->SetBranchAddress("partonId", 			&partonId);
    qgMiniTuple->SetBranchAddress("jetIdLevel",	 		&jetIdLevel);
    qgMiniTuple->SetBranchAddress("balanced",			&balanced);
    qgMiniTuple->SetBranchAddress("matchedJet",			&matchedJet);
    qgMiniTuple->SetBranchAddress("nGenJetsInCone",		&nGenJetsInCone);
    qgMiniTuple->SetBranchAddress("nGenJetsForGenParticle",	&nGenJetsForGenParticle);
    qgMiniTuple->SetBranchAddress("nJetsForGenParticle",	&nJetsForGenParticle);

    // Creation of the pdfs
    std::map<TString, TH1D*> pdfs;
    for(int etaBin = 0; etaBin < getNBins(etaBins); ++etaBin){
      for(int ptBin = 0; ptBin < getNBins(ptBins); ++ptBin){
        for(int rhoBin = 0; rhoBin < getNBins(rhoBins); ++rhoBin){
          for(TString type : {"quark","gluon"}){
            TString histName = "_" + type + TString::Format("_eta%d_pt%d_rho%d", etaBin, ptBin, rhoBin);
            pdfs["axis2" + histName] = new TH1D("axis2" + histName, "axis2" + histName, fineBinning ? 1000 : 200, 0, 8);
            pdfs["ptD"   + histName] = new TH1D("ptD"   + histName, "ptD"   + histName, fineBinning ? 1000 : 200, 0, 1);
            pdfs["mult"  + histName] = new TH1D("mult"  + histName, "mult"  + histName, 100, 0.5, 100.5);
          }
        }
      }
    }

    // Fill pdfs
//  std::map<int, std::map<int, std::map<int, std::vector<int>*>>> eventsWithRho;
//  int countDoubles = 0;
//  int lastEvent = 0;
//  bool skipEvent = false;
    for(int i = 0; i < qgMiniTuple->GetEntries(); ++i){
      qgMiniTuple->GetEntry(i);
      int rhoBin, etaBin, ptBin;
      if(!getBinNumber(rhoBins, rho, rhoBin)) 		continue;
      if(!getBinNumber(etaBins, fabs(eta), etaBin)) 	continue;
      if(!getBinNumber(ptBins, pt, ptBin)) 		continue;

//    if(lastEvent != event){												// A bit more complicated and less efficient because of our switch to flat trees
//      lastEvent = event;
//      if(eventAlreadyInRhoBin(eventsWithRho, rhoBin, event)){ ++countDoubles; skipEvent = true;}
//      else skipEvent = false;
//    }
//    if(skipEvent) continue;

      if(jetIdLevel < 3) continue;											// Select tight jets
      if(!matchedJet || nGenJetsInCone != 1 || nJetsForGenParticle != 1 || nGenJetsForGenParticle != 1) continue;	// Use only jets matched to exactly one gen jet and gen particle, and no other jet candidates
      if((fabs(partonId) > 3 && partonId != 21) || partonId == 0) continue;						// Keep only udsg
      if(jetType.Contains("antib") && bTag >  0.244) continue;								// Anti-b tagging
      if(!balanced) continue;												// Take only two leading jets with pt3 < 0.15*(pt1+pt2)
      TString type = (partonId == 21? "gluon" : "quark");								// Define q/g

      TString histName = "_" + type + TString::Format("_eta%d_pt%d_rho%d", etaBin, ptBin, rhoBin);

      pdfs["axis2" + histName]->Fill(axis2);										// "axis2" already contains the log
      pdfs["ptD"   + histName]->Fill(ptD);
      pdfs["mult"  + histName]->Fill(mult);
    }
//  for(auto& i : eventsWithRho) for(auto& j : i.second) for(auto& k : j.second) delete k.second;
//  std::cout << countDoubles << " doubles on a total of " << qgMiniTuple->GetEntries() << std::endl;

    // Merging bins as defined above and storing a copy in each of them (could be checked by eye using ./plots/mergedBins)
    for(auto i : associatedBins){
      std::cout << "Merging ";
      TCanvas c;
      TH1D* temp = nullptr;
      for(auto j : i.second){
        TString histName = i.first + TString::Format("_eta%d_pt%d_rho%d", j[0], j[1], j[2]);
        std::cout << histName << ", ";
        if(!pdfs[histName]) continue;
        if(!temp) temp = (TH1D*) pdfs[histName]->Clone();
        else temp->Add(pdfs[histName]);
      }
      temp->Scale(1./temp->Integral(0, temp->GetNbinsX() + 1));
      temp->SetLineColor(kRed);
      temp->Draw("L");
      std::cout << "..." << std::endl;
      for(auto j : i.second){
        TString histName = i.first + TString::Format("_eta%d_pt%d_rho%d", j[0], j[1], j[2]);
        if(pdfs[histName]){
          pdfs[histName]->Scale(1./pdfs[histName]->Integral(0, pdfs[histName]->GetNbinsX() + 1));
          pdfs[histName]->DrawCopy("same L");
          delete pdfs[histName];
        }
        pdfs[histName] = (TH1D*) temp->Clone(histName);
      }
      temp->Draw("same L");
      TString pdfDir = "./plots/mergedBins/" + jetType + "/";
      TString pdfName = pdfDir + temp->GetName() + ".pdf";
      system("mkdir -p " + pdfDir);
      c.SaveAs(pdfName);
      delete temp;
    }

    // Make file and write binnings
    TFile *pdfFile = new TFile("../data/pdfQG_"+jetType + (fineBinning ? "_fineBinning":"") + "_13TeV.root","RECREATE");
    pdfFile->cd();
    writeBinsToFile(rhoBins, "rhoBins");
    writeBinsToFile(etaBins, "etaBins");
    writeBinsToFile(ptBins, "ptBins");

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
    }

    for(auto& pdf : pdfs) delete pdf.second;
    for(auto& file : {pdfFile, qgMiniTupleFile}){ file->Close(); delete file;}
  }
  return 0;
}

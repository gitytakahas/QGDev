#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "binFunctions.h"


int main(int argc, char**argv){
  bool fineBinning = false;

  // Define binning for pdfs
  std::vector<float> etaBins = {0,2.5,4.7};
  std::vector<float> ptBinsC; getBins(ptBinsC, 20, 20, 2000, true); ptBinsC.push_back(4000);
  std::vector<float> ptBinsF; getBins(ptBinsF, 20, 20, 2000, true); ptBinsF.erase(ptBinsF.end() - 12, ptBinsF.end()); ptBinsF.push_back(4000);
  std::vector<float> rhoBins; getBins(rhoBins, 40, 0, 40, false); rhoBins.push_back(42); rhoBins.push_back(44); rhoBins.push_back(50);

  printBins("eta", etaBins);
  printBins("pt (central)", ptBinsC);
  printBins("pt (forward)", ptBinsF);
  printBins("rho", rhoBins); std::cout << std::endl;

  // For different jet types
  for(TString jetType : {"AK4","AK4chs","AK4chs_antib","AK5","AK5chs","AK7","AK7chs"}){
    std::cout << "Building pdf's for " << jetType << "..." << std::endl;
    TString treePath = "qgMiniTuple"+jetType+"/qgMiniTuple";
    treePath.ReplaceAll("_antib","");

    // Init qgMiniTuple
    TFile *qgMiniTupleFile = new TFile("~tomc/public/merged/QGMiniTuple/qgMiniTuple_QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_AllPU.root");
    TTree *qgMiniTuple; qgMiniTupleFile->GetObject(treePath, qgMiniTuple);
    float rho; int event;
    std::vector<float> *pt 		= nullptr;
    std::vector<float> *eta 		= nullptr;
    std::vector<float> *axis2 		= nullptr;
    std::vector<float> *ptD 		= nullptr;
    std::vector<float> *bTag 		= nullptr;
    std::vector<int> *mult 		= nullptr;
    std::vector<int> *partonId 		= nullptr;
    std::vector<bool> *jetIdLoose 	= nullptr;
    qgMiniTuple->SetBranchAddress("rho", 	&rho);
    qgMiniTuple->SetBranchAddress("nEvent", 	&event);
    qgMiniTuple->SetBranchAddress("pt", 	&pt);
    qgMiniTuple->SetBranchAddress("eta", 	&eta);
    qgMiniTuple->SetBranchAddress("axis2", 	&axis2);
    qgMiniTuple->SetBranchAddress("ptD", 	&ptD);
    qgMiniTuple->SetBranchAddress("mult", 	&mult);
    qgMiniTuple->SetBranchAddress("bTag", 	&bTag);
    qgMiniTuple->SetBranchAddress("partonId", 	&partonId);
    qgMiniTuple->SetBranchAddress("jetIdLoose", &jetIdLoose);

    // Creation of the pdfs
    std::map<TString, TH1F*> pdfs;
    for(int etaBin = 0; etaBin < getNBins(etaBins); ++etaBin){
      for(int ptBin = 0; ptBin < getNBins(etaBin == 0? ptBinsC : ptBinsF); ++ptBin){
        for(int rhoBin = 0; rhoBin < getNBins(rhoBins); ++rhoBin){
          for(TString type : {"quark","gluon"}){
            TString histName = "_" + type + TString::Format("_eta-%d_pt-%d_rho-%d", etaBin, ptBin, rhoBin);
            pdfs["axis2" + histName] = new TH1F("axis2" + histName, "axis2" + histName, fineBinning ? 1000 : 200, 0, 8);
            pdfs["ptD"   + histName] = new TH1F("ptD"   + histName, "ptD"   + histName, fineBinning ? 1000 : 200, 0, 1);
            pdfs["mult"  + histName] = new TH1F("mult"  + histName, "mult"  + histName, 100, 0.5, 100.5);
          }
        }
      }
    }

    // Fill pdfs
    std::map<int, std::map<int, std::map<int, std::vector<int>*>>> eventsWithRho;
    int countDoubles = 0;
    for(int i = 0; i < qgMiniTuple->GetEntries(); ++i){
      qgMiniTuple->GetEntry(i);
      int rhoBin;
      if(!getBinNumber(rhoBins, rho, rhoBin)) continue;
      if(eventAlreadyInRhoBin(eventsWithRho, rhoBin, event)){ ++countDoubles; continue;}
      for(int j = 0; j < pt->size(); ++j){
        if(!jetIdLoose->at(j)) continue;
        if(fabs(partonId->at(j)) > 3 && partonId->at(j) != 21) continue;						// Keep only udsg
        if(jetType.Contains("antib") && bTag->at(j) > 0.423) continue;
        TString type = (partonId->at(j) == 21? "gluon" : "quark");							// Define q/g

        int etaBin, ptBin;												// Calculate bin numbers
        if(!getBinNumber(etaBins, fabs(eta->at(j)), etaBin)) 			continue;
        if(!getBinNumber(etaBin == 0? ptBinsC : ptBinsF, pt->at(j), ptBin)) 	continue;

        if(fabs(eta->at(j)) > 2 && fabs(eta->at(j)) < 3) 			continue;				// Don't use 2 < |eta| < 3
        TString histName = "_" + type + TString::Format("_eta-%d_pt-%d_rho-%d", etaBin, ptBin, rhoBin);
        pdfs["axis2" + histName]->Fill(-std::log(axis2->at(j)));							// QGTagger uses -log(axis2) as pdf
        pdfs["ptD"   + histName]->Fill(ptD->at(j));
        pdfs["mult"  + histName]->Fill(mult->at(j));
      }
    }
    for(auto& i : eventsWithRho) for(auto& j : i.second) for(auto& k : j.second) delete k.second;
    std::cout << countDoubles << " doubles on a total of " << qgMiniTuple->GetEntries() << std::endl;

    // Write pdfs and binning to file
    TFile *pdfFile = new TFile("../data/pdfQG_"+jetType + (fineBinning ? "_fineBinning":"") + "_13TeV.root","RECREATE");
    pdfFile->cd();
    writeBinsToFile(etaBins, "etaBins");
    writeBinsToFile(ptBinsC, "ptBinsC");
    writeBinsToFile(ptBinsF, "ptBinsF");
    writeBinsToFile(rhoBins, "rhoBins");

    for(TString var : {"axis2","ptD","mult"}) pdfFile->mkdir(var);
    for(auto& pdf : pdfs){
      for(TString var: {"axis2","ptD","mult"}) if(pdf.first.Contains(var)) pdfFile->cd(var);
      if(pdf.second->GetEntries() == 0) 	std::cout << "Warning: no entries in " << pdf.first << std::endl;	// Give warning for empty pdfs

      if(!fineBinning){													// Make smooth until no empty bins between low_mean - RMS <--> high_mean + RMS
        TString theOther = pdf.first;
        if(theOther.Contains("gluon")) theOther.ReplaceAll("gluon","quark");
        else theOther.ReplaceAll("quark","gluon");
        float thisMean 	= pdf.second->GetMean();
        float otherMean = pdfs[theOther]->GetMean();
        float thisRMS 	= pdf.second->GetRMS();
        float otherRMS 	= pdfs[theOther]->GetRMS();
        int leftBin 	= pdf.second->FindBin(thisMean < otherMean ? thisMean - thisRMS : otherMean - otherRMS) - 1;
        int rightBin 	= pdf.second->FindBin(thisMean > otherMean ? thisMean + thisRMS : otherMean + otherRMS) + 1;
        int emptyBins = 0;
        int maxEmptyBins = 0;
        for(int bin = leftBin; bin <= rightBin; ++bin){
          if(pdf.second->GetBinContent(bin) <= 0) emptyBins++;
          else {
            if(emptyBins > maxEmptyBins) maxEmptyBins = emptyBins;
            emptyBins = 0;
          }
        }
        if(emptyBins > maxEmptyBins) maxEmptyBins = emptyBins;
        if(maxEmptyBins > 24) std::cout << "This bin is really empty: " << pdf.first << std::endl;
        else if(maxEmptyBins > 19) pdf.second->Rebin(25);
        else if(maxEmptyBins > 9) pdf.second->Rebin(20);
        else if(maxEmptyBins > 4) pdf.second->Rebin(10);
        else if(maxEmptyBins > 1) pdf.second->Rebin(5);
        else if(maxEmptyBins > 0) pdf.second->Rebin(2);
      }

      pdf.second->Scale(1./pdf.second->Integral(0, pdf.second->GetNbinsX() + 1));					// Scale to integral=1 (also include underflow/overflow)
      pdf.second->Write();
    }

    for(auto& pdf : pdfs) delete pdf.second;
    for(auto& file : {pdfFile, qgMiniTupleFile}){ file->Close(); delete file;}
  }
  return 0;
}

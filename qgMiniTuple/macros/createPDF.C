#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "binFunctions.h"


int main(int argc, char**argv){
  // Define binning for pdfs
  std::vector<float> etaBins = {0,2.5,4.7};
  std::vector<float> ptBinsC; getBins(ptBinsC, 20, 20, 2000, true); ptBinsC.push_back(4000);
  std::vector<float> ptBinsF; getBins(ptBinsF, 20, 20, 2000, true); ptBinsF.erase(ptBinsF.end() - 12, ptBinsF.end()); ptBinsF.push_back(4000);
  std::vector<float> rhoBins; getBins(rhoBins, 42, 4, 46, false);

  printBins("eta", etaBins);
  printBins("pt (central)", ptBinsC);
  printBins("pt (forward)", ptBinsF);
  printBins("rho", rhoBins); std::cout << std::endl;

  // For different jet types
  for(TString jetType : {"AK4","AK4chs","AK5","AK5chs","AK7","AK7chs"}){
    std::cout << "Building pdf's for " << jetType << "..." << std::endl;

    // Init qgMiniTuple
    TFile *qgMiniTupleFile = new TFile("~tomc/public/merged/QGMiniTuple/qgMiniTuple_QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8.root");
    TTree *qgMiniTuple; qgMiniTupleFile->GetObject("qgMiniTuple"+jetType+"/qgMiniTuple",qgMiniTuple);
    float rho = 0;
    std::vector<float> *pt 	= nullptr;
    std::vector<float> *eta 	= nullptr;
    std::vector<float> *axis2 	= nullptr;
    std::vector<float> *ptD 	= nullptr;
    std::vector<int> *mult 	= nullptr;
    std::vector<int> *partonId 	= nullptr;
    qgMiniTuple->SetBranchAddress("rho", 	&rho);
    qgMiniTuple->SetBranchAddress("pt", 	&pt);
    qgMiniTuple->SetBranchAddress("eta", 	&eta);
    qgMiniTuple->SetBranchAddress("axis2", 	&axis2);
    qgMiniTuple->SetBranchAddress("ptD", 	&ptD);
    qgMiniTuple->SetBranchAddress("mult", 	&mult);
    qgMiniTuple->SetBranchAddress("partonId", 	&partonId);

    // Creation of the pdfs
    std::map<TString, TH1F*> pdfs;
    for(int etaBin = 0; etaBin < getNBins(etaBins); ++etaBin){
      for(int ptBin = 0; ptBin < getNBins(etaBin == 0? ptBinsC : ptBinsF); ++ptBin){
        for(int rhoBin = 0; rhoBin < getNBins(rhoBins); ++rhoBin){
          for(TString type : {"quark","gluon"}){
            TString histName = "_" + type + TString::Format("_eta-%d_pt-%d_rho-%d", etaBin, ptBin, rhoBin);
            pdfs["axis2" + histName] = new TH1F("axis2" + histName, "axis2" + histName, 100, 0, 8);
            pdfs["ptD"   + histName] = new TH1F("ptD"   + histName, "ptD"   + histName, 100, 0, 1);
            pdfs["mult"  + histName] = new TH1F("mult"  + histName, "mult"  + histName, 100, 0.5, 100.5);
          }
        }
      }
    }

    // Fill pdfs
    for(int i = 0; i < qgMiniTuple->GetEntries(); ++i){
      qgMiniTuple->GetEntry(i);
      for(int j = 0; j < pt->size(); ++j){
        if(fabs(partonId->at(j)) > 3 && partonId->at(j) != 21) continue;						// Keep only udsg
        TString type = (partonId->at(j) == 21? "gluon" : "quark");							// Define q/g

        int etaBin, ptBin, rhoBin;											// Calculate bin numbers
        if(!getBinNumber(etaBins, fabs(eta->at(j)), etaBin)) 			continue;
        if(!getBinNumber(etaBin == 0? ptBinsC : ptBinsF, pt->at(j), ptBin)) 	continue;
        if(!getBinNumber(rhoBins, rho, rhoBin)) 				continue;
        if(fabs(eta->at(j)) > 2 && fabs(eta->at(j)) < 3) 			continue;				// Don't use 2 < |eta| < 3
        TString histName = "_" + type + TString::Format("_eta-%d_pt-%d_rho-%d", etaBin, ptBin, rhoBin);
        pdfs["axis2" + histName]->Fill(-std::log(axis2->at(j)));							// QGTagger uses -log(axis2) as pdf
        pdfs["ptD"   + histName]->Fill(ptD->at(j));
        pdfs["mult"  + histName]->Fill(mult->at(j));
      }
    }

    // Write pdfs and binning to file
    TFile *pdfFile = new TFile("../data/pdfQG_"+jetType +"_13TeV.root","RECREATE");
    pdfFile->cd();
    writeBinsToFile(etaBins, "etaBins");
    writeBinsToFile(ptBinsC, "ptBinsC");
    writeBinsToFile(ptBinsF, "ptBinsF");
    writeBinsToFile(rhoBins, "rhoBins");

    for(TString var : {"axis2","ptD","mult"}) pdfFile->mkdir(var);
    for(auto& pdf : pdfs){
      for(TString var: {"axis2","ptD","mult"}) if(pdf.first.Contains(var)) pdfFile->cd(var);
      if(pdf.second->GetEntries() == 0) 	std::cout << "Warning: no entries in " << pdf.first << std::endl;	// Give warning for empty pdfs
      else if(pdf.second->GetEntries() < 50)   	pdf.second->Rebin(5);							// Make the pdf more stable
      else if(pdf.second->GetEntries() < 500)  	pdf.second->Rebin(2);
      pdf.second->Scale(1./pdf.second->Integral(0, pdf.second->GetNbinsX() + 1, "width"));				// Scale to integral (also include underflow/overflow)
      pdf.second->Write();
    }

    for(auto& pdf : pdfs) delete pdf.second;
    for(auto& file : {pdfFile, qgMiniTupleFile}){ file->Close(); delete file;}
  }
  return 0;
}

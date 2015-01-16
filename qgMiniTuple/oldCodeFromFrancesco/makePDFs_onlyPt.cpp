#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TLine.h"
#include "TVectorT.h"

#include "../macros/binFunctions.h"



bool debug = false;


void normalizeHistos( std::vector<TH1F*> vh1 );
//int findBin( float x, std::vector<float> v);


int main() {


  TFile* file = TFile::Open("/afs/cern.ch/work/t/tomc/public/qgMiniTuples/qgMiniTuple_QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14.root");
  std::cout << "-> Opened file: " << file->GetName() << std::endl;
  TTree* tree = (TTree*)file->Get("qgMiniTupleAK4chs/qgMiniTuple");


  float rho;
  tree->SetBranchAddress("rho", &rho);
  float pt;
  tree->SetBranchAddress("pt", &pt);
  float eta;
  tree->SetBranchAddress("eta", &eta);
  float ptD;
  tree->SetBranchAddress("ptD", &ptD);
  float axis2;
  tree->SetBranchAddress("axis2", &axis2);
  int mult;
  tree->SetBranchAddress("mult", &mult);
  int partonId;
  tree->SetBranchAddress("partonId", &partonId);
  float bTag;
  tree->SetBranchAddress("bTag", &bTag);
  int jetIdLevel;
  tree->SetBranchAddress("jetIdLevel", &jetIdLevel);
  int nGenJetsInCone;
  tree->SetBranchAddress("nGenJetsInCone", &nGenJetsInCone);
  bool balanced;
  tree->SetBranchAddress("balanced", &balanced);


  std::vector<float> etaBins = {0,2.5};
  std::vector<float> ptBinsC; getBins(ptBinsC, 20, 20, 2000, true); ptBinsC.push_back(4000);

  if( debug ) {

    std::cout << "ptBinsC: " << std::endl;
    for( unsigned i=0; i<ptBinsC.size(); ++i ) std::cout << ptBinsC[i] << std::endl;
    std::cout << std::endl;

  }


  int nbins = 100;


  std::vector< TH1F* > vh1_mult_quark;
  std::vector< TH1F* > vh1_mult_gluon;

  std::vector< TH1F* > vh1_ptD_quark;
  std::vector< TH1F* > vh1_ptD_gluon;

  std::vector< TH1F* > vh1_axis2_quark;
  std::vector< TH1F* > vh1_axis2_gluon;



  for( unsigned i=0; i<ptBinsC.size()-1; ++i ) { 

    TH1F* h1_mult_quark = new TH1F(Form("mult_quark_eta0_pt%d", i), "", 100, 0, 100);
    TH1F* h1_mult_gluon = new TH1F(Form("mult_gluon_eta0_pt%d", i), "", 100, 0, 100);

    vh1_mult_quark.push_back(h1_mult_quark);
    vh1_mult_gluon.push_back(h1_mult_gluon);
    
    TH1F* h1_ptD_quark = new TH1F(Form("ptD_quark_eta0_pt%d", i), "", nbins, 0., 1.0001);
    TH1F* h1_ptD_gluon = new TH1F(Form("ptD_gluon_eta0_pt%d", i), "", nbins, 0., 1.0001);

    vh1_ptD_quark.push_back(h1_ptD_quark);
    vh1_ptD_gluon.push_back(h1_ptD_gluon);
    
    TH1F* h1_axis2_quark = new TH1F(Form("axis2_quark_eta0_pt%d", i), "", nbins, 0., 10.);
    TH1F* h1_axis2_gluon = new TH1F(Form("axis2_gluon_eta0_pt%d", i), "", nbins, 0., 10.);

    vh1_axis2_quark.push_back(h1_axis2_quark);
    vh1_axis2_gluon.push_back(h1_axis2_gluon);
      
      
  } // for pt bins






  int nentries = tree->GetEntries();
  if( debug ) 
    nentries=100;



  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    if( iEntry % 500000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;

    if( pt<20. ) continue;
    if( fabs(eta)>2.5 ) continue;
    if( !balanced ) continue;
    if( nGenJetsInCone!=1 ) continue;
    if( bTag > 0.244 ) continue; // CSVL
    
    bool is_uds = (partonId!=0 && abs(partonId)<4);
    bool is_g   = partonId==21;
    if( !is_uds && !is_g ) continue;	// Keep only udsg

    std::vector<float> vars;
    vars.push_back( mult );
    vars.push_back( ptD );
    vars.push_back( axis2 );
    
    int ibin = -1;
    if( !getBinNumber( ptBinsC, pt, ibin) ) continue;

    if( debug ) {
      std::cout << "pt: " << pt << " ptBin: " << ibin << std::endl;
    }

    if( ibin<0 ) {
      std::cout << "CONTINUING (pt: " << pt << " )" << std::endl;
      continue;
    }

    if( is_uds ) {

      vh1_mult_quark.at(ibin)->Fill( mult );
      vh1_ptD_quark.at(ibin)->Fill( ptD );
      vh1_axis2_quark.at(ibin)->Fill( axis2 );

    } else if( is_g ) {

      vh1_mult_gluon.at(ibin)->Fill( mult );
      vh1_ptD_gluon.at(ibin)->Fill( ptD );
      vh1_axis2_gluon.at(ibin)->Fill( axis2 );

    }

  
  } // for entries


  //normalizeHistos( vh1_mult_quark );
  //normalizeHistos( vh1_mult_gluon );
  //normalizeHistos( vh1_ptD_quark );
  //normalizeHistos( vh1_ptD_gluon );
  //normalizeHistos( vh1_axis2_quark );
  //normalizeHistos( vh1_axis2_gluon );
  

  std::string outfilename = (debug) ? "pdfsOnlyPt_tmp.root" : "pdfsOnlyPt.root";

  TFile* outfile = TFile::Open(outfilename.c_str(), "recreate");
  outfile->cd();
  outfile->mkdir("mult");
  outfile->cd("mult");
  for( unsigned i=0; i<vh1_mult_quark.size(); ++i ) {
    vh1_mult_quark[i]->Write();
    vh1_mult_gluon[i]->Write();
  }
  outfile->cd();
  outfile->mkdir("ptD");
  outfile->cd("ptD");
  for( unsigned i=0; i<vh1_ptD_quark.size(); ++i ) {
    vh1_ptD_quark[i]->Write();
    vh1_ptD_gluon[i]->Write();
  }
  outfile->cd();
  outfile->mkdir("axis2");
  outfile->cd("axis2");
  for( unsigned i=0; i<vh1_axis2_quark.size(); ++i ) {
    vh1_axis2_quark[i]->Write();
    vh1_axis2_gluon[i]->Write();
  }

  outfile->cd();

  TVectorT<float> etaBins_t(etaBins.size(), &etaBins[0]);
  etaBins_t.Write("etaBins");
  TVectorT<float> ptBinsC_t(ptBinsC.size(), &ptBinsC[0]);
  ptBinsC_t.Write("ptBinsC");


  outfile->Close();
    
  return 0;

}



void normalizeHistos( std::vector<TH1F*> vh1 ) {

  for( unsigned i=0; i<vh1.size(); ++i ) 
    for( unsigned ibin=1; ibin<vh1[i]->GetNbinsX()+1; ++ibin ) 
      vh1[i]->SetBinContent(ibin, vh1[i]->GetBinContent(ibin)/vh1[i]->Integral());


}



//int findBin( float x, std::vector<float> v) {
//
//  int foundBin=-1;
//
//  if( x<v[0] ) foundBin=0;
//  else if( x>v[v.size()-1] ) foundBin=v.size()-2;
//  else {
//    for( unsigned i=0; i<v.size()-1; ++i ) {
//      if( x>=v[i] && x<v[i+1]) {
//        foundBin = i;
//        break;
//      }
//    }
//  }
//  
//
//  return foundBin;
//
//}




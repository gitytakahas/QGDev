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

#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator2.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.h"
#include "../macros/binFunctions.h"
#include "DrawTools.h"



void drawSingleRoC( const std::string& outputdir, const std::string& saveName, float ptMin, float ptMax, const std::string& name1, TH1F* h1_qgl_quark, TH1F* h1_qgl_gluon, const std::string& name2, TH1F* h1_qgl2_quark, TH1F* h1_qgl2_gluon );
TGraph* getSingleRoC( const std::string& name, TH1F* h1_quark, TH1F* h1_gluon );
void drawQuarkVsGluon( const std::string& outputdir, const std::string& savename, const std::string& axisName, float ptMin, float ptMax, TH1F* h1_quark, TH1F* h1_gluon );



int main() {

  DrawTools::setStyle();

  std::string qgfileName = "pdfsOnlyPt.root";
  QGLikelihoodCalculator2* qglc = new QGLikelihoodCalculator2(qgfileName); 

  std::string qgfileName_old = "/afs/cern.ch/user/t/tomc/public/QG_pdfs_13TeV_2014-10-12/pdfQG_AK4chs_antib_NoQC_13TeV.root";
  QGLikelihoodCalculator* qglc_old = new QGLikelihoodCalculator(qgfileName_old); 


  TFile* file = TFile::Open("/afs/cern.ch/work/t/tomc/public/qgMiniTuples/qgMiniTuple_QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14.root");
  std::cout << "-> Opened file: " << file->GetName() << std::endl;
  TTree* tree = (TTree*)file->Get("qgMiniTupleAK4chs/qgMiniTuple");

  //TFile* file = TFile::Open("/afs/cern.ch/work/t/tomc/public/qgMiniTuples/qgMiniTupleForMiniAOD_TTJets.root");
  //std::cout << "-> Opened file: " << file->GetName() << std::endl;
  //TTree* tree = (TTree*)file->Get("qgMiniTupleMiniAOD/qgMiniTupleForMiniAOD");


  int event;
  tree->SetBranchAddress("nEvent", &event);
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
  //int jetIdLevel;
  //tree->SetBranchAddress("jetIdLevel", &jetIdLevel);
  bool hasGenJet;
  tree->SetBranchAddress("hasGenJet", &hasGenJet);
  bool balanced;
  tree->SetBranchAddress("balanced", &balanced);
  int nGenJetsInCone;
  tree->SetBranchAddress("nGenJetsInCone", &nGenJetsInCone);

  std::vector<float> etaBins = {0,2.5};
  std::vector<float> ptBinsC; getBins(ptBinsC, 100, 20, 2000, true); ptBinsC.push_back(4000);


  int nbins = 50;

  std::vector< TH1F* > vh1_mult_quark;
  std::vector< TH1F* > vh1_mult_gluon;

  std::vector< TH1F* > vh1_ptD_quark;
  std::vector< TH1F* > vh1_ptD_gluon;

  std::vector< TH1F* > vh1_axis2_quark;
  std::vector< TH1F* > vh1_axis2_gluon;

  std::vector< TH1F* > vh1_qgl_quark;
  std::vector< TH1F* > vh1_qgl_gluon;

  std::vector< TH1F* > vh1_qgl_old_quark;
  std::vector< TH1F* > vh1_qgl_old_gluon;


  for( unsigned i=0; i<ptBinsC.size()-1; ++i ) { // integrated in rho

    TH1F* h1_mult_quark = new TH1F(Form("mult_quark_pt%d", i), "", 100, 0., 100);
    TH1F* h1_mult_gluon = new TH1F(Form("mult_gluon_pt%d", i), "", 100, 0., 100);

    vh1_mult_quark.push_back(h1_mult_quark);
    vh1_mult_gluon.push_back(h1_mult_gluon);


    TH1F* h1_ptD_quark = new TH1F(Form("ptD_quark_pt%d", i), "", nbins, 0., 1.0001);
    TH1F* h1_ptD_gluon = new TH1F(Form("ptD_gluon_pt%d", i), "", nbins, 0., 1.0001);

    vh1_ptD_quark.push_back(h1_ptD_quark);
    vh1_ptD_gluon.push_back(h1_ptD_gluon);


    TH1F* h1_axis2_quark = new TH1F(Form("axis2_quark_pt%d", i), "", nbins, 0., 10);
    TH1F* h1_axis2_gluon = new TH1F(Form("axis2_gluon_pt%d", i), "", nbins, 0., 10);

    vh1_axis2_quark.push_back(h1_axis2_quark);
    vh1_axis2_gluon.push_back(h1_axis2_gluon);


    
    TH1F* h1_qgl_quark = new TH1F(Form("qgl_quark_pt%d", i), "", nbins, 0., 1.0001);
    TH1F* h1_qgl_gluon = new TH1F(Form("qgl_gluon_pt%d", i), "", nbins, 0., 1.0001);

    vh1_qgl_quark.push_back(h1_qgl_quark);
    vh1_qgl_gluon.push_back(h1_qgl_gluon);
    

    TH1F* h1_qgl_old_quark = new TH1F(Form("qgl_old_quark_pt%d", i), "", nbins, 0., 1.0001);
    TH1F* h1_qgl_old_gluon = new TH1F(Form("qgl_old_gluon_pt%d", i), "", nbins, 0., 1.0001);

    vh1_qgl_old_quark.push_back(h1_qgl_old_quark);
    vh1_qgl_old_gluon.push_back(h1_qgl_old_gluon);
    
  }








  std::string outputdir = "RoCs_onlyPt";
  system( Form("mkdir -p %s", outputdir.c_str()));

  int nentries = tree->GetEntries();
  int nentries_max = 10000000;
  if( nentries_max>0 && nentries>nentries_max ) nentries = nentries_max;

  int njets_tot = 0;
  int njets_uds = 0;
  int njets_g = 0;

  int cachedEvent = -1;

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    if( iEntry % 500000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;

    if( pt<20. ) continue;
    if( fabs(eta)>2.5 ) continue;
    //if( !hasGenJet ) continue;
    if( nGenJetsInCone!=1 ) continue;
    if( !balanced ) continue;
    if( bTag > 0.244 ) continue; // CSVL
    
    //if( event!=cachedEvent ) {
    //  cachedEvent = event;
    //} else { 
    //  continue;
    //}
      
    bool is_uds = (partonId!=0 && abs(partonId)<4);
    bool is_g   = partonId==21;
    njets_tot += 1;
    if( is_uds ) njets_uds += 1;
    if( is_g ) njets_g += 1;
    if( !is_uds && !is_g ) continue;	// Keep only udsg

    std::vector<float> vars;
    vars.push_back( mult );
    vars.push_back( ptD );
    vars.push_back( axis2 );

    float qgl = qglc->computeQGLikelihood     ( pt, eta, vars );
    float qgl_old = qglc_old->computeQGLikelihood     ( pt, eta, rho, vars );
    
    int ibin = -1;
    if( !getBinNumber( ptBinsC, pt, ibin) ) continue;


    if( ibin<0 ) {
      std::cout << "CONTINUING (pt: " << pt << " )" << std::endl;
      continue;
    }


    if( is_uds ) {

      vh1_mult_quark[ibin]->Fill( mult );
      vh1_ptD_quark[ibin]->Fill( ptD );
      vh1_axis2_quark[ibin]->Fill( axis2 );

      vh1_qgl_quark[ibin]->Fill( qgl );
      vh1_qgl_old_quark[ibin]->Fill( qgl_old );

    } else if( is_g ) {

      vh1_mult_gluon[ibin]->Fill( mult );
      vh1_ptD_gluon[ibin]->Fill( ptD );
      vh1_axis2_gluon[ibin]->Fill( axis2 );

      vh1_qgl_gluon[ibin]->Fill( qgl );
      vh1_qgl_old_gluon[ibin]->Fill( qgl_old );

    }

  
  } // for entries


  std::cout << "quark fraction: " << (float)njets_uds/njets_tot << std::endl;
  std::cout << "gluon fraction: " << (float)njets_g/njets_tot << std::endl;

  for( unsigned i=0; i<vh1_qgl_quark.size(); ++i ) {

    drawSingleRoC( outputdir, "vsRhoBins", ptBinsC[i], ptBinsC[i+1], "QG Likelihood", vh1_qgl_quark[i], vh1_qgl_gluon[i], "QG Likelihood (#rho bins)", vh1_qgl_old_quark[i], vh1_qgl_old_gluon[i] );
    //drawSingleRoC( outputdir, "vsNoQC", ptBinsC[i], ptBinsC[i+1], "QG Likelihood", vh1_qgl_quark[i], vh1_qgl_gluon[i], "QG Likelihood (no QC)", vh1_qgl_noQC_quark[i], vh1_qgl_noQC_gluon[i] );

    drawQuarkVsGluon( outputdir, "mult", "Jet Multiplicity", ptBinsC[i], ptBinsC[i+1], vh1_mult_quark[i], vh1_mult_gluon[i] );
    drawQuarkVsGluon( outputdir, "ptD", "p_{T}D", ptBinsC[i], ptBinsC[i+1], vh1_ptD_quark[i], vh1_ptD_gluon[i] );
    drawQuarkVsGluon( outputdir, "axis2", "-log( #sigma_{2} )", ptBinsC[i], ptBinsC[i+1], vh1_axis2_quark[i], vh1_axis2_gluon[i] );
    drawQuarkVsGluon( outputdir, "qgl", "Quark-Gluon Likelihood", ptBinsC[i], ptBinsC[i+1], vh1_qgl_quark[i], vh1_qgl_gluon[i] );
    drawQuarkVsGluon( outputdir, "qgl_old", "Quark-Gluon Likelihood (#rho bins)", ptBinsC[i], ptBinsC[i+1], vh1_qgl_old_quark[i], vh1_qgl_old_gluon[i] );
  }


  std::string outfilename = "pdfs_plus_qgl_onlyPt.root";

  TFile* outfile = TFile::Open(outfilename.c_str(), "recreate");
  outfile->cd();
  outfile->mkdir("qgl");
  outfile->cd("qgl");
  for( unsigned i=0; i<vh1_qgl_quark.size(); ++i ) {
    vh1_qgl_quark.at(i)->Write();
    vh1_qgl_gluon.at(i)->Write();
  }
  //outfile->cd();
  //outfile->mkdir("mult");
  //outfile->cd("mult");
  //for( unsigned i=0; i<vvh1_mult_quark.size(); ++i ) {
  //  for( unsigned j=0; j<vvh1_mult_quark[i].size(); ++j ) {
  //    vvh1_mult_quark.at(i).at(j)->Write();
  //    vvh1_mult_gluon.at(i).at(j)->Write();
  //  }
  //}
  //outfile->cd();
  //outfile->mkdir("ptD");
  //outfile->cd("ptD");
  //for( unsigned i=0; i<vvh1_ptD_quark.size(); ++i ) {
  //  for( unsigned j=0; j<vvh1_ptD_quark[i].size(); ++j ) {
  //    vvh1_ptD_quark.at(i).at(j)->Write();
  //    vvh1_ptD_gluon.at(i).at(j)->Write();
  //  }
  //}
  //outfile->cd();
  //outfile->mkdir("axis2");
  //outfile->cd("axis2");
  //for( unsigned i=0; i<vvh1_axis2_quark.size(); ++i ) {
  //  for( unsigned j=0; j<vvh1_axis2_quark[i].size(); ++j ) {
  //    vvh1_axis2_quark.at(i).at(j)->Write();
  //    vvh1_axis2_gluon.at(i).at(j)->Write();
  //  }
  //}

  //outfile->cd();

  TVectorT<float> etaBins_t(etaBins.size(), &etaBins[0]);
  etaBins_t.Write("etaBins");
  TVectorT<float> ptBinsC_t(ptBinsC.size(), &ptBinsC[0]);
  ptBinsC_t.Write("ptBinsC");

  outfile->Close();
    
  return 0;

}






void drawSingleRoC( const std::string& outputdir, const std::string& saveName, float ptMin, float ptMax, const std::string& name1, TH1F* h1_qgl_quark, TH1F* h1_qgl_gluon, const std::string& name2, TH1F* h1_qgl2_quark, TH1F* h1_qgl2_gluon ) {

  TGraph* gr_qgl = getSingleRoC( "qgl", h1_qgl_quark, h1_qgl_gluon );
  TGraph* gr_qgl2 = getSingleRoC( "qgl2", h1_qgl2_quark, h1_qgl2_gluon );

  TCanvas* c1 = new TCanvas("c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, 0., 1.0001, 10, 0., 1.0001);
  h2_axes->SetYTitle( "Quark Efficiency" );
  h2_axes->SetXTitle( "Gluon Rejection" );
  h2_axes->Draw();

  TLine* diag = new TLine( 0., 1., 1., 0.);
  diag->Draw("same");

  gr_qgl->SetMarkerSize(1.5);
  gr_qgl2->SetMarkerSize(1.5);

  gr_qgl->SetMarkerColor(kOrange);
  gr_qgl2->SetMarkerColor(kRed+2);

  gr_qgl->SetMarkerStyle(20);
  gr_qgl2->SetMarkerStyle(24);


  gr_qgl->Draw("P same");
  gr_qgl2->Draw("P same");

 
  
  TLegend* legend = new TLegend( 0.2, 0.2, 0.55, 0.45, Form("%.0f < p_{T} < %.0f GeV", ptMin, ptMax) );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( gr_qgl, name1.c_str(), "P" );
  legend->AddEntry( gr_qgl2, name2.c_str(), "P" );
  legend->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/RoC_%s_pt%.0f.eps", outputdir.c_str(), saveName.c_str(), ptMin) );
  c1->SaveAs( Form("%s/RoC_%s_pt%.0f.png", outputdir.c_str(), saveName.c_str(), ptMin) );

  delete c1;
  delete h2_axes;

}



TGraph* getSingleRoC( const std::string& name, TH1F* h1_quark, TH1F* h1_gluon ) {

  TGraph* gr_RoC = new TGraph(0);
  gr_RoC->SetName(Form("roc_%s", name.c_str()));

  int nbins = h1_quark->GetNbinsX();

  for( unsigned ibin=1; ibin<nbins+1; ++ibin ) {

    float eff_q = h1_quark->Integral( nbins-ibin, nbins )/h1_quark->Integral( 1, nbins );
    float eff_g = h1_gluon->Integral( nbins-ibin, nbins )/h1_gluon->Integral( 1, nbins );
    
    gr_RoC->SetPoint( ibin-1, 1.-eff_g, eff_q );

  }

  return gr_RoC;

}




void drawQuarkVsGluon( const std::string& outputdir, const std::string& savename, const std::string& axisName, float ptMin, float ptMax, TH1F* h1_quark, TH1F* h1_gluon ) {

  TCanvas* c2 = new TCanvas( "c2", "", 600, 600 );
  c2->cd();

  float xMax = h1_quark->GetXaxis()->GetXmax();

  float ymax_q = h1_quark->GetMaximum()/h1_quark->Integral();
  float ymax_g = h1_gluon->GetMaximum()/h1_gluon->Integral();

  float yMax = (ymax_q>ymax_g) ? ymax_q : ymax_g;
  yMax *= 1.1;

  TH2D* h2_axes = new TH2D("axes", "", 10, 0., xMax, 10, 0., yMax );
  h2_axes->SetXTitle( axisName.c_str() );
  h2_axes->SetYTitle( "Normalized To Unity" );
  h2_axes->Draw();

  h1_quark->SetFillStyle(3004);
  h1_gluon->SetFillStyle(3005);

  h1_quark->SetFillColor(38);
  h1_gluon->SetFillColor(46);

  h1_quark->SetLineColor(38);
  h1_gluon->SetLineColor(46);

  h1_quark->SetLineWidth(2);
  h1_gluon->SetLineWidth(2);

  float xMin_legend = (savename=="qgl" || savename=="qglcdf" || savename=="qgl_noQC") ? 0.3 : 0.6;
  float xMax_legend = xMin_legend + 0.35;

  TLegend* legend = new TLegend( xMin_legend, 0.7, xMax_legend, 0.9, Form("%.0f < p_{T} < %.0f GeV", ptMin, ptMax) );
  legend->SetTextSize(0.038);
  legend->SetFillColor(0);
  legend->AddEntry( h1_quark, "Quarks", "F" );
  legend->AddEntry( h1_gluon, "Gluons", "F" );
  legend->Draw("same");

  h1_quark->DrawNormalized("same");
  h1_gluon->DrawNormalized("same");


  gPad->RedrawAxis();

  c2->SaveAs( Form("%s/%s_pt%.0f.eps", outputdir.c_str(), savename.c_str(), ptMin) );
  c2->SaveAs( Form("%s/%s_pt%.0f.png", outputdir.c_str(), savename.c_str(), ptMin) );

  delete h2_axes;
  delete c2;


}




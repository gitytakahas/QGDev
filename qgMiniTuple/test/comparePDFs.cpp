#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"


void compareSinglePDF( const std::string& outputdir, TFile* file_tom, TFile* file_fp, const std::string& var, int etaBin, int ptBin, int rhoBin );
void compareVsRho( const std::string& outputdir, TFile* file_tom, const std::string& var, int etaBin, int ptBin );
void compareVsPt( const std::string& outputdir, TFile* file, const std::string& var );


int main() {

  

  TFile* file_tom = TFile::Open("/afs/cern.ch/user/t/tomc/public/QG_pdfs_13TeV_2014-10-12/pdfQG_AK4chs_antib_NoQC_13TeV.root");
  TFile* file_fp  = TFile::Open("pdfs.root");
  
  std::string outputdir = "comparePDFsPlots";
  system( Form("mkdir -p %s", outputdir.c_str()) );
  compareSinglePDF( outputdir, file_tom, file_fp, "axis2", 0, 10, 10 );
  compareSinglePDF( outputdir, file_tom, file_fp, "ptD", 0, 10, 10 );
  compareSinglePDF( outputdir, file_tom, file_fp, "mult", 0, 10, 10 );

  compareVsRho( outputdir, file_tom, "mult", 0, 10 );
  compareVsRho( outputdir, file_tom, "ptD", 0, 10 );
  compareVsRho( outputdir, file_tom, "axis2", 0, 10 );

  TFile* file_noRho  = TFile::Open("pdfsOnlyPt.root");
  compareVsPt( outputdir, file_noRho, "mult" );
  compareVsPt( outputdir, file_noRho, "ptD" );
  compareVsPt( outputdir, file_noRho, "axis2" );


  return 0;

}



void compareSinglePDF( const std::string& outputdir, TFile* file_tom, TFile* file_fp, const std::string& var, int etaBin, int ptBin, int rhoBin ) {

  std::string name_quark( Form("%s/%s_quark_eta-%d_pt-%d_rho-%d", var.c_str(), var.c_str(), etaBin, ptBin, rhoBin) );
  std::string name_gluon( Form("%s/%s_gluon_eta-%d_pt-%d_rho-%d", var.c_str(), var.c_str(), etaBin, ptBin, rhoBin) );

  std::string name_quark_p1( Form("%s/%s_quark_eta-%d_pt-%d_rho-%d", var.c_str(), var.c_str(), etaBin, ptBin, rhoBin-4) );
  std::string name_gluon_p1( Form("%s/%s_gluon_eta-%d_pt-%d_rho-%d", var.c_str(), var.c_str(), etaBin, ptBin, rhoBin-4) );

  TH1F* h1_gluon_fp = (TH1F*)file_fp->Get(name_gluon.c_str());
  TH1F* h1_quark_fp = (TH1F*)file_fp->Get(name_quark.c_str());

  TH1F* h1_gluon_tom = (TH1F*)file_tom->Get(name_gluon_p1.c_str());
  TH1F* h1_quark_tom = (TH1F*)file_tom->Get(name_quark_p1.c_str());

  h1_gluon_fp->SetMarkerSize(1.5);
  h1_gluon_fp->SetMarkerColor(kRed);
  h1_gluon_fp->SetMarkerStyle(20);

  h1_quark_fp->SetMarkerSize(1.5);
  h1_quark_fp->SetMarkerColor(kBlue);
  h1_quark_fp->SetMarkerStyle(20);

  h1_gluon_tom->SetLineWidth(2);
  h1_gluon_tom->SetLineColor(kRed);

  h1_quark_tom->SetLineWidth(2);
  h1_quark_tom->SetLineColor(kBlue);

  if( var!="mult" ) {
    h1_gluon_tom->Rebin(4);
    h1_quark_tom->Rebin(4);
  }


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  std::cout << h1_gluon_fp->GetNbinsX() << std::endl;;
  std::cout << h1_quark_fp->GetNbinsX() << std::endl;;
  std::cout << h1_gluon_tom->GetNbinsX() << std::endl;;
  std::cout << h1_quark_tom->GetNbinsX() << std::endl;;

  h1_gluon_fp->DrawNormalized("P");
  h1_quark_fp->DrawNormalized("P same");
  h1_gluon_tom->DrawNormalized("same");
  h1_quark_tom->DrawNormalized("same");

  c1->SaveAs(Form("%s/%s_eta-%d_pt-%d_rho-%d.eps", outputdir.c_str(), var.c_str(), etaBin, ptBin, rhoBin));
  c1->SaveAs(Form("%s/%s_eta-%d_pt-%d_rho-%d.png", outputdir.c_str(), var.c_str(), etaBin, ptBin, rhoBin));

  delete c1;

}



void compareVsRho( const std::string& outputdir, TFile* file, const std::string& var, int etaBin, int ptBin ) {


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  for( unsigned i=0; i<4; ++i ) {

    int rhoBin = 10*i;
    std::string name_quark( Form("%s/%s_quark_eta-%d_pt-%d_rho-%d", var.c_str(), var.c_str(), etaBin, ptBin, rhoBin) );
    std::string name_gluon( Form("%s/%s_gluon_eta-%d_pt-%d_rho-%d", var.c_str(), var.c_str(), etaBin, ptBin, rhoBin) );

    TH1F* h1_quark = (TH1F*)file->Get(name_quark.c_str());
    TH1F* h1_gluon = (TH1F*)file->Get(name_gluon.c_str());
  
    h1_quark->SetLineColor(kBlue+i);
    h1_gluon->SetLineColor(kRed+i);

    h1_quark->SetLineWidth(2);
    h1_gluon->SetLineWidth(2);

    if( i==0 )
      h1_quark->DrawNormalized("");
    else
      h1_quark->DrawNormalized("same");

    h1_gluon->DrawNormalized("same");

  } 

  gPad->RedrawAxis();

  c1->SaveAs(Form("%s/%s_eta-%d_pt-%d_vsRho.eps", outputdir.c_str(), var.c_str(), etaBin, ptBin));
  c1->SaveAs(Form("%s/%s_eta-%d_pt-%d_vsRho.png", outputdir.c_str(), var.c_str(), etaBin, ptBin));

  delete c1;


}


void compareVsPt( const std::string& outputdir, TFile* file, const std::string& var ) {

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  for( unsigned i=1; i<5; ++i ) {

    int ptBin = 3*i + 4;
    std::string name_quark( Form("%s/%s_quark_eta0_pt%d", var.c_str(), var.c_str(), ptBin) );
    std::string name_gluon( Form("%s/%s_gluon_eta0_pt%d", var.c_str(), var.c_str(), ptBin) );

    TH1F* h1_quark = (TH1F*)file->Get(name_quark.c_str());
    TH1F* h1_gluon = (TH1F*)file->Get(name_gluon.c_str());
  
    h1_quark->SetLineColor(kBlue+i);
    h1_gluon->SetLineColor(kRed+i);

    h1_quark->SetLineWidth(2);
    h1_gluon->SetLineWidth(2);

    if( i==0 )
      h1_quark->DrawNormalized("");
    else
      h1_quark->DrawNormalized("same");

    h1_gluon->DrawNormalized("same");

  } 

  gPad->RedrawAxis();

  c1->SaveAs(Form("%s/%s_vsPt.eps", outputdir.c_str(), var.c_str()));
  c1->SaveAs(Form("%s/%s_vsPt.png", outputdir.c_str(), var.c_str()));

  delete c1;

  


}

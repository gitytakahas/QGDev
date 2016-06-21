#include <map>
#include "TH1D.h"
#include "TCanvas.h"
#include "binClass.h"
#include "binningConfigurations.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.h"
#include "../localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc"

/*
 * Plot the single variable likelihood for each bin of a pdf
 * (useful to inspect if pdf's are smooth enough, even though fluctuations on the edges or low statistic bins are unavoidable)
 */
int main(int argc, char**argv){
  system("rm -r ./plots/likelihoodValues/");

  QGLikelihoodCalculator localQG("pdfQG_AK4chs_13TeV_76X.root");						// Init localQGLikelihoodCalculator
  binClass bins = get76XBinning();									// Define binning for the plots (i.e. the same as used to create the pdf's)

  std::map<TString, TH1D*> plots;
  for(TString binName : bins.getAllBinNames()){
    plots["axis2"]	= new TH1D("axis2_" + binName, "axis2_" + binName, 200, 1, 9);
    plots["ptD"]	= new TH1D("ptD_"   + binName, "ptD_"   + binName, 200, 0, 1);
    plots["mult"]	= new TH1D("mult_"  + binName, "mult_"  + binName, 140, 2.5, 144.5);

    for(int bin = 0; bin <= plots["mult"]->GetNbinsX() + 1; ++bin){
      float var     		= plots["mult"]->GetBinCenter(bin);
      float var_l   		= localQG.computeQGLikelihood(binName, {var, -1, -1});
      plots["mult"]->SetBinContent(bin, var_l);
    }
    for(int bin = 0; bin < plots["ptD"]->GetNbinsX() + 1; ++bin){
      float var     		= plots["ptD"]->GetBinCenter(bin);
      float var_l   		= localQG.computeQGLikelihood(binName, {-1, var, -1});
      plots["ptD"]->SetBinContent(bin, var_l);
    }
    for(int bin = 0; bin < plots["axis2"]->GetNbinsX() + 1; ++bin){
      float var     		= plots["axis2"]->GetBinCenter(bin);
      float var_l   		= localQG.computeQGLikelihood(binName, {-1, -1, var});
      plots["axis2"]->SetBinContent(bin, var_l);
    }

    for(auto& plot : plots){
      TCanvas c;
      plot.second->SetLineWidth(2);
      plot.second->SetStats(0);
      plot.second->GetXaxis()->SetTitle(plot.first);
      plot.second->GetYaxis()->SetTitle("likelihood");
      plot.second->Draw();
      TString pdfDir = "./plots/likelihoodValues/" + plot.first + "/";
      TString pdfName = binName + ".pdf";
      system("mkdir -p " + pdfDir);
      c.SaveAs(pdfDir + pdfName);
    }
  }
  return 0;
}

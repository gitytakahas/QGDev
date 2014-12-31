#include "binClass.h"

binClass getDefaultBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.3,1.5,2,2.5,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));				// i.e. 20 bins from 20 to 2000 with log=true and with an additional bin up to 6500
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.printBinRanges();

  for(int i=10; i < bins.getNBins("pt"); ++i) bins.setLink("eta5_pt9_rho0",  TString("eta5_pt") + i + "_rho0");
  for(int i=14; i < bins.getNBins("pt"); ++i) bins.setLink("eta4_pt13_rho0", TString("eta4_pt") + i + "_rho0");
  for(int i=15; i < bins.getNBins("pt"); ++i) bins.setLink("eta3_pt14_rho0", TString("eta3_pt") + i + "_rho0");
  for(int i=17; i < bins.getNBins("pt"); ++i) bins.setLink("eta2_pt16_rho0", TString("eta2_pt") + i + "_rho0");
  for(int i=20; i < bins.getNBins("pt"); ++i) bins.setLink("eta1_pt19_rho0", TString("eta1_pt") + i + "_rho0");

  return bins;
}



binClass get8TeVBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,2.5,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.printBinRanges();

  for(int i=10; i < bins.getNBins("pt"); ++i) bins.setLink("eta1_pt9_rho0", TString("eta1_pt") + i + "_rho0");

  return bins;
}



binClass getSmallEtaBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.3,1.5,1.8,2.1,2.4,2.7,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.printBinRanges();

  for(int i=10; i < bins.getNBins("pt"); ++i) bins.setLink("eta7_pt9_rho0",  TString("eta7_pt") + i + "_rho0");
  for(int i=14; i < bins.getNBins("pt"); ++i) bins.setLink("eta6_pt13_rho0", TString("eta6_pt") + i + "_rho0");
  for(int i=15; i < bins.getNBins("pt"); ++i) bins.setLink("eta5_pt14_rho0", TString("eta5_pt") + i + "_rho0");
  for(int i=17; i < bins.getNBins("pt"); ++i) bins.setLink("eta4_pt16_rho0", TString("eta4_pt") + i + "_rho0");
  for(int i=18; i < bins.getNBins("pt"); ++i) bins.setLink("eta3_pt17_rho0", TString("eta3_pt") + i + "_rho0");
  for(int i=19; i < bins.getNBins("pt"); ++i) bins.setLink("eta2_pt18_rho0", TString("eta2_pt") + i + "_rho0");
  for(int i=20; i < bins.getNBins("pt"); ++i) bins.setLink("eta1_pt19_rho0", TString("eta1_pt") + i + "_rho0");

  return bins;
}



binClass getV1Binning(TString weightOption){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.3,1.5,1.8,2.1,2.4,2.7,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));
  bins.setBinRange("rho",	"#rho",			{0,15,20,25,30,9999});
  bins.printBinRanges();

  for(int j=0; j < 5; ++j){														// No rho binning needed in the central (eta < 2.4)
    for(int i=0; i < bins.getNBins("pt"); ++i){
      for(int k=1; k < bins.getNBins("rho"); ++k){
        bins.setLink(TString("eta") + j + "_pt" + i + "_rho0", TString("eta") + j + "_pt" + i + "_rho" + k);
      }
    }
  }
  for(int j=5; j < bins.getNBins("eta"); ++j){												// Unfortunately, need to merge some rho bins in lowest pt bins of forward because of very low statistics
    for(int i=0; i < 3; ++i){														// If we get a forward jets QCD sample, try to go with finer binning, it could help a lot to boost performance
      bins.setLink(TString("eta") + j + "_pt" + i + "_rho2", TString("eta") + j + "_pt" + i + "_rho3");
      bins.setLink(TString("eta") + j + "_pt" + i + "_rho2", TString("eta") + j + "_pt" + i + "_rho4");
    }
  }

  for(int j=0; j < bins.getNBins("rho"); ++j){
    for(int i=9;  i < bins.getNBins("pt"); ++i) bins.setLink(TString("eta7_pt8_rho")  + j, TString("eta7_pt") + i + "_rho" + j);	// Merge high pt bins in forward (no or low statistics)
    for(int i=13; i < bins.getNBins("pt"); ++i) bins.setLink(TString("eta6_pt12_rho") + j, TString("eta6_pt") + i + "_rho" + j);
    for(int i=14; i < bins.getNBins("pt"); ++i) bins.setLink(TString("eta5_pt13_rho") + j, TString("eta5_pt") + i + "_rho" + j);
    for(int i=16; i < bins.getNBins("pt"); ++i) bins.setLink(TString("eta4_pt15_rho") + j, TString("eta4_pt") + i + "_rho" + j);
    for(int i=17; i < bins.getNBins("pt"); ++i) bins.setLink(TString("eta3_pt16_rho") + j, TString("eta3_pt") + i + "_rho" + j);
    for(int i=18; i < bins.getNBins("pt"); ++i) bins.setLink(TString("eta2_pt17_rho") + j, TString("eta2_pt") + i + "_rho" + j);
    for(int i=19; i < bins.getNBins("pt"); ++i) bins.setLink(TString("eta1_pt18_rho") + j, TString("eta1_pt") + i + "_rho" + j);
   
    for(int i = 3; i < bins.getNBins("eta"); ++i){
      bins.setLink(TString("eta") + i + "_pt0_rho" + j, TString("eta") + i + "_pt1_rho" + j);						// Merge lowest two pt bins in forward (also low statistics, especially in non-CHS mode)
    }
  }

/*
  for(int j=0; j < 5; ++j){														// Set weights in the central region, high pt
    for(int i=10; i < bins.getNBins("pt"); ++i) bins.setWeights(TString("eta") + j + "_pt" + i + "_rho0", {5./3.,2./3.,2./3.});
  }
*/

  for(int j=0; j < bins.getNBins("eta"); ++j){														// Set weights in the central region, high pt
    for(int i=0; i < bins.getNBins("pt"); ++i){
      for(int k=0; k < bins.getNBins("rho"); ++k){
        if(weightOption == "a") bins.setWeights(TString("eta") + j + "_pt" + i + "_rho" + k, {1..,1..,1.});
        if(weightOption == "b") bins.setWeights(TString("eta") + j + "_pt" + i + "_rho" + k, {3.5/3.,2.75/3.,2.75/3.});
        if(weightOption == "c") bins.setWeights(TString("eta") + j + "_pt" + i + "_rho" + k, {4./3.,2.5/3.,2.5/3.});
        if(weightOption == "d") bins.setWeights(TString("eta") + j + "_pt" + i + "_rho" + k, {4.5/3.,2.25/3.,2.25/3.});
        if(weightOption == "e") bins.setWeights(TString("eta") + j + "_pt" + i + "_rho" + k, {5./3.,2./3.,2./3.});
        if(weightOption == "f") bins.setWeights(TString("eta") + j + "_pt" + i + "_rho" + k, {5.5/3.,1.75/3.,1.75/3.});
        if(weightOption == "g") bins.setWeights(TString("eta") + j + "_pt" + i + "_rho" + k, {6./3.,1.5/3.,1.5/3.});
      }
    }
  }
  return bins;
}

#include "binClass.h"

binClass getDefaultBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.3,1.5,2,2.5,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));				// i.e. 20 bins from 20 to 2000 with log=true and with an additional bin up to 6500
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.setBinRange("cbj",	"N_{#DeltaR<0.8}",	{-0.5,9999});
  bins.printBinRanges();

  for(int i=10; i < 21; ++i) bins.setLinks("eta5_pt9_rho0_cbj0", {TString::Format("eta5_pt%d_rho0_cbj0",i)});
  for(int i=14; i < 21; ++i) bins.setLinks("eta4_pt13_rho0_cbj0", {TString::Format("eta4_pt%d_rho0_cbj0",i)});
  for(int i=15; i < 21; ++i) bins.setLinks("eta3_pt14_rho0_cbj0", {TString::Format("eta3_pt%d_rho0_cbj0",i)});
  for(int i=17; i < 21; ++i) bins.setLinks("eta2_pt16_rho0_cbj0", {TString::Format("eta2_pt%d_rho0_cbj0",i)});
  for(int i=20; i < 21; ++i) bins.setLinks("eta1_pt19_rho0_cbj0", {TString::Format("eta1_pt%d_rho0_cbj0",i)});

  return bins;
}



binClass get8TeVBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,2.5,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));				// i.e. 20 bins from 20 to 2000 with log=true and with an additional bin up to 6500
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.setBinRange("cbj",	"N_{#DeltaR<0.8}",	{-0.5,9999});
  bins.printBinRanges();

  for(int i=10; i < 21; ++i) bins.setLinks("eta1_pt9_rho0_cbj0", {TString::Format("eta1_pt%d_rho0_cbj0",i)});

  return bins;
}



binClass getSmallEtaBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.3,1.5,1.8,2.1,2.4,2.7,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));				// i.e. 20 bins from 20 to 2000 with log=true and with an additional bin up to 6500
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.setBinRange("cbj",	"N_{#DeltaR<0.8}",	{-0.5,9999});
  bins.printBinRanges();

  for(int i=10; i < 21; ++i) bins.setLinks("eta7_pt9_rho0_cbj0", {TString::Format("eta7_pt%d_rho0_cbj0",i)});
  for(int i=14; i < 21; ++i) bins.setLinks("eta6_pt13_rho0_cbj0", {TString::Format("eta6_pt%d_rho0_cbj0",i)});
  for(int i=15; i < 21; ++i) bins.setLinks("eta5_pt14_rho0_cbj0", {TString::Format("eta5_pt%d_rho0_cbj0",i)});
  for(int i=17; i < 21; ++i) bins.setLinks("eta4_pt16_rho0_cbj0", {TString::Format("eta4_pt%d_rho0_cbj0",i)});
  for(int i=18; i < 21; ++i) bins.setLinks("eta3_pt17_rho0_cbj0", {TString::Format("eta3_pt%d_rho0_cbj0",i)});
  for(int i=19; i < 21; ++i) bins.setLinks("eta2_pt18_rho0_cbj0", {TString::Format("eta2_pt%d_rho0_cbj0",i)});
  for(int i=20; i < 21; ++i) bins.setLinks("eta1_pt19_rho0_cbj0", {TString::Format("eta1_pt%d_rho0_cbj0",i)});

  return bins;
}



binClass getClosebyJetBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.3,1.5,2.0,2.5,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));				// i.e. 20 bins from 20 to 2000 with log=true and with an additional bin up to 6500
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.setBinRange("cbj",	"N_{#DeltaR<0.8}",	bins.getBins(5, -0.5, 4.5, false, {6.5, 15.5}));			// if only cbj > 10 GeV
  bins.printBinRanges();

  for(int j=1; j < 7; ++j){															// No close-by jet binning in forward
    for(int i=0; i < 21; ++i) bins.setLinks(TString::Format("eta4_pt%d_rho0_cbj0", i), {TString::Format("eta4_pt%d_rho0_cbj%d",i,j)});
    for(int i=0; i < 21; ++i) bins.setLinks(TString::Format("eta5_pt%d_rho0_cbj0", i), {TString::Format("eta5_pt%d_rho0_cbj%d",i,j)});
  }

  for(int j=0; j < 7; ++j){
    for(int i=10; i < 21; ++i) bins.setLinks(TString::Format("eta5_pt9_rho0_cbj%d",j), {TString::Format("eta5_pt%d_rho0_cbj%d",i,j)});		// Merge high pt bins in forward
    for(int i=14; i < 21; ++i) bins.setLinks(TString::Format("eta4_pt13_rho0_cbj%d",j), {TString::Format("eta4_pt%d_rho0_cbj%d",i,j)});
    for(int i=15; i < 21; ++i) bins.setLinks(TString::Format("eta3_pt14_rho0_cbj%d",j), {TString::Format("eta3_pt%d_rho0_cbj%d",i,j)});
    for(int i=17; i < 21; ++i) bins.setLinks(TString::Format("eta2_pt16_rho0_cbj%d",j), {TString::Format("eta2_pt%d_rho0_cbj%d",i,j)});
    for(int i=18; i < 21; ++i) bins.setLinks(TString::Format("eta1_pt17_rho0_cbj%d",j), {TString::Format("eta1_pt%d_rho0_cbj%d",i,j)});
  }
  bins.setLinks("eta0_pt18_rho0_cbj0", {"eta0_pt19_rho0_cbj0", "eta0_pt20_rho0_cbj0"});
  bins.setLinks("eta1_pt17_rho0_cbj0", {"eta1_pt17_rho0_cbj1"});
  bins.setLinks("eta1_pt15_rho0_cbj0", {"eta1_pt16_rho0_cbj0", "eta1_pt15_rho0_cbj1", "eta1_pt16_rho0_cbj1"});
  bins.setLinks("eta1_pt13_rho0_cbj0", {"eta1_pt14_rho0_cbj0", "eta1_pt13_rho0_cbj1", "eta1_pt14_rho0_cbj1"});
  bins.setLinks("eta1_pt11_rho0_cbj0", {"eta1_pt12_rho0_cbj0", "eta1_pt11_rho0_cbj1", "eta1_pt12_rho0_cbj1"});
  bins.setLinks("eta1_pt9_rho0_cbj0",  {"eta1_pt10_rho0_cbj0", "eta1_pt9_rho0_cbj1",  "eta1_pt10_rho0_cbj1"});
  bins.setLinks("eta2_pt16_rho0_cbj0", {"eta2_pt16_rho0_cbj1"});
  bins.setLinks("eta2_pt14_rho0_cbj0", {"eta2_pt15_rho0_cbj0", "eta2_pt14_rho0_cbj1", "eta2_pt15_rho0_cbj1"});
  bins.setLinks("eta2_pt12_rho0_cbj0", {"eta2_pt13_rho0_cbj0", "eta2_pt12_rho0_cbj1", "eta3_pt13_rho0_cbj1"});
  bins.setLinks("eta2_pt10_rho0_cbj0", {"eta2_pt11_rho0_cbj0"});
  bins.setLinks("eta2_pt8_rho0_cbj0",  {"eta2_pt9_rho0_cbj0"});
  bins.setLinks("eta3_pt14_rho0_cbj0", {"eta3_pt14_rho0_cbj1"});
  bins.setLinks("eta3_pt11_rho0_cbj0", {"eta3_pt12_rho0_cbj0", "eta3_pt13_rho0_cbj0", "eta3_pt11_rho0_cbj1", "eta3_pt12_rho0_cbj1", "eta3_pt13_rho0_cbj1"});
  bins.setLinks("eta3_pt8_rho0_cbj0",  {"eta3_pt9_rho0_cbj0",  "eta3_pt10_rho0_cbj0", "eta3_pt8_rho0_cbj1",  "eta3_pt9_rho0_cbj1",  "eta3_pt10_rho0_cbj1"});
  bins.setLinks("eta3_pt12_rho0_cbj1", {"eta3_pt13_rho0_cbj1"});

  for(int j=0; j < 4; ++j){
    for(int i=2; i < 7; ++i) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj5", j,i), {TString::Format("eta%d_pt%d_rho0_cbj6",j,i)});	// At low pt low statistics in cbj binning, also not as important as at high pt
    for(int i=2; i < 5; ++i) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj2", j,i), {TString::Format("eta%d_pt%d_rho0_cbj4",j,i)});
    for(int i=2; i < 5; ++i) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj2", j,i), {TString::Format("eta%d_pt%d_rho0_cbj3",j,i)});
    for(int i=0; i < 5; ++i) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj0", j,i), {TString::Format("eta%d_pt%d_rho0_cbj1",j,i)});
    for(int i=0; i < 2; ++i) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj2", j,i), {TString::Format("eta%d_pt%d_rho0_cbj6",j,i)});
    for(int i=0; i < 2; ++i) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj2", j,i), {TString::Format("eta%d_pt%d_rho0_cbj5",j,i)});
    for(int i=0; i < 2; ++i) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj2", j,i), {TString::Format("eta%d_pt%d_rho0_cbj4",j,i)});
    for(int i=0; i < 2; ++i) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj2", j,i), {TString::Format("eta%d_pt%d_rho0_cbj3",j,i)});
  }
  bins.setLinks("eta1_pt0_rho0_cbj0", {"eta1_pt1_rho0_cbj0"});
  return bins;
}




binClass getAllClosebyJetBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.3,1.5,2,2.5,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));				// i.e. 20 bins from 20 to 2000 with log=true and with an additional bin up to 6500
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.setBinRange("cbj",	"N_{#DeltaR<0.8}",	{-0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 15.5});				// if all cbj
  bins.printBinRanges();

  for(int j=1; j < 7; ++j){															// No close-by jet binning in forward
    for(int i=0; i < 21; ++i) bins.setLinks(TString::Format("eta4_pt%d_rho0_cbj0", i), {TString::Format("eta4_pt%d_rho0_cbj%d",i,j)});
    for(int i=0; i < 21; ++i) bins.setLinks(TString::Format("eta5_pt%d_rho0_cbj0", i), {TString::Format("eta5_pt%d_rho0_cbj%d",i,j)});

    for(int i=10; i < 21; ++i) bins.setLinks(TString::Format("eta5_pt9_rho0_cbj%d",j), {TString::Format("eta5_pt%d_rho0_cbj%d",i,j)});		// Merge high pt bins in forward
    for(int i=14; i < 21; ++i) bins.setLinks(TString::Format("eta4_pt13_rho0_cbj%d",j), {TString::Format("eta4_pt%d_rho0_cbj%d",i,j)});
    for(int i=15; i < 21; ++i) bins.setLinks(TString::Format("eta3_pt14_rho0_cbj%d",j), {TString::Format("eta3_pt%d_rho0_cbj%d",i,j)});
    for(int i=17; i < 21; ++i) bins.setLinks(TString::Format("eta2_pt16_rho0_cbj%d",j), {TString::Format("eta2_pt%d_rho0_cbj%d",i,j)});
    for(int i=18; i < 21; ++i) bins.setLinks(TString::Format("eta1_pt17_rho0_cbj%d",j), {TString::Format("eta1_pt%d_rho0_cbj%d",i,j)});
  }
  for(int j=0; j < 4; ++j){
    for(int i=2; i < 7; ++i) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj5", j,i), {TString::Format("eta%d_pt%d_rho0_cbj6",j,i)});	// At low pt low statistics in cbj binning, also not as important as at high pt
    for(int i=2; i < 5; ++i) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj2", j,i), {TString::Format("eta%d_pt%d_rho0_cbj4",j,i)});
    for(int i=2; i < 5; ++i) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj2", j,i), {TString::Format("eta%d_pt%d_rho0_cbj3",j,i)});
    for(int i=0; i < 5; ++i) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj0", j,i), {TString::Format("eta%d_pt%d_rho0_cbj1",j,i)});
    for(int i=0; i < 2; ++i) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj2", j,i), {TString::Format("eta%d_pt%d_rho0_cbj6",j,i)});
    for(int i=0; i < 2; ++i) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj2", j,i), {TString::Format("eta%d_pt%d_rho0_cbj5",j,i)});
    for(int i=0; i < 2; ++i) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj2", j,i), {TString::Format("eta%d_pt%d_rho0_cbj4",j,i)});
    for(int i=0; i < 2; ++i) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj2", j,i), {TString::Format("eta%d_pt%d_rho0_cbj3",j,i)});
  }
  return bins;
}



binClass getPtDoubleConeBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.3,1.5,2.0,2.5,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));				// i.e. 20 bins from 20 to 2000 with log=true and with an additional bin up to 6500
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.setBinRange("cbj",	"pt(doubleCone)",	bins.getBins(20, 20, 2000, true, {6500}));				// if pt(doubleCone)
  bins.printBinRanges();

  // Link some bins to be merged because of low statistics (for example higher pT bins at large eta)
  for(int j=1; j < 21; ++j){															// No close-by jet binning in forward
    for(int i=0; i < 21; ++i) bins.setLinks(TString::Format("eta4_pt%d_rho0_cbj0", i), {TString::Format("eta4_pt%d_rho0_cbj%d",i,j)});
    for(int i=0; i < 21; ++i) bins.setLinks(TString::Format("eta5_pt%d_rho0_cbj0", i), {TString::Format("eta5_pt%d_rho0_cbj%d",i,j)});
  }
  for(int j=0; j < 21; ++j){
    for(int i=10; i < 21; ++i) bins.setLinks(TString::Format("eta5_pt9_rho0_cbj%d",j), {TString::Format("eta5_pt%d_rho0_cbj%d",i,j)});
    for(int i=14; i < 21; ++i) bins.setLinks(TString::Format("eta4_pt13_rho0_cbj%d",j), {TString::Format("eta4_pt%d_rho0_cbj%d",i,j)});
    for(int i=15; i < 21; ++i) bins.setLinks(TString::Format("eta3_pt14_rho0_cbj%d",j), {TString::Format("eta3_pt%d_rho0_cbj%d",i,j)});
    for(int i=17; i < 21; ++i) bins.setLinks(TString::Format("eta2_pt16_rho0_cbj%d",j), {TString::Format("eta2_pt%d_rho0_cbj%d",i,j)});
    for(int i=18; i < 21; ++i) bins.setLinks(TString::Format("eta1_pt17_rho0_cbj%d",j), {TString::Format("eta1_pt%d_rho0_cbj%d",i,j)});
  }
  for(int k=0; k < 4; ++k){
    for(int i=0; i < 21; ++i){
      for(int j=0; j < i; ++j) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj%d",k,i,i), {TString::Format("eta%d_pt%d_rho0_cbj%d",k,i,j)});
      if(i > 5 && i+6-k < 21) for(int j = i+6-k; j < 21; ++j) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj%d",k,i,i+5-k), {TString::Format("eta%d_pt%d_rho0_cbj%d",k,i,j)});
      else if(i < 6) for(int j=9; j < 21; ++j) bins.setLinks(TString::Format("eta%d_pt%d_rho0_cbj%d",k,i,i+5), {TString::Format("eta%d_pt%d_rho0_cbj%d",k,i,j)});
    }
    bins.setLinks(TString::Format("eta%d_pt0_rho0_cbj0",k), {TString::Format("eta%d_pt0_rho0_cbj1",k), TString::Format("eta%d_pt0_rho0_cbj2",k), TString::Format("eta%d_pt0_rho0_cbj3",k)});
    bins.setLinks(TString::Format("eta%d_pt1_rho0_cbj1",k), {TString::Format("eta%d_pt1_rho0_cbj2",k), TString::Format("eta%d_pt1_rho0_cbj3",k)});
    bins.setLinks(TString::Format("eta%d_pt2_rho0_cbj2",k), {TString::Format("eta%d_pt2_rho0_cbj3",k)});
    bins.setLinks("eta3_pt3_rho0_cbj3", {"eta3_pt3_rho0_cbj4"});
  }
  return bins;
}

#include "binClass.h"

binClass getDefaultBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.3,1.5,2,2.5,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));				// i.e. 20 bins from 20 to 2000 with log=true and with an additional bin up to 6500
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.printBinRanges();

  for(int i=10; i < 21; ++i) bins.setLinks("eta5_pt9_rho0", {TString::Format("eta5_pt%d_rho0",i)});
  for(int i=14; i < 21; ++i) bins.setLinks("eta4_pt13_rho0", {TString::Format("eta4_pt%d_rho0",i)});
  for(int i=15; i < 21; ++i) bins.setLinks("eta3_pt14_rho0", {TString::Format("eta3_pt%d_rho0",i)});
  for(int i=17; i < 21; ++i) bins.setLinks("eta2_pt16_rho0", {TString::Format("eta2_pt%d_rho0",i)});
  for(int i=20; i < 21; ++i) bins.setLinks("eta1_pt19_rho0", {TString::Format("eta1_pt%d_rho0",i)});

  return bins;
}



binClass get8TeVBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,2.5,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));				// i.e. 20 bins from 20 to 2000 with log=true and with an additional bin up to 6500
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.printBinRanges();

  for(int i=10; i < 21; ++i) bins.setLinks("eta1_pt9_rho0", {TString::Format("eta1_pt%d_rho0",i)});

  return bins;
}



binClass getSmallEtaBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.3,1.5,1.8,2.1,2.4,2.7,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));				// i.e. 20 bins from 20 to 2000 with log=true and with an additional bin up to 6500
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.printBinRanges();

  for(int i=10; i < 21; ++i) bins.setLinks("eta7_pt9_rho0", {TString::Format("eta7_pt%d_rho0",i)});
  for(int i=14; i < 21; ++i) bins.setLinks("eta6_pt13_rho0", {TString::Format("eta6_pt%d_rho0",i)});
  for(int i=15; i < 21; ++i) bins.setLinks("eta5_pt14_rho0", {TString::Format("eta5_pt%d_rho0",i)});
  for(int i=17; i < 21; ++i) bins.setLinks("eta4_pt16_rho0", {TString::Format("eta4_pt%d_rho0",i)});
  for(int i=18; i < 21; ++i) bins.setLinks("eta3_pt17_rho0", {TString::Format("eta3_pt%d_rho0",i)});
  for(int i=19; i < 21; ++i) bins.setLinks("eta2_pt18_rho0", {TString::Format("eta2_pt%d_rho0",i)});
  for(int i=20; i < 21; ++i) bins.setLinks("eta1_pt19_rho0", {TString::Format("eta1_pt%d_rho0",i)});

  return bins;
}

#include "binClass.h"

/*
 * Extended binning configurations, for testing additional jet categories
 */
binClass getSmallEtaBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.3,1.5,1.8,2.1,2.4,2.7,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.setBinRange("aj",	"additional jets",	{-0.5,9999});
  bins.printBinRanges();

  for(int j=10; j < bins.getNBins("pt"); ++j) bins.setLink("eta7_pt9_rho0_aj0",  TString("eta7_pt") + j + "_rho0_aj0");
  for(int j=14; j < bins.getNBins("pt"); ++j) bins.setLink("eta6_pt13_rho0_aj0", TString("eta6_pt") + j + "_rho0_aj0");
  for(int j=15; j < bins.getNBins("pt"); ++j) bins.setLink("eta5_pt14_rho0_aj0", TString("eta5_pt") + j + "_rho0_aj0");
  for(int j=17; j < bins.getNBins("pt"); ++j) bins.setLink("eta4_pt16_rho0_aj0", TString("eta4_pt") + j + "_rho0_aj0");
  for(int j=18; j < bins.getNBins("pt"); ++j) bins.setLink("eta3_pt17_rho0_aj0", TString("eta3_pt") + j + "_rho0_aj0");
  for(int j=19; j < bins.getNBins("pt"); ++j) bins.setLink("eta2_pt18_rho0_aj0", TString("eta2_pt") + j + "_rho0_aj0");
  for(int j=20; j < bins.getNBins("pt"); ++j) bins.setLink("eta1_pt19_rho0_aj0", TString("eta1_pt") + j + "_rho0_aj0");

  return bins;
}


binClass getNoBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		{20,6500});
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.setBinRange("aj",	"additional jets",	{-0.5,9999});
  bins.printBinRanges();
  return bins;
}


binClass getAdditionalJetBinning(){													// Same as V1 + additional jets binning
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.3,1.5,1.8,2.1,2.4,2.7,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));
  bins.setBinRange("rho",	"#rho",			{0,15,20,25,30,9999});
  bins.setBinRange("aj",	"additional jets",	{-0.5,0.5,9999});								// Adding categories n = 0 and n > 0, with n the number of additional jets	
  bins.printBinRanges();

  for(int i=0; i < bins.getNBins("eta"); ++i){												// Use only additional jet binnings in the central bins (eta < 2.4), where enough statistics is available
   for(int j=0; j < bins.getNBins("pt"); ++j){												// Low pt and high pt bins also get merged because of too low statistics
    for(int k=0; k < bins.getNBins("rho"); ++k){
     if(i == 0 && j > 2 && j < 20) continue;
     if(i == 1 && j > 2 && j < 18) continue;
     if(i == 2 && j > 3 && j < 17) continue;
     if(i == 3 && j > 3 && j < 15) continue;
     if(i == 4 && j > 3 && j < 15) continue;
     bins.setLink(TString("eta")+i+"_pt"+j+"_rho"+k+"_aj0", TString("eta")+i+"_pt"+j+"_rho"+k+"_aj1");
    }
   }
  }


  for(int i=0; i < 5; ++i){														// No rho binning needed in the central (eta < 2.4) as rho dependencies are very small
   for(int j=0; j < bins.getNBins("pt"); ++j){												// On the other hand some bins (medium pt) have a lot of statistics, so it would probably be possible to
    for(int k=1; k < bins.getNBins("rho"); ++k){											// use also some rho binning in the medium pt (and correcting for even the smallest dependency) without
     for(int l=0; l < bins.getNBins("aj"); ++l){											// giving up on smooth pdf's
      bins.setLink(TString("eta")+i+"_pt"+j+"_rho0_aj"+l, TString("eta")+i+"_pt"+j+"_rho"+k+"_aj"+l);
     }
    }
   }
  }

  for(int i=5; i < bins.getNBins("eta"); ++i){												// Unfortunately, need to merge some rho bins in lowest pt bins of forward because of very low statistics
   for(int j=0; j < 3; ++j){														// If we get a forward jets QCD sample, try to go with finer binning, it could help a lot to boost performance
    for(int k=3; k < 5; ++k){														// Also there is some overtraining issues for the lowest and highest rho bins in some pt bins which
     for(int l=0; l < bins.getNBins("aj"); ++l){
      bins.setLink(TString("eta")+i+"_pt"+j+"_rho2_aj"+l, TString("eta")+i+"_pt"+j+"_rho"+k+"_aj"+l);					//is again something which can be solved by more statistics
     }
    }
   }
  }
  for(int l=0; l < bins.getNBins("aj"); ++l){
   bins.setLink(TString("eta5_pt13_rho3_aj")+l, TString("eta5_pt13_rho4_aj")+l);							// A couple of other low statistics bins
   bins.setLink(TString("eta6_pt3_rho3_aj")+l,  TString("eta6_pt3_rho4_aj")+l);
   bins.setLink(TString("eta6_pt12_rho3_aj")+l, TString("eta6_pt12_rho4_aj")+l);
  }

  for(int k=0; k < bins.getNBins("rho"); ++k){
   for(int l=0; l < bins.getNBins("aj"); ++l){
    for(int j=9;  j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta7_pt8_rho") +k+"_aj"+l, TString("eta7_pt")+j+"_rho"+k+"_aj"+l);// Merge high pt bins in forward (no or low statistics)
    for(int j=13; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta6_pt12_rho")+k+"_aj"+l, TString("eta6_pt")+j+"_rho"+k+"_aj"+l);
    for(int j=14; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta5_pt13_rho")+k+"_aj"+l, TString("eta5_pt")+j+"_rho"+k+"_aj"+l);
    for(int j=16; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta4_pt15_rho")+k+"_aj"+l, TString("eta4_pt")+j+"_rho"+k+"_aj"+l);
    for(int j=17; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta3_pt16_rho")+k+"_aj"+l, TString("eta3_pt")+j+"_rho"+k+"_aj"+l);
    for(int j=18; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta2_pt17_rho")+k+"_aj"+l, TString("eta2_pt")+j+"_rho"+k+"_aj"+l);
    for(int j=19; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta1_pt18_rho")+k+"_aj"+l, TString("eta1_pt")+j+"_rho"+k+"_aj"+l);
   
    for(int i = 3; i < bins.getNBins("eta"); ++i){
     bins.setLink(TString("eta")+i+"_pt0_rho"+k+"_aj"+l, TString("eta")+i+"_pt1_rho"+k+"_aj"+l);					// Merge lowest two pt bins in forward (also low statistics, especially in non-CHS mode)
    }
   }
  }

  // The weights
  // The multiplicity is less correlated with the other two variables and has better discrimination power than those two at high pt
  // Therefore we give a higher weight to the multiplicity contribution to increase the performance
  // The modified likelihood algorithm has L/(1-L) = Mult((Qi/Gi)^wi), the original algorithm has wi = 1 for all three variables
  // Note 1: We take the sum of the three weights always as 3, this ensures we keep the same final likelihood value if all three
  //         variables have the same Qi/Gi regardeless of weights (i.e. if all three variables agree on their individual probabilities
  //         the weights have not much impact, if they disagree the variable with higher weight has priority could pull it more in
  //         his direction)
  // Note 2: The multiplicity is a discrete variable, and therefore its contribution to the likelihood is also discrete. If the
  //         multiplicity weight is to high, this discrete structure could show up in the likelihood shape. However, at the high pt
  //         bins this effect is only observable if you go higher than {6/3,1.5/3,1.5/3}. For the lower pt bins, this is threshold
  //         is lower, because we have less multiplicity bins (but here there is also no advantage of using a higher weight for
  //         the multiplicity)
  // Note 3: We use the same weights for the 5 most central eta bins (<2.4). It turns out we can give the same weights to these 5
  //         neighbouring eta bins. In the forward we use always {1,1,1}, as the multiplicity is in general less performing and due
  //         to the low statistics it is difficult to find the optimal weights. (Though if there would be more statistics in the forward,
  //         it could be useful to re-check for the higher pt bins)
  for(int i=0; i < 5; ++i){
   for(int k=0; k < bins.getNBins("rho"); ++k){
    for(int l = 0; l < bins.getNBins("aj"); ++l){
     for(int j=12; j < 21; ++j) bins.setWeights(TString("eta")+i+"_pt"+j+"_rho"+k+"_aj"+l, {6./3.,1.5/3.,1.5/3.});
     for(int j=10; j < 12; ++j) bins.setWeights(TString("eta")+i+"_pt"+j+"_rho"+k+"_aj"+l, {5.5/3.,1.75/3.,1.75/3.});
     for(int j=8 ; j < 10; ++j) bins.setWeights(TString("eta")+i+"_pt"+j+"_rho"+k+"_aj"+l, {5./3.,2./3.,2./3.});
     for(int j=7 ; j < 8;  ++j) bins.setWeights(TString("eta")+i+"_pt"+j+"_rho"+k+"_aj"+l, {4.5/3.,2.25/3.,2.25/3.});
     for(int j=6 ; j < 7;  ++j) bins.setWeights(TString("eta")+i+"_pt"+j+"_rho"+k+"_aj"+l, {4./3.,2.5/3.,2.5/3.});
     for(int j=5 ; j < 6;  ++j) bins.setWeights(TString("eta")+i+"_pt"+j+"_rho"+k+"_aj"+l, {3.5/3.,2.75/3.,2.75/3.});
    }
   }
  }

  return bins;
}



binClass getV1Binning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.3,1.5,1.8,2.1,2.4,2.7,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));
  bins.setBinRange("rho",	"#rho",			{0,15,20,25,30,9999});
  bins.setBinRange("aj",	"additional jets",	{-0.5,9999});									// Adding categories n = 0 and n > 0, with n the number of additional jets	
  bins.printBinRanges();

  for(int i=0; i < 5; ++i){														// No rho binning needed in the central (eta < 2.4) as rho dependencies are very small
   for(int j=0; j < bins.getNBins("pt"); ++j){												// On the other hand some bins (medium pt) have a lot of statistics, so it would probably be possible to
    for(int k=1; k < bins.getNBins("rho"); ++k){											// use also some rho binning in the medium pt (and correcting for even the smallest dependency) without
     for(int l=0; l < bins.getNBins("aj"); ++l){											// giving up on smooth pdf's
      bins.setLink(TString("eta")+i+"_pt"+j+"_rho0_aj"+l, TString("eta")+i+"_pt"+j+"_rho"+k+"_aj"+l);
     }
    }
   }
  }

  for(int i=5; i < bins.getNBins("eta"); ++i){												// Unfortunately, need to merge some rho bins in lowest pt bins of forward because of very low statistics
   for(int j=0; j < 3; ++j){														// If we get a forward jets QCD sample, try to go with finer binning, it could help a lot to boost performance
    for(int k=3; k < 5; ++k){														// Also there is some overtraining issues for the lowest and highest rho bins in some pt bins which
     for(int l=0; l < bins.getNBins("aj"); ++l){
      bins.setLink(TString("eta")+i+"_pt"+j+"_rho2_aj"+l, TString("eta")+i+"_pt"+j+"_rho"+k+"_aj"+l);					//is again something which can be solved by more statistics
     }
    }
   }
  }
  for(int l=0; l < bins.getNBins("aj"); ++l){
   bins.setLink(TString("eta5_pt13_rho3_aj")+l, TString("eta5_pt13_rho4_aj")+l);							// A couple of other low statistics bins
   bins.setLink(TString("eta6_pt3_rho3_aj")+l,  TString("eta6_pt3_rho4_aj")+l);
   bins.setLink(TString("eta6_pt12_rho3_aj")+l, TString("eta6_pt12_rho4_aj")+l);
  }

  for(int k=0; k < bins.getNBins("rho"); ++k){
   for(int l=0; l < bins.getNBins("aj"); ++l){
    for(int j=9;  j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta7_pt8_rho") +k+"_aj"+l, TString("eta7_pt")+j+"_rho"+k+"_aj"+l);// Merge high pt bins in forward (no or low statistics)
    for(int j=13; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta6_pt12_rho")+k+"_aj"+l, TString("eta6_pt")+j+"_rho"+k+"_aj"+l);
    for(int j=14; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta5_pt13_rho")+k+"_aj"+l, TString("eta5_pt")+j+"_rho"+k+"_aj"+l);
    for(int j=16; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta4_pt15_rho")+k+"_aj"+l, TString("eta4_pt")+j+"_rho"+k+"_aj"+l);
    for(int j=17; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta3_pt16_rho")+k+"_aj"+l, TString("eta3_pt")+j+"_rho"+k+"_aj"+l);
    for(int j=18; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta2_pt17_rho")+k+"_aj"+l, TString("eta2_pt")+j+"_rho"+k+"_aj"+l);
    for(int j=19; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta1_pt18_rho")+k+"_aj"+l, TString("eta1_pt")+j+"_rho"+k+"_aj"+l);
   
    for(int i = 3; i < bins.getNBins("eta"); ++i){
     bins.setLink(TString("eta")+i+"_pt0_rho"+k+"_aj"+l, TString("eta")+i+"_pt1_rho"+k+"_aj"+l);					// Merge lowest two pt bins in forward (also low statistics, especially in non-CHS mode)
    }
   }
  }

  // The weights
  // The multiplicity is less correlated with the other two variables and has better discrimination power than those two at high pt
  // Therefore we give a higher weight to the multiplicity contribution to increase the performance
  // The modified likelihood algorithm has L/(1-L) = Mult((Qi/Gi)^wi), the original algorithm has wi = 1 for all three variables
  // Note 1: We take the sum of the three weights always as 3, this ensures we keep the same final likelihood value if all three
  //         variables have the same Qi/Gi regardeless of weights (i.e. if all three variables agree on their individual probabilities
  //         the weights have not much impact, if they disagree the variable with higher weight has priority could pull it more in
  //         his direction)
  // Note 2: The multiplicity is a discrete variable, and therefore its contribution to the likelihood is also discrete. If the
  //         multiplicity weight is to high, this discrete structure could show up in the likelihood shape. However, at the high pt
  //         bins this effect is only observable if you go higher than {6/3,1.5/3,1.5/3}. For the lower pt bins, this is threshold
  //         is lower, because we have less multiplicity bins (but here there is also no advantage of using a higher weight for
  //         the multiplicity)
  // Note 3: We use the same weights for the 5 most central eta bins (<2.4). It turns out we can give the same weights to these 5
  //         neighbouring eta bins. In the forward we use always {1,1,1}, as the multiplicity is in general less performing and due
  //         to the low statistics it is difficult to find the optimal weights. (Though if there would be more statistics in the forward,
  //         it could be useful to re-check for the higher pt bins)
  for(int i=0; i < 5; ++i){
   for(int k=0; k < bins.getNBins("rho"); ++k){
    for(int l = 0; l < bins.getNBins("aj"); ++l){
     for(int j=12; j < 21; ++j) bins.setWeights(TString("eta")+i+"_pt"+j+"_rho"+k+"_aj"+l, {6./3.,1.5/3.,1.5/3.});
     for(int j=10; j < 12; ++j) bins.setWeights(TString("eta")+i+"_pt"+j+"_rho"+k+"_aj"+l, {5.5/3.,1.75/3.,1.75/3.});
     for(int j=8 ; j < 10; ++j) bins.setWeights(TString("eta")+i+"_pt"+j+"_rho"+k+"_aj"+l, {5./3.,2./3.,2./3.});
     for(int j=7 ; j < 8;  ++j) bins.setWeights(TString("eta")+i+"_pt"+j+"_rho"+k+"_aj"+l, {4.5/3.,2.25/3.,2.25/3.});
     for(int j=6 ; j < 7;  ++j) bins.setWeights(TString("eta")+i+"_pt"+j+"_rho"+k+"_aj"+l, {4./3.,2.5/3.,2.5/3.});
     for(int j=5 ; j < 6;  ++j) bins.setWeights(TString("eta")+i+"_pt"+j+"_rho"+k+"_aj"+l, {3.5/3.,2.75/3.,2.75/3.});
    }
   }
  }

  return bins;
}

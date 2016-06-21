#include "binClass.h"

/*
 * Finer eta binning
 */
binClass getSmallEtaBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.3,1.5,1.8,2.1,2.4,2.7,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.printBinRanges();

  for(int j=10; j < bins.getNBins("pt"); ++j) bins.setLink("eta7_pt9_rho0",  TString("eta7_pt") + j + "_rho0");
  for(int j=14; j < bins.getNBins("pt"); ++j) bins.setLink("eta6_pt13_rho0", TString("eta6_pt") + j + "_rho0");
  for(int j=15; j < bins.getNBins("pt"); ++j) bins.setLink("eta5_pt14_rho0", TString("eta5_pt") + j + "_rho0");
  for(int j=17; j < bins.getNBins("pt"); ++j) bins.setLink("eta4_pt16_rho0", TString("eta4_pt") + j + "_rho0");
  for(int j=18; j < bins.getNBins("pt"); ++j) bins.setLink("eta3_pt17_rho0", TString("eta3_pt") + j + "_rho0");
  for(int j=19; j < bins.getNBins("pt"); ++j) bins.setLink("eta2_pt18_rho0", TString("eta2_pt") + j + "_rho0");
  for(int j=20; j < bins.getNBins("pt"); ++j) bins.setLink("eta1_pt19_rho0", TString("eta1_pt") + j + "_rho0");

  return bins;
}

/*
 * For plot
 */
binClass getTestBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,2.0,3.0,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		{30,40,80,100,10e4});
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.printBinRanges();
  return bins;
}

/*
 * Similar to old V1 binning (as used for CSA14 training), but adapted for Spring15dr asympt25ns flat QCD training 
 */
binClass getV2Binning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.3,1.5,1.8,2.1,2.4,2.7,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));
  bins.setBinRange("rho",	"#rho",			{0,8,11,14,17,9999});								// Very preliminary binning for spring15dr asympt25ns rho distribution (peaks around ~12)
  bins.printBinRanges();

  for(int i=0; i < 5; ++i){														// No rho binning needed in the central (eta < 2.4) as rho dependencies are very small
    for(int j=0; j < bins.getNBins("pt"); ++j){
      for(int k=1; k < bins.getNBins("rho"); ++k){
        bins.setLink(TString("eta") + i + "_pt" + j + "_rho0", TString("eta") + i + "_pt" + j + "_rho" + k);
      }
    }
  }

  for(int i=5; i < bins.getNBins("eta"); ++i){												// Unfortunately, need to merge some rho bins in lowest pt bins of forward because of very low statistics
    for(int j=0; j < 3; ++j){														// If we get a forward jets QCD sample, try to go with finer binning, it could help a lot to boost performance
      bins.setLink(TString("eta") + i + "_pt" + j + "_rho2", TString("eta") + i + "_pt" + j + "_rho3");					// Also there are some overtraining issues for the lowest and highest rho bins in some pt bins which
      bins.setLink(TString("eta") + i + "_pt" + j + "_rho2", TString("eta") + i + "_pt" + j + "_rho4");					// is again something which can be solved by more statistics
    }
  }

  for(int k=0; k < bins.getNBins("rho"); ++k){
    for(int j=9;  j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta7_pt8_rho")  + k, TString("eta7_pt") + j + "_rho" + k);	// Merge high pt bins jn forward (no or low statistics)
    for(int j=13; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta6_pt12_rho") + k, TString("eta6_pt") + j + "_rho" + k);
    for(int j=14; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta5_pt13_rho") + k, TString("eta5_pt") + j + "_rho" + k);
    for(int j=16; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta4_pt15_rho") + k, TString("eta4_pt") + j + "_rho" + k);
    for(int j=17; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta3_pt16_rho") + k, TString("eta3_pt") + j + "_rho" + k);
    for(int j=18; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta2_pt17_rho") + k, TString("eta2_pt") + j + "_rho" + k);
    for(int j=19; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta1_pt18_rho") + k, TString("eta1_pt") + j + "_rho" + k);
   
    for(int i = 3; i < bins.getNBins("eta"); ++i){
      bins.setLink(TString("eta") + i + "_pt0_rho" + k, TString("eta") + i + "_pt1_rho" + k);						// Merge lowest two pt bins in forward (also low statistics, especially in non-CHS mode)
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
  for(int k=0; k < bins.getNBins("rho"); ++k){
    for(int i=0; i < 5; ++i){
      for(int j=12; j < 21; ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {6./3.,1.5/3.,1.5/3.});
      for(int j=10; j < 12; ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {5.5/3.,1.75/3.,1.75/3.});
      for(int j=8 ; j < 10; ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {5./3.,2./3.,2./3.});
      for(int j=7 ; j < 8;  ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {4.5/3.,2.25/3.,2.25/3.});
      for(int j=6 ; j < 7;  ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {4./3.,2.5/3.,2.5/3.});
      for(int j=5 ; j < 6;  ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {3.5/3.,2.75/3.,2.75/3.});
    }
  }
  return bins;
}


/*
 * Experimental binning (formerly called v2): try to introduce overlapping bins to get smoother transitions
 */
binClass getExperimentalBinning(bool PU40bx50 = false){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3,4.7});						// Small bins in the HF region in combination with the neighbour method result in smooth transitions
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));
  if(PU40bx50) bins.setBinRange("rho",	"#rho",		{0,21,22.5,24,25.5,27,9999});								// For spring14dr PU40bx50, should be different in Phys14
  else         bins.setBinRange("rho",	"#rho",		{0,16,17,18,19,20,9999});								// For spring14dr PU20bx25, should be different in Phys14
  bins.printBinRanges();

  for(int i=0; i < 8; ++i){															// No rho binning needed in the central (eta < 2.4) as rho dependencies are very small
    for(int j=0; j < bins.getNBins("pt"); ++j){
      for(int k=1; k < bins.getNBins("rho"); ++k){
        bins.setLink(TString("eta") + i + "_pt" + j + "_rho0", TString("eta") + i + "_pt" + j + "_rho" + k);
      }
    }
  }


  for(int k=0; k < bins.getNBins("rho"); ++k){
    for(int j=9;  j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta10_pt8_rho")  + k, TString("eta10_pt") + j + "_rho" + k);		// Merge higher pt bins in forward
    for(int j=13; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta9_pt12_rho") + k, TString("eta9_pt") + j + "_rho" + k);
    for(int j=14; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta8_pt13_rho") + k, TString("eta8_pt") + j + "_rho" + k);
    for(int j=15; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta7_pt14_rho") + k, TString("eta7_pt") + j + "_rho" + k);
    for(int j=16; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta6_pt15_rho") + k, TString("eta6_pt") + j + "_rho" + k);
    for(int j=17; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta5_pt16_rho") + k, TString("eta5_pt") + j + "_rho" + k);
    for(int j=18; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta4_pt17_rho") + k, TString("eta4_pt") + j + "_rho" + k);
    for(int j=19; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta3_pt18_rho") + k, TString("eta3_pt") + j + "_rho" + k);
    for(int j=20; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta2_pt19_rho") + k, TString("eta2_pt") + j + "_rho" + k);

    for(int i = 8; i < bins.getNBins("eta"); ++i){
      bins.setLink(TString("eta") + i + "_pt0_rho" + k, TString("eta") + i + "_pt1_rho" + k);							// Merge lowest two pt bins
    }
  }


  // Neighbour bins: those are a bit different than the "linked" (a.k.a "merged") bins,
  // they can be used to add the statistics of other bins, which allow to go for a finer binning, with smoother transitions,
  // but still keeping enough statistics
  for(int i=0; i < bins.getNBins("eta"); ++i){
    for(int j=0; j < bins.getNBins("pt"); ++j){
      for(int k=0; k < bins.getNBins("rho"); ++k){
        for(int ii : {-1,0,1}){															// Loop over neighbours
          for(int jj : {-1,0,1}){
            for(int kk : {-1,0,1}){
              if(ii == 0 && jj == 0 && kk == 0) continue;											// Exclude the bin itself
              if((i+ii < 0) || (j+jj < 0) || (k+kk < 0)) continue;										// Range check
              if((i+ii >= bins.getNBins("eta")) || (j+jj >= bins.getNBins("pt")) || (k+kk >= bins.getNBins("rho"))) continue;
              if((i+ii >= bins.getNBins("eta")-2) && ii != 0) continue;										// For eta: don't add statistics from neighbours to the 2.7-3 bin, as this one is too different from its neighbours
              bins.setNeighbour(TString("eta") + i + "_pt" + j + "_rho" + k, TString("eta") + (i+ii) + "_pt" + (j+jj) + "_rho" + (k+kk));
            }
          }
        }
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
  // Note 3: We use the same weights for the 8 most central eta bins (<2.5). It turns out we can give the same weights to these 5
  //         neighbouring eta bins. In the forward we use always {1,1,1}, as the multiplicity is in general less performing and due
  //         to the low statistics it is difficult to find the optimal weights. (Though if there would be more statistics in the forward,
  //         it could be useful to re-check for the higher pt bins)
  for(int k=0; k < bins.getNBins("rho"); ++k){
    for(int i=0; i < 8; ++i){
      for(int j=12; j < 21; ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {6./3.,1.5/3.,1.5/3.});
      for(int j=10; j < 12; ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {5.5/3.,1.75/3.,1.75/3.});
      for(int j=8 ; j < 10; ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {5./3.,2./3.,2./3.});
      for(int j=7 ; j < 8;  ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {4.5/3.,2.25/3.,2.25/3.});
      for(int j=6 ; j < 7;  ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {4./3.,2.5/3.,2.5/3.});
      for(int j=5 ; j < 6;  ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {3.5/3.,2.75/3.,2.75/3.});
    }
  }
  return bins;
}


/*
 * To use no binning at all
 */
binClass getNoBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		{20,6500});
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.printBinRanges();
  return bins;
}


/*
 * No eta binning, but using pt slices
 */
binClass getCentralPtSlices(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,2.4});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));
  bins.setBinRange("rho",	"#rho",			{0,9999});
  bins.printBinRanges();
  return bins;
}


/*
 * Similar to old V1 binning (as used for CSA14 training), but adapted for 76X asympt25ns pthat bin training
 */
binClass get76XBinning(){
  binClass bins;
  bins.setBinRange("eta", 	"#eta",			{0,1.3,1.5,1.8,2.1,2.4,2.7,3,4.7});
  bins.setBinRange("pt" , 	"p_{T}",		bins.getBins(20, 20, 2000, true, {6500}));
  bins.setBinRange("rho",	"#rho",			{0,8,11,14,17,9999});								// Very preliminary binning for spring15dr asympt25ns rho distribution (peaks around ~12)
  bins.printBinRanges();

  for(int i=0; i < 5; ++i){														// No rho binning needed in the central (eta < 2.4) as rho dependencies are very small
    for(int j=0; j < bins.getNBins("pt"); ++j){
      for(int k=1; k < bins.getNBins("rho"); ++k){
        bins.setLink(TString("eta") + i + "_pt" + j + "_rho0", TString("eta") + i + "_pt" + j + "_rho" + k);
      }
    }
  }

  for(int i=5; i < bins.getNBins("eta"); ++i){												// Unfortunately, need to merge some rho bins in many pt bins of forward because of very low statistics
    for(int j=0; j < bins.getNBins("pt"); ++j){														// If we get a forward jets QCD sample, try to go with finer binning, it could help a lot to boost performance
      bins.setLink(TString("eta") + i + "_pt" + j + "_rho2", TString("eta") + i + "_pt" + j + "_rho3");					// Also there are some overtraining issues for the lowest and highest rho bins in some pt bins which
      bins.setLink(TString("eta") + i + "_pt" + j + "_rho2", TString("eta") + i + "_pt" + j + "_rho4");					// is again something which can be solved by more statistics
    }
  }

  for(int k=0; k < bins.getNBins("rho"); ++k){
    for(int j=9;  j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta7_pt8_rho")  + k, TString("eta7_pt") + j + "_rho" + k);	// Merge high pt bins jn forward (no or low statistics)
    for(int j=13; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta6_pt12_rho") + k, TString("eta6_pt") + j + "_rho" + k);
    for(int j=14; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta5_pt13_rho") + k, TString("eta5_pt") + j + "_rho" + k);
    for(int j=16; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta4_pt15_rho") + k, TString("eta4_pt") + j + "_rho" + k);
    for(int j=17; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta3_pt16_rho") + k, TString("eta3_pt") + j + "_rho" + k);
    for(int j=18; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta2_pt17_rho") + k, TString("eta2_pt") + j + "_rho" + k);
    for(int j=19; j < bins.getNBins("pt"); ++j) bins.setLink(TString("eta1_pt18_rho") + k, TString("eta1_pt") + j + "_rho" + k);
   
    for(int i = 3; i < bins.getNBins("eta"); ++i){
      bins.setLink(TString("eta") + i + "_pt0_rho" + k, TString("eta") + i + "_pt1_rho" + k);						// Merge lowest two pt bins in forward (also low statistics, especially in non-CHS mode)
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
  for(int k=0; k < bins.getNBins("rho"); ++k){
    for(int i=0; i < 5; ++i){
      for(int j=12; j < 21; ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {6./3.,1.5/3.,1.5/3.});
      for(int j=10; j < 12; ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {5.5/3.,1.75/3.,1.75/3.});
      for(int j=8 ; j < 10; ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {5./3.,2./3.,2./3.});
      for(int j=7 ; j < 8;  ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {4.5/3.,2.25/3.,2.25/3.});
      for(int j=6 ; j < 7;  ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {4./3.,2.5/3.,2.5/3.});
      for(int j=5 ; j < 6;  ++j) bins.setWeights(TString("eta") + i + "_pt" + j + "_rho" + k, {3.5/3.,2.75/3.,2.75/3.});
    }
  }
  return bins;
}

#ifndef treeLooper_h
#define treeLooper_h

#include <vector>
#include <map>
#include <iostream>
#include "TFile.h"
#include "TChain.h"
#include "TString.h"

/*
 * Class holding the tree variables, and taking care of the weights in case of binned QCD
 */
class treeLooper{
  private:
    TChain *qgMiniTuple;
    int eventNumber, maxEntries;
    bool usePtHatBins, useBTagging;
    std::vector<TString> ptHatBins;
    std::vector<float> ptHatMin, nEvents, xsec;
    void setWeight();

  public:
    treeLooper(TString, TString, TString);
    ~treeLooper(){ delete qgMiniTuple;};
    bool next();
    void setMaxEntries(int i){ maxEntries = i;};
 
    float rho, pt, eta, axis2, ptD, bTagValue, ptDoubleCone, ratioDoubleCone, weight, motherMass;
    int nEvent, mult, partonId, jetIdLevel, nGenJetsInCone, nJetsForGenParticle, nGenJetsForGenParticle, nPriVtxs, nPileUp, motherId;
    bool balanced, matchedJet, bTag;
    float additionalJets;
    std::vector<float> *closebyJetdR, *closebyJetPt;
};


/*
 * Initialise tree with sample, jet type and (optional) location of tuples
 */
treeLooper::treeLooper(TString file, TString jetType, TString qgMiniTuplesDir = "~tomc/public/merged/QGMiniTupleSpring15/"){
  useBTagging = jetType.Contains("antib");
  jetType.ReplaceAll("_antib","");

  qgMiniTuple = new TChain("qgMiniTuple"+jetType+"/qgMiniTuple");

  usePtHatBins = (file == "QCD_AllPtBins");
  if(usePtHatBins){
    ptHatBins = {"15to30",    "30to50",   "50to80",  "80to120", "120to170", "170to300", "300to470", "470to600", "600to800", "800to1000", "1000to1400", "1400to1800", "1800to2400", "2400to3200", "3200"};
    ptHatMin  = { 15,          30,         50,        80,        120,        170,        300,        470,        600,        800,         1000,         1400,         1800,         2400,         3200};
    nEvents   = { 4942232,     4957245,    4978425,   3424782,   3452896,    3364368,    2933611,    1936832,    1964128,    1937216,     1487136,      197959,       193608,       194456,       192944};
    xsec      = { 1837410000., 140932000., 19204300., 2762530.,  471100.,    117276.,    7823.,      648.,       186.9 ,     32.293 ,     9.4183,       0.84265,      0.114943,     0.00682981,   0.000165445};
    for(TString ptHatBin : ptHatBins) 	qgMiniTuple->Add(qgMiniTuplesDir + "/qgMiniTuple_QCD_Pt-"+ ptHatBin +"_TuneCUETP8M1_13TeV-pythia8_asympt25ns.root", -1);
  } else if(file == "test"){	 	qgMiniTuple->Add("../test/qgMiniTuple.root", -1);
  } else 				qgMiniTuple->Add(qgMiniTuplesDir + "qgMiniTuple_" + file + ".root", -1);

  qgMiniTuple->SetBranchAddress("nEvent", 			&nEvent);
  qgMiniTuple->SetBranchAddress("rho", 				&rho);
  qgMiniTuple->SetBranchAddress("pt", 				&pt);
  qgMiniTuple->SetBranchAddress("eta", 				&eta);
  qgMiniTuple->SetBranchAddress("axis2", 			&axis2);
  qgMiniTuple->SetBranchAddress("ptD", 				&ptD);
  qgMiniTuple->SetBranchAddress("mult",	 			&mult);
  qgMiniTuple->SetBranchAddress("partonId", 			&partonId);
  qgMiniTuple->SetBranchAddress("motherId", 			&motherId);
  qgMiniTuple->SetBranchAddress("motherMass", 			&motherMass);
  qgMiniTuple->SetBranchAddress("bTag", 			&bTagValue);
  qgMiniTuple->SetBranchAddress("jetIdLevel",			&jetIdLevel);
  qgMiniTuple->SetBranchAddress("balanced",			&balanced);
  qgMiniTuple->SetBranchAddress("matchedJet",			&matchedJet);
  qgMiniTuple->SetBranchAddress("ptDoubleCone", 		&ptDoubleCone);
  qgMiniTuple->SetBranchAddress("nGenJetsInCone",		&nGenJetsInCone);
  qgMiniTuple->SetBranchAddress("nGenJetsForGenParticle",	&nGenJetsForGenParticle);
  qgMiniTuple->SetBranchAddress("nJetsForGenParticle",		&nJetsForGenParticle);

  closebyJetdR = nullptr;
  closebyJetPt = nullptr;
  qgMiniTuple->SetBranchAddress("closebyJetdR",			&closebyJetdR);
  qgMiniTuple->SetBranchAddress("closebyJetPt",			&closebyJetPt);

  eventNumber = 0;
  maxEntries = qgMiniTuple->GetEntries();
}


/*
 * Get next entry (and return true), set back to start at end (and return false)
 */
bool treeLooper::next(){
  if(eventNumber < maxEntries){
    qgMiniTuple->GetEntry(eventNumber++);
    eta = fabs(eta);
    bTag = (useBTagging && bTagValue > 0.244);
    additionalJets = 0;
    for(int i = 0; i < closebyJetdR->size(); ++i){
      if(closebyJetPt->at(i) < 20) continue;
      if(closebyJetdR->at(i) < 0.8) ++additionalJets;
    }
    ratioDoubleCone = ptDoubleCone/pt;
    setWeight();
    return true;
  } else {
    eventNumber = 0;
    return false;
  }
}


/*
 * Set weight in case of binned QCD
 */
void treeLooper::setWeight(){
  if(usePtHatBins){
    int ptIndex = 0;
    while(pt > ptHatMin[ptIndex]) ++ptIndex;
    int treeIndex = qgMiniTuple->GetTreeNumber();
    if(pt < ptHatMin[treeIndex])    weight = xsec[treeIndex]/nEvents[treeIndex];						// Apply weight to avoid too high contribution of jets with pt < pthat, those are radiated jets which could be wrongly matched
    else	                    weight = xsec[ptIndex]/nEvents[ptIndex];							// Avoid high weights from jets with pt > ptHat, those could destroy the smoothness of the histrograms
  } else			    weight = 1.;
}
#endif

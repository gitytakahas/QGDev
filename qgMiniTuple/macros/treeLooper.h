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
treeLooper::treeLooper(TString file, TString jetType, TString qgMiniTuplesDir = "~tomc/public/merged/QGMiniTuple4/"){
  useBTagging = jetType.Contains("antib");
  jetType.ReplaceAll("_antib","");

  if(file == "TTJets" || file.Contains("VBFHbb")) qgMiniTuple = new TChain("qgMiniTuple/qgMiniTuple");
  else qgMiniTuple = new TChain("qgMiniTuple"+jetType+"/qgMiniTuple");

  usePtHatBins = (file == "QCD_AllPtBins");
  if(usePtHatBins){
    ptHatBins = {"15to30","30to50","50to80","80to120","120to170","170to300","300to470","470to600","600to800","800to1000","1000to1400","1400to1800","1800to2400", "2400to3200","3200"};
    ptHatMin  = { 15,      30,      50,      80,       120,       170,       300,       470,       600,       800,        1000,        1400,        1800,         2400,        3200};
    nEvents   = { 2498841, 2449363, 2500315, 2500098,  2491398,   1490834,   1498032,   1498417,   1465278,   1500369,    1500642,     1500040,     2953210 ,     2958105,     2953431};
    xsec      = { 2237e6,  1615e5,  2211e4,  3000114,  493200,    120300,    7475,      587.1,     167,       28.25,      8.195,       0.7346,      0.102,        0.00644,     0.000163};
    for(TString ptHatBin : ptHatBins) 	qgMiniTuple->Add(qgMiniTuplesDir + "/qgMiniTuple_QCD_Pt-"+ ptHatBin +"_Tune4C_13TeV_pythia8.root", -1);
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
    if(pt < 10*ptHatMin[treeIndex]) weight = xsec[treeIndex]/nEvents[treeIndex];
    else	                    weight = xsec[ptIndex]/nEvents[ptIndex];							// Avoid high weights from jets with pt >>> ptHat
  } else			    weight = 1.;
}
#endif

#ifndef treeLooper_h
#define treeLooper_h

#include <vector>
#include <map>
#include <iostream>
#include "TFile.h"
#include "TChain.h"
#include "TString.h"

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
    ~treeLooper();
    bool next();
    void setMaxEntries(int i){ maxEntries = i;};
 
    float rho, pt, eta, axis2, axis2_dR2, axis2_dR3, ptD, ptD_dR2, ptD_dR3, bTagValue, weight; 
    int nEvent, nChg, nChg_dR2, nChg_dR3, mult, mult_dR2, mult_dR3, partonId, jetIdLevel, nGenJetsInCone, nJetsForGenParticle, nGenJetsForGenParticle, nPriVtxs, nPileUp; 
    bool balanced, matchedJet, bTag;
    float closestJetdR, closestJetPt, closebyJetsInCone;
    std::vector<float> *closebyJetdR, *closebyJetPt;
    std::vector<int> *closebyJetGenJetsInCone;
};


treeLooper::treeLooper(TString file, TString jetType, TString qgMiniTuplesDir = "~tomc/public/merged/QGMiniTuple3/"){
  useBTagging = jetType.Contains("antib");
  jetType.ReplaceAll("_antib","");
  qgMiniTuple = new TChain("qgMiniTuple"+jetType+"/qgMiniTuple");

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
  qgMiniTuple->SetBranchAddress("nPriVtxs", 			&nPriVtxs);
  qgMiniTuple->SetBranchAddress("nPileUp", 			&nPileUp);
  qgMiniTuple->SetBranchAddress("rho", 				&rho);
  qgMiniTuple->SetBranchAddress("pt", 				&pt);
  qgMiniTuple->SetBranchAddress("eta", 				&eta);
  qgMiniTuple->SetBranchAddress("axis2", 			&axis2);
  qgMiniTuple->SetBranchAddress("axis2_dR2", 			&axis2_dR2);
  qgMiniTuple->SetBranchAddress("axis2_dR3", 			&axis2_dR3);
  qgMiniTuple->SetBranchAddress("ptD", 				&ptD);
  qgMiniTuple->SetBranchAddress("ptD_dR2", 			&ptD_dR2);
  qgMiniTuple->SetBranchAddress("ptD_dR3", 			&ptD_dR3);
  qgMiniTuple->SetBranchAddress("mult",	 			&mult);
  qgMiniTuple->SetBranchAddress("mult_dR2",	 		&mult_dR2);
  qgMiniTuple->SetBranchAddress("mult_dR3",	 		&mult_dR3);
  qgMiniTuple->SetBranchAddress("nChg",	 			&nChg);
  qgMiniTuple->SetBranchAddress("nChg_dR2",	 		&nChg_dR2);
  qgMiniTuple->SetBranchAddress("nChg_dR3",	 		&nChg_dR3);
  qgMiniTuple->SetBranchAddress("partonId", 			&partonId);
  qgMiniTuple->SetBranchAddress("bTag", 			&bTagValue);
  qgMiniTuple->SetBranchAddress("jetIdLevel",			&jetIdLevel);
  qgMiniTuple->SetBranchAddress("balanced",			&balanced);
  qgMiniTuple->SetBranchAddress("matchedJet",			&matchedJet);
  qgMiniTuple->SetBranchAddress("nGenJetsInCone",		&nGenJetsInCone);
  qgMiniTuple->SetBranchAddress("nGenJetsForGenParticle",	&nGenJetsForGenParticle);
  qgMiniTuple->SetBranchAddress("nJetsForGenParticle",		&nJetsForGenParticle);

  closebyJetdR = nullptr;
  closebyJetPt = nullptr;
  closebyJetGenJetsInCone = nullptr;
  qgMiniTuple->SetBranchAddress("closebyJetdR",			&closebyJetdR);
  qgMiniTuple->SetBranchAddress("closebyJetPt",			&closebyJetPt);
  qgMiniTuple->SetBranchAddress("closebyJetGenJetsInCone",	&closebyJetGenJetsInCone);

  eventNumber = 0;
  maxEntries = qgMiniTuple->GetEntries();
}

treeLooper::~treeLooper(){
  delete qgMiniTuple;
}

bool treeLooper::next(){
  if(eventNumber < maxEntries){
    qgMiniTuple->GetEntry(eventNumber++);
    eta = fabs(eta);
    bTag = (useBTagging && bTagValue > 0.244);
    closebyJetsInCone = 0;
    closestJetdR = 99999;
    for(int i = 0; i < closebyJetdR->size(); ++i){
      if(closebyJetdR->at(i) < 0.8) ++closebyJetsInCone;
      if(closebyJetdR->at(i) < closestJetdR) closestJetdR = closebyJetdR->at(i);
    }
    setWeight();
    return true;
  } else {
    eventNumber = 0;
    return false;
  }
}

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

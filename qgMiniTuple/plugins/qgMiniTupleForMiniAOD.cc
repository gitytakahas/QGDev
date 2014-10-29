/*
  Package:    		QGDev/qgMiniTupleForMiniAOD
  Class:     		qgMiniTupleForMiniAOD
  Original Author:  	Tom Cornelis
 
  Description: 		Create small ntuple of QG-Likelihood variables and binning variables

*/

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "localQGLikelihoodCalculator.h"

#include "TFile.h"
#include "TTree.h"


class qgMiniTupleForMiniAOD : public edm::EDAnalyzer{
   public:
      explicit qgMiniTupleForMiniAOD(const edm::ParameterSet&);
      ~qgMiniTupleForMiniAOD(){};
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override {};
      template <class jetClass> void calcVariables(const jetClass *jet, float& axis2_, float& ptD_, int& mult_);

      edm::EDGetTokenT<double> rhoToken;
      edm::EDGetTokenT<pat::JetCollection> jetsToken;
      edm::EDGetTokenT<reco::GenJetCollection> genJetsToken;
      edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken;
      const double minJetPt;
      const double deltaRcut;
      const bool pythia6;

      edm::Service<TFileService> fs;
      TTree *tree;

      float rho, pt, eta, axis2, ptD, bTag, qg, closestJetdR;
      int nRun, nLumi, nEvent, mult, partonId, partonFlavour, nGenJetsInCone, nGenJetsForGenParticle, nJetsForGenParticle, nOtherJetsInCone;
      bool matchedJet, balanced;

      QGLikelihoodCalculator *qglcalc;
};


qgMiniTupleForMiniAOD::qgMiniTupleForMiniAOD(const edm::ParameterSet& iConfig) :
  rhoToken( 		consumes<double>(			iConfig.getParameter<edm::InputTag>("rhoInputTag"))),
  jetsToken(    	consumes<pat::JetCollection>(		iConfig.getParameter<edm::InputTag>("jetsInputTag"))),
  genJetsToken(    	consumes<reco::GenJetCollection>(	iConfig.getParameter<edm::InputTag>("genJetsInputTag"))),
  genParticlesToken(    consumes<reco::GenParticleCollection>(	iConfig.getParameter<edm::InputTag>("genParticlesInputTag"))),
  minJetPt(							iConfig.getUntrackedParameter<double>("minJetPt", 10.)),
  deltaRcut(							iConfig.getUntrackedParameter<double>("deltaRcut", 0.3)),
  pythia6(							iConfig.getUntrackedParameter<bool>("pythia6", false))
{
  qglcalc = new QGLikelihoodCalculator("/user/tomc/QGTagger/CMSSW_7_0_9_patch1/src/QGDev/qgMiniTuple/data/pdfQG_AK4chs_antib_13TeV.root");
}


void qgMiniTupleForMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  nRun 		= (int) iEvent.id().run();
  nLumi 	= (int) iEvent.id().luminosityBlock();
  nEvent	= (int) iEvent.id().event();

  edm::Handle<double> rhoHandle;					iEvent.getByToken(rhoToken, 		rhoHandle);
  edm::Handle<pat::JetCollection> jets;					iEvent.getByToken(jetsToken, 		jets);
  edm::Handle<reco::GenJetCollection> genJets;				iEvent.getByToken(genJetsToken, 	genJets);
  edm::Handle<reco::GenParticleCollection> genParticles;		iEvent.getByToken(genParticlesToken, 	genParticles);

  rho = (float) *rhoHandle;

  if(genJets->size() > 2){
    auto jet1 = genJets->begin();
    auto jet2 = genJets->begin() + 1;
    auto jet3 = genJets->begin() + 2;
    balanced = (jet3->pt() < 0.15*(jet1->pt()+jet2->pt()));
  } else balanced = true;

  for(auto jet = jets->begin();  jet != jets->end(); ++jet){
    if(jet == jets->begin() + 2) balanced = false;
    if(jet->pt() < minJetPt) continue;

    nGenJetsInCone = 0;
    for(auto genJet = genJets->begin(); genJet != genJets->end(); ++genJet){
      if(reco::deltaR(*jet, *genJet) < deltaRcut) ++nGenJetsInCone;
    }

    nOtherJetsInCone = 0;
    closestJetdR = 999;
    for(auto otherJet = jets->begin(); otherJet != jets->end(); ++otherJet){
      if(otherJet == jet) continue;
      if(reco::deltaR(*jet, *otherJet) < 0.8) ++nOtherJetsInCone;
      if(reco::deltaR(*jet, *otherJet) < closestJetdR) closestJetdR = reco::deltaR(*jet, *otherJet);
    }


    partonId = 0; matchedJet = false;
    float deltaRmin = 999;
    auto matchedGenParticle = genParticles->end();
    for(auto genParticle = genParticles->begin(); genParticle != genParticles->end(); ++genParticle){
      if(genParticle->status() != (pythia6? 3 : 23)) continue; 							//status 3 (pythia6) / status 23 (pythia8) for outgoing particles from the hardest subprocess
      if(abs(genParticle->pdgId()) > 5 && abs(genParticle->pdgId()) != 21) continue;				//Only keep quarks and gluons
      float thisDeltaR = reco::deltaR(genParticle->eta(), genParticle->phi(), jet->eta(), jet->phi());
      if(thisDeltaR < deltaRmin){
        deltaRmin = thisDeltaR;
        matchedGenParticle = genParticle;
      }
    }
    if(deltaRmin < deltaRcut){
      partonId 		= matchedGenParticle->pdgId();
      matchedJet	= true;

      nJetsForGenParticle = 0;
      for(auto otherJet = jets->begin(); otherJet != jets->end(); ++otherJet){
        if(reco::deltaR(*matchedGenParticle, *otherJet) < deltaRcut) ++nJetsForGenParticle;
      }
      nGenJetsForGenParticle = 0;
      for(auto genJet = genJets->begin(); genJet != genJets->end(); ++genJet){
        if(reco::deltaR(*matchedGenParticle, *genJet) < deltaRcut) ++nGenJetsForGenParticle;
      }
    }

    calcVariables(&*jet, axis2, ptD, mult);
    axis2 		= -std::log(axis2);
    partonFlavour	= jet->partonFlavour();
    pt			= jet->pt();
    eta			= jet->eta();
    bTag		= jet->bDiscriminator("combinedSecondaryVertexBJetTags");

    qg 			= qglcalc->computeQGLikelihood(pt, eta, rho, {(float) mult, ptD, axis2});

    tree->Fill();
  }
}


/// Calculation of axis2_, mult_ and ptD_
template <class jetClass> void qgMiniTupleForMiniAOD::calcVariables(const jetClass *jet, float& axis2_, float& ptD_, int& mult_){
  float sum_weight = 0., sum_deta = 0., sum_dphi = 0., sum_deta2 = 0., sum_dphi2 = 0., sum_detadphi = 0., sum_pt = 0.;
  int nChg_QC = 0, nNeutral_ptCut = 0;

  //Loop over the jet constituents
  for(int i = 0; i < jet->numberOfDaughters(); ++i){
    auto part = dynamic_cast<const pat::PackedCandidate*> (jet->daughter(i));
    if(part->charge()){
      if(part->fromPV() > 1 && part->trackHighPurity()) nChg_QC++;
      else continue;
    } else {
      if(part->pt() > 1.0) nNeutral_ptCut++;
      else continue;
    }
	  
    float deta = part->eta() - jet->eta();
    float dphi = reco::deltaPhi(part->phi(), jet->phi());
    float partPt = part->pt(); 
    float weight = partPt*partPt;

    sum_weight += weight;
    sum_pt += partPt;
    sum_deta += deta*weight;                  
    sum_dphi += dphi*weight;                                                                                             
    sum_deta2 += deta*deta*weight;                    
    sum_detadphi += deta*dphi*weight;                               
    sum_dphi2 += dphi*dphi*weight;
  }

  //Calculate axis2_ and ptD_
  float a = 0., b = 0., c = 0.;
  float ave_deta = 0., ave_dphi = 0., ave_deta2 = 0., ave_dphi2 = 0.;
  if(sum_weight > 0){
    ptD_ = sqrt(sum_weight)/sum_pt;
    ave_deta = sum_deta/sum_weight;
    ave_dphi = sum_dphi/sum_weight;
    ave_deta2 = sum_deta2/sum_weight;
    ave_dphi2 = sum_dphi2/sum_weight;
    a = ave_deta2 - ave_deta*ave_deta;                          
    b = ave_dphi2 - ave_dphi*ave_dphi;                          
    c = -(sum_detadphi/sum_weight - ave_deta*ave_dphi);                
  } else ptD_ = 0;
  float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
  if(a+b-delta > 0) axis2_ = sqrt(0.5*(a+b-delta));
  else axis2_ = 0.;

  mult_ = (nChg_QC + nNeutral_ptCut);
}


void qgMiniTupleForMiniAOD::beginJob(){
  tree = fs->make<TTree>("qgMiniTupleForMiniAOD","qgMiniTuple");
  tree->Branch("nRun" ,			&nRun, 			"nRun/I");
  tree->Branch("nLumi" ,		&nLumi, 		"nLumi/I");
  tree->Branch("nEvent" ,		&nEvent, 		"nEvent/I");
  tree->Branch("rho" ,			&rho, 			"rho/F");
  tree->Branch("pt" ,			&pt,			"pt/F");
  tree->Branch("eta",			&eta,			"eta/F");
  tree->Branch("axis2",			&axis2,			"axis2/F");
  tree->Branch("ptD",			&ptD,			"ptD/F");
  tree->Branch("mult",			&mult,			"mult/I");
  tree->Branch("qg",			&qg,			"qg/F");
  tree->Branch("bTag",			&bTag,			"bTag/F");
  tree->Branch("partonId",		&partonId,		"partonId/I");
  tree->Branch("partonFlavour",		&partonFlavour,		"partonFlavour/I");
  tree->Branch("nGenJetsInCone",	&nGenJetsInCone,	"nGenJetsInCone/I");
  tree->Branch("closestJetdR",		&closestJetdR,		"closestJetdR/F");
  tree->Branch("nOtherJetsInCone",	&nOtherJetsInCone,	"nOtherJetsInCone/I");
  tree->Branch("matchedJet",		&matchedJet,		"matchedJet/O");
  tree->Branch("balanced",		&balanced,		"balanced/O");
  tree->Branch("nGenJetsForGenParticle",&nGenJetsForGenParticle,"nGenJetsForGenParticle/I");
  tree->Branch("nJetsForGenParticle",   &nJetsForGenParticle,   "nJetsForGenParticle/I");
}


void qgMiniTupleForMiniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
  edm::ParameterSetDescription desc;
  desc.addUntracked<std::string>("fileName","qgMiniTupleForMiniAOD.root");
  desc.add<edm::InputTag>("rhoInputTag");
  desc.add<edm::InputTag>("jetsInputTag");
  desc.add<edm::InputTag>("genJetsInputTag");
  desc.add<edm::InputTag>("genParticlesInputTag");
  desc.addUntracked<double>("minJetPt", 10.);
  desc.addUntracked<double>("deltaRcut", 0.3);
  desc.addUntracked<bool>("pythia6", false);
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(qgMiniTupleForMiniAOD);

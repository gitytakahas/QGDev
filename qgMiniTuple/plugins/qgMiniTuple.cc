/*
  Package:    		QGDev/qgMiniTuple
  Class:     		qgMiniTuple
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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TFile.h"
#include "TTree.h"


class qgMiniTuple : public edm::EDAnalyzer{
   public:
      explicit qgMiniTuple(const edm::ParameterSet&);
      ~qgMiniTuple(){};
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      bool jetId(const reco::PFJet *jet, bool tight = false, bool medium = false);
      virtual void endJob() override {};

      edm::EDGetTokenT<double> rhoToken;
      edm::EDGetTokenT<reco::PFJetCollection> jetsToken;
      edm::InputTag qgVariablesInputTag;
      edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken;
      edm::EDGetTokenT<edm::ValueMap<float>> qgToken, axis2Token, ptDToken;
      edm::EDGetTokenT<edm::ValueMap<int>> multToken;
      edm::EDGetTokenT<reco::JetTagCollection> bTagToken;
      const double minJetPt, deltaRcut;
      const bool pythia6, bTagUsed;

      edm::Service<TFileService> fs;
      TTree *tree;

      bool jetIdLoose, jetIdMedium, jetIdTight;
      float rho, qg, pt, eta, axis2, ptD, deltaR, bTag;
      int nRun, nLumi, nEvent, mult, partonId;
};


qgMiniTuple::qgMiniTuple(const edm::ParameterSet& iConfig) :
  rhoToken( 		consumes<double>(			iConfig.getParameter<edm::InputTag>("rhoInputTag"))),
  jetsToken(    	consumes<reco::PFJetCollection>(	iConfig.getParameter<edm::InputTag>("jetsInputTag"))),
  qgVariablesInputTag(  					iConfig.getParameter<edm::InputTag>("qgVariablesInputTag")),
  genParticlesToken(    consumes<reco::GenParticleCollection>(	iConfig.getParameter<edm::InputTag>("genParticlesInputTag"))),
  minJetPt(							iConfig.getUntrackedParameter<double>("minJetPt", 10.)),
  deltaRcut(							iConfig.getUntrackedParameter<double>("deltaRcut", 0.3)),
  pythia6(							iConfig.getUntrackedParameter<bool>("pythia6", false)),
  bTagUsed(							iConfig.getUntrackedParameter<bool>("bTagUsed", false))
{
  qgToken	= 	consumes<edm::ValueMap<float>>(		edm::InputTag(qgVariablesInputTag.label(), "qgLikelihood"));
  axis2Token	= 	consumes<edm::ValueMap<float>>(		edm::InputTag(qgVariablesInputTag.label(), "axis2Likelihood"));
  multToken	= 	consumes<edm::ValueMap<int>>(		edm::InputTag(qgVariablesInputTag.label(), "multLikelihood"));
  ptDToken	= 	consumes<edm::ValueMap<float>>(		edm::InputTag(qgVariablesInputTag.label(), "ptDLikelihood"));
  bTagToken	=	consumes<reco::JetTagCollection>(	edm::InputTag("combinedSecondaryVertexV2BJetTags"));
}


void qgMiniTuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  nRun 		= (int) iEvent.id().run();
  nLumi 	= (int) iEvent.id().luminosityBlock();
  nEvent	= (int) iEvent.id().event();

  edm::Handle<double> rho_;
  iEvent.getByToken(rhoToken, rho_);
  rho = (float) *rho_;

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetsToken, jets);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesToken, genParticles);

  edm::Handle<reco::JetTagCollection> bTagHandle;
  iEvent.getByToken(bTagToken, bTagHandle);

  edm::Handle<edm::ValueMap<float>> qgHandle, axis2Handle, ptDHandle;
  edm::Handle<edm::ValueMap<int>> multHandle;
  iEvent.getByToken(qgToken, qgHandle);
  iEvent.getByToken(axis2Token, axis2Handle);
  iEvent.getByToken(multToken, multHandle);
  iEvent.getByToken(ptDToken, ptDHandle);

  for(auto jet = jets->begin();  jet != jets->end(); ++jet){
    if(jet->pt() < minJetPt) continue;
    edm::RefToBase<reco::Jet> jetRef(edm::Ref<reco::PFJetCollection>(jets, (jet - jets->begin())));

    float deltaRmin = 999;
    auto matchedGenParticle = genParticles->end();
    for(auto genParticle = genParticles->begin(); genParticle != genParticles->end(); ++genParticle){
      if(pythia6 && genParticle->status() != 3 ) continue; 							//Pythia6: status 3  for outgoing particles from the hardest subprocess
      if(!pythia6 && genParticle->status() != 23 ) continue; 							//Pythia8: status 23 for outgoing particles from the hardest subprocess
      if(abs(genParticle->pdgId()) > 5 && abs(genParticle->pdgId()) !=21) continue;				//Only keep quarks and gluons
      float thisDeltaR = reco::deltaR(genParticle->eta(), genParticle->phi(), jet->eta(), jet->phi());
      if(thisDeltaR < deltaRmin){
        deltaRmin = thisDeltaR;
        matchedGenParticle = genParticle;
      }
    }

    if(deltaRmin > deltaRcut) continue;
    partonId	= matchedGenParticle->pdgId();
    deltaR	= deltaRmin;

    pt		= jet->pt();
    eta		= jet->eta();
    qg		= (*qgHandle)[jetRef];
    axis2	= (*axis2Handle)[jetRef];
    mult	= (*multHandle)[jetRef];
    ptD		= (*ptDHandle)[jetRef];
    bTag	= (bTagUsed ? (*bTagHandle)[jetRef] : -1);

    jetIdLoose	= jetId(&*jet); 
    jetIdMedium	= jetId(&*jet, false, true); 
    jetIdTight  = jetId(&*jet, true); 
    tree->Fill();
  }
}

void qgMiniTuple::beginJob(){
  tree = fs->make<TTree>("qgMiniTuple","qgMiniTuple");
  tree->Branch("nRun" ,		&nRun, 		"nRun/I");
  tree->Branch("nLumi" ,	&nLumi, 	"nLumi/I");
  tree->Branch("nEvent" ,	&nEvent, 	"nEvent/I");
  tree->Branch("rho" ,		&rho, 		"rho/F");
  tree->Branch("pt" ,		&pt,		"pt/F");
  tree->Branch("eta",		&eta,		"eta/F");
  tree->Branch("qg",		&qg,		"qg/F");
  tree->Branch("axis2",		&axis2,		"axis2/F");
  tree->Branch("ptD",		&ptD,		"ptD/F");
  tree->Branch("mult",		&mult,		"mult/I");
  tree->Branch("bTag",		&bTag,		"bTag/F");
  tree->Branch("partonId",	&partonId,	"partonId/I");
  tree->Branch("deltaR",	&deltaR,	"deltaR/F");
  tree->Branch("jetIdLoose",	&jetIdLoose,	"jetIdLoose/O");
  tree->Branch("jetIdMedium",	&jetIdMedium,	"jetIdMedium/O");
  tree->Branch("jetIdTight",	&jetIdTight,	"jetIdTight/O");
}


bool qgMiniTuple::jetId(const reco::PFJet *jet, bool tight, bool medium){
  float jetEnergyUncorrected 		= jet->chargedHadronEnergy() + jet->neutralHadronEnergy() + jet->photonEnergy() +
  					  jet->electronEnergy() + jet->muonEnergy() + jet->HFHadronEnergy() + jet->HFEMEnergy();
  float neutralHadronEnergyFraction 	= (jet->neutralHadronEnergy() + jet->HFHadronEnergy())/jetEnergyUncorrected;
  float neutralEmEnergyFraction 	= (jet->neutralEmEnergy())/jetEnergyUncorrected;
  float chargedHadronEnergyFraction 	= (jet->chargedHadronEnergy())/jetEnergyUncorrected;
  float chargedEmEnergyFraction 	= (jet->chargedEmEnergy())/jetEnergyUncorrected;

  if(!(neutralHadronEnergyFraction 	< (tight ? 0.90 : (medium ? 0.95 : .99)))) 	return false;
  if(!(neutralEmEnergyFraction 		< (tight ? 0.90 : (medium ? 0.95 : .99)))) 	return false;
  if(!((jet->chargedMultiplicity() + jet->neutralMultiplicity()) > 1)) 			return false;
  if(fabs(jet->eta()) < 2.4){
    if(!(chargedHadronEnergyFraction > 0)) 						return false;
    if(!(chargedEmEnergyFraction < .99)) 						return false;
    if(!(jet->chargedMultiplicity() > 0)) 						return false;
  }
  return true;
}


void qgMiniTuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
  edm::ParameterSetDescription desc;
  desc.addUntracked<std::string>("fileName","qgMiniTuple.root");
  desc.add<edm::InputTag>("rhoInputTag");
  desc.add<edm::InputTag>("jetsInputTag");
  desc.add<edm::InputTag>("qgVariablesInputTag");
  desc.add<edm::InputTag>("genParticlesInputTag");
  desc.addUntracked<double>("minJetPt", 10.);
  desc.addUntracked<double>("deltaRcut", 0.3);
  desc.addUntracked<bool>("pythia6", false);
  desc.addUntracked<bool>("bTagUsed", false);
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(qgMiniTuple);

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
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "localQGLikelihoodCalculator.h"
#include "TFile.h"
#include "TTree.h"


/*
 * qgMiniTuple class
 */
class qgMiniTuple : public edm::EDAnalyzer{
   public:
      explicit qgMiniTuple(const edm::ParameterSet&);
      ~qgMiniTuple(){};

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      template <class jetClass> bool jetId(const jetClass *jet, bool tight = false, bool medium = false);
      template <class jetClass> void calcVariables(const jetClass *jet, float& axis2_, float& ptD_, int& mult_, int& nChg_, edm::Handle<reco::VertexCollection> vC);
      template <class jetCollection, class candidateCollection> void analyzeEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, const jetCollection& jets, const jetCollection& jetsAK8, const candidateCollection& pfCandidates);
      template <class jetCollection> edm::RefToBase<reco::Jet> jetRef(const jetCollection& jets, typename jetCollection::element_type::const_iterator& jet);
      virtual void endJob() override;

      edm::EDGetTokenT<double> rhoToken;
      edm::EDGetTokenT<reco::VertexCollection> vertexToken;
      edm::EDGetTokenT<reco::GenJetCollection> genJetsToken;
      edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken;
      edm::EDGetTokenT<reco::JetTagCollection> bTagToken;
      edm::EDGetTokenT<reco::PFJetCollection> jetsToken, jetsTokenAK8;
      edm::EDGetTokenT<pat::JetCollection> patJetsToken, patJetsTokenAK8;
      edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidatesToken;
      edm::EDGetTokenT<pat::PackedCandidateCollection> patCandidatesToken;
      edm::InputTag jetsInputTag, jetsInputTagAK8, pfCandidatesInputTag;
      edm::InputTag csvInputTag;
      std::string jecService;
      std::string jecServiceAK8;
/*    edm::InputTag qgVariablesInputTag;
      edm::EDGetTokenT<edm::ValueMap<float>> qgToken, axis2Token, ptDToken;
      edm::EDGetTokenT<edm::ValueMap<int>> multToken;*/
      const double minJetPt;
      const double deltaRcut;
      const bool pythia6;
      const bool usePatJets;

      const JetCorrector *JEC, *JECAK8;
      edm::Service<TFileService> fs;
      TTree *tree;
//    QGLikelihoodCalculator *qglcalc;

      float rho, pt, eta, axis2, ptD, bTag, deltaRAK8, ptAK8, ptDoubleCone;
      int nRun, nLumi, nEvent, nPileUp, nPriVtxs, mult, nChg, partonId, jetIdLevel, nGenJetsInCone, nGenJetsForGenParticle, nJetsForGenParticle;
      bool matchedJet, balanced;
      std::vector<float> *closebyJetdR, *closebyJetPt;
      std::vector<int> *closebyJetGenJetsInCone;
};


/*
 * Constructor
 */
qgMiniTuple::qgMiniTuple(const edm::ParameterSet& iConfig) :
  rhoToken( 		consumes<double>(					iConfig.getParameter<edm::InputTag>("rhoInputTag"))),
  vertexToken(    	consumes<reco::VertexCollection>(			iConfig.getParameter<edm::InputTag>("vertexInputTag"))),
  genJetsToken(    	consumes<reco::GenJetCollection>(			iConfig.getParameter<edm::InputTag>("genJetsInputTag"))),
  genParticlesToken(    consumes<reco::GenParticleCollection>(			iConfig.getParameter<edm::InputTag>("genParticlesInputTag"))),
  jetsInputTag(    								iConfig.getParameter<edm::InputTag>("jetsInputTag")),
  jetsInputTagAK8( 	 							iConfig.getParameter<edm::InputTag>("jetsInputTagAK8")),
  pfCandidatesInputTag(								iConfig.getParameter<edm::InputTag>("pfCandidatesInputTag")),
  csvInputTag(    								iConfig.getParameter<edm::InputTag>("csvInputTag")),
  jecService( 									iConfig.getParameter<std::string>("jec")),
  jecServiceAK8(								iConfig.getParameter<std::string>("jecAK8")),
//qgVariablesInputTag(  							iConfig.getParameter<edm::InputTag>("qgVariablesInputTag")),
  minJetPt(									iConfig.getUntrackedParameter<double>("minJetPt", 20.)),
  deltaRcut(									iConfig.getUntrackedParameter<double>("deltaRcut", 0.3)),
  pythia6(									iConfig.getUntrackedParameter<bool>("pythia6", false)),
  usePatJets(									iConfig.getUntrackedParameter<bool>("usePatJets", false))
{
  jetsToken	=	consumes<reco::PFJetCollection>(	edm::InputTag(jetsInputTag));
  jetsTokenAK8	=	consumes<reco::PFJetCollection>(	edm::InputTag(jetsInputTagAK8));
  pfCandidatesToken =	consumes<reco::PFCandidateCollection>(  edm::InputTag(pfCandidatesInputTag));
  patJetsToken	=	consumes<pat::JetCollection>(		edm::InputTag(jetsInputTag));
  patJetsTokenAK8 =	consumes<pat::JetCollection>(		edm::InputTag(jetsInputTagAK8));
  patCandidatesToken =	consumes<pat::PackedCandidateCollection>(edm::InputTag(pfCandidatesInputTag));
  bTagToken	= 	consumes<reco::JetTagCollection>(       edm::InputTag(csvInputTag));
/*qgToken	= 	consumes<edm::ValueMap<float>>(		edm::InputTag(qgVariablesInputTag.label(), "qgLikelihood"));
  axis2Token	= 	consumes<edm::ValueMap<float>>(		edm::InputTag(qgVariablesInputTag.label(), "axis2Likelihood"));
  multToken	= 	consumes<edm::ValueMap<int>>(		edm::InputTag(qgVariablesInputTag.label(), "multLikelihood"));
  ptDToken	= 	consumes<edm::ValueMap<float>>(		edm::InputTag(qgVariablesInputTag.label(), "ptDLikelihood"));*/
//qglcalc 	= 	new QGLikelihoodCalculator("/user/tomc/QGTagger/CMSSW_7_0_9_patch1/src/QGDev/qgMiniTuple/data/pdfQG_AK4chs_antib_13TeV.root");
}


/*
 * Prepare for analyzing the event: choose pat or reco jets
 */
void qgMiniTuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  for(auto v : {closebyJetdR, closebyJetPt}) v->clear();
  for(auto v : {closebyJetGenJetsInCone})    v->clear();

  nRun 		= (int) iEvent.id().run();
  nLumi 	= (int) iEvent.id().luminosityBlock();
  nEvent	= (int) iEvent.id().event();

  if(usePatJets){
    edm::Handle<pat::JetCollection> jets;
    edm::Handle<pat::JetCollection> jetsAK8;
    edm::Handle<pat::PackedCandidateCollection> pfCandidates;
    iEvent.getByToken(patJetsToken, jets);
    iEvent.getByToken(patJetsTokenAK8, jetsAK8);
    iEvent.getByToken(patCandidatesToken, pfCandidates);
    analyzeEvent(iEvent, iSetup, jets, jetsAK8, pfCandidates);
  } else {
    JEC = JetCorrector::getJetCorrector(jecService, iSetup);
    JECAK8 = JetCorrector::getJetCorrector(jecServiceAK8, iSetup);
    edm::Handle<reco::PFJetCollection> jets;
    edm::Handle<reco::PFJetCollection> jetsAK8;
    edm::Handle<reco::PFCandidateCollection> pfCandidates;
    iEvent.getByToken(jetsToken, jets);
    iEvent.getByToken(jetsTokenAK8, jetsAK8);
    iEvent.getByToken(pfCandidatesToken, pfCandidates);
    analyzeEvent(iEvent, iSetup, jets, jetsAK8, pfCandidates);
  }
}


/*
 * Analyze event
 */
template <class jetCollection, class candidateCollection> void qgMiniTuple::analyzeEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, const jetCollection& jets, const jetCollection& jetsAK8, const candidateCollection& pfCandidates){
  edm::Handle<std::vector<PileupSummaryInfo>> PupInfo;					iEvent.getByLabel("addPileupInfo", 	PupInfo);
  edm::Handle<reco::VertexCollection> vertexCollection;					iEvent.getByToken(vertexToken, 		vertexCollection);
  edm::Handle<double> rhoHandle;							iEvent.getByToken(rhoToken, 		rhoHandle);
  edm::Handle<reco::GenJetCollection> genJets;						iEvent.getByToken(genJetsToken, 	genJets);
  edm::Handle<reco::GenParticleCollection> genParticles;				iEvent.getByToken(genParticlesToken, 	genParticles);
//  edm::Handle<reco::PFCandidateCollection> pfCandidates;				iEvent.getByToken(pfCandidatesToken, 	pfCandidates);
  edm::Handle<reco::JetTagCollection> bTagHandle;			if(!usePatJets) iEvent.getByToken(bTagToken, 		bTagHandle);
/*edm::Handle<edm::ValueMap<float>> qgHandle;				if(!usePatJets)	iEvent.getByToken(qgToken, 		qgHandle);
  edm::Handle<edm::ValueMap<float>> axis2Handle; 			if(!usePatJets)	iEvent.getByToken(axis2Token, 		axis2Handle);
  edm::Handle<edm::ValueMap<float>> ptDHandle;  			if(!usePatJets)	iEvent.getByToken(multToken, 		multHandle);
  edm::Handle<edm::ValueMap<int>> multHandle;   			if(!usePatJets)	iEvent.getByToken(ptDToken,		ptDHandle);
*/

  nPriVtxs 	= vertexCollection->size();
  rho 		= (float) *rhoHandle;
  nPileUp 	= -1;
  if(PupInfo.isValid()){
    auto PVI = PupInfo->begin();
    while(PVI->getBunchCrossing() != 0 && PVI != PupInfo->end()) ++PVI;
    if(PVI != PupInfo->end()) nPileUp = PVI->getPU_NumInteractions();
  }

 
  if(genJets->size() > 2){
    auto jet1 = genJets->begin();
    auto jet2 = genJets->begin() + 1;
    auto jet3 = genJets->begin() + 2;
    balanced = (jet3->pt() < 0.15*(jet1->pt()+jet2->pt()));
  } else balanced = true;


  for(auto jet = jets->begin();  jet != jets->end(); ++jet){
    if(jet == jets->begin() + 2) balanced = false;

    pt 			= jet->pt()*(usePatJets? 1. : JEC->correction(*jet, iEvent, iSetup));
    if(pt < minJetPt) continue;

    eta			= jet->eta();
    jetIdLevel		= jetId(&*jet) + jetId(&*jet, false, true) + jetId(&*jet, true); 

    nGenJetsInCone 	= 0;
    for(auto genJet = genJets->begin(); genJet != genJets->end(); ++genJet){
      if(reco::deltaR(*jet, *genJet) < deltaRcut) ++nGenJetsInCone;
    }

    // Closeby jet study variables
    for(auto otherJet = jets->begin(); otherJet != jets->end(); ++otherJet){
      if(otherJet == jet) continue;
      float dR = reco::deltaR(*jet, *otherJet);
      if(dR > 1.2) continue;
      float nGenJetsInConeOtherJet = 0;
      for(auto genJet = genJets->begin(); genJet != genJets->end(); ++genJet){
        if(reco::deltaR(*otherJet, *genJet) < deltaRcut) ++nGenJetsInConeOtherJet;
      }
      closebyJetdR->push_back(dR);
      closebyJetPt->push_back(otherJet->pt()*(usePatJets? 1. : JEC->correction(*jet, iEvent, iSetup)));
      closebyJetGenJetsInCone->push_back(nGenJetsInConeOtherJet);
    }

    // Parton Id matching
    partonId = 0; matchedJet = false;
    float deltaRmin = 999;
    auto matchedGenParticle = genParticles->end();
    for(auto genParticle = genParticles->begin(); genParticle != genParticles->end(); ++genParticle){
      if(genParticle->status() != (pythia6? 3 : 23)) continue; 							//status 3 (pythia6) / status 23 (pythia8) for outgoing particles from the hardest subprocess
      if(abs(genParticle->pdgId()) > 5 && abs(genParticle->pdgId()) != 21) continue;				//Only keep quarks and gluons
      float thisDeltaR = reco::deltaR(*genParticle, *jet);
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
    } else continue;												//Keep only matched jets for the moment

    if(usePatJets){
      auto patjet = dynamic_cast<const pat::Jet*> (&*jet);
      bTag		= patjet->bDiscriminator(csvInputTag.label());
    } else {
      bTag		= (*bTagHandle)[jetRef(jets, jet)];
/*    qg		= (*qgHandle)[jetRef(jets, jet)];
      axis2		= (*axis2Handle)[jetRef(jets, jet)];
      mult		= (*multHandle)[jetRef(jets, jet)];
      ptD		= (*ptDHandle)[jetRef(jets, jet)];*/
    }

    calcVariables(&*jet, axis2,     ptD,     mult,     nChg,     vertexCollection);
    axis2 		= -std::log(axis2);

//  qg 			= qglcalc->computeQGLikelihood(pt, eta, rho, {(float) mult, ptD, axis2});

    deltaRAK8 = 9999;
    ptAK8 = 0;
    for(auto jetAK8 = jetsAK8->begin(); jetAK8 != jetsAK8->end(); ++jetAK8){
      if(reco::deltaR(*jet, *jetAK8) < deltaRAK8){
        deltaRAK8 = reco::deltaR(*jet, *jetAK8);
        ptAK8 = jetAK8->pt()*(usePatJets? 1. : JECAK8->correction(*jet, iEvent, iSetup));
      }
    }

    ptDoubleCone = 0;
    for(auto pfCandidate = pfCandidates->begin(); pfCandidate != pfCandidates->end(); ++pfCandidate){
      if(reco::deltaR(*pfCandidate, *jet) < 0.8) ptDoubleCone += pfCandidate->pt();
    }
 
    tree->Fill();
  }
}


// Some dirty C++ stuff to get it compiled for both pat::Jet and reco::PFJet
template <class jetCollection> edm::RefToBase<reco::Jet> qgMiniTuple::jetRef(const jetCollection& jets, typename jetCollection::element_type::const_iterator& jet){
  return edm::RefToBase<reco::Jet>();
}
template <> edm::RefToBase<reco::Jet> qgMiniTuple::jetRef<edm::Handle<reco::PFJetCollection>>(const edm::Handle<reco::PFJetCollection>& jets, reco::PFJetCollection::const_iterator& jet){
  edm::RefToBase<reco::Jet> thisJetRef(edm::Ref<reco::PFJetCollection>(jets, (jet- jets->begin())));
  return thisJetRef;
}

/*
 * Begin job: create vectors and set up tree
 */
void qgMiniTuple::beginJob(){
  for(auto v : {&closebyJetdR, &closebyJetPt}) *v = new std::vector<float>();
  for(auto v : {&closebyJetGenJetsInCone})     *v = new std::vector<int>();

  tree = fs->make<TTree>("qgMiniTuple","qgMiniTuple");
  tree->Branch("nRun" ,			&nRun, 			"nRun/I");
  tree->Branch("nLumi" ,		&nLumi, 		"nLumi/I");
  tree->Branch("nEvent" ,		&nEvent, 		"nEvent/I");
  tree->Branch("nPileUp",		&nPileUp, 		"nPileUp/I");
  tree->Branch("nPriVtxs",		&nPriVtxs, 		"nPriVtxs/I");
  tree->Branch("rho" ,			&rho, 			"rho/F");
  tree->Branch("pt" ,			&pt,			"pt/F");
  tree->Branch("eta",			&eta,			"eta/F");
  tree->Branch("axis2",			&axis2,			"axis2/F");
  tree->Branch("ptD",			&ptD,			"ptD/F");
  tree->Branch("mult",			&mult,			"mult/I");
  tree->Branch("nChg",			&nChg,			"nChg/I");
  tree->Branch("bTag",			&bTag,			"bTag/F");
  tree->Branch("partonId",		&partonId,		"partonId/I");
  tree->Branch("jetIdLevel",		&jetIdLevel,		"jetIdLevel/I");
  tree->Branch("nGenJetsInCone",	&nGenJetsInCone,	"nGenJetsInCone/I");
  tree->Branch("matchedJet",		&matchedJet,		"matchedJet/O");
  tree->Branch("balanced",		&balanced,		"balanced/O");
  tree->Branch("deltaRAK8",		&deltaRAK8,		"deltaRAK8/F");
  tree->Branch("ptAK8",			&ptAK8,			"ptAK8/F");
  tree->Branch("ptDoubleCone",		&ptDoubleCone,		"ptDoubleCone/F");
  tree->Branch("nGenJetsForGenParticle",&nGenJetsForGenParticle,"nGenJetsForGenParticle/I");
  tree->Branch("nJetsForGenParticle",   &nJetsForGenParticle,   "nJetsForGenParticle/I");
  tree->Branch("closebyJetdR",			"vector<float>",	&closebyJetdR);
  tree->Branch("closebyJetPt",			"vector<float>",	&closebyJetPt);
  tree->Branch("closebyJetGenJetsInCone",	"vector<int>",		&closebyJetGenJetsInCone);
}


/*
 * End of jobs: delete vectors
 */
void qgMiniTuple::endJob(){
 for(auto v : {closebyJetdR, closebyJetPt}) delete v;
 for(auto v : {closebyJetGenJetsInCone})    delete v;
}


/*
 * Calculation for of axis2, ptD and mult (works both on reco and pat)
 */
template <class jetClass> void qgMiniTuple::calcVariables(const jetClass *jet, float& axis2_, float& ptD_, int& mult_, int& nChg_, edm::Handle<reco::VertexCollection> vC){
  auto vtxLead = vC->begin();

  float sum_weight = 0., sum_deta = 0., sum_dphi = 0., sum_deta2 = 0., sum_dphi2 = 0., sum_detadphi = 0., sum_pt = 0.;
  mult_ = 0; nChg_ = 0;

  //Loop over the jet constituents
  for(auto daughter = jet->begin(); daughter < jet->end(); ++daughter){
    if(usePatJets){
      auto part = dynamic_cast<const pat::PackedCandidate*> (&*daughter);
      if(!part) continue;
      if(part->charge()){
        if(!(part->fromPV() > 1 && part->trackHighPurity())) continue;
        else ++nChg_;
      } else if(part->pt() < 1.0) continue;
    } else {
      auto part = dynamic_cast<const reco::PFCandidate*> (&*daughter);
      if(!part) continue;
      reco::TrackRef itrk = part->trackRef();
      if(itrk.isNonnull()){;
        auto vtxClose = vC->begin();
        for(auto vtx = vC->begin(); vtx != vC->end(); ++vtx){
          if(fabs(itrk->dz(vtx->position())) < fabs(itrk->dz(vtxClose->position()))) vtxClose = vtx;
        }
        if(!(vtxClose == vtxLead && itrk->quality(reco::TrackBase::qualityByName("highPurity")))) continue;
        else ++nChg_;
      } else if(part->pt() < 1.0) continue;
    }

    float deta 	 = daughter->eta() - jet->eta();
    float dphi 	 = reco::deltaPhi(*daughter, *jet);
    float partPt = daughter->pt();
    float weight = partPt*partPt;

    ++mult_;

    sum_weight 	 += weight;
    sum_pt 	 += partPt;
    sum_deta     += deta*weight;
    sum_dphi 	 += dphi*weight;
    sum_deta2  	 += deta*deta*weight;
    sum_detadphi += deta*dphi*weight;
    sum_dphi2 	 += dphi*dphi*weight;
  }

  //Calculate axis2 and ptD
  float a = 0., b = 0., c = 0.;
  float ave_deta = 0., ave_dphi = 0., ave_deta2 = 0., ave_dphi2 = 0.;
  if(sum_weight > 0){
    ptD_ 	= sqrt(sum_weight)/sum_pt;
    ave_deta 	= sum_deta/sum_weight;
    ave_dphi 	= sum_dphi/sum_weight;
    ave_deta2 	= sum_deta2/sum_weight;
    ave_dphi2 	= sum_dphi2/sum_weight;
    a 		= ave_deta2 - ave_deta*ave_deta;
    b 		= ave_dphi2 - ave_dphi*ave_dphi;
    c 		= -(sum_detadphi/sum_weight - ave_deta*ave_dphi);
  } else ptD_ 	= 0;
  float delta 	= sqrt(fabs((a-b)*(a-b)+4*c*c));
  if(a+b-delta > 0) axis2_ = sqrt(0.5*(a+b-delta));
  else 		    axis2_ = 0.;
}


/*
 * Calculate jetId for levels loose, medium and tight
 */
template<class jetClass> bool qgMiniTuple::jetId(const jetClass *jet, bool tight, bool medium){
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


DEFINE_FWK_MODULE(qgMiniTuple);

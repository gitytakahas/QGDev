/*
  Package:    		QGDev/qgMiniTuple
  Class:     		qgMiniTuple
  Original Author:  	Tom Cornelis
 
  Description: 		Create small ntuple of QG-Likelihood variables and binning variables

*/

#include <memory>
#include <tuple>

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
      virtual void 			beginJob() override;
      virtual void 			endJob() override;
      virtual void 			analyze(const edm::Event&, const edm::EventSetup&) override;
      template <class jetClass> bool 	jetId(const jetClass *jet, bool tight = false, bool medium = false);
      std::tuple<int, float, float> 	calcVariables(const reco::Jet *jet, edm::Handle<reco::VertexCollection>& vC);
      bool 				isPatJetCollection(const edm::Handle<edm::View<reco::Jet>>& jets);
      bool 				isPackedCandidate(const reco::Candidate* candidate);
      template<class a, class b> int 	countInCone(a center, b objectsToCount);

      edm::EDGetTokenT<edm::View<reco::Jet>> 		jetsToken;
      edm::EDGetTokenT<edm::View<reco::Candidate>> 	candidatesToken;
      edm::EDGetTokenT<double> 				rhoToken;
      edm::EDGetTokenT<reco::VertexCollection> 		vertexToken;
      edm::EDGetTokenT<reco::GenJetCollection> 		genJetsToken;
      edm::EDGetTokenT<reco::GenParticleCollection> 	genParticlesToken;
      edm::EDGetTokenT<reco::JetTagCollection> 		bTagToken;
      edm::InputTag					csvInputTag;
      std::string 					jecService;
/*    edm::InputTag 					qgVariablesInputTag;
      edm::EDGetTokenT<edm::ValueMap<float>> 		qgToken, axis2Token, ptDToken;
      edm::EDGetTokenT<edm::ValueMap<int>> 		multToken;*/
      const double 					minJetPt, deltaRcut;
      const bool 					pythia6, useQC;

      const JetCorrector 				*JEC;
      edm::Service<TFileService> 			fs;
      TTree 						*tree;
//    QGLikelihoodCalculator 				*qglcalc;

      float rho, pt, eta, axis2, ptD, bTag, ptDoubleCone, motherMass;
      int nEvent, nPileUp, nPriVtxs, mult, partonId, jetIdLevel, nGenJetsInCone, nGenJetsForGenParticle, nJetsForGenParticle, motherId;
      bool matchedJet, balanced;
      std::vector<float> *closebyJetdR, *closebyJetPt;

      bool weStillNeedToCheckJets, weAreUsingPatJets;
      bool weStillNeedToCheckJetCandidates, weAreUsingPackedCandidates;
};


/*
 * Constructor
 */
qgMiniTuple::qgMiniTuple(const edm::ParameterSet& iConfig) :
  jetsToken( 		consumes<edm::View<reco::Jet>>(				iConfig.getParameter<edm::InputTag>("jetsInputTag"))),
  candidatesToken(	consumes<edm::View<reco::Candidate>>(			iConfig.getParameter<edm::InputTag>("pfCandidatesInputTag"))),
  rhoToken( 		consumes<double>(					iConfig.getParameter<edm::InputTag>("rhoInputTag"))),
  vertexToken(    	consumes<reco::VertexCollection>(			iConfig.getParameter<edm::InputTag>("vertexInputTag"))),
  genJetsToken(    	consumes<reco::GenJetCollection>(			iConfig.getParameter<edm::InputTag>("genJetsInputTag"))),
  genParticlesToken(    consumes<reco::GenParticleCollection>(			iConfig.getParameter<edm::InputTag>("genParticlesInputTag"))),
  csvInputTag(									iConfig.getParameter<edm::InputTag>("csvInputTag")),
//qgVariablesInputTag(  							iConfig.getParameter<edm::InputTag>("qgVariablesInputTag")),
  jecService( 									iConfig.getParameter<std::string>("jec")),
  minJetPt(									iConfig.getUntrackedParameter<double>("minJetPt", 20.)),
  deltaRcut(									iConfig.getUntrackedParameter<double>("deltaRcut", 0.3)),
  pythia6(									iConfig.getUntrackedParameter<bool>("pythia6", false)),
  useQC(									iConfig.getUntrackedParameter<bool>("useQualityCuts", false))
{
  weStillNeedToCheckJets	  = true;
  weStillNeedToCheckJetCandidates = true;
  bTagToken		=	consumes<reco::JetTagCollection>(		edm::InputTag(csvInputTag));
/*qgToken		= 	consumes<edm::ValueMap<float>>(			edm::InputTag(qgVariablesInputTag.label(), "qgLikelihood"));
  axis2Token		= 	consumes<edm::ValueMap<float>>(			edm::InputTag(qgVariablesInputTag.label(), "axis2Likelihood"));
  multToken		= 	consumes<edm::ValueMap<int>>(			edm::InputTag(qgVariablesInputTag.label(), "multLikelihood"));
  ptDToken		= 	consumes<edm::ValueMap<float>>(			edm::InputTag(qgVariablesInputTag.label(), "ptDLikelihood"));*/
//qglcalc 		= 	new QGLikelihoodCalculator("/user/tomc/QGTagger/CMSSW_7_0_9_patch1/src/QGDev/qgMiniTuple/data/pdfQG_AK4chs_antib_13TeV.root");
}


template<class a, class b> int qgMiniTuple::countInCone(a center, b objectsToCount){
  int counter = 0;
  for(auto object = objectsToCount->begin(); object != objectsToCount->end(); ++object){
    if(reco::deltaR(*center, *object) < deltaRcut) ++counter;
  }
  return counter;
}


/*
 * Prepare for analyzing the event: choose pat or reco jets
 */
void qgMiniTuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  for(auto v : {closebyJetdR, closebyJetPt}) v->clear();

  nEvent	= (int) iEvent.id().event();

  edm::Handle<edm::View<reco::Jet>> 		jets;							iEvent.getByToken(jetsToken, 		jets);
  edm::Handle<edm::View<reco::Candidate>> 	candidates;						iEvent.getByToken(candidatesToken, 	candidates);
  edm::Handle<std::vector<PileupSummaryInfo>> 	pupInfo;						iEvent.getByLabel("addPileupInfo", 	pupInfo);
  edm::Handle<reco::VertexCollection> 		vertexCollection;					iEvent.getByToken(vertexToken, 		vertexCollection);
  edm::Handle<double> 				rhoHandle;						iEvent.getByToken(rhoToken, 		rhoHandle);
  edm::Handle<reco::GenJetCollection> 		genJets;						iEvent.getByToken(genJetsToken, 	genJets);
  edm::Handle<reco::GenParticleCollection> 	genParticles;						iEvent.getByToken(genParticlesToken, 	genParticles);
  edm::Handle<reco::JetTagCollection> 		bTagHandle;		if(!isPatJetCollection(jets))
                                                                        if(!csvInputTag.label().empty())iEvent.getByToken(bTagToken, 		bTagHandle);
/*edm::Handle<edm::ValueMap<float>> 		qgHandle;		if(!isPatJetCollection(jets))	iEvent.getByToken(qgToken, 		qgHandle);
  edm::Handle<edm::ValueMap<float>> 		axis2Handle; 		if(!isPatJetCollection(jets))	iEvent.getByToken(axis2Token, 		axis2Handle);
  edm::Handle<edm::ValueMap<float>> 		ptDHandle;  		if(!isPatJetCollection(jets))	iEvent.getByToken(multToken, 		multHandle);
  edm::Handle<edm::ValueMap<int>> 		multHandle;   		if(!isPatJetCollection(jets))	iEvent.getByToken(ptDToken,		ptDHandle);
*/

  if(!isPatJetCollection(jets)) JEC = JetCorrector::getJetCorrector(jecService, iSetup);

  nPriVtxs 	= vertexCollection->size();
  rho 		= (float) *rhoHandle;
  nPileUp 	= -1;
  if(pupInfo.isValid()){
    auto PVI = pupInfo->begin();
    while(PVI->getBunchCrossing() != 0 && PVI != pupInfo->end()) ++PVI;
    if(PVI != pupInfo->end()) nPileUp = PVI->getPU_NumInteractions();
  }

 
  if(genJets->size() > 2){
    auto jet1 = genJets->begin();
    auto jet2 = genJets->begin() + 1;
    auto jet3 = genJets->begin() + 2;
    balanced = (jet3->pt() < 0.15*(jet1->pt()+jet2->pt()));
  } else balanced = true;

  for(auto jet = jets->begin();  jet != jets->end(); ++jet){
    if(jet == jets->begin() + 2) balanced = false;

    if(isPatJetCollection(jets)) pt = jet->pt();
    else			 pt = jet->pt()*JEC->correction(*jet, iEvent, iSetup);
    if(pt < minJetPt) continue;


    // Closeby jet study variables
    for(auto otherJet = jets->begin(); otherJet != jets->end(); ++otherJet){
      if(otherJet == jet) continue;
      float dR = reco::deltaR(*jet, *otherJet);
      if(dR > 1.2) continue;
      closebyJetdR->push_back(dR);
      if(isPatJetCollection(jets)) closebyJetPt->push_back(otherJet->pt());
      else			   closebyJetPt->push_back(otherJet->pt()*JEC->correction(*otherJet, iEvent, iSetup));
    }

    // Parton Id matching
    float deltaRmin = 999;
    auto matchedGenParticle = genParticles->end();
    for(auto genParticle = genParticles->begin(); genParticle != genParticles->end(); ++genParticle){
      if(genParticle->fromHardProcessBeforeFSR()) continue; 							//New status flag, supposed to be similar to status 3 (pythia6) / status 23 (pythia8)
      if(abs(genParticle->pdgId()) > 5 && abs(genParticle->pdgId()) != 21) continue;				//Only keep quarks and gluons
      float thisDeltaR = reco::deltaR(*genParticle, *jet);
      if(thisDeltaR < deltaRmin){
        deltaRmin = thisDeltaR;
        matchedGenParticle = genParticle;
      }
    }
    if(deltaRmin < deltaRcut){
      matchedJet		= true;
      partonId 			= matchedGenParticle->pdgId();
      nJetsForGenParticle 	= countInCone(matchedGenParticle, jets);
      nGenJetsForGenParticle 	= countInCone(matchedGenParticle, genJets);
      if(matchedGenParticle->numberOfMothers() == 1){								//Very experimental, but first tests shows it's good at finding W's and t's
        motherId		= matchedGenParticle->mother()->pdgId();					//A bit more difficult for QCD, where it's sometimes a quark, decaying into
        motherMass		= matchedGenParticle->mother()->mass();						//quark+gluon, and sometimes just a proton with a lot of other QCD mess and 
      } else {													//sometimes 2 mothers (mostly two quarks recoiling each other, but sometimes
        motherId		= 0;										//also two quarks going into two gluons etc...)
        motherMass		= 0;
      }
    } else {
      matchedJet 		= false;
      partonId 			= 0; 
      nJetsForGenParticle 	= 0;
      nGenJetsForGenParticle 	= 0;
      motherId		= 0;
      motherMass		= 0;
    }
    if(!matchedJet) continue;
    nGenJetsInCone 		= countInCone(jet, genJets);

    if(isPatJetCollection(jets)){
      auto patJet 	= static_cast<const pat::Jet*> (&*jet);
      jetIdLevel	= jetId(patJet) + jetId(patJet, false, true) + jetId(patJet, true); 
      bTag		= patJet->bDiscriminator(csvInputTag.label());
    } else {
      edm::RefToBase<reco::Jet> jetRef(edm::Ref<edm::View<reco::Jet>>(jets, (jet - jets->begin())));
      auto recoJet 	= static_cast<const reco::PFJet*>(&*jet);
      jetIdLevel	= jetId(recoJet) + jetId(recoJet, false, true) + jetId(recoJet, true); 
      bTag		= csvInputTag.label().empty() ? 0 : (*bTagHandle)[jetRef];
/*    qg		= (*qgHandle)[jetRef];
      axis2		= (*axis2Handle)[jetRef];
      mult		= (*multHandle)[jetRef];
      ptD		= (*ptDHandle)[jetRef];*/
    }

    std::tie(mult, ptD, axis2) 	= calcVariables(&*jet, vertexCollection);
    axis2 			= -std::log(axis2);
    eta				= jet->eta();
//  qg 				= qglcalc->computeQGLikelihood(pt, eta, rho, {(float) mult, ptD, axis2});

    ptDoubleCone = 0;
    for(auto candidate = candidates->begin(); candidate != candidates->end(); ++candidate){
      if(reco::deltaR(*candidate, *jet) < 0.8) ptDoubleCone += candidate->pt();
    }
 
    tree->Fill();
  }
}


/*
 * Begin job: create vectors and set up tree
 */
void qgMiniTuple::beginJob(){
  for(auto v : {&closebyJetdR, &closebyJetPt}) *v = new std::vector<float>();

  tree = fs->make<TTree>("qgMiniTuple","qgMiniTuple");
  tree->Branch("nEvent" ,			&nEvent, 			"nEvent/I");
  tree->Branch("nPileUp",			&nPileUp, 			"nPileUp/I");
  tree->Branch("nPriVtxs",			&nPriVtxs, 			"nPriVtxs/I");
  tree->Branch("rho" ,				&rho, 				"rho/F");
  tree->Branch("pt" ,				&pt,				"pt/F");
  tree->Branch("eta",				&eta,				"eta/F");
  tree->Branch("axis2",				&axis2,				"axis2/F");
  tree->Branch("ptD",				&ptD,				"ptD/F");
  tree->Branch("mult",				&mult,				"mult/I");
  tree->Branch("bTag",				&bTag,				"bTag/F");
  tree->Branch("partonId",			&partonId,			"partonId/I");
  tree->Branch("motherId",			&motherId,			"motherId/I");
  tree->Branch("motherMass",			&motherMass,			"motherMass/F");
  tree->Branch("jetIdLevel",			&jetIdLevel,			"jetIdLevel/I");
  tree->Branch("nGenJetsInCone",		&nGenJetsInCone,		"nGenJetsInCone/I");
  tree->Branch("matchedJet",			&matchedJet,			"matchedJet/O");
  tree->Branch("balanced",			&balanced,			"balanced/O");
  tree->Branch("ptDoubleCone",			&ptDoubleCone,			"ptDoubleCone/F");
  tree->Branch("nGenJetsForGenParticle",	&nGenJetsForGenParticle,	"nGenJetsForGenParticle/I");
  tree->Branch("nJetsForGenParticle",   	&nJetsForGenParticle,   	"nJetsForGenParticle/I");
  tree->Branch("closebyJetdR",			"vector<float>",		&closebyJetdR);
  tree->Branch("closebyJetPt",			"vector<float>",		&closebyJetPt);
}


/*
 * End of jobs: delete vectors
 */
void qgMiniTuple::endJob(){
 for(auto v : {closebyJetdR, closebyJetPt}) delete v;
}


/*
 * Function to tell us if we are using pat::Jet or reco::PFJet
 */
bool qgMiniTuple::isPatJetCollection(const edm::Handle<edm::View<reco::Jet>>& jets){
  if(weStillNeedToCheckJets){
    if(typeid(pat::Jet)==typeid(*(jets->begin())))         weAreUsingPatJets = true;
    else if(typeid(reco::PFJet)==typeid(*(jets->begin()))) weAreUsingPatJets = false;
    else throw cms::Exception("WrongJetCollection", "Expecting pat::Jet or reco::PFJet");
    weStillNeedToCheckJets = false;
  }
  return weAreUsingPatJets;
}


/*
 * Function to tell us if we are using packedCandidates, only test for first candidate
 */
bool qgMiniTuple::isPackedCandidate(const reco::Candidate* candidate){
  if(weStillNeedToCheckJetCandidates){
    if(typeid(pat::PackedCandidate)==typeid(*candidate))   weAreUsingPackedCandidates = true;
    else if(typeid(reco::PFCandidate)==typeid(*candidate)) weAreUsingPackedCandidates = false;
    else throw cms::Exception("WrongJetCollection", "Jet constituents are not particle flow candidates");
    weStillNeedToCheckJetCandidates = false;
  }
  return weAreUsingPackedCandidates;
}


/* 
 * Calculation of axis2, mult and ptD
 */
std::tuple<int, float, float> qgMiniTuple::calcVariables(const reco::Jet *jet, edm::Handle<reco::VertexCollection>& vC){
  float sum_weight = 0., sum_deta = 0., sum_dphi = 0., sum_deta2 = 0., sum_dphi2 = 0., sum_detadphi = 0., sum_pt = 0.;
  int mult = 0;

  //Loop over the jet constituents
  for(auto daughter : jet->getJetConstituentsQuick()){
    if(isPackedCandidate(daughter)){											//packed candidate situation
      auto part = static_cast<const pat::PackedCandidate*>(daughter);

      if(part->charge()){
        if(!(part->fromPV() > 1 && part->trackHighPurity())) continue;
        if(useQC){
          if((part->dz()*part->dz())/(part->dzError()*part->dzError()) > 25.) continue;
          if((part->dxy()*part->dxy())/(part->dxyError()*part->dxyError()) < 25.) ++mult;
        } else ++mult;
      } else {
        if(part->pt() < 1.0) continue;
        ++mult;
      }
    } else {
      auto part = static_cast<const reco::PFCandidate*>(daughter);

      reco::TrackRef itrk = part->trackRef();
      if(itrk.isNonnull()){												//Track exists --> charged particle
        auto vtxLead  = vC->begin();
        auto vtxClose = vC->begin();											//Search for closest vertex to track
        for(auto vtx = vC->begin(); vtx != vC->end(); ++vtx){
          if(fabs(itrk->dz(vtx->position())) < fabs(itrk->dz(vtxClose->position()))) vtxClose = vtx;
        }
        if(!(vtxClose == vtxLead && itrk->quality(reco::TrackBase::qualityByName("highPurity")))) continue;

        if(useQC){													//If useQC, require dz and d0 cuts
          float dz = itrk->dz(vtxClose->position());
          float d0 = itrk->dxy(vtxClose->position());
          float dz_sigma_square = pow(itrk->dzError(),2) + pow(vtxClose->zError(),2);
          float d0_sigma_square = pow(itrk->d0Error(),2) + pow(vtxClose->xError(),2) + pow(vtxClose->yError(),2);
          if(dz*dz/dz_sigma_square > 25.) continue;
          if(d0*d0/d0_sigma_square < 25.) ++mult;
        } else ++mult;
      } else {														//No track --> neutral particle
        if(part->pt() < 1.0) continue;											//Only use neutrals with pt > 1 GeV
        ++mult;
      }
    }

    float deta   = daughter->eta() - jet->eta();
    float dphi   = reco::deltaPhi(daughter->phi(), jet->phi());
    float partPt = daughter->pt();
    float weight = partPt*partPt;

    sum_weight   += weight;
    sum_pt       += partPt;
    sum_deta     += deta*weight;
    sum_dphi     += dphi*weight;
    sum_deta2    += deta*deta*weight;
    sum_detadphi += deta*dphi*weight;
    sum_dphi2    += dphi*dphi*weight;
  }

  //Calculate axis2 and ptD
  float a = 0., b = 0., c = 0.;
  float ave_deta = 0., ave_dphi = 0., ave_deta2 = 0., ave_dphi2 = 0.;
  if(sum_weight > 0){
    ave_deta  = sum_deta/sum_weight;
    ave_dphi  = sum_dphi/sum_weight;
    ave_deta2 = sum_deta2/sum_weight;
    ave_dphi2 = sum_dphi2/sum_weight;
    a         = ave_deta2 - ave_deta*ave_deta;                          
    b         = ave_dphi2 - ave_dphi*ave_dphi;                          
    c         = -(sum_detadphi/sum_weight - ave_deta*ave_dphi);                
  }
  float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
  float axis2 = (a+b-delta > 0 ?  sqrt(0.5*(a+b-delta)) : 0);
  float ptD   = (sum_weight > 0 ? sqrt(sum_weight)/sum_pt : 0);
  return std::make_tuple(mult, ptD, axis2);
}


/*
 * Calculate jetId for levels loose, medium and tight
 */
template <class jetClass> bool qgMiniTuple::jetId(const jetClass *jet, bool tight, bool medium){
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

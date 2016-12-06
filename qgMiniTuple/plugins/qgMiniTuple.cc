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
#include "TH2F.h"


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
  std::tuple<int, int, int, int, int, float, float, float, float, float, float, float, float> 	calcVariables(const reco::Jet *jet, edm::Handle<reco::VertexCollection>& vC);
  bool 				isPatJetCollection(const edm::Handle<edm::View<reco::Jet>>& jets);
  bool 				isPackedCandidate(const reco::Candidate* candidate);
  template<class a, class b> int 	countInCone(a center, b objectsToCount);
  template<class a> bool		isBalanced(a objects);
  int 				getPileUp(edm::Handle<std::vector<PileupSummaryInfo>>& pupInfo);
  //  float deltaPhi(Float_t p1, Float_t p2);
  //  float deltaR(Float_t deta, Float_t dphi);
  reco::GenParticleCollection::const_iterator getMatchedGenParticle(const reco::Jet *jet, edm::Handle<reco::GenParticleCollection>& genParticles);

  edm::EDGetTokenT<edm::View<reco::Jet>> 		jetsToken;
  edm::EDGetTokenT<edm::View<reco::Candidate>> 	candidatesToken;
  edm::EDGetTokenT<double> 				rhoToken;
  edm::EDGetTokenT<reco::VertexCollection> 		vertexToken;
  edm::EDGetTokenT<reco::GenJetCollection> 		genJetsToken;
  edm::EDGetTokenT<reco::GenParticleCollection> 	genParticlesToken;
  edm::EDGetTokenT<reco::JetTagCollection> 		bTagToken;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puToken;
  edm::InputTag					csvInputTag;
  std::string 					jecService;
  /*    edm::InputTag 					qgVariablesInputTag;
	edm::EDGetTokenT<edm::ValueMap<float>> 		qgToken, axis2Token, ptDToken;
	edm::EDGetTokenT<edm::ValueMap<int>> 		multToken;*/
  const double 					minJetPt, deltaRcut;
  const bool 					useQC;

  const JetCorrector 				*JEC;
  edm::Service<TFileService> 			fs;
  TTree 						*tree;
  TTree 						*jtree;
  //    QGLikelihoodCalculator 				*qglcalc;

  float rho, pt, eta, axis2, ptD, bTag, ptDoubleCone, motherMass, mass, pt_dr, pt_dr2, pt_dr_log, sumpt_dr, sumpt_dr2, sumpt_dr_log, gluon_ll, quark_ll;
  int nEvent, nPileUp, nPriVtxs, mult, mult_new, pmult, cmult, nmult, chadmult, nhadmult, partonId, jetIdLevel, nGenJetsInCone, nGenJetsForGenParticle, nJetsForGenParticle, motherId, charge, jetarea;
  bool matchedJet, balanced;
  std::vector<float> *closebyJetdR, *closebyJetPt;


  std::vector<int> *cnEvent, *cpdgId, *ccharge, *cfromPV, *ctrackHighPurity, *cjetid, *cjetpdg;
  std::vector<float> *crho, *cpt, *ceta, *cphi, *cmass, *cdz, *cdxy, *cjetpt, *cjeteta, *cjetphi;


  TH2F* gluon_pdf;
  TH2F* quark_pdf;

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
  puToken ( 		consumes<std::vector<PileupSummaryInfo> > (edm::InputTag("slimmedAddPileupInfo")) ),
  csvInputTag(									iConfig.getParameter<edm::InputTag>("csvInputTag")),
  //qgVariablesInputTag(  							iConfig.getParameter<edm::InputTag>("qgVariablesInputTag")),
  jecService( 									iConfig.getParameter<std::string>("jec")),
  minJetPt(									iConfig.getUntrackedParameter<double>("minJetPt", 20.)),
  deltaRcut(									iConfig.getUntrackedParameter<double>("deltaRcut", 0.3)),
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
  
  TFile *input = new TFile("/mnt/t3nfs01/data01/shome/ytakahas/work/QG_disc/CMSSW_8_1_0_pre11/src/QGDev/qgMiniTuple/test/pdf_prod.root");

  gluon_pdf = (TH2F*) input->Get("gluon");
  quark_pdf = (TH2F*) input->Get("quark");

  input->Close();

}


/*
 * Prepare for analyzing the event: choose pat or reco jets
 */
void qgMiniTuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  for(auto v : {closebyJetdR, closebyJetPt}) v->clear();


  nEvent	= (int) iEvent.id().event();

  edm::Handle<edm::View<reco::Jet>> 		jets;							iEvent.getByToken(jetsToken, 		jets);
  edm::Handle<edm::View<reco::Candidate>> 	candidates;						iEvent.getByToken(candidatesToken, 	candidates);
  edm::Handle<std::vector<PileupSummaryInfo>> 	pupInfo;						iEvent.getByToken(puToken, 	pupInfo);
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

  // Get number of primary vertices, pile-up and rho
  nPriVtxs 	= vertexCollection->size();
  rho 		= (float) *rhoHandle;
  nPileUp 	= getPileUp(pupInfo);

  // Start jet loop (the tree is filled for each jet separately)
  balanced = isBalanced(genJets);										// Check if first two generator jets are balanced

  int counter = 0;

  for(auto jet = jets->begin();  jet != jets->end(); ++jet){
    
    if(jet == jets->begin() + 2) balanced = false;								// 3rd jet and higher are never balanced

    if(isPatJetCollection(jets)) pt = jet->pt();								// If miniAOD, jets are already corrected
    else			 pt = jet->pt()*JEC->correction(*jet, iEvent, iSetup);				// If RECO, we correct them on the fly
    if(pt < minJetPt) continue;

    // Closeby jet study variables (i.e. dR and pt of other jets within 1.2)
    for(auto otherJet = jets->begin(); otherJet != jets->end(); ++otherJet){
      if(otherJet == jet) continue;
      float dR = reco::deltaR(*jet, *otherJet);
      if(dR > 1.2) continue;
      closebyJetdR->push_back(dR);
      if(isPatJetCollection(jets)) closebyJetPt->push_back(otherJet->pt());
      else			   closebyJetPt->push_back(otherJet->pt()*JEC->correction(*otherJet, iEvent, iSetup));
    }

    // Parton Id matching
    auto matchedGenParticle = getMatchedGenParticle(&*jet, genParticles);
    matchedJet = (matchedGenParticle != genParticles->end());
    if(matchedJet){
      partonId 			= matchedGenParticle->pdgId();
      nJetsForGenParticle 	= countInCone(matchedGenParticle, jets);
      nGenJetsForGenParticle 	= countInCone(matchedGenParticle, genJets);
      if(matchedGenParticle->numberOfMothers() == 1){								// Very experimental, but first tests shows it's good at finding W's and t's
        motherId		= matchedGenParticle->mother()->pdgId();					// A bit more difficult for QCD, where it's sometimes a quark, decaying into
        motherMass		= matchedGenParticle->mother()->mass();						// quark+gluon, and sometimes just a proton with a lot of other QCD mess and 
      } else {													// sometimes 2 mothers (mostly two quarks recoiling each other, but sometimes
        motherId		= 0;										// also two quarks going into two gluons etc...)
        motherMass		= 0;
      }
    } else {
      partonId 			= 0; 
      nJetsForGenParticle 	= 0;
      nGenJetsForGenParticle 	= 0;
      motherId			= 0;
      motherMass		= 0;
      continue;													// To keep the tuples small, we only save matched jets
    }
    nGenJetsInCone 		= countInCone(jet, genJets);



    for(auto v : {cnEvent, cpdgId, ccharge, cfromPV, ctrackHighPurity, cjetid, cjetpdg}) v->clear();
    for(auto v : {crho, cpt, ceta, cphi, cmass, cdz, cdxy, cjetpt, cjeteta, cjetphi}) v->clear();

    TH2F* h_map = new TH2F("h_map", "h_map", 120,-0.5,0.5,120,-0.5,0.5);
    h_map->Sumw2();

    for(auto daughter : jet->getJetConstituentsQuick()){
      if(isPackedCandidate(daughter)){								    	//packed candidate situation
	auto part = static_cast<const pat::PackedCandidate*>(daughter);

	//	if(part->pt() > 1.){
	  cnEvent->push_back(nEvent);
	  crho->push_back(rho);
	  ccharge->push_back(part->charge());
	  cfromPV->push_back(part->fromPV());
	  ctrackHighPurity->push_back(part->trackHighPurity());
	  cdz->push_back(part->dz());
	  cdxy->push_back(part->dxy());
	  cpt->push_back(part->pt());
	  ceta->push_back(part->eta());
	  cphi->push_back(part->phi());
	  cmass->push_back(part->mass());
	  cpdgId->push_back(part->pdgId());
	  cjetid->push_back(counter);
	  cjetpt->push_back(pt);
	  cjeteta->push_back(jet->eta());
	  cjetphi->push_back(jet->phi());
	  cjetpdg->push_back(partonId);

	  float dphi = reco::deltaPhi(part->phi(), jet->phi());
	  float deta = part->eta() - jet->eta();

	  h_map->Fill(deta, dphi, part->pt()/jet->pt());
	  
	  //	}
      }
    }

    jtree->Fill();

    counter++;

    // perform 2D fitting, manually 

    float_t chi2_gluon = 0;
    float_t chi2_quark = 0;
    int npoints = 0;

    for(int ix=1; ix<h_map->GetXaxis()->GetNbins()+1; ix++){
      for(int iy=1; iy<h_map->GetYaxis()->GetNbins()+1; iy++){
	float val = h_map->GetBinContent(ix, iy);
	if(val == 0) continue;
	float val_err = h_map->GetBinError(ix, iy);	
	float diff_gluon = val - gluon_pdf->GetBinContent(ix,iy);
	float diff_quark = val - quark_pdf->GetBinContent(ix,iy);

	chi2_quark += (diff_quark*diff_quark);
	chi2_gluon += (diff_gluon*diff_gluon);
	npoints += 1;
      }
    }

    gluon_ll = sqrt(chi2_gluon)/npoints;
    quark_ll = sqrt(chi2_quark)/npoints;


    if(isPatJetCollection(jets)){
      auto patJet 	= static_cast<const pat::Jet*> (&*jet);
      jetIdLevel	= jetId(patJet) + jetId(patJet, false, true) + jetId(patJet, true); 
      bTag		= patJet->bDiscriminator(csvInputTag.label());
      chadmult                    = patJet->chargedMultiplicity();
      nhadmult                    = patJet->neutralMultiplicity();

    } else {
      edm::RefToBase<reco::Jet> jetRef(edm::Ref<edm::View<reco::Jet>>(jets, (jet - jets->begin())));
      auto recoJet 	= static_cast<const reco::PFJet*>(&*jet);
      jetIdLevel	= jetId(recoJet) + jetId(recoJet, false, true) + jetId(recoJet, true); 
      bTag		= csvInputTag.label().empty() ? 0 : (*bTagHandle)[jetRef];
      chadmult                    = recoJet->chargedMultiplicity();
      nhadmult                    = recoJet->neutralMultiplicity();
      /*    qg		= (*qgHandle)[jetRef];
	    axis2		= (*axis2Handle)[jetRef];
	    mult		= (*multHandle)[jetRef];
	    ptD		= (*ptDHandle)[jetRef];*/
    }

    //    std::cout << "(chad_multiplicity, nhad_multiplicity) = " << chadmult << " " << nhadmult << std::endl;


    std::tie(mult, mult_new, pmult, cmult, nmult, ptD, axis2, pt_dr, pt_dr2, pt_dr_log, sumpt_dr, sumpt_dr2, sumpt_dr_log) = calcVariables(&*jet, vertexCollection);
    axis2 			= -std::log(axis2);
    eta				= jet->eta();
    charge                      = jet->charge();
    mass                        = jet->mass();
    jetarea                     = jet->jetArea();
    //  qg 				= qglcalc->computeQGLikelihood(pt, eta, rho, {(float) mult, ptD, axis2});
    if(mult < 2) continue;  

    ptDoubleCone = 0;
    for(auto candidate = candidates->begin(); candidate != candidates->end(); ++candidate){
      if(reco::deltaR(*candidate, *jet) < 0.8) ptDoubleCone += candidate->pt();
    }
 


    tree->Fill();
  } // jet loop
}


/*
 * Begin job: create vectors and set up tree
 */
void qgMiniTuple::beginJob(){
  for(auto v : {&closebyJetdR, &closebyJetPt}) *v = new std::vector<float>();
  for(auto v : {&cnEvent, &cpdgId, &ccharge, &cfromPV, &ctrackHighPurity, &cjetid, &cjetpdg}) *v = new std::vector<int>();
  for(auto v : {&crho, &cpt, &ceta, &cphi, &cmass, &cdz, &cdxy, &cjetpt, &cjeteta, &cjetphi}) *v = new std::vector<float>();


  tree = fs->make<TTree>("qgMiniTuple","qgMiniTuple");
  tree->Branch("nEvent" ,			&nEvent, 			"nEvent/I");
  tree->Branch("nPileUp",			&nPileUp, 			"nPileUp/I");
  tree->Branch("nPriVtxs",			&nPriVtxs, 			"nPriVtxs/I");
  tree->Branch("rho" ,				&rho, 				"rho/F");
  tree->Branch("pt" ,				&pt,				"pt/F");
  tree->Branch("eta",				&eta,				"eta/F");
  tree->Branch("axis2",				&axis2,				"axis2/F");
  tree->Branch("charge",		        &charge,			"charge/I");
  tree->Branch("mass",				&mass,				"mass/F");
  tree->Branch("pt_dr",				&pt_dr,				"pt_dr/F");
  tree->Branch("pt_dr2",			&pt_dr2,			"pt_dr2/F");
  tree->Branch("pt_dr_log",			&pt_dr_log,			"pt_dr_log/F");
  tree->Branch("sumpt_dr",			&sumpt_dr,			"sumpt_dr/F");
  tree->Branch("sumpt_dr2",			&sumpt_dr2,			"sumpt_dr2/F");
  tree->Branch("sumpt_dr_log",			&sumpt_dr_log,			"sumpt_dr_log/F");
  tree->Branch("jetarea",			&jetarea,			"jetarea/F");
  tree->Branch("gluon_ll",			&gluon_ll,			"gluon_ll/F");
  tree->Branch("quark_ll",			&quark_ll,			"quark_ll/F");
  tree->Branch("ptD",				&ptD,				"ptD/F");
  tree->Branch("mult",				&mult,				"mult/I");
  tree->Branch("mult_new",			&mult_new,	      	        "mult_new/I");
  tree->Branch("pmult",				&pmult,				"pmult/I");
  tree->Branch("cmult",				&cmult,				"cmult/I");
  tree->Branch("nmult",				&nmult,				"nmult/I");
  tree->Branch("chadmult",			&chadmult,			"chadmult/I");
  tree->Branch("nhadmult",			&nhadmult,			"nhadmult/I");
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


  jtree = fs->make<TTree>("perJet","perJet");
  jtree->Branch("cnEvent",                      "vector<int>",			&cnEvent);
  jtree->Branch("crho",                         "vector<float>",		&crho);
  jtree->Branch("cpt",                          "vector<float>",		&cpt);
  jtree->Branch("ceta",                         "vector<float>",		&ceta);
  jtree->Branch("cphi",                         "vector<float>",		&cphi);
  jtree->Branch("cmass",                        "vector<float>",		&cmass);
  jtree->Branch("cpdgId",                       "vector<int>",			&cpdgId);
  jtree->Branch("ccharge",                      "vector<int>",			&ccharge);
  jtree->Branch("cfromPV",                      "vector<int>",			&cfromPV);
  jtree->Branch("ctrackHighPurity",             "vector<int>",		        &ctrackHighPurity);
  jtree->Branch("cdz",                          "vector<float>",		&cdz);
  jtree->Branch("cdxy",                         "vector<float>",		&cdxy);
  jtree->Branch("cjetid",                       "vector<int>",			&cjetid);
  jtree->Branch("cjetpt",                       "vector<float>",		&cjetpt);
  jtree->Branch("cjeteta",                      "vector<float>",		&cjeteta);
  jtree->Branch("cjetphi",                      "vector<float>",		&cjetphi);
  jtree->Branch("cjetpdg",                      "vector<int>",		        &cjetpdg);


}


/*
 * End of jobs: delete vectors
 */
void qgMiniTuple::endJob(){
  for(auto v : {closebyJetdR, closebyJetPt}) delete v;
  for(auto v : {&cnEvent, &cpdgId, &ccharge, &cfromPV, &ctrackHighPurity, &cjetid, &cjetpdg}) delete v;
  for(auto v : {&crho, &cpt, &ceta, &cphi, &cmass, &cdz, &cdxy, &cjetpt, &cjeteta, &cjetphi}) delete v;
  
  //  fs->Write();
  //  fs->Close();
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
std::tuple<int, int, int, int, int, float, float, float, float, float, float, float, float> qgMiniTuple::calcVariables(const reco::Jet *jet, edm::Handle<reco::VertexCollection>& vC){
  float sum_weight = 0., sum_deta = 0., sum_dphi = 0., sum_deta2 = 0., sum_dphi2 = 0., sum_detadphi = 0., sum_pt = 0.;
  int mult = 0, mult_new = 0;
  int pmult = 0, cmult = 0, nmult = 0;
  float pt_dr = 0., pt_dr2 = 0., pt_dr_log = 0;
  float sumpt_dr = 0., sumpt_dr2 = 0., sumpt_dr_log = 0;

  //  std::cout << "enter calcVariable !" << std::endl;

  //Loop over the jet constituents
  for(auto daughter : jet->getJetConstituentsQuick()){
    if(isPackedCandidate(daughter)){								    	//packed candidate situation
      auto part = static_cast<const pat::PackedCandidate*>(daughter);

      if(part->charge()){
	//	std::cout << "charge = " << part->charge() << std::endl;

        if(!(part->fromPV() > 1 && part->trackHighPurity())) continue;
        if(useQC){
          if((part->dz()*part->dz())/(part->dzError()*part->dzError()) > 25.) continue;
          if((part->dxy()*part->dxy())/(part->dxyError()*part->dxyError()) < 25.) ++mult;
        } else{ ++mult; 
	  //	  std::cout << "added charged particle : " << part->pt() << " " << part->pdgId() << std::endl;
	};

      } else {
        if(part->pt() < 1.0) continue;
	//	std::cout << "add multi : pdg = " << part->pdgId() << " " << part->pt() << std::endl;
        ++mult;
	//	std::cout << "added neutral particle : " << part->pt() << " " << part->pdgId() << std::endl;
      }


      Int_t pdg = abs(part->pdgId());
      
      //      std::cout << "isPackedCandidate : pdg ID = " << pdg << " " << part->pt() << std::endl;
      if(part->pt() > 1.0){
	mult_new ++;
	if(pdg<=9 || pdg==21) continue;
	if(pdg>=11 && pdg<=16) continue;
	if(pdg==211){ cmult++;}
	else if(pdg==310 || pdg==130){ nmult++;}
	else if(pdg==22){pmult++;}
	else std::cout << "What !!!!!!!!!!!!!!!!!!!!!!!! " << pdg << std::endl;
      }



      if(part->pt() > 1.0){

	//	float min_dr = 100;
	
//	for(auto daughter : jet->getJetConstituentsQuick()){
//	  auto part2 = static_cast<const pat::PackedCandidate*>(daughter);
//	  if(part2->pt()==part->pt() && part2->eta()==part->eta()) continue;
//	  //	  float dphi = reco::deltaPhi(part->phi(), part2->phi());
//	  //	  float deta = part->eta(), part2->eta();
	float dr = reco::deltaR(*jet, *part);
	
	pt_dr += (part->pt()/dr);
	pt_dr2 += (part->pt()/(dr*dr)); 
	pt_dr_log += log(part->pt()/dr);
	sumpt_dr += 1/dr;
	sumpt_dr2 += 1/(dr*dr);
	sumpt_dr_log += log(1/dr);
//
//	if(dr < min_dr){
//	  min_dr = dr;
//	}	  
//	}

//	if(min_dr !=100) pt_dr += 

      }
      
    } else {

      std::cout << "!!!!!!!!!!!!!!!!!!!!!! PF Candidate Type !!" << std::endl;
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

      std::cout << "isPFCandidate : pdg ID = " <<  part->pdgId() << std::endl;
      if(part->pt() > 1.0){

	float dr = reco::deltaR(*jet, *part);
	std::cout << "dr = " << dr << std::endl;
	pt_dr += (part->pt()/dr);
	pt_dr2 += (part->pt()/(dr*dr)); 
	pt_dr_log += log(part->pt()/dr);
	sumpt_dr += 1/dr;
	sumpt_dr2 += 1/(dr*dr);
	sumpt_dr_log += log(1/dr);
//
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
  return std::make_tuple(mult, mult_new, pmult, cmult, nmult, ptD, axis2, pt_dr, pt_dr2, pt_dr_log, sumpt_dr, sumpt_dr2, sumpt_dr_log);
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


/*
 * Count objects around another object within dR < deltaRcut
 */
template<class a, class b> int qgMiniTuple::countInCone(a center, b objectsToCount){
  int counter = 0;
  for(auto object = objectsToCount->begin(); object != objectsToCount->end(); ++object){
    if(reco::deltaR(*center, *object) < deltaRcut) ++counter;
  }
  return counter;
}


/*
 * Are two leading objects in the collection balanced ?
 */
template<class a> bool qgMiniTuple::isBalanced(a objects){
  if(objects->size() > 2){
    auto object1 = objects->begin();
    auto object2 = objects->begin() + 1;
    auto object3 = objects->begin() + 2;
    return (object3->pt() < 0.15*(object1->pt()+object2->pt()));
  } else return true;
}


/*
 * Get nPileUp
 */
int qgMiniTuple::getPileUp(edm::Handle<std::vector<PileupSummaryInfo>>& pupInfo){
  if(!pupInfo.isValid()) return -1;
  auto PVI = pupInfo->begin();
  while(PVI->getBunchCrossing() != 0 && PVI != pupInfo->end()) ++PVI;
  if(PVI != pupInfo->end()) return PVI->getPU_NumInteractions();
  else return -1;
}


/*
 * Parton Id matching (physics definition)
 */
reco::GenParticleCollection::const_iterator qgMiniTuple::getMatchedGenParticle(const reco::Jet *jet, edm::Handle<reco::GenParticleCollection>& genParticles){
  float deltaRmin = 999;
  auto matchedGenParticle = genParticles->end();
  for(auto genParticle = genParticles->begin(); genParticle != genParticles->end(); ++genParticle){
    if(!genParticle->isHardProcess()) continue;								// This status flag is exactly the pythia8 status-23 we need (i.e. the same as genParticles->status() == 23), probably also ok to use for other generators
    if(abs(genParticle->pdgId()) > 5 && abs(genParticle->pdgId() != 21)) continue;			// Only consider udscb quarks and gluons
    float thisDeltaR = reco::deltaR(*genParticle, *jet);
    if(thisDeltaR < deltaRmin && thisDeltaR < deltaRcut){
      deltaRmin = thisDeltaR;
      matchedGenParticle = genParticle;
    }
  }
  return matchedGenParticle;
}

//float qgMiniTuple::deltaPhi(Float_t p1, Float_t p2){
//
//  Float_t res = p1 - p2;
//  while(res > TMath::Pi()){
//    res -= 2*TMath::Pi();
//  }
//  while(res < -TMath::Pi()){
//    res += 2*TMath::Pi();
//  }
//  
//  return res;
//
//}

//float qgMiniTuple::deltaR(Float_t deta, Float_t dphi){
//  
//  return TMath::Sqrt(TMath::Power(deta,2) + TMath::Power(dphi,2));
//
//}


DEFINE_FWK_MODULE(qgMiniTuple);

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
      template <class jetClass> void calcVariables(const jetClass *jet, float& axis2_NoQC, float& ptD_NoQC, int& mult_NoQC);

      edm::EDGetTokenT<double> rhoToken;
      edm::EDGetTokenT<pat::JetCollection> jetsToken;
      const double minJetPt;

      edm::Service<TFileService> fs;
      TTree *tree;

      float rho, pt, eta, axis2_NoQC, laxis2_NoQC, ptD_NoQC, bTag;
      int nRun, nLumi, nEvent, mult_NoQC, partonId;
};


qgMiniTupleForMiniAOD::qgMiniTupleForMiniAOD(const edm::ParameterSet& iConfig) :
  rhoToken( 		consumes<double>(			iConfig.getParameter<edm::InputTag>("rhoInputTag"))),
  jetsToken(    	consumes<pat::JetCollection>(		iConfig.getParameter<edm::InputTag>("jetsInputTag"))),
  minJetPt(							iConfig.getUntrackedParameter<double>("minJetPt", 10.))
{
}


void qgMiniTupleForMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  nRun 		= (int) iEvent.id().run();
  nLumi 	= (int) iEvent.id().luminosityBlock();
  nEvent	= (int) iEvent.id().event();

  edm::Handle<double> rho_;
  iEvent.getByToken(rhoToken, rho_);
  rho = (float) *rho_;

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetsToken, jets);

  for(auto jet = jets->begin();  jet != jets->end(); ++jet){
    if(jet->pt() < minJetPt) continue;
    if(jet->partonFlavour() == 0) continue;

    calcVariables(&*jet, axis2_NoQC, ptD_NoQC, mult_NoQC);
    axis2_NoQC 	= -std::log(axis2_NoQC);
    partonId	= jet->partonFlavour();
    pt		= jet->pt();
    eta		= jet->eta();
    bTag	= jet->bDiscriminator("combinedInclusiveSecondaryVertexBJetTags");
    tree->Fill();
  }
}


/// Calculation of axis2_NoQC, mult_NoQC and ptD_NoQC
template <class jetClass> void qgMiniTupleForMiniAOD::calcVariables(const jetClass *jet, float& axis2_NoQC_, float& ptD_NoQC_, int& mult_NoQC_){
  float sum_weight = 0., sum_deta = 0., sum_dphi = 0., sum_deta2 = 0., sum_dphi2 = 0., sum_detadphi = 0., sum_pt = 0.;
  int nChg_QC = 0, nNeutral_ptCut = 0;

  //Loop over the jet constituents
  for(int i = 0; i < jet->numberOfDaughters(); ++i){
    auto part = dynamic_cast<const pat::PackedCandidate*> (jet->daughter(i));
    if(part->charge()){
      if(part->fromPV() > 0 && part->trackHighPurity()){								// fromPV > 0 should always be the case for CHS
        nChg_QC++;
      } else continue;
    } else {
      if(part->pt() > 1.0) nNeutral_ptCut++;
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

  //Calculate axis2_NoQC and ptD_NoQC
  float a = 0., b = 0., c = 0.;
  float ave_deta = 0., ave_dphi = 0., ave_deta2 = 0., ave_dphi2 = 0.;
  if(sum_weight > 0){
    ptD_NoQC_ = sqrt(sum_weight)/sum_pt;
    ave_deta = sum_deta/sum_weight;
    ave_dphi = sum_dphi/sum_weight;
    ave_deta2 = sum_deta2/sum_weight;
    ave_dphi2 = sum_dphi2/sum_weight;
    a = ave_deta2 - ave_deta*ave_deta;                          
    b = ave_dphi2 - ave_dphi*ave_dphi;                          
    c = -(sum_detadphi/sum_weight - ave_deta*ave_dphi);                
  } else ptD_NoQC_ = 0;
  float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
  if(a+b-delta > 0) axis2_NoQC_ = sqrt(0.5*(a+b-delta));
  else axis2_NoQC_ = 0.;

  mult_NoQC_ = (nChg_QC + nNeutral_ptCut);
}


void qgMiniTupleForMiniAOD::beginJob(){
  tree = fs->make<TTree>("qgMiniTupleForMiniAOD","qgMiniTuple");
  tree->Branch("nRun" ,		&nRun, 		"nRun/I");
  tree->Branch("nLumi" ,	&nLumi, 	"nLumi/I");
  tree->Branch("nEvent" ,	&nEvent, 	"nEvent/I");
  tree->Branch("rho" ,		&rho, 		"rho/F");
  tree->Branch("pt" ,		&pt,		"pt/F");
  tree->Branch("eta",		&eta,		"eta/F");
  tree->Branch("axis2_NoQC",	&axis2_NoQC,	"axis2_NoQC/F");
  tree->Branch("laxis2_NoQC",	&laxis2_NoQC,	"laxis2_NoQC/F");
  tree->Branch("ptD_NoQC",	&ptD_NoQC,	"ptD_NoQC/F");
  tree->Branch("mult_NoQC",	&mult_NoQC,	"mult_NoQC/I");
  tree->Branch("bTag",		&bTag,		"bTag/F");
  tree->Branch("partonId",	&partonId,	"partonId/I");
}


void qgMiniTupleForMiniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
  edm::ParameterSetDescription desc;
  desc.addUntracked<std::string>("fileName","qgMiniTupleForMiniAOD.root");
  desc.add<edm::InputTag>("rhoInputTag");
  desc.add<edm::InputTag>("jetsInputTag");
  desc.addUntracked<double>("minJetPt", 10.);
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(qgMiniTupleForMiniAOD);

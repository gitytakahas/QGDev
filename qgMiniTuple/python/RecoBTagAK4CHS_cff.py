import FWCore.ParameterSet.Config as cms
from RecoJets.JetAssociationProducers.ak4JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

ak4CHSJetTracksAssociator                    = ak4JetTracksAssociatorAtVertexPF.clone(jets = "ak4PFJetsCHS")
ak4CHSImpactParameterTagInfos                = impactParameterTagInfos.clone(jetTracks = 'ak4CHSJetTracksAssociator')
ak4CHSInclusiveSecondaryVertexFinderTagInfos = inclusiveSecondaryVertexFinderTagInfos.clone(trackIPTagInfos = 'ak4CHSImpactParameterTagInfos')

ak4CHSCombinedInclusiveSecondaryVertexV2BJetTags = combinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("ak4CHSImpactParameterTagInfos"),
			     cms.InputTag("ak4CHSInclusiveSecondaryVertexFinderTagInfos"))
)


ak4CHSBTagging = cms.Sequence(ak4CHSJetTracksAssociator *
                              ak4CHSImpactParameterTagInfos *
                              ak4CHSInclusiveSecondaryVertexFinderTagInfos *
                              ak4CHSCombinedInclusiveSecondaryVertexV2BJetTags
)

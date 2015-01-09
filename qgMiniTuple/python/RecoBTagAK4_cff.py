import FWCore.ParameterSet.Config as cms
from RecoJets.JetAssociationProducers.ak4JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

ak4JetTracksAssociator                    = ak4JetTracksAssociatorAtVertexPF.clone(jets = "ak4PFJets")
ak4ImpactParameterTagInfos                = impactParameterTagInfos.clone(jetTracks = 'ak4JetTracksAssociator')
ak4InclusiveSecondaryVertexFinderTagInfos = inclusiveSecondaryVertexFinderTagInfos.clone(trackIPTagInfos = 'ak4ImpactParameterTagInfos')

ak4CombinedInclusiveSecondaryVertexV2BJetTags = combinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("ak4ImpactParameterTagInfos"),
			     cms.InputTag("ak4InclusiveSecondaryVertexFinderTagInfos"))
)


ak4BTagging = cms.Sequence(ak4JetTracksAssociator * 
                           ak4ImpactParameterTagInfos *
                           ak4InclusiveSecondaryVertexFinderTagInfos *
                           ak4CombinedInclusiveSecondaryVertexV2BJetTags
)

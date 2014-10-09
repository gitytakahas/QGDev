import FWCore.ParameterSet.Config as cms
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

ak7JetTracksAssociatorAtVertex 		= ak5JetTracksAssociatorAtVertexPF.clone(jets = "ak7PFJets")
ak7ImpactParameterTagInfos 		= impactParameterTagInfos.clone(jetTracks = 'ak7JetTracksAssociatorAtVertex')
ak7TrackCountingHighEffBJetTags 	= trackCountingHighEffBJetTags.clone(tagInfos = cms.VInputTag('ak7ImpactParameterTagInfos'))
ak7TrackCountingHighPurBJetTags 	= trackCountingHighPurBJetTags.clone(tagInfos = cms.VInputTag('ak7ImpactParameterTagInfos'))
ak7JetProbabilityBJetTags 		= jetProbabilityBJetTags.clone(tagInfos = cms.VInputTag('ak7ImpactParameterTagInfos'))
ak7JetBProbabilityBJetTags 		= jetBProbabilityBJetTags.clone(tagInfos = cms.VInputTag('ak7ImpactParameterTagInfos'))

ak7SecondaryVertexTagInfos 		= secondaryVertexTagInfos.clone(trackIPTagInfos = 'ak7ImpactParameterTagInfos')
ak7SimpleSecondaryVertexBJetTags 	= simpleSecondaryVertexBJetTags.clone(tagInfos = cms.VInputTag('ak7SecondaryVertexTagInfos'))
ak7CombinedSecondaryVertexBJetTags 	= combinedSecondaryVertexBJetTags.clone(tagInfos = cms.VInputTag('ak7ImpactParameterTagInfos', 'ak7SecondaryVertexTagInfos'))
ak7CombinedSecondaryVertexMVABJetTags 	= combinedSecondaryVertexMVABJetTags.clone(tagInfos = cms.VInputTag('ak7ImpactParameterTagInfos', 'ak7SecondaryVertexTagInfos'))

ak7JetBtaggingIP = cms.Sequence(
    ak7ImpactParameterTagInfos * (
        ak7TrackCountingHighEffBJetTags +
        ak7TrackCountingHighPurBJetTags +
        ak7JetProbabilityBJetTags +
        ak7JetBProbabilityBJetTags
    )
)

ak7JetBtaggingSV = cms.Sequence(
    ak7ImpactParameterTagInfos *
    ak7SecondaryVertexTagInfos * (
        ak7SimpleSecondaryVertexBJetTags +
        ak7CombinedSecondaryVertexBJetTags +
        ak7CombinedSecondaryVertexMVABJetTags
    )
)

ak7JetBtagging = cms.Sequence(
    ak7JetBtaggingIP +
    ak7JetBtaggingSV
)

ak7BTagging = cms.Sequence(ak7JetTracksAssociatorAtVertex * ak7JetBtagging)

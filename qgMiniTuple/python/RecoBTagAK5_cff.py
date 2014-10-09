import FWCore.ParameterSet.Config as cms
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

ak5JetTracksAssociatorAtVertex2		= ak5JetTracksAssociatorAtVertexPF.clone(jets = "ak5PFJets")
ak5ImpactParameterTagInfos 		= impactParameterTagInfos.clone(jetTracks = 'ak5JetTracksAssociatorAtVertex2')
ak5TrackCountingHighEffBJetTags 	= trackCountingHighEffBJetTags.clone(tagInfos = cms.VInputTag('ak5ImpactParameterTagInfos'))
ak5TrackCountingHighPurBJetTags 	= trackCountingHighPurBJetTags.clone(tagInfos = cms.VInputTag('ak5ImpactParameterTagInfos'))
ak5JetProbabilityBJetTags 		= jetProbabilityBJetTags.clone(tagInfos = cms.VInputTag('ak5ImpactParameterTagInfos'))
ak5JetBProbabilityBJetTags 		= jetBProbabilityBJetTags.clone(tagInfos = cms.VInputTag('ak5ImpactParameterTagInfos'))

ak5SecondaryVertexTagInfos 		= secondaryVertexTagInfos.clone(trackIPTagInfos = 'ak5ImpactParameterTagInfos')
ak5SimpleSecondaryVertexBJetTags 	= simpleSecondaryVertexBJetTags.clone(tagInfos = cms.VInputTag('ak5SecondaryVertexTagInfos'))
ak5CombinedSecondaryVertexBJetTags 	= combinedSecondaryVertexBJetTags.clone(tagInfos = cms.VInputTag('ak5ImpactParameterTagInfos', 'ak5SecondaryVertexTagInfos'))
ak5CombinedSecondaryVertexMVABJetTags 	= combinedSecondaryVertexMVABJetTags.clone(tagInfos = cms.VInputTag('ak5ImpactParameterTagInfos', 'ak5SecondaryVertexTagInfos'))

ak5JetBtaggingIP = cms.Sequence(
    ak5ImpactParameterTagInfos * (
        ak5TrackCountingHighEffBJetTags +
        ak5TrackCountingHighPurBJetTags +
        ak5JetProbabilityBJetTags +
        ak5JetBProbabilityBJetTags
    )
)

ak5JetBtaggingSV = cms.Sequence(
    ak5ImpactParameterTagInfos *
    ak5SecondaryVertexTagInfos * (
        ak5SimpleSecondaryVertexBJetTags +
        ak5CombinedSecondaryVertexBJetTags +
        ak5CombinedSecondaryVertexMVABJetTags
    )
)

ak5JetBtagging = cms.Sequence(
    ak5JetBtaggingIP +
    ak5JetBtaggingSV
)

ak5BTagging = cms.Sequence(ak5JetTracksAssociatorAtVertex2 * ak5JetBtagging)

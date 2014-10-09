import FWCore.ParameterSet.Config as cms
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

ak4JetTracksAssociatorAtVertex 		= ak5JetTracksAssociatorAtVertexPF.clone(jets = "ak4PFJets")
ak4ImpactParameterTagInfos 		= impactParameterTagInfos.clone(jetTracks = 'ak4JetTracksAssociatorAtVertex')
ak4TrackCountingHighEffBJetTags 	= trackCountingHighEffBJetTags.clone(tagInfos = cms.VInputTag('ak4ImpactParameterTagInfos'))
ak4TrackCountingHighPurBJetTags 	= trackCountingHighPurBJetTags.clone(tagInfos = cms.VInputTag('ak4ImpactParameterTagInfos'))
ak4JetProbabilityBJetTags 		= jetProbabilityBJetTags.clone(tagInfos = cms.VInputTag('ak4ImpactParameterTagInfos'))
ak4JetBProbabilityBJetTags 		= jetBProbabilityBJetTags.clone(tagInfos = cms.VInputTag('ak4ImpactParameterTagInfos'))

ak4SecondaryVertexTagInfos 		= secondaryVertexTagInfos.clone(trackIPTagInfos = 'ak4ImpactParameterTagInfos')
ak4SimpleSecondaryVertexBJetTags 	= simpleSecondaryVertexBJetTags.clone(tagInfos = cms.VInputTag('ak4SecondaryVertexTagInfos'))
ak4CombinedSecondaryVertexBJetTags 	= combinedSecondaryVertexBJetTags.clone(tagInfos = cms.VInputTag('ak4ImpactParameterTagInfos', 'ak4SecondaryVertexTagInfos'))
ak4CombinedSecondaryVertexMVABJetTags 	= combinedSecondaryVertexMVABJetTags.clone(tagInfos = cms.VInputTag('ak4ImpactParameterTagInfos', 'ak4SecondaryVertexTagInfos'))

ak4JetBtaggingIP = cms.Sequence(
    ak4ImpactParameterTagInfos * (
        ak4TrackCountingHighEffBJetTags +
        ak4TrackCountingHighPurBJetTags +
        ak4JetProbabilityBJetTags +
        ak4JetBProbabilityBJetTags
    )
)

ak4JetBtaggingSV = cms.Sequence(
    ak4ImpactParameterTagInfos *
    ak4SecondaryVertexTagInfos * (
        ak4SimpleSecondaryVertexBJetTags +
        ak4CombinedSecondaryVertexBJetTags +
        ak4CombinedSecondaryVertexMVABJetTags
    )
)

ak4JetBtagging = cms.Sequence(
    ak4JetBtaggingIP +
    ak4JetBtaggingSV
)

ak4BTagging = cms.Sequence(ak4JetTracksAssociatorAtVertex * ak4JetBtagging)

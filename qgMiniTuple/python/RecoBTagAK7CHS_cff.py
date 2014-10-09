import FWCore.ParameterSet.Config as cms
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

ak7CHSJetTracksAssociatorAtVertex 		= ak5JetTracksAssociatorAtVertexPF.clone(jets = "ak7PFJetsCHS")
ak7CHSImpactParameterTagInfos 			= impactParameterTagInfos.clone(jetTracks = 'ak7CHSJetTracksAssociatorAtVertex')
ak7CHSTrackCountingHighEffBJetTags 		= trackCountingHighEffBJetTags.clone(tagInfos = cms.VInputTag('ak7CHSImpactParameterTagInfos'))
ak7CHSTrackCountingHighPurBJetTags 		= trackCountingHighPurBJetTags.clone(tagInfos = cms.VInputTag('ak7CHSImpactParameterTagInfos'))
ak7CHSJetProbabilityBJetTags 			= jetProbabilityBJetTags.clone(tagInfos = cms.VInputTag('ak7CHSImpactParameterTagInfos'))
ak7CHSJetBProbabilityBJetTags 			= jetBProbabilityBJetTags.clone(tagInfos = cms.VInputTag('ak7CHSImpactParameterTagInfos'))

ak7CHSSecondaryVertexTagInfos 			= secondaryVertexTagInfos.clone(trackIPTagInfos = 'ak7CHSImpactParameterTagInfos')
ak7CHSSimpleSecondaryVertexBJetTags 		= simpleSecondaryVertexBJetTags.clone(tagInfos = cms.VInputTag('ak7CHSSecondaryVertexTagInfos'))
ak7CHSCombinedSecondaryVertexBJetTags 		= combinedSecondaryVertexBJetTags.clone(tagInfos = cms.VInputTag('ak7CHSImpactParameterTagInfos', 'ak7CHSSecondaryVertexTagInfos'))
ak7CHSCombinedSecondaryVertexMVABJetTags 	= combinedSecondaryVertexMVABJetTags.clone(tagInfos = cms.VInputTag('ak7CHSImpactParameterTagInfos', 'ak7CHSSecondaryVertexTagInfos'))

ak7CHSJetBtaggingIP = cms.Sequence(
    ak7CHSImpactParameterTagInfos * (
        ak7CHSTrackCountingHighEffBJetTags +
        ak7CHSTrackCountingHighPurBJetTags +
        ak7CHSJetProbabilityBJetTags +
        ak7CHSJetBProbabilityBJetTags
    )
)

ak7CHSJetBtaggingSV = cms.Sequence(
    ak7CHSImpactParameterTagInfos *
    ak7CHSSecondaryVertexTagInfos * (
        ak7CHSSimpleSecondaryVertexBJetTags +
        ak7CHSCombinedSecondaryVertexBJetTags +
        ak7CHSCombinedSecondaryVertexMVABJetTags
    )
)

ak7CHSJetBtagging = cms.Sequence(
    ak7CHSJetBtaggingIP +
    ak7CHSJetBtaggingSV
)

ak7CHSBTagging = cms.Sequence(ak7CHSJetTracksAssociatorAtVertex * ak7CHSJetBtagging)

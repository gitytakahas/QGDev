import FWCore.ParameterSet.Config as cms
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

ak5CHSJetTracksAssociatorAtVertex 		= ak5JetTracksAssociatorAtVertexPF.clone(jets = "ak5PFJetsCHS")
ak5CHSImpactParameterTagInfos 			= impactParameterTagInfos.clone(jetTracks = 'ak5CHSJetTracksAssociatorAtVertex')
ak5CHSTrackCountingHighEffBJetTags 		= trackCountingHighEffBJetTags.clone(tagInfos = cms.VInputTag('ak5CHSImpactParameterTagInfos'))
ak5CHSTrackCountingHighPurBJetTags 		= trackCountingHighPurBJetTags.clone(tagInfos = cms.VInputTag('ak5CHSImpactParameterTagInfos'))
ak5CHSJetProbabilityBJetTags 			= jetProbabilityBJetTags.clone(tagInfos = cms.VInputTag('ak5CHSImpactParameterTagInfos'))
ak5CHSJetBProbabilityBJetTags 			= jetBProbabilityBJetTags.clone(tagInfos = cms.VInputTag('ak5CHSImpactParameterTagInfos'))

ak5CHSSecondaryVertexTagInfos 			= secondaryVertexTagInfos.clone(trackIPTagInfos = 'ak5CHSImpactParameterTagInfos')
ak5CHSSimpleSecondaryVertexBJetTags 		= simpleSecondaryVertexBJetTags.clone(tagInfos = cms.VInputTag('ak5CHSSecondaryVertexTagInfos'))
ak5CHSCombinedSecondaryVertexBJetTags 		= combinedSecondaryVertexBJetTags.clone(tagInfos = cms.VInputTag('ak5CHSImpactParameterTagInfos', 'ak5CHSSecondaryVertexTagInfos'))
ak5CHSCombinedSecondaryVertexMVABJetTags 	= combinedSecondaryVertexMVABJetTags.clone(tagInfos = cms.VInputTag('ak5CHSImpactParameterTagInfos', 'ak5CHSSecondaryVertexTagInfos'))

ak5CHSJetBtaggingIP = cms.Sequence(
    ak5CHSImpactParameterTagInfos * (
        ak5CHSTrackCountingHighEffBJetTags +
        ak5CHSTrackCountingHighPurBJetTags +
        ak5CHSJetProbabilityBJetTags +
        ak5CHSJetBProbabilityBJetTags
    )
)

ak5CHSJetBtaggingSV = cms.Sequence(
    ak5CHSImpactParameterTagInfos *
    ak5CHSSecondaryVertexTagInfos * (
        ak5CHSSimpleSecondaryVertexBJetTags +
        ak5CHSCombinedSecondaryVertexBJetTags +
        ak5CHSCombinedSecondaryVertexMVABJetTags
    )
)

ak5CHSJetBtagging = cms.Sequence(
    ak5CHSJetBtaggingIP +
    ak5CHSJetBtaggingSV
)

ak5CHSBTagging = cms.Sequence(ak5CHSJetTracksAssociatorAtVertex * ak5CHSJetBtagging)

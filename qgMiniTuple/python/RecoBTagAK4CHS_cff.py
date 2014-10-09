import FWCore.ParameterSet.Config as cms
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

ak4CHSJetTracksAssociatorAtVertex 		= ak5JetTracksAssociatorAtVertexPF.clone(jets = "ak4PFJetsCHS")
ak4CHSImpactParameterTagInfos 			= impactParameterTagInfos.clone(jetTracks = 'ak4CHSJetTracksAssociatorAtVertex')
ak4CHSTrackCountingHighEffBJetTags 		= trackCountingHighEffBJetTags.clone(tagInfos = cms.VInputTag('ak4CHSImpactParameterTagInfos'))
ak4CHSTrackCountingHighPurBJetTags 		= trackCountingHighPurBJetTags.clone(tagInfos = cms.VInputTag('ak4CHSImpactParameterTagInfos'))
ak4CHSJetProbabilityBJetTags 			= jetProbabilityBJetTags.clone(tagInfos = cms.VInputTag('ak4CHSImpactParameterTagInfos'))
ak4CHSJetBProbabilityBJetTags 			= jetBProbabilityBJetTags.clone(tagInfos = cms.VInputTag('ak4CHSImpactParameterTagInfos'))

ak4CHSSecondaryVertexTagInfos 			= secondaryVertexTagInfos.clone(trackIPTagInfos = 'ak4CHSImpactParameterTagInfos')
ak4CHSSimpleSecondaryVertexBJetTags 		= simpleSecondaryVertexBJetTags.clone(tagInfos = cms.VInputTag('ak4CHSSecondaryVertexTagInfos'))
ak4CHSCombinedSecondaryVertexBJetTags 		= combinedSecondaryVertexBJetTags.clone(tagInfos = cms.VInputTag('ak4CHSImpactParameterTagInfos', 'ak4CHSSecondaryVertexTagInfos'))
ak4CHSCombinedSecondaryVertexMVABJetTags 	= combinedSecondaryVertexMVABJetTags.clone(tagInfos = cms.VInputTag('ak4CHSImpactParameterTagInfos', 'ak4CHSSecondaryVertexTagInfos'))

ak4CHSJetBtaggingIP = cms.Sequence(
    ak4CHSImpactParameterTagInfos * (
        ak4CHSTrackCountingHighEffBJetTags +
        ak4CHSTrackCountingHighPurBJetTags +
        ak4CHSJetProbabilityBJetTags +
        ak4CHSJetBProbabilityBJetTags
    )
)

ak4CHSJetBtaggingSV = cms.Sequence(
    ak4CHSImpactParameterTagInfos *
    ak4CHSSecondaryVertexTagInfos * (
        ak4CHSSimpleSecondaryVertexBJetTags +
        ak4CHSCombinedSecondaryVertexBJetTags +
        ak4CHSCombinedSecondaryVertexMVABJetTags
    )
)

ak4CHSJetBtagging = cms.Sequence(
    ak4CHSJetBtaggingIP +
    ak4CHSJetBtaggingSV
)

ak4CHSBTagging = cms.Sequence(ak4CHSJetTracksAssociatorAtVertex * ak4CHSJetBtagging)

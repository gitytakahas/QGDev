import FWCore.ParameterSet.Config as cms
from RecoBTag.Configuration.RecoBTag_cff import *

ak4CHSImpactParameterTagInfos                = pfImpactParameterTagInfos.clone(jets = 'ak4PFJetsCHS')
ak4CHSInclusiveSecondaryVertexFinderTagInfos = pfInclusiveSecondaryVertexFinderTagInfos.clone(trackIPTagInfos = 'ak4CHSImpactParameterTagInfos')

ak4CHSCombinedInclusiveSecondaryVertexV2BJetTags = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("ak4CHSImpactParameterTagInfos"),
			     cms.InputTag("ak4CHSInclusiveSecondaryVertexFinderTagInfos"))
)


ak4CHSBTagging = cms.Sequence(ak4CHSImpactParameterTagInfos *
                              ak4CHSInclusiveSecondaryVertexFinderTagInfos *
                              ak4CHSCombinedInclusiveSecondaryVertexV2BJetTags
)

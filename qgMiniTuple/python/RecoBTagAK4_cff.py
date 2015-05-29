import FWCore.ParameterSet.Config as cms
from RecoBTag.Configuration.RecoBTag_cff import *

ak4ImpactParameterTagInfos                = pfImpactParameterTagInfos.clone(jets = 'ak4PFJets')
ak4InclusiveSecondaryVertexFinderTagInfos = pfInclusiveSecondaryVertexFinderTagInfos.clone(trackIPTagInfos = 'ak4ImpactParameterTagInfos')

ak4CombinedInclusiveSecondaryVertexV2BJetTags = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("ak4ImpactParameterTagInfos"),
			     cms.InputTag("ak4InclusiveSecondaryVertexFinderTagInfos"))
)


ak4BTagging = cms.Sequence(ak4ImpactParameterTagInfos *
                           ak4InclusiveSecondaryVertexFinderTagInfos *
                           ak4CombinedInclusiveSecondaryVertexV2BJetTags
)

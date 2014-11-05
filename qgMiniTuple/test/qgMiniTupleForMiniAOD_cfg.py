import FWCore.ParameterSet.Config as cms

process = cms.Process("qgMiniTupleProducer")

# Settings for local tests
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('file:/user/tomc/public/TTJets_forSynch_1.root')
)

# Standard configurations
process.load('Configuration.StandardSequences.Services_cff')

# Use TFileService to put trees from different analyzers in one file
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("qgMiniTupleForMiniAOD.root"),
    closeFileFast = cms.untracked.bool(True)
)


process.qgMiniTuple = cms.EDAnalyzer("qgMiniTuple",
    usePatJets			= cms.untracked.bool(True),
    fileName 			= cms.untracked.string('qgMiniTuple.root'),
    rhoInputTag			= cms.InputTag('fixedGridRhoFastjetAll'),
    csvInputTag			= cms.InputTag('combinedSecondaryVertexBJetTags'),
    vertexInputTag		= cms.InputTag('offlineSlimmedPrimaryVertices'),
    jetsInputTag		= cms.InputTag('slimmedJets'),
    genJetsInputTag		= cms.InputTag('slimmedGenJets'),
    genParticlesInputTag	= cms.InputTag('prunedGenParticles'),
    minJetPt			= cms.untracked.double(20.),
    deltaRcut			= cms.untracked.double(0.3),
    jec				= cms.string(''),						# Ignored when using pat mode
    jetFlavourInputTag		= cms.InputTag(''),						# Ignored when using pat mode
)


process.p = cms.Path(process.qgMiniTuple)

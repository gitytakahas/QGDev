import FWCore.ParameterSet.Config as cms

process = cms.Process("qgMiniTupleProducer")

# Settings for local tests
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('file:/user/tomc/public/TTJets_forSynch_1.root')
)

# Standard configurations
process.load('Configuration.StandardSequences.Services_cff')

# Use TFileService to put trees from different analyzers in one file
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("qgMiniTuple.root"),
    closeFileFast = cms.untracked.bool(True)
)


process.qgMiniTupleMiniAOD = cms.EDAnalyzer("qgMiniTupleForMiniAOD",
    fileName 			= cms.untracked.string('qgMiniTuple.root'),
    rhoInputTag			= cms.InputTag('fixedGridRhoFastjetAll'),
    jetsInputTag		= cms.InputTag('slimmedJets'),
    minJetPt			= cms.untracked.double(20.),
)

process.p = cms.Path(process.qgMiniTupleMiniAOD)

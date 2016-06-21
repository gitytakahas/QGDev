import FWCore.ParameterSet.Config as cms
import glob

process = cms.Process("qgMiniTupleProducer")

# Settings for local tests
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))

fileList=[""]

### do not remove the line below!
###FILELIST###

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(fileList) # Test file available at T2B
)

# Standard configurations
process.load('Configuration.StandardSequences.Services_cff')

# Use TFileService to put trees from different analyzers in one file
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("qgMiniTuple.root"),
    closeFileFast = cms.untracked.bool(True)
)


process.qgMiniTupleAK4chs = cms.EDAnalyzer("qgMiniTuple",
    rhoInputTag			= cms.InputTag('fixedGridRhoFastjetAll'),
    csvInputTag			= cms.InputTag('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
    vertexInputTag		= cms.InputTag('offlineSlimmedPrimaryVertices'),
    jetsInputTag		= cms.InputTag('slimmedJets'),
    genJetsInputTag		= cms.InputTag('slimmedGenJets'),
    genParticlesInputTag	= cms.InputTag('prunedGenParticles'),
    pfCandidatesInputTag	= cms.InputTag('packedPFCandidates'),
    minJetPt			= cms.untracked.double(20.),
    deltaRcut			= cms.untracked.double(0.3),
    jec				= cms.string(''),						# Ignored when using pat mode
)


process.p = cms.Path(process.qgMiniTupleAK4chs)

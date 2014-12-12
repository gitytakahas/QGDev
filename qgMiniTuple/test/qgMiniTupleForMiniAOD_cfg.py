import FWCore.ParameterSet.Config as cms
import glob

process = cms.Process("qgMiniTupleProducer")

files = glob.glob("/pnfs/iihe/cms/ph/sc4/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/*")
files = ["/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/" + file.split('/')[-1] for file in files]

# Settings for local tests
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.source = cms.Source("PoolSource", 
#    fileNames = cms.untracked.vstring('file:/user/tomc/public/TTJets_forSynch_1.root')
    fileNames = cms.untracked.vstring(files),
    skipEvents=cms.untracked.uint32(25000000)
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
    pfCandidatesInputTag	= cms.InputTag('packedPFCandidates'),
    minJetPt			= cms.untracked.double(20.),
    deltaRcut			= cms.untracked.double(0.3),
    jec				= cms.string(''),						# Ignored when using pat mode
)


process.p = cms.Path(process.qgMiniTuple)

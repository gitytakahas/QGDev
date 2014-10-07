import FWCore.ParameterSet.Config as cms

process = cms.Process("qgMiniTupleProducer")

# Settings for local tests
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('/store/mc/Spring14miniaod/QCD_Pt-30to50_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/02EF5246-7F06-E411-A4A2-002590596486.root') #Test file available at T2B.
)

# Standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START71_V8A::All'

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

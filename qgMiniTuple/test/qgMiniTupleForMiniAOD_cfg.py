import FWCore.ParameterSet.Config as cms
import glob

process = cms.Process("qgMiniTupleProducer")

# Settings for local tests
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))
process.source = cms.Source("PoolSource", 
#    fileNames = cms.untracked.vstring('/store/mc/Phys14DR/QCD_Pt-50to80_Tune4C_13TeV_pythia8/MINIAODSIM/AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/00000/0C31D9A0-548B-E411-ABBB-0025905A60F4.root') # Test file available at T2B
#    fileNames = cms.untracked.vstring('/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root') # Test file available at T2B
    fileNames = cms.untracked.vstring('/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/0066F143-F8FD-E411-9A0B-D4AE526A0D2E.root') # Test file available at T2B
#    fileNames = cms.untracked.vstring('file:miniAOD.root') # Test file available at T2B
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

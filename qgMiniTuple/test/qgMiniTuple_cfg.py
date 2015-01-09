import FWCore.ParameterSet.Config as cms

process = cms.Process("qgMiniTupleProducer")

# Settings for local tests
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('/store/mc/Phys14DR/QCD_Pt-50to80_Tune4C_13TeV_pythia8/AODSIM/AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/00000/00E7A97D-118B-E411-9CA8-002618943920.root') #Test file available at T2B (PHYS14DR)
)

# Standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.GlobalTag.globaltag = 'PHYS14_50_V1::All'

# Use TFileService to put trees from different analyzers in one file
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("qgMiniTuple.root"),
    closeFileFast = cms.untracked.bool(True)
)


# Sequence for making all jet collections
# To do: add puppi jets when CMSSW recipe is available
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.myRecoPFJets = cms.Sequence(process.fixedGridRhoFastjetAll + process.pfNoPileUpJMESequence + process.genParticlesForJets + 
                                    process.ak4GenJets + process.ak4PFJets + process.ak4PFJetsCHS
)

# Jet energy corrections and b-tagging
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('QGDev.qgMiniTuple.RecoBTagAK4_cff')
process.load('QGDev.qgMiniTuple.RecoBTagAK4CHS_cff')

# Use clones of the QGTagger process for the different jet collections available [CMSSW_7_4_X and higher]
# process.load('RecoJets.JetProducers.QGTagger_cfi')
# process.QGTaggerAK4 		= process.QGTagger.clone(srcJets = 'ak4PFJets',		jetsLabel = 'QGL_AK4PF')
# process.QGTaggerAK4chs 	= process.QGTagger.clone(srcJets = 'ak4PFJetsCHS',	jetsLabel = 'QGL_AK4PFchs')

# Init qgMiniTuple analyzer for AK4, make clones for the other ones
process.qgMiniTupleAK4 = cms.EDAnalyzer("qgMiniTuple",
    rhoInputTag			= cms.InputTag('fixedGridRhoFastjetAll'),
    vertexInputTag		= cms.InputTag('offlinePrimaryVerticesWithBS'),
    jetsInputTag		= cms.InputTag('ak4PFJets'),
    pfCandidatesInputTag	= cms.InputTag('particleFlow'),
    genJetsInputTag		= cms.InputTag('ak4GenJets'),
    jec				= cms.string('ak4PFL1FastL2L3'),
    genParticlesInputTag	= cms.InputTag('genParticles'),
    csvInputTag			= cms.InputTag('ak4CombinedInclusiveSecondaryVertexV2BJetTags'),
#   qgVariablesInputTag		= cms.InputTag('QGTaggerAK4'),
    minJetPt			= cms.untracked.double(20.),
    deltaRcut			= cms.untracked.double(0.3),
)
process.qgMiniTupleAK4chs 	= process.qgMiniTupleAK4.clone(jetsInputTag = 'ak4PFJetsCHS', jec = 'ak4PFCHSL1FastL2L3', csvInputTag = 'ak4CHSCombinedInclusiveSecondaryVertexV2BJetTags')#, qgVariablesInputTag = 'QGTaggerAK4chs')

# The path: jet sequence + b-tagging + (QGTagger) + qgMiniTuple for every jet collection
process.p = cms.Path(process.myRecoPFJets *
                     process.ak4BTagging    * process.ak4CHSBTagging *
                     process.qgMiniTupleAK4 * process.qgMiniTupleAK4chs
)

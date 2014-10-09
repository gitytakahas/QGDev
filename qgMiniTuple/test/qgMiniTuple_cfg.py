import FWCore.ParameterSet.Config as cms

process = cms.Process("qgMiniTupleProducer")

# Settings for local tests
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('/store/mc/Fall13dr/EWKZjj_mqq120_mll50_13TeV_madgraph-pythia8/AODSIM/tsg_PU40bx50_POSTLS162_V2-v1/00000/16C169C2-C475-E311-8BEE-7845C4FC3A2B.root') #Test file available at T2B.
)

# Standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.GlobalTag.globaltag = 'PLS170_V6AN1::All'

# Use TFileService to put trees from different analyzers in one file
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("qgMiniTuple.root"),
    closeFileFast = cms.untracked.bool(True)
)


# Sequence for making all jet collections
# To do: add puppi jets when CMSSW recipe is available
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.ak7PFJetsCHS = process.ak4PFJetsCHS.clone(rParam = 0.7) 
process.myRecoPFJets = cms.Sequence(process.fixedGridRhoAll +
                                    process.fixedGridRhoFastjetAll +
                                    process.fixedGridRhoFastjetCentralChargedPileUp +
                                    process.fixedGridRhoFastjetCentralNeutral +
                                    process.pfNoPileUpJMESequence +
                                    process.ak4PFJets + process.ak4PFJetsCHS + 
                                    process.ak5PFJets + process.ak5PFJetsCHS +
                                    process.ak7PFJets + process.ak7PFJetsCHS
)


# Jet energy corrections [do not waste time on these corrections, they are not available with the PLS170_V6AN1 global tag]
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.ak4PFJetsL1FastL2L3 	= process.ak5PFJetsL2L3.clone(src = 'ak4PFJets',    correctors = ['ak4PFL1FastL2L3'])
process.ak7PFJetsL1FastL2L3    	= process.ak5PFJetsL2L3.clone(src = 'ak7PFJets',    correctors = ['ak7PFL1FastL2L3'])
process.ak4PFchsJetsL1FastL2L3 	= process.ak5PFJetsL2L3.clone(src = 'ak4PFJetsCHS', correctors = ['ak4PFCHSL1FastL2L3'])
process.ak5PFchsJetsL1FastL2L3 	= process.ak5PFJetsL2L3.clone(src = 'ak5PFJetsCHS', correctors = ['ak5PFCHSL1FastL2L3'])
process.ak7PFchsJetsL1FastL2L3 	= process.ak5PFJetsL2L3.clone(src = 'ak7PFJetsCHS', correctors = ['ak7PFCHSL1FastL2L3'])


# b-tagging
process.load('QGDev.qgMiniTuple.RecoBTagAK4_cff')
process.load('QGDev.qgMiniTuple.RecoBTagAK5_cff')
process.load('QGDev.qgMiniTuple.RecoBTagAK7_cff')
process.load('QGDev.qgMiniTuple.RecoBTagAK4CHS_cff')
process.load('QGDev.qgMiniTuple.RecoBTagAK5CHS_cff')
process.load('QGDev.qgMiniTuple.RecoBTagAK7CHS_cff')

# Flavour matching
process.selectedHadronsAndPartons = cms.EDProducer('HadronAndPartonSelector',
    src = cms.InputTag("generator"),
    particles = cms.InputTag("genParticles"),
    partonMode = cms.string("Auto")
)
process.jetFlavourInfosAK4 = cms.EDProducer("JetFlavourClustering",
    jets = cms.InputTag("ak4PFJets"),
    bHadrons = cms.InputTag("selectedHadronsAndPartons","bHadrons"),
    cHadrons = cms.InputTag("selectedHadronsAndPartons","cHadrons"),
    partons = cms.InputTag("selectedHadronsAndPartons","partons"),
    jetAlgorithm = cms.string("AntiKt"),
    rParam = cms.double(0.4),
    ghostRescaling = cms.double(1e-18),
    hadronFlavourHasPriority = cms.bool(True)
)
process.jetFlavourInfosAK5	= process.jetFlavourInfosAK4.clone(jets = 'ak5PFJets', rParam = 0.5)
process.jetFlavourInfosAK7	= process.jetFlavourInfosAK4.clone(jets = 'ak7PFJets', rParam = 0.7)
process.jetFlavourInfosAK4CHS	= process.jetFlavourInfosAK4.clone(jets = 'ak4PFJetsCHS')
process.jetFlavourInfosAK5CHS	= process.jetFlavourInfosAK4.clone(jets = 'ak5PFJetsCHS', rParam = 0.5)
process.jetFlavourInfosAK7CHS	= process.jetFlavourInfosAK4.clone(jets = 'ak7PFJetsCHS', rParam = 0.7)


# Use clones of the QGTagger process for the different jet collections available [CMSSW_7_1_X and higher]
# process.load('QGDev.qgMiniTuple.QGTagger_v1-preliminary_cfi') # v1-preliminary includes: AK4,AK4chs,AK5,AK5chs
# process.QGTaggerAK4 		= process.QGTagger.clone(srcJets = 'ak4PFJets',		jetsLabel = 'QGL_AK4PF')
# process.QGTaggerAK4chs 	= process.QGTagger.clone(srcJets = 'ak4PFJetsCHS',	jetsLabel = 'QGL_AK4PFchs')

# Init qgMiniTuple analyzer for AK4, make clones for the other ones
process.qgMiniTupleAK4 = cms.EDAnalyzer("qgMiniTuple",
    fileName 			= cms.untracked.string('qgMiniTuple.root'),
    rhoInputTag			= cms.InputTag('fixedGridRhoFastjetAll'),
    vertexInputTag		= cms.InputTag('offlinePrimaryVerticesWithBS'),
    jetsInputTag		= cms.InputTag('ak4PFJets'),
    jetFlavourInputTag		= cms.InputTag('jetFlavourInfosAK4'),
    csvInputTag			= cms.InputTag('ak4CombinedSecondaryVertexBJetTags'),
#   qgVariablesInputTag		= cms.InputTag('QGTaggerAK4'),
    minJetPt			= cms.untracked.double(20.),
)
process.qgMiniTupleAK5 		= process.qgMiniTupleAK4.clone(jetsInputTag = 'ak5PFJets',	jetFlavourInputTag = 'jetFlavourInfosAK5',	csvInputTag = 'ak5CombinedSecondaryVertexBJetTags')#,    qgVariablesInputTag = 'QGTaggerAK5')
process.qgMiniTupleAK7 		= process.qgMiniTupleAK4.clone(jetsInputTag = 'ak7PFJets',	jetFlavourInputTag = 'jetFlavourInfosAK7',	csvInputTag = 'ak7CombinedSecondaryVertexBJetTags')#,    qgVariablesInputTag = 'QGTaggerAK7')
process.qgMiniTupleAK4chs 	= process.qgMiniTupleAK4.clone(jetsInputTag = 'ak4PFJetsCHS',	jetFlavourInputTag = 'jetFlavourInfosAK4CHS',	csvInputTag = 'ak4CHSCombinedSecondaryVertexBJetTags')#, qgVariablesInputTag = 'QGTaggerAK4chs')
process.qgMiniTupleAK5chs 	= process.qgMiniTupleAK4.clone(jetsInputTag = 'ak5PFJetsCHS',	jetFlavourInputTag = 'jetFlavourInfosAK5CHS',	csvInputTag = 'ak5CHSCombinedSecondaryVertexBJetTags')#, qgVariablesInputTag = 'QGTaggerAK5chs')
process.qgMiniTupleAK7chs 	= process.qgMiniTupleAK4.clone(jetsInputTag = 'ak7PFJetsCHS',	jetFlavourInputTag = 'jetFlavourInfosAK7CHS',	csvInputTag = 'ak7CHSCombinedSecondaryVertexBJetTags')#, qgVariablesInputTag = 'QGTaggerAK7chs')

#process.qgMiniTupleAK4chs.bTagUsed = cms.untracked.bool(True)


# The path: jet sequence + b tagging + (QGTagger + qgMiniTuple) for every jet collection
process.p = cms.Path(process.myRecoPFJets *
#                     process.ak4PFJetsL1FastL2L3 * process.ak4PFchsJetsL1FastL2L3 *
#                     process.ak5PFJetsL1FastL2L3 * process.ak5PFchsJetsL1FastL2L3 *
#                     process.ak7PFJetsL1FastL2L3 * process.ak7PFJetsL1FastL2L3 *
                     process.ak4BTagging * process.ak4CHSBTagging *
#                     process.ak5BTagging * process.ak5CHSBTagging *
#                     process.ak7BTagging * process.ak7CHSBTagging *
                     process.selectedHadronsAndPartons *
                     process.jetFlavourInfosAK4 * process.jetFlavourInfosAK4CHS *
#                     process.jetFlavourInfosAK5 * process.jetFlavourInfosAK5CHS *
#                     process.jetFlavourInfosAK7 * process.jetFlavourInfosAK7CHS *
                     process.qgMiniTupleAK4 * process.qgMiniTupleAK4chs *
#                     process.qgMiniTupleAK5 * process.qgMiniTupleAK5chs *
#                     process.qgMiniTupleAK7 * process.qgMiniTupleAK7chs)


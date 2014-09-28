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
process.GlobalTag.globaltag = 'START71_V8A::All'

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


# Jet energy corrections [do not waste time on these corrections, payloads are not available through the database, it probably gets fixed in a future release]
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.ak5PFJetsL1FastL2L3 	= process.ak4PFJetsL2L3.clone(src = 'ak5PFJets',    correctors = ['ak5PFL1FastJet',   'ak5PFL2Relative',   'ak5PFL3Absolute'])
process.ak7PFJetsL1FastL2L3    	= process.ak4PFJetsL2L3.clone(src = 'ak7PFJets',    correctors = ['ak7PFL1FastJet',   'ak7PFL2Relative',   'ak7PFL3Absolute'])
process.ak4PFchsJetsL1FastL2L3 	= process.ak4PFJetsL2L3.clone(src = 'ak4PFchsJets', correctors = ['ak4PFchsL1FastJet','ak4PFchsL2Relative','ak4PFchsL3Absolute'])
process.ak5PFchsJetsL1FastL2L3 	= process.ak4PFJetsL2L3.clone(src = 'ak5PFchsJets', correctors = ['ak5PFchsL1FastJet','ak5PFchsL2Relative','ak5PFchsL3Absolute'])
process.ak7PFchsJetsL1FastL2L3 	= process.ak4PFJetsL2L3.clone(src = 'ak7PFchsJets', correctors = ['ak7PFchsL1FastJet','ak7PFchsL2Relative','ak7PFchsL3Absolute'])


# b-tagging
process.load("RecoJets.JetAssociationProducers.ak4JTA_cff")
process.load('RecoBTag.Configuration.RecoBTag_cff')
process.load('RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff')
process.load('RecoBTag.SecondaryVertex.combinedSecondaryVertexES_cfi')
process.load('RecoBTag.SecondaryVertex.combinedSecondaryVertexBJetTags_cfi')
process.combinedSecondaryVertexV2.calibrationRecords = cms.vstring(
   'CombinedSVIVFV2RecoVertex',
   'CombinedSVIVFV2PseudoVertex',
   'CombinedSVIVFV2NoVertex'
)
process.combinedSecondaryVertexV2BJetTags.tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfos"), cms.InputTag("inclusiveSecondaryVertexFinderTagInfos"))
process.CSVIVFV2ak4chs = cms.Sequence(process.ak4JetTracksAssociatorAtVertexPF * process.btagging *
                                      process.inclusiveVertexing * process.inclusiveSecondaryVertexFinderTagInfos * process.combinedSecondaryVertexV2BJetTags)

# Use clones of the QGTagger process for the different jet collections available
process.load('QGDev.qgMiniTuple.QGTagger_v1-preliminary_cfi') # v1-preliminary includes: AK4,AK4chs,AK5,AK5chs
process.QGTaggerAK4 		= process.QGTagger.clone(srcJets = 'ak4PFJets',		jetsLabel = 'QGL_AK4PF')
process.QGTaggerAK5 		= process.QGTagger.clone(srcJets = 'ak5PFJets',		jetsLabel = 'QGL_AK5PF')
process.QGTaggerAK7 		= process.QGTagger.clone(srcJets = 'ak7PFJets',		jetsLabel = 'QGL_AK5PF') 			# Change to AK7 on next .db file
process.QGTaggerAK4chs 		= process.QGTagger.clone(srcJets = 'ak4PFJetsCHS',	jetsLabel = 'QGL_AK4PFchs')
process.QGTaggerAK5chs 		= process.QGTagger.clone(srcJets = 'ak5PFJetsCHS',	jetsLabel = 'QGL_AK5PFchs')
process.QGTaggerAK7chs 		= process.QGTagger.clone(srcJets = 'ak7PFJetsCHS',	jetsLabel = 'QGL_AK5PFchs')			# Change to AK7 on next .db file

# Init qgMiniTuple analyzer for AK4, make clones for the other ones
process.qgMiniTupleAK4 = cms.EDAnalyzer("qgMiniTuple",
    fileName 			= cms.untracked.string('qgMiniTuple.root'),
    rhoInputTag			= cms.InputTag('fixedGridRhoFastjetAll'),
    jetsInputTag		= cms.InputTag('ak4PFJets'),
    qgVariablesInputTag		= cms.InputTag('QGTaggerAK4'),
    genParticlesInputTag	= cms.InputTag('genParticles'),
    minJetPt			= cms.untracked.double(20.),
    deltaRcut			= cms.untracked.double(0.3),
    pythia6			= cms.untracked.bool(False)
)
process.qgMiniTupleAK5 		= process.qgMiniTupleAK4.clone(jetsInputTag = 'ak5PFJets',    qgVariablesInputTag = 'QGTaggerAK5')
process.qgMiniTupleAK7 		= process.qgMiniTupleAK4.clone(jetsInputTag = 'ak7PFJets',    qgVariablesInputTag = 'QGTaggerAK7')
process.qgMiniTupleAK4chs 	= process.qgMiniTupleAK4.clone(jetsInputTag = 'ak4PFJetsCHS', qgVariablesInputTag = 'QGTaggerAK4chs')
process.qgMiniTupleAK5chs 	= process.qgMiniTupleAK4.clone(jetsInputTag = 'ak5PFJetsCHS', qgVariablesInputTag = 'QGTaggerAK5chs')
process.qgMiniTupleAK7chs 	= process.qgMiniTupleAK4.clone(jetsInputTag = 'ak7PFJetsCHS', qgVariablesInputTag = 'QGTaggerAK7chs')

process.qgMiniTupleAK4chs.bTagUsed = cms.untracked.bool(True)


# The path: jet sequence + b tagging + (QGTagger + qgMiniTuple) for every jet collection
process.p = cms.Path(process.myRecoPFJets *
#                     process.ak4PFJetsL1FastL2L3 * process.ak4PFchsJetsL1FastL2L3 * 
#                     process.ak5PFJetsL1FastL2L3 * process.ak5PFchsJetsL1FastL2L3 * 
#                     process.ak7PFJetsL1FastL2L3 * process.ak7PFJetsL1FastL2L3 *
                     process.CSVIVFV2ak4chs *
                     process.QGTaggerAK4 * process.qgMiniTupleAK4 * process.QGTaggerAK4chs * process.qgMiniTupleAK4chs *
                     process.QGTaggerAK5 * process.qgMiniTupleAK5 * process.QGTaggerAK5chs * process.qgMiniTupleAK5chs *
                     process.QGTaggerAK7 * process.qgMiniTupleAK7 * process.QGTaggerAK7chs * process.qgMiniTupleAK7chs)


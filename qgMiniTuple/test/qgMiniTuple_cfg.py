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
process.load('RecoJets.Configuration.RecoGenJets_cff')
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.ak7GenJets 	= process.ak4GenJets.clone(rParam = 0.7)
process.ak7PFJetsCHS 	= process.ak4PFJetsCHS.clone(rParam = 0.7) 
process.myRecoPFJets = cms.Sequence(process.fixedGridRhoFastjetAll + process.pfNoPileUpJMESequence + process.genParticlesForJets + 
#                                   process.ak5GenJets + process.ak5PFJets + process.ak5PFJetsCHS +
#                                   process.ak7GenJets + process.ak7PFJets + process.ak7PFJetsCHS +
                                    process.ak4GenJets + process.ak4PFJets + process.ak4PFJetsCHS
)

# Jet energy corrections
process.load('QGDev.qgMiniTuple.correctionSources_cff')

# b-tagging
for jetCollection in ['AK4','AK5','AK7','AK4CHS','AK5CHS','AK7CHS']: process.load('QGDev.qgMiniTuple.RecoBTag' + jetCollection + '_cff')


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
    pfCandidatesInputTag	= cms.InputTag('particleFlow'),
    genJetsInputTag		= cms.InputTag('ak4GenJets'),
    jec				= cms.string('ak4PFL1FastL2L3'),
    genParticlesInputTag	= cms.InputTag('genParticles'),
    csvInputTag			= cms.InputTag('ak4CombinedSecondaryVertexBJetTags'),
#   qgVariablesInputTag		= cms.InputTag('QGTaggerAK4'),
    minJetPt			= cms.untracked.double(20.),
    deltaRcut			= cms.untracked.double(0.3),
)
process.qgMiniTupleAK5 		= process.qgMiniTupleAK4.clone(jetsInputTag = 'ak5PFJets',	genJetsInputTag = 'ak5GenJets',	jec = 'ak5PFL1FastL2L3',	csvInputTag = 'ak5CombinedSecondaryVertexBJetTags')#,    qgVariablesInputTag = 'QGTaggerAK5')
process.qgMiniTupleAK7 		= process.qgMiniTupleAK4.clone(jetsInputTag = 'ak7PFJets',	genJetsInputTag = 'ak7GenJets',	jec = 'ak7PFL1FastL2L3',	csvInputTag = 'ak7CombinedSecondaryVertexBJetTags')#,    qgVariablesInputTag = 'QGTaggerAK7')
process.qgMiniTupleAK4chs 	= process.qgMiniTupleAK4.clone(jetsInputTag = 'ak4PFJetsCHS',	genJetsInputTag = 'ak4GenJets',	jec = 'ak4PFCHSL1FastL2L3',	csvInputTag = 'ak4CHSCombinedSecondaryVertexBJetTags')#, qgVariablesInputTag = 'QGTaggerAK4chs')
process.qgMiniTupleAK5chs 	= process.qgMiniTupleAK4.clone(jetsInputTag = 'ak4PFJetsCHS',	genJetsInputTag = 'ak5GenJets',	jec = 'ak5PFCHSL1FastL2L3',	csvInputTag = 'ak5CHSCombinedSecondaryVertexBJetTags')#, qgVariablesInputTag = 'QGTaggerAK5chs')
process.qgMiniTupleAK7chs 	= process.qgMiniTupleAK4.clone(jetsInputTag = 'ak4PFJetsCHS',	genJetsInputTag = 'ak7GenJets',	jec = 'ak7PFCHSL1FastL2L3',	csvInputTag = 'ak7CHSCombinedSecondaryVertexBJetTags')#, qgVariablesInputTag = 'QGTaggerAK7chs')

# The path: jet sequence + b tagging + (QGTagger + qgMiniTuple) for every jet collection
process.p = cms.Path(process.myRecoPFJets *
                     process.ak4BTagging * process.ak4CHSBTagging *
#                     process.ak5BTagging * process.ak5CHSBTagging *
#                     process.ak7BTagging * process.ak7CHSBTagging *
                     process.qgMiniTupleAK4 * process.qgMiniTupleAK4chs #*
#                     process.qgMiniTupleAK5 * process.qgMiniTupleAK5chs *
#                     process.qgMiniTupleAK7 * process.qgMiniTupleAK7chs
)

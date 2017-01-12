import FWCore.ParameterSet.Config as cms
import glob

process = cms.Process("qgMiniTupleProducer")

# Settings for local tests
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))

isData=False

#######################################
###              GT                 ###
#######################################
# used by photon id and jets
process.load("Configuration.Geometry.GeometryIdeal_cff") 
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")

#mc https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Global_Tags_for_Run2_MC_Producti
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

if (isData):
    process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v3'
else:
    ## tranch IV v6
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'

############## END GT #################

#######################################
###           DB FOR JEC            ###
#######################################
# comment the all block if they are in GT
#### Load from a sqlite db, if not read from the global tag
#process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *

if isData:
    connectString = cms.string('sqlite:jec/Spring16_23Sep2016AllV1_DATA.db')
    tagName = 'Spring16_23Sep2016AllV1_DATA_AK4PFchs'
else:
    connectString = cms.string('sqlite:jec/Spring16_23Sep2016V1_MC.db')
    tagName = 'Spring16_23Sep2016V1_MC_AK4PFchs'

process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_%s'%tagName),
            label  = cms.untracked.string('AK4PFchs')
            ),
      ## here you add as many jet types as you need
      ## note that the tag name is specific for the particular sqlite file 
      ), 
      connect = connectString
     # uncomment above tag lines and this comment to use MC JEC
)
## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')
## 

#######################################
###            Input                ###
#######################################
fileList=['/store/mc/RunIISummer16MiniAODv2/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/349F7699-A5B6-E611-BE96-0CC47A5FBE35.root']


### do not remove the line below!
###FILELIST###

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(fileList)
)

#######################################
###            Output               ###
#######################################

# Use TFileService to put trees from different analyzers in one file
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("qgMiniTuple.root"),
    closeFileFast = cms.untracked.bool(True)
)
#######################################
###             JEC                 ###
#######################################

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

jecLevels= ['L1FastJet',  'L2Relative', 'L3Absolute']
if isData:
	jecLevels.append( 'L2L3Residuals')

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(jecLevels), 'None')  # Do not forget 'L2L3Residual' on data!
)
print "-> Updating the jets collection to run on to 'updatedPatJetsUpdatedJEC' with the new jec in the GT"
process.jecSequence = cms.Sequence( process.patJetCorrFactorsUpdatedJEC* process.updatedPatJetsUpdatedJEC )
############### END JEC ###############

#######################################
###            Config               ###
#######################################
process.qgMiniTupleAK4chs = cms.EDAnalyzer("qgMiniTuple",
    rhoInputTag			= cms.InputTag('fixedGridRhoFastjetAll'),
    csvInputTag			= cms.InputTag('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
    vertexInputTag		= cms.InputTag('offlineSlimmedPrimaryVertices'),
    #jetsInputTag		= cms.InputTag('slimmedJets'),
    jetsInputTag		= cms.InputTag('updatedPatJetsUpdatedJEC'),
    genJetsInputTag		= cms.InputTag('slimmedGenJets'),
    genParticlesInputTag	= cms.InputTag('prunedGenParticles'),
    pfCandidatesInputTag	= cms.InputTag('packedPFCandidates'),
    minJetPt			= cms.untracked.double(20.),
    deltaRcut			= cms.untracked.double(0.3),
    jec				= cms.string(''),						# Ignored when using pat mode
)


#######################################
###              Path               ###
#######################################

#process.p = cms.Path(process.qgMiniTupleAK4chs)
process.p = cms.Path(process.jecSequence * process.qgMiniTupleAK4chs)

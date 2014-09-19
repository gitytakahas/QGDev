import FWCore.ParameterSet.Config as cms

qgDatabaseVersion = 'v1-preliminary'

from CondCore.DBCommon.CondDBSetup_cfi import *
QGPoolDBESSource = cms.ESSource("PoolDBESSource",
      CondDBSetup,
      toGet = cms.VPSet(
        cms.PSet(
            record = cms.string('QGLikelihoodRcd'),
            tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_AK5PFchs'),
            label  = cms.untracked.string('QGL_AK5PFchs')
        ),
        cms.PSet(
            record = cms.string('QGLikelihoodRcd'),
            tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_AK5PF'),
            label  = cms.untracked.string('QGL_AK5PF')
        ),
        cms.PSet(
            record = cms.string('QGLikelihoodRcd'),
            tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_AK4PFchs'),
            label  = cms.untracked.string('QGL_AK4PFchs')
        ),
        cms.PSet(
            record = cms.string('QGLikelihoodRcd'),
            tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_AK4PF'),
            label  = cms.untracked.string('QGL_AK4PF')
        ),
      ),
      connect = cms.string('sqlite:QGL_'+qgDatabaseVersion+'.db')
)

QGTagger = cms.EDProducer('QGTagger',
  srcRho 		= cms.InputTag('fixedGridRhoFastjetAll'),		
  srcVertexCollection	= cms.InputTag('offlinePrimaryVerticesWithBS'),
  srcJets		= cms.InputTag('ak4PFJetsCHS'),
  jetsLabel		= cms.string('QGL_AK4PF'),
  jec			= cms.string(''),
  systematicsLabel	= cms.string('')
)

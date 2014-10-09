# Provide ESSources for JEC missing in CMSSW_7_0_X

import FWCore.ParameterSet.Config as cms
from JetMETCorrections.Configuration.DefaultJEC_cff import *

ak4PFL1Fastjet			= ak5PFL1Fastjet.clone(algorithm = 'AK4PF')
ak4PFL2Relative 		= ak5PFL2Relative.clone(algorithm = 'AK4PF')
ak4PFL3Absolute 		= ak5PFL3Absolute.clone(algorithm = 'AK4PF')

ak4PFCHSL1Fastjet		= ak5PFL1Fastjet.clone(algorithm = 'AK4PFchs')
ak4PFCHSL2Relative 		= ak5PFL2Relative.clone(algorithm = 'AK4PFchs')
ak4PFCHSL3Absolute 		= ak5PFL3Absolute.clone(algorithm = 'AK4PFchs')

ak5PFCHSL1Fastjet		= ak5PFL1Fastjet.clone(algorithm = 'AK5PFchs')
ak5PFCHSL2Relative 		= ak5PFL2Relative.clone(algorithm = 'AK5PFchs')
ak5PFCHSL3Absolute 		= ak5PFL3Absolute.clone(algorithm = 'AK5PFchs')

ak7PFCHSL1Fastjet		= ak5PFL1Fastjet.clone(algorithm = 'AK7PFchs')
ak7PFCHSL2Relative 		= ak5PFL2Relative.clone(algorithm = 'AK7PFchs')
ak7PFCHSL3Absolute 		= ak5PFL3Absolute.clone(algorithm = 'AK7PFchs')

ak4PFL1FastL2L3 		= cms.ESProducer('JetCorrectionESChain', correctors = cms.vstring('ak4PFL1Fastjet','ak4PFL2Relative','ak4PFL3Absolute'))
ak4PFCHSL1FastL2L3 		= cms.ESProducer('JetCorrectionESChain', correctors = cms.vstring('ak4PFCHSL1Fastjet','ak4PFCHSL2Relative','ak4PFCHSL3Absolute'))
ak5PFCHSL1FastL2L3 		= cms.ESProducer('JetCorrectionESChain', correctors = cms.vstring('ak5PFCHSL1Fastjet','ak5PFCHSL2Relative','ak5PFCHSL3Absolute'))
ak7PFCHSL1FastL2L3 		= cms.ESProducer('JetCorrectionESChain', correctors = cms.vstring('ak7PFCHSL1Fastjet','ak7PFCHSL2Relative','ak7PFCHSL3Absolute'))

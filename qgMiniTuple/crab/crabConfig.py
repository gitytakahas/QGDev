from WMCore.Configuration import Configuration
import shutil, os

shutil.copyfile('../../qgMiniTuple/test/qgMiniTupleForMiniAOD_cfg.py', './qgMiniTupleForMiniAOD_cfg.py')

config = Configuration()
config.section_('General')

config.section_('JobType')
config.JobType.psetName    = 'qgMiniTupleForMiniAOD_cfg.py'
config.JobType.pluginName  = 'analysis'
config.JobType.outputFiles = ['qgMiniTuple.root']
config.JobType.allowUndistributedCMSSW = True 

config.section_('Data')
config.Data.inputDataset  = '/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'
config.Data.unitsPerJob   = 5
config.Data.splitting     = 'LumiBased'
config.Data.outLFNDirBase = '/store/user/' + os.environ['USER'] + '/qgMiniTuple_80X/'
config.Data.publication = False

config.section_('Site')
config.Site.storageSite = 'T2_BE_IIHE'

config.section_('User')
config.User.voGroup = 'becms'

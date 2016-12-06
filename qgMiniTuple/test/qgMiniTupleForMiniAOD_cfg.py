import FWCore.ParameterSet.Config as cms
import glob

process = cms.Process("qgMiniTupleProducer")

# Settings for local tests
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


#fileList=[
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/02DA6F10-0B1A-E611-BE2F-7845C4FC3B8A.root", 
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/02DA6F10-0B1A-E611-BE2F-7845C4FC3B8A.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/5285D64A-0A1A-E611-A678-3417EBE65F65.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/9056503A-0A1A-E611-A3FD-7845C4FC3B1B.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/CADE229C-091A-E611-AFCA-3417EBE646E4.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/0AD7EC6D-0D1A-E611-8812-008CFA000B68.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/56152442-091A-E611-8CD5-008CFA0E9E28.root",
#]


#fileList=[
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/0694F42A-031A-E611-87D7-0025905A48D0.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/523EFF33-031A-E611-AA13-24BE05C63721.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/D043E934-031A-E611-A56D-0025905A60CE.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/F8B61C34-031A-E611-A3E5-0025905A60B6.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/1078336C-031A-E611-B8C6-90B11C2CB7A9.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/8AF36C1F-031A-E611-926E-5065F381B211.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/D602842A-031A-E611-A6AF-0CC47A78A3D8.root"
#    ]

#fileList=[
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/00499750-DC19-E611-8340-001517EA7068.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/5A3101BA-DA19-E611-AE81-0CC47A123F98.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/AAA0F4FB-D619-E611-BA8B-0025905AC806.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/028A1AB9-DA19-E611-9715-001EC9ADDBFF.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/5AB025C1-DA19-E611-ACAE-001EC9ADE794.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/AAF3EF89-D719-E611-B701-A0369F310374.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/0428C060-DE19-E611-9E46-0025900EB230.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/5C0DBCD4-DB19-E611-924F-0025905AC806.root", 
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/AC4A6A75-D719-E611-BC71-24BE05CE2D41.root",    
#    ]

#fileList=[
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/0223FB02-DF1A-E611-985E-0CC47A1DF64A.root", 
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/5C714A01-DF1A-E611-9421-0CC47A1DF7E6.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/B8AC48C6-471B-E611-901E-008CFA152104.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/088A2C2E-411C-E611-BD74-0CC47A1E0472.root",  
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/5CEF82AE-DF1A-E611-B378-A0369F7FC770.root",  
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/B8BC0981-421C-E611-9604-008CFA05EA18.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/08D4FE15-DF1A-E611-BBEC-549F35AC7DFB.root",  
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/5EA6CFB7-DF1A-E611-BFD1-A0369F7FC8E8.root", 
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/B8D754B6-421C-E611-BE50-00266CFBE88C.root"
#]


fileList=[
    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/3031A43C-5A1C-E611-A856-14187741136B.root",
    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/6629E139-5B1C-E611-8070-7845C4FBBD07.root",
    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/9C7D540C-591C-E611-A017-002590FD5A3A.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/CC4AE144-5A1C-E611-9D99-B083FED43141.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/04036652-5A1C-E611-958B-0019B9CB030F.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/32F093FA-591C-E611-9F3A-6CC2173D9860.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/68BB5179-621C-E611-9DC6-782BCB20F910.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/A0137EF3-581C-E611-9890-1418774117C7.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/CCF816BC-591C-E611-B1B1-782BCB21110D.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/065EC87A-591C-E611-A1DB-000E1E878860.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/3453A438-5B1C-E611-85A4-C81F66C0F057.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/6CD50D62-591C-E611-9996-0090FAA57EA4.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/A069F40D-591C-E611-8FC2-141877410512.root",
#    "/store/mc/RunIISpring16MiniAODv2/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/CEC77B4A-591C-E611-959C-002590D9D9DA.root"
]

### do not remove the line below!
###FILELIST###

process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(fileList), # Test file available at T2B
#fileNames = cms.untracked.vstring("/store/mc/RunIISpring16MiniAODv2/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/3031A43C-5A1C-E611-A856-14187741136B.root"),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

# Standard configurations
process.load('Configuration.StandardSequences.Services_cff')

# Use TFileService to put trees from different analyzers in one file
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("qgMiniTuple_1.root"),
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

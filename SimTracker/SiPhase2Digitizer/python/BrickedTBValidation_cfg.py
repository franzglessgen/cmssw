
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2_cff import Phase2
process = cms.Process('digiTest',Phase2)
process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(10000)
 )
 #process.MessageLogger = cms.Service("MessageLogger",
#    debugModules = cms.untracked.vstring('siPixelRawData'),
#    destinations = cms.untracked.vstring("cout"),
 #    cout = cms.untracked.PSet(
 #        threshold = cms.untracked.string('ERROR')
 #    )
 #)
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D57Reco_cff')

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
 
# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
# list of files

process.source = cms.Source("PoolSource",
     fileNames =  cms.untracked.vstring(
         'file:Digis_10000.root'
       )
 )
 # Production Info
process.configurationMetadata = cms.untracked.PSet(
     version = cms.untracked.string('$Revision: 1.19 $'),
     annotation = cms.untracked.string('step1 nevts:1'), 
     name = cms.untracked.string('Applications')
 )
# Output definition


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo_TB_10000.root')
)
 
process.DQMoutput = cms.OutputModule("PoolOutputModule",
     splitLevel = cms.untracked.int32(0),
     outputCommands = process.DQMEventContent.outputCommands,
     fileName = cms.untracked.string('step1_TBVal.root'),
     dataset = cms.untracked.PSet(
     filterName = cms.untracked.string(''),
     dataTier = cms.untracked.string('')))

#process.load('SimTracker.SiPhase2Digitizer.Phase2TrackerMonitorDigi_cff')
process.load('SimTracker.SiPhase2Digitizer.BrickedTBValidation_cff')

#process.digiana_seq = cms.Sequence(process.pixDigiMon * process.otDigiMon * process.pixDigiValid * process.otDigiValid)
 
process.digiana_seq = cms.Sequence(process.pixelcells)

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)
#process.digi_step = cms.Sequence(process.siPixelRawData*process.siPixelDigis)
#process.p = cms.Path(process.digiana_seq * process.dqm_comm )
process.p = cms.Path(process.digiana_seq)

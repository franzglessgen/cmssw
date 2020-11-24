# Imports
import FWCore.ParameterSet.Config as cms

# Create a new CMS process
process = cms.Process('HitTest')

# Import all the necessary files
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D57Reco_cff') 
#commented for bricked
#process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# clusterizer 
process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")

# for raw
#process.load("EventFilter.SiPixelRawToDigi.SiPixelDigiToRaw_cfi")
#process.load("EventFilter.SiPixelRawToDigi.SiPixelRawToDigi_cfi")
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')


# needed for pixel RecHits (templates?)
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')


# Number of events (-1 = all)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input file
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring('file:RecoHits_1000.root')
)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('file:rechits_validation.root')
)

#process.load('RecoLocalTracker.SiPhase2Clusterizer.phase2TrackerClusterizer_cfi')
#process.load('RecoLocalTracker.Phase2TrackerRecHits.Phase2StripCPEESProducer_cfi')
#process.load('RecoLocalTracker.Phase2TrackerRecHits.Phase2StripCPEGeometricESProducer_cfi')

#process.load('RecoLocalTracker.Phase2TrackerRecHits.Phase2TrackerRecHits_cfi')
#process.load('RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi')


#process.siPhase2RecHits.Phase2StripCPE = cms.ESInputTag("phase2StripCPEESProducer", "Phase2StripCPE")
#process.siPhase2RecHits.Phase2StripCPE = cms.ESInputTag("phase2StripCPEGeometricESProducer", "Phase2StripCPEGeometric")


# Analyzer
process.analysis = cms.EDAnalyzer('BrickedRecHits',
    #src = cms.InputTag("siPhase2RecHits"),
    #clusters = cms.InputTag("siPhase2Clusters"),
    #links = cms.InputTag("simSiPixelDigis", "Tracker"),
   
    src = cms.InputTag("siPixelRecHits"),
    clusters = cms.InputTag("siPixelClusters"),
    links = cms.InputTag("simSiPixelDigis", "Pixel"),
    
    PSimHitSourceB  = cms.VInputTag('g4SimHits:TrackerHitsPixelBarrelLowTof',
                                   'g4SimHits:TrackerHitsPixelBarrelHighTof',
                                   'g4SimHits:TrackerHitsPixelEndcapLowTof',
                                   'g4SimHits:TrackerHitsPixelEndcapHighTof'),

    PSimHitSourceE  = cms.VInputTag('g4SimHits:TrackerHitsPixelBarrelLowTof',
                                   'g4SimHits:TrackerHitsPixelBarrelHighTof',
                                   'g4SimHits:TrackerHitsPixelEndcapLowTof',
                                   'g4SimHits:TrackerHitsPixelEndcapHighTof'),
    
    #simhitsbarrel = cms.InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof"),
    #simhitsendcap = cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof"),
    
    

    simtracks = cms.InputTag("g4SimHits"),
    ECasRings = cms.bool(True),
    SimTrackMinPt = cms.double(0.0),
    MakeEtaPlots = cms.bool(True),
    MinEta = cms.double(0.),
    MaxEta = cms.double(10.)
)

# Processes to run
#process.rechits_step = cms.Path(process.siPhase2Clusters + process.siPhase2RecHits)
#process.rechits_step = cms.Path(process.siPhase2RecHits)
process.rechits_step = cms.Path(process.siPixelRecHits)
process.validation_step = cms.Path(process.analysis)

process.schedule = cms.Schedule(process.rechits_step, process.validation_step)


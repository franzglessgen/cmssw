# Imports
import FWCore.ParameterSet.Config as cms

# Create a new CMS process
process = cms.Process('cluTest')

# Import all the necessary files
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.EventContent.EventContent_cff')

#bricked geometry
process.load('Configuration.Geometry.GeometryExtended2026D57Reco_cff') 
process.load('Configuration.StandardSequences.MagneticField_cff')

#process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')


# Number of events (-1 = all)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

# Input file
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring('file:RecoHits_NB_10000.root')
)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('file:cluster_NB_validation.root')
)

# DEBUG
#process.MessageLogger = cms.Service('MessageLogger',
#	debugModules = cms.untracked.vstring('siPhase2Clusters'),
#	destinations = cms.untracked.vstring('cout'),
#	cout = cms.untracked.PSet(
#		threshold = cms.untracked.string('ERROR')
#	)
#)

# Analyzer
process.analysis = cms.EDAnalyzer('NBClust',
    #src = cms.InputTag("siPhase2Clusters"),
    #src = cms.InputTag("simSiPixelClusters", "Pixel"),
    src = cms.InputTag("siPixelClusters"),
    #links = cms.InputTag("simSiPixelDigis", "Tracker"),
    links = cms.InputTag("simSiPixelDigis", "Pixel"),
    
    PSimHitSourceB  = cms.VInputTag('g4SimHits:TrackerHitsPixelBarrelLowTof',
                                   'g4SimHits:TrackerHitsPixelBarrelHighTof',
                                   'g4SimHits:TrackerHitsPixelEndcapLowTof',
                                   'g4SimHits:TrackerHitsPixelEndcapHighTof'),
                                   #'g4SimHits:TrackerHitsTIBLowTof',
                                   #'g4SimHits:TrackerHitsTIBHighTof',
                                   #'g4SimHits:TrackerHitsTIDLowTof',
                                   #'g4SimHits:TrackerHitsTIDHighTof',
                                   #'g4SimHits:TrackerHitsTOBLowTof',
                                   #'g4SimHits:TrackerHitsTOBHighTof',
                                   #'g4SimHits:TrackerHitsTECLowTof',
                                   #'g4SimHits:TrackerHitsTECHighTof'),

    PSimHitSourceE  = cms.VInputTag('g4SimHits:TrackerHitsPixelBarrelLowTof',
                                   'g4SimHits:TrackerHitsPixelBarrelHighTof',
                                   'g4SimHits:TrackerHitsPixelEndcapLowTof',
                                   'g4SimHits:TrackerHitsPixelEndcapHighTof'),
                                   #'g4SimHits:TrackerHitsTIBLowTof',
                                   #'g4SimHits:TrackerHitsTIBHighTof',
                                   #'g4SimHits:TrackerHitsTIDLowTof',
                                   #'g4SimHits:TrackerHitsTIDHighTof',
                                   #'g4SimHits:TrackerHitsTOBLowTof',
                                   #'g4SimHits:TrackerHitsTOBHighTof',
                                   #'g4SimHits:TrackerHitsTECLowTof',
                                   #'g4SimHits:TrackerHitsTECHighTof'),


    #simhitsbarrel = cms.InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof"),
    #simhitsendcap = cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof"),
    simtracks = cms.InputTag("g4SimHits"),
    ECasRings = cms.bool(True),
    SimTrackMinPt = cms.double(0.0) #2.
)

# Processes to run
process.p = cms.Path(process.analysis)

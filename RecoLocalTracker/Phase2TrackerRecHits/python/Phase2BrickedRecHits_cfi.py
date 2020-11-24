
import FWCore.ParameterSet.Config as cms

# RecHits options
siPhase2BrickedRecHits = cms.EDAnalyzer("BrickedRecHits",
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

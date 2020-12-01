import FWCore.ParameterSet.Config as cms

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
rechitValidIT = DQMEDAnalyzer('Phase2ITValidateRecHit',#Phase2TrackerValidateRecHit,
	    Verbosity = cms.bool(False),
	    TopFolderName = cms.string("Phase2TrackerRecHitV"),
	    ITPlotFillingFlag = cms.bool(False),
	    rechitsSrc = cms.InputTag("siPixelRecHits"),#rechitsIT
	    InnerPixelDigiSource   = cms.InputTag("simSiPixelDigis","Pixel"),                          
	    InnerPixelDigiSimLinkSource = cms.InputTag("simSiPixelDigis", "Pixel"), 
	    PSimHitSource  = cms.VInputTag('g4SimHits:TrackerHitsPixelBarrelLowTof',
					   'g4SimHits:TrackerHitsPixelBarrelHighTof',
					   'g4SimHits:TrackerHitsPixelEndcapLowTof',
					   'g4SimHits:TrackerHitsPixelEndcapHighTof',
					   'g4SimHits:TrackerHitsTIBLowTof',
					   'g4SimHits:TrackerHitsTIBHighTof',
					   'g4SimHits:TrackerHitsTIDLowTof',
					   'g4SimHits:TrackerHitsTIDHighTof',
					   'g4SimHits:TrackerHitsTOBLowTof',
					   'g4SimHits:TrackerHitsTOBHighTof',
					   'g4SimHits:TrackerHitsTECLowTof',
					   'g4SimHits:TrackerHitsTECHighTof'),
	      simTracksSrc = cms.InputTag("g4SimHits"), #SimTrackSource 
	      SimVertexSource = cms.InputTag("g4SimHits"),
	      usePhase2Tracker = cms.bool(True),#these are used by simHit assoc.
	      associatePixel = cms.bool(True),
	      SimTrackMinPt = cms.double(0.0),
              associateRecoTracks = cms.bool(False),
	      associateStrip = cms.bool(False),
	      associateHitbySimTrack = cms.bool(True),
	      pixelSimLinkSrc = cms.InputTag("simSiPixelDigis","Pixel"),
	      ROUList  =  cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof',
				      'g4SimHitsTrackerHitsPixelBarrelHighTof',
				      'g4SimHitsTrackerHitsPixelEndcapLowTof',
				      'g4SimHitsTrackerHitsPixelEndcapHighTof',
				  ),

) 


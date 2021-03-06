import FWCore.ParameterSet.Config as cms

#
# Geometry master configuration
#
# Ideal geometry, needed for simulation
DDDetectorESProducer = cms.ESSource("DDDetectorESProducer",
                                    confGeomXMLFiles = cms.FileInPath('Geometry/CMSCommonData/data/dd4hep/cmsExtendedGeometry2026D41.xml'),
                                    appendToDataLabel = cms.string('')
)

DDSpecParRegistryESProducer = cms.ESProducer("DDSpecParRegistryESProducer",
                                             appendToDataLabel = cms.string('')
)

DDVectorRegistryESProducer = cms.ESProducer("DDVectorRegistryESProducer",
                                            appendToDataLabel = cms.string(''))

DDCompactViewESProducer = cms.ESProducer("DDCompactViewESProducer",
                                         appendToDataLabel = cms.string('')
)

import FWCore.ParameterSet.Config as cms

skimAllPathsFilter = cms.EDFilter(
    "SkimAllPathsFilter",
    #electronCollection = cms.InputTag("gedGsfElectrons"),
    electronCollection = cms.InputTag("gsfElectrons"),
    muonCollection = cms.InputTag("muons"),
    mode = cms.string("ML"),
    eleID = cms.string("VBTF95")
)


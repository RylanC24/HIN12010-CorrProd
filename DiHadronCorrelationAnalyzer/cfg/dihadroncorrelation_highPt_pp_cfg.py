import FWCore.ParameterSet.Config as cms

process = cms.Process("corr")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_P_V43D::All'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)
process.HeavyIonGlobalParameters = cms.PSet(
        centralityVariable = cms.string("HFtowersTrunc"), #or HFtowersPlusTrunc
        nonDefaultGlauberModel = cms.string(""),
        centralitySrc = cms.InputTag("pACentrality"),
        pPbRunFlip = cms.untracked.uint32(211313)
        )

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'/store/user/davidlw/PPJet/PP2013_FlowCorr_PromptReco_HighPtTrk_NoGplus_v2/2a01db98f9c5d38ba9c3802e841ac721/pp_HighPtTrk_100_1_oc4.root'
                )
                            )
process.load("FlowCorrAna.DiHadronCorrelationAnalyzer.dihadroncorrelation_cff")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

# Additional output definition
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('dihadroncorrelation_pp.root')
                                   )

process.ana_step = cms.Path(process.corr_ana_pp_highPt)

import FWCore.ParameterSet.Config as cms

from FlowCorrAna.DiHadronCorrelationAnalyzer.dihadroncorrelation_cfi import *

corr_ana_HI = corr_ana.clone(
#  TrgTrackCollection = cms.string('hiLowPtPixelTracks'),
  TrgTrackCollection = cms.string('hiGeneralTracks'),
  VertexCollection = cms.string('hiSelectedVertex'),
  GenParticleCollection = cms.string('hiGenParticles'),

  rhomin = cms.double(0.0),
  rhomax = cms.double(0.02),
  xvtxcenter = cms.double(0.077),
  yvtxcenter = cms.double(0.037),
  zvtxcenter = cms.double(-0.54),

  IsHI = cms.bool(True),
  IsHITrkQuality = cms.bool(True),
  IsPPTrkQuality = cms.bool(False),
)

corr_ana_HI_highPt = corr_ana_HI.clone(
  IsDebug = cms.bool(False),
  IsHarmonics = cms.bool(False),
  pttrgmin = cms.vdouble(12.0, 14.0, 16.0, 19.2, 24.0, 28.8, 35.2),
  pttrgmax = cms.vdouble(14.0, 16.0, 19.2, 24.0, 28.8, 35.2, 48.0),
  ptassmin = cms.vdouble(0.5,1.0,2.0,3.0,4.0,6.0,8.0),
  ptassmax = cms.vdouble(1.0,2.0,3.0,4.0,6.0,8.0,12.0)
)

corr_ana_pp_highPt = corr_ana.clone(
  IsDebug = cms.bool(False),
  IsHarmonics = cms.bool(False),
  pttrgmin = cms.vdouble(12.0, 14.0, 16.0, 19.2, 24.0, 28.8, 35.2),
  pttrgmax = cms.vdouble(14.0, 16.0, 19.2, 24.0, 28.8, 35.2, 48.0),
  ptassmin = cms.vdouble(0.5,1.0,2.0,3.0,4.0,6.0,8.0),
  ptassmax = cms.vdouble(1.0,2.0,3.0,4.0,6.0,8.0,12.0),
  EffFileName = cms.string('../data/TrkCorr/TrackCorrection_PP_pTmax50_Oct28.root')
)

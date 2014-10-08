import FWCore.ParameterSet.Config as cms

#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
hltHM220 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
hltHM220.HLTPaths = ['HLT_PAPixelTracks_Multiplicity100_v*','HLT_PAPixelTracks_Multiplicity130_v*','HLT_PAPixelTracks_Multiplicity160_v*','HLT_PAPixelTracks_Multiplicity190_v*']
hltHM220.andOr = cms.bool(True)
hltHM220.throw = cms.bool(False)

hltHM185220 = hltHM220.clone()
hltHM185220.HLTPaths = ['HLT_PAPixelTracks_Multiplicity100_v*','HLT_PAPixelTracks_Multiplicity130_v*','HLT_PAPixelTracks_Multiplicity160_v*']

hltHM150185 = hltHM220.clone()
hltHM150185.HLTPaths = ['HLT_PAPixelTracks_Multiplicity100_v*','HLT_PAPixelTracks_Multiplicity130_v*']

hltHM120150 = hltHM220.clone()
hltHM120150.HLTPaths = ['HLT_PAPixelTracks_Multiplicity100_v*']

hltMB = hltHM220.clone()
hltMB.HLTPaths = ['HLT_PAZeroBiasPixel_SingleTrack_v*']

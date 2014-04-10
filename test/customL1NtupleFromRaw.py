import FWCore.ParameterSet.Config as cms

from L1TriggerDPG.L1Menu.customL1Ntuple_cfg import *

process.p.remove(process.l1RecoTreeProducer)
process.p.remove(process.l1MuonRecoTreeProducer)

# edit here
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

readFiles.extend( ['file:///data2/battilan/L1Trigger/62X_RAW_RECO.root'] )

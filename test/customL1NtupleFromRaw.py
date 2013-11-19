import FWCore.ParameterSet.Config as cms

customDTTF  = True
customCSCTF = False
customPACT  = True
useUCT2015  = False

# make ntuples from RAW (ie. remove RECO)
from L1TriggerDPG.L1Ntuples.l1Ntuple_cfg import *

process.p.remove(process.l1RecoTreeProducer)
process.p.remove(process.l1MuonRecoTreeProducer)
process.p.remove(process.l1MenuTreeProducer)

if customDTTF or customCSCTF or customPACT :
    from L1TriggerDPG.L1TMenu.customiseL1Muons_cff import *
    customiseL1Muons(process,customDTTF,customCSCTF,customPACT)

# edit here
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.GlobalTag.globaltag = 'GR_P_V41_AN1::All'

readFiles.extend(["file:///data2/battilan/L1Trigger/202299_RAW-RECO.root"])

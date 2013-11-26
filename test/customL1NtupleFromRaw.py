import FWCore.ParameterSet.Config as cms

reEmulation = True
reEmulMuons = True
reEmulCalos = True
patchNtuple = True

runOnMC = True

customDTTF  = True
customCSCTF = False
customPACT  = False

useUCT2015  = True

# make ntuples from RAW (ie. remove RECO)
from L1TriggerDPG.L1Ntuples.l1Ntuple_cfg import *

process.p.remove(process.l1RecoTreeProducer)
process.p.remove(process.l1MuonRecoTreeProducer)
process.p.remove(process.l1MenuTreeProducer)

if reEmulation :
    from L1TriggerDPG.L1TMenu.reEmulation_cff import *
    reEmulation(process, reEmulMuons, reEmulCalos, patchNtuple)
    process.p.replace(process.l1NtupleProducer, process.reEmul + process.l1NtupleProducer)

if reEmulation and (customDTTF or customCSCTF or customPACT) :
    from L1TriggerDPG.L1TMenu.customiseL1Muons_cff import *
    customiseL1Muons(process, customDTTF, customCSCTF, customPACT)

if reEmulation and useUCT2015 :
    from L1TriggerDPG.L1TMenu.customiseL1Calos_cff import *
    customiseUCT2015(process, runOnMC)

# edit here
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.GlobalTag.globaltag = 'GR_P_V41_AN1::All'

readFiles.extend(["file:///data2/battilan/L1Trigger/202299_RAW-RECO.root"])

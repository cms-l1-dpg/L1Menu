import FWCore.ParameterSet.Config as cms

from L1TriggerDPG.L1Ntuples.l1Ntuple_cfg import *

# configuration parameters
readFiles.extend(["file:///data2/battilan/L1Trigger/202299_RAW-RECO.root"])
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.GlobalTag.globaltag = 'GR_P_V41_AN1::All'
process.GlobalTag.globaltag = 'START53_V19D::All'
process.GlobalTag.toGet = cms.VPSet()

reEmulation     = False
reEmulMuons     = True
reEmulCalos     = True
patchNtuple     = True
force2012Config = True

runOnMC       = True
keepEDMOutput = True

customDTTF  = False
customCSCTF = False
customPACT  = False
useUCT2015  = False

# make ntuples from RAW (ie. remove RECO)

process.p.remove(process.l1RecoTreeProducer)
process.p.remove(process.l1MuonRecoTreeProducer)
process.p.remove(process.l1MenuTreeProducer)

if reEmulation :
    from L1TriggerDPG.L1TMenu.reEmulation_cff import *
    reEmulation(process, reEmulMuons, reEmulCalos, patchNtuple)
    process.p.replace(process.l1NtupleProducer, process.reEmul + process.l1NtupleProducer)
    if force2012Config :
         run2012CConfiguration(process)

if reEmulation and (customDTTF or customCSCTF or customPACT) :
    from L1TriggerDPG.L1TMenu.customiseL1Muons_cff import *
    customiseL1Muons(process, customDTTF, customCSCTF, customPACT)

if reEmulation and useUCT2015 :
    from L1TriggerDPG.L1TMenu.customiseL1Calos_cff import *
    customiseUCT2015(process, runOnMC)

if keepEDMOutput :
    
    process.output = cms.OutputModule("PoolOutputModule",
                                      fileName = cms.untracked.string('L1GmtGt.root'),
                                      outputCommands = cms.untracked.vstring('drop *',
                                                                             'keep *_gtDigis_*_*',
                                                                             'keep *_gtReEmulDigis_*_*',
                                                                             'keep *_gmtReEmulDigis_*_*',
                                                                             'keep *_rpcTriggerReEmulDigis_*_*',
                                                                             'keep *_csctfReEmulDigis_*_*',
                                                                             'keep *_dttfReEmulDigis_*_*')
                                  )

    process.out = cms.EndPath(process.output)

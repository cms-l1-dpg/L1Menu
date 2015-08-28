import FWCore.ParameterSet.Config as cms

# General config options
import FWCore.ParameterSet.VarParsing as VarParsing
import sys

options = VarParsing.VarParsing()

options.register('globalTag',
                 'GR_P_V41_AN1', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Global Tag")

options.register('reEmulation',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Run re-emulation")

options.register('reEmulMuons',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Run re-emulation of L1 muons")

options.register('reEmulCalos',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Run re-emulation of L1 calos")

options.register('reEmulRCT',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Run re-emulation of L1 RCT")

options.register('patchNtuple',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Patch ntuple inputs to use re-emulation ones")

options.register('jetSeedThr10GeV',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Switches on 10 GeV jet Seed Thresholds for 2012 GCT")

options.register('runOnMC',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Set to True when running on MC")

options.register('runOnPostLS1',
                 True, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Set to True when running on MC and this postLS1")

options.register('whichPU',
                 40, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of average PU interactions for UCT PU subtractions")

options.register('keepEDMOutput',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "When True keeps also EDM GMT/GT skimmmed collections")

options.register('customDTTF',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Enables usage of new DTTF LUTs")

options.register('dttfLutsFile',
                 'sqlite:../data/dttf_config.db', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "DTTF LUTs sqlite input file")

options.register('customCSCTF',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Enables usage of new CSCTF FW and LUTs (actually does nothing right now)")

options.register('customPACT',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Enables usage of new RPC PACT patterns")

options.register('customGMT',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Switches to minPT for the GMT")

options.register('useStage1Layer2',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Enables new Stage1Layer2 emulation for calos")

options.register('puReweightingFile',
                 'none', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "InputFile to be used for PU reweighting (for example to scale luminosity)")

options.parseArguments()


#L1 ntuple
from L1TriggerDPG.L1Ntuples.l1Ntuple_cfg import *

if options.runOnMC and hasattr(process,'l1NtupleProducer') :
    print "[L1Menu]: Running on MC reading also PileUpSummary info"
    ntuple                  = getattr(process,'l1NtupleProducer')
    ntuple.simulationSource = cms.InputTag("addPileupInfo")


print "[L1Menu]: Using GlobalTag", options.globalTag
process.GlobalTag.globaltag = options.globalTag
process.GlobalTag.toGet     = cms.VPSet()

# make ntuples from RAW (ie. remove RECO)
if options.reEmulMuons :
    process.p.remove(process.muonDTDigis)

# re-emulation customisations

if options.reEmulation :
    from L1TriggerDPG.L1Menu.reEmulation_cff import *
    reEmulation(process, options.reEmulMuons, options.reEmulCalos, options.patchNtuple, options.runOnPostLS1, options.useStage1Layer2, options.reEmulRCT)
    process.p.replace(process.l1NtupleProducer, process.reEmul + process.l1NtupleProducer)

if options.reEmulCalos and options.jetSeedThr10GeV :
    from L1TriggerDPG.L1Menu.customiseL1Calos_cff import *
    set10GCTtreshold(process)

if options.reEmulation and (options.customDTTF or options.customCSCTF or options.customPACT or options.customGMT ) :
    from L1TriggerDPG.L1Menu.customiseL1Muons_cff import *
    customiseL1Muons(process, options.customDTTF, options.customCSCTF, options.customPACT, options.customGMT, options.dttfLutsFile)

if options.reEmulation and options.useStage1Layer2 :
    from L1TriggerDPG.L1Menu.customiseL1Calos_cff import *
    if options.useStage1Layer2:
        customiseStage1(process, options.runOnMC, options.runOnPostLS1, options.whichPU)

if options.puReweightingFile != "none" :
    from L1TriggerDPG.L1Menu.pileUpReweighting_cff import *
    pileUpReweighting(process,options.puReweightingFile, "productionPileUpHisto", "targetPileUpHisto")


# EDM keep statement

if options.keepEDMOutput :
    
    process.output = cms.OutputModule("PoolOutputModule",
                                      fileName = cms.untracked.string('L1GmtGt.root'),
                                      outputCommands = cms.untracked.vstring('drop *',
                                                                             'keep *_gtDigis_*_*',
                                                                             'keep *_gtReEmulDigis_*_*',
                                                                             'keep *_gmtReEmulDigis_*_*',
                                                                             'keep *_rpcTriggerReEmulDigis_*_*',
                                                                             'keep *_csctfReEmulDigis_*_*',
                                                                             'keep *_dttfReEmulDigis_*_*',
                                                                             'keep *_caloStage1LegacyFormatDigis_*_*',
                                                                             'keep *BXVector_*__L1TEMULATION',
                                                                             'keep *_gctDigis_*_*')
                                  )

    process.out = cms.EndPath(process.output)

import FWCore.ParameterSet.Config as cms

from L1TriggerDPG.L1Menu.customL1Ntuple_cfg import *
from Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff import *

process.p.remove(process.l1RecoTreeProducer)
process.p.remove(process.l1MuonRecoTreeProducer)
process.p.remove(process.l1MenuTreeProducer)
process.p.remove(process.csctfDigis)
process.p.remove(process.gtEvmDigis)

# uncomment the following lines to override the L1RCT configuration parameters in the GlobalTag
if options.reEmulRCT :
    recordOverrides = { ('L1RCTParametersRcd', None) :
                          ##('L1RCTParametersRcd_L1TDevelCollisions_ExtendedScaleFactors_NewTau_FullEGTransparency_v1', None) }   #This is like 50 ns but with transparency corrections at RCT
                          ##('L1RCTParametersRcd_L1TDevelCollisions_ExtendedScaleFactors_EGOnly_v1', None) }   #This is was we use in 50 ns, new RCT for EG, old RCT for hadronic
                          ('L1RCTParametersRcd_L1TDevelCollisions_ExtendedScaleFactorsV4', None) }    #This is what we use in 25 ns
    process.GlobalTag = GlobalTag(process.GlobalTag,'MCRUN2_74_V9A', recordOverrides) #for MC emulation
    process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_HLT_v1', recordOverrides)  # for data emulation

#This is to use the new JEC
if options.reEmulCalos and options.useStage1Layer2:
    process.caloStage1Params.jetCalibrationLUTFile = cms.FileInPath("L1TriggerDPG/L1Menu/data/Jet_Stage1_2015_v2.txt")

# Get the ECAL transparency corrections
# DO NOT USE UNTIL ECAL UNDERSTANDS WHAT'S GOING ON
#process.GlobalTag.toGet = cms.VPSet(
#    cms.PSet(record = cms.string("EcalTPGLinearizationConstRcd"),
#             tag = cms.string("EcalTPGLinearizationConst_weekly_test2_hlt"),
#             connect =cms.untracked.string('frontier://FrontierPrep/CMS_CONDITIONS')
#    )
#)
#process.GlobalTag.toGet = cms.VPSet(
#    cms.PSet(record = cms.string("EcalTPGTowerStatusRcd"),
#             tag = cms.string("EcalTPGTowerStatus_confid722_plus_ebm10_tt64"),
#             connect =cms.untracked.string('frontier://FrontierPrep/CMS_CONDITIONS')
#             )
#)

##To reemulate on data you need to remove Randoms for the ZeroBias dataset
if not options.runOnMC :
    import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt
    process.ZBfilter = hlt.triggerResultsFilter.clone()
    process.ZBfilter.triggerConditions = cms.vstring('HLT_ZeroBias_*',)
    process.ZBfilter.hltResults = cms.InputTag( "TriggerResults", "", "HLT" )
    process.ZBfilter.l1tResults = cms.InputTag("")
    process.ZBfilter.throw = cms.bool( False )
    process.FilterPath = cms.Path(process.ZBfilter)

process.TFileService.fileName = "L1Tree.root"
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

## readFiles.extend( ['file:///data2/battilan/L1Trigger/62X_RAW_RECO.root'] )
## readFiles.extend( ['file:/uscmst1b_scratch/lpc1/lpctrig/apana/L1Upgrade/262AA156-744A-E311-9829-002618943945.root'] )
## readFiles.extend( ['root://eoscms.cern.ch//eos/cms/store/group/comm_trigger/L1Trigger/apana/DYJetsToLL_M-50_13TeV-pythia6_Fall13dr-tsg_PU40bx25__skim_150_1_6UN.root'] )
## readFiles.extend( ['root://xrootd.unl.edu//store/data/Commissioning2014/Cosmics/RAW/v3/000/228/929/00000/509C6C80-3164-E411-92E0-02163E00FFE1.root'] )
## readFiles.extend( ['root://lxcms02//data2/p/pellicci/L1DPG/root/RelValTTbar_730_GENSIMRECO.root'] )
## readFiles.extend( ['/store/data/Run2015C/HLTPhysicspart0/RAW/v1/000/254/790/00000/FE944ECF-5148-E511-96F3-02163E012BBB.root'] )

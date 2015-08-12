import FWCore.ParameterSet.Config as cms

from L1TriggerDPG.L1Menu.customL1Ntuple_cfg import *

process.p.remove(process.l1RecoTreeProducer)
process.p.remove(process.l1MuonRecoTreeProducer)
process.p.remove(process.l1MenuTreeProducer)
process.p.remove(process.csctfDigis)

# uncomment the following lines to override the L1RCT configuration parameters in the GlobalTag
recordOverrides = { ('L1RCTParametersRcd', None) :
                    ('L1RCTParametersRcd_L1TDevelCollisions_ExtendedScaleFactors_EGOnly_v1', None) }
                    ## ('L1RCTParametersRcd_L1TDevelCollisions_ExtendedScaleFactorsV4', None) }
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V9', recordOverrides)

### Get the ECAL transparency corrections
#process.GlobalTag.toGet = cms.VPSet(
#    cms.PSet(record = cms.string("EcalTPGLinearizationConstRcd"),
#             tag = cms.string("EcalTPGLinearizationConst_weekly_test2_hlt"),
#             connect =cms.untracked.string('frontier://FrontierPrep/CMS_CONDITIONS')
#    )
#)

OUTFILE="L1Tree.root"
NEVTS=200

process.TFileService.fileName=OUTFILE
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(NEVTS) )

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

## readFiles.extend( ['file:///data2/battilan/L1Trigger/62X_RAW_RECO.root'] )
## readFiles.extend( ['file:/uscmst1b_scratch/lpc1/lpctrig/apana/L1Upgrade/262AA156-744A-E311-9829-002618943945.root'] )
## readFiles.extend( ['root://eoscms.cern.ch//eos/cms/store/group/comm_trigger/L1Trigger/apana/DYJetsToLL_M-50_13TeV-pythia6_Fall13dr-tsg_PU40bx25__skim_150_1_6UN.root'] )
## readFiles.extend( ['root://xrootd.unl.edu//store/data/Commissioning2014/Cosmics/RAW/v3/000/228/929/00000/509C6C80-3164-E411-92E0-02163E00FFE1.root'] )
## readFiles.extend( ['root://lxcms02//data2/p/pellicci/L1DPG/root/RelValTTbar_730_GENSIMRECO.root'] )


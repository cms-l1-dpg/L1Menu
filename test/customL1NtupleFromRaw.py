import FWCore.ParameterSet.Config as cms

from L1TriggerDPG.L1Menu.customL1Ntuple_cfg import *

process.p.remove(process.l1RecoTreeProducer)
process.p.remove(process.l1MuonRecoTreeProducer)

# edit here

OUTFILE="L1Tree.root"
NEVTS=-1

process.TFileService.fileName=OUTFILE
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(NEVTS) )

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

## readFiles.extend( ['file:///data2/battilan/L1Trigger/62X_RAW_RECO.root'] )
## readFiles.extend( ['file:/uscmst1b_scratch/lpc1/lpctrig/apana/L1Upgrade/262AA156-744A-E311-9829-002618943945.root'] )
readFiles.extend( ['root://eoscms.cern.ch//eos/cms/store/group/comm_trigger/L1Trigger/apana/DYJetsToLL_M-50_13TeV-pythia6_Fall13dr-tsg_PU40bx25__skim_150_1_6UN.root'] )

## processDumpFile = open('config.dump', 'w')
## print >> processDumpFile, process.dumpPython()

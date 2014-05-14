import FWCore.ParameterSet.Config as cms

from L1TriggerDPG.L1Menu.customL1Ntuple_cfg import *

process.p.remove(process.l1RecoTreeProducer)
process.p.remove(process.l1MuonRecoTreeProducer)

if options.useStage1Layer2:
    process.p *= process.Layer2
    process.p *= process.l1NtupleProducerStage1Layer2

# edit here

OUTFILE="L1Tree_both.root"
NEVTS=-1

process.TFileService.fileName=OUTFILE
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(NEVTS) )

## readFiles.extend( ['file:///data2/battilan/L1Trigger/62X_RAW_RECO.root'] )
readFiles.extend( ['root://eoscms.cern.ch//eos/cms/store/group/comm_trigger/L1Trigger/apana/262AA156-744A-E311-9829-002618943945.root'] )

## processDumpFile = open('config.dump', 'w')
## print >> processDumpFile, process.dumpPython()

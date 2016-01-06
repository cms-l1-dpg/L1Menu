import FWCore.ParameterSet.Config as cms

process = cms.Process("L1NTUPLE")

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.load("L1TriggerDPG.L1Ntuples.l1ExtraTreeProducer_cfi")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

# output file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('L1Tree.root')
)

# L1 ntuple producers
## process.load("L1TriggerDPG.L1Ntuples.l1NtupleProducer_cfi")
import L1TriggerDPG.L1Ntuples.l1NtupleProducerNano_cfi 
process.l1NtupleProducer = L1TriggerDPG.L1Ntuples.l1NtupleProducerNano_cfi.l1NtupleProducer.clone()

process.p = cms.Path(
    process.l1NtupleProducer
    +process.l1extraParticles
    +process.l1ExtraTreeProducer
    +process.l1GtTriggerMenuLite
)

# job options
process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v1'

SkipEvent = cms.untracked.vstring('ProductNotFound')

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

process.source = cms.Source ("PoolSource",
                             fileNames = readFiles,
                             secondaryFileNames = secFiles
                             )

# edit here
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

## readFiles.extend( ['file:///data2/battilan/L1Trigger/62X_RAW_RECO.root'] )
## readFiles.extend( ['file:/uscmst1b_scratch/lpc1/lpctrig/apana/L1Upgrade/262AA156-744A-E311-9829-002618943945.root'] )
## readFiles.extend( ['root://eoscms.cern.ch//eos/cms/store/group/comm_trigger/L1Trigger/apana/DYJetsToLL_M-50_13TeV-pythia6_Fall13dr-tsg_PU40bx25__skim_150_1_6UN.root'] )
## readFiles.extend( ['root://xrootd.unl.edu//store/data/Commissioning2014/Cosmics/RAW/v3/000/228/929/00000/509C6C80-3164-E411-92E0-02163E00FFE1.root'] )

readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2015C/L1Accept/RAW/v1/000/254/790/00000/16005522-3448-E511-AAFE-02163E011E5B.root'] )

## processDumpFile = open('config.dump', 'w')
## print >> processDumpFile, process.dumpPython()

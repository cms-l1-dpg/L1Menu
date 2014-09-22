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
## readFiles.extend( ['root://eoscms.cern.ch//eos/cms/store/group/comm_trigger/L1Trigger/apana/DYJetsToLL_M-50_13TeV-pythia6_Fall13dr-tsg_PU40bx25__skim_150_1_6UN.root'] )

readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/003B1379-5E00-E211-8542-001D09F24024.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/0209ABF5-7900-E211-8F66-BCAEC518FF44.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/067305B3-6200-E211-9F2B-001D09F2AD4D.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/083E64A6-4300-E211-BA44-BCAEC518FF30.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/0E941746-4900-E211-A8B8-0025901D631E.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/20A5C502-3E00-E211-812D-003048F11C58.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/244682E2-4700-E211-9B81-5404A63886C5.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/281C6EE7-7200-E211-B605-001D09F295FB.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/32B1FD3E-4200-E211-814E-E0CB4E55367F.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/3A310F0F-5100-E211-A62A-0030486730C6.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/3A6BEA28-7700-E211-B182-485B3962633D.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/3ACEA6F3-5A00-E211-A3E1-003048D37694.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/465B3737-6600-E211-B28D-001D09F27067.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/50BAA9A4-4F00-E211-9652-0015C5FDE067.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/50D2F56A-5700-E211-9636-5404A63886AE.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/58CD5882-5200-E211-BCE3-5404A6388698.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/5E94CDB0-6E00-E211-9D87-001D09F2924F.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/64F450D2-4000-E211-AB1E-BCAEC532971D.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/66F40371-3300-E211-ADFB-003048D37538.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/68FDABD0-7000-E211-A677-001D09F2B30B.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/763FCA73-4600-E211-84B9-BCAEC518FF65.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/78825086-5900-E211-A07D-001D09F24682.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/78A0D6A9-3700-E211-97C7-BCAEC5329717.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/7C17DE85-4C00-E211-A439-001D09F27003.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/7CE72DDD-3400-E211-A6ED-003048D374CA.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/8812B896-3C00-E211-BF97-001D09F24353.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/88B4757F-3A00-E211-894A-003048CF9B28.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/8E9ABC4B-3100-E211-88F7-BCAEC53296F3.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/9059595C-6800-E211-B97F-001D09F24399.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/90C57231-3B00-E211-ABFC-0025901D627C.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/9CD61C10-4500-E211-AD76-BCAEC518FF80.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/A4751591-6000-E211-A2F2-001D09F242EF.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/A6929C89-3500-E211-B7F0-5404A638868F.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/AE02DFD1-6400-E211-B12C-001D09F242EF.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/AEF03DBD-4A00-E211-8110-001D09F2905B.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/B4B8609B-3000-E211-8EF2-0030486780B8.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/C4BEE117-7C00-E211-8434-001D09F2932B.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/CA1A0684-4D00-E211-8A83-003048F117EC.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/CA4B0995-6C00-E211-9FBC-001D09F241B9.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/D0F6A9DE-5300-E211-A99C-003048D2BD66.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/D62CBA74-3F00-E211-BA10-00215AEDFD74.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/E07B8F34-7E00-E211-92C3-001D09F2A690.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/EECBEC00-5600-E211-8537-003048D2BB58.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/F0102C56-5C00-E211-A062-001D09F28EA3.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/F2FD4E63-3800-E211-97E2-E0CB4E553673.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/F863BD14-3900-E211-9EC8-001D09F2447F.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/F8914E89-8200-E211-A3A7-003048D3733E.root'] )
readFiles.extend( ['root://xrootd.unl.edu//store/data/Run2012C/MinimumBias/RAW/v1/000/203/002/FEB5BC0D-7500-E211-89FE-001D09F26509.root'] )



## processDumpFile = open('config.dump', 'w')
## print >> processDumpFile, process.dumpPython()

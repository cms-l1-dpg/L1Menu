import FWCore.ParameterSet.Config as cms

def pileUpReweighting(process, fileName, origHisto, targetHisto ):

    if hasattr(process,'l1NtupleProducer'):
        print "[L1Menu]: Customising PU reweighting into ntuple with file" , fileName
 
        ntuple = getattr(process,'l1NtupleProducer')
        ntuple.puMCFile   = cms.untracked.string(fileName)
        ntuple.puDataFile = cms.untracked.string(fileName)
        ntuple.puMCHist   = cms.untracked.string(origHisto)
        ntuple.puDataHist = cms.untracked.string(targetHisto)

    else :
        print "[L1Menu]: Ntuple configuration not found. Can't customise PU reweighting!"

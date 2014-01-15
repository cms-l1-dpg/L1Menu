import FWCore.ParameterSet.Config as cms

def customiseUCT2015(process, runOnMC):

    if hasattr(process,'reEmulCaloChain') :
        print "[L1Menu]: Customising calo chain with UCT2015"

        if runOnMC :
            print "[L1Menu]:\tUsing MC configuration"
            process.load("L1Trigger.UCT2015.emulationMC_cfi")
        else :
            print "[L1Menu]:\tUsing DATA configuration"
            process.load("L1Trigger.UCT2015.emulation_cfi") # For running on data

        process.uctGctDigis =cms.EDProducer("UCT2015GctCandsProducer",
                                           egRelaxed = cms.InputTag("UCT2015Producer","RelaxedEGUnpacked"),
                                           egIsolated = cms.InputTag("UCT2015Producer","IsolatedEGUnpacked"),
                                           tauRelaxed = cms.InputTag("UCT2015Producer","RelaxedTauUnpacked"), # this collection is ignored in the final output, GT constraints
                                           tauIsolated = cms.InputTag("UCT2015Producer","IsolatedTauUnpacked"),
                                           jetSource = cms.InputTag("UCT2015Producer","CorrJetUnpacked"), # default are corrected jets
                                           # jetSource = cms.InputTag("UCT2015Producer","JetUnpacked"),
                                           setSource = cms.InputTag("UCT2015Producer","SETUnpacked"),
                                           metSource = cms.InputTag("UCT2015Producer","METUnpacked"),
                                           shtSource = cms.InputTag("UCT2015Producer","SHTUnpacked"),
                                           mhtSource = cms.InputTag("UCT2015Producer","MHTUnpacked")
        )

        process.reEmulUctChain = cms.Sequence(process.emulationSequence + process.uctGctDigis)

        getattr(process,'reEmul').replace(process.reEmulCaloChain, process.reEmulUctChain)
    
        getattr(process,'gtReEmulDigis').GctInputTag = cms.InputTag("uctGctDigis")
        getattr(process,'gtReEmulDigis').EmulateBxInEvent = cms.int32(1)

    else :
       print "[L1Menu]: ERROR: Can't customise calo chain with UCT2015, reEmulCaloChain is missing!"

    if hasattr(process,'l1NtupleProducer') :
        print "[L1Menu]:\tConfiguring Ntuple to use UCT2015 information"
 
        ntuple = getattr(process,'l1NtupleProducer')
        ntuple.gctCentralJetsSource = cms.InputTag("uctGctDigis","cenJets")
        ntuple.gctNonIsoEmSource    = cms.InputTag("uctGctDigis","rlxEm")
        ntuple.gctForwardJetsSource = cms.InputTag("uctGctDigis","forJets")
        ntuple.gctIsoEmSource       = cms.InputTag("uctGctDigis","isoEm")
        ntuple.gctEnergySumsSource  = cms.InputTag("uctGctDigis","")
        ntuple.gctTauJetsSource     = cms.InputTag("uctGctDigis","isoTau")
        ntuple.rctSource            = cms.InputTag("none")

def customiseL1Calos(process, customGCT=True):

    print "[L1Menu]: Customising legacy calo chain with possible 2015 improvements"

    if customGCT and hasattr(process,"gctReEmulDigis") :
        
        print "[L1Menu]:\tCustomising GCT configuration to use 10 GeV jet Seeds"
        
        process.load("L1TriggerConfig.GctConfigProducers.L1GctConfig_cff")

        process.L1GctConfigProducers.JetFinderCentralJetSeed = 10.0
        process.L1GctConfigProducers.JetFinderForwardJetSeed = 10.0
        
        process.es_prefer_gct = cms.ESPrefer("L1GctConfigProducers")
        

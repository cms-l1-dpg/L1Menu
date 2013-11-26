import FWCore.ParameterSet.Config as cms

def customiseUCT2015(process, runOnMC):

    if hasattr(process,'reEmulCaloChain') :
        print "[L1TMenu]: Customising calo chain with UCT2015"

        if runOnMC :
            print "[L1TMenu]:\tUsing MC configuration"
            process.load("L1Trigger.UCT2015.emulationMC_cfi")
        else :
            print "[L1TMenu]:\tUsing DATA configuration"
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

        process.uctReEmulDigis = cms.Sequence(process.emulationSequence + process.uctGctDigis)

        getattr(process,'reEmul').replace(process.reEmulCaloChain, process.uctReEmulDigis)
    
        getattr(process,'gtReEmulDigis').GctInputTag = cms.InputTag("uctGctDigis")
        getattr(process,'gtReEmulDigis').EmulateBxInEvent = cms.int32(1)

    if hasattr(process,'l1NtupleProducer') :
        print "[L1TMenu]: Configuring Ntuple to use UCT2015 information"
 
        ntuple = getattr(process,'l1NtupleProducer')
        ntuple.gctCentralJetsSource = cms.InputTag("none")
        ntuple.gctNonIsoEmSource    = cms.InputTag("none")
        ntuple.gctForwardJetsSource = cms.InputTag("none")
        ntuple.gctIsoEmSource       = cms.InputTag("none")
        ntuple.gctEnergySumsSource  = cms.InputTag("none")
        ntuple.gctTauJetsSource     = cms.InputTag("none")
        ntuple.rctSource            = cms.InputTag("none")

    else :
       print "[L1TMenu]: ERROR: Can't customise calo chain with UCT2015, reEmulCaloChain is missing!"


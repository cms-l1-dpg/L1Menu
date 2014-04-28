import FWCore.ParameterSet.Config as cms

def customiseUCT2015(process, runOnMC, runOnPostLS1, whichPU ):

    if hasattr(process,'reEmulCaloChain') :
        print "[L1Menu]: Customising calo chain with UCT2015"

        if runOnMC and runOnPostLS1 :
            print "[L1Menu]:\tUsing MC configuration for post LS1"
            process.load("L1Trigger.UCT2015.emulationMC_cfi")
            from L1Trigger.UCT2015.regionSF_cfi import *
            if whichPU == 20 :
                process.CorrectedDigis.regionSubtraction = regionSubtraction_PU20_MC13TeV
        elif not runOnMC : 
            print "[L1Menu]:\tUsing DATA configuration"
            process.load("L1Trigger.UCT2015.emulation_cfi") # For running on data
        else :
            # If run on MC and on 53X needs a patch different w.r.t. the UCT "standard" one
            # in this case there are no hcal digis stored so need to run on unpacked (Zero Suppressed) ones
            print "[L1Menu]:\tUsing MC configuration for 53X"

            process.load("L1Trigger.UCT2015.emulation_cfi") # For running on data patch by hand
            
            process.load("SimCalorimetry.HcalSimProducers.hcalUnsuppressedDigis_cfi")
            process.load("SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff")

            process.hcalReEmulDigis = process.simHcalTriggerPrimitiveDigis.clone()
            
            process.hcalReEmulDigis.inputLabel = cms.VInputTag(cms.InputTag('hcalDigis'), cms.InputTag('hcalDigis'))
            process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(True)

            process.hackHCALMIPs.src = cms.InputTag("hcalReEmulDigis")


        process.uctGctDigis =cms.EDProducer("UCT2015GctCandsProducer",
                                           egRelaxed = cms.InputTag("UCT2015Producer","RelaxedEGUnpacked"),
                                           egIsolated = cms.InputTag("UCT2015Producer","IsolatedEGUnpacked"),
                                           tauRelaxed = cms.InputTag("UCT2015Producer","RelaxedTauUnpacked"), # this collection is ignored in the final output, GT constraints
                                           tauIsolated = cms.InputTag("UCT2015Producer","IsolatedTauUnpacked"),
                                           jetSource = cms.InputTag("UCT2015Producer","CorrJetUnpacked"), # default are corrected jets
                                           #jetSource = cms.InputTag("UCT2015Producer","JetUnpacked"),
                                           setSource = cms.InputTag("UCT2015Producer","SETUnpacked"),
                                           metSource = cms.InputTag("UCT2015Producer","METUnpacked"),
                                           shtSource = cms.InputTag("UCT2015Producer","SHTUnpacked"),
                                           mhtSource = cms.InputTag("UCT2015Producer","MHTUnpacked")
        )
        

        if runOnMC and not runOnPostLS1 :
            # If run on MC and on 53X needs a patch different w.r.t. the UCT "standard" one
            # in this case there are no hcal digis stored so need to run on unpacked (Zero Suppressed) ones
            
            process.reEmulUctChain = cms.Sequence(process.hcalReEmulDigis + process.emulationSequence + process.uctGctDigis)
        else :
            process.reEmulUctChain = cms.Sequence(process.emulationSequence + process.uctGctDigis)


        getattr(process,'reEmul').replace(process.reEmulCaloChain, process.reEmulUctChain)

        l1ExtraReEmul = getattr(process,'l1ExtraReEmul') 

        l1ExtraReEmul.isolatedEmSource    = cms.InputTag("uctGctDigis","isoEm")
        l1ExtraReEmul.nonIsolatedEmSource = cms.InputTag("uctGctDigis","nonIsoEm")

        l1ExtraReEmul.forwardJetSource = cms.InputTag("uctGctDigis","forJets")
        l1ExtraReEmul.centralJetSource = cms.InputTag("uctGctDigis","cenJets")
        l1ExtraReEmul.tauJetSource     = cms.InputTag("uctGctDigis","tauJets")
            
        l1ExtraReEmul.etTotalSource = cms.InputTag("uctGctDigis")
        l1ExtraReEmul.etHadSource   = cms.InputTag("uctGctDigis")
        l1ExtraReEmul.etMissSource  = cms.InputTag("uctGctDigis")
        l1ExtraReEmul.htMissSource  = cms.InputTag("uctGctDigis")

        l1ExtraReEmul.hfRingEtSumsSource    = cms.InputTag("uctGctDigis")
        l1ExtraReEmul.hfRingBitCountsSource = cms.InputTag("uctGctDigis")

        getattr(process,'gtReEmulDigis').GctInputTag = cms.InputTag("uctGctDigis")
        getattr(process,'gtReEmulDigis').EmulateBxInEvent = cms.int32(1)

    else :
       print "[L1Menu]: ERROR: Can't customise calo chain with UCT2015, reEmulCaloChain is missing!"

    if hasattr(process,'l1NtupleProducer') and hasattr(process,'l1ExtraTreeProducer') :
        print "[L1Menu]:\tConfiguring Ntuple to use UCT2015 information"
 
        ntuple = getattr(process,'l1NtupleProducer')
        ntuple.gctCentralJetsSource = cms.InputTag("uctGctDigis","cenJets")
        ntuple.gctNonIsoEmSource    = cms.InputTag("uctGctDigis","nonIsoEm")
        ntuple.gctForwardJetsSource = cms.InputTag("uctGctDigis","forJets")
        ntuple.gctIsoEmSource       = cms.InputTag("uctGctDigis","isoEm")
        ntuple.gctTauJetsSource     = cms.InputTag("uctGctDigis","tauJets")
        ntuple.gctEnergySumsSource  = cms.InputTag("uctGctDigis","")
        ntuple.rctSource            = cms.InputTag("none")



def customiseL1Calos(process, customGCT=True):

    print "[L1Menu]: Customising legacy calo chain with possible 2015 improvements"

    if customGCT and hasattr(process,"gctReEmulDigis") :

        print "[L1Menu]:\tCustomising GCT configuration to use 10 GeV jet Seeds"

        process.load("L1TriggerConfig.GctConfigProducers.L1GctConfig_cff")
    
        process.L1GctConfigProducers.JetFinderCentralJetSeed = 10.0
        process.L1GctConfigProducers.JetFinderForwardJetSeed = 10.0
    
        process.es_prefer_gct = cms.ESPrefer("L1GctConfigProducers")
        

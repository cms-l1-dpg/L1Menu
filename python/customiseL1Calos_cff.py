import FWCore.ParameterSet.Config as cms

def updatel1ExtraReEmulTag(process,inputTag):

    l1ExtraReEmul = getattr(process,'l1ExtraReEmul') 

    l1ExtraReEmul.isolatedEmSource    = cms.InputTag(inputTag,"isoEm")
    l1ExtraReEmul.nonIsolatedEmSource = cms.InputTag(inputTag,"nonIsoEm")

    l1ExtraReEmul.forwardJetSource = cms.InputTag(inputTag,"forJets")
    l1ExtraReEmul.centralJetSource = cms.InputTag(inputTag,"cenJets")
    l1ExtraReEmul.tauJetSource     = cms.InputTag(inputTag,"tauJets")
    l1ExtraReEmul.isoTauJetSource  = cms.InputTag(inputTag,"isoTauJets")
        
    l1ExtraReEmul.etTotalSource = cms.InputTag(inputTag)
    l1ExtraReEmul.etHadSource   = cms.InputTag(inputTag)
    l1ExtraReEmul.etMissSource  = cms.InputTag(inputTag)
    l1ExtraReEmul.htMissSource  = cms.InputTag(inputTag)
    
    l1ExtraReEmul.hfRingEtSumsSource    = cms.InputTag(inputTag)
    l1ExtraReEmul.hfRingBitCountsSource = cms.InputTag(inputTag)

def updategtReEmulTag(process,inputTag):

    getattr(process,'gtReEmulDigis').GctInputTag = cms.InputTag(inputTag)
    getattr(process,'gtReEmulDigis').EmulateBxInEvent = cms.int32(1)

def updatel1ntupleTag(process,inputTag):

    ntuple = getattr(process,'l1NtupleProducer')
    ntuple.gctCentralJetsSource = cms.InputTag(inputTag,"cenJets")
    ntuple.gctNonIsoEmSource    = cms.InputTag(inputTag,"nonIsoEm")
    ntuple.gctForwardJetsSource = cms.InputTag(inputTag,"forJets")
    ntuple.gctIsoEmSource       = cms.InputTag(inputTag,"isoEm")
    ntuple.gctTauJetsSource     = cms.InputTag(inputTag,"tauJets")
    ntuple.gctIsoTauJetsSource  = cms.InputTag(inputTag,"isoTauJets")
    ntuple.gctEnergySumsSource  = cms.InputTag(inputTag,"")
    ntuple.rctSource            = cms.InputTag("none")

def set10GCTtreshold(process):

    if hasattr(process,"gctReEmulDigis") :

        print "[L1Menu]:\tCustomising GCT configuration to use 10 GeV jet Seeds"

        process.load("L1TriggerConfig.GctConfigProducers.L1GctConfig_cff")
        process.L1GctConfigProducers.JetFinderCentralJetSeed = 10.0
        process.L1GctConfigProducers.JetFinderForwardJetSeed = 10.0
        process.es_prefer_gct = cms.ESPrefer("L1GctConfigProducers")

def customiseStage1(process, runOnMC, runOnPostLS1, whichPU ):

    if hasattr(process,'reEmulCaloChain') :
        print "[L1Menu]: Customising calo chain with new L1 Stage1 Emulator"

        process.load('L1Trigger.L1TCalorimeter.L1TCaloStage1_PPFromRaw_cff')
        process.load('L1Trigger/L1TCalorimeter/caloStage1RegionSF_cfi')
        process.caloStage1Params.jetSeedThreshold = 5.0
        from L1Trigger.L1TCalorimeter.caloStage1RegionSF_cfi import regionSubtraction_PU40_MC13TeV
        from L1Trigger.L1TCalorimeter.caloStage1RegionSF_cfi import regionSubtraction_PU20_MC13TeV
        if runOnMC and whichPU == 20 :
            process.caloStage1Params.regionPUSParams = regionSubtraction_PU20_MC13TeV


        getattr(process,'reEmul').replace(process.reEmulCaloChain, process.L1TCaloStage1_PPFromRaw)

        ## unpack stage1 digis as well as gct digis
        from EventFilter.L1TRawToDigi.caloStage1Digis_cfi import caloStage1Digis
        process.caloStage1Digis = caloStage1Digis.clone()
        process.p.replace(process.gctDigis, process.gctDigis + process.caloStage1Digis)

        
        l1ExtraReEmul = getattr(process,'l1ExtraReEmul') 

        updatel1ExtraReEmulTag(process,"simCaloStage1LegacyFormatDigis")
        updategtReEmulTag(process,"simCaloStage1LegacyFormatDigis")

    else :
       print "[L1Menu]: ERROR: Can't customise calo chain with Stage1 emulator, reEmulCaloChain is missing!"

    if hasattr(process,'l1NtupleProducer') and hasattr(process,'l1ExtraTreeProducer') :
        print "[L1Menu]:\tConfiguring Ntuple to use Stage1 emulator information"
 
        updatel1ntupleTag(process,"simCaloStage1LegacyFormatDigis")

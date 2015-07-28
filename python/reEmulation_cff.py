import FWCore.ParameterSet.Config as cms

def reEmulation(process, reEmulMuons=True, reEmulCalos=True, patchNtuple=True, runOnPostLS1 = False, useStage1Layer2=False):

    print "[L1Menu]: Setting up overall re-emulation"        

    if patchNtuple and hasattr(process,'l1NtupleProducer') and hasattr(process,'l1ExtraTreeProducer') :
        print "[L1Menu]:\tConfiguring Ntuple to use re-emulated information"        
        ntuple        = getattr(process,'l1NtupleProducer')
        l1ExtraNtuple = getattr(process,'l1ExtraTreeProducer')
    elif patchNtuple :
        print "[L1Menu]:\tERROR: FAILED to find ntuple! switching patchNtuple to False"        
        patchNtuple=False

    process.l1ExtraReEmul = cms.EDProducer(
        "L1ExtraParticlesProd",
        muonSource = cms.InputTag("gtDigis"),
        isolatedEmSource    = cms.InputTag("gctDigis","isoEm"),
        nonIsolatedEmSource = cms.InputTag("gctDigis","nonIsoEm"),
        
        forwardJetSource = cms.InputTag("gctDigis","forJets"),
        centralJetSource = cms.InputTag("gctDigis","cenJets"),
        tauJetSource     = cms.InputTag("gctDigis","tauJets"),
        isoTauJetSource  = cms.InputTag("gctDigis","isoTauJets"),
        
        etTotalSource = cms.InputTag("gctDigis"),
        etHadSource   = cms.InputTag("gctDigis"),
        etMissSource  = cms.InputTag("gctDigis"),
        htMissSource  = cms.InputTag("gctDigis"),
        
        hfRingEtSumsSource    = cms.InputTag("gctDigis"),
        hfRingBitCountsSource = cms.InputTag("gctDigis"),
        
        produceMuonParticles = cms.bool(True),
        produceCaloParticles = cms.bool(True),
        centralBxOnly = cms.bool(True),
        ignoreHtMiss = cms.bool(False)
        )
        
    if reEmulMuons :
        print "[L1Menu]:\tSetting up muon re-emulation"        
        
        from L1Trigger.DTTrackFinder.dttfDigis_cfi import dttfDigis
        process.dttfReEmulDigis       = dttfDigis.clone()
        process.dttfReEmulDigis.DTDigi_Source  = cms.InputTag("dttfDigis")
        process.dttfReEmulDigis.CSCStub_Source = cms.InputTag("csctfReEmulTrackDigis")

        from L1Trigger.RPCTrigger.rpcTriggerDigis_cfi import rpcTriggerDigis
        process.rpcTriggerReEmulDigis = rpcTriggerDigis.clone()

        if not runOnPostLS1 :
            from L1Trigger.CSCTrackFinder.csctfTrackDigis_cfi import csctfTrackDigis
            from L1Trigger.CSCTrackFinder.csctfDigis_cfi import csctfDigis

        
            process.csctfReEmulTrackDigis = csctfTrackDigis.clone()
            process.csctfReEmulDigis      = csctfDigis.clone()
        
            process.csctfReEmulTrackDigis.readDtDirect        = True
            process.csctfReEmulTrackDigis.SectorReceiverInput = cms.untracked.InputTag("csctfDigis")
            process.csctfReEmulTrackDigis.DtDirectProd        = cms.untracked.InputTag("csctfDigis","DT")
            process.csctfReEmulDigis.CSCTrackProducer         = cms.untracked.InputTag("csctfReEmulTrackDigis")
            #process.csctfReEmulDigis.SectorProcessor.initializeFromPSet = True
            
            process.csctfReEmulSequence = cms.Sequence(
                process.csctfReEmulTrackDigis
                * process.csctfReEmulDigis
            )
        else :
            from SLHCUpgradeSimulations.Configuration.muonCustoms import customise_csc_L1Emulator_sim
            from L1Trigger.CSCTrackFinder.csctfDigis_cfi import csctfDigis

            customise_csc_L1Emulator_sim(process) 

            process.csctfReEmulTrackDigis = process.simCsctfTrackDigis.clone()
            process.csctfReEmulDigis      = csctfDigis.clone()

            process.csctfReEmulTrackDigis.DTproducer  = cms.untracked.InputTag("dttfDigis")
            process.csctfReEmulDigis.CSCTrackProducer = cms.untracked.InputTag("csctfReEmulTrackDigis")

            process.csctfReEmulTrackDigis.SectorProcessor.PTLUT.PtMethod = cms.untracked.uint32(34) # no triple ganging in ME11a
            process.csctfReEmulTrackDigis.SectorProcessor.gangedME1a = cms.untracked.bool(False)
            process.csctfReEmulTrackDigis.SectorProcessor.firmwareSP = cms.uint32(20140515) #core 20120730
            process.csctfReEmulTrackDigis.SectorProcessor.initializeFromPSet = cms.bool(True) 

            process.csctfReEmulSequence = cms.Sequence(
                process.simCscTriggerPrimitiveDigis
                * process.csctfReEmulTrackDigis
                * process.csctfReEmulDigis
            )

            process.load('L1TriggerConfig.GMTConfigProducers.L1MuGMTParameters_cfi')
            process.L1MuGMTParameters.MergeMethodPtBrl=cms.string("byCombi")
            process.L1MuGMTParameters.MergeMethodPtFwd=cms.string("byCombi")
            process.L1MuGMTParameters.VersionSortRankEtaQLUT = cms.uint32(1043)
            process.L1MuGMTParameters.VersionLUTs = cms.uint32(1)


        from L1Trigger.GlobalMuonTrigger.gmtDigis_cfi import gmtDigis
        process.gmtReEmulDigis  = gmtDigis.clone()

        process.gmtReEmulDigis.DTCandidates   = cms.InputTag("dttfReEmulDigis","DT")
        process.gmtReEmulDigis.CSCCandidates  = cms.InputTag("csctfReEmulDigis","CSC")
        process.gmtReEmulDigis.RPCbCandidates = cms.InputTag("rpcTriggerReEmulDigis","RPCb")
        process.gmtReEmulDigis.RPCfCandidates = cms.InputTag("rpcTriggerReEmulDigis","RPCf")
        process.gmtReEmulDigis.MipIsoData     = cms.InputTag("none")
        
        process.l1ExtraReEmul.muonSource = cms.InputTag("gmtReEmulDigis")            

        if patchNtuple :
            ntuple.gmtSource          = cms.InputTag("gmtReEmulDigis")
            ntuple.dttfSource         = cms.InputTag("dttfReEmulDigis")
            ntuple.csctfTrkSource     = cms.InputTag("csctfReEmulTrackDigis")
            ntuple.csctfStatusSource  = cms.InputTag("csctfReEmulTrackDigis")

            l1ExtraNtuple.muonLabel = cms.untracked.InputTag("l1ExtraReEmul")

        process.reEmulMuonChain = cms.Sequence(
            process.rpcTriggerReEmulDigis
            *process.csctfReEmulSequence
            *process.dttfReEmulDigis
            *process.gmtReEmulDigis
            )

    if reEmulCalos :
        print "[L1Menu]:\tSetting up calo re-emulation"        

        # Need to have RCT emulator configurable and not UCT 2015 patches
        # in order to run 2012 RCT emulator correctly
        
        #	# In MC HCAL need to be re-run as there is no TPG information stored
        #	process.load("SimCalorimetry.HcalSimProducers.hcalUnsuppressedDigis_cfi")
        #	process.load("SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff")
        #	
        #	from L1Trigger.RegionalCaloTrigger.rctDigis_cfi import rctDigis
        #	from L1Trigger.GlobalCaloTrigger.gctDigis_cfi import gctDigis
        #	
        #	process.hcalReEmulDigis = process.simHcalTriggerPrimitiveDigis.clone()
        #	process.rctReEmulDigis  = rctDigis.clone()
        #	process.gctReEmulDigis  = gctDigis.clone()
        #	
        #	process.hcalReEmulDigis.inputLabel = cms.VInputTag(cms.InputTag('hcalDigis'), cms.InputTag('hcalDigis'))
        #	#process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(False)
        #	
        #	process.rctReEmulDigis.ecalDigis = cms.VInputTag( cms.InputTag( 'ecalDigis:EcalTriggerPrimitives' ) )
        #	process.rctReEmulDigis.hcalDigis = cms.VInputTag( cms.InputTag( 'hcalReEmulDigis' ) )
        #	
        #	process.gctReEmulDigis.inputLabel  = cms.InputTag("rctReEmulDigis")

        from L1Trigger.Configuration.SimL1Emulator_cff import simRctDigis
        process.rctReEmulDigis = process.simRctDigis.clone()
        process.rctReEmulDigis.ecalDigis = cms.VInputTag( cms.InputTag( 'ecalDigis:EcalTriggerPrimitives' ) )
        process.rctReEmulDigis.hcalDigis = cms.VInputTag( cms.InputTag( 'hcalDigis' ) )

        from L1Trigger.GlobalCaloTrigger.gctDigis_cfi import gctDigis

        process.gctReEmulDigis = gctDigis.clone()        
        process.gctReEmulDigis.inputLabel = cms.InputTag("gctDigis")

        process.l1ExtraReEmul.isolatedEmSource    = cms.InputTag("gctReEmulDigis","isoEm")
        process.l1ExtraReEmul.nonIsolatedEmSource = cms.InputTag("gctReEmulDigis","nonIsoEm")

        process.l1ExtraReEmul.forwardJetSource = cms.InputTag("gctReEmulDigis","forJets")
        process.l1ExtraReEmul.centralJetSource = cms.InputTag("gctReEmulDigis","cenJets")
        process.l1ExtraReEmul.tauJetSource     = cms.InputTag("gctReEmulDigis","tauJets")
        process.l1ExtraReEmul.isoTauJetSource  = cms.InputTag("gctReEmulDigis","isoTauJets")            

        process.l1ExtraReEmul.etTotalSource = cms.InputTag("gctDigis")
        process.l1ExtraReEmul.etHadSource   = cms.InputTag("gctReEmulDigis")
        process.l1ExtraReEmul.etMissSource  = cms.InputTag("gctReEmulDigis")
        process.l1ExtraReEmul.htMissSource  = cms.InputTag("gctReEmulDigis")

        process.l1ExtraReEmul.hfRingEtSumsSource    = cms.InputTag("gctReEmulDigis")
        process.l1ExtraReEmul.hfRingBitCountsSource = cms.InputTag("gctReEmulDigis")

        if patchNtuple :
            ntuple.gctCentralJetsSource = cms.InputTag("gctReEmulDigis","cenJets")
            ntuple.gctNonIsoEmSource    = cms.InputTag("gctReEmulDigis","nonIsoEm")
            ntuple.gctForwardJetsSource = cms.InputTag("gctReEmulDigis","forJets")
            ntuple.gctIsoEmSource       = cms.InputTag("gctReEmulDigis","isoEm")
            ntuple.gctEnergySumsSource  = cms.InputTag("gctReEmulDigis","")
            ntuple.gctTauJetsSource     = cms.InputTag("gctReEmulDigis","tauJets")
            if useStage1Layer2:
                ntuple.gctIsoTauJetsSource  = cms.InputTag("gctReEmulDigis","isoTauJets")
            else:
                ntuple.gctIsoTauJetsSource  = cms.InputTag("none","isoTauJets")

            l1ExtraNtuple.nonIsoEmLabel = cms.untracked.InputTag("l1ExtraReEmul:NonIsolated")
            l1ExtraNtuple.isoEmLabel    = cms.untracked.InputTag("l1ExtraReEmul:Isolated")
            l1ExtraNtuple.tauJetLabel   = cms.untracked.InputTag("l1ExtraReEmul:Tau")
            l1ExtraNtuple.isoTauJetLabel   = cms.untracked.InputTag("l1ExtraReEmul:IsoTau")
            l1ExtraNtuple.cenJetLabel   = cms.untracked.InputTag("l1ExtraReEmul:Central")
            l1ExtraNtuple.fwdJetLabel   = cms.untracked.InputTag("l1ExtraReEmul:Forward")
            l1ExtraNtuple.metLabel      = cms.untracked.InputTag("l1ExtraReEmul:MET")
            l1ExtraNtuple.mhtLabel      = cms.untracked.InputTag("l1ExtraReEmul:MHT")

            # Need to have RCT emulator configurable and not UCT 2015 patches
            # in order to run 2012 RCT emulator correctly
            #    ntuple.rctSource            = cms.InputTag("rctReEmulDigis")

        process.reEmulCaloChain = cms.Sequence(
            #process.hcalReEmulDigis
            process.rctReEmulDigis
            +process.gctReEmulDigis
        )


    from L1Trigger.GlobalTrigger.gtDigis_cfi import gtDigis
    process.gtReEmulDigis   = gtDigis.clone()

    if reEmulMuons :
        process.gtReEmulDigis.GmtInputTag  = cms.InputTag("gmtReEmulDigis")
    if reEmulCalos :
        process.gctReEmulDigis.inputLabel  = cms.InputTag("rctReEmulDigis")
        process.gtReEmulDigis.GctInputTag  = cms.InputTag("gctReEmulDigis")

    if patchNtuple :
        ntuple.gtSource = cms.InputTag("gtReEmulDigis")
        
    if reEmulMuons and reEmulCalos :
        process.reEmul = cms.Sequence(process.reEmulCaloChain + process.reEmulMuonChain + process.gtReEmulDigis + process.l1ExtraReEmul)
    elif reEmulMuons :
        process.reEmul = cms.Sequence(process.reEmulMuonChain + process.gtReEmulDigis + process.l1ExtraReEmul)
    elif reEmulCalos :
        process.reEmul = cms.Sequence(process.reEmulCaloChain + process.gtReEmulDigis + process.l1ExtraReEmul)
    else :
        process.reEmul = cms.Sequence(process.gtReEmulDigis + process.l1ExtraReEmul)


def run2012CConfiguration(process):

        print "[L1Menu]: Setting up muon/calo configuration to correspond to RUN2012C"
        print "[L1Menu]: Forcing RPC muon trigger to switch off HSCP BX extension"

        process.GlobalTag.toGet.extend(
            cms.VPSet(cms.PSet(record = cms.string("L1GtTriggerMenuRcd"),
                               tag = cms.string("L1GtTriggerMenu_L1Menu_Collisions2012_v2_mc"),
                               connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T")
                               ),
                      cms.PSet(record = cms.string("L1GctChannelMaskRcd"),
                               tag = cms.string("L1GctChannelMask_AllEnergySumsMaskedFromHF_jetCentresToEta3Allowed_mc"),
                               connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T")
                               ),
                      cms.PSet(record = cms.string("L1GctJetFinderParamsRcd"),
                               tag = cms.string("L1GctJetFinderParams_GCTPhysics_2012_04_27_JetSeedThresh5GeV_mc"),
                               connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T")
                               ),
                      cms.PSet(record = cms.string("L1HfRingEtScaleRcd"),
                               tag = cms.string("L1HfRingEtScale_GCTPhysics_2012_04_27_JetSeedThresh5GeV_mc"),
                               connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T")
                               ),
                      cms.PSet(record = cms.string("L1JetEtScaleRcd"),
                               tag = cms.string("L1JetEtScale_GCTPhysics_2012_04_27_JetSeedThresh5GeV_mc"),
                               connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T")
                               ),
                      cms.PSet(record = cms.string("L1HtMissScaleRcd"),
                               tag = cms.string("L1HtMissScale_GCTPhysics_2012_04_27_JetSeedThresh5GeV_mc"),
                               connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T")
                               ),
                      cms.PSet(record = cms.string("L1MuCSCPtLutRcd"),
                               tag = cms.string("L1MuCSCPtLut_key-11_mc"),
                               connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T")
                               ),
                      cms.PSet(record = cms.string("L1MuDTTFParametersRcd"),
                               tag = cms.string("L1MuDTTFParameters_dttf12_TSC_03_csc_col_mc"),
                               connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T")
                               ),
                      # Forcing RPCs without HSCP (also in 50 ns as it has minor impact in rates there)
                      cms.PSet(record = cms.string("L1RPCBxOrConfigRcd"),
                               tag = cms.string("L1RPCBxOrConfig_LHC8_mc"),
                               connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T")
                               ),
                      cms.PSet(record = cms.string("L1RPCConeDefinitionRcd"),
                               tag = cms.string("L1RPCConeDefinition_LHC8_mc"),
                               connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T")
                               ),
                      cms.PSet(record = cms.string("L1RPCConfigRcd"),
                               tag = cms.string("L1RPCConfig_LHC8_mc"),
                               connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T")
                               ),
                      cms.PSet(record = cms.string("L1RPCHsbConfigRcd"),
                               tag = cms.string("L1RPCHsbConfig_LHC8_mc"),
                               connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T")
                               )
                      )
            )


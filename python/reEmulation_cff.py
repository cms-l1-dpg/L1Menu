import FWCore.ParameterSet.Config as cms

def reEmulation(process, reEmulMuons=True, reEmulCalos=True, patchNtuple=True):

    print "[L1Menu]: Setting up overall re-emulation"        

    if patchNtuple and hasattr(process,'l1NtupleProducer') :
        print "[L1Menu]:\tConfiguring Ntuple to use re-emulated information"        
        ntuple = getattr(process,'l1NtupleProducer')
    elif patchNtuple :
        print "[L1Menu]:\tERROR: FAILED to find ntuple! switching patchNtuple to False"        
        patchNtuple=False

    if reEmulMuons :
        print "[L1Menu]:\tSetting up muon re-emulation"        
        
        from L1Trigger.DTTrackFinder.dttfDigis_cfi import dttfDigis
        from L1Trigger.CSCTrackFinder.csctfTrackDigis_cfi import csctfTrackDigis
        from L1Trigger.CSCTrackFinder.csctfDigis_cfi import csctfDigis
        from L1Trigger.RPCTrigger.rpcTriggerDigis_cfi import rpcTriggerDigis
        from L1Trigger.GlobalMuonTrigger.gmtDigis_cfi import gmtDigis


        process.dttfReEmulDigis       = dttfDigis.clone()
        process.csctfReEmulTrackDigis = csctfTrackDigis.clone()
        process.rpcTriggerReEmulDigis = rpcTriggerDigis.clone()
        process.csctfReEmulDigis      = csctfDigis.clone()
        
        process.gmtReEmulDigis  = gmtDigis.clone()

        process.dttfReEmulDigis.DTDigi_Source  = cms.InputTag("dttfDigis")
        process.dttfReEmulDigis.CSCStub_Source = cms.InputTag("csctfReEmulTrackDigis")

        process.csctfReEmulTrackDigis.readDtDirect        = True
        process.csctfReEmulTrackDigis.SectorReceiverInput = cms.untracked.InputTag("csctfDigis")
        process.csctfReEmulTrackDigis.DtDirectProd        = cms.untracked.InputTag("csctfDigis","DT")
        process.csctfReEmulDigis.CSCTrackProducer         = cms.untracked.InputTag("csctfReEmulTrackDigis")
        #process.csctfReEmulDigis.SectorProcessor.initializeFromPSet = True

        process.gmtReEmulDigis.DTCandidates   = cms.InputTag("dttfReEmulDigis","DT")
        process.gmtReEmulDigis.CSCCandidates  = cms.InputTag("csctfReEmulDigis","CSC")
        process.gmtReEmulDigis.RPCbCandidates = cms.InputTag("rpcTriggerReEmulDigis","RPCb")
        process.gmtReEmulDigis.RPCfCandidates = cms.InputTag("rpcTriggerReEmulDigis","RPCf")
        process.gmtReEmulDigis.MipIsoData     = cms.InputTag("none")
        

        if patchNtuple :
            ntuple.gmtSource          = cms.InputTag("gmtReEmulDigis")
            ntuple.dttfSource         = cms.InputTag("dttfReEmulDigis")
            ntuple.csctfTrkSource     = cms.InputTag("csctfReEmulTrackDigis")
            ntuple.csctfStatusSource  = cms.InputTag("csctfReEmulTrackDigis")

        process.reEmulMuonChain = cms.Sequence(
            process.rpcTriggerReEmulDigis
            *process.csctfReEmulTrackDigis
            *process.csctfReEmulDigis
            *process.dttfReEmulDigis
            *process.gmtReEmulDigis
            )

    if reEmulCalos :
        print "[L1Menu]:\tSetting up calo re-emulation"        

        # In MC HCAL need to be re-run as there is no TPG information stored
        process.load("SimCalorimetry.HcalSimProducers.hcalUnsuppressedDigis_cfi")
        process.load("SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff")
        
        from L1Trigger.RegionalCaloTrigger.rctDigis_cfi import rctDigis
        from L1Trigger.GlobalCaloTrigger.gctDigis_cfi import gctDigis
        from L1Trigger.GlobalTrigger.gtDigis_cfi import gtDigis    

        process.hcalReEmulDigis = process.simHcalTriggerPrimitiveDigis.clone()
        process.rctReEmulDigis  = rctDigis.clone()
        process.gctReEmulDigis  = gctDigis.clone()
    
        process.hcalReEmulDigis.inputLabel = cms.VInputTag(cms.InputTag('hcalDigis'), cms.InputTag('hcalDigis'))
        #process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(False)

        process.rctReEmulDigis.ecalDigis = cms.VInputTag( cms.InputTag( 'ecalDigis:EcalTriggerPrimitives' ) )
        process.rctReEmulDigis.hcalDigis = cms.VInputTag( cms.InputTag( 'hcalReEmulDigis' ) )

        process.gctReEmulDigis.inputLabel  = cms.InputTag("rctReEmulDigis")
    
        if patchNtuple :
            ntuple.gctCentralJetsSource = cms.InputTag("gctReEmulDigis","cenJets")
            ntuple.gctNonIsoEmSource    = cms.InputTag("gctReEmulDigis","nonIsoEm")
            ntuple.gctForwardJetsSource = cms.InputTag("gctReEmulDigis","forJets")
            ntuple.gctIsoEmSource       = cms.InputTag("gctReEmulDigis","isoEm")
            ntuple.gctEnergySumsSource  = cms.InputTag("gctReEmulDigis","")
            ntuple.gctTauJetsSource     = cms.InputTag("gctReEmulDigis","tauJets")
            ntuple.rctSource            = cms.InputTag("rctReEmulDigis")

        process.reEmulCaloChain = cms.Sequence(
            process.hcalReEmulDigis
            + process.rctReEmulDigis
            + process.gctReEmulDigis
        )


    from L1Trigger.GlobalTrigger.gtDigis_cfi import gtDigis
    process.gtReEmulDigis   = gtDigis.clone()


    if reEmulMuons :
        process.gtReEmulDigis.GmtInputTag  = cms.InputTag("gmtReEmulDigis")
    if reEmulCalos :
        process.gtReEmulDigis.GctInputTag  = cms.InputTag("gctReEmulDigis")

    if patchNtuple :
        ntuple.gtSource = cms.InputTag("gtReEmulDigis")

    if reEmulMuons and reEmulCalos :
        process.reEmul = cms.Sequence(process.reEmulCaloChain + process.reEmulMuonChain + process.gtReEmulDigis)
    elif reEmulMuons :
        process.reEmul = cms.Sequence(process.reEmulMuonChain + process.gtReEmulDigis)
    elif reEmulCalos :
        process.reEmul = cms.Sequence(process.reEmulCaloChain + process.gtReEmulDigis)
    else :
        process.reEmul = cms.Sequence(process.gtReEmulDigis)


def run2012CConfiguration(process):

        print "[L1Menu]: Setting up muon/calo configuration to correspond to RUN2012C"
        print "[L1Menu]: Forcing RPC muon to switch off HSCP BX extension"

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


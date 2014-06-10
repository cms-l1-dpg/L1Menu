import FWCore.ParameterSet.Config as cms

def customiseL1Muons(process, customDTTF=True, customCSCTF=True, customPACT=True, customGMT=True, dttfFile = "sqlite_file:crab/dttf_config.db"):

    print "[L1Menu]: Customising muon chain with 2015 improvements"

    if customDTTF and hasattr(process,"dttfReEmulDigis") :
        
        print "[L1Menu]:\tCustomising DTTF LUTs"
        

        process.GlobalTag.toGet.extend(
            cms.VPSet(cms.PSet(record = cms.string("L1MuDTEtaPatternLutRcd"),
                               tag = cms.string("L1MuDTEtaPatternLut_CRAFT09_hlt"),
                               connect = cms.untracked.string(dttfFile)
                           ),
                      cms.PSet(record = cms.string("L1MuDTExtLutRcd"),
                               tag = cms.string("L1MuDTExtLut_CRAFT09_hlt"),
                               connect = cms.untracked.string(dttfFile)
                           ),
                      cms.PSet(record = cms.string("L1MuDTPhiLutRcd"),
                               tag = cms.string("L1MuDTPhiLut_CRAFT09_hlt"),
                               connect = cms.untracked.string(dttfFile)
                           ),
                      cms.PSet(record = cms.string("L1MuDTPtaLutRcd"),
                               tag = cms.string("L1MuDTPtaLut_CRAFT09_hlt"),
                               connect = cms.untracked.string(dttfFile)
                           ),
                      cms.PSet(record = cms.string("L1MuDTQualPatternLutRcd"),
                               tag = cms.string("L1MuDTQualPatternLut_CRAFT09_hlt"),
                               connect = cms.untracked.string(dttfFile)
                           )
                 )
        )

    if customPACT and hasattr(process,"rpcTriggerReEmulDigis") :

        print "[L1Menu]:\tCustomising PACT patterns"

        patternDirectory = "L1TriggerDPG/L1Menu/data/rpc_patterns/xml/"
        
        process.load("L1TriggerConfig.RPCTriggerConfig.RPCConeDefinition_cff")
        process.load("L1TriggerConfig.RPCTriggerConfig.L1RPCConfig_cff")
        process.load("L1Trigger.RPCTrigger.RPCConeConfig_cff")
        process.rpcconf.filedir = cms.untracked.string(patternDirectory)
        process.es_prefer_rpcPats = cms.ESPrefer("RPCTriggerConfig","rpcconf")

    if customGMT and hasattr(process,"gmtReEmulDigis") :

        print "[L1Menu]:\tCustomising GMT to use min-pt"

        process.load('L1TriggerConfig.GMTConfigProducers.L1MuGMTParameters_cfi')
        process.L1MuGMTParameters.MergeMethodPtBrl=cms.string("byMinPt")
        process.L1MuGMTParameters.MergeMethodPtFwd=cms.string("byMinPt")
        process.L1MuGMTParameters.VersionSortRankEtaQLUT = cms.uint32(275)
        process.L1MuGMTParameters.VersionLUTs = cms.uint32(1) 
        process.es_prefer_gmtConfig = cms.ESPrefer("L1MuGMTParametersProducer","L1MuGMTParameters")

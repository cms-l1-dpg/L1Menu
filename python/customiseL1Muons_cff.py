import FWCore.ParameterSet.Config as cms

def customiseL1Muons(process, customDTTF=True, customCSCTF=True, customPACT=True):

    if customDTTF and hasattr(process,"dttfReEmulDigis") :
        
        print "[L1TMenu]: Customising DTTF LUTs"
        
        connectFile = "sqlite_file:../data/dttf_config.db"

        process.GlobalTag.toGet.extend(
            cms.VPSet(cms.PSet(record = cms.string("L1MuDTEtaPatternLutRcd"),
                               tag = cms.string("L1MuDTEtaPatternLut_CRAFT09_hlt"),
                               connect = cms.untracked.string(connectFile)
                           ),
                      cms.PSet(record = cms.string("L1MuDTExtLutRcd"),
                               tag = cms.string("L1MuDTExtLut_CRAFT09_hlt"),
                               connect = cms.untracked.string(connectFile)
                           ),
                      cms.PSet(record = cms.string("L1MuDTPhiLutRcd"),
                               tag = cms.string("L1MuDTPhiLut_CRAFT09_hlt"),
                               connect = cms.untracked.string(connectFile)
                           ),
                      cms.PSet(record = cms.string("L1MuDTPtaLutRcd"),
                               tag = cms.string("L1MuDTPtaLut_CRAFT09_hlt"),
                               connect = cms.untracked.string(connectFile)
                           ),
                      cms.PSet(record = cms.string("L1MuDTQualPatternLutRcd"),
                               tag = cms.string("L1MuDTQualPatternLut_CRAFT09_hlt"),
                               connect = cms.untracked.string(connectFile)
                           )
                 )
        )

    if customPACT and hasattr(process,"rpcTriggerReEmulDigis") :

        print "[L1TMenu]: Customising PACT patterns"

        patternDirectory = "/L1TriggerDPG/L1TMenu/data/rpc_patterns/xml/"
        
        process.load("L1TriggerConfig.RPCTriggerConfig.RPCConeDefinition_cff")
        process.load("L1TriggerConfig.RPCTriggerConfig.L1RPCConfig_cff")
        process.load("L1Trigger.RPCTrigger.RPCConeConfig_cff")
        process.rpcconf.filedir = cms.untracked.string(patternDirectory)
        process.es_prefer_rpcPats = cms.ESPrefer("RPCTriggerConfig","rpcconf")

from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = '../customL1NtupleFromRaw.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['L1Tree.root']

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    from multiprocessing import Process

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    config.General.requestName = '13TeV_20PU_25ns_ReEmul2015_v16'
    config.Data.inputDataset = '/SingleNeutrino/RunIISpring15Digi74-AVE_20_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW'
    config.Data.outLFNDirBase = '/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2015/v16/'
    config.JobType.pyCfgParams = ['reEmulation=True', 'reEmulMuons=True', 'reEmulCalos=True', 'patchNtuple=True', 'useStage1Layer2=True', 'globalTag=MCRUN2_74_V9', 'runOnMC=True', 'runOnPostLS1=True', 'whichPU=20']
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = '13TeV_40PU_25ns_ReEmul2015_v16'
    config.Data.inputDataset = '/SingleNeutrino/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v2/GEN-SIM-RAW'
    config.Data.outLFNDirBase = '/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2015/v16/'
    config.JobType.pyCfgParams = ['reEmulation=True', 'reEmulMuons=True', 'reEmulCalos=True', 'patchNtuple=True', 'useStage1Layer2=True', 'globalTag=MCRUN2_74_V9', 'runOnMC=True', 'runOnPostLS1=True', 'whichPU=40']
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = '13TeV_30PU_50ns_ReEmul2012Gct10GeV_v16'
    config.Data.inputDataset = '/SingleNeutrino/RunIISpring15Digi74-AVE_30_BX_50ns_tsg_MCRUN2_74_V6-v1/GEN-SIM-RAW'
    config.Data.outLFNDirBase = '/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2015/v16/'
    config.JobType.pyCfgParams = ['reEmulation=True', 'reEmulMuons=True', 'reEmulCalos=True', 'patchNtuple=True', 'globalTag=MCRUN2_74_V8', 'runOnMC=True', 'runOnPostLS1=True', 'jetSeedThr10GeV=True']
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()


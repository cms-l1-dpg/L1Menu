from WMCore.Configuration import Configuration

config = Configuration()

RunOnMC = True

requestName = '251244'
pyCfg       = ['runOnMC=False', 'globalTag=GR_P_V56', 'reEmulation=True', 'patchNtuple=True', 'reEmulCalos=True', 'jetSeedThr10GeV=True', 'runOnPostLS1=True']
dataset     = '/ZeroBias/Run2015B-v1/RAW'
splitting   = 'LumiBased'
output      = '/store/group/dpg_trigger/comm_trigger/L1Trigger/Data/Collisions/'

if RunOnMC :
    requestName = 'Spring15_Flat10_50_50ns'
    pyCfg       = ['runOnMC=True', 'globalTag=MCRUN2_74_V8', 'reEmulation=True', 'patchNtuple=True', 'reEmulCalos=True', 'runOnPostLS1=True', 'jetSeedThr10GeV=True']
    dataset     = '/SingleNeutrino/RunIISpring15Digi74-Flat_10_50_50ns_tsg_MCRUN2_74_V6-v1/GEN-SIM-RAW'
    splitting   = 'FileBased'
    output      = '/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2015/'

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = requestName
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = '../customL1NtupleFromRaw.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['L1Tree.root']
config.JobType.pyCfgParams = pyCfg

config.section_('Data')
config.Data.inputDataset = dataset
config.Data.inputDBS = 'global'
if not(RunOnMC) : config.Data.lumiMask = 'run_mask.txt'
config.Data.splitting = splitting
config.Data.unitsPerJob = 25
config.Data.outLFNDirBase = output

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

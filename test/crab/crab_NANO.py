from WMCore.Configuration import Configuration

config = Configuration()

requestName = '254790_L1Accept'
pyCfg       = ['runOnMC=False', 'globalTag=74X_dataRun2_Prompt_v1']
dataset     = '/L1Accept/Run2015C-v1/RAW'
splitting   = 'LumiBased'
output      = '/store/group/dpg_trigger/comm_trigger/L1Trigger/Data/Collisions/'

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = requestName
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = '../customL1NtupleFromNanoDST.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['L1Tree.root']
config.JobType.pyCfgParams = pyCfg

config.section_('Data')
config.Data.inputDataset = dataset
config.Data.inputDBS = 'global'
config.Data.lumiMask = 'run_mask.txt'
config.Data.splitting = splitting
config.Data.unitsPerJob = 25
config.Data.outLFNDirBase = output

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

from WMCore.Configuration import Configuration

config = Configuration()

RunOnMC = True

requestName = 'RECO'
pyCfg       = ['runOnMC=False', 'globalTag=GR_P_V56']
dataset     = '/SingleMuon/Run2015B-PromptReco-v1/RECO'
splitting   = 'LumiBased'
output      = '/store/group/dpg_trigger/comm_trigger/L1Trigger/Data/Collisions/'

if RunOnMC :
    requestName = '50nsStart'
    pyCfg       = ['runOnMC=True', 'globalTag=MCRUN2_74_V8']
    dataset     = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/AODSIM'
    splitting   = 'FileBased'
    output      = '/store/group/dpg_trigger/comm_trigger/L1Trigger/efficiency/'

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = requestName
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = '../customL1NtupleFromFevt.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['L1Tree.root']
config.JobType.pyCfgParams = pyCfg

config.section_('Data')
config.Data.inputDataset = dataset
config.Data.inputDBS = 'global'
#if not(RunOnMC) : config.Data.lumiMask = 'run_mask.txt'
config.Data.splitting = splitting
config.Data.unitsPerJob = 25
config.Data.outLFNDirBase = output

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

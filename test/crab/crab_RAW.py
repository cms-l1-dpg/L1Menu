from WMCore.Configuration import Configuration

config = Configuration()

RunOnMC = False

requestName = '254790_ZeroBias_Reemul_2'
pyCfg       = ['runOnMC=False', 'globalTag=GR_P_V56','reEmulation=True', 'patchNtuple=True', 'reEmulCalos=True', 'runOnPostLS1=True', 'useStage1Layer2=True', 'whichPU=20']
dataset     = '/ZeroBias/Run2015C-v1/RAW'
splitting   = 'LumiBased'
output      = '/store/group/dpg_trigger/comm_trigger/L1Trigger/Data/Collisions/'

if RunOnMC :
    requestName = 'MC_50ns_Flat_for25'
    pyCfg       = ['runOnMC=True', 'globalTag=MCRUN2_74_V9A', 'reEmulation=True', 'patchNtuple=True', 'reEmulCalos=True', 'reEmulMuons=True', 'runOnPostLS1=True', 'useStage1Layer2=True', 'whichPU=20']
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
config.JobType.inputFiles = ['../../data/Jet_Stage1_2015_v2.txt']

config.section_('Data')
config.Data.inputDataset = dataset
config.Data.inputDBS = 'global'
if not(RunOnMC) : config.Data.lumiMask = 'run_mask.txt'
config.Data.splitting = splitting
config.Data.unitsPerJob = 25
config.Data.outLFNDirBase = output

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

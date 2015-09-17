from WMCore.Configuration import Configuration

config = Configuration()

RunOnMC = False

requestName = 'RECO_Run2015C'
pyCfg       = ['runOnMC=False', 'globalTag=74X_dataRun2_Prompt_v1']
dataset     = '/SingleMuon/Run2015C-PromptReco-v1/RECO'
splitting   = 'LumiBased'
output      = '/store/group/dpg_trigger/comm_trigger/L1Trigger/Data/Collisions/'

if RunOnMC :
    requestName = '25nsStart'
    pyCfg       = ['runOnMC=True', 'globalTag=MCRUN2_74_V9']
    dataset     = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-AsymptFlat10to50bx25Raw_MCRUN2_74_V9-v1/AODSIM'
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
if not(RunOnMC) : config.Data.lumiMask = 'Cert_246908-254879_13TeV_PromptReco_Collisions15_JSON_MuonPhys.txt'
config.Data.splitting = splitting
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = output
config.Data.useParent = True

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

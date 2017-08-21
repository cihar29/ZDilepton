from WMCore.Configuration import Configuration
config = Configuration()

isMC = True

if isMC:
  inputFiles = []
  outputFiles = ["analysis_MC.root"]
  splitting = 'FileBased'
  lumiMask = ''
  unitsPerJob = 3

else:
  inputFiles = ["pileup_12_6_16.txt"]
  outputFiles = ["analysis_Data.root"]
  splitting = 'LumiBased'
  lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
            #'/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
            #'/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
  unitsPerJob = 100

config.section_("General")
config.General.requestName = 'ST_tAug2017'
config.General.workArea = 'crab_projects'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_zdilepton.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = inputFiles
config.JobType.outputFiles = outputFiles

config.section_("Data")
config.Data.inputDataset = '/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
  #'/SingleMuon/Run2016B-03Feb2017_ver2-v2/MINIAOD'
  #'/SingleMuon/Run2016G-03Feb2017-v1/MINIAOD'
  #'/SingleMuon/Run2016H-03Feb2017_ver2-v1/MINIAOD'
  #'/DoubleEG/Run2016B-03Feb2017_ver2-v2/MINIAOD'
  #'/DoubleEG/Run2016G-03Feb2017-v1/MINIAOD'
  #'/DoubleEG/Run2016H-03Feb2017_ver2-v1/MINIAOD'
  #'/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
  #'/ST_s-channel_4f_InclusiveDecays_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
  #'/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
  #'/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
  #'/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'

config.Data.splitting = splitting
config.Data.lumiMask = lumiMask

config.Data.unitsPerJob = unitsPerJob
config.Data.outLFNDirBase = '/store/user/charring/AnalysisAug2017'
config.Data.publication = False
#config.Data.ignoreLocality = True
#config.Data.publishDataName = 'analysis_tree'

config.section_("Site")
#config.Site.blacklist = ['T1_US_FNAL']
#config.Site.whitelist = ['T2_CH_CERN']
config.Site.storageSite = "T3_US_FNALLPC"

# source /cvmfs/cms.cern.ch/crab3/crab.csh

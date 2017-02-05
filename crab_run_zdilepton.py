from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'MuonF'
config.General.workArea = 'crab_projects'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_zdilepton.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ["pileup_12_6_16.txt"]
config.JobType.outputFiles = ["analysis_Data.root"]

config.section_("Data")
config.Data.inputDataset = '/SingleMuon/Run2016F-23Sep2016-v1/MINIAOD'
  #'/SingleMuon/Run2016E-23Sep2016-v1/MINIAOD'
  #'/SingleMuon/Run2016B-23Sep2016-v1/MINIAOD'
  #'/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/MINIAODSIM'
  #'/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'
  #'/SingleMuon/Run2016B-23Sep2016-v1/MINIAOD'

config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'

config.Data.unitsPerJob = 25
config.Data.outLFNDirBase = '/store/user/charring/analysis'
config.Data.publication = False
#config.Data.ignoreLocality = True
#config.Data.publishDataName = 'analysis_tree'

config.section_("Site")
#config.Site.blacklist = ['T1_US_FNAL']
#config.Site.whitelist = ['T2_CH_CERN']
config.Site.storageSite = "T3_US_FNALLPC"

# source /cvmfs/cms.cern.ch/crab3/crab.csh

# PYTHON configuration file for class: ZDilepton
# Author: C. Harrington
# Date:  5 - October - 2016

import FWCore.ParameterSet.Config as cms

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

process = cms.Process("Ana")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

readFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource", fileNames = readFiles)
readFiles.extend( [

  #'/store/data/Run2016H/SingleMuon/MINIAOD/PromptReco-v2/000/281/280/00000/008681A1-7382-E611-B223-02163E01190F.root'

  #'/store/data/Run2016E/SingleMuon/MINIAOD/23Sep2016-v1/50000/00CFC689-8D8D-E611-9F90-0CC47A13D16E.root'

  #'/store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/1110000/A2C0F697-B19C-E611-A4D8-F04DA275BFF2.root'

  #'/store/data/Run2016B/SingleMuon/MINIAOD/23Sep2016-v1/70000/029A7F31-B187-E611-B7D5-0CC47A13CECE.root'

  #'/store/data/Run2016G/SingleElectron/MINIAOD/23Sep2016-v1/100000/004A7893-A990-E611-B29F-002590E7DE36.root'

  #'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root'

  #'/store/mc/RunIISummer16MiniAODv2/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/001C9040-67B9-E611-AD6E-0CC47A7EEE32.root'

  '/store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_InclusiveDecays_13TeV-amcatnlo-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/0057C751-D8F6-E611-9151-68B59972C578.root'

] );

gt = '80X_dataRun2_2016SeptRepro_v5'
#'80X_dataRun2_Prompt_v15'
#'80X_dataRun2_2016SeptRepro_v5'

process.load( "Configuration.Geometry.GeometryIdeal_cff" )
process.load( "Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff" )
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, gt)

isMC = cms.bool(True)

if isMC:
  OutputName = "_MC"  
  metLabel = "SIM"

else:
  OutputName = "_Data"
  metLabel = "RECO"

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)

my_eid_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff']
for idmod in my_eid_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.analysis = cms.EDAnalyzer("ZDilepton",
    fileName = cms.string("analysis" + OutputName + ".root"),
    btag = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    isMC = isMC,
    patTrgLabel = cms.InputTag("TriggerResults", "", "metLabel"),
    BadPFMuonFilter = cms.InputTag("BadPFMuonFilter",""),
    BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter",""),
    rhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
    pvTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
    genParticleTag = cms.InputTag("prunedGenParticles"),
    genJetTag = cms.InputTag("slimmedGenJets"),
    muonTag = cms.InputTag("slimmedMuons"),
    electronTag = cms.InputTag("slimmedElectrons"),
    effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt"),
    eleVetoIdMapToken = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
    eleLooseIdMapToken = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
    eleMediumIdMapToken = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
    eleTightIdMapToken = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
    convLabel = cms.InputTag("reducedEgamma:reducedConversions"),
    jetTag = cms.InputTag("slimmedJets"),
    metTag = cms.InputTag("slimmedMETs"),
    minLepPt = cms.double(45.),
    minSubLepPt = cms.double(25.),
    minDiLepMass = cms.double(20.),
    minLeadJetPt = cms.double(80.),
    triggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
    prescalesTag = cms.InputTag("patTrigger"),
    genEventTag = cms.InputTag("generator"),
    muTag = cms.InputTag("slimmedAddPileupInfo")
    #convTag = cms.InputTag("reducedConversions"),
    #bsTag = cms.InputTag("offlineBeamSpot")
)

process.myseq = cms.Sequence( process.BadPFMuonFilter *
                              process.BadChargedCandidateFilter *
                              process.egmGsfElectronIDSequence *
                              process.analysis )

process.p = cms.Path( process.myseq )

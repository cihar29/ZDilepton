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

  #RSGluon M-3000
  #'/store/mc/RunIISummer16MiniAODv2/RSGluonToTT_M-3000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/220B9096-25B7-E611-9E16-141877344D39.root',
  #'/store/mc/RunIISummer16MiniAODv2/RSGluonToTT_M-3000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/2CE514B5-2AB7-E611-BAAC-008CFA165F18.root',
  #'/store/mc/RunIISummer16MiniAODv2/RSGluonToTT_M-3000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/AABE4A83-29B7-E611-8747-0025907B4EC8.root',
  #'/store/mc/RunIISummer16MiniAODv2/RSGluonToTT_M-3000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/F697AF2E-26B7-E611-BA71-0090FAA57400.root'

  #Z' M-3000 W-300
  #'/store/mc/RunIISummer16MiniAODv2/ZprimeToTT_M-3000_W-300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/68589978-7CBE-E611-9AE2-0CC47A1E0DC2.root',
  #'/store/mc/RunIISummer16MiniAODv2/ZprimeToTT_M-3000_W-300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/68728876-7ABE-E611-97A7-0CC47A1E046E.root',
  #'/store/mc/RunIISummer16MiniAODv2/ZprimeToTT_M-3000_W-300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/A0510FBD-63BE-E611-A4C4-484D7E8DF0C6.root',
  #'/store/mc/RunIISummer16MiniAODv2/ZprimeToTT_M-3000_W-300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/DA381FBB-75BE-E611-BC99-0CC47A1E0DC8.root',
  #'/store/mc/RunIISummer16MiniAODv2/ZprimeToTT_M-3000_W-300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/E6F9045C-77BE-E611-B2B5-0CC47A1E046E.root'

  '/store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver2-v2/100000/000C6E52-8BEC-E611-B3FF-0025905C42FE.root'

  #'/store/data/Run2016G/SingleMuon/MINIAOD/03Feb2017-v1/100000/00E6DF50-70EA-E611-ACC4-0CC47A1E089C.root'

  #'/store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/00933E2A-A0D5-E611-B2CD-00266CF89130.root'

  #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/02A210D6-F5C3-E611-B570-008CFA197BD4.root'

  #'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root'

  #'/store/mc/RunIISummer16MiniAODv2/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/001C9040-67B9-E611-AD6E-0CC47A7EEE32.root'

  #'/store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_InclusiveDecays_13TeV-amcatnlo-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/0057C751-D8F6-E611-9151-68B59972C578.root'

] );

isMC = cms.bool(False)
gts = {'BCDEFG':'80X_dataRun2_2016SeptRepro_v7', 'H':'80X_dataRun2_Prompt_v16'}

if isMC:
  OutputName = "_MC"  
  metLabel = "SIM"
  gt = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'

else:
  OutputName = "_Data"
  metLabel = "RECO"
  gt = gts['BCDEFG']

process.load( "Configuration.Geometry.GeometryIdeal_cff" )
process.load( "Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff" )
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, gt)

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
    extLHETag = cms.InputTag("externalLHEProducer"),
    muTag = cms.InputTag("slimmedAddPileupInfo")
    #convTag = cms.InputTag("reducedConversions"),
    #bsTag = cms.InputTag("offlineBeamSpot")
)

process.myseq = cms.Sequence( process.BadPFMuonFilter *
                              process.BadChargedCandidateFilter *
                              process.egmGsfElectronIDSequence *
                              process.analysis )

process.p = cms.Path( process.myseq )

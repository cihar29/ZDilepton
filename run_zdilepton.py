# PYTHON configuration file for class: ZDilepton
# Author: C. Harrington
# Date:  5 - October - 2016

import FWCore.ParameterSet.Config as cms

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

process = cms.Process("Ana")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

readFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource", fileNames = readFiles)
readFiles.extend( [

  #'/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/150/00000/34A57FB8-D819-E611-B0A4-02163E0144EE.root'

  #'/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/1CCC1100-0E1A-E611-98C7-02163E014332.root'

  #'file:/uscms_data/d3/broozbah/Analysis_Zprime/CMSSW_8_0_19/src/Analysis_Zprime/ZDilepton/singleElectron.root'

  '/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00000/004A0552-3929-E611-BD44-0025905A48F0.root'

] );

gt = "80X_dataRun2_Prompt_v12"

process.load( "Configuration.Geometry.GeometryIdeal_cff" )
process.load( "Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff" )
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, gt)

isMC = cms.bool(True)

if isMC:
  OutputName = "_MC"  

else:
  OutputName = "_Data"

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)

my_eid_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']
                  #'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff'
for idmod in my_eid_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.analysis = cms.EDAnalyzer("ZDilepton",
    fileName = cms.string("analysis" + OutputName + ".root"),
    btag = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    isMC = isMC,
    metFilters = cms.bool(True),
    patTrgLabel = cms.InputTag("TriggerResults", "", "RECO"),
    BadPFMuonFilter = cms.InputTag("BadPFMuonFilter",""),
    BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter",""),
    rhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
    pvTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
    genParticleTag = cms.InputTag("prunedGenParticles"),
    genJetTag = cms.InputTag("slimmedGenJets"),
    muonTag = cms.InputTag("slimmedMuons"),
    electronTag = cms.InputTag("slimmedElectrons"),
    effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt"),
    eleVetoIdMapToken = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    eleLooseIdMapToken = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    eleMediumIdMapToken = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    eleTightIdMapToken = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    convLabel = cms.InputTag("reducedEgamma:reducedConversions"),
    jetTag = cms.InputTag("slimmedJets"),
    metTag = cms.InputTag("slimmedMETs"),
    metPuppiTag = cms.InputTag("slimmedMETsPuppi"),
    minLepPt = cms.double(45.),
    minSubLepPt = cms.double(25.),
    triggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
    prescalesTag = cms.InputTag("patTrigger"),
    genEventTag = cms.InputTag("generator"),
    muTag = cms.InputTag("slimmedAddPileupInfo")
)

process.myseq = cms.Sequence( process.BadPFMuonFilter *
                              process.BadChargedCandidateFilter *
                              process.egmGsfElectronIDSequence *
                              process.analysis )

process.p = cms.Path( process.myseq )

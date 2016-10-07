# PYTHON configuration file for class: ZDilepton
# Author: C. Harrington
# Date:  5 - October - 2016

import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

readFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource", fileNames = readFiles)
readFiles.extend( [

  '/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/150/00000/34A57FB8-D819-E611-B0A4-02163E0144EE.root'

] );

gt = "80X_dataRun2_Prompt_v12"

process.load( "Configuration.Geometry.GeometryIdeal_cff" )
process.load( "Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff" )
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, gt)

#Beam Halo
process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')

#HCAL HBHE
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
  inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResultRun2Tight'),
  reverseDecision = cms.bool(False)
)

#Bad EE Supercrystal filter
#process.load('RecoMET.METFilters.eeBadScFilter_cfi')

process.analysis = cms.EDAnalyzer("ZDilepton",
    RootFileName = cms.string("analysis.root"),
    rhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
    pvTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muonTag = cms.InputTag("slimmedMuons"),
    electronTag = cms.InputTag("slimmedElectrons")
)

process.myseq = cms.Sequence( process.analysis )

process.p = cms.Path( process.CSCTightHaloFilter *
                      process.HBHENoiseFilterResultProducer *
                      process.ApplyBaselineHBHENoiseFilter *
                      # process.eeBadScFilter *
                      process.myseq )

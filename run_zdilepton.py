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

era = "Spring16_25nsV6_DATA"
tag = "80X_dataRun2_Prompt_v12"

process.load( "Configuration.Geometry.GeometryIdeal_cff" )
process.load( "Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff" )
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, tag)

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
process.load('RecoMET.METFilters.eeBadScFilter_cfi')

#Jet Corrections
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
    connect = cms.string('sqlite_file:'+era+'.db'),
    toGet =  cms.VPSet(
        cms.PSet(
            record = cms.string("JetCorrectionsRecord"),
            tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
            label= cms.untracked.string("AK4PFchs")
        )
    )
)
process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")

process.load('JetMETCorrections.Configuration.JetCorrectors_cff')

process.ak4PFJetsCHSl1l2l3 = cms.EDProducer('CorrectedPFJetProducer',
    src         = cms.InputTag("ak4PFJetsCHS"),
    correctors  = cms.VInputTag("ak4PFCHSL1FastL2L3ResidualCorrector")
)

process.analysis = cms.EDAnalyzer("ZDilepton",
    RootFileName = cms.string("analysis.root"),
    muTag = cms.InputTag("addPileupInfo"),
    rhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
    pvTag = cms.InputTag("offlinePrimaryVertices"),
    muonTag = cms.InputTag("slimmedMuons")
)

process.myseq = cms.Sequence( process.ak4PFJetsCHSl1l2l3 * process.analysis )

process.p = cms.Path( process.CSCTightHaloFilter *
                      process.HBHENoiseFilterResultProducer *
                      process.ApplyBaselineHBHENoiseFilter *
                      process.eeBadScFilter *
                      process.myseq )

// -*- C++ -*-
//
// Package:    analysis/ZDilepton
// Class:      ZDilepton
// 
/**\class ZDilepton ZDilepton.cc ZDilepton/plugins/ZDilepton.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  charles harrington and bahareh roozbahani
//         Created:  Wed, 05 Oct 2016 17:09:43 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
//#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "parsePileUpJSON2.h"
#include <vector>
#include <algorithm>
#include <utility>

//root files
#include <TFile.h>
#include <TTree.h>
#include <TVector.h>
#include <TLorentzVector.h>

using namespace std;
using namespace edm;

const int MAXLEP = 20;
const int MAXGEN = 20;
const int MAXJET = 50;
const int nFilters = 6;
const int nTriggers = 8;
const int METUNCERT = 4;

bool sortLepPt(const pair<reco::CandidatePtr, char>& lep1, const pair<reco::CandidatePtr, char>& lep2){ return lep1.first->pt() > lep2.first->pt(); }

class ZDilepton : public edm::EDAnalyzer {
  public:
    explicit ZDilepton(const edm::ParameterSet&);

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    TFile* root_file;
    TTree* tree;

    string filters[nFilters] = {
      "Flag_HBHENoiseFilter",
      "Flag_HBHENoiseIsoFilter", 
      "Flag_EcalDeadCellTriggerPrimitiveFilter",
      "Flag_goodVertices",
      "Flag_eeBadScFilter",
      "Flag_globalTightHalo2016Filter"
    };
    vector<int> totalEvts, filter_failed, dilep_cut, leppt_cut, jetpteta_cut, met_cut, dilepmass_cut;

    string triggers[nTriggers] = {
      "HLT_Mu45_eta2p1_v",
      "HLT_Mu50_v",
      "HLT_TkMu50_v",
      "HLT_Mu30_TkMu11_v",
      "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v",
      "HLT_Ele105_CaloIdVT_GsfTrkIdT_v",
      "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
      "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v"
    };
    vector<int> trig_prescale; vector<bool> trig_passed; vector<string> trig_name;

    ULong64_t event;
    int run, lumi, bx;

    float rho, mu, genweight;
    int nPV;

    int nGen;
    int gen_status[MAXGEN], gen_PID[MAXGEN], gen_mother0[MAXGEN], gen_mother1[MAXGEN], gen_index[MAXGEN]; //gen_nMothers[MAXGEN], gen_nDaughters[MAXGEN];
    float gen_pt[MAXGEN], gen_mass[MAXGEN], gen_eta[MAXGEN], gen_phi[MAXGEN];
    int nGenJet;
    float genJet_pt[MAXJET], genJet_eta[MAXJET], genJet_phi[MAXJET], genJet_mass[MAXJET], genJet_area[MAXJET], genJet_nDaught[MAXJET];

    char lep0flavor, lep1flavor;

    int nMuon;
    bool muon_isGlob[MAXLEP], muon_IsMediumID[MAXLEP], muon_IsTightID[MAXLEP];
    int muon_charge[MAXLEP];
    float muon_pt[MAXLEP], muon_eta[MAXLEP], muon_phi[MAXLEP], muon_D0[MAXLEP], muon_Dz[MAXLEP];
    float muon_chi2[MAXLEP], muon_tspm[MAXLEP], muon_kinkf[MAXLEP], muon_segcom[MAXLEP], muon_ftrackhits[MAXLEP];
    
    int nEle;
    int ele_charge[MAXLEP], ele_missinghits[MAXLEP];
    bool ele_LooseID[MAXLEP], ele_MediumID[MAXLEP], ele_TightID[MAXLEP]; //ele_passConv[MAXLEP];
    float ele_pt[MAXLEP], ele_eta[MAXLEP], ele_phi[MAXLEP], ele_D0[MAXLEP], ele_Dz[MAXLEP], ele_etaSupClust[MAXLEP];
    float ele_dPhiIn[MAXLEP], ele_sigmaIetaIeta[MAXLEP], ele_dEtaSeed[MAXLEP], ele_HE[MAXLEP], ele_rcpiwec[MAXLEP], ele_overEoverP[MAXLEP];

    int nJet;
    float jet_pt[MAXJET], jet_eta[MAXJET], jet_phi[MAXJET], jet_mass[MAXJET], jet_area[MAXJET], jet_jec[MAXJET], jet_btag[MAXJET];
    float jet_nhf[MAXJET], jet_nef[MAXJET], jet_chf[MAXJET], jet_muf[MAXJET];
    float jet_elef[MAXJET], jet_numneutral[MAXJET], jet_chmult[MAXJET];
    char jet_clean[MAXJET];

    float genmet_pt, genmet_px, genmet_py, genmet_sumet, genmet_phi;

    float met_pt, met_px, met_py, met_sumet, met_phi;
    int nMETUncert;
    float met_shiftedpx[METUNCERT], met_shiftedpy[METUNCERT];

    TString fileName_;
    string btag_;
    bool isMC_;
    double minLepPt_, minSubLepPt_, minDiLepMass_, minLeadJetPt_;

    EffectiveAreas ele_areas_;
    edm::EDGetTokenT<edm::TriggerResults> patTrgLabel_;
    edm::EDGetTokenT<double> rhoTag_;
    edm::EDGetTokenT< edm::View<reco::Vertex> > pvTag_;
    edm::EDGetTokenT< edm::View<reco::GenParticle> > genParticleTag_;
    edm::EDGetTokenT< edm::View<reco::GenJet> > genJetTag_;
    edm::EDGetTokenT< edm::View<pat::Muon> > muonTag_;
    edm::EDGetTokenT< edm::View<pat::Electron> > electronTag_;
    edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleVetoIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleLooseIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleMediumIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleTightIdMapToken_;
    edm::EDGetTokenT< edm::View<pat::Jet> > jetTag_;
    edm::EDGetTokenT< edm::View<pat::MET> > metTag_;
    edm::EDGetTokenT< edm::View<pat::MET> > metPuppiTag_;
    edm::EDGetTokenT<bool> BadChCandFilterToken_;
    edm::EDGetTokenT<bool> BadPFMuonFilterToken_;
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsTag_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> prescalesTag_;
    edm::EDGetTokenT<GenEventInfoProduct> genEventTag_;
    edm::EDGetTokenT< edm::View<PileupSummaryInfo> > muTag_;
    //edm::EDGetTokenT< vector<reco::Conversion> > convTag_;
    //edm::EDGetTokenT<reco::BeamSpot> bsTag_;
};

ZDilepton::ZDilepton(const edm::ParameterSet& iConfig):
  ele_areas_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath() )
{
  fileName_ = iConfig.getParameter<string>("fileName");
  btag_ = iConfig.getParameter<string>("btag");
  isMC_ = iConfig.getParameter<bool>("isMC");
  patTrgLabel_ = consumes<edm::TriggerResults>( iConfig.getParameter<edm::InputTag>("patTrgLabel") );
  rhoTag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoTag") );
  pvTag_ = consumes< edm::View<reco::Vertex> >( iConfig.getParameter<edm::InputTag>("pvTag") );
  genParticleTag_ =  consumes< edm::View<reco::GenParticle> >( iConfig.getParameter<edm::InputTag>("genParticleTag") );
  genJetTag_ =  consumes< edm::View<reco::GenJet> >( iConfig.getParameter<edm::InputTag>("genJetTag") );
  muonTag_ = consumes< edm::View<pat::Muon> >( iConfig.getParameter<edm::InputTag>("muonTag") );
  electronTag_ = consumes< edm::View<pat::Electron> >( iConfig.getParameter<edm::InputTag>("electronTag") );
  eleVetoIdMapToken_ = consumes<edm::ValueMap<vid::CutFlowResult> >( iConfig.getParameter<edm::InputTag>("eleVetoIdMapToken") );
  eleLooseIdMapToken_ = consumes<edm::ValueMap<vid::CutFlowResult> >( iConfig.getParameter<edm::InputTag>("eleLooseIdMapToken") );
  eleMediumIdMapToken_ = consumes<edm::ValueMap<vid::CutFlowResult> >( iConfig.getParameter<edm::InputTag>("eleMediumIdMapToken") );
  eleTightIdMapToken_ = consumes<edm::ValueMap<vid::CutFlowResult> >( iConfig.getParameter<edm::InputTag>("eleTightIdMapToken") );
  jetTag_ = consumes< edm::View<pat::Jet> >( iConfig.getParameter<edm::InputTag>("jetTag") );
  metTag_ = consumes< edm::View<pat::MET> >( iConfig.getParameter<edm::InputTag>("metTag") );
  BadChCandFilterToken_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter"));
  BadPFMuonFilterToken_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadPFMuonFilter"));
  minLepPt_ = iConfig.getParameter<double>("minLepPt");
  minSubLepPt_ = iConfig.getParameter<double>("minSubLepPt");
  minDiLepMass_ = iConfig.getParameter<double>("minDiLepMass");
  minLeadJetPt_ = iConfig.getParameter<double>("minLeadJetPt");
  triggerResultsTag_ = consumes<edm::TriggerResults>( iConfig.getParameter<edm::InputTag>("triggerResultsTag") );
  prescalesTag_ = consumes<pat::PackedTriggerPrescales>( iConfig.getParameter<edm::InputTag>("prescalesTag") );
  genEventTag_ = consumes<GenEventInfoProduct>( iConfig.getParameter<edm::InputTag>("genEventTag") );
  muTag_ = consumes< edm::View<PileupSummaryInfo> >( iConfig.getParameter<edm::InputTag>("muTag") );
  //convTag_ = consumes< vector<reco::Conversion> >( iConfig.getParameter<edm::InputTag>("convTag") );
  //bsTag_ = consumes<reco::BeamSpot>( iConfig.getParameter<edm::InputTag>("bsTag") );
}

// ------------ method called once each job just before starting event loop  ------------
void  ZDilepton::beginJob() {

  root_file = new TFile(fileName_, "RECREATE");
  tree = new TTree("T", "Analysis Tree");

  filter_failed.assign(nFilters+2, 0);
  totalEvts.assign(1, 0); dilep_cut.assign(1, 0); leppt_cut.assign(1, 0); jetpteta_cut.assign(1, 0); met_cut.assign(1, 0); dilepmass_cut.assign(1, 0);

  tree->Branch("trig_prescale", "std::vector<int>", &trig_prescale);
  tree->Branch("trig_passed", "std::vector<bool>", &trig_passed);
  tree->Branch("trig_name", "std::vector<string>", &trig_name);

  tree->Branch("run", &run, "run/I");
  tree->Branch("lumi", &lumi, "lumi/I");
  tree->Branch("bx", &bx, "bx/I");
  tree->Branch("event", &event, "event/l");

  if(isMC_){
    tree->Branch("nGen", &nGen, "nGen/I");
    tree->Branch("gen_status",  gen_status, "gen_status[nGen]/I");
    tree->Branch("gen_PID",  gen_PID, "gen_PID[nGen]/I");
    tree->Branch("gen_pt", gen_pt, "gen_pt[nGen]/F");
    tree->Branch("gen_mass", gen_mass, "gen_mass[nGen]/F");
    tree->Branch("gen_eta", gen_eta, "gen_eta[nGen]/F");
    tree->Branch("gen_phi", gen_phi, "gen_phi[nGen]/F");
    tree->Branch("gen_index",  gen_index, "gen_index[nGen]/I");
    //tree->Branch("gen_nMothers", gen_nMothers, "gen_nMothers[nGen]/I");
    //tree->Branch("gen_nDaughters", gen_nDaughters, "gen_nDaughters[nGen]/I");
    tree->Branch("gen_mother0", gen_mother0, "gen_mother0[nGen]/I");
    tree->Branch("gen_mother1", gen_mother1, "gen_mother1[nGen]/I");

    tree->Branch("nGenJet", &nGenJet, "nGenJet/I");
    tree->Branch("genJet_pt", genJet_pt, "genJet_pt[nGenJet]/F");
    tree->Branch("genJet_eta", genJet_eta, "genJet_eta[nGenJet]/F");
    tree->Branch("genJet_phi", genJet_phi, "genJet_phi[nGenJet]/F");
    tree->Branch("genJet_mass", genJet_mass, "genJet_mass[nGenJet]/F");
    tree->Branch("genJet_area",  genJet_area, "genJet_area[nGenJet]/F");
    tree->Branch("genJet_nDaught",  genJet_nDaught, "genJet_nDaught[nGenJet]/F");

    tree->Branch("genmet_pt", &genmet_pt, "genmet_pt/F");
    tree->Branch("genmet_px", &genmet_px, "genmet_px/F");
    tree->Branch("genmet_py", &genmet_py, "genmet_py/F");
    tree->Branch("genmet_sumet", &genmet_sumet, "genmet_sumet/F");
    tree->Branch("genmet_phi", &genmet_phi, "genmet_phi/F");

    tree->Branch("genweight", &genweight, "genweight/F");
  }
  else{
    parsePileUpJSON2();
  }

  tree->Branch("rho", &rho, "rho/F");
  tree->Branch("mu", &mu, "mu/F");
  tree->Branch("nPV", &nPV, "nPV/I");

  tree->Branch("lep0flavor", &lep0flavor, "lep0flavor/B");
  tree->Branch("lep1flavor", &lep1flavor, "lep1flavor/B");

  tree->Branch("nMuon", &nMuon, "nMuon/I");
  tree->Branch("muon_charge", muon_charge, "muon_charge[nMuon]/I");
  tree->Branch("muon_isGlob", muon_isGlob, "muon_isGlob[nMuon]/O");
  tree->Branch("muon_pt", muon_pt, "muon_pt[nMuon]/F");
  tree->Branch("muon_eta", muon_eta, "muon_eta[nMuon]/F");
  tree->Branch("muon_phi", muon_phi, "muon_phi[nMuon]/F");
  tree->Branch("muon_D0", muon_D0, "muon_D0[nMuon]/F");
  tree->Branch("muon_Dz", muon_Dz, "muon_Dz[nMuon]/F");
  tree->Branch("muon_chi2", muon_chi2, "muon_chi2[nMuon]/F");
  tree->Branch("muon_tspm", muon_tspm, "muon_tspm[nMuon]/F");
  tree->Branch("muon_kinkf", muon_kinkf, "muon_kinkf[nMuon]/F");
  tree->Branch("muon_ftrackhits", muon_ftrackhits, "muon_ftrackhits[nMuon]/F");
  tree->Branch("muon_segcom", muon_segcom, "muon_segcom[nMuon]/F");
  tree->Branch("muon_IsMediumID", muon_IsMediumID, "muon_IsMediumID[nMuon]/O");
  tree->Branch("muon_IsTightID", muon_IsTightID, "muon_IsTightID[nMuon]/O");

  tree->Branch("nEle", &nEle, "nEle/I");
  tree->Branch("ele_charge", ele_charge, "ele_charge[nEle]/I");
  tree->Branch("ele_pt", ele_pt, "ele_pt[nEle]/F");
  tree->Branch("ele_eta", ele_eta, "ele_eta[nEle]/F");
  tree->Branch("ele_phi", ele_phi, "ele_phi[nEle]/F");
  tree->Branch("ele_LooseID", ele_LooseID , "ele_LooseID[nEle]/O");
  tree->Branch("ele_MediumID", ele_MediumID , "ele_MediumID[nEle]/O");
  tree->Branch("ele_TightID", ele_TightID , "ele_TightID[nEle]/O");
  //tree->Branch("ele_passConv", ele_passConv , "ele_passConv[nEle]/O");
  tree->Branch("ele_D0", ele_D0, "ele_D0[nEle]/F");
  tree->Branch("ele_Dz", ele_Dz, "ele_Dz[nEle]/F");
  tree->Branch("ele_sigmaIetaIeta", ele_sigmaIetaIeta, "ele_sigmaIetaIeta[nEle]/F");
  tree->Branch("ele_dEtaSeed", ele_dEtaSeed, "ele_dEtaSeed[nEle]/F");
  tree->Branch("ele_dPhiIn", ele_dPhiIn, "ele_dPhiIn[nEle]/F");
  tree->Branch("ele_etaSupClust", ele_etaSupClust, "ele_etaSupClust[nEle]/F");
  tree->Branch("ele_overEoverP", ele_overEoverP, "ele_overEoverP[nEle]/F");
  tree->Branch("ele_HE", ele_HE, "ele_HE[nEle]/F");
  tree->Branch("ele_rcpiwec", ele_rcpiwec, "ele_rcpiwec[nEle]/F");
  tree->Branch("ele_missinghits", ele_missinghits, "ele_missinghits[nEle]/I");

  tree->Branch("nJet", &nJet, "nJet/I");
  tree->Branch("jet_pt", jet_pt, "jet_pt[nJet]/F");
  tree->Branch("jet_eta", jet_eta, "jet_eta[nJet]/F");
  tree->Branch("jet_phi", jet_phi, "jet_phi[nJet]/F");
  tree->Branch("jet_mass", jet_mass, "jet_mass[nJet]/F");
  tree->Branch("jet_area", jet_area, "jet_area[nJet]/F");
  tree->Branch("jet_jec", jet_jec, "jet_jec[nJet]/F");
  tree->Branch("jet_btag", jet_btag, "jet_btag[nJet]/F");
  tree->Branch("jet_clean", jet_clean, "jet_clean[nJet]/B");

  tree->Branch("jet_nhf", jet_nhf, "jet_nhf[nJet]/F");
  tree->Branch("jet_nef", jet_nef, "jet_nef[nJet]/F");
  tree->Branch("jet_chf", jet_chf, "jet_chf[nJet]/F");
  tree->Branch("jet_muf", jet_muf, "jet_muf[nJet]/F");
  tree->Branch("jet_elef", jet_elef, "jet_elef[nJet]/F");
  tree->Branch("jet_numneutral", jet_numneutral, "jet_numneutral[nJet]/F");
  tree->Branch("jet_chmult", jet_chmult, "jet_chmult[nJet]/F");

  tree->Branch("met_pt", &met_pt, "met_pt/F");
  tree->Branch("met_px", &met_px, "met_px/F");
  tree->Branch("met_py", &met_py, "met_py/F");
  tree->Branch("met_sumet", &met_sumet, "met_sumet/F");
  tree->Branch("met_phi", &met_phi, "met_phi/F");

  tree->Branch("nMETUncert", &nMETUncert, "nMETUncert/I");
  tree->Branch("met_shiftedpx", met_shiftedpx, "met_shiftedpx[nMETUncert]/F");
  tree->Branch("met_shiftedpy", met_shiftedpy, "met_shiftedpy[nMETUncert]/F");
}

// ------------ method called for each event  ------------
void ZDilepton::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  totalEvts[0]++;

  //------------ Lepton Pt Filter ------------//

  edm::Handle< edm::View<pat::Muon> > muons;
  iEvent.getByToken(muonTag_, muons);

  edm::Handle< edm::View<pat::Electron> > electrons;
  iEvent.getByToken(electronTag_, electrons);

  edm::Handle<edm::ValueMap<vid::CutFlowResult> > Veto_id_decisions;
  iEvent.getByToken(eleVetoIdMapToken_, Veto_id_decisions);

  //edm::Handle< vector<reco::Conversion> > conversions;
  //iEvent.getByToken(convTag_, conversions);
  //edm::Handle<reco::BeamSpot> beamspot;
  //iEvent.getByToken(bsTag_, beamspot);

  vector<pair<reco::CandidatePtr, char> > leps;

  for (int i=0, n=muons->size(); i<n; i++) {
    if (!muons->at(i).isLooseMuon()) continue;
    leps.push_back( make_pair(muons->ptrAt(i),'m') );
  }
  const int isoCut = 7;
  for (int i=0, n=electrons->size(); i<n; i++) {
    const Ptr<pat::Electron> elPtr(electrons, i);
    if (!(*Veto_id_decisions)[elPtr].getCutFlowResultMasking(isoCut).cutFlowPassed()) continue;
    leps.push_back( make_pair(electrons->ptrAt(i),'e') );
  }

  if (leps.size() < 2) return;
  dilep_cut[0]++;

  sort(leps.begin(), leps.end(), sortLepPt);

  if (leps[0].first->pt() < minLepPt_ || leps[1].first->pt() < minSubLepPt_) return;
  leppt_cut[0]++;

  bool dilepmass_flag = false;
  for (int i=0, n=leps.size(); i<n && !dilepmass_flag; i++) {

    char flavi = leps[i].second;
    reco::CandidatePtr lepi = leps[i].first;
    TLorentzVector vlepi;
    vlepi.SetPtEtaPhiM( lepi->pt(), lepi->eta(), lepi->phi(), lepi->mass() );

    for (int j=i+1; j<n && !dilepmass_flag; j++) {
      char flavj = leps[j].second;

      if (flavi==flavj) {
        reco::CandidatePtr lepj = leps[j].first;
        TLorentzVector vlepj;
        vlepj.SetPtEtaPhiM( lepj->pt(), lepj->eta(), lepj->phi(), lepj->mass() );

        if ( (vlepi+vlepj).M() > minDiLepMass_ ) dilepmass_flag = true;
      }
    }
  }
  if (!dilepmass_flag) return;
  dilepmass_cut[0]++;

  lep0flavor = leps[0].second;
  lep1flavor = leps[1].second;

  vector<reco::CandidatePtr> lep0Sources, lep1Sources;
  for (unsigned int i=0, n=leps[0].first->numberOfSourceCandidatePtrs(); i<n; i++) lep0Sources.push_back(leps[0].first->sourceCandidatePtr(i));
  for (unsigned int i=0, n=leps[1].first->numberOfSourceCandidatePtrs(); i<n; i++) lep1Sources.push_back(leps[1].first->sourceCandidatePtr(i));

  //------------ Jets ------------//

  edm::Handle< edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetTag_, jets);

  nJet = jets->size();
  int leadJetPt_flag = false;

  for (int i=0; i<nJet; i++){

    const pat::Jet& jet = jets->at(i).correctedJet(0);
    reco::Candidate::LorentzVector jet_p4 = jet.p4();

    jet_clean[i] = 'n';
    if ( reco::deltaR(jet, *leps[0].first) < 0.4 || reco::deltaR(jet, *leps[1].first) < 0.4 ){

      const vector<reco::CandidatePtr> & dvec = jet.daughterPtrVector();
      for (vector<reco::CandidatePtr>::const_iterator i_d = dvec.begin(); i_d != dvec.end(); ++i_d){

        if ( find(lep0Sources.begin(), lep0Sources.end(), *i_d ) != lep0Sources.end() ) {
          jet_p4 -= (*i_d)->p4();
          if (jet_clean[i] == 'n') jet_clean[i] = 'l';
          else if (jet_clean[i] == 's') jet_clean[i] = 'b';
        }
        else if ( find(lep1Sources.begin(), lep1Sources.end(), *i_d ) != lep1Sources.end() ) {
          jet_p4 -= (*i_d)->p4();
          if (jet_clean[i] == 'n') jet_clean[i] = 's';
          else if (jet_clean[i] == 'l') jet_clean[i] = 'b';
        }
      }
    }

    jet_pt[i] = jet_p4.Pt();
    jet_eta[i] = jet_p4.Eta();
    jet_phi[i] = jet_p4.Phi();
    jet_mass[i] = jet_p4.M();
    jet_area[i] = jet.jetArea();
    jet_jec[i] = jets->at(i).jecFactor(0);

    jet_nhf[i] = jet.neutralHadronEnergyFraction();
    jet_nef[i] = jet.neutralEmEnergyFraction();
    jet_chf[i] = jet.chargedHadronEnergyFraction();
    jet_muf[i] = jet.muonEnergyFraction();
    jet_elef[i] = jet.chargedEmEnergyFraction();
    jet_numneutral[i] =jet.neutralMultiplicity();
    jet_chmult[i] = jet.chargedMultiplicity();
    jet_btag[i] = jet.bDiscriminator(btag_);

    if ( jet_pt[i]>minLeadJetPt_ && fabs(jet_eta[i])<2.5 ) leadJetPt_flag = true;
  }
  if (!leadJetPt_flag) return;
  jetpteta_cut[0]++;

  //------------ MET Filters ------------//

  Handle<edm::TriggerResults> patFilterHandle;
  iEvent.getByToken(patTrgLabel_, patFilterHandle);

  if (patFilterHandle.isValid()){
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*patFilterHandle);

    //for (unsigned int i=0; i<trigNames.size(); i++) cout << trigNames.triggerName(i) << endl;

    for (int i=0; i<nFilters; i++){

      if( filters[i] == "Flag_eeBadScFilter" && isMC_ ) continue;

      const unsigned int trig = trigNames.triggerIndex(filters[i]);
      if (trig != trigNames.size()){
        if (!patFilterHandle->accept(trig)){
          filter_failed[i]++;
          return;
        }
      }
      else cout << filters[i] << " not found." << endl;
    }
  }

  edm::Handle<bool> ifilterbadChCand;
  iEvent.getByToken(BadChCandFilterToken_, ifilterbadChCand);
  bool filterbadChCandidate = *ifilterbadChCand;

  if (!filterbadChCandidate){
    filter_failed[nFilters]++;
    return;
  }

  edm::Handle<bool> ifilterbadPFMuon;
  iEvent.getByToken(BadPFMuonFilterToken_, ifilterbadPFMuon);
  bool filterbadPFMuon = *ifilterbadPFMuon;

  if (!filterbadPFMuon){
    filter_failed[nFilters+1]++;
    return;
  }
  met_cut[0]++;

  //------------ Event Info ------------//

  run = int(iEvent.id().run());
  lumi = int(iEvent.getLuminosityBlock().luminosityBlock());
  bx = iEvent.bunchCrossing();
  event = iEvent.id().event();

  if (isMC_){
    edm::Handle<GenEventInfoProduct> genEventHandle;
    iEvent.getByToken(genEventTag_, genEventHandle);
    genweight = genEventHandle->weight();

    edm::Handle< edm::View<PileupSummaryInfo> > pileups;
    iEvent.getByToken(muTag_, pileups);

    mu = pileups->at(1).getTrueNumInteractions();
  }
  else{
    mu = getAvgPU( run, lumi );
  }

  //------------ Rho ------------//

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoTag_, rhoHandle);
  rho = *rhoHandle;

  //------------ Primary Vertices ------------//

  edm::Handle< edm::View<reco::Vertex> > primaryVertices;
  iEvent.getByToken(pvTag_, primaryVertices);

  nPV = primaryVertices->size();
  reco::Vertex pvtx = primaryVertices->at(0);

  //--------------Generated Particles-------------//

  if(isMC_){
    edm::Handle< edm::View<reco::GenParticle> > genParticles;
    iEvent.getByToken(genParticleTag_, genParticles);

    vector<pair<reco::GenParticle, int> > reducedGens;

    //cout << "#\tID\tstatus\td1\td2\tm1\tm2\t4 momentum" << endl;
    //cout << "--------------------------------------------------------------------"<<endl;

    for (int i=0, n=genParticles->size(); i<n; i++) {
      const reco::GenParticle& p = genParticles->at(i);

      /*cout << i << "\t" << p.pdgId() << "\t" << p.status() << "\t";
      if (p.numberOfDaughters() > 0) cout << p.daughterRef(0).key() << "\t";
      else cout << -1 << "\t";
      if (p.numberOfDaughters() > 1) cout << p.daughterRef(1).key() << "\t"; 
      else cout << -1 << "\t";
      if (p.numberOfMothers() > 0) cout << p.motherRef(0).key() << "\t"; 
      else cout << -1 << "\t";
      if (p.numberOfMothers() > 1) cout << p.motherRef(1).key() << "\t";
      else cout << -1 << "\t";
      cout << "( " << p.mass() << ", " << p.pt() << ", " << p.eta() << ", " << p.phi() << " )" << endl;*/

      int id = p.pdgId();
      int status = p.status();
      int nDaught = p.numberOfDaughters();

      if (id>1000000)  //Z'
        reducedGens.push_back(make_pair(p,i));

      if (fabs(id)==6 && 20<=status && status<30) //first t's
        reducedGens.push_back(make_pair(p,i));

      else if (fabs(id)==6 && nDaught==2){   //last t's
        reducedGens.push_back(make_pair(p,i));
        reducedGens.push_back(make_pair( genParticles->at(p.daughterRef(0).key()),p.daughterRef(0).key()) ); //b or W (first)
        reducedGens.push_back(make_pair( genParticles->at(p.daughterRef(1).key()),p.daughterRef(1).key()) ); //b or W (first)
      }

      else if (fabs(id)==24 && nDaught==2){   //W's
        reducedGens.push_back(make_pair( genParticles->at(p.daughterRef(0).key()),p.daughterRef(0).key()) ); //q or lep
        reducedGens.push_back(make_pair( genParticles->at(p.daughterRef(1).key()),p.daughterRef(1).key()) ); //q or lep
      }
    }
    //cout << endl << "STORED PARTICLES" << endl;
      
    nGen = reducedGens.size();

    for (int i=0; i<nGen; i++) {
      const reco::GenParticle& p = reducedGens[i].first;

      /*cout << reducedGens[i].second << "\t" << p.pdgId() << "\t" << p.status() << "\t";
      if (p.numberOfDaughters() > 0) cout << p.daughterRef(0).key() << "\t";
      else cout << -1 << "\t";
      if (p.numberOfDaughters() > 1) cout << p.daughterRef(1).key() << "\t";
      else cout << -1 << "\t";
      if (p.numberOfMothers() > 0) cout << p.motherRef(0).key() << "\t";
      else cout << -1 << "\t";
      if (p.numberOfMothers() > 1) cout << p.motherRef(1).key() << "\t";
      else cout << -1 << "\t";
      cout << "( " << p.mass() << ", " << p.pt() << ", " << p.eta() << ", " << p.phi() << " )" << endl;*/

      gen_status[i] = p.status();
      gen_PID[i] = p.pdgId();
      gen_pt[i]  = p.pt();
      gen_mass[i] = p.mass();
      gen_eta[i] = p.eta();
      gen_phi[i] = p.phi();
      gen_index[i] = reducedGens[i].second;
      //gen_nMothers[i] = p.numberOfMothers();
      //gen_nDaughters[i] = p.numberOfDaughters();

      if (p.numberOfMothers() > 0) gen_mother0[i] = p.motherRef(0).key();
      else gen_mother0[i] = -1;
      if (p.numberOfMothers() > 1) gen_mother1[i] = p.motherRef(1).key();
      else gen_mother1[i] = -1;
    }
    //cout << endl;

    edm::Handle< edm::View<reco::GenJet> > genJets;
    iEvent.getByToken(genJetTag_, genJets);

    nGenJet = 0;

    for (int i=0, n=genJets->size(); i<n; i++){
      if (genJets->at(i).pt() < 15) continue;

      const reco::GenJet& jet = genJets->at(i);

      genJet_pt[nGenJet]  = jet.pt();
      genJet_eta[nGenJet] = jet.eta();
      genJet_phi[nGenJet] = jet.phi();
      genJet_mass[nGenJet] = jet.mass();
      genJet_area[nGenJet] = jet.jetArea();
      genJet_nDaught[nGenJet] = jet.numberOfDaughters();

      nGenJet++;
    }
  }

  //--------------Muons-------------//

  nMuon = 0;

  for (int i=0, n=muons->size(); i<n; i++){
    //save only loose muons above 15 gev
    if ( muons->at(i).pt() < 15 || !muons->at(i).isLooseMuon() ) continue;
    const pat::Muon& muon = muons->at(i);

    muon_isGlob[nMuon] = muon.isGlobalMuon();
    muon_charge[nMuon] = muon.charge();
    muon_pt[nMuon] = muon.pt();
    muon_eta[nMuon] = muon.eta();
    muon_phi[nMuon] = muon.phi();

    muon_D0[nMuon] = muon.muonBestTrack()->dxy(pvtx.position());
    muon_Dz[nMuon] = muon.muonBestTrack()->dz(pvtx.position());
    muon_tspm[nMuon] = muon.combinedQuality().chi2LocalPosition;
    muon_kinkf[nMuon] = muon.combinedQuality().trkKink;
    muon_ftrackhits[nMuon] = muon.innerTrack()->validFraction();

    muon_segcom[nMuon] = muon::segmentCompatibility(muon);

    if (muon_isGlob[nMuon]) muon_chi2[nMuon] = muon.globalTrack()->normalizedChi2();
    else                muon_chi2[nMuon] = 99;

    muon_IsMediumID[nMuon] = muon.isMediumMuon();
    muon_IsTightID[nMuon] = muon.isTightMuon(pvtx);   

    nMuon++;
  }

  //------------ Electrons ------------//

  edm::Handle<edm::ValueMap<vid::CutFlowResult> > Loose_id_decisions;
  iEvent.getByToken(eleLooseIdMapToken_, Loose_id_decisions);

  edm::Handle<edm::ValueMap<vid::CutFlowResult> > Medium_id_decisions;
  iEvent.getByToken(eleMediumIdMapToken_, Medium_id_decisions);

  edm::Handle<edm::ValueMap<vid::CutFlowResult> > Tight_id_decisions;
  iEvent.getByToken(eleTightIdMapToken_, Tight_id_decisions);

  nEle = 0;

  for (int i=0, n=electrons->size(); i<n; i++){
    const pat::Electron& ele = electrons->at(i);
    const Ptr<pat::Electron> elPtr(electrons, i);

    //save only veto electrons (without isocut) above 15 gev
    if (ele.pt() < 15 || !(*Veto_id_decisions)[elPtr].getCutFlowResultMasking(isoCut).cutFlowPassed()) continue;


    ele_charge[nEle] = ele.charge();
    ele_pt[nEle] = ele.pt();
    ele_eta[nEle] = ele.eta();
    ele_phi[nEle] = ele.phi();

    ele_D0[nEle] = ele.gsfTrack()->dxy(pvtx.position());
    ele_Dz[nEle] = ele.gsfTrack()->dz(pvtx.position());
    ele_sigmaIetaIeta[nEle] = ele.full5x5_sigmaIetaIeta();

    if (ele.superCluster().isNonnull()) ele_etaSupClust[nEle] = ele.superCluster()->eta();
    else ele_etaSupClust[nEle] = -99;
    if ( ele.superCluster().isNonnull() && ele.superCluster()->seed().isNonnull() )
      ele_dEtaSeed[nEle] = ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta();
    else ele_dEtaSeed[nEle] = -99;

    ele_dPhiIn[nEle] = ele.deltaPhiSuperClusterTrackAtVtx();
    ele_overEoverP[nEle] = fabs(1.0 - ele.eSuperClusterOverP()) / ele.ecalEnergy();

    ele_HE[nEle] = ele.hadronicOverEm();
    ele_missinghits[nEle] = ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

    reco::GsfElectron::PflowIsolationVariables pfIso = ele.pfIsolationVariables();
    float eA = ele_areas_.getEffectiveArea( fabs(ele_etaSupClust[nEle]) );
    ele_rcpiwec[nEle] = ( pfIso.sumChargedHadronPt + max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho) ) / ele_pt[nEle];

    //if ( beamspot.isValid() && conversions.isValid() )
    //  ele_passConv[nEle] = !ConversionTools::hasMatchedConversion(ele, conversions, beamspot->position());
    //else
    //  ele_passConv[nEle] = true;
    //Tried implementing the passConversion but didn't agree with official passConversion cut, so use below
    //const int passConvCut = 8;
    //ele_passConv[nEle] = (*Veto_id_decisions)[elPtr].getCutResultByIndex(passConvCut);

    ele_LooseID[nEle]  = (*Loose_id_decisions)[elPtr].getCutFlowResultMasking(isoCut).cutFlowPassed();
    ele_MediumID[nEle] = (*Medium_id_decisions)[elPtr].getCutFlowResultMasking(isoCut).cutFlowPassed();
    ele_TightID[nEle]  = (*Tight_id_decisions)[elPtr].getCutFlowResultMasking(isoCut).cutFlowPassed();

/*
0: MinPtCut_0
1: GsfEleSCEtaMultiRangeCut_0
2: GsfEleDEtaInSeedCut_0
3: GsfEleDPhiInCut_0
4: GsfEleFull5x5SigmaIEtaIEtaCut_0
5: GsfEleHadronicOverEMCut_0
6: GsfEleEInverseMinusPInverseCut_0
7: GsfEleEffAreaPFIsoCut_0
8: GsfEleConversionVetoCut_0
9: GsfEleMissingHitsCut_0    
*/
    nEle++;
  }

  //------------ MET ------------//

  edm::Handle< edm::View<pat::MET> > mets;
  iEvent.getByToken(metTag_, mets);

  const pat::MET& met = mets->at(0);

  if (isMC_){
    const reco::GenMET *genmet = met.genMET();

    genmet_pt = genmet->pt();
    genmet_px = genmet->px();
    genmet_py = genmet->py();
    genmet_sumet = genmet->sumEt();
    genmet_phi = genmet->phi();
  }

  met_pt = met.uncorPt();
  met_px = met.uncorPx();
  met_py = met.uncorPy();
  met_sumet = met.uncorSumEt();
  met_phi = met.uncorPhi();

  nMETUncert = METUNCERT;

  for (int i=0; i<nMETUncert; i++){

    pat::MET::METUncertainty uncert = static_cast<pat::MET::METUncertainty>(i);
    pat::MET::METCorrectionLevel level = pat::MET::Type1;

    met_shiftedpx[i] = met.shiftedPx(uncert, level);
    met_shiftedpy[i] = met.shiftedPy(uncert, level);
  }

  //------------ Triggers ------------//

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerResultsTag_, triggerResults);
  
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(prescalesTag_, triggerPrescales);
    
  if (triggerResults.isValid()) {
    const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults); 

    trig_prescale.clear(); trig_passed.clear(); trig_name.clear();

    for (int i=0; i<nTriggers; i++) {
      string myTrigger = triggers[i], name;
      int len = myTrigger.length(), index = -1;

      for (int i=0, n=triggerResults->size(); i<n; i++){
        name = triggerNames.triggerName(i);
        if (name.substr(0, len) == myTrigger) {index = i; break;}
      }
      if (index == -1) continue;

      trig_prescale.push_back( triggerPrescales->getPrescaleForIndex(index) );  
      trig_passed.push_back( triggerResults->accept(index) );
      trig_name.push_back( name );
    }
  }

  //------------ Fill Tree ------------//

  tree->Fill();
}


// ------------ method called once each job just after ending the event loop  ------------
void ZDilepton::endJob() {

  if (root_file !=0) {

    root_file->WriteObject(&totalEvts, "totalEvts");
    root_file->WriteObject(&filter_failed, "filter_failed");
    root_file->WriteObject(&dilep_cut, "dilep_cut");
    root_file->WriteObject(&leppt_cut, "leppt_cut");
    root_file->WriteObject(&jetpteta_cut, "jetpteta_cut");
    root_file->WriteObject(&met_cut, "met_cut");
    root_file->WriteObject(&dilepmass_cut, "dilepmass_cut");

    root_file->Write();
    delete root_file;
    root_file = 0;
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(ZDilepton);

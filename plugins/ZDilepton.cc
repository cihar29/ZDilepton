// -*- C++ -*-
//
// Package:    treemaker/ZDilepton
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
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <vector>

//root files
#include <TFile.h>
#include <TTree.h>

using namespace std;
using namespace edm;

const int MAXLEP = 10;
const int MAXGEN = 50;
const int MAXJET = 30;
const int MAXNPV = 50;

class ZDilepton : public edm::EDAnalyzer {
  public:
    explicit ZDilepton(const edm::ParameterSet&);

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    TFile* root_file;
    TTree* tree;

    ULong64_t event;
    int run, lumi, bx;

    float rho;
    int nPV;

    int nGen;
    int gen_status[MAXGEN], gen_PID[MAXGEN];
    float gen_pt[MAXGEN], gen_mass[MAXGEN], gen_eta[MAXGEN], gen_phi[MAXGEN];
    int nGenJet;
    float genJet_pt[MAXJET], genJet_eta[MAXJET], genJet_phi[MAXJET], genJet_area[MAXJET], genJet_nDaught[MAXJET];

    int nMuon;
    bool muon_isGlob[MAXLEP], muon_IsLooseID[MAXLEP], muon_IsMediumID[MAXLEP], muon_IsTightID[MAXLEP];
    int muon_type[MAXLEP], muon_charge[MAXLEP];
    float muon_pt[MAXLEP], muon_eta[MAXLEP], muon_phi[MAXLEP], muon_D0[MAXLEP], muon_Dz[MAXLEP];
    float muon_chi2[MAXLEP], muon_tspm[MAXLEP], muon_kinkf[MAXLEP], muon_segcom[MAXLEP];
    
    int nEle;
    int ele_charge[MAXLEP];
    bool ele_VetoID[MAXLEP], ele_LooseID[MAXLEP], ele_MediumID[MAXLEP], ele_TightID[MAXLEP];
    float ele_dEtaIn[MAXLEP], ele_dPhiIn[MAXLEP];
    float ele_pt[MAXLEP], ele_eta[MAXLEP], ele_phi[MAXLEP], ele_D0[MAXLEP], ele_Dz[MAXLEP];
    float ele_sigmaIetaIeta[MAXLEP], ele_dEtaSeed[MAXLEP], ele_dPhiSeed[MAXLEP], ele_HE[MAXLEP], ele_rcpiwec[MAXLEP], ele_overEoverP[MAXLEP];
    float ele_missinghits[MAXLEP];

    int nJet;
    float jet_pt[MAXJET], jet_eta[MAXJET], jet_phi[MAXJET], jet_area[MAXJET], jet_jec[MAXJET];
    float jet_nhf[MAXJET], jet_nef[MAXJET], jet_chf[MAXJET], jet_muf[MAXJET];
    float jet_elef[MAXJET], jet_numconst[MAXJET], jet_numneutral[MAXJET], jet_chmult[MAXJET];

    int nMET;
    float met_pt[MAXJET], met_eta[MAXJET], met_phi[MAXJET];

    TString RootFileName_;
    bool isMC_;

    EffectiveAreas ele_areas_;
    edm::EDGetTokenT< vector<PileupSummaryInfo> > muTag_;
    edm::EDGetTokenT<double> rhoTag_;
    edm::EDGetTokenT< vector<reco::Vertex> > pvTag_;
    edm::EDGetTokenT< vector<reco::GenParticle> > genParticleTag_;
    edm::EDGetTokenT< vector<reco::GenJet> > genJetTag_;
    edm::EDGetTokenT< vector<pat::Muon> > muonTag_;
    edm::EDGetTokenT< vector<pat::Electron> > electronTag_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
    edm::EDGetTokenT< vector<pat::Jet> > jetTag_;
    edm::EDGetTokenT< vector<pat::MET> > metTag_;
};

ZDilepton::ZDilepton(const edm::ParameterSet& iConfig):
  ele_areas_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath() )
{
  RootFileName_ = iConfig.getParameter<string>("RootFileName");
  isMC_ = iConfig.getParameter<bool>("isMC");
  rhoTag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoTag") );
  pvTag_ = consumes< vector<reco::Vertex> >( iConfig.getParameter<edm::InputTag>("pvTag") );
  genParticleTag_ =  consumes< vector<reco::GenParticle> >( iConfig.getParameter<edm::InputTag>("genParticleTag") );
  genJetTag_ =  consumes< vector<reco::GenJet> >( iConfig.getParameter<edm::InputTag>("genJetTag") );
  muonTag_ = consumes< vector<pat::Muon> >( iConfig.getParameter<edm::InputTag>("muonTag") );
  electronTag_ = consumes< vector<pat::Electron> >( iConfig.getParameter<edm::InputTag>("electronTag") );
  eleVetoIdMapToken_ = consumes<edm::ValueMap<bool> >( iConfig.getParameter<edm::InputTag>("eleVetoIdMapToken") );
  eleLooseIdMapToken_ = consumes<edm::ValueMap<bool> >( iConfig.getParameter<edm::InputTag>("eleLooseIdMapToken") );
  eleMediumIdMapToken_ = consumes<edm::ValueMap<bool> >( iConfig.getParameter<edm::InputTag>("eleMediumIdMapToken") );
  eleTightIdMapToken_ = consumes<edm::ValueMap<bool> >( iConfig.getParameter<edm::InputTag>("eleTightIdMapToken") );
  jetTag_ = consumes< vector<pat::Jet> >( iConfig.getParameter<edm::InputTag>("jetTag") );
  metTag_ = consumes< vector<pat::MET> >( iConfig.getParameter<edm::InputTag>("metTag") );
}

// ------------ method called once each job just before starting event loop  ------------
void  ZDilepton::beginJob() {

  root_file = new TFile(RootFileName_,"RECREATE");
  tree = new TTree("T","RECREATE");

  if(isMC_){
    tree->Branch("nGen", &nGen, "nGen/I");
    tree->Branch("gen_status",  gen_status, "gen_status[nGen]/I");
    tree->Branch("gen_PID",  gen_PID, "gen_PID[nGen]/I");
    tree->Branch("gen_pt", gen_pt, "gen_pt[nGen]/F");
    tree->Branch("gen_mass", gen_mass, "gen_mass[nGen]/F");
    tree->Branch("gen_eta", gen_eta, "gen_eta[nGen]/F");
    tree->Branch("gen_phi", gen_phi, "gen_phi[nGen]/F");

    tree->Branch("nGenJet", &nGenJet, "nGenJet/I");
    tree->Branch("genJet_pt", genJet_pt, "genJet_pt[nGenJet]/F");
    tree->Branch("genJet_eta", genJet_eta, "genJet_eta[nGenJet]/F");
    tree->Branch("genJet_phi", genJet_phi, "genJet_phi[nGenJet]/F");
    tree->Branch("genJet_area",  genJet_area, "genJet_area[nGenJet]/F");
    tree->Branch("genJet_nDaught",  genJet_nDaught, "genJet_nDaught[nGenJet]/F");
  }
  else{
    tree->Branch("run", &run, "run/I");
    tree->Branch("lumi", &lumi, "lumi/I");
    tree->Branch("bx", &bx, "bx/I");
    tree->Branch("event", &event, "event/l");
  }

  tree->Branch("rho",   &rho,   "rho/F");
  tree->Branch("nPV",     &nPV,    "nPV/I");

  tree->Branch("nMuon", &nMuon, "nMuon/I");
  tree->Branch("muon_charge", muon_charge, "muon_charge[nMuon]/I");
  tree->Branch("muon_type", muon_type, "muon_type[nMuon]/I");
  tree->Branch("muon_isGlob", muon_isGlob, "muon_isGlob[nMuon]/O");
  tree->Branch("muon_pt", muon_pt, "muon_pt[nMuon]/F");
  tree->Branch("muon_eta", muon_eta, "muon_eta[nMuon]/F");
  tree->Branch("muon_phi", muon_phi, "muon_phi[nMuon]/F");
  tree->Branch("muon_D0", muon_D0, "muon_D0[nMuon]/F");
  tree->Branch("muon_Dz", muon_Dz, "muon_Dz[nMuon]/F");
  tree->Branch("muon_chi2", muon_chi2, "muon_chi2[nMuon]/F");
  tree->Branch("muon_tspm", muon_tspm, "muon_tspm[nMuon]/F");
  tree->Branch("muon_kinkf", muon_kinkf, "muon_kinkf[nMuon]/F");
  if (!isMC_) tree->Branch("muon_segcom", muon_segcom, "muon_segcom[nMuon]/F");
  tree->Branch("muon_IsLooseID", muon_IsLooseID, "muon_IsLooseID[nMuon]/O");
  tree->Branch("muon_IsMediumID", muon_IsMediumID, "muon_IsMediumID[nMuon]/O");
  tree->Branch("muon_IsTightID", muon_IsTightID, "muon_IsTightID[nMuon]/O");

  tree->Branch("nEle", &nEle, "nEle/I");
  tree->Branch("ele_charge", ele_charge, "ele_charge[nEle]/I");
  tree->Branch("ele_pt", ele_pt, "ele_pt[nEle]/F");
  tree->Branch("ele_eta", ele_eta, "ele_eta[nEle]/F");
  tree->Branch("ele_phi", ele_phi, "ele_phi[nEle]/F");
  tree->Branch("ele_VetoID", ele_VetoID , "ele_VetoID/O");
  tree->Branch("ele_LooseID", ele_LooseID , "ele_LooseID/O");
  tree->Branch("ele_MediumID", ele_MediumID , "ele_MediumID/O");
  tree->Branch("ele_TightID", ele_TightID , "ele_TightID/O");
  tree->Branch("ele_D0", ele_D0, "ele_D0[nEle]/F");
  tree->Branch("ele_Dz", ele_Dz, "ele_Dz[nEle]/F");
  tree->Branch("ele_sigmaIetaIeta", ele_sigmaIetaIeta, "ele_sigmaIetaIeta[nEle]/F");
  tree->Branch("ele_dEtaSeed", ele_dEtaSeed, "ele_dEtaSeed[nEle]/F");
  tree->Branch("ele_dPhiSeed", ele_dPhiSeed, "ele_dPhiSeed[nEle]/F");
  tree->Branch("ele_dEtaIn", ele_dEtaIn, "ele_dEtaIn[nEle]/F");
  tree->Branch("ele_dPhiIn", ele_dPhiIn, "ele_dPhiIn[nEle]/F");
  tree->Branch("ele_overEoverP", ele_overEoverP, "ele_overEoverP[nEle]/F");
  if (!isMC_){
    tree->Branch("ele_HE", ele_HE, "ele_HE[nEle]/F");
    tree->Branch("ele_rcpiwec", ele_rcpiwec, "ele_rcpiwec[nEle]/F");
    tree->Branch("ele_missinghits", ele_missinghits, "ele_missinghits[nEle]/F");
  }

  tree->Branch("nJet", &nJet, "nJet/I");
  tree->Branch("jet_pt", jet_pt, "jet_pt[nJet]/F");
  tree->Branch("jet_eta", jet_eta, "jet_eta[nJet]/F");
  tree->Branch("jet_phi", jet_phi, "jet_phi[nJet]/F");
  tree->Branch("jet_area", jet_area, "jet_area[nJet]/F");
  tree->Branch("jet_jec", jet_jec, "jet_jec[nJet]/F");

  tree->Branch("jet_nhf", jet_nhf, "jet_nhf[nJet]/F");
  tree->Branch("jet_nef", jet_nef, "jet_nef[nJet]/F");
  tree->Branch("jet_chf", jet_chf, "jet_chf[nJet]/F");
  tree->Branch("jet_muf", jet_muf, "jet_muf[nJet]/F");
  tree->Branch("jet_elef", jet_elef, "jet_elef[nJet]/F");
  tree->Branch("jet_numconst", jet_numconst, "jet_numconst[nJet]/F");
  tree->Branch("jet_numneutral", jet_numneutral, "jet_numneutral[nJet]/F");
  tree->Branch("jet_chmult", jet_chmult, "jet_chmult[nJet]/F");

  tree->Branch("nMET", &nMET, "nMET/I");
  tree->Branch("met_pt", met_pt, "met_pt[nMET]/F");
  tree->Branch("met_eta", met_eta, "met_eta[nMET]/F");
  tree->Branch("met_phi", met_phi, "met_phi[nMET]/F");
}

// ------------ method called for each event  ------------
void ZDilepton::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //------------ Event Info ------------//

  if(!isMC_){
    run = int(iEvent.id().run());
    lumi = int(iEvent.getLuminosityBlock().luminosityBlock());
    bx = iEvent.bunchCrossing();
    event = iEvent.id().event();
  }

  //------------ Rho ------------//

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoTag_, rhoHandle);
  rho = *rhoHandle;

  //------------ Primary Vertices ------------//

  edm::Handle< vector<reco::Vertex> > primaryVertices;
  iEvent.getByToken(pvTag_, primaryVertices);

  nPV = primaryVertices->size();
  reco::Vertex pvtx = primaryVertices->at(0);

  //--------------Generated Particles-------------//

  if(isMC_){
    edm::Handle< vector<reco::GenParticle> > genParticles;
    iEvent.getByToken(genParticleTag_, genParticles);

    nGen = genParticles->size();

    for (int i=0; i<nGen; i++){
      reco::GenParticle p = genParticles->at(i);

      gen_status[i] = p.status();
      gen_PID[i] = p.pdgId();
      gen_pt[i]  = p.pt();
      gen_mass[i] = p.mass();
      gen_eta[i] = p.eta();
      gen_phi[i] = p.phi();
    }

    edm::Handle< vector<reco::GenJet> > genJets;
    iEvent.getByToken(genJetTag_, genJets);

    nGenJet = genJets->size();

    for (int i=0; i<nGenJet; i++){
      reco::GenJet jet = genJets->at(i);

      genJet_pt[i]  = jet.pt();
      genJet_eta[i] = jet.eta();
      genJet_phi[i] = jet.phi();
      genJet_area[i] = jet.jetArea();
      genJet_nDaught[i] = jet.numberOfDaughters();
    }
  }

  //--------------Muons-------------//

  edm::Handle< vector<pat::Muon> > muons;
  iEvent.getByToken(muonTag_, muons);

  nMuon = muons->size();

  for (int i=0; i<nMuon; i++){
    pat::Muon muon = muons->at(i);

    muon_isGlob[i] = muon.isGlobalMuon();
    muon_charge[i] = muon.charge();
    muon_type[i] = muon.type();
    muon_pt[i] = muon.pt();
    muon_eta[i] = muon.eta();
    muon_phi[i] = muon.phi();

    muon_D0[i] = muon.muonBestTrack()->dxy(pvtx.position());
    muon_Dz[i] = muon.muonBestTrack()->dz(pvtx.position());
    muon_tspm[i] = muon.combinedQuality().chi2LocalPosition;
    muon_kinkf[i] = muon.combinedQuality().trkKink;

    if (!isMC_) muon_segcom[i] = muon::segmentCompatibility(muon);

    if (muon_isGlob[i]) muon_chi2[i] = muon.globalTrack()->normalizedChi2();
    else                muon_chi2[i] = -1;

    muon_IsLooseID[i] = muon.isLooseMuon();
    muon_IsMediumID[i] = muon.isMediumMuon();
    muon_IsTightID[i] = muon.isTightMuon(pvtx);   
  }

  //------------ Electrons ------------//

  edm::Handle< vector<pat::Electron> > electrons;
  iEvent.getByToken(electronTag_, electrons);

  edm::Handle<edm::ValueMap<bool> > Veto_id_decisions;
  iEvent.getByToken(eleVetoIdMapToken_, Veto_id_decisions);

  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  iEvent.getByToken(eleLooseIdMapToken_, loose_id_decisions);

  edm::Handle<edm::ValueMap<bool> > Medium_id_decisions;
  iEvent.getByToken(eleMediumIdMapToken_, Medium_id_decisions);

  edm::Handle<edm::ValueMap<bool> > Tight_id_decisions;
  iEvent.getByToken(eleTightIdMapToken_, Tight_id_decisions);

  nEle = electrons->size();

  for (int i=0; i<nEle; i++){
    pat::Electron ele = electrons->at(i);
    const Ptr<pat::Electron> elPtr(electrons, i);

    ele_charge[i] = ele.charge();
    ele_pt[i] = ele.pt();
    ele_eta[i] = ele.eta();
    ele_phi[i] = ele.phi();

    ele_D0[i] = ele.gsfTrack()->dxy(pvtx.position());
    ele_Dz[i] = ele.gsfTrack()->dz(pvtx.position());
    ele_sigmaIetaIeta[i] = ele.full5x5_sigmaIetaIeta();

    if ( ele.superCluster().isNonnull() && ele.superCluster()->seed().isNonnull() )
      ele_dEtaSeed[i] = ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta();
    else ele_dEtaSeed[i] = -99;

    ele_dPhiSeed[i] = abs(ele.deltaPhiSuperClusterTrackAtVtx());

    ele_dEtaIn[i] = ele.deltaEtaSuperClusterTrackAtVtx();
    ele_dPhiIn[i] = ele.deltaPhiSuperClusterTrackAtVtx();
    ele_overEoverP[i] = abs(1.0 - ele.eSuperClusterOverP()) / ele.ecalEnergy();

    if (!isMC_){
      ele_HE[i] = ele.hadronicOverEm();
      ele_missinghits[i] = ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

      reco::GsfElectron::PflowIsolationVariables pfIso = ele.pfIsolationVariables();
      float abseta =  abs(ele.superCluster()->eta());
      float eA = ele_areas_.getEffectiveArea(abseta);
      ele_rcpiwec[i] = ( pfIso.sumChargedHadronPt + max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho) ) / ele_pt[i];
    }

    ele_VetoID[i]   = (*Veto_id_decisions)[ elPtr ];
    ele_LooseID[i]  = (*loose_id_decisions)[ elPtr ];
    ele_MediumID[i]  = (*Medium_id_decisions)[ elPtr ];
    ele_TightID[i]  = (*Tight_id_decisions)[ elPtr ];
  }

  //------------ Jets ------------//

  edm::Handle< vector<pat::Jet> > jets;
  iEvent.getByToken(jetTag_, jets);

  nJet = jets->size();

  for (int i=0; i<nJet; i++){
    pat::Jet jet = jets->at(i);

    jet_pt[i] = jet.pt();
    jet_eta[i] = jet.eta();
    jet_phi[i] = jet.phi();
    jet_area[i] = jet.jetArea();
    jet_jec[i] = jet.jecFactor(0);

    jet_nhf[i] = jet.neutralHadronEnergyFraction();
    jet_nef[i] = jet.neutralEmEnergyFraction();
    jet_chf[i] = jet.chargedHadronEnergyFraction();
    jet_muf[i] = jet.muonEnergyFraction();
    jet_elef[i] = jet.chargedEmEnergyFraction();
    jet_numconst[i] = jet.chargedMultiplicity()+jet.neutralMultiplicity();
    jet_numneutral[i] =jet.neutralMultiplicity();
    jet_chmult[i] = jet.chargedMultiplicity();
  }

  //------------ MET ------------//

  edm::Handle< vector<pat::MET> > mets;
  iEvent.getByToken(metTag_, mets);

  nMET = mets->size();

  for (int i=0; i<nMET; i++){
    pat::MET met = mets->at(i);

    met_pt[i] = met.pt();
    met_eta[i] = met.eta();
    met_phi[i] = met.phi();
  }

  //------------ Fill Tree ------------//

  tree->Fill();
}


// ------------ method called once each job just after ending the event loop  ------------
void ZDilepton::endJob() {

  if (root_file !=0) {

    root_file->Write();
    delete root_file;
    root_file = 0;
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(ZDilepton);


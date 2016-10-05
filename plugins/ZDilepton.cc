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
// Original Author:  charles harrington
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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include <vector>

//root files
#include <TFile.h>
#include <TTree.h>

using namespace std;

const int MAXLEP = 10;
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

    int nPVall, nPV;
    float pv_ndof[MAXNPV], pv_z[MAXNPV], pv_rho[MAXNPV];

    int nMuon;
    bool muon_isGlob[MAXLEP];
    float muon_pt[MAXLEP], muon_eta[MAXLEP], muon_phi[MAXLEP], muon_chi2[MAXLEP], muon_tspm[MAXLEP], muon_kinkf[MAXLEP], muon_segcom[MAXLEP];

    //int nElectron;
    //float muon_pt[MAXLEP], muon_eta[MAXLEP], muon_phi[MAXLEP];

    TString RootFileName_;

    edm::EDGetTokenT< vector<PileupSummaryInfo> > muTag_;
    edm::EDGetTokenT<double> rhoTag_;
    edm::EDGetTokenT< vector<reco::Vertex> > pvTag_;
    edm::EDGetTokenT< vector<pat::Muon> > muonTag_;
};

ZDilepton::ZDilepton(const edm::ParameterSet& iConfig)
{
  RootFileName_ = iConfig.getParameter<string>("RootFileName");
  rhoTag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoTag") );
  pvTag_ = consumes< vector<reco::Vertex> >( iConfig.getParameter<edm::InputTag>("pvTag") );
  muonTag_ = consumes< vector<pat::Muon> >( iConfig.getParameter<edm::InputTag>("muonTag") );
}

// ------------ method called once each job just before starting event loop  ------------
void  ZDilepton::beginJob() {

  root_file = new TFile(RootFileName_,"RECREATE");
  tree = new TTree("T","Offset Tree");

  tree->Branch("run", &run, "run/I");
  tree->Branch("lumi", &lumi, "lumi/I");
  tree->Branch("bx", &bx, "bx/I");
  tree->Branch("event", &event, "event/l");

  tree->Branch("rho",   &rho,   "rho/F");

  tree->Branch("nPVall",  &nPVall, "nPVall/I");
  tree->Branch("nPV",     &nPV,    "nPV/I");
  tree->Branch("pv_ndof", pv_ndof, "pv_ndof[nPVall]/F");
  tree->Branch("pv_z",    pv_z,    "pv_z[nPVall]/F");
  tree->Branch("pv_rho",  pv_rho,  "pv_rho[nPVall]/F");

  tree->Branch("nMuon", &nMuon, "nMuon/I");
  tree->Branch("muon_isGlob", muon_isGlob, "muon_isGlob[nMuon]/O");
  tree->Branch("muon_pt", muon_pt, "muon_pt[nMuon]/F");
  tree->Branch("muon_eta", muon_eta, "muon_eta[nMuon]/F");
  tree->Branch("muon_phi", muon_phi, "muon_phi[nMuon]/F");
  tree->Branch("muon_chi2", muon_chi2, "muon_chi2[nMuon]/F");
  tree->Branch("muon_tspm", muon_tspm, "muon_tspm[nMuon]/F");
  tree->Branch("muon_kinkf", muon_kinkf, "muon_kinkf[nMuon]/F");
  tree->Branch("muon_segcom", muon_segcom, "muon_segcom[nMuon]/F");
}

// ------------ method called for each event  ------------
void ZDilepton::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

//------------ Event Info ------------//

  run = int(iEvent.id().run());
  lumi = int(iEvent.getLuminosityBlock().luminosityBlock());
  bx = iEvent.bunchCrossing();
  event = iEvent.id().event();

//------------ Rho ------------//

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoTag_, rhoHandle);
  rho = *rhoHandle;

//------------ Primary Vertices ------------//

  edm::Handle< vector<reco::Vertex> > primaryVertices;
  iEvent.getByToken(pvTag_, primaryVertices);

  nPVall = primaryVertices->size();
  nPV = 0;

  for (int i=0; i<nPVall; i++){
    reco::Vertex pv = primaryVertices->at(i);

    pv_ndof[i] = pv.ndof();
    pv_z[i] = pv.z();
    pv_rho[i] = pv.position().rho();

    if( !pv.isFake() && pv_ndof[i] > 4 && pv_z[i] <= 24 && pv_rho[i] <= 2 )
      nPV++;
  }

//------------ Muons ------------//

  edm::Handle< vector<pat::Muon> > muons;
  iEvent.getByToken(muonTag_, muons);

  nMuon = muons->size();

  for (int i=0; i<nMuon; i++){
    pat::Muon muon = muons->at(i);

    muon_isGlob[i] = muon.isGlobalMuon();
    muon_pt[i] = muon.pt();
    muon_eta[i] = muon.eta();
    muon_phi[i] = muon.phi();
    muon_tspm[i] = muon.combinedQuality().chi2LocalPosition;
    muon_kinkf[i] = muon.combinedQuality().trkKink;
    muon_segcom[i] = muon::segmentCompatibility(muon);

    if (muon_isGlob[i]) muon_chi2[i] = muon.globalTrack()->normalizedChi2();
    else                muon_chi2[i] = -1;
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

//charles harrington and bahareh roozbahani
//execute as analyze pars.txt mc_weights.txt

#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"

#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <time.h>

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

using namespace std;

void FillHist1D(const TString& histName, const Double_t& value, const double& weight);
void setPars(const string& parFile); 
void setWeight(const string& parFile); 
bool isMediumMuonBCDEF(const bool& isGlob, const float& chi2, const float& tspm, const float& kinkf, const float& segcom, const float& ftrackhits);
bool sortJetPt(const pair<int, float>& jet1, const pair<int, float>& jet2){ return jet1.second > jet2.second; }
void FillHists(const TString& prefix, const int& nEle, const int& nGoodEle, const int& nMuon, const int& nGoodMuon, const int& nJet, const int& nGoodJet,
               const TLorentzVector& lep0, const TLorentzVector& lep1, const float& dilepmass, const float& lepept, const float& lepmpt,
               const float& rmin0, const float& rmin1, const float& rl0l1, const float& rl0cleanj, const float& rl1cleanj, const float& lep0perp, const float& lep1perp,
               const TLorentzVector& jet0, const TLorentzVector& jet1, const float& jet0btag, const float& jet1btag, const int& nbtag,
               const float& hT, const float& met_pt, const float& met_corrpt, const float& sT, const float& sT_met,
               const float& rbal, const float& rabl, const float& minjet0pt, const float& minjet1pt, const float& cleanjet0pt, const float& cleanjet1pt,
               const float& masslmin0, const float& masslmin1, const float& masslljjm, const float& deta_lep, const float& deta_lepJet);

map<TString, TH1*> m_Histos1D;

//parameters- edit in pars.txt
bool ISMC;
TString inName, outName, muTrigSfName, muIdSfName, muTrackSfName, eRecoSfName, eIdSfName;
string channel, jet_type;
vector<string> eras;
double weight0, weight;

const int MAXJET = 50;
const int MAXLEP = 20;
const int MAXGEN = 20;
const float MUONMASS = 0.10566;
const float ELEMASS = 0.;

int main(int argc, char* argv[]){

  if (argc == 1) { cout << "Please provide a parameter file." << endl; return -1; }

  string parFile = argv[1];
  setPars(parFile);

  weight0 = -1;
  if (argc == 3)    { string wFile = argv[2]; setWeight(wFile); }
  if (weight0 == -1) { cout << "Weight set to 1" << endl; weight0 = 1.; }
  else                cout << "Weight set to " << weight0 << endl;

  //Jet Corrections//

  map<string, JetCorrectorParameters*> ResJetPars, L3JetPars, L2JetPars, L1JetPars;
  map<string, vector<JetCorrectorParameters> > jetPars, jetL1Pars;
  map<string, FactorizedJetCorrector*> jetCorrectors, jetL1Correctors;
  map<string, pair<int, int> > m_IOV;

  cout << endl << "Using eras: " << endl;
  for(vector<string>::iterator i_era = eras.begin(); i_era != eras.end(); ++i_era) {
    string era = *i_era;

    ResJetPars[era] = new JetCorrectorParameters(era + "/" + era + "_L2L3Residual_" + jet_type + ".txt");
    L3JetPars[era] = new JetCorrectorParameters(era + "/" + era + "_L3Absolute_" + jet_type + ".txt");
    L2JetPars[era] = new JetCorrectorParameters(era + "/" + era + "_L2Relative_" + jet_type + ".txt");
    L1JetPars[era] = new JetCorrectorParameters(era + "/" + era + "_L1FastJet_" + jet_type + ".txt");

    jetPars[era].push_back( *L1JetPars[era] );
    jetPars[era].push_back( *L2JetPars[era] );
    jetPars[era].push_back( *L3JetPars[era] );
    jetPars[era].push_back( *ResJetPars[era] );

    jetCorrectors[era] = new FactorizedJetCorrector( jetPars[era] );

    jetL1Pars[era].push_back( *L1JetPars[era] );
    jetL1Correctors[era] = new FactorizedJetCorrector( jetL1Pars[era] );

    if (ISMC) cout << era << endl;
    else {
      TString tera = era.data();

      if ( tera.Contains("BCDV", TString::kIgnoreCase) ) { cout << "BCD "; m_IOV[era] = make_pair(1, 276811); }
      else if ( tera.Contains("EFV", TString::kIgnoreCase) ) { cout << "EF "; m_IOV[era] = make_pair(276831, 278801); }
      else if ( tera.Contains("GV", TString::kIgnoreCase) ) { cout << "G "; m_IOV[era] = make_pair(278802, 280385); }
      else if ( tera.Contains("HV", TString::kIgnoreCase) ) { cout << "H "; m_IOV[era] = make_pair(280919, 300000); }

      cout << "[" << m_IOV[era].first << ", " << m_IOV[era].second << "]" << endl;
    }   
  }

  //Open Files//

  TFile* inFile = TFile::Open(inName);

  TTree* T = (TTree*) inFile->Get("T");
  Long64_t nEntries = T->GetEntries();
  cout << endl << nEntries << " Events" << endl;
  cout << "Processing " + inName << endl;
  cout << "Channel: " + channel << endl;

  //SF Files//

  TH2F* muTrigSfHist=0, *muIdSfHist=0, *eRecoSfHist=0, *eIdSfHist=0;
  TGraphAsymmErrors* muTrackSfGraph=0;
  
  int muTrig_pT=0, muId_pT=0, eReco_pT=0, eId_pT=0;

  if (ISMC) {
    TFile* muTrigSfFile = TFile::Open(muTrigSfName);
    muTrigSfHist = (TH2F*) muTrigSfFile->Get("Mu50_OR_TkMu50_PtEtaBins/abseta_pt_ratio");

    TFile* muIdSfFile = TFile::Open(muIdSfName);
    muIdSfHist = (TH2F*) muIdSfFile->Get("MC_NUM_MediumID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");

    TFile* muTrackSfFile = TFile::Open(muTrackSfName);
    muTrackSfGraph = (TGraphAsymmErrors*) muTrackSfFile->Get("ratio_eff_eta3_dr030e030_corr");

    TFile* eRecoSfFile = TFile::Open(eRecoSfName);
    eRecoSfHist = (TH2F*) eRecoSfFile->Get("EGamma_SF2D");

    TFile* eIdSfFile = TFile::Open(eIdSfName);
    eIdSfHist = (TH2F*) eIdSfFile->Get("EGamma_SF2D");

    muTrig_pT = muTrigSfHist->GetYaxis()->GetBinCenter(muTrigSfHist->GetYaxis()->GetNbins());
    muId_pT = muIdSfHist->GetYaxis()->GetBinCenter(muIdSfHist->GetYaxis()->GetNbins());
    eReco_pT = eRecoSfHist->GetYaxis()->GetBinCenter(eRecoSfHist->GetYaxis()->GetNbins());
    eId_pT = eIdSfHist->GetYaxis()->GetBinCenter(eIdSfHist->GetYaxis()->GetNbins());
  }

  //Skims//

  enum Cuts{
    countEvts, countDilep, countLeppt, countDilepmass, countJetpteta, countMet,
    channelCut, trigCut, lepkinCut, signCut, thirdLepCut, dilepmassCut, dilepVetoCut, ptrelCut, jetCut1, onebtagCut1jet, jetCut2, onebtagCut2jets, twobtagsCut2jets, numCuts
  };
  vector<pair<string, double> > v_cuts(numCuts);

  v_cuts[countEvts]=make_pair("Initial",0.); v_cuts[countDilep]=make_pair("Dilepton selection",0.); v_cuts[countLeppt]=make_pair("Lepton Pt Cut",0.);
  v_cuts[countDilepmass]=make_pair("Dilepton Mass Cut",0.); v_cuts[countJetpteta]=make_pair("Leading Jet Pt/eta cut",0.);
  v_cuts[countMet]=make_pair("MET Filters",0.); v_cuts[channelCut]=make_pair("Correct Channel",0.); v_cuts[signCut]=make_pair("Opposite Lepton Sign",0.);
  v_cuts[trigCut]=make_pair("HLT Trigger",0.); v_cuts[lepkinCut]=make_pair("Lepton kinematics cut",0.); v_cuts[thirdLepCut]=make_pair("Third lepton cut",0.);
  v_cuts[dilepmassCut]=make_pair("Dilepton mass cut",0.); v_cuts[jetCut1]=make_pair("1 Jet, 0 btags",0.); v_cuts[ptrelCut]=make_pair("pTrel cut",0.);
  v_cuts[jetCut2]=make_pair("2 Jets, 0 btags",0.); v_cuts[dilepVetoCut]=make_pair("Z-mass veto",0.);
  v_cuts[onebtagCut1jet]=make_pair("1 Jet, 1 btag",0.); v_cuts[onebtagCut2jets]=make_pair("2 Jets, 1 btag",0.);
  v_cuts[twobtagsCut2jets]=make_pair("2 Jets, 2 btags",0.);

  TIter nextkey(inFile->GetListOfKeys());
  TKey* key;
  while ( (key = (TKey*)nextkey()) ) {
    TString keyname = key->GetName();

    if (keyname.EqualTo("totalEvts"))          v_cuts[countEvts].second += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname.EqualTo("dilep_cut"))     v_cuts[countDilep].second += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname.EqualTo("leppt_cut"))     v_cuts[countLeppt].second += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname.EqualTo("dilepmass_cut")) v_cuts[countDilepmass].second += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname.EqualTo("jetpteta_cut"))  v_cuts[countJetpteta].second += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname.EqualTo("met_cut"))       v_cuts[countMet].second += weight0 * (*(vector<int>*)key->ReadObj())[0];
    //else if (keyname.EqualTo("filter_failed"))
  }
  if ( (int) (v_cuts[countMet].second + 0.5) != (int) (weight0 * nEntries + 0.5) ) { cout << "hadd added incorrectly." << endl; return -1; }

  //Histograms//

  int nDirs = 5;
  for (int i=0; i<nDirs; i++) {
    TString hname = Form("%i_nJet",i);
    m_Histos1D[hname] = new TH1F(hname,hname,MAXJET,0,MAXJET);
    hname = Form("%i_nGoodJet",i);
    m_Histos1D[hname] = new TH1F(hname,hname,MAXJET,0,MAXJET);
    hname = Form("%i_nJetDiff",i);
    m_Histos1D[hname] = new TH1F(hname,hname,MAXJET,0,MAXJET);
    hname = Form("%i_jet0pt",i);
    m_Histos1D[hname] = new TH1F(hname,hname,200,0,2000);
    hname = Form("%i_jet1pt",i);
    m_Histos1D[hname] = new TH1F(hname,hname,200,0,2000);
    hname = Form("%i_jet0eta",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,-5,5);
    hname = Form("%i_jet1eta",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,-5,5);
    hname = Form("%i_jet0phi",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,-4,4);
    hname = Form("%i_jet1phi",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,-4,4);
    hname = Form("%i_jet0btag",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,1);
    hname = Form("%i_jet1btag",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,1);
    hname = Form("%i_nbtag",i);
    m_Histos1D[hname] = new TH1F(hname,hname,5,0,5);
    hname = Form("%i_jethT",i);
    m_Histos1D[hname] = new TH1F(hname,hname,200,0,2000);

    hname = Form("%i_minjet0pt",i);
    m_Histos1D[hname] = new TH1F(hname,hname,200,0,2000);
    hname = Form("%i_minjet1pt",i);
    m_Histos1D[hname] = new TH1F(hname,hname,200,0,2000);
    hname = Form("%i_cleanjet0pt",i);
    m_Histos1D[hname] = new TH1F(hname,hname,200,0,2000);
    hname = Form("%i_cleanjet1pt",i);
    m_Histos1D[hname] = new TH1F(hname,hname,200,0,2000);
    hname = Form("%i_masslmin0",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,1000);
    hname = Form("%i_masslmin1",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,1000);
    hname = Form("%i_masslljjm",i);
    m_Histos1D[hname] = new TH1F(hname,hname,500,0,10000);

    hname = Form("%i_nEle",i);
    m_Histos1D[hname] = new TH1F(hname,hname,MAXLEP,0,MAXLEP);
    hname = Form("%i_nEleDiff",i);
    m_Histos1D[hname] = new TH1F(hname,hname,MAXLEP,0,MAXLEP);
    hname = Form("%i_nMuon",i);
    m_Histos1D[hname] = new TH1F(hname,hname,MAXLEP,0,MAXLEP);
    hname = Form("%i_nMuonDiff",i);
    m_Histos1D[hname] = new TH1F(hname,hname,MAXLEP,0,MAXLEP);
    hname = Form("%i_nGoodEle",i);
    m_Histos1D[hname] = new TH1F(hname,hname,MAXLEP,0,MAXLEP);
    hname = Form("%i_nGoodMuon",i);
    m_Histos1D[hname] = new TH1F(hname,hname,MAXLEP,0,MAXLEP);
    hname = Form("%i_lep0pt",i);
    m_Histos1D[hname] = new TH1F(hname,hname,500,0,5000);
    hname = Form("%i_lep0eta",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,-5,5);
    hname = Form("%i_lep1pt",i);
    m_Histos1D[hname] = new TH1F(hname,hname,500,0,5000);
    hname = Form("%i_lep1eta",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,-5,5);
    hname = Form("%i_lep0phi",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,-4,4);
    hname = Form("%i_lep1phi",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,-4,4);
    hname = Form("%i_dilepmass",i);
    m_Histos1D[hname] = new TH1F(hname,hname,500,0,2000);
    hname = Form("%i_lepept",i);
    m_Histos1D[hname] = new TH1F(hname,hname,500,0,1000);
    hname = Form("%i_lepmpt",i);
    m_Histos1D[hname] = new TH1F(hname,hname,500,0,1000);
    hname = Form("%i_deta_lep",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,-5,5);
    hname = Form("%i_deta_lepJet",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,-5,5);

    hname = Form("%i_muonD0",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,-1,1);
    hname = Form("%i_muonDz",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,-1,1);
    hname = Form("%i_rmin0",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,5);
    hname = Form("%i_rmin1",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,5);
    hname = Form("%i_rbl",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,5);
    hname = Form("%i_rl0l1",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,5);
    hname = Form("%i_rl0cleanj",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,5);
    hname = Form("%i_rl1cleanj",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,5);
    hname = Form("%i_lep0perp",i);
    m_Histos1D[hname] = new TH1F(hname,hname,400,0,2000);
    hname = Form("%i_lep1perp",i);
    m_Histos1D[hname] = new TH1F(hname,hname,400,0,2000);

    hname = Form("%i_metpt",i);
    m_Histos1D[hname] = new TH1F(hname,hname,400,0,2000);
    hname = Form("%i_metcorrpt",i);
    m_Histos1D[hname] = new TH1F(hname,hname,400,0,2000);
    hname = Form("%i_sT",i);
    m_Histos1D[hname] = new TH1F(hname,hname,300,0,3000);
    hname = Form("%i_sT_met",i);
    m_Histos1D[hname] = new TH1F(hname,hname,300,0,3000);
  }

  TString hname = "muTrigSf";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,2);
  hname = "muIdSf";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,2);
  hname = "muTrackSf";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,2);
  hname = "eRecoSf";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,2);
  hname = "eIdSf";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,2);

  //Set Branches//

  //ULong64_t event;
  int run; //, lumi, bx;

  if (!ISMC){
    //T->SetBranchAddress("event", &event);
    T->SetBranchAddress("run", &run);
    //T->SetBranchAddress("lumi", &lumi);
    //T->SetBranchAddress("bx", &bx);
  }

 /* string triggers[nTriggers] = {
    "HLT_Mu45_eta2p1_v",
    "HLT_Mu50_v",
    "HLT_TkMu50_v",
    "HLT_Mu30_TkMu11_v",
    "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v",
    "HLT_Ele105_CaloIdVT_GsfTrkIdT_v",
    "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v"
  };*/

  vector<bool> *trig_passed = 0;
  //vector<int> *trig_prescale = 0;
  //vector<string> *trig_name = 0;

  T->SetBranchAddress("trig_passed", &trig_passed);
  //T->SetBranchAddress("trig_prescale", &trig_prescale);
  //T->SetBranchAddress("trig_name", &trig_name);

  int nGen=MAXGEN, gen_PID[nGen], gen_index[nGen], gen_mother0[nGen], gen_mother1[nGen];
  float gen_pt[nGen], gen_mass[nGen], gen_eta[nGen], gen_phi[nGen];

  if (inName.EqualTo("TTbar.root")) {
    T->SetBranchAddress("nGen", &nGen);
    T->SetBranchAddress("gen_PID", gen_PID);
    T->SetBranchAddress("gen_pt", gen_pt);
    T->SetBranchAddress("gen_mass", gen_mass);
    T->SetBranchAddress("gen_eta", gen_eta);
    T->SetBranchAddress("gen_phi", gen_phi);
    T->SetBranchAddress("gen_index",  gen_index);
    T->SetBranchAddress("gen_mother0", gen_mother0);
    T->SetBranchAddress("gen_mother1", gen_mother1);
  }

  char lep0flavor, lep1flavor;
  T->SetBranchAddress("lep0flavor", &lep0flavor);
  T->SetBranchAddress("lep1flavor", &lep1flavor);

  int nMuon=MAXLEP, muon_charge[nMuon];
  float muon_eta[nMuon], muon_pt[nMuon], muon_phi[nMuon], muon_D0[nMuon], muon_Dz[nMuon];
  float muon_chi2[nMuon], muon_tspm[nMuon], muon_kinkf[nMuon], muon_segcom[nMuon], muon_ftrackhits[nMuon];
  bool muon_isGlob[nMuon], muon_IsMediumID[nMuon];

  T->SetBranchAddress("nMuon", &nMuon);
  T->SetBranchAddress("muon_charge", muon_charge);
  T->SetBranchAddress("muon_eta", muon_eta);
  T->SetBranchAddress("muon_pt", muon_pt);
  T->SetBranchAddress("muon_phi", muon_phi);
  T->SetBranchAddress("muon_D0", muon_D0);
  T->SetBranchAddress("muon_Dz", muon_Dz);
  T->SetBranchAddress("muon_IsMediumID", muon_IsMediumID);

  T->SetBranchAddress("muon_isGlob", muon_isGlob);
  T->SetBranchAddress("muon_chi2", muon_chi2);
  T->SetBranchAddress("muon_tspm", muon_tspm);
  T->SetBranchAddress("muon_kinkf", muon_kinkf);
  T->SetBranchAddress("muon_segcom", muon_segcom);
  T->SetBranchAddress("muon_ftrackhits", muon_ftrackhits);

  int nEle=MAXLEP, ele_charge[nEle];
  float ele_eta[nEle], ele_pt[nEle], ele_phi[nEle], ele_etaSupClust[nEle];
  bool ele_MediumID[nEle];

  T->SetBranchAddress("nEle", &nEle);
  T->SetBranchAddress("ele_charge", ele_charge);
  T->SetBranchAddress("ele_eta", ele_eta);
  T->SetBranchAddress("ele_pt", ele_pt);
  T->SetBranchAddress("ele_phi", ele_phi);
  T->SetBranchAddress("ele_etaSupClust", ele_etaSupClust);
  T->SetBranchAddress("ele_MediumID", ele_MediumID);

  int nJet=MAXJET;
  float jet_eta[nJet], jet_phi[nJet], jet_pt[nJet], jet_mass[nJet], jet_area[nJet];
  float jet_btag[nJet], jet_nhf[nJet], jet_nef[nJet], jet_chf[nJet], jet_muf[nJet], jet_elef[nJet], jet_numneutral[nJet], jet_chmult[nJet];
  char jet_clean[nJet];

  T->SetBranchAddress("nJet", &nJet);
  T->SetBranchAddress("jet_eta", jet_eta);
  T->SetBranchAddress("jet_phi", jet_phi);
  T->SetBranchAddress("jet_pt", jet_pt);
  T->SetBranchAddress("jet_mass", jet_mass);
  T->SetBranchAddress("jet_area", jet_area);
  T->SetBranchAddress("jet_clean", jet_clean);

  T->SetBranchAddress("jet_btag", jet_btag);
  T->SetBranchAddress("jet_nhf", jet_nhf);
  T->SetBranchAddress("jet_nef", jet_nef);
  T->SetBranchAddress("jet_chf", jet_chf);
  T->SetBranchAddress("jet_muf", jet_muf);
  T->SetBranchAddress("jet_elef", jet_elef);
  T->SetBranchAddress("jet_numneutral", jet_numneutral);
  T->SetBranchAddress("jet_chmult", jet_chmult);

  float met_pt, met_px, met_py, met_phi;
  T->SetBranchAddress("met_pt", &met_pt);
  T->SetBranchAddress("met_px", &met_px);
  T->SetBranchAddress("met_py", &met_py);
  T->SetBranchAddress("met_phi", &met_phi);

  float rho;
  T->SetBranchAddress("rho", &rho);

  //Loop Over Entries//
  int sameRlepjet=0;
  time_t start = time(NULL);

  for (Long64_t n=0; n<nEntries; n++){
    T->GetEntry(n);

    TLorentzVector lep0, lep1;
    weight = weight0;

    if (channel == "mm") {
      if (lep0flavor == 'm' && lep1flavor == 'm'){
        v_cuts[channelCut].second += weight;

        //HLT_Mu50 or HLT_TkMu50 triggers
        if ( !(*trig_passed)[1] && !(*trig_passed)[2] ) continue;

        if (ISMC) {
          double muTrackSf = muTrackSfGraph->Eval(muon_eta[0]) * muTrackSfGraph->Eval(muon_eta[1]);
          weight *= muTrackSf;

          double trigSf0=0, trigSf1=0;
          if (muon_pt[0] > 53) trigSf0 = muTrigSfHist->GetBinContent( muTrigSfHist->FindBin( fabs(muon_eta[0]), muon_pt[0]>muTrig_pT?muTrig_pT:muon_pt[0] ) );
          if (muon_pt[1] > 53) trigSf1 = muTrigSfHist->GetBinContent( muTrigSfHist->FindBin( fabs(muon_eta[1]), muon_pt[1]>muTrig_pT?muTrig_pT:muon_pt[1] ) );

          double muTrigSf = ( 1. - (1. - (trigSf0>1. ? 1. : trigSf0) )*(1. - (trigSf1>1. ? 1. : trigSf1) ) );
          weight *= muTrigSf;

          FillHist1D("muTrigSf", muTrigSf, 1.);
          FillHist1D("muTrackSf", muTrackSf, 1.);
        }
        v_cuts[trigCut].second += weight;

        if ( ISMC || inName.Contains("GH", TString::kIgnoreCase) ) {
          if ( !muon_IsMediumID[0] || !muon_IsMediumID[1] ) continue;
        }
        else {
          if ( !isMediumMuonBCDEF(muon_isGlob[0], muon_chi2[0], muon_tspm[0], muon_kinkf[0], muon_segcom[0], muon_ftrackhits[0]) ||
               !isMediumMuonBCDEF(muon_isGlob[1], muon_chi2[1], muon_tspm[1], muon_kinkf[1], muon_segcom[1], muon_ftrackhits[1]) ) continue;
        }
        if ( muon_pt[0] < 53 || muon_pt[1] < 25 ) continue;
        if (fabs(muon_eta[0]) > 2.4 || fabs(muon_eta[1]) > 2.4) continue;

        if (ISMC) {
          double muIdSf = muIdSfHist->GetBinContent( muIdSfHist->FindBin( fabs(muon_eta[0]), muon_pt[0]>muId_pT?muId_pT:muon_pt[0] ) )
                        * muIdSfHist->GetBinContent( muIdSfHist->FindBin( fabs(muon_eta[1]), muon_pt[1]>muId_pT?muId_pT:muon_pt[1] ) );

          weight *= muIdSf;
          FillHist1D("muIdSf", muIdSf, 1.);
        }
        v_cuts[lepkinCut].second += weight;

        if (muon_charge[0]*muon_charge[1] > 0) continue;
        v_cuts[signCut].second += weight;

        //use these events for em channel
        if ( nEle>0 && ele_MediumID[0] && ele_pt[0] > 25 ) continue;
        v_cuts[thirdLepCut].second += weight;

        lep0.SetPtEtaPhiM(muon_pt[0], muon_eta[0], muon_phi[0], MUONMASS);
        lep1.SetPtEtaPhiM(muon_pt[1], muon_eta[1], muon_phi[1], MUONMASS);

        if ((lep0+lep1).M() < 20) continue;
      }
      else continue;
    }
    else if (channel == "ee"){
      if (lep0flavor == 'e' && lep1flavor == 'e'){
        v_cuts[channelCut].second += weight;

        if (fabs(ele_eta[0]) > 2.5 || fabs(ele_eta[1]) > 2.5) continue;
        v_cuts[lepkinCut].second += weight;

        if (ele_charge[0]*ele_charge[1] > 0) continue;
        v_cuts[signCut].second += weight;

        lep0.SetPtEtaPhiM(ele_pt[0], ele_eta[0], ele_phi[0], ELEMASS);
        lep1.SetPtEtaPhiM(ele_pt[1], ele_eta[1], ele_phi[1], ELEMASS);

        if ((lep0+lep1).M() < 20) continue;
      }
      else continue;
    }
    else{
      if (lep0flavor != lep1flavor) {
        v_cuts[channelCut].second += weight;

        //HLT_Mu50 or HLT_TkMu50 triggers
        if ( !(*trig_passed)[1] && !(*trig_passed)[2] ) continue;

        if (ISMC) {
          double muTrackSf = muTrackSfGraph->Eval(muon_eta[0]);
          weight *= muTrackSf;

          double muTrigSf = 0;
          if (muon_pt[0] > 53) muTrigSf = muTrigSfHist->GetBinContent( muTrigSfHist->FindBin( fabs(muon_eta[0]), muon_pt[0]>muTrig_pT?muTrig_pT:muon_pt[0] ) );
          weight *= muTrigSf;

          FillHist1D("muTrigSf", muTrigSf, 1.);
          FillHist1D("muTrackSf", muTrackSf, 1.);
        }
        v_cuts[trigCut].second += weight;

        if ( !ele_MediumID[0] ) continue;

        if ( ISMC || inName.Contains("GH", TString::kIgnoreCase) ) {
          if ( !muon_IsMediumID[0] ) continue;
        }
        else {
          if ( !isMediumMuonBCDEF(muon_isGlob[0], muon_chi2[0], muon_tspm[0], muon_kinkf[0], muon_segcom[0], muon_ftrackhits[0]) ) continue;
        }
        if ( muon_pt[0] < 53 || ele_pt[0] < 25 ) continue;
        if (fabs(muon_eta[0]) > 2.4 || fabs(ele_eta[0]) > 2.5) continue;

        if (ISMC) {
          double muIdSf = muIdSfHist->GetBinContent( muIdSfHist->FindBin( fabs(muon_eta[0]), muon_pt[0]>muId_pT?muId_pT:muon_pt[0] ) );
          double eRecoSf = eRecoSfHist->GetBinContent( eRecoSfHist->FindBin( ele_etaSupClust[0], ele_pt[0]>eReco_pT?eReco_pT:ele_pt[0] ) );
          double eIdSf = eIdSfHist->GetBinContent( eIdSfHist->FindBin( ele_etaSupClust[0], ele_pt[0]>eId_pT?eId_pT:ele_pt[0] ) );

          weight *= muIdSf * eRecoSf * eIdSf;

          FillHist1D("muIdSf", muIdSf, 1.);
          FillHist1D("eRecoSf", eRecoSf, 1.);
          FillHist1D("eIdSf", eIdSf, 1.);
        }

        v_cuts[lepkinCut].second += weight;

        if (ele_charge[0]*muon_charge[0] > 0) continue;
        v_cuts[signCut].second += weight;

        v_cuts[thirdLepCut].second += weight;

        if (lep0flavor=='e') {
          lep0.SetPtEtaPhiM(ele_pt[0], ele_eta[0], ele_phi[0], ELEMASS);
          lep1.SetPtEtaPhiM(muon_pt[0], muon_eta[0], muon_phi[0], MUONMASS);
        }
        else {
          lep0.SetPtEtaPhiM(muon_pt[0], muon_eta[0], muon_phi[0], MUONMASS);
          lep1.SetPtEtaPhiM(ele_pt[0], ele_eta[0], ele_phi[0], ELEMASS);
        }
      }
      else continue;
    }
    v_cuts[dilepmassCut].second += weight;
    double dilepmass = (lep0+lep1).M();

    if ( lep0flavor==lep1flavor && (76<dilepmass && dilepmass<106) ) continue;
    v_cuts[dilepVetoCut].second += weight;

    string era = eras[0];
    for( map<string, pair<int, int> >::const_iterator it = m_IOV.begin(); it != m_IOV.end(); ++it ) {
      const pair<int, int>& interval = it->second;
      if (interval.first <= run && run <= interval.second) { era = it->first; break; }
    }

    vector<pair<int, double> > jet_index_corrpt;
    TLorentzVector minjet0, minjet1;
    double rmin0=99, rmin1=99;
    double ctype1_x=0, ctype1_y=0;
    double rl0cleanj=-1, rl1cleanj=-1, cleanjet0pt=-1, cleanjet1pt=-1;

    int nGoodJet=0;
    for (int i=0; i<nJet; i++){

      //loose jet cut
      if (fabs(jet_eta[i]) <= 2.7) {
        if (jet_nhf[i]>=0.99 || jet_nef[i]>=0.99 || (jet_numneutral[i]+jet_chmult[i])<=1) continue;
        if (fabs(jet_eta[i]) <= 2.4 && ( jet_chf[i]<=0 || jet_chmult[i]<=0 || jet_elef[i]>=0.99 )) continue;
      }
      else if (2.7 < fabs(jet_eta[i]) && fabs(jet_eta[i]) <= 3.0) {
        if (jet_nhf[i]>=0.98 || jet_nef[i]<=0.01 || jet_numneutral[i]<=2) continue;
      }
      else {
        if (jet_nef[i]>=0.9 || jet_numneutral[i]<=10) continue;
      }

      jetCorrectors[era]->setJetEta( jet_eta[i] );
      jetCorrectors[era]->setJetPt( jet_pt[i] );
      jetCorrectors[era]->setJetA( jet_area[i] );
      jetCorrectors[era]->setRho(rho);
      double corr_pt = jetCorrectors[era]->getCorrection() * jet_pt[i];
      jet_index_corrpt.push_back( make_pair(i, corr_pt) );

      TLorentzVector jet;
      jet.SetPtEtaPhiM(corr_pt, jet_eta[i], jet_phi[i], jet_mass[i]);

      if (corr_pt>30 && fabs(jet_eta[i])<2.5) {
        nGoodJet++;

        if (lep0.DeltaR(jet) < rmin0) {
          rmin0 = lep0.DeltaR(jet);
          minjet0 = jet;
        }
        if (lep1.DeltaR(jet) < rmin1) {
          rmin1 = lep1.DeltaR(jet);
          minjet1 = jet;
        }
      }

      if (jet_clean[i] == 'l' || jet_clean[i] == 'b') { rl0cleanj = lep0.DeltaR(jet); cleanjet0pt = corr_pt; }
      if (jet_clean[i] == 's' || jet_clean[i] == 'b') { rl1cleanj = lep1.DeltaR(jet); cleanjet1pt = corr_pt; }

      //corrected MET
      if ( corr_pt>15 && (jet_elef[i]+jet_nef[i])<0.9 ) {
        jetL1Correctors[era]->setJetEta( jet_eta[i] );
        jetL1Correctors[era]->setJetPt( jet_pt[i] );
        jetL1Correctors[era]->setJetA( jet_area[i] );
        jetL1Correctors[era]->setRho(rho);

        TLorentzVector jetL1;
        jetL1.SetPtEtaPhiM( jetL1Correctors[era]->getCorrection()*jet_pt[i], jet_eta[i], jet_phi[i], jet_mass[i] );

        ctype1_x += (jet.Px()-jetL1.Px());
        ctype1_y += (jet.Py()-jetL1.Py());
      }
    }

    if (minjet0 == minjet1) sameRlepjet++;
    double lep0perp = lep0.Perp( minjet0.Vect() );
    double lep1perp = lep1.Perp( minjet1.Vect() );

    if ( (lep0perp<10 && rmin0<0.4) || (lep1perp<10 && rmin1<0.4) ) continue;
    v_cuts[ptrelCut].second += weight;

    if (nGoodJet < 2) continue;
    sort(jet_index_corrpt.begin(), jet_index_corrpt.end(), sortJetPt);

    int jet0index = jet_index_corrpt[0].first, jet1index = jet_index_corrpt[1].first;
    double jet0pt = jet_index_corrpt[0].second, jet1pt = jet_index_corrpt[1].second;

   //at least one jet 
   if ( jet0pt < 100 || fabs(jet_eta[jet0index]) > 2.5 ) continue;

    double minjet0pt = minjet0.Pt();
    double minjet1pt = minjet1.Pt();
    double masslmin0 = (lep0+minjet0).M();
    double masslmin1 = (lep1+minjet1).M();

    double deta_lep = lep0.Eta() - lep1.Eta();
    double deta_lepJet = (lep0+minjet0).Eta() - (lep1+minjet1).Eta();

    double met_corrpx = met_px - ctype1_x;
    double met_corrpy = met_py - ctype1_y;
    double met_corrpt = sqrt(met_corrpx*met_corrpx + met_corrpy*met_corrpy);

    TLorentzVector met;
    met.SetPtEtaPhiE(met_corrpt, 0, met_phi, met_corrpt);
    TLorentzVector jet0, jet1;
    jet0.SetPtEtaPhiM(jet0pt, jet_eta[jet0index], jet_phi[jet0index], jet_mass[jet0index]);
    jet1.SetPtEtaPhiM(jet1pt, jet_eta[jet1index], jet_phi[jet1index], jet_mass[jet1index]);

    double masslljjm = (lep0+lep1+jet0+jet1+met).M();

    int nGoodMuon=0;
    for (int i=0; i<nMuon; i++) {
      if ( ISMC || inName.Contains("GH", TString::kIgnoreCase) ) {
        if (muon_IsMediumID[i]) nGoodMuon++;
      }
      else {
        if (isMediumMuonBCDEF(muon_isGlob[i], muon_chi2[i], muon_tspm[i], muon_kinkf[i], muon_segcom[i], muon_ftrackhits[i])) nGoodMuon++;
      }
    }

    int nGoodEle=0;
    for (int i=0; i<nEle; i++) {
      if (ele_MediumID[i]) nGoodEle++;
    }

    double hT=0;
    for (int i=0; i<nGoodJet; i++) {
      int index = jet_index_corrpt[i].first;
      double corr_pt = jet_index_corrpt[i].second;

      if ( corr_pt>30 && fabs(jet_eta[index])<2.5 ) hT+=corr_pt;
    }
    double sT = hT+lep0.Pt()+lep1.Pt();
    double sT_met = sT + met_corrpt;

    bool jet0btag = jet_btag[jet0index] > 0.8484 && fabs(jet_eta[jet0index]) < 2.4;
    bool jet1btag = jet_btag[jet1index] > 0.8484 && fabs(jet_eta[jet1index]) < 2.4;

    double rl0l1 = lep0.DeltaR(lep1);
    double lepept=0, lepmpt=0;
    if (lep0flavor == 'm') lepmpt += lep0.Pt();
    else lepept += lep0.Pt();
    if (lep1flavor == 'm') lepmpt += lep1.Pt();
    else lepept += lep1.Pt();

    double rbal=-1, rabl=-1;
    if (inName.EqualTo("TTbar.root")) {
      //lepton, anti-lepton, b quark, anti-b quark
      TLorentzVector glep, galep, gb, gab;

      for (int i=0; i<nGen; i++) {

        if (gen_PID[i] == 5) gb.SetPtEtaPhiM(gen_pt[i], gen_eta[i], gen_phi[i], gen_mass[i]);
        else if (gen_PID[i] == -5) gab.SetPtEtaPhiM(gen_pt[i], gen_eta[i], gen_phi[i], gen_mass[i]);
        else if (gen_PID[i] == 11 || gen_PID[i] == 13) glep.SetPtEtaPhiM(gen_pt[i], gen_eta[i], gen_phi[i], gen_mass[i]);
        else if (gen_PID[i] == -11 || gen_PID[i] == -13) galep.SetPtEtaPhiM(gen_pt[i], gen_eta[i], gen_phi[i], gen_mass[i]);
      }

      //b goes with anti-lepton and anti-b goes with lepton
      if (gb.Pt() != 0 && galep.Pt() != 0) rbal = gb.DeltaR(galep);
      if (gab.Pt() != 0 && glep.Pt() != 0) rabl = gab.DeltaR(glep);
    }

    //exactly one jet
    if ( jet1pt < 50 || fabs(jet_eta[jet1index]) > 2.5 ) {

      //zero btags
      if (!jet0btag) {
        v_cuts[jetCut1].second += weight;

        FillHists("0_", nEle, nGoodEle, nMuon, nGoodMuon, nJet, nGoodJet, lep0, lep1, dilepmass, lepept, lepmpt,
                  rmin0, rmin1, rl0l1, rl0cleanj, rl1cleanj, lep0perp, lep1perp, jet0, jet1, jet_btag[jet0index], jet_btag[jet1index], int(jet0btag)+int(jet1btag),
                  hT, met_pt, met_corrpt, sT, sT_met, rbal, rabl, minjet0pt, minjet1pt, cleanjet0pt, cleanjet1pt, masslmin0, masslmin1, masslljjm,
                  deta_lep, deta_lepJet);
      }
      //one btag
      else {
        v_cuts[onebtagCut1jet].second += weight;

        FillHists("1_", nEle, nGoodEle, nMuon, nGoodMuon, nJet, nGoodJet, lep0, lep1, dilepmass, lepept, lepmpt,
                  rmin0, rmin1, rl0l1, rl0cleanj, rl1cleanj, lep0perp, lep1perp, jet0, jet1, jet_btag[jet0index], jet_btag[jet1index], int(jet0btag)+int(jet1btag),
                  hT, met_pt, met_corrpt, sT, sT_met, rbal, rabl, minjet0pt, minjet1pt, cleanjet0pt, cleanjet1pt, masslmin0, masslmin1, masslljjm,
                  deta_lep, deta_lepJet);
      }
    }
    //two jets
    else {

      //exactly zero btags
      if ( !jet0btag && !jet1btag ) {
        v_cuts[jetCut2].second += weight;

        FillHists("2_", nEle, nGoodEle, nMuon, nGoodMuon, nJet, nGoodJet, lep0, lep1, dilepmass, lepept, lepmpt,
                  rmin0, rmin1, rl0l1, rl0cleanj, rl1cleanj, lep0perp, lep1perp, jet0, jet1, jet_btag[jet0index], jet_btag[jet1index], int(jet0btag)+int(jet1btag),
                  hT, met_pt, met_corrpt, sT, sT_met, rbal, rabl, minjet0pt, minjet1pt, cleanjet0pt, cleanjet1pt, masslmin0, masslmin1, masslljjm,
                  deta_lep, deta_lepJet);
      }
      //exactly two btags
      else if ( jet0btag && jet1btag ) {
        v_cuts[twobtagsCut2jets].second += weight;

        FillHists("4_", nEle, nGoodEle, nMuon, nGoodMuon, nJet, nGoodJet, lep0, lep1, dilepmass, lepept, lepmpt,
                  rmin0, rmin1, rl0l1, rl0cleanj, rl1cleanj, lep0perp, lep1perp, jet0, jet1, jet_btag[jet0index], jet_btag[jet1index], int(jet0btag)+int(jet1btag),
                  hT, met_pt, met_corrpt, sT, sT_met, rbal, rabl, minjet0pt, minjet1pt, cleanjet0pt, cleanjet1pt, masslmin0, masslmin1, masslljjm,
                  deta_lep, deta_lepJet);
      }
      //exactly one btag
      else {
        v_cuts[onebtagCut2jets].second += weight;

        FillHists("3_", nEle, nGoodEle, nMuon, nGoodMuon, nJet, nGoodJet, lep0, lep1, dilepmass, lepept, lepmpt,
                  rmin0, rmin1, rl0l1, rl0cleanj, rl1cleanj, lep0perp, lep1perp, jet0, jet1, jet_btag[jet0index], jet_btag[jet1index], int(jet0btag)+int(jet1btag),
                  hT, met_pt, met_corrpt, sT, sT_met, rbal, rabl, minjet0pt, minjet1pt, cleanjet0pt, cleanjet1pt, masslmin0, masslmin1, masslljjm,
                  deta_lep, deta_lepJet);
      }
    }
  }
  cout << difftime(time(NULL), start) << " s" << endl;
  cout << "Min_jet0 = Min_jet1: " << sameRlepjet << endl;

  TH1D* cuts = new TH1D("cuts","cuts",numCuts,-0.5,float(numCuts)-0.5);

  //Cutflow Table//
  cout<<"===================================================================================================\n";
  cout<<"                                     Cut Flow Table: " + inName + "\n";
  cout<<"===================================================================================================\n";

  cout<<      "                          |||          Nevent          |||     Efficiency (Relative Efficiency)\n";

  for (int i=0; i<numCuts; i++) {
    if (i == 0)
      cout << Form("%-25s |||       %12.1f       |||       %1.6f (%1.4f)",
                   v_cuts[i].first.data(), v_cuts[i].second, v_cuts[i].second/v_cuts[0].second, v_cuts[i].second/v_cuts[0].second) << endl;

    else if (i >= jetCut1)
      cout << Form("%-25s |||       %12.1f       |||       %1.6f (%1.4f)",
                   v_cuts[i].first.data(), v_cuts[i].second, v_cuts[i].second/v_cuts[0].second, v_cuts[i].second/v_cuts[jetCut1-1].second) << endl;

    else
      cout << Form("%-25s |||       %12.1f       |||       %1.6f (%1.4f)",
                   v_cuts[i].first.data(), v_cuts[i].second, v_cuts[i].second/v_cuts[0].second, v_cuts[i].second/v_cuts[i-1].second) << endl;

    if (i==countMet)
      cout << "---------------------------------------------------------------------------------------------------" << endl;

    cuts->SetBinContent(i+1, v_cuts[i].second);
    cuts->GetXaxis()->SetBinLabel(i+1,Form("%i",i+1));
  }
  cout << endl;

  //Write Histograms//

  TFile* outFile = new TFile(outName,"RECREATE");
  outFile->cd();

  cuts->Write();

  for (int i=0; i<nDirs; i++) outFile->mkdir( Form("%i/", i) );

  for (map<TString, TH1*>::iterator hid = m_Histos1D.begin(); hid != m_Histos1D.end(); hid++) {
    outFile->cd();
    TString prefix = hid->first(0, 1);

    if ( prefix.IsDigit() ) outFile->cd(outName + ":/" + prefix);

    hid->second->Write();
  }

  outFile->Write();
  delete outFile;
  outFile = 0;
}

void FillHist1D(const TString& histName, const Double_t& value, const double& weight) {
  map<TString, TH1*>::iterator hid=m_Histos1D.find(histName);
  if (hid==m_Histos1D.end())
    cout << "%FillHist1D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value, weight);
}

void FillHists(const TString& prefix, const int& nEle, const int& nGoodEle, const int& nMuon, const int& nGoodMuon, const int& nJet, const int& nGoodJet,
               const TLorentzVector& lep0, const TLorentzVector& lep1, const float& dilepmass, const float& lepept, const float& lepmpt,
               const float& rmin0, const float& rmin1, const float& rl0l1, const float& rl0cleanj, const float& rl1cleanj, const float& lep0perp, const float& lep1perp,
               const TLorentzVector& jet0, const TLorentzVector& jet1, const float& jet0btag, const float& jet1btag, const int& nbtag,
               const float& hT, const float& met_pt, const float& met_corrpt, const float& sT, const float& sT_met,
               const float& rbal, const float& rabl, const float& minjet0pt, const float& minjet1pt, const float& cleanjet0pt, const float& cleanjet1pt,
               const float& masslmin0, const float& masslmin1, const float& masslljjm, const float& deta_lep, const float& deta_lepJet) {

    FillHist1D(prefix+"nEleDiff", nEle-nGoodEle, weight);
    FillHist1D(prefix+"nMuonDiff", nMuon-nGoodMuon, weight);
    FillHist1D(prefix+"nJetDiff", nJet-nGoodJet, weight);
    FillHist1D(prefix+"nEle", nEle, weight);
    FillHist1D(prefix+"nMuon", nMuon, weight);
    FillHist1D(prefix+"nJet", nJet, weight);
    FillHist1D(prefix+"nGoodEle", nGoodEle, weight);
    FillHist1D(prefix+"nGoodMuon", nGoodMuon, weight);
    FillHist1D(prefix+"nGoodJet", nGoodJet, weight);

    FillHist1D(prefix+"lep0pt", lep0.Pt(), weight);
    FillHist1D(prefix+"lep0eta", lep0.Eta(), weight);
    FillHist1D(prefix+"lep0phi", lep0.Phi(), weight);
    FillHist1D(prefix+"lep1pt", lep1.Pt(), weight);
    FillHist1D(prefix+"lep1eta", lep1.Eta(), weight);
    FillHist1D(prefix+"lep1phi", lep1.Phi(), weight);
    FillHist1D(prefix+"dilepmass", dilepmass, weight);
    FillHist1D(prefix+"lepept", lepept, weight);
    FillHist1D(prefix+"lepmpt", lepmpt, weight);

    FillHist1D(prefix+"rmin0", rmin0, weight);
    FillHist1D(prefix+"rmin1", rmin1, weight);
    FillHist1D(prefix+"rl0l1", rl0l1, weight);
    FillHist1D(prefix+"rl0cleanj", rl0cleanj, weight);
    FillHist1D(prefix+"rl1cleanj", rl1cleanj, weight);
    FillHist1D(prefix+"lep0perp", lep0perp, weight);
    FillHist1D(prefix+"lep1perp", lep1perp, weight);

    FillHist1D(prefix+"jet0pt", jet0.Pt(), weight);
    FillHist1D(prefix+"jet0eta", jet0.Eta(), weight);
    FillHist1D(prefix+"jet0phi", jet0.Phi(), weight);
    FillHist1D(prefix+"jet0btag", jet0btag, weight);
    FillHist1D(prefix+"jet1pt", jet1.Pt(), weight);
    FillHist1D(prefix+"jet1eta", jet1.Eta(), weight);
    FillHist1D(prefix+"jet1phi", jet1.Phi(), weight);
    FillHist1D(prefix+"jet1btag", jet1btag, weight);
    FillHist1D(prefix+"nbtag", nbtag, weight);
    FillHist1D(prefix+"jethT", hT, weight);

    FillHist1D(prefix+"metpt", met_pt, weight);
    FillHist1D(prefix+"metcorrpt", met_corrpt, weight);
    FillHist1D(prefix+"sT", sT, weight);
    FillHist1D(prefix+"sT_met", sT_met, weight);

    if (rbal != -1) FillHist1D(prefix+"rbl", rbal, weight);
    if (rabl != -1) FillHist1D(prefix+"rbl", rabl, weight);

    FillHist1D(prefix+"minjet0pt", minjet0pt, weight);
    FillHist1D(prefix+"minjet1pt", minjet1pt, weight);
    FillHist1D(prefix+"cleanjet0pt", cleanjet0pt, weight);
    FillHist1D(prefix+"cleanjet1pt", cleanjet1pt, weight);
    FillHist1D(prefix+"masslmin0", masslmin0, weight);
    FillHist1D(prefix+"masslmin1", masslmin1, weight);
    FillHist1D(prefix+"masslljjm", masslljjm, weight);

    FillHist1D(prefix+"deta_lep", deta_lep, weight);
    FillHist1D(prefix+"deta_lepJet", deta_lepJet, weight);
}

bool isMediumMuonBCDEF(const bool& isGlob, const float& chi2, const float& tspm, const float& kinkf, const float& segcom, const float& ftrackhits) {

  bool goodGlob = isGlob && chi2<3 && tspm<12 && kinkf<20; 
  
  return ( ftrackhits>0.49 && segcom>(goodGlob ? 0.303 : 0.451) );
}

void setWeight(const string& wFile) {

  ifstream file(wFile);
  string line;

  while (getline(file, line)){

    if (line.length() > 0) {
      while (line.at(0) == ' ') line.erase(0, 1);
    }

    int delim_pos = line.find(' ');
    if (delim_pos == -1) continue;

    TString dataset = line.substr(0, delim_pos).data();
    if ( inName.Contains(dataset, TString::kIgnoreCase) ) {

      while (line.at(line.length()-1) == ' ') line.erase(line.length()-1, line.length());

      //weight is found in the last column
      delim_pos = line.length()-1;
      while (line.at(delim_pos) != ' ') delim_pos--;

      line.erase(0, delim_pos+1);
      weight0 = stod(line);
      break;
    }
  }
  file.close();
}

void setPars(const string& parFile) {

  ifstream file(parFile);
  string line;

  while (getline(file, line)){

    if (line.length() > 0){
      while (line.at(0) == ' ') line.erase(0, 1);
    }

    int delim_pos = line.find(' ');
    if (delim_pos == -1) continue;

    string var = line.substr(0, delim_pos);
    line.erase(0, delim_pos + 1);

    while (line.at(0) == ' ') line.erase(0, 1);
    while (line.at(line.length()-1) == ' ') line.erase(line.length()-1, line.length());

    if (var == "ISMC"){
      if (line == "true") ISMC = true;
      else ISMC = false;
    }
    else if (var == "inName") inName = line.data();
    else if (var == "outName") outName = line.data();
    else if (var == "muTrigSfName") muTrigSfName = line.data();
    else if (var == "muIdSfName") muIdSfName = line.data();
    else if (var == "muTrackSfName") muTrackSfName = line.data();
    else if (var == "eRecoSfName") eRecoSfName = line.data();
    else if (var == "eIdSfName") eIdSfName = line.data();
    else if (var == "channel") channel = line;
    else if (var == "eras") {
      while ( (delim_pos = line.find(' ')) != -1) {
        eras.push_back( line.substr(0, delim_pos) );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      eras.push_back( line );
    }
    else if (var == "jet_type") jet_type = line;
  }
  file.close();
}

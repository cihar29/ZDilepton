//charles harrington and bahareh roozbahani
//execute as analyze pars.txt mc_weights.txt

#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TTree.h"
#include "TLorentzVector.h"

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

map<TString, TH1*> m_Histos1D;

//parameters- edit in pars.txt
bool ISMC;
TString inName, outName;
string channel, era, jet_type;
double weight;

const int MAXJET = 50;
const int MAXLEP = 20;
const float MUONMASS = 0.10566;
const float ELEMASS = 0.;

int main(int argc, char* argv[]){

  if (argc == 1) { cout << "Please provide a parameter file." << endl; return -1; }

  string parFile = argv[1];
  setPars(parFile);

  weight = -1;
  if (argc == 3)    { string wFile = argv[2]; setWeight(wFile); }
  if (weight == -1) { cout << "Weight set to 1" << endl; weight = 1.; }
  else                cout << "Weight set to " << weight << endl;

  //Open Files//

  TFile* inFile = TFile::Open(inName);

  TTree* T = (TTree*) inFile->Get("T");
  Long64_t nEntries = T->GetEntries();
  cout << nEntries << " Events" << endl;
  cout << "Processing " + inName << endl;

  //Skims//

  int countEvts=0, countDilep=0, countLeppt=0, countDilepmass=0, countJetpteta=0, countMet=0;

  TIter nextkey(inFile->GetListOfKeys());
  TKey* key;
  while ( (key = (TKey*)nextkey()) ) {
    TString keyname = key->GetName();

    if (keyname.EqualTo("totalEvts"))          countEvts += (*(vector<int>*)key->ReadObj())[0];
    else if (keyname.EqualTo("dilep_cut"))     countDilep += (*(vector<int>*)key->ReadObj())[0];
    else if (keyname.EqualTo("leppt_cut"))     countLeppt += (*(vector<int>*)key->ReadObj())[0];
    else if (keyname.EqualTo("dilepmass_cut")) countDilepmass += (*(vector<int>*)key->ReadObj())[0];
    else if (keyname.EqualTo("jetpteta_cut"))  countJetpteta += (*(vector<int>*)key->ReadObj())[0];
    else if (keyname.EqualTo("met_cut"))       countMet += (*(vector<int>*)key->ReadObj())[0];
    //else if (keyname.EqualTo("filter_failed"))
  }
  if (countMet != nEntries) { cout << "hadd added incorrectly." << endl; return -1; }

  //Histograms//

  int nDirs = 3;
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
    hname = Form("%i_jet0btag",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,1);
    hname = Form("%i_jet1btag",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,1);
    hname = Form("%i_nbtag",i);
    m_Histos1D[hname] = new TH1F(hname,hname,5,0,5);
    hname = Form("%i_jethT",i);
    m_Histos1D[hname] = new TH1F(hname,hname,200,0,2000);

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
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,500);
    hname = Form("%i_lep0eta",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,-5,5);
    hname = Form("%i_lep1pt",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,500);
    hname = Form("%i_lep1eta",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,-5,5);
    hname = Form("%i_dilepmass",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,500);

    hname = Form("%i_muonD0",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,-5,5);
    hname = Form("%i_muonDz",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,-10,10);
    hname = Form("%i_rmin0",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,5);
    hname = Form("%i_rmin1",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,5);
    hname = Form("%i_lep0perp",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,500);
    hname = Form("%i_lep1perp",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,500);

    hname = Form("%i_metpt",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,500);
    hname = Form("%i_metcorrpt",i);
    m_Histos1D[hname] = new TH1F(hname,hname,100,0,500);
    hname = Form("%i_sT",i);
    m_Histos1D[hname] = new TH1F(hname,hname,200,0,2000);
  }

  //Jet Corrections//

  string data = "DATA";
  if (ISMC) data = "MC";
  era += "_" + data;

  JetCorrectorParameters *ResJetPar = new JetCorrectorParameters(era + "/" + era + "_L2L3Residual_" + jet_type + ".txt"); 
  JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(era + "/" + era + "_L3Absolute_" + jet_type + ".txt");
  JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(era + "/" + era + "_L2Relative_" + jet_type + ".txt");
  JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(era + "/" + era + "_L1FastJet_" + jet_type + ".txt");

  vector<JetCorrectorParameters> jetPars;
  jetPars.push_back(*L1JetPar);
  jetPars.push_back(*L2JetPar);
  jetPars.push_back(*L3JetPar);
  jetPars.push_back(*ResJetPar);

  FactorizedJetCorrector *jetCorrector = new FactorizedJetCorrector(jetPars);

  vector<JetCorrectorParameters> jetL1Pars;
  jetL1Pars.push_back(*L1JetPar);
  FactorizedJetCorrector *jetL1Corrector = new FactorizedJetCorrector(jetL1Pars);

  //Set Branches//

  //ULong64_t event;
  //int run, lumi, bx;

  if (!ISMC){
    //T->SetBranchAddress("event", &event);
    //T->SetBranchAddress("run", &run);
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

  char lep0flavor, lep1flavor;
  T->SetBranchAddress("lep0flavor", &lep0flavor);
  T->SetBranchAddress("lep1flavor", &lep1flavor);

  int nMuon=MAXLEP, muon_charge[nMuon];
  float muon_eta[nMuon], muon_pt[nMuon], muon_phi[nMuon], muon_D0[nMuon], muon_Dz[nMuon];
  float muon_chi2[nMuon], muon_tspm[nMuon], muon_kinkf[nMuon], muon_segcom[nMuon], muon_ftrackhits[nMuon];
  bool muon_isGlob[nMuon], muon_IsMediumID[nMuon], muon_IsTightID[nMuon];

  T->SetBranchAddress("nMuon", &nMuon);
  T->SetBranchAddress("muon_charge", muon_charge);
  T->SetBranchAddress("muon_eta", muon_eta);
  T->SetBranchAddress("muon_pt", muon_pt);
  T->SetBranchAddress("muon_phi", muon_phi);
  T->SetBranchAddress("muon_D0", muon_D0);
  T->SetBranchAddress("muon_Dz", muon_Dz);
  T->SetBranchAddress("muon_IsMediumID", muon_IsMediumID);
  T->SetBranchAddress("muon_IsTightID", muon_IsTightID);

  T->SetBranchAddress("muon_isGlob", muon_isGlob);
  T->SetBranchAddress("muon_chi2", muon_chi2);
  T->SetBranchAddress("muon_tspm", muon_tspm);
  T->SetBranchAddress("muon_kinkf", muon_kinkf);
  T->SetBranchAddress("muon_segcom", muon_segcom);
  T->SetBranchAddress("muon_ftrackhits", muon_ftrackhits);

  int nEle=MAXLEP, ele_charge[nEle];
  float ele_eta[nEle], ele_pt[nEle], ele_phi[nEle];
  bool ele_MediumID[nEle];

  T->SetBranchAddress("nEle", &nEle);
  T->SetBranchAddress("ele_charge", ele_charge);
  T->SetBranchAddress("ele_eta", ele_eta);
  T->SetBranchAddress("ele_pt", ele_pt);
  T->SetBranchAddress("ele_phi", ele_phi);
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

  float met_pt, met_px, met_py;
  T->SetBranchAddress("met_pt", &met_pt);
  T->SetBranchAddress("met_px", &met_px);
  T->SetBranchAddress("met_py", &met_py);

  float rho;
  T->SetBranchAddress("rho", &rho);

  //Loop Over Entries//
  int channelCut=0, signCut=0, trigCut=0, lepkinCut=0, thirdLepCut=0, dilepmassCut=0, jetCut=0, btagCut=0 , DilepVetoCut=0;
  int sameRlepjet=0;
  time_t start = time(NULL);

  for (Long64_t n=0; n<nEntries; n++){
    T->GetEntry(n);

    TLorentzVector lep0, lep1;

   if (channel == "mm"){
      if (lep0flavor == 'm' && lep1flavor == 'm'){
        channelCut++;

        if (muon_charge[0]*muon_charge[1] > 0) continue;
        signCut++;

        //HLT_Mu50 or HLT_TkMu50 triggers
        //TString trig1 = (*trig_name)[1].data(), trig2 = (*trig_name)[2].data();
        //if ( !trig1.Contains("HLT_Mu50_v", TString::kIgnoreCase) || !trig2.Contains("HLT_TkMu50_v", TString::kIgnoreCase) ) {
        //  cout << "wrong trigger\t" << trig1 << "\t" << trig2 << endl;
        //  return -1;
        //}
        if ( !(*trig_passed)[1] && !(*trig_passed)[2] ) continue;
        trigCut++;

        //if ( ISMC || inName.Contains("GH", TString::kIgnoreCase) ) {
        //  if ( !muon_IsMediumID[0] || !muon_IsMediumID[1] ) continue;
        //}
        //else {
        //  if ( !isMediumMuonBCDEF(muon_isGlob[0], muon_chi2[0], muon_tspm[0], muon_kinkf[0], muon_segcom[0], muon_ftrackhits[0]) ||
        //       !isMediumMuonBCDEF(muon_isGlob[1], muon_chi2[1], muon_tspm[1], muon_kinkf[1], muon_segcom[1], muon_ftrackhits[1]) ) continue;
        //}

        if ( !muon_IsTightID[0] || !muon_IsTightID[1] ) continue;

        if ( muon_pt[0] < 53 || muon_pt[1] < 25 ) continue;
        if (fabs(muon_eta[0]) > 2.4 || fabs(muon_eta[1]) > 2.4) continue;
        lepkinCut++;

        //use these events for em channel
        if ( nEle>0 && ele_MediumID[0] && ele_pt[0] > 25 ) continue;
        thirdLepCut++;

        lep0.SetPtEtaPhiM(muon_pt[0], muon_eta[0], muon_phi[0], MUONMASS);
        lep1.SetPtEtaPhiM(muon_pt[1], muon_eta[1], muon_phi[1], MUONMASS);

        if ((lep0+lep1).M() < 20) continue;
      }
      else continue;
    }
    else if (channel == "ee"){
      if (lep0flavor == 'e' && lep1flavor == 'e'){
        channelCut++;

        if (ele_charge[0]*ele_charge[1] > 0) continue;
        signCut++;

        if (fabs(ele_eta[0]) > 2.5 || fabs(ele_eta[1]) > 2.5) continue;
        lepkinCut++;

        lep0.SetPtEtaPhiM(ele_pt[0], ele_eta[0], ele_phi[0], ELEMASS);
        lep1.SetPtEtaPhiM(ele_pt[1], ele_eta[1], ele_phi[1], ELEMASS);

        if ((lep0+lep1).M() < 20) continue;
      }
      else continue;
    }
    else{
      if (lep0flavor != lep1flavor) {
        channelCut++;

        if (ele_charge[0]*muon_charge[0] > 0) continue;
        signCut++;

        //HLT_Mu50 or HLT_TkMu50 triggers
        //TString trig1 = (*trig_name)[1].data(), trig2 = (*trig_name)[2].data();
        //if ( !trig1.Contains("HLT_Mu50_v", TString::kIgnoreCase) || !trig2.Contains("HLT_TkMu50_v", TString::kIgnoreCase) ) {
        //  cout << "wrong trigger\t" << trig1 << "\t" << trig2 << endl;
        //  return -1;
        //}
        if ( !(*trig_passed)[1] && !(*trig_passed)[2] ) continue;
        trigCut++;

        if ( !ele_MediumID[0] ) continue;

        //if ( ISMC || inName.Contains("GH", TString::kIgnoreCase) ) {
        //  if ( !muon_IsMediumID[0] ) continue;
        //}
        //else {
        //  //cout<< "custom made medium muons"<<endl;
        //  if ( !isMediumMuonBCDEF(muon_isGlob[0], muon_chi2[0], muon_tspm[0], muon_kinkf[0], muon_segcom[0], muon_ftrackhits[0]) ) continue;
        //}

        if ( !muon_IsTightID[0] ) continue;

        if ( muon_pt[0] < 53 || ele_pt[0] < 25 ) continue;
        if (fabs(muon_eta[0]) > 2.4 || fabs(ele_eta[0]) > 2.5) continue;
        lepkinCut++;
        thirdLepCut++;

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
    dilepmassCut++;
    double dilepmass = (lep0+lep1).M();

    vector<pair<int, float> > jet_index_corrpt;
    TLorentzVector minjet0, minjet1;
    float rmin0=99, rmin1=99;
    double ctype1_x=0, ctype1_y=0;

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

      jetCorrector->setJetEta( jet_eta[i] );
      jetCorrector->setJetPt( jet_pt[i] );
      jetCorrector->setJetA( jet_area[i] );
      jetCorrector->setRho(rho);
      double corr_pt = jetCorrector->getCorrection() * jet_pt[i];
      jet_index_corrpt.push_back( make_pair(i, corr_pt) );

      TLorentzVector jet;
      jet.SetPtEtaPhiM(corr_pt, jet_eta[i], jet_phi[i], jet_mass[i]);

      if (lep0.DeltaR(jet) < rmin0) {
        rmin0 = lep0.DeltaR(jet);
        minjet0 = jet;
      }
      if (lep1.DeltaR(jet) < rmin1) {
        rmin1 = lep1.DeltaR(jet);
        minjet1 = jet;
      }

      //corrected MET
      if ( corr_pt>15 && (jet_elef[i]+jet_nef[i])<0.9 ) {
        jetL1Corrector->setJetEta( jet_eta[i] );
        jetL1Corrector->setJetPt( jet_pt[i] );
        jetL1Corrector->setJetA( jet_area[i] );
        jetL1Corrector->setRho(rho);

        TLorentzVector jetL1;
        jetL1.SetPtEtaPhiM( jetL1Corrector->getCorrection()*jet_pt[i], jet_eta[i], jet_phi[i], jet_mass[i] );

        ctype1_x += (jet.Px()-jetL1.Px());
        ctype1_y += (jet.Py()-jetL1.Py());
      }
    }
    sort(jet_index_corrpt.begin(), jet_index_corrpt.end(), sortJetPt);

    int nGoodJet = jet_index_corrpt.size();
    if (nGoodJet < 2) continue;

    int jet0index = jet_index_corrpt[0].first, jet1index = jet_index_corrpt[1].first;
    float jet0pt = jet_index_corrpt[0].second, jet1pt = jet_index_corrpt[1].second;

    //if ( jet0pt < 100 || jet1pt < 50 ) continue;
    //if ( fabs(jet_eta[jet0index]) > 2.5 || fabs(jet_eta[jet1index]) > 2.5 ) continue;
    if ( jet0pt < 100 ) continue;
    if ( fabs(jet_eta[jet0index]) > 2.5 ) continue;
    jetCut++;

    if (minjet0 == minjet1) sameRlepjet++;
    float perp0 = lep0.Perp( minjet0.Vect() );
    float perp1 = lep1.Perp( minjet1.Vect() );

    double met_corrpx = met_px - ctype1_x;
    double met_corrpy = met_py - ctype1_y;
    double met_corrpt = sqrt(met_corrpx*met_corrpx + met_corrpy*met_corrpy);

    int nGoodMuon=0;
    for (int i=0; i<nMuon; i++) {
      if (muon_IsTightID[i]) nGoodMuon++;

      //if ( ISMC || inName.Contains("GH", TString::kIgnoreCase) ) {
      //  if (muon_IsMediumID[i]) nGoodMuon++;
      //}
      //else {
      //  if (isMediumMuonBCDEF(muon_isGlob[i], muon_chi2[i], muon_tspm[i], muon_kinkf[i], muon_segcom[i], muon_ftrackhits[i])) nGoodMuon++;
      //}
    }

    int nGoodEle=0;
    for (int i=0; i<nEle; i++) {
      if (ele_MediumID[i]) nGoodEle++;
    }

    float hT=0;
    for (int i=0; i<nGoodJet; i++) {
      int index = jet_index_corrpt[i].first;
      float corr_pt = jet_index_corrpt[i].second;

      if ( corr_pt>30 && fabs(jet_eta[index])<2.5 ) hT+=corr_pt;
    }
    float sT = hT+lep0.Pt()+lep1.Pt();

    bool jet0btag = jet_btag[jet0index] > 0.8484 && fabs(jet_eta[jet0index]) < 2.4;
    bool jet1btag = jet_btag[jet1index] > 0.8484 && fabs(jet_eta[jet1index]) < 2.4;

    TString prefix = "0_";
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
    FillHist1D(prefix+"lep1pt", lep1.Pt(), weight);
    FillHist1D(prefix+"lep1eta", lep1.Eta(), weight);
    FillHist1D(prefix+"dilepmass", dilepmass, weight);

    if (lep0flavor == 'm') {
      FillHist1D(prefix+"muonD0", muon_D0[0], weight);
      FillHist1D(prefix+"muonDz", muon_Dz[0], weight);
    }
    if (lep1flavor == 'm') {
      FillHist1D(prefix+"muonD0", muon_D0[1], weight);
      FillHist1D(prefix+"muonDz", muon_Dz[1], weight);
    }
    FillHist1D(prefix+"rmin0", rmin0, weight);
    FillHist1D(prefix+"rmin1", rmin1, weight);
    FillHist1D(prefix+"lep0perp", perp0, weight);
    FillHist1D(prefix+"lep1perp", perp1, weight);

    FillHist1D(prefix+"jet0pt", jet0pt, weight);
    FillHist1D(prefix+"jet0eta", jet_eta[jet0index], weight);
    FillHist1D(prefix+"jet0btag", jet_btag[jet0index], weight);
    FillHist1D(prefix+"jet1pt", jet1pt, weight);
    FillHist1D(prefix+"jet1eta", jet_eta[jet1index], weight);
    FillHist1D(prefix+"jet1btag", jet_btag[jet1index], weight);
    FillHist1D(prefix+"nbtag", int(jet0btag)+int(jet1btag), weight);
    FillHist1D(prefix+"jethT", hT, weight);

    FillHist1D(prefix+"metpt", met_pt, weight);
    FillHist1D(prefix+"metcorrpt", met_corrpt, weight);
    FillHist1D(prefix+"sT", sT, weight);

    if ( !jet0btag && !jet1btag ) continue;
    btagCut++;

    prefix = "1_";
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
    FillHist1D(prefix+"lep1pt", lep1.Pt(), weight);
    FillHist1D(prefix+"lep1eta", lep1.Eta(), weight);
    FillHist1D(prefix+"dilepmass", dilepmass, weight);

    if (lep0flavor == 'm') {
      FillHist1D(prefix+"muonD0", muon_D0[0], weight);
      FillHist1D(prefix+"muonDz", muon_Dz[0], weight);
    }
    if (lep1flavor == 'm') {
      FillHist1D(prefix+"muonD0", muon_D0[1], weight);
      FillHist1D(prefix+"muonDz", muon_Dz[1], weight);
    }
    FillHist1D(prefix+"rmin0", rmin0, weight);
    FillHist1D(prefix+"rmin1", rmin1, weight);
    FillHist1D(prefix+"lep0perp", perp0, weight);
    FillHist1D(prefix+"lep1perp", perp1, weight);

    FillHist1D(prefix+"jet0pt", jet0pt, weight);
    FillHist1D(prefix+"jet0eta", jet_eta[jet0index], weight);
    FillHist1D(prefix+"jet0btag", jet_btag[jet0index], weight);
    FillHist1D(prefix+"jet1pt", jet1pt, weight);
    FillHist1D(prefix+"jet1eta", jet_eta[jet1index], weight);
    FillHist1D(prefix+"jet1btag", jet_btag[jet1index], weight);
    FillHist1D(prefix+"nbtag", int(jet0btag)+int(jet1btag), weight);
    FillHist1D(prefix+"jethT", hT, weight);

    FillHist1D(prefix+"metpt", met_pt, weight);
    FillHist1D(prefix+"metcorrpt", met_corrpt, weight);
    FillHist1D(prefix+"sT", sT, weight);


    if (lep0flavor == lep1flavor) {
      if ( 76<dilepmass && dilepmass<106 ) continue;

      DilepVetoCut++;
      prefix = "2_";
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
      FillHist1D(prefix+"lep1pt", lep1.Pt(), weight);
      FillHist1D(prefix+"lep1eta", lep1.Eta(), weight);
      FillHist1D(prefix+"dilepmass", dilepmass, weight);

      if (lep0flavor == 'm') {
        FillHist1D(prefix+"muonD0", muon_D0[0], weight);
        FillHist1D(prefix+"muonDz", muon_Dz[0], weight);
      }
      if (lep1flavor == 'm') {
        FillHist1D(prefix+"muonD0", muon_D0[1], weight);
        FillHist1D(prefix+"muonDz", muon_Dz[1], weight);
      }
      FillHist1D(prefix+"rmin0", rmin0, weight);
      FillHist1D(prefix+"rmin1", rmin1, weight);
      FillHist1D(prefix+"lep0perp", perp0, weight);
      FillHist1D(prefix+"lep1perp", perp1, weight);

      FillHist1D(prefix+"jet0pt", jet0pt, weight);
      FillHist1D(prefix+"jet0eta", jet_eta[jet0index], weight);
      FillHist1D(prefix+"jet0btag", jet_btag[jet0index], weight);
      FillHist1D(prefix+"jet1pt", jet1pt, weight);
      FillHist1D(prefix+"jet1eta", jet_eta[jet1index], weight);
      FillHist1D(prefix+"jet1btag", jet_btag[jet1index], weight);
      FillHist1D(prefix+"nbtag", int(jet0btag)+int(jet1btag), weight);
      FillHist1D(prefix+"jethT", hT, weight);

      FillHist1D(prefix+"metpt", met_pt, weight);
      FillHist1D(prefix+"metcorrpt", met_corrpt, weight);
      FillHist1D(prefix+"sT", sT, weight);
    }
  }
  cout << difftime(time(NULL), start) << " s" << endl;
  cout << "Min_jet0 = Min_jet1: " << sameRlepjet << endl;

  //Cutflow Table//

  cout<<"====================================================================================================================="<< "\n" ;
  cout<<"                                     Cut Flow Table: " + inName + "\n" ;
  cout<<"====================================================================================================================="<< "\n" ;

  cout<<      "                                  |||         Nevent        |||   Relative Efficiency    |||     Efficiency      " << "\n" ;
  cout<< Form("        Initial                   |||         %10i          |||           %1.3f          |||       %1.6f         ",countEvts,float(countEvts)/countEvts,float(countEvts)/countEvts) << "\n";
  cout<< Form("  Passed Dilepton selection       |||         %10i          |||           %1.3f          |||       %1.6f         ",countDilep,float(countDilep)/countEvts,float(countDilep)/countEvts) << "\n";
  cout<< Form("  Passed lepton Pt Cut            |||         %10i          |||           %1.3f          |||       %1.6f         ",countLeppt,float(countLeppt)/countDilep,float(countLeppt)/countEvts) << "\n";
  cout<< Form("  Passed Dilepton Mass Cut        |||         %10i          |||           %1.3f          |||       %1.6f         ",countDilepmass,float(countDilepmass)/countLeppt,float(countDilepmass)/countEvts) << "\n";
  cout<< Form("  Passed Leading Jet Pt_eta cut   |||         %10i          |||           %1.3f          |||       %1.6f         ",countJetpteta,float(countJetpteta)/countDilepmass,float(countJetpteta)/countEvts) << "\n";
  cout<< Form("  Passed MET Filters              |||         %10i          |||           %1.3f          |||       %1.6f         ",countMet,float(countMet)/countJetpteta,float(countMet)/countEvts) << "\n";
  cout<< Form("  Passed Correct Channel          |||         %10i          |||           %1.3f          |||       %1.6f         ",channelCut,float(channelCut)/countMet,float(channelCut)/countEvts) << "\n";
  cout<< Form("  Passed Opposite Lepton Sign     |||         %10i          |||           %1.3f          |||       %1.6f         ",signCut,float(signCut)/channelCut,float(signCut)/countEvts) << "\n";
  cout<< Form("  Passed HLT Trigger              |||         %10i          |||           %1.3f          |||       %1.6f         ",trigCut,float(trigCut)/signCut,float(trigCut)/countEvts) << "\n";
  cout<< Form("  Passed lepton kinematics cut    |||         %10i          |||           %1.3f          |||       %1.6f         ",lepkinCut,float(lepkinCut)/trigCut,float(lepkinCut)/countEvts) << "\n";
  cout<< Form("  Passed third Lepton Cut         |||         %10i          |||           %1.3f          |||       %1.6f         ",thirdLepCut,float(thirdLepCut)/lepkinCut,float(thirdLepCut)/countEvts) << "\n";
  cout<< Form("  Passed dilepton mass Cut        |||         %10i          |||           %1.3f          |||       %1.6f         ",dilepmassCut,float(dilepmassCut)/thirdLepCut,float(dilepmassCut)/countEvts) << "\n";
  cout<< Form("  Passed jet Cut                  |||         %10i          |||           %1.3f          |||       %1.6f         ",jetCut,float(jetCut)/dilepmassCut,float(jetCut)/countEvts) << "\n";
  cout<< Form("  Passed btag Cut                 |||         %10i          |||           %1.3f          |||       %1.6f         ",btagCut,float(btagCut)/jetCut,float(btagCut)/countEvts) << "\n";
  cout<< Form("  Passed Dilepmass window Cut     |||         %10i          |||           %1.3f          |||       %1.6f         ",DilepVetoCut,float(DilepVetoCut)/btagCut,float(DilepVetoCut)/countEvts) << "\n";

  //Write Histograms//

  TFile* outFile = new TFile(outName,"RECREATE");
  outFile->cd();

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
      weight = stod(line);
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
    if ((delim_pos = line.find(' ')) != -1) line.erase(delim_pos, line.length());

    if (var == "ISMC"){
      if (line == "true") ISMC = true;
      else ISMC = false;
    }
    else if (var == "inName") inName = line.data();
    else if (var == "outName") outName = line.data();
    else if (var == "channel") channel = line;
    else if (var == "era") era = line;
    else if (var == "jet_type") jet_type = line;
  }
  file.close();
}

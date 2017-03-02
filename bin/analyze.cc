//charles harrington and bahareh roozbahani
//execute as analyze pars.txt

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
bool sortJetPt(pair<int, float> jet1, pair<int, float> jet2){ return jet1.second > jet2.second; }

map<TString, TH1*> m_Histos1D;

//default parameters- edit in pars.txt
bool ISMC = false;
string inName = "analysis_Data.root";
string outName = "hist.root";
string channel = "mm";
string era = "Spring16_25nsV10BCD";
string jet_type = "AK4PFchs";
double weight = 1.;

const int MAXJET = 50;
const int MAXLEP = 20;
const float MUONMASS = 0.10566;
const float ELEMASS = 0;

int main(int argc, char* argv[]){

  string parFile = argv[1];
  setPars(parFile);

  //Open Files//

  TFile* inFile = TFile::Open(inName.data());

  TTree* T = (TTree*) inFile->Get("T");
  Long64_t nEntries = T->GetEntries();
  cout << nEntries << " Events" << endl;

  //Skims//

  int countEvts=0, countDilep=0, countLeppt=0, countDilepmass=0, countJetpteta=0, countMet=0;

  TIter nextkey(inFile->GetListOfKeys());
  TKey *key;
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
  if (countMet != nEntries) { cout << "hadd added incorrectly" << endl; return 0; }

  //Histograms//

  TString hname = "nJet";
  m_Histos1D[hname] = new TH1F(hname,hname,MAXJET,0,MAXJET);
  hname = "jet0pt";
  m_Histos1D[hname] = new TH1F(hname,hname,200,0,2000);
  hname = "jet1pt";
  m_Histos1D[hname] = new TH1F(hname,hname,200,0,2000);
  hname = "jet0eta";
  m_Histos1D[hname] = new TH1F(hname,hname,100,-5,5);
  hname = "jet1eta";
  m_Histos1D[hname] = new TH1F(hname,hname,100,-5,5);

  hname = "nEle";
  m_Histos1D[hname] = new TH1F(hname,hname,MAXLEP,0,MAXLEP);
  hname = "nMuon";
  m_Histos1D[hname] = new TH1F(hname,hname,MAXLEP,0,MAXLEP);
  hname = "lep0pt";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,500);
  hname = "lep0eta";
  m_Histos1D[hname] = new TH1F(hname,hname,100,-5,5);
  hname = "lep1pt";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,500);
  hname = "lep1eta";
  m_Histos1D[hname] = new TH1F(hname,hname,100,-5,5);
  hname = "dilepmass";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,500);

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
  float muon_eta[nMuon], muon_pt[nMuon], muon_phi[nMuon];
  bool muon_IsMediumID[nMuon];

  T->SetBranchAddress("nMuon", &nMuon);
  T->SetBranchAddress("muon_charge", muon_charge);
  T->SetBranchAddress("muon_eta", muon_eta);
  T->SetBranchAddress("muon_pt", muon_pt);
  T->SetBranchAddress("muon_phi", muon_phi);
  T->SetBranchAddress("muon_IsMediumID", muon_IsMediumID);

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
  float jet_eta[nJet], jet_pt[nJet], jet_area[nJet];

  T->SetBranchAddress("nJet", &nJet);
  T->SetBranchAddress("jet_eta", jet_eta);
  T->SetBranchAddress("jet_pt", jet_pt);
  T->SetBranchAddress("jet_area", jet_area);

  float rho;
  T->SetBranchAddress("rho", &rho);

  //Loop Over Entries//
  int channelCut=0, signCut=0, trigCut=0, lepptCut=0, thirdLepCut=0, etaCut=0, dilepmassCut=0, jetCut=0;
  time_t start = time(NULL);

  for (Long64_t n=0; n<nEntries; n++){
    if (n % 10000 == 0)
      cout << "Processing Event " << n+1 << endl;
    T->GetEntry(n);

    TLorentzVector lep0, lep1;

   if (channel == "mm"){
      if (lep0flavor == 'm' && lep1flavor == 'm'){
        channelCut++;

        if (muon_charge[0]*muon_charge[1] > 0) continue;
        signCut++;

        //HLT_Mu50 or HLT_TkMu50 triggers
        if ( !(*trig_passed)[1] && !(*trig_passed)[2] ) continue;
        trigCut++;

        if ( !muon_IsMediumID[0] || !muon_IsMediumID[1] ) continue;
        if ( muon_pt[0] < 52 || muon_pt[1] < 25 ) continue;
        lepptCut++;

        //use these events for em channel
        if ( nEle>0 && ele_MediumID[0] && ele_pt[0] > 25 ) continue;
        thirdLepCut++;

        if (fabs(muon_eta[0]) > 2.4 || fabs(muon_eta[1]) > 2.4) continue;
        etaCut++;

        lep0.SetPtEtaPhiM(muon_pt[0], muon_eta[0], muon_phi[0], MUONMASS);
        lep1.SetPtEtaPhiM(muon_pt[1], muon_eta[1], muon_phi[1], MUONMASS);
      }
      else continue;
    }
    else if (channel == "ee"){
      if (lep0flavor == 'e' && lep1flavor == 'e'){
        channelCut++;

        if (ele_charge[0]*ele_charge[1] > 0) continue;
        signCut++;

        if (fabs(ele_eta[0]) > 2.5 || fabs(ele_eta[1]) > 2.5) continue;
        etaCut++;

        lep0.SetPtEtaPhiM(ele_pt[0], ele_eta[0], ele_phi[0], ELEMASS);
        lep1.SetPtEtaPhiM(ele_pt[1], ele_eta[1], ele_phi[1], ELEMASS);
      }
      else continue;
    }
    else{
      if (lep0flavor == 'e' && lep1flavor == 'm'){
        channelCut++;

        if (ele_charge[0]*muon_charge[0] > 0) continue;
        signCut++;

        if (fabs(ele_eta[0]) > 2.5 || fabs(muon_eta[0]) > 2.4) continue;
        etaCut++;

        lep0.SetPtEtaPhiM(ele_pt[0], ele_eta[0], ele_phi[0], ELEMASS);
        lep1.SetPtEtaPhiM(muon_pt[0], muon_eta[0], muon_phi[0], MUONMASS);
      }
      else if (lep0flavor == 'm' && lep1flavor == 'e'){
        channelCut++;

        if (muon_charge[0]*ele_charge[0] > 0) continue;
        signCut++;

        if (fabs(muon_eta[0]) > 2.4 || fabs(ele_eta[0]) > 2.5) continue;
        etaCut++;

        lep0.SetPtEtaPhiM(muon_pt[0], muon_eta[0], muon_phi[0], MUONMASS);
        lep1.SetPtEtaPhiM(ele_pt[0], ele_eta[0], ele_phi[0], ELEMASS);
      }
      else continue;
    }
    if ((lep0+lep1).M() < 20) continue;
    dilepmassCut++;

    if (nJet < 2) continue;
    vector<pair<int, float> > jet_index_corrpt;
    for (int i=0; i<nJet; i++){

      float eta = jet_eta[i];
      float pt = jet_pt[i];
      float area = jet_area[i];

      jetCorrector->setJetEta(eta);
      jetCorrector->setJetPt(pt);
      jetCorrector->setJetA(area);
      jetCorrector->setRho(rho); 
      double corr_pt = jetCorrector->getCorrection() * pt;
      jet_index_corrpt.push_back( make_pair(i, corr_pt) );
    }
    sort(jet_index_corrpt.begin(), jet_index_corrpt.end(), sortJetPt);

    int jet0index = jet_index_corrpt[0].first, jet1index = jet_index_corrpt[1].first;
    float jet0pt = jet_index_corrpt[0].second, jet1pt = jet_index_corrpt[1].second;

    if ( fabs(jet_eta[jet0index]) > 2.5 ) continue;
    if ( jet0pt < 100 || jet1pt < 50 ) continue;
    jetCut++;

    FillHist1D("nEle", nEle, weight);
    FillHist1D("nMuon", nMuon, weight);
    FillHist1D("nJet", nJet, weight);

    FillHist1D("lep0pt", lep0.Pt(), weight);
    FillHist1D("lep0eta", lep0.Eta(), weight);
    FillHist1D("lep1pt", lep1.Pt(), weight);
    FillHist1D("lep1eta", lep1.Eta(), weight);
    FillHist1D("dilepmass", (lep0+lep1).M(), weight);

    FillHist1D("jet0pt", jet0pt, weight);
    FillHist1D("jet0eta", jet_eta[jet0index], weight);
    FillHist1D("jet1pt", jet1pt, weight);
    FillHist1D("jet1eta", jet_eta[jet1index], weight);
  }
  cout << difftime(time(NULL), start) << " s" << endl;

  //Cutflow Table//

  cout<<"====================================================================================================================="<< "\n" ;
  cout<<"                                         Cut Flow Table                                                              "<< "\n" ;
  cout<<"====================================================================================================================="<< "\n" ;

  cout<<      "                                  |||          Nevent          |||   Relative Efficiency    |||     Efficiency      " << "\n" ;
  cout<< Form("        Initial                   |||         %10i          |||           %1.3f          |||       %1.6f         ",countEvts,float(countEvts)/countEvts,float(countEvts)/countEvts) << "\n";
  cout<< Form("  Passed Dilepton selection       |||         %10i          |||           %1.3f          |||       %1.6f         ",countDilep,float(countDilep)/countEvts,float(countDilep)/countEvts) << "\n";
  cout<< Form("  Passed lepton Pt Cut            |||         %10i          |||           %1.3f          |||       %1.6f         ",countLeppt,float(countLeppt)/countDilep,float(countLeppt)/countEvts) << "\n";
  cout<< Form("  Passed Dilepton Mass Cut        |||         %10i          |||           %1.3f          |||       %1.6f         ",countDilepmass,float(countDilepmass)/countLeppt,float(countDilepmass)/countEvts) << "\n";
  cout<< Form("  Passed Leading Jet Pt_eta cut   |||         %10i          |||           %1.3f          |||       %1.6f         ",countJetpteta,float(countJetpteta)/countDilepmass,float(countJetpteta)/countEvts) << "\n";
  cout<< Form("  Passed MET Filters              |||         %10i          |||           %1.3f          |||       %1.6f         ",countMet,float(countMet)/countJetpteta,float(countMet)/countEvts) << "\n";
  cout<< Form("  Passed Correct Channel          |||         %10i          |||           %1.3f          |||       %1.6f         ",channelCut,float(channelCut)/countMet,float(channelCut)/countEvts) << "\n";
  cout<< Form("  Passed Opposite Lepton Sign     |||         %10i          |||           %1.3f          |||       %1.6f         ",signCut,float(signCut)/channelCut,float(signCut)/countEvts) << "\n";
  cout<< Form("  Passed HLT Trigger              |||         %10i          |||           %1.3f          |||       %1.6f         ",trigCut,float(trigCut)/signCut,float(trigCut)/countEvts) << "\n";
  cout<< Form("  Passed lepton pt cut            |||         %10i          |||           %1.3f          |||       %1.6f         ",lepptCut,float(lepptCut)/trigCut,float(lepptCut)/countEvts) << "\n";
  cout<< Form("  Passed third Lepton Cut         |||         %10i          |||           %1.3f          |||       %1.6f         ",thirdLepCut,float(thirdLepCut)/lepptCut,float(thirdLepCut)/countEvts) << "\n";
  cout<< Form("  Passed Eta Cut                  |||         %10i          |||           %1.3f          |||       %1.6f         ",etaCut,float(etaCut)/thirdLepCut,float(etaCut)/countEvts) << "\n";
  cout<< Form("  Passed dilepton mass Cut        |||         %10i          |||           %1.3f          |||       %1.6f         ",dilepmassCut,float(dilepmassCut)/etaCut,float(dilepmassCut)/countEvts) << "\n";
  cout<< Form("  Passed jet Cut                  |||         %10i          |||           %1.3f          |||       %1.6f         ",jetCut,float(jetCut)/dilepmassCut,float(jetCut)/countEvts) << "\n";

  //Write Histograms//

  TFile* outFile = new TFile(outName.data(),"RECREATE");
  outFile->cd();

  for (map<TString, TH1*>::iterator hid = m_Histos1D.begin(); hid != m_Histos1D.end(); hid++)
    hid->second->Write();

  outFile->Write();
  delete outFile;
  outFile = 0;
}

void FillHist1D(const TString& histName, const Double_t& value, const double& weight) 
{
  map<TString, TH1*>::iterator hid=m_Histos1D.find(histName);
  if (hid==m_Histos1D.end())
    cout << "%FillHist1D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value, weight);
}

void setPars(const string& parFile){

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
    else if (var == "inName") inName = line;
    else if (var == "outName") outName = line;
    else if (var == "channel") channel = line;
    else if (var == "era") era = line;
    else if (var == "jet_type") jet_type = line;
    else if (var == "weight") weight = stod(line);
  }
  file.close();
}

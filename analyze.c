#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <cmath>
#include <iomanip>
#include <vector>
#include <map>
#include <utility>
#include <time.h>

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
using namespace std;

void FillHist1D(const TString& histName, const Double_t& value, double weight);
void setPars(string parFile); 

map<TString, TH1*> m_Histos1D;

bool ISMC = false;
string inName = "analysis_Data.root";
string channel = "mm";
string era = "Spring16_25nsV10BCD";
string jet_type = "AK4PFchs";

const int MAXJET = 50;
const int MAXLEP = 20;

void analyze(string parFile = "pars.txt"){

  setPars(parFile);

  //Open Files//

  TString outName = "hist_" + inName + "_" + channel + ".root";
  TFile* inFile = TFile::Open(inName.data());

  TTree* T = (TTree*) inFile->Get("T");
  Long64_t nEntries = T->GetEntries();
  cout << nEntries << " Events" << endl;

  //MET Filters//

  int metCut = 0;
  for (int j=1; j<=6; j++){
    vector<int> *vFilter = (vector<int>*) inFile->Get( Form("filter_failed;%i", j) );
    for (int i=0; i<vFilter->size(); i++) metCut += vFilter->at(i);
  }
  cout << "MET cut: " << metCut << endl;

  //Histograms//

  TString hname = "nJet";
  m_Histos1D[hname] = new TH1F(hname,hname,MAXJET,0,MAXJET);
  hname = "jet_pt";
  m_Histos1D[hname] = new TH1F(hname,hname,200,0,2000);
  hname = "jet_eta";
  m_Histos1D[hname] = new TH1F(hname,hname,100,-5,5);

  hname = "lep0pt";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,500);
  hname = "lep0eta";
  m_Histos1D[hname] = new TH1F(hname,hname,100,-5,5);
  hname = "lep1pt";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,500);
  hname = "lep1eta";
  m_Histos1D[hname] = new TH1F(hname,hname,100,-5,5);
  hname = "lepmass";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,500);

  //Jet Corrections//

  string data = "DATA";
  if (ISMC) data = "MC";
  era += "_" + data;

  JetCorrectorParameters *ResJetPar = new JetCorrectorParameters("../../" + era + "/" + era + "_L2L3Residual_" + jet_type + ".txt"); 
  JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("../../" + era + "/" + era + "_L3Absolute_" + jet_type + ".txt");
  JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("../../" + era + "/" + era + "_L2Relative_" + jet_type + ".txt");
  JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("../../" + era + "/" + era + "_L1FastJet_" + jet_type + ".txt");

  vector<JetCorrectorParameters> jetPars;
  jetPars.push_back(*L1JetPar);
  jetPars.push_back(*L2JetPar);
  jetPars.push_back(*L3JetPar);
  jetPars.push_back(*ResJetPar);

  FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(jetPars);

  //Set Branches//

  double weight = 1;

  //ULong64_t event;
  //int run, lumi, bx;

  if (!ISMC){
    //T->SetBranchAddress("event", &event);
    //T->SetBranchAddress("run", &run);
    //T->SetBranchAddress("lumi", &lumi);
    //T->SetBranchAddress("bx", &bx);
  }

  //vector<bool> *trig_failed = 0;
  //vector<float> *trig_prescale = 0;
  //vector<string> *trig_name = 0;

  //T->SetBranchAddress("trig_failed", &trig_failed);
  //T->SetBranchAddress("trig_prescale", &trig_prescale);
  //T->SetBranchAddress("trig_name", &trig_name);

  char lep0flavor, lep1flavor;
  T->SetBranchAddress("lep0flavor", &lep0flavor);
  T->SetBranchAddress("lep1flavor", &lep1flavor);

  int nMuon=MAXLEP, muon_charge[nMuon];
  float muon_eta[nMuon], muon_pt[nMuon], muon_phi[nMuon], muon_mass[nMuon];

  T->SetBranchAddress("nMuon", &nMuon);
  T->SetBranchAddress("muon_charge", muon_charge);
  T->SetBranchAddress("muon_eta", muon_eta);
  T->SetBranchAddress("muon_pt", muon_pt);
  T->SetBranchAddress("muon_phi", muon_phi);
  T->SetBranchAddress("muon_mass", muon_mass);

  int nEle=MAXLEP, ele_charge[nEle];
  float ele_eta[nEle], ele_pt[nEle], ele_phi[nEle], ele_mass[nEle];

  T->SetBranchAddress("nEle", &nEle);
  T->SetBranchAddress("ele_charge", ele_charge);
  T->SetBranchAddress("ele_eta", ele_eta);
  T->SetBranchAddress("ele_pt", ele_pt);
  T->SetBranchAddress("ele_phi", ele_phi);
  T->SetBranchAddress("ele_mass", ele_mass);

  int nJet=MAXJET;
  float jet_eta[nJet], jet_pt[nJet], jet_area[nJet];

  T->SetBranchAddress("nJet", &nJet);
  T->SetBranchAddress("jet_eta", jet_eta);
  T->SetBranchAddress("jet_pt", jet_pt);
  T->SetBranchAddress("jet_area", jet_area);

  float rho;
  T->SetBranchAddress("rho", &rho);

  //Loop Over Entries//
  int etaCut = 0;
  int signCut = 0;
  time_t start = time(NULL);

  for (Long64_t n=0; n<nEntries; n++){
    if (n % 1000 == 0)
      cout << "Processing Event " << n+1 << endl;
    T->GetEntry(n);

    TLorentzVector lep0, lep1;
    int lep0charge, lep1charge;

    if (channel == "mm"){
      if (lep0flavor == 'm' && lep1flavor == 'm'){
        if (muon_eta[0] > 2.4 || muon_eta[1] > 2.4) { etaCut++; continue;}

        lep0.SetPtEtaPhiM(muon_pt[0], muon_eta[0], muon_phi[0], muon_mass[0]);
        lep1.SetPtEtaPhiM(muon_pt[1], muon_eta[1], muon_phi[1], muon_mass[1]);
        lep0charge = muon_charge[0]; lep1charge = muon_charge[1];
      }
      else continue;
    }
    else if (channel == "ee"){
      if (lep0flavor == 'e' && lep1flavor == 'e'){
        if (ele_eta[0] > 2.5 || ele_eta[1] > 2.5) { etaCut++; continue;}

        lep0.SetPtEtaPhiM(ele_pt[0], ele_eta[0], ele_phi[0], ele_mass[0]);
        lep1.SetPtEtaPhiM(ele_pt[1], ele_eta[1], ele_phi[1], ele_mass[1]);
        lep0charge = ele_charge[0]; lep1charge = ele_charge[1];
      }
      else continue;
    }
    else{
      if (lep0flavor == 'e' && lep1flavor == 'm'){
        if (ele_eta[0] > 2.5 || muon_eta[0] > 2.4) { etaCut++; continue;}

        lep0.SetPtEtaPhiM(ele_pt[0], ele_eta[0], ele_phi[0], ele_mass[0]);
        lep1.SetPtEtaPhiM(muon_pt[0], muon_eta[0], muon_phi[0], muon_mass[0]);
        lep0charge = ele_charge[0]; lep1charge = muon_charge[0];
      }
      else if (lep0flavor == 'm' && lep1flavor == 'e'){
        if (muon_eta[0] > 2.4 || ele_eta[0] > 2.5) { etaCut++; continue;}

        lep0.SetPtEtaPhiM(muon_pt[0], muon_eta[0], muon_phi[0], muon_mass[0]);
        lep1.SetPtEtaPhiM(ele_pt[0], ele_eta[0], ele_phi[0], ele_mass[0]);
        lep0charge = muon_charge[0]; lep1charge = ele_charge[0];
      }
      else continue;
    }
    //cout << lep0flavor << lep0charge << "\t" << lep0.Pt() << "\t" << endl;
    //cout << lep1flavor << lep1charge << "\t" << lep1.Pt() << "\t" << endl << endl;

    if (lep0charge*lep1charge > 0) { signCut++; continue; }

    FillHist1D("lep0pt", lep0.Pt(), weight);
    FillHist1D("lep0eta", lep0.Eta(), weight);
    FillHist1D("lep1pt", lep1.Pt(), weight);
    FillHist1D("lep1eta", lep1.Eta(), weight);
    FillHist1D("lepmass", (lep0+lep1).M(), weight);

    FillHist1D("nJet", nJet, weight);
    for (int i=0; i<nJet; i++){

      float eta = jet_eta[i];
      float pt = jet_pt[i];
      float area = jet_area[i];

      JetCorrector->setJetEta(eta);
      JetCorrector->setJetPt(pt);
      JetCorrector->setJetA(area);
      JetCorrector->setRho(rho); 
      double corr_pt = JetCorrector->getCorrection() * pt;

      FillHist1D("jet_pt", corr_pt, weight);
      FillHist1D("jet_eta", eta, weight);
    }

  }
  cout << difftime(time(NULL), start) << " s" << endl;
  cout << "eta cut: " << etaCut << endl;
  cout << "sign cut: " << signCut << endl;

  //Write Histograms//

  TFile* outFile = new TFile(outName,"RECREATE");
  outFile->cd();

  for (map<TString, TH1*>::iterator hid = m_Histos1D.begin(); hid != m_Histos1D.end(); hid++)
    hid->second->Write();

  outFile->Write();
  delete outFile;
  outFile = 0;
}

void FillHist1D(const TString& histName, const Double_t& value, double weight) 
{
  map<TString, TH1*>::iterator hid=m_Histos1D.find(histName);
  if (hid==m_Histos1D.end())
    cout << "%FillHist1D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value, weight);
}

void setPars(string parFile){

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
    else if (var == "channel") channel = line;
    else if (var == "era") era = line;
    else if (var == "jet_type") jet_type = line;
  }
  file.close();
}

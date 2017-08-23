//charles harrington and bahareh roozbahani
//execute as plot plot_pars.txt

#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

#include <vector>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

double delta(double delta1, double delta2, string pos_neg);
void setPars(const string& parFile);
void setStyle();
void drawText();

//parameters- edit in plot_pars.txt
vector<TString> mcFileNames, sigFileNames, systematics;
vector<double> mcScales, sigScales, rebin;
map<TString, float> sys_norm;
string subplot, dataName;
TString dataFileName, outName;
TString hname, leftText, rightText;
float xmin, xmax, ymin, ymax, subymin, subymax;
bool logx, logy, plotData, plotImpact;

int main(int argc, char* argv[]) {

  if (argc == 1) { cout << "Please provide a parameter file." << endl; return -1; }

  string parFile = argv[1];
  setPars(parFile);

  setStyle();

  enum SetEnum { wjet=0, vv, st, dy, ttbar, gluon, zprime, bkg, undefined };
  TString labels[] = { "W+Jets", "VV", "Single-Top", "Z/#gamma*#rightarrowl^{+}l^{-}", "t#bar{t}", "g_{kk} 3 TeV(#sigma=10 pb)", "Z'   3 TeV(#sigma=10 pb)", "Background", "Undefined" };

  TFile* dataFile = TFile::Open(dataFileName);
  TH1F* h_Data = (TH1F*) dataFile->FindObjectAny(hname);

  if (rebin.size() == 1) h_Data->Rebin(rebin[0]);
  else h_Data = (TH1F*) h_Data->Rebin(rebin.size()-1, "h_Data", &rebin[0]);
  int nBins = h_Data->GetNbinsX();

  h_Data->SetMarkerStyle(20);
  h_Data->SetLineColor(kBlack);
  h_Data->SetMarkerSize(0.65);

  map<SetEnum, TH1F*> m_MCs, m_bkg;
  map<pair<SetEnum, TString>, TH1F*> m_MCUPs, m_MCDNs, m_bkgUPs, m_bkgDNs;
  for (int i=0,n=mcFileNames.size(); i<n; i++) {
    TFile* mcFile = TFile::Open(mcFileNames[i]);

    TH1F* h_MC = (TH1F*) mcFile->FindObjectAny(hname);
    //h_MC->Scale(mcScales[i]);

    if (rebin.size() == 1) h_MC->Rebin(rebin[0]);
    else h_MC = (TH1F*) h_MC->Rebin(rebin.size()-1, "h_MC", &rebin[0]);

    SetEnum dataset = undefined;
    if ( mcFileNames[i].Contains("ttbar", TString::kIgnoreCase) )     dataset = ttbar;
    else if ( mcFileNames[i].Contains("dy", TString::kIgnoreCase) )   dataset = dy;
    else if ( mcFileNames[i].Contains("wjet", TString::kIgnoreCase) ) dataset = wjet;
    else if ( mcFileNames[i].Contains("st", TString::kIgnoreCase)
           || mcFileNames[i].Contains("sat", TString::kIgnoreCase) )  dataset = st;
    else if ( mcFileNames[i].Contains("ww", TString::kIgnoreCase)
           || mcFileNames[i].Contains("wz", TString::kIgnoreCase)
           || mcFileNames[i].Contains("zz", TString::kIgnoreCase) )   dataset = vv;

    if ( m_MCs.find(dataset) == m_MCs.end() ) m_MCs[dataset] = h_MC;
    else m_MCs[dataset]->Add(h_MC);

    if ( m_bkg.find(bkg) == m_bkg.end() ) m_bkg[bkg] = (TH1F*) h_MC->Clone("h_bkg"); //second map to start with h_MC (needs new object)
    else m_bkg[bkg]->Add(h_MC);

    for (unsigned int i_sys = 0; i_sys != systematics.size(); ++i_sys) {
      TString sys = systematics[i_sys];
      TFile* mcFileUP = TFile::Open( mcFileNames[i]( 0, mcFileNames[i].Index(".root") ) + "_" + sys + "UP.root" );
      TFile* mcFileDN = TFile::Open( mcFileNames[i]( 0, mcFileNames[i].Index(".root") ) + "_" + sys + "DOWN.root" );

      TH1F* h_MCUP = (TH1F*) mcFileUP->FindObjectAny(hname);
      //h_MCUP->Scale(mcScales[i]);
      TH1F* h_MCDN = (TH1F*) mcFileDN->FindObjectAny(hname);
      //h_MCDN->Scale(mcScales[i]);

      if (rebin.size() == 1) h_MCUP->Rebin(rebin[0]);
      else h_MCUP = (TH1F*) h_MCUP->Rebin(rebin.size()-1, "h_MCUP", &rebin[0]);
      if (rebin.size() == 1) h_MCDN->Rebin(rebin[0]);
      else h_MCDN = (TH1F*) h_MCDN->Rebin(rebin.size()-1, "h_MCDN", &rebin[0]);

      pair<SetEnum, TString> keypair = make_pair(dataset, sys);
      if ( m_MCUPs.find(keypair) == m_MCUPs.end() ) m_MCUPs[keypair] = h_MCUP;
      else m_MCUPs[keypair]->Add(h_MCUP);
      if ( m_MCDNs.find(keypair) == m_MCDNs.end() ) m_MCDNs[keypair] = h_MCDN;
      else m_MCDNs[keypair]->Add(h_MCDN);

      TH1F* h_bkgUP = 0, *h_bkgDN = 0;
      if (sys == "topPt_weight" && dataset != ttbar) { h_bkgUP = h_MC; h_bkgDN = h_MC; } // topPt_weight only affects ttbar
      else { h_bkgUP = h_MCUP; h_bkgDN = h_MCDN; }

      keypair = make_pair(bkg, sys);
      if ( m_bkgUPs.find(keypair) == m_bkgUPs.end() ) m_bkgUPs[keypair] = (TH1F*) h_bkgUP->Clone("h_bkgUP");
      else m_bkgUPs[keypair]->Add(h_bkgUP);
      if ( m_bkgDNs.find(keypair) == m_bkgDNs.end() ) m_bkgDNs[keypair] = (TH1F*) h_bkgDN->Clone("h_bkgDN");
      else m_bkgDNs[keypair]->Add(h_bkgDN);
    }
  }

  THStack* mcStack = new THStack();
  for (map<SetEnum, TH1F*>::const_iterator i_MC=m_MCs.begin(); i_MC != m_MCs.end(); ++i_MC) mcStack->Add( i_MC->second );

  if ( !m_MCs.empty() ) {
    m_MCs[ttbar]->SetLineColor(2);
    m_MCs[ttbar]->SetFillColor(2);

    m_MCs[dy]->SetLineColor(8);
    m_MCs[dy]->SetFillColor(8);

    m_MCs[wjet]->SetLineColor(4);
    m_MCs[wjet]->SetFillColor(4);

    m_MCs[vv]->SetLineColor(6);
    m_MCs[vv]->SetFillColor(6);

    m_MCs[st]->SetLineColor(28);
    m_MCs[st]->SetFillColor(28);
  }

  map<SetEnum, TH1F*> m_sigs;
  map<pair<SetEnum, TString>, TH1F*> m_sigUPs, m_sigDNs;
  for (int i=0,n=sigFileNames.size(); i<n; i++) {

    TFile* sigFile = TFile::Open(sigFileNames[i]);
    TH1F* h_sig = (TH1F*) sigFile->FindObjectAny(hname);
    h_sig->Scale(sigScales[i]);

    if (rebin.size() == 1) h_sig->Rebin(rebin[0]);
    else h_sig = (TH1F*) h_sig->Rebin(rebin.size()-1, "h_sig", &rebin[0]);

    SetEnum dataset = undefined;
    if ( sigFileNames[i].Contains("zprime", TString::kIgnoreCase) )     dataset = zprime;
    else if ( sigFileNames[i].Contains("gluon", TString::kIgnoreCase) ) dataset = gluon;
    m_sigs[dataset] = h_sig;

    for (unsigned int i_sys = 0; i_sys != systematics.size(); ++i_sys) {
      TString sys = systematics[i_sys];
      TFile* sigFileUP = TFile::Open( sigFileNames[i]( 0, sigFileNames[i].Index(".root") ) + "_" + sys + "UP.root" );
      TFile* sigFileDN = TFile::Open( sigFileNames[i]( 0, sigFileNames[i].Index(".root") ) + "_" + sys + "DOWN.root" );

      TH1F* h_sigUP = (TH1F*) sigFileUP->FindObjectAny(hname);
      h_sigUP->Scale(sigScales[i]);
      TH1F* h_sigDN = (TH1F*) sigFileDN->FindObjectAny(hname);
      h_sigDN->Scale(sigScales[i]);

      if (rebin.size() == 1) h_sigUP->Rebin(rebin[0]);
      else h_sigUP = (TH1F*) h_sigUP->Rebin(rebin.size()-1, "h_sigUP", &rebin[0]);
      if (rebin.size() == 1) h_sigDN->Rebin(rebin[0]);
      else h_sigDN = (TH1F*) h_sigDN->Rebin(rebin.size()-1, "h_sigDN", &rebin[0]);

      pair<SetEnum, TString> keypair = make_pair(dataset, sys);
      m_sigUPs[keypair] = h_sigUP;
      m_sigDNs[keypair] = h_sigDN;
    }
  }

  if ( !m_sigs.empty() ) {
    m_sigs[zprime]->SetLineColor(12);
    m_sigs[zprime]->SetLineWidth(2);

    m_sigs[gluon]->SetLineColor(9);
    m_sigs[gluon]->SetLineWidth(2);
    m_sigs[gluon]->SetLineStyle(2);
  }

  string channel = "em";
  if      (dataFileName.Contains("mm", TString::kIgnoreCase)) channel = "mm";
  else if (dataFileName.Contains("ee", TString::kIgnoreCase)) channel = "ee";

  //Combine sys_norm and systematics
  for (map<TString, float>::const_iterator it = sys_norm.begin(); it != sys_norm.end(); ++it) {
    TString sys = it->first;

    if( channel == "mm" && (sys == "eltrig" || sys == "elid" || sys == "eliso") ) continue ;
    if( channel == "ee" && (sys == "mutrig" || sys == "muid" || sys == "muiso") ) continue ;
    if( channel == "em" &&  sys == "eltrig" )                                     continue ;
    systematics.push_back(sys);
  }

  TGraphAsymmErrors* background = new TGraphAsymmErrors(); //graph for background uncertainty
  background->SetFillStyle(3013);
  background->SetFillColor(kBlack);

  for (int pt=0; pt<nBins; pt++) {
    int bin = pt+1;
    double nom = m_bkg[bkg]->GetBinContent(bin), errorUP=0, errorDN=0;

    for (unsigned int i_sys = 0; i_sys != systematics.size(); ++i_sys) {
      TString sys = systematics[i_sys];
      double deltaUP, deltaDN;

      if (sys_norm.find(sys) != sys_norm.end()) {  // normalization-only systematics
        double perEvent_sys = sys_norm[sys] ; // per event systematics

        if (sys == "muid" || sys == "muiso") {
          if (channel == "mm") perEvent_sys *= 2.;
          perEvent_sys *= nom;
        }
        else if (sys == "eltrig" || sys == "elid" || sys == "eliso") {
          if (channel == "ee") perEvent_sys *= 2.;
          perEvent_sys *= nom;
        }
        else if (sys == "lumi" || sys == "mutrig") perEvent_sys *= nom;

        else if (sys=="sig_tt") perEvent_sys *= m_MCs[ttbar]->GetBinContent(bin);

        else if (sys=="sig_dy") perEvent_sys *= m_MCs[dy]->GetBinContent(bin);

        else if (sys=="sig_st") perEvent_sys *= m_MCs[st]->GetBinContent(bin);

        else if (sys=="sig_db") perEvent_sys *= m_MCs[vv]->GetBinContent(bin);

        else perEvent_sys = 0.;

        deltaUP = perEvent_sys;
        deltaDN = -1 * perEvent_sys;
      }
      else {
        deltaUP = m_bkgUPs[make_pair(bkg, sys)]->GetBinContent(bin) - nom;
        deltaDN = m_bkgDNs[make_pair(bkg, sys)]->GetBinContent(bin) - nom;
      }
      double d1 = delta(deltaUP, deltaDN, "+");
      double d2 = delta(deltaUP, deltaDN, "-");

      errorUP += d1*d1;
      errorDN += d2*d2;
    }
    background->SetPoint(pt, h_Data->GetBinCenter(bin), nom);
    if (nom != 0) {
      background->SetPointEYhigh(pt, sqrt(errorUP) );
      background->SetPointEYlow(pt, sqrt(errorDN) );
      background->SetPointEXhigh(pt, h_Data->GetBinWidth(bin)/2);
      background->SetPointEXlow(pt, h_Data->GetBinWidth(bin)/2);
    }
  }

  TCanvas* c = new TCanvas("c", "c", 600, 600);

  float b_scale = 0.3, t_scale = 1 - b_scale;
  TPad* top = new TPad("top", "top", 0, b_scale, 1, 1);
  TPad* bottom = new TPad("bottom", "bottom", 0, 0, 1, b_scale);

  if ( plotData && (subplot=="ratio" || subplot=="pull") ) {
    top->SetTopMargin(0.05);
    top->SetBottomMargin(0.05);
    top->Draw();
    bottom->SetBottomMargin(0.35);
    bottom->SetGridy();
    bottom->Draw();
    top->cd();
  }
  if (logx) gPad->SetLogx();
  if (logy) gPad->SetLogy();

  TString xtitle = "No Title";
  string keytitle = hname(2, hname.Length()).Data();
  map<string, string> xtitles = {{"dilepmass","M_{ll} (Gev)"},{"lep0eta","#eta_{Leading Lepton}"},
  {"lep1eta","#eta_{Subleading Lepton}"},{"lep0pt","Leading Lepton p_{T} (GeV)"},{"lep1pt","Subleading Lepton p_{T} (GeV)"},
  {"jet0eta","#eta_{Leading Jet}"},{"jet1eta","#eta_{Subleading Jet}"},{"jet0pt","Leading Jet p_{T} (GeV)"},
  {"jet1pt","Subleading Jet p_{T} (GeV)"},{"nEle","Number of Electrons"},{"nMuon","Number of Muons"},
  {"nJet","Number of Jets"},{"nEleDiff","N_{Electrons}-N_{Good Electrons}"},{"nMuonDiff","N_{Muons}-N_{Good Muons} "},
  {"nJetDiff","N_{Jets}-N_{Good Jets} "},{"nGoodEle","N_{Good Electrons}"},{"nGoodMuon","N_{Good Muons}"},
  {"nGoodJet","N_{Good Jets}"},{"jethT","H_{T} (GeV)"},{"sT","H_{T}^{L} (GeV)"},{"sT_met","S_{T} (GeV)"},{"metpt","MET p_{T} (GeV)"},
  {"jet0btag","btag_{Leading Jet}"},{"jet1btag","btag_{Subeading Jet}"},{"nbtag","Number of btagged Jets"},{"sumrmin","#DeltaR_{min0} + #DeltaR_{min1}"},
  {"metcorrpt","MET p_{T} (GeV)"},{"muonD0","Muon D_{0} (cm)"},{"muonDz","Muon D_{z} (cm)"},{"rmin0","#DeltaR_{min}(leading lepton, jet)"},
  {"rmin1","#DeltaR_{min}(subleading lepton, jet)"},{"lep0perp","Leading Lepton p_{T}^{rel} (GeV)"},{"lep1perp","Subleading Lepton p_{T}^{rel} (GeV)"},
  {"rl0cleanj","#DeltaR(leading lepton, cleaned jet)"},{"rl1cleanj","#DeltaR(subleading lepton, cleaned jet)"},{"rl0l1","#DeltaR(leading lepton, subleading lepton)"},
  {"jet0phi","#phi_{Leading Jet}"}, {"jet1phi","#phi_{Subleading Jet}"}, {"lep0phi","#phi_{Leading Lepton}"}, {"lep1phi","#phi_{Subleading Lepton}"},
  {"lepept","electron p_{T} (GeV)"}, {"lepmpt","muon p_{T} (GeV)"},{"rbl","#DeltaR(b quark, lepton)"},{"minjet0pt","Jet p_{T}^{rmin leading lepton} (GeV)"},
  {"minjet1pt","Jet p_{T}^{rmin subleading lepton} (GeV)"},{"cleanjet0pt","Jet p_{T}^{cleaned from leading lepton} (GeV)"},
  {"cleanjet1pt","Jet p_{T}^{cleaned from subleading lepton} (GeV)"},{"masslmin0","M_{leading lep,rmin jet} (Gev)"},{"masslmin1","M_{subleading lep,rmin jet} (Gev)"},
  {"masslljjm","M_{lljjmet} (Gev)"},{"dphi_jet0met","#Delta #phi_{Leading Jet,MET}"},{"dphi_jet1met","#Delta #phi_{Subleading Jet,MET}"}};
  if (xtitles.find(keytitle) != xtitles.end()) xtitle = xtitles[keytitle];

  TH1F* hist = 0;
  if (rebin.size() == 1) hist = new TH1F("hist", "hist", nBins, h_Data->GetBinLowEdge(1), h_Data->GetBinLowEdge(nBins+1));
  else hist = new TH1F("hist", "hist", rebin.size()-1, &rebin[0]);

  if ( plotData && (subplot=="ratio" || subplot=="pull") ) {
    hist->GetXaxis()->SetTickLength(0.03/t_scale);
    hist->GetXaxis()->SetLabelSize(0);
    hist->GetYaxis()->SetTitleSize(0.03/t_scale);
    hist->GetYaxis()->SetTitleOffset(0.95);
    hist->GetYaxis()->SetLabelSize(0.03/t_scale);
  }
  else {
    hist->GetXaxis()->SetTitle(xtitle);
    hist->GetXaxis()->SetTitleSize(0.03/t_scale);
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetXaxis()->SetLabelSize(0.03/t_scale);

    hist->GetYaxis()->SetTitleSize(0.03/t_scale);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetLabelSize(0.03/t_scale);

    if (logx) { hist->GetXaxis()->SetNoExponent(); hist->GetXaxis()->SetMoreLogLabels(); }
  }

  hist->GetXaxis()->SetNdivisions(5, 5, 0);
  hist->GetXaxis()->SetRangeUser(xmin, xmax);

  if (logy) hist->GetYaxis()->SetRangeUser(0.1, int(h_Data->GetMaximum()*100) );
  else      hist->GetYaxis()->SetRangeUser(ymin, int(h_Data->GetMaximum()*1.5) );

  hist->Draw();
  int legEntries = m_MCs.size() + m_sigs.size();
  if (plotData) legEntries++;

  TLegend* leg = 0;
  if ( plotData && (subplot=="ratio" || subplot=="pull") ) {
    leg = new TLegend(.65,.9-.06*legEntries,.85,.9);
    leg->SetTextSize(0.04);
  }
  else {
    leg = new TLegend(.65,.9-.04*legEntries,.85,.9);
    leg->SetTextSize(0.03);
  }
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  for (map<SetEnum, TH1F*>::const_reverse_iterator i_MC=m_MCs.rbegin(); i_MC != m_MCs.rend(); ++i_MC)
    leg->AddEntry(i_MC->second, labels[i_MC->first].Data(), "F");
  for (map<SetEnum, TH1F*>::const_iterator i_sig=m_sigs.begin(); i_sig != m_sigs.end(); ++i_sig)
    leg->AddEntry(i_sig->second, labels[i_sig->first].Data(), "L");

  mcStack->Draw("samehist");
  for (map<SetEnum, TH1F*>::const_iterator i_sig=m_sigs.begin(); i_sig != m_sigs.end(); ++i_sig)
    i_sig->second->Draw("samehist");

  leg->AddEntry(background, "Bkg Uncertainty", "F");
  background->Draw("sameE2");

  if (plotData) {
    h_Data->Draw("samePE");
    leg->AddEntry(h_Data, Form( "%s", dataName.data() ), "PLE");
  }

  gPad->RedrawAxis();
  leg->Draw();

  drawText();

  TH1F* bhist = 0;
  if ( plotData && (subplot=="ratio" || subplot=="pull") ) {
    bottom->cd();
    if (logx) gPad->SetLogx();

    bhist = (TH1F*) hist->Clone("bhist");

    TH1F* hsubplot = (TH1F*) h_Data->Clone("hsubplot");

    TGraphAsymmErrors* bkg_env = (TGraphAsymmErrors*) background->Clone("bkg_env");
    double *yarray = background->GetY(), *eyhigh = background->GetEYhigh(), *eylow = background->GetEYlow();

    TString subytitle = "#frac{Data}{Background}";
    if (subplot=="ratio") {
      bhist->GetYaxis()->SetTitleOffset(0.53);
      for (int pt=0; pt<nBins; pt++) {
        if (yarray[pt] == 0) continue;
        int bin=pt+1;

        hsubplot->SetBinContent(bin, h_Data->GetBinContent(bin)/yarray[pt]);
        hsubplot->SetBinError(bin, h_Data->GetBinError(bin)/yarray[pt]);

        bkg_env->SetPoint(pt, h_Data->GetBinCenter(bin), 1.);
        bkg_env->SetPointEYhigh( pt, eyhigh[pt]/yarray[pt] );
        bkg_env->SetPointEYlow( pt, eylow[pt]/yarray[pt] );
      }
    }
    else {
      subytitle = "Pull";
      bhist->GetYaxis()->SetTitleOffset(0.43);
      for (int pt=0; pt<nBins; pt++) {
        if (yarray[pt] == 0) continue;
        int bin=pt+1;

        double diff = h_Data->GetBinContent(bin)-yarray[pt];
        double sigma = diff>0 ? sqrt( h_Data->GetBinError(bin)*h_Data->GetBinError(bin) + eyhigh[pt]*eyhigh[pt] )
                              : sqrt( h_Data->GetBinError(bin)*h_Data->GetBinError(bin) + eylow[pt]*eylow[pt] );

        if (sigma == 0) continue;
        hsubplot->SetBinContent(bin, diff/sigma);
      }
    }

    bhist->GetXaxis()->SetLabelSize(0.03/b_scale);
    bhist->GetXaxis()->SetTickLength(0.03/b_scale);
    bhist->GetXaxis()->SetTitleSize(0.03/b_scale);
    bhist->GetXaxis()->SetTitleOffset(1.3);
    bhist->GetXaxis()->SetTitle(xtitle);
    bhist->GetXaxis()->SetRangeUser(xmin, xmax);
    if (logx) { bhist->GetXaxis()->SetNoExponent(); bhist->GetXaxis()->SetMoreLogLabels(); }

    bhist->GetYaxis()->SetRangeUser(subymin, subymax);
    bhist->GetYaxis()->SetNdivisions(5, 3, 0);
    bhist->GetYaxis()->SetLabelSize(0.03/b_scale);
    bhist->GetYaxis()->SetTitle(subytitle);
    bhist->GetYaxis()->SetTitleSize(0.03/b_scale);
    bhist->GetYaxis()->CenterTitle(true);

    bhist->Draw();

    gStyle->SetHatchesLineWidth(1);
    gStyle->SetHatchesSpacing(2);

    if (subplot=="ratio") {
      bkg_env->SetFillColor(kGray+2);
      bkg_env->SetFillStyle(3144);
      bkg_env->Draw("samePE2");
      hsubplot->Draw("samePE");
    }
    else {
      hsubplot->SetFillColor(kGray+2);
      hsubplot->SetFillStyle(3144);
      hsubplot->Draw("samehist");
    }
  }

  c->Print("./plots/" + outName + ".pdf");

  if (plotImpact) {

    TString outNames[]  = { "wjet", "vv", "st", "dy", "ttbar", "gluon", "zprime", "bkg", "undefined" };

    delete leg;
    TLegend* leg = new TLegend(.65,.9-.06*4,.85,.9);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);

    bhist->GetYaxis()->SetTitleOffset(0.53);
    if (subplot=="ratio") bhist->GetYaxis()->SetTitle("#frac{Up/Down}{Nom}");
    else                  bhist->GetYaxis()->SetTitle("#frac{Up/Down - Nom}{Nom}");

    m_MCs.insert(m_sigs.begin(), m_sigs.end()); //merge signals to mc map to take care of plots in one loop
    m_MCUPs.insert(m_sigUPs.begin(), m_sigUPs.end());
    m_MCDNs.insert(m_sigDNs.begin(), m_sigDNs.end());

    m_MCs[bkg] = m_bkg[bkg]; //merge background to mc map to take care of plots in one loop

    for (unsigned int i_sys = 0; i_sys != systematics.size(); ++i_sys) {
      TString sys = systematics[i_sys];
      if (sys_norm.find(sys) != sys_norm.end()) continue;

      m_MCUPs[make_pair(bkg, sys)] = m_bkgUPs[make_pair(bkg, sys)]; //merge background to mc map
      m_MCDNs[make_pair(bkg, sys)] = m_bkgDNs[make_pair(bkg, sys)];

      for (map<SetEnum, TH1F*>::const_iterator i_MC=m_MCs.begin(); i_MC != m_MCs.end(); ++i_MC) {

        c->Clear("D");  leg->Clear();

        SetEnum dataset = i_MC->first;
        if (sys == "topPt_weight" && dataset != ttbar && dataset != bkg) continue;

        top->cd();  hist->Draw();  drawText();

        TH1F* h_NM = i_MC->second;
        TH1F* h_UP = m_MCUPs[make_pair(dataset, sys)];
        TH1F* h_DN = m_MCDNs[make_pair(dataset, sys)];

        if (logy) hist->GetYaxis()->SetRangeUser(0.1, int(h_NM->GetMaximum()*100) );
        else      hist->GetYaxis()->SetRangeUser(ymin, int(h_NM->GetMaximum()*1.5) );

        h_NM->ResetAttFill();  h_NM->ResetAttLine();

        h_NM->SetLineColor(kBlack);  h_NM->SetLineWidth(2);  h_NM->Draw("histsame");
        h_UP->SetLineColor(8);       h_UP->SetLineWidth(2);  h_UP->Draw("histsame");
        h_DN->SetLineColor(kRed+1);  h_DN->SetLineWidth(2);  h_DN->Draw("histsame");

        leg->SetHeader( labels[dataset].Data() );  leg->AddEntry(h_NM, "Nominal", "L");
        leg->AddEntry(h_UP, sys+" Up", "L");       leg->AddEntry(h_DN, sys+" Down", "L");
        leg->Draw();

        bottom->cd();  bhist->Draw();

        TH1F* hsubUP = (TH1F*) h_UP->Clone("hsubUP");
        TH1F* hsubDN = (TH1F*) h_DN->Clone("hsubDN");
        if (subplot=="ratio") { hsubUP->Divide(h_NM); hsubDN->Divide(h_NM); }
        else                  { hsubUP->Add(h_NM, -1); hsubUP->Divide(h_NM); hsubDN->Add(h_NM, -1); hsubDN->Divide(h_NM); }
        hsubUP->Draw("histsame");  hsubDN->Draw("histsame");

        c->Print("./plots/impact/" + outName + "_" + outNames[dataset] + "_" + sys + ".pdf");
      }
    }
  }

}

void drawText() {
  TLatex text;
  text.SetNDC();

  text.SetTextSize(0.035);
  text.SetTextFont(42);
  if ( plotData && (subplot=="ratio" || subplot=="pull") ) text.DrawLatex(1-rightText.Length()/95., 0.96, rightText);
  else text.DrawLatex(1-rightText.Length()/68., 0.96, rightText);

  text.SetTextSize(0.05);
  text.SetTextFont(61);
  text.DrawLatex(0.18, 0.96, "CMS");
  //text.DrawLatex(0.18, 0.96, leftText);
  text.SetTextSize(0.03);

  //text.SetTextSize(0.04);
  //text.SetTextFont(52);
  //text.DrawLatex(0.29, 0.96, "Simulation Preliminary"); //make bool

  text.SetTextSize(0.035);
  text.SetTextFont(42);
  float textposx = 0.2, textposy = 0.9;

  if      (dataFileName.Contains("mm", TString::kIgnoreCase))  text.DrawLatex(textposx,textposy,"#bf{#mu#mu}");
  else if (dataFileName.Contains("ee", TString::kIgnoreCase))  text.DrawLatex(textposx,textposy,"#bf{ee}");
  else                                                         text.DrawLatex(textposx,textposy,"#bf{e#mu}");

  if      (hname.Contains("0_")) { text.DrawLatex(textposx,textposy-0.05,"#bf{= 0 btags}");    text.DrawLatex(textposx,textposy-0.1,"#bf{p_{T}^{j0}>100 GeV}"); }
  else if (hname.Contains("1_")) { text.DrawLatex(textposx,textposy-0.05,"#bf{#geq 1 btag}");  text.DrawLatex(textposx,textposy-0.1,"#bf{p_{T}^{j0}>100 GeV}"); }
  else if (hname.Contains("2_")) { text.DrawLatex(textposx,textposy-0.05,"#bf{= 0 btags}");    text.DrawLatex(textposx,textposy-0.1,"#bf{p_{T}^{j0}>100 GeV, p_{T}^{j1}>50 GeV}"); }
  else if (hname.Contains("3_")) { text.DrawLatex(textposx,textposy-0.05,"#bf{= 1 btag}");     text.DrawLatex(textposx,textposy-0.1,"#bf{p_{T}^{j0}>100 GeV, p_{T}^{j1}>50 GeV}"); }
  else if (hname.Contains("4_")) { text.DrawLatex(textposx,textposy-0.05,"#bf{#geq 2 btags}"); text.DrawLatex(textposx,textposy-0.1,"#bf{p_{T}^{j0}>100 GeV, p_{T}^{j1}>50 GeV}"); }
  else if (hname.Contains("5_")) { text.DrawLatex(textposx,textposy-0.05,"#bf{#geq 1 btag}");  text.DrawLatex(textposx,textposy-0.1,"#bf{p_{T}^{j0}>100 GeV, p_{T}^{j1}>50 GeV}"); }
}

void setPars(const string& parFile) {

  ifstream file(parFile);
  string line;

  while (getline(file, line)) {

    if (line.length() > 0){
      while (line.at(0) == ' ') line.erase(0, 1);
    }

    int delim_pos = line.find(' ');
    if (delim_pos == -1) continue;

    string var = line.substr(0, delim_pos);
    line.erase(0, delim_pos + 1);

    while (line.at(0) == ' ') line.erase(0, 1);
    while (line.at(line.length()-1) == ' ') line.erase(line.length()-1, line.length());

    if (var == "dataFileName")   dataFileName = line.data();
    else if (var == "dataName")  dataName = line;
    else if (var == "mcFileNames") {
      while ( (delim_pos = line.find(' ')) != -1) {
        mcFileNames.push_back( line.substr(0, delim_pos).data() );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      mcFileNames.push_back( line.data() );
    }
    else if (var == "sigFileNames") {
      while ( (delim_pos = line.find(' ')) != -1) {
        sigFileNames.push_back( line.substr(0, delim_pos).data() );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      sigFileNames.push_back( line.data() );
    }
    else if (var == "systematics") {
      while ( (delim_pos = line.find(' ')) != -1) {
        systematics.push_back( line.substr(0, delim_pos).data() );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      systematics.push_back( line.data() );
    }
    else if (var == "sys_norm") {
      while ( (delim_pos = line.find(' ')) != -1) {
        int col = line.find(':');
        sys_norm[line.substr(0, col).data()] = stof( line.substr(col+1, delim_pos-col-1) );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      int col = line.find(':');
      sys_norm[line.substr(0, col).data()] = stof( line.substr(col+1, delim_pos-col-1) );
    }
    else if (var == "mcScales") {
      while ( (delim_pos = line.find(' ')) != -1) {
        mcScales.push_back( stod( line.substr(0, delim_pos) ) );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      mcScales.push_back( stod(line) );
    }
    else if (var == "sigScales") {
      while ( (delim_pos = line.find(' ')) != -1) {
        sigScales.push_back( stod( line.substr(0, delim_pos) ) );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      sigScales.push_back( stod(line) );
    }
    else if (var == "rebin") {
      while ( (delim_pos = line.find(' ')) != -1) {
        rebin.push_back( stod( line.substr(0, delim_pos) ) );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      rebin.push_back( stod(line) );
    }
    else if (var == "outName")   outName = line.data();
    else if (var == "hname")     hname = line.data();
    else if (var == "leftText")  leftText = line.data();
    else if (var == "rightText") rightText = line.data();
    else if (var == "xmin")      xmin = stof(line);
    else if (var == "xmax")      xmax = stof(line);
    else if (var == "ymin")      ymin = stof(line);
    else if (var == "ymax")      ymax = stof(line);
    else if (var == "subymin")   subymin = stof(line);
    else if (var == "subymax")   subymax = stof(line);
    else if (var == "logx") {
      if (line == "true") logx = true;
      else logx = false;
    }
    else if (var == "logy") {
      if (line == "true") logy = true;
      else logy = false;
    }
    else if (var == "plotData") {
      if (line == "true") plotData = true;
      else plotData = false;
    }
    else if (var == "plotImpact") {
      if (line == "true") plotImpact = true;
      else plotImpact = false;
    }
    else if (var == "subplot")   subplot = line;
  }
  file.close();
}

void setStyle(){

//Style//

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(0.9);

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  //tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  tdrStyle->SetPaperSize(20.,20.);

  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->SetOptStat(0);

  tdrStyle->cd();
}

double delta(double delta1, double delta2, string plus_minus) {
// pos_neg should be either "+" or "-"

  double del = 0. ;
  if (plus_minus=="+") {
    double dd = (delta1 > delta2) ? delta1 : delta2 ;
    if (dd > 0.) del = dd ;
  }
  else if (plus_minus=="-") {
    double dd = (delta1 < delta2) ? delta1 : delta2 ;
    if (dd < 0.) del = dd ;
  }
  return del;
}

//charles harrington and bahareh roozbahani
//execute as plot plot_pars.txt

#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

#include <vector>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void setPars(const string& parFile);
void setStyle();

//parameters- edit in plot_pars.txt
vector<TString> mcFileNames, sigFileNames;
vector<double> mcScales, sigScales;
string subplot, dataName;
TString dataFileName;
TString hname, outName, leftText, rightText;
float xmin, xmax, ymin, ymax, subymin, subymax;
bool logx, logy;

int main(int argc, char* argv[]) {

  if (argc == 1) { cout << "Please provide a parameter file." << endl; return -1; }

  string parFile = argv[1];
  setPars(parFile);

  setStyle();

  TFile* dataFile = TFile::Open(dataFileName);

  //TString className = dataFile->FindObjectAny(hname)->ClassName();

  TH1F* h_Data = (TH1F*) dataFile->FindObjectAny(hname);
  //h_Data->Scale( 1 / h_Data->Integral() );  //make this bool, as well as exponents?

  h_Data->SetMarkerStyle(20);
  h_Data->SetLineColor(kBlack);
  h_Data->SetMarkerSize(0.65);

  unordered_map<string, TH1F*> m_MCs;
  for (int i=0,n=mcFileNames.size(); i<n; i++) {

    TFile* mcFile = TFile::Open(mcFileNames[i]);
    TH1F* h_MC = (TH1F*) mcFile->FindObjectAny(hname);
    //h_MC->Scale(mcScales[i]);
    h_MC->Sumw2();

    string key;
    if ( mcFileNames[i].Contains("ttbar", TString::kIgnoreCase) )     key = "t#bar{t}";
    else if ( mcFileNames[i].Contains("dy", TString::kIgnoreCase) )   key = "Z/#gamma*#rightarrowl^{+}l^{-}";
    else if ( mcFileNames[i].Contains("wjet", TString::kIgnoreCase) ) key = "W+Jets";
    else if ( mcFileNames[i].Contains("st", TString::kIgnoreCase)
           || mcFileNames[i].Contains("sat", TString::kIgnoreCase) )  key = "Single-Top";

    if ( m_MCs.find(key) == m_MCs.end() ) m_MCs[key] = h_MC;
    else m_MCs[key]->Add(h_MC);
  }

  THStack* mcStack = new THStack();
  //int color = 4;
  for (unordered_map<string, TH1F*>::const_iterator i_MC=m_MCs.begin(); i_MC != m_MCs.end(); ++i_MC) {

    mcStack->Add( i_MC->second );
    //i_MC->second->SetLineColor(color+2);
    //i_MC->second->SetFillColor(color+2);
    //i_MC->second->Scale( 1 / i_MC->second->Integral() );
    //color++;
  }

  if ( !m_MCs.empty() ) {
    m_MCs["t#bar{t}"]->SetLineColor(2);
    m_MCs["t#bar{t}"]->SetFillColor(2);

    m_MCs["Z/#gamma*#rightarrowl^{+}l^{-}"]->SetLineColor(8);      
    m_MCs["Z/#gamma*#rightarrowl^{+}l^{-}"]->SetFillColor(8);

    m_MCs["W+Jets"]->SetLineColor(4);      
    m_MCs["W+Jets"]->SetFillColor(4);

    m_MCs["Single-Top"]->SetLineColor(28);      
    m_MCs["Single-Top"]->SetFillColor(28);
  }

  unordered_map<string, TH1F*> m_sigs;
  for (int i=0,n=sigFileNames.size(); i<n; i++) {

    TFile* sigFile = TFile::Open(sigFileNames[i]);
    TH1F* h_sig = (TH1F*) sigFile->FindObjectAny(hname);
    h_sig->Scale(sigScales[i]);
    //h_sig->Sumw2();

    string key="";
    if ( sigFileNames[i].Contains("zprime", TString::kIgnoreCase) )     key = "Z' 3 TeV (x300)";
    else if ( sigFileNames[i].Contains("gluon", TString::kIgnoreCase) ) key = "g_{kk} 3 TeV";
    m_sigs[key] = h_sig;
  }

  if ( !m_sigs.empty() ) {
    m_sigs["Z' 3 TeV (x300)"]->SetLineColor(30);
    //m_sigs["Z' 3 TeV (x300)"]->SetFillColor(9);

    m_sigs["g_{kk} 3 TeV"]->SetLineColor(6);      
    //m_sigs["g_{kk} 3 TeV"]->SetFillColor(12);
  }

  TCanvas* c = new TCanvas("c", "c", 600, 600);
  if (logx) c->SetLogx();

  float b_scale = 0.3, t_scale = 1 - b_scale;
  TPad* top = new TPad("top", "top", 0, b_scale, 1, 1);
  TPad* bottom = new TPad("bottom", "bottom", 0, 0, 1, b_scale);

  if ( subplot=="ratio" || subplot=="diff" ) {
    top->SetTopMargin(0.05);
    top->SetBottomMargin(0.05);
    top->Draw();
    bottom->SetBottomMargin(0.35);
    bottom->SetGridy();
    bottom->Draw();
    top->cd();
  }
  TH1F* hist = new TH1F("hist", "hist", h_Data->GetNbinsX(), h_Data->GetBinLowEdge(1), h_Data->GetBinLowEdge(h_Data->GetNbinsX()+1));

  TString xtitle = "";
  string keytitle = hname(2, hname.Length()).Data();
  unordered_map<string, string> xtitles = {{"dilepmass","M_{ll} (Gev)"},{"lep0eta","#eta_{Leading Lepton}"},
  {"lep1eta","#eta_{Subleading Lepton}"},{"lep0pt","Leading Lepton p_{T} (GeV)"},{"lep1pt","Subleading Lepton p_{T} (GeV)"},
  {"jet0eta","#eta_{Leading Jet}"},{"jet1eta","#eta_{Subleading Jet}"},{"jet0pt","Leading Jet p_{T} (GeV)"},
  {"jet1pt","Subleading Jet p_{T} (GeV)"},{"nEle","Number of Electrons"},{"nMuon","Number of Muons"},
  {"nJet","Number of Jets"},{"nEleDiff","N_{Electrons}-N_{Good Electrons}"},{"nMuonDiff","N_{Muons}-N_{Good Muons} "},
  {"nJetDiff","N_{Jets}-N_{Good Jets} "},{"nGoodEle","N_{Good Electrons}"},{"nGoodMuon","N_{Good Muons}"},
  {"nGoodJet","N_{Good Jets}"},{"jethT","H_{T} (GeV)"},{"sT","H_{T}^{L} (GeV)"},{"sT_met","S_{T} (GeV)"},{"metpt","MET p_{T} (GeV)"},
  {"jet0btag","btag_{Leading Jet}"},{"jet1btag","btag_{Subeading Jet}"},{"nbtag","Number of btagged Jets"},
  {"metcorrpt","Corrected MET p_{T} (GeV)"},{"muonD0","Muon D_{0} (cm)"},{"muonDz","Muon D_{z} (cm)"},{"rmin0","#DeltaR_{min}(leading lepton, jet)"},
  {"rmin1","#DeltaR_{min}(subleading lepton, jet)"},{"lep0perp","Leading Lepton p_{T}^{rel} (GeV)"},{"lep1perp","Subleading Lepton p_{T}^{rel} (GeV)"},
  {"rl0cleanj","#DeltaR(leading lepton, cleaned jet)"},{"rl1cleanj","#DeltaR(subleading lepton, cleaned jet)"},{"rl0l1","#DeltaR(leading lepton, subleading lepton)"},
  {"jet0phi","#phi_{Leading Jet}"}, {"jet1phi","#phi_{Subleading Jet}"}, {"lep0phi","#phi_{Leading Lepton}"}, {"lep1phi","#phi_{Subleading Lepton}"},
  {"lepept","electron p_{T} (GeV)"}, {"lepmpt","muon p_{T} (GeV)"},{"rbl","#DeltaR(b quark, lepton)"},{"minjet0pt","Jet p_{T}^{rmin leading lepton} (GeV)"},
  {"minjet1pt","Jet p_{T}^{rmin subleading lepton} (GeV)"},{"cleanjet0pt","Jet p_{T}^{cleaned from leading lepton} (GeV)"},
  {"cleanjet1pt","Jet p_{T}^{cleaned from subleading lepton} (GeV)"},{"masslmin0","M_{leading lep,rmin jet} (Gev)"},{"masslmin1","M_{subleading lep,rmin jet} (Gev)"},
  {"masslljjm","M_{lljjmet} (Gev)"}};
  if (xtitles.find(keytitle) != xtitles.end()) xtitle = xtitles[keytitle];

  if ( subplot=="ratio" || subplot=="diff" ) {
    hist->GetXaxis()->SetTickLength(0.03/t_scale);
    hist->GetXaxis()->SetLabelSize(0);
    hist->GetYaxis()->SetTitleSize(0.04/t_scale);
    hist->GetYaxis()->SetTitleOffset(0.95);
    hist->GetYaxis()->SetLabelSize(0.03/t_scale);
  }
  else {
    hist->GetXaxis()->SetTitle(xtitle);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetXaxis()->SetNdivisions(5, 5, 0);
  }

  hist->GetXaxis()->SetRangeUser(xmin, xmax);
  //hist->GetXaxis()->SetNoExponent(false);
  //hist->GetXaxis()->SetMoreLogLabels();
  //hist->GetYaxis()->SetTitle("Events");
  hist->GetYaxis()->SetRangeUser(ymin, int(h_Data->GetMaximum()*1.2) );
  //hist->GetYaxis()->SetRangeUser(ymin, ymax);
  hist->Draw();
  mcStack->Draw("samehist");
  h_Data->Draw("samePE");

  for (unordered_map<string, TH1F*>::const_iterator i_sig=m_sigs.begin(); i_sig != m_sigs.end(); ++i_sig)
    i_sig->second->Draw("samehist");

  TLegend* leg = new TLegend(.7,.9-.05*(1.+m_MCs.size()+m_sigs.size()),.9,.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_Data, Form( "#bf{%s}", dataName.data() ), "P");
  for (unordered_map<string, TH1F*>::const_iterator i_MC=m_MCs.begin(); i_MC != m_MCs.end(); ++i_MC)
    leg->AddEntry(i_MC->second, Form( "#bf{%s}", i_MC->first.data() ), "L");
  for (unordered_map<string, TH1F*>::const_iterator i_sig=m_sigs.begin(); i_sig != m_sigs.end(); ++i_sig)
    leg->AddEntry(i_sig->second, Form( "#bf{%s}", i_sig->first.data() ), "L");

  leg->Draw();

  TLatex text;
  text.SetNDC();

  text.SetTextSize(0.035);
  text.SetTextFont(42);
  if (subplot=="ratio" || subplot=="diff" ) text.DrawLatex(1-rightText.Length()/90., 0.96, rightText);
  else text.DrawLatex(1-rightText.Length()/68., 0.96, rightText);

  text.SetTextSize(0.05);
  text.SetTextFont(61);
  text.DrawLatex(0.18, 0.96, "CMS");
  //text.DrawLatex(0.18, 0.96, leftText);
  text.SetTextSize(0.03);

  //text.SetTextSize(0.04);
  //text.SetTextFont(52);
  //text.DrawLatex(0.29, 0.96, "Simulation Preliminary"); //make bool

  text.SetTextSize(0.04);
  text.SetTextFont(42);
  if (dataFileName.Contains("mm", TString::kIgnoreCase)) text.DrawLatex(0.22,0.85,"#bf{#mu#mu}");
  else if (dataFileName.Contains("ee", TString::kIgnoreCase)) text.DrawLatex(0.22,0.85,"#bf{ee}");
  else text.DrawLatex(0.22,0.85,"#bf{e#mu}");

  if (hname.Contains("0_")) { text.DrawLatex(0.22,0.8,"#geq 0 btags"); text.DrawLatex(0.22,0.75,"p_{T}^{Leading Jet}>100 GeV"); }
  else if (hname.Contains("1_")) { text.DrawLatex(0.22,0.8,"#geq 1 btag"); text.DrawLatex(0.22,0.75,"p_{T}^{Leading Jet}>100 GeV"); }
  else if (hname.Contains("2_")) { text.DrawLatex(0.22,0.8,"#geq 0 btags"); text.DrawLatex(0.22,0.75,"p_{T}^{Leading Jet}>100 GeV, p_{T}^{Subleading Jet}>50 GeV"); }
  else if (hname.Contains("3_")) { text.DrawLatex(0.22,0.8,"#geq 1 btag"); text.DrawLatex(0.22,0.75,"p_{T}^{Leading Jet}>100 GeV, p_{T}^{Subleading Jet}>50 GeV"); }
  else if (hname.Contains("4_")) { text.DrawLatex(0.22,0.8,"#geq 2 btags"); text.DrawLatex(0.22,0.75,"p_{T}^{Leading Jet}>100 GeV, p_{T}^{Subleading Jet}>50 GeV"); }

  if (subplot=="ratio" || subplot=="diff" ) {
    bottom->cd();
    TH1F* bhist = new TH1F("bhist", "bhist", h_Data->GetNbinsX(), h_Data->GetBinLowEdge(1), h_Data->GetBinLowEdge(h_Data->GetNbinsX()+1));

    TH1F* hsubplot = (TH1F*) h_Data->Clone("hsubplot");
    TH1F* h_MC = new TH1F("h_MC", "h_MC", h_Data->GetNbinsX(), h_Data->GetBinLowEdge(1), h_Data->GetBinLowEdge(h_Data->GetNbinsX()+1));

    for (unordered_map<string, TH1F*>::const_iterator i_MC=m_MCs.begin(); i_MC != m_MCs.end(); ++i_MC)
      h_MC->Add(i_MC->second);

    TString subytitle = "Data/MC";
    if (subplot=="ratio") hsubplot->Divide(h_MC);
    else {
      subytitle = "(Data-MC)/#sigma";
      hsubplot->Add(h_MC, -1);
//      for (int i=1; i<hsubplot->GetNbinsX(); i++)
//        hsubplot->SetBinContent( i, hsubplot->GetBinContent(i) / hsubplot->GetBinError(i) );
    }

    bhist->GetXaxis()->SetNdivisions(5, 5, 0);
    bhist->GetXaxis()->SetLabelSize(0.02/b_scale);
    bhist->GetXaxis()->SetTickLength(0.03/b_scale);
    bhist->GetXaxis()->SetTitleSize(0.04/b_scale);
    bhist->GetXaxis()->SetTitleOffset(0.75);
    bhist->GetXaxis()->SetTitle(xtitle);
    bhist->GetXaxis()->SetRangeUser(xmin, xmax);
    bhist->GetYaxis()->SetRangeUser(subymin, subymax);
    bhist->GetYaxis()->SetNdivisions(5, 3, 0);
    bhist->GetYaxis()->SetLabelSize(0.03/b_scale);
    bhist->GetYaxis()->SetTitle(subytitle);
    bhist->GetYaxis()->SetTitleSize(0.035/b_scale);
    bhist->GetYaxis()->SetTitleOffset(0.43);

    bhist->Draw();
    hsubplot->Draw("sameP");
  }

  c->Print("./plots/" + outName + ".pdf");
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

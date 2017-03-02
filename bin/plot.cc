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
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void setPars(const string& parFile);
void setStyle();
void plotHist(const TString& hname);

//parameters- edit in plot_pars.txt
vector<TString> mcFileNames, mcNames;
vector<int> mcScales;
string subplot;
TString dataFileName, dataName;
TString hname, outName, leftText, rightText;
float xmax, ymax;
bool logx, logy;

int main(int argc, char* argv[]) {

  string parFile = argv[1];
  setPars(parFile);

  setStyle();

  TFile* dataFile = TFile::Open(dataFileName);

  //TString className = dataFile->Get(hname)->ClassName();

  TH1F* h_Data = (TH1F*) dataFile->Get(hname);
  //h_Data->Scale( 1 / h_Data->Integral() );  //make this bool, as well as exponents?

  h_Data->SetMarkerStyle(20);
  h_Data->SetLineColor(kBlack);
  h_Data->SetMarkerSize(0.65);

  vector<TH1F*> h_MCs;
  THStack* mcStack = new THStack();

  for (int i=0,n=mcFileNames.size(); i<n; i++) {

    TFile* mcFile = TFile::Open(mcFileNames[i]);
    h_MCs.push_back( (TH1F*) mcFile->Get(hname) );
    mcStack->Add( h_MCs.back() );
    h_MCs.back()->Scale(mcScales[i]);
    //h_MCs.back()->Sumw2();
    h_MCs.back()->SetLineColor(i+2);
    h_MCs.back()->SetFillColor(i+2);
    //h_MCs.back()->Scale( 1 / h_MCs.back()->Integral() );
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
  if (hname.EqualTo("jet0pt")) xtitle = "Leading Jet p_{T} (Gev)";
  else if (hname.EqualTo("dilepmass")) xtitle = "M_{ll} (Gev)";

  if ( subplot=="ratio" || subplot=="diff" ) {
    hist->GetXaxis()->SetTickLength(0.03/t_scale);
    hist->GetXaxis()->SetLabelSize(0);
    hist->GetYaxis()->SetTitleSize(0.06/t_scale);
    hist->GetYaxis()->SetTitleOffset(1);
    hist->GetYaxis()->SetLabelSize(0.05/t_scale);
  }
  else {
    hist->GetXaxis()->SetTitle(xtitle);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetXaxis()->SetNdivisions(5, 5, 0);
  }

  hist->GetXaxis()->SetRangeUser(0, xmax);
  //hist->GetXaxis()->SetNoExponent(false);
  //hist->GetXaxis()->SetMoreLogLabels();
  hist->GetYaxis()->SetTitle("Events");
  hist->GetYaxis()->SetRangeUser(0, ymax);
  hist->Draw();
  mcStack->Draw("samehist");
  h_Data->Draw("sameP");

  TLegend* leg = new TLegend(.7,.9-.05*(1.+mcNames.size()),.9,.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_Data, Form( "#bf{%s}", dataName.Data() ), "P");
  for (int i=0, n=mcNames.size(); i<n; i++) leg->AddEntry(h_MCs[i], Form( "#bf{%s}", mcNames[i].Data() ), "L");

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

  //text.SetTextSize(0.04);
  //text.SetTextFont(52);
  //text.DrawLatex(0.29, 0.96, "Simulation Preliminary"); //make bool

  if (subplot=="ratio" || subplot=="diff" ) {
    bottom->cd();
    TH1F* bhist = new TH1F("bhist", "bhist", h_Data->GetNbinsX(), h_Data->GetBinLowEdge(1), h_Data->GetBinLowEdge(h_Data->GetNbinsX()+1));

    TH1F* hsubplot = (TH1F*) h_Data->Clone("hsubplot");
    TH1F* h_MC = (TH1F*) h_MCs[0]->Clone("h_MC");
    for (int i=1, n=h_MCs.size(); i<n; i++) h_MC->Add(h_MCs[i]);

    TString subytitle = "Data/MC";
    if (subplot=="ratio") hsubplot->Divide(h_MC);
    else {
      subytitle = "(Data-MC)/#sigma";
      hsubplot->Add(h_MC, -1);
//      for (int i=1; i<hsubplot->GetNbinsX(); i++)
//        hsubplot->SetBinContent( i, hsubplot->GetBinContent(i) / hsubplot->GetBinError(i) );
    }

    int subymax = hsubplot->GetMaximum()+0.5;

    bhist->GetXaxis()->SetNdivisions(5, 5, 0);
    bhist->GetXaxis()->SetLabelSize(0.05/b_scale);
    bhist->GetXaxis()->SetTickLength(0.03/b_scale);
    bhist->GetXaxis()->SetTitleSize(0.06/b_scale);
    bhist->GetXaxis()->SetTitleOffset(0.75);
    bhist->GetXaxis()->SetTitle(xtitle);
    bhist->GetYaxis()->SetRangeUser(0, subymax);
    bhist->GetYaxis()->SetNdivisions(5, 3, 0);
    bhist->GetYaxis()->SetLabelSize(0.05/b_scale);
    bhist->GetYaxis()->SetTitle(subytitle);
    bhist->GetYaxis()->SetTitleSize(0.055/b_scale);
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
    else if (var == "dataName")  dataName = line.data();
    else if (var == "mcFileNames") {
      while ( (delim_pos = line.find(' ')) != -1) {
        mcFileNames.push_back( line.substr(0, delim_pos).data() );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      mcFileNames.push_back( line.data() );
    }
    else if (var == "mcNames") {
      while ( (delim_pos = line.find(' ')) != -1) {
        mcNames.push_back( line.substr(0, delim_pos).data() );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      mcNames.push_back( line.data() );
    }
    else if (var == "mcScales") {
      while ( (delim_pos = line.find(' ')) != -1) {
        mcScales.push_back( stoi( line.substr(0, delim_pos) ) );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      mcScales.push_back( stoi(line) );
    }
    else if (var == "outName")   outName = line.data();
    else if (var == "hname")     hname = line.data();
    else if (var == "leftText")  leftText = line.data();
    else if (var == "rightText") rightText = line.data();
    else if (var == "xmax")      xmax = stof(line);
    else if (var == "ymax")      ymax = stof(line);
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

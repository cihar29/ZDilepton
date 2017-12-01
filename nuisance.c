// Chad Harrington 10/31/2017 - create nuisance/correlation plots

void setStyle();
void readFile(const string& dir, const string& parFile, vector< vector<string> >& vec);

TString rightText = "Run 2016 - 35.9 fb^{-1} (13 TeV)";

void nuisance( string folder = "all_st/gkk", string filename = "gkk3000.txt", bool asimov = false ) {
  string dir = "/uscms_data/d3/cihar29/Analysis/CMSSW_8_1_0/src/theta/utils2/2017/";
  if (asimov) filename =  "asv_" + filename;

  vector< vector<string> > v_pars;
  readFile( dir+folder+"/", "nuisances_"+filename, v_pars );

  setStyle();

  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.25);

  TCanvas* c = new TCanvas("c", "c", 600, 900);
  TGraphErrors* g_band0 = new TGraphErrors(), *g_band1 = new TGraphErrors(), *g_pars = new TGraphErrors();

  int nPars = v_pars.size();
  for (int i=0; i<nPars; i++) {
    g_band0->SetPoint(i, 0, i+1);
    g_band1->SetPoint(i, 0, i+1);
    g_pars->SetPoint(i, stof(v_pars[i][1]), i+1);

    g_band0->SetPointError(i, 2, 1);
    g_band1->SetPointError(i, 1, 1);
    g_pars->SetPointError(i, stof(v_pars[i][2]), 0);
  }

  g_band0->SetFillColor(kYellow+1);
  g_band1->SetFillColor(kGreen+1);

  g_pars->SetMarkerStyle(20);
  g_pars->SetMarkerColor(kBlack);
  g_pars->SetLineColor(kBlack);
  g_pars->SetMarkerSize(0.65);

  g_band0->GetXaxis()->SetNdivisions(5, 5, 0);
  g_band0->GetXaxis()->SetRangeUser(-2.5, 2.5);
  g_band0->GetXaxis()->SetTitle("Post-Fit Values");
  g_band0->GetXaxis()->SetLabelSize(0.04);
  g_band0->GetXaxis()->SetTitleSize(0.05);
  g_band0->GetXaxis()->SetTitleOffset(0.85);

  g_band0->GetYaxis()->SetNdivisions(nPars+1, 0, 0);
  g_band0->GetYaxis()->SetRangeUser(0, nPars+1);
  g_band0->GetYaxis()->SetLabelSize(0);

  g_band0->Draw("AE2");
  g_band1->Draw("sameE2");
  g_pars->Draw("Psame");

  gPad->RedrawAxis();

  TLatex text;
  text.SetNDC();

  text.SetTextSize(0.05);
  text.SetTextFont(61);
  text.DrawLatex(0.26, 0.96, "CMS");

  text.SetTextSize(0.035);
  text.SetTextFont(42);
  text.DrawLatex(1-rightText.Length()/72., 0.96, rightText);

  for (int i=0; i<nPars; i++) text.DrawLatex(0.01, 0.095+(i+1)*0.00225*nPars, Form("#bf{%s}", v_pars[i][0].data()));

  c->Print( ( folder + (asimov ? "/nuisance_asv.pdf" : "/nuisance.pdf") ).data() );

  /// Correlation Matrix ///

  vector< vector<string> > v_cors;
  readFile( dir+folder+"/", "covariance_"+filename, v_cors );

  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPaintTextFormat("4.2f");

  TCanvas* c2 = new TCanvas("c2", "c2", 900, 600);
  TH2F* h = new TH2F("h", "h", nPars, 0, nPars, nPars, 0, nPars);

  for (int i=0; i<nPars; i++) {
    h->GetXaxis()->SetBinLabel(i+1, v_cors[i][0].data());
    h->GetYaxis()->SetBinLabel(i+1, v_cors[i][0].data());

    for (int j=0; j<nPars; j++) {
      double cor = stof(v_cors[i][j+1]) / sqrt( stof(v_cors[i][i+1]) * stof(v_cors[j][j+1]) );
      cor = cor==0 ? 0.000001 : cor;
      h->SetBinContent(i+1, j+1, cor);
    }
  }

  h->GetXaxis()->SetTickLength(0);
  h->GetXaxis()->SetLabelSize(0.04);

  h->GetYaxis()->SetTickLength(0);
  h->GetYaxis()->SetLabelSize(0.04);

  h->GetZaxis()->SetRangeUser(-1, 1);
  h->GetZaxis()->SetLabelSize(0.03);
  h->Draw("textcolz");

  text.SetTextSize(0.05);
  text.SetTextFont(61);
  text.DrawLatex(0.15, 0.96, "CMS");

  text.SetTextSize(0.035);
  text.SetTextFont(42);
  text.DrawLatex(1-rightText.Length()/85., 0.96, rightText);

  c2->Print( ( folder + (asimov ? "/correlation_asv.pdf" : "/correlation.pdf") ).data() );
}

void readFile(const string& dir, const string& parFile, vector< vector<string> >& vec) {

  ifstream file(dir + parFile);
  string line;

  while (getline(file, line)) {

    while (line.length() > 0 && line.at(0) == ' ') line.erase(0, 1);
    if (line.length() > 0 && line.at(0) == '#') continue;

    vector<string> subvec;
    int delim_pos;
    while ( (delim_pos = line.find(' ')) != -1) {
      subvec.push_back( line.substr(0, delim_pos) );

      line.erase(0, delim_pos+1);
      while (line.length() > 0 && line.at(0) == ' ') line.erase(0, 1);
    }
    if (line.length() > 0) subvec.push_back( line.substr(0, line.length()) );
    vec.push_back(subvec);
  }
}

void setStyle() {

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

//End Style//
}

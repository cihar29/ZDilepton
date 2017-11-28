// Chad Harrington 10/8/2017 - create brazilian plots

void setStyle();
void readFile(const string& dir, const string& parFile, vector< vector<double> >& vec);
double intersection(const TGraph* g1, const TGraph* g2, double start, double end);

TString rightText = "Run 2016 - 35.9 fb^{-1} (13 TeV)";

void brazilian( string folder = "all_st/gkk", bool plotObs=false) {
  string dir = "/uscms_data/d3/cihar29/Analysis/CMSSW_8_1_0/src/theta/utils2/2017/";

  vector< vector<double> > v_exp, v_obs, v_theory;
  readFile( dir+folder+"/", "bayesian_limits_expected.txt", v_exp );
  readFile( dir+folder+"/", "bayesian_limits_observed.txt", v_obs );
  readFile( "", "theory.txt", v_theory );

  TString signame = "", tfolder = folder.data(), ytitle="#sigma(Z'#rightarrowt#bar{t}) (pb)";
  int sig = -1;
  if      (tfolder.Contains("zp10", TString::kIgnoreCase)) { signame = "Z' (10% width)"; sig = 3; }
  else if (tfolder.Contains("zp1", TString::kIgnoreCase))  { signame = "Z' (1% width)";  sig = 2; }
  else if (tfolder.Contains("zp30", TString::kIgnoreCase)) { signame = "Z' (30% width)"; sig = 4; }
  else if (tfolder.Contains("gkk", TString::kIgnoreCase))  { signame = "g_{kk}";         sig = 1; ytitle = "#sigma(g_{kk}#rightarrowt#bar{t}) (pb)";}

  setStyle();

  TCanvas* c = new TCanvas("c", "c", 600, 600);
  gPad->SetLogy();
  TGraphAsymmErrors* g_band0 = new TGraphAsymmErrors(), *g_band1 = new TGraphAsymmErrors();
  TGraphErrors* g_obs = new TGraphErrors();
  TGraph* g_exp = new TGraph(), *g_theory = new TGraph();

  for (int i=0,n=v_exp.size(); i<n; i++) {
    g_exp->SetPoint(i, v_exp[i][0], v_exp[i][1]);
    g_band0->SetPoint(i, v_exp[i][0], 0);
    g_band1->SetPoint(i, v_exp[i][0], 0);
    g_obs->SetPoint(i, v_obs[i][0], v_obs[i][1]);

    g_band0->SetPointEYlow( i, -1*v_exp[i][2] );
    g_band0->SetPointEYhigh( i, v_exp[i][3] );
    g_band1->SetPointEYlow( i, -1*v_exp[i][4] );
    g_band1->SetPointEYhigh( i, v_exp[i][5] );

    g_obs->SetPointError( i, 0, v_obs[i][2] );
  }
  for (int i=0,n=v_theory.size(); i<n; i++) {
    g_theory->SetPoint(i, v_theory[i][0], v_theory[i][sig]);
  }

  g_band0->SetFillColor(kYellow+1);
  g_band1->SetFillColor(kGreen+1);
  g_exp->SetLineColor(kBlack);
  g_exp->SetLineWidth(2);
  g_exp->SetLineStyle(2);

  g_obs->SetMarkerStyle(20);
  g_obs->SetMarkerColor(kGray+2);
  g_obs->SetLineColor(kGray+2);
  g_obs->SetMarkerSize(0.65);
  g_theory->SetLineWidth(3);
  g_theory->SetLineColor(kRed);

  int legEntries = plotObs ? 5 : 4;
  TLegend* leg = new TLegend(.6,.9-.06*legEntries,.85,.9);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetHeader("95% CL upper limits");

  leg->AddEntry(g_exp,   "#bf{Expected}", "L");
  leg->AddEntry(g_band1, "#bf{68% expected}", "F");
  leg->AddEntry(g_band0, "#bf{95% expected}", "F");
  if (plotObs) leg->AddEntry(g_obs,   "#bf{Observed}", "PLE");
  leg->AddEntry(g_theory, Form( "#bf{%s}", signame.Data() ), "L");

  g_band0->GetXaxis()->SetNdivisions(5, 5, 0);
  g_band0->GetXaxis()->SetRangeUser(500, 5000);
  g_band0->GetXaxis()->SetTitle("Mass (GeV)");
  g_band0->GetXaxis()->SetLabelSize(0.04);
  g_band0->GetXaxis()->SetTitleSize(0.05);
  g_band0->GetXaxis()->SetTitleOffset(1.1);

  g_band0->GetYaxis()->SetRangeUser(0.01, 100);
  g_band0->GetYaxis()->SetTitle(ytitle);
  g_band0->GetYaxis()->SetLabelSize(0.04);
  g_band0->GetYaxis()->SetTitleSize(0.05);
  g_band0->GetYaxis()->SetTitleOffset(1.3);

  g_band0->Draw("AE3");
  g_band1->Draw("sameE3");
  g_exp->Draw("Lsame");
  if (plotObs) g_obs->Draw("PLsame");
  g_theory->Draw("Lsame");

  leg->Draw();
  gPad->RedrawAxis();

  TLatex text;
  text.SetNDC();

  text.SetTextSize(0.05);
  text.SetTextFont(61);
  text.DrawLatex(0.18, 0.96, "CMS");

  text.SetTextSize(0.035);
  text.SetTextFont(42);
  text.DrawLatex(1-rightText.Length()/70., 0.96, rightText);

  //text.DrawLatex(0.2,0.87,"#bf{#geq 1 btag}");
  //text.DrawLatex(0.2,0.82,"#bf{p_{T}^{j0}>100 GeV, p_{T}^{j1}>50 GeV}");

  c->Print( tfolder + (plotObs ? "/brazilian_obs.pdf" : "/brazilian.pdf") );

  double minx = g_exp->GetX()[0];
  double maxx = g_exp->GetX()[g_exp->GetN()-1];

  cout << "Expected: " << intersection(g_theory, g_exp, minx, maxx) << endl;
  cout << "Observed: " << intersection(g_theory, g_obs, minx, maxx) << endl;
}

//linear intersection
double intersection(const TGraph* g1, const TGraph* g2, double start, double end) {
  int n = 500;
  double width, min_diff = 10., intersect = 0.;

  do {
    width = (end-start)/n;
    for (int i=0; i<n; i++) {
      double x = start + i*width;
      double diff = fabs(g1->Eval(x) - g2->Eval(x));

      if (diff < min_diff) { min_diff = diff; intersect = x; }
    }
    start = intersect-10; end = intersect+10;
  } while (width > 0.1);

  return intersect;
}

void readFile(const string& dir, const string& parFile, vector< vector<double> >& vec) {

  ifstream file(dir + parFile);
  string line;

  while (getline(file, line)) {

    while (line.length() > 0 && line.at(0) == ' ') line.erase(0, 1);
    if (line.length() > 0 && line.at(0) == '#') continue;

    vector<double> subvec;
    int delim_pos;
    while ( (delim_pos = line.find(' ')) != -1) {
      subvec.push_back( stof( line.substr(0, delim_pos) ) );

      line.erase(0, delim_pos+1);
      while (line.length() > 0 && line.at(0) == ' ') line.erase(0, 1);
    }
    if (line.length() > 0) subvec.push_back( stof( line.substr(0, line.length()) ) );
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

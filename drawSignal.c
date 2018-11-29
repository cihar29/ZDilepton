// Chad Harrington 12/8/2017 - Draw signal masses

void setStyle();

void drawSignal( TString sigstring = "gkk" ) {

  TString signame = "";
  vector<TString> files;

  if      (sigstring == "zp1") {
    signame = "Z' (1% width)";
    files = { "zprime_M-1000_W-10.root", "zprime_M-2000_W-20.root", "zprime_M-3000_W-30.root", "zprime_M-4000_W-40.root", "zprime_M-4500_W-45.root", "zprime_M-5000_W-50.root" };
  }
  else if (sigstring == "zp10") {
    signame = "Z' (10% width)";
    files = { "zprime_M-1000_W-100.root", "zprime_M-2000_W-200.root", "zprime_M-3000_W-300.root", "zprime_M-4000_W-400.root", "zprime_M-4500_W-450.root", "zprime_M-5000_W-500.root" };
  }
  else if (sigstring == "zp30") {
    signame = "Z' (30% width)";
    files = { "zprime_M-1000_W-300.root", "zprime_M-2000_W-600.root", "zprime_M-3000_W-900.root", "zprime_M-4000_W-1200.root", "zprime_M-5000_W-1500.root" };
  }
  else if (sigstring == "gkk") {
    signame = "g_{kk}";
    files = { "gluon_M-1000.root", "gluon_M-2000.root", "gluon_M-3000.root", "gluon_M-4000.root", "gluon_M-4500.root", "gluon_M-5000.root" };
  }

  TString dir = "root_trees/";
  int nFiles = files.size();

  setStyle();
  TCanvas* c = new TCanvas("c", "c", 600, 600);

  TLegend* leg = new TLegend(.6,.9-.05*nFiles,.75,.9);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetHeader(signame);

  TH1F* h = new TH1F("h", "h", 300, 0, 6000);
  h->GetXaxis()->SetNdivisions(6, 6, 0);
  h->GetXaxis()->SetTitle("Mass (GeV)");
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(1.1);

  h->GetYaxis()->SetLabelSize(0.04);
  h->Draw();

  float max=0;
  for (int i=0; i<nFiles; i++) {

    TFile* file = TFile::Open(dir + files[i]);
    TTree* tree = (TTree*) file->Get("T");

    TH1F* h_mass = new TH1F("mass", "mass", 300, 0, 6000);
    tree->Draw("gen_mass>>mass", "gen_PID>1000000 && gen_status==22", "histsame");

    if (sigstring == "zp30" && i==4) h_mass->SetLineColor(i+2);
    else                             h_mass->SetLineColor(i+1);
    h_mass->Scale( 1 / h_mass->Integral() );
    if (i==0) max = h_mass->GetMaximum();

    int index = sigstring == "gkk" ? files[i].Last('.') : files[i].Last('_');
    TString mass = files[i](files[i].Index("M-")+2, index-files[i].Index("M-")-2);
    leg->AddEntry(h_mass, Form( "#bf{%s (RMS %.0f)}", mass.Data(), h_mass->GetRMS() ), "L");

  }
  leg->Draw();
  h->GetYaxis()->SetRangeUser(0, max*1.2);

  TLatex text;
  text.SetNDC();

  text.SetTextSize(0.05);
  text.SetTextFont(61);
  text.DrawLatex(0.18, 0.96, "CMS");

  c->Print(sigstring + ".pdf");
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

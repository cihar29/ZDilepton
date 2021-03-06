// Chad Harrington 5/11/2017 - create histograms for pileup reweighting

void setStyle();
void drawText();

void pileup_weights() {

  TString dir = "/uscms_data/d3/cihar29/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/root_trees/";
  TString files[] = {
    "mu_DataUP",
    "mu_DataNOMINAL",
    "mu_DataDOWN",
    "DYhigh",
    "DYlow",
    "STschannel",
    "STtWchannel",
    "STtchannel",
    "SaTtWchannel",
    "SaTtchannel",
    "TTbar0-700",
    "TTbar700-1000",
    "TTbar1000-inf",
    "WJets",
    "WW",
    "WZ",
    "ZZ",
    "gluon_M-1000",
    "gluon_M-1250",
    "gluon_M-1500",
    "gluon_M-2000",
    "gluon_M-2500",
    "gluon_M-3000",
    "gluon_M-3500",
    "gluon_M-4000",
    "gluon_M-4500",
    "gluon_M-5000",
    "gluon_M-500",
    "gluon_M-750",
    "zprime_M-1000_W-10",
    "zprime_M-1000_W-100",
    "zprime_M-1000_W-300",
    "zprime_M-1250_W-125",
    "zprime_M-1250_W-12p5",
    "zprime_M-1500_W-15",
    "zprime_M-1500_W-150",
    "zprime_M-3000_W-300",
    "zprime_M-3000_W-900",
    "zprime_M-500_W-5",
    "zprime_M-500_W-50",
    "zprime_M-750_W-7p5",
    "zprime_M-750_W-75",
    "zprime_M-2000_W-20",
    "zprime_M-2000_W-200",
    "zprime_M-2000_W-600",
    "zprime_M-2500_W-25",
    "zprime_M-2500_W-250",
    "zprime_M-3000_W-30",
    "zprime_M-3500_W-35",
    "zprime_M-3500_W-350",
    "zprime_M-4000_W-40",
    "zprime_M-4000_W-400",
    "zprime_M-4000_W-1200",
    "zprime_M-4500_W-45",
    "zprime_M-4500_W-450",
    "zprime_M-5000_W-50",
    "zprime_M-5000_W-500",
    "zprime_M-5000_W-1500",
    "qcd30-50_Mu",
    "qcd50-80_Mu",
    "qcd80-120_Mu",
    "qcd120-170_Mu",
    "qcd170-300_Mu",
    "qcd300-470_Mu",
    "qcd470-600_Mu",
    "qcd600-800_Mu",
    "qcd800-1000_Mu",
    "qcd1000-inf_Mu",
    "qcd30-50_EM",
    "qcd50-80_EM",
    "qcd80-120_EM",
    "qcd120-170_EM",
    "qcd170-300_EM",
    "qcd300-inf_EM",
    "qcd30-80_bcToE",
    "qcd80-170_bcToE",
    "qcd170-250_bcToE",
    "qcd250-inf_bcToE"
  };
  int file_size = sizeof(files)/sizeof(files[0]);

  setStyle();
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  TH1F* h = new TH1F("h", "h", 100, 0, 100);

  h->GetXaxis()->SetTitle("#mu");
  h->GetYaxis()->SetRangeUser(0, 0.07);
  h->Draw();

  TLegend* leg = new TLegend(.5,.9-file_size*0.03,.7,.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.025);
  leg->SetTextFont(42);

  map<TString, TH1F*> hists;
  for (int i=0; i<file_size; i++) {

    TString name = files[i];
    TFile* file = TFile::Open( dir + name + ".root" );

    hists[name] = (TH1F*) file->Get( "fullMu" );
    hists[name]->Scale( 1 / hists[name]->Integral() );

    if ( name.Contains("Data") ) hists[name]->SetLineStyle(2);
    hists[name]->SetLineColor(i+1);
    hists[name]->Draw("histsame");

    leg->AddEntry(hists[name],  Form( "%-15s (Mean %4.1f, RMS %4.1f)" , name.Data(), hists[name]->GetMean(), hists[name]->GetRMS() ), "L");
  }
  leg->Draw();
  drawText();
  c->Print("mu.pdf");

  c->Clear();
  leg->Clear();
  h->GetYaxis()->SetTitle("Weight");
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetRangeUser(0, 2.5);
  h->Draw();

  TString wgt_type[] = { "UP", "NOMINAL", "DOWN" };

  TFile* outFile = new TFile("mu_weights.root","RECREATE");
  outFile->cd();

  for (int i=0, n=sizeof(wgt_type)/sizeof(wgt_type[0]); i<n; i++) outFile->mkdir( wgt_type[i] + "/" );

  map<TString, TH1F*> weights;
  for (int i=0, n=sizeof(wgt_type)/sizeof(wgt_type[0]); i<n; i++) {
    outFile->cd( "mu_weights.root:/" + wgt_type[i] );

    for (map<TString, TH1F*>::iterator it = hists.begin(); it != hists.end(); it++) {

      TString name = it->first + "_" + wgt_type[i];
      if ( name.Contains("Data") ) continue;

      TH1F* mcHist = it->second;

      weights[name] = (TH1F*) hists[ "mu_Data" + wgt_type[i] ]->Clone(name);
      weights[name]->Divide( mcHist );

      weights[name]->SetLineColor(i+1);
      weights[name]->SetLineStyle(1);
      weights[name]->Draw("histsame");

      leg->AddEntry(weights[name], name, "L");
      weights[name]->Write(name);
    }
    //outFile->Write();
  }

  leg->SetX1NDC(0.6);
  leg->SetX2NDC(0.8);
  leg->Draw();
  drawText();
  c->Print("mu_weights.pdf");

  delete outFile;
  outFile = 0;
}

void drawText() {

  TLatex text;
  text.SetNDC();

  text.SetTextSize(0.035);
  text.SetTextFont(42);
  TString rightText = "Run 2016 - 35.9 fb^{-1} (13 TeV)";
  text.DrawLatex(1-rightText.Length()/68., 0.96, rightText);

  text.SetTextSize(0.05);
  text.SetTextFont(61);
  text.DrawLatex(0.18, 0.96, "CMS");
  text.SetTextSize(0.03);
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

//End Style//
}

// Chad Harrington and Bahareh Roozbahani 5/11/2017 - btagging efficiency

void setStyle();

void btag_efficiency(TString channel = "mm", TString flavor = "b") {

  if ( !channel.EqualTo("mm") && !channel.EqualTo("ee") ) channel = "em";
  if ( !flavor.EqualTo("b") && !flavor.EqualTo("c") ) flavor = "udsg";

  const int nBins = 15;
  const double bins[nBins+1] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1000, 1200, 2000};

  map<TString, TString> files;
  files["t#bar{t}"]             = "TTbar";
  //files["lowDY"]             = "lowDY";
  //files["highDY"]            = "highDY";
  //files["STtchannel"]        = "STtchannel";
  //files["SaTtchannel"]       = "SaTtchannel";
  //files["STschannel"]        = "STschannel";
  //files["STtWchannel"]       = "STtWchannel";
  //files["SaTtWchannel"]      = "STtWchannel";
  //files["Wjet"]              = "Wjet";
  //files["gluonkk-M3000"]     = "gluonkk-M3000";
  files["Z' 3 TeV (width=0.1)"] = "zprime-M3000-W300";

  setStyle();
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  c->SetGrid();
  TH1F* h = new TH1F("h", "h", 200, 0, 2000);

  h->GetXaxis()->SetTitle("Jet p_{T} (GeV)");
  h->GetXaxis()->SetRangeUser(0, 1200);
  h->GetYaxis()->SetTitle("b-tagging #varepsilon");
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetRangeUser(0, 1.5);
  h->Draw();

  TLegend* leg = new TLegend(.5,.9-files.size()*2*0.04,.7,.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);

  map<TString, TH1F*> hists;
  for (map<TString, TString>::iterator it = files.begin(); it != files.end(); it++){

    TFile* file = TFile::Open( it->second + "_" + channel + ".root" );
    TH1F* jetPt = (TH1F*) file->Get( "jetPt_" + flavor );
    jetPt->Sumw2();
    jetPt = (TH1F*) jetPt->Rebin(nBins, "jetPt", bins);

    int j = distance( files.begin(), it );

    TString name = "CSVL, " + it->first;
    hists[name] = (TH1F*) file->Get( "jetPt_bTagL_" + flavor );
    hists[name]->Sumw2();
    hists[name] = (TH1F*) hists[name]->Rebin(nBins, name, bins);
    hists[name]->Divide(jetPt);

    hists[name]->SetMarkerStyle(j+20);
    hists[name]->SetMarkerColor(j+1);
    hists[name]->SetLineColor(j+1);
    hists[name]->Draw("pesame");

    leg->AddEntry(hists[name], name, "PL");

    name = "CSVM, " + it->first;
    hists[name] = (TH1F*) file->Get( "jetPt_bTagM_" + flavor );
    hists[name]->Sumw2();
    hists[name] = (TH1F*) hists[name]->Rebin(nBins, name, bins);
    hists[name]->Divide(jetPt);

    hists[name]->SetMarkerStyle(j+24);
    hists[name]->SetMarkerColor(j+1);
    hists[name]->SetLineColor(j+1);
    hists[name]->Draw("pesame");

    leg->AddEntry(hists[name], name, "PL");
  }
  leg->Draw();

  TLatex text;
  text.SetNDC();

  text.SetTextSize(0.05);
  text.SetTextFont(61);
  text.DrawLatex(0.18, 0.96, "CMS");

  text.SetTextSize(0.04);
  text.DrawLatex(0.2, 0.87, flavor + "-Jets");
  text.DrawLatex(0.2, 0.83, channel + " channel");

  text.SetTextFont(52);
  text.DrawLatex(0.29, 0.96, "Simulation");

  text.SetTextSize(0.035);
  text.SetTextFont(42);
  text.DrawLatex(0.85, 0.96, "(13 TeV)");

  c->Print("btag_efficiency_" + flavor + "_" + channel + ".pdf");
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

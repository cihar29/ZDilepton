// Chad Harrington and Bahareh Roozbahani 5/11/2017 - btagging efficiency

void setStyle();

void btag_efficiency() {

  const int nBins = 15;
  const double bins[nBins+1] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1000, 1400, 2000};

  TString files[] = {
    "DYhigh",
    "DYlow",
    "STschannel",
    "STtWchannel",
    "STtchannel",
    "SaTtWchannel",
    "SaTtchannel",
    "TTbar",
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
    "gluon_M-500",
    "gluon_M-750",
    "zprime_M-1000_W-10",
    "zprime_M-1000_W-100",
    "zprime_M-1000_W-300",
    "zprime_M-1250_W-125",
    "zprime_M-1250_W-12p5",
    "zprime_M-1500_W-15",
    "zprime_M-1500_W-150",
    "zprime_M-3000_W-300"
  };
  int file_size = sizeof(files)/sizeof(files[0]);

  setStyle();
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  c->SetGrid();
  TH1F* h = new TH1F("h", "h", nBins, bins);

  h->GetXaxis()->SetTitle("Jet p_{T} (GeV)");
  h->GetXaxis()->SetRangeUser(0, 1400);
  h->GetYaxis()->SetTitle("b-tagging #varepsilon");
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetRangeUser(0, 1.4);

  TLegend* leg = new TLegend(.5,.9-2*0.04,.7,.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);

  TLatex text;
  text.SetNDC();

  map<TString, TH1F*> effs;

  TString channels[] = {"mm", "em", "ee"}, flavors[] = {"b", "c", "udsg"};
  int nChannels = sizeof(channels)/sizeof(channels[0]);
  int nFlavors = sizeof(flavors)/sizeof(flavors[0]);

  TString clabels[] = {"#mu#mu", "e#mu", "ee"};
  float ymax[] = {1.4, 1.4, 0.5};

  for (int i=0; i<nChannels; i++) {

    for (int j=0; j<nFlavors; j++) {

      TString LWP_num = flavors[j] + "_LWP_" + channels[i];
      effs[LWP_num] = new TH1F("LWP_num", "LWP_num", nBins, bins);
      TString MWP_num = flavors[j] + "_MWP_" + channels[i];
      effs[MWP_num] = new TH1F("MWP_num", "MWP_num", nBins, bins);
      TString denom = flavors[j] + "_denom_" + channels[i];
      effs[denom] = new TH1F("denom", "denom", nBins, bins);

      for (int k=0; k<file_size; k++) {

        if ( files[k].Contains("zprime") || files[k].Contains("gluon") ) continue;
        TFile* file = TFile::Open( channels[i] + "/" + files[k] + "_" + channels[i] + ".root" );

        TH1F* hist = (TH1F*) file->Get( "jetPt_" + flavors[j] );
        hist->Sumw2();
        hist = (TH1F*) hist->Rebin(nBins, "hist", bins);
        effs[denom]->Add(hist);

        hist = (TH1F*) file->Get( "jetPt_bTagL_" + flavors[j] );
        hist->Sumw2();
        hist = (TH1F*) hist->Rebin(nBins, "hist", bins);
        effs[LWP_num]->Add(hist);

        hist = (TH1F*) file->Get( "jetPt_bTagM_" + flavors[j] );
        hist->Sumw2();
        hist = (TH1F*) hist->Rebin(nBins, "hist", bins);
        effs[MWP_num]->Add(hist);
      }
      effs[LWP_num]->Divide(effs[denom]);
      effs[MWP_num]->Divide(effs[denom]);

      h->GetYaxis()->SetRangeUser(0, ymax[j]);
      h->Draw();

      effs[LWP_num]->SetMarkerStyle(20);
      effs[LWP_num]->SetMarkerColor(1);
      effs[LWP_num]->SetLineColor(1);
      effs[LWP_num]->Draw("pesame");

      effs[MWP_num]->SetMarkerStyle(24);
      effs[MWP_num]->SetMarkerColor(1);
      effs[MWP_num]->SetLineColor(1);
      effs[MWP_num]->Draw("pesame");

      leg->AddEntry(effs[LWP_num], "CSVL", "PL");
      leg->AddEntry(effs[MWP_num], "CSVM", "PL");
      leg->Draw();

      text.SetTextSize(0.05);  text.SetTextFont(61); text.DrawLatex(0.18, 0.96, "CMS");
      text.SetTextSize(0.04);  text.SetTextFont(52); text.DrawLatex(0.29, 0.96, "Simulation");
      text.SetTextSize(0.035); text.SetTextFont(42); text.DrawLatex(0.85, 0.96, "(13 TeV)");

      text.SetTextSize(0.04);  text.SetTextFont(62);
      text.DrawLatex(0.2, 0.87, flavors[j] + "-Jets");
      text.DrawLatex(0.2, 0.83, clabels[i] + " channel");

      c->Print("plots/btag_eff_" + flavors[j] + "_" + channels[i] + ".pdf");
      leg->Clear();
      c->Clear();
    }
  }

  TFile* outFile = new TFile("btag_eff.root","RECREATE");
  outFile->cd();

  for (int i=0; i<nChannels; i++) outFile->mkdir( channels[i] + "/" );

  for (map<TString, TH1F*>::iterator it = effs.begin(); it != effs.end(); it++) {
    outFile->cd();
    TString name = it->first;
    TString postfix = name(name.Last('_')+1, 2);

    outFile->cd("btag_eff.root:/" + postfix);

    if (!name.Contains("_denom_")) it->second->Write(name);
  }

  outFile->Write();
  delete outFile;
  outFile = 0;

  //////////////////////////////////////////////////
  //Plot individual efficiencies on same histogram//
  //////////////////////////////////////////////////

  TString filesToPlot[] = { "TTbar", "zprime_M-3000_W-300" };
  file_size = sizeof(filesToPlot)/sizeof(filesToPlot[0]);

  TString channel = "mm", flavor = "b";

  h->GetYaxis()->SetRangeUser(0, 1.4);
  h->Draw();

  leg = new TLegend(.5,.9-2*file_size*0.04,.7,.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);

  text.SetTextSize(0.05);  text.SetTextFont(61); text.DrawLatex(0.18, 0.96, "CMS");
  text.SetTextSize(0.04);  text.SetTextFont(52); text.DrawLatex(0.29, 0.96, "Simulation");
  text.SetTextSize(0.035); text.SetTextFont(42); text.DrawLatex(0.85, 0.96, "(13 TeV)");

  text.SetTextSize(0.04);  text.SetTextFont(62);
  text.DrawLatex(0.2, 0.87, flavor + "-Jets");
  text.DrawLatex(0.2, 0.83, clabels[0] + " channel");

  map<TString, TH1F*> hists;
  for (int i=0; i<file_size; i++) {

    TString dataset = filesToPlot[i];

    TFile* file = TFile::Open( channel + "/" + dataset + "_" + channel + ".root" );
    TH1F* jetPt = (TH1F*) file->Get( "jetPt_" + flavor );
    jetPt->Sumw2();
    jetPt = (TH1F*) jetPt->Rebin(nBins, "jetPt", bins);

    TString name = "CSVL,  " + dataset;
    hists[name] = (TH1F*) file->Get( "jetPt_bTagL_" + flavor );
    hists[name]->Sumw2();
    hists[name] = (TH1F*) hists[name]->Rebin(nBins, name, bins);
    hists[name]->Divide(jetPt);

    hists[name]->SetMarkerStyle(i+20);
    hists[name]->SetMarkerColor(i+1);
    hists[name]->SetLineColor(i+1);
    hists[name]->Draw("pesame");

    leg->AddEntry(hists[name], name, "PL");

    name = "CSVM, " + dataset;
    hists[name] = (TH1F*) file->Get( "jetPt_bTagM_" + flavor );
    hists[name]->Sumw2();
    hists[name] = (TH1F*) hists[name]->Rebin(nBins, name, bins);
    hists[name]->Divide(jetPt);

    hists[name]->SetMarkerStyle(i+24);
    hists[name]->SetMarkerColor(i+1);
    hists[name]->SetLineColor(i+1);
    hists[name]->Draw("pesame");

    leg->AddEntry(hists[name], name, "PL");
  }
  leg->Draw();
  c->Print("plots/btag_eff_" + flavor + "_" + channel + "_individual.pdf");
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

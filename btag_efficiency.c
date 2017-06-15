// Chad Harrington and Bahareh Roozbahani 5/11/2017 - btagging efficiency

void setStyle();

void btag_efficiency(TString channel = "mm", TString flavor = "b", bool createFile = false) {

  if ( !channel.EqualTo("mm") && !channel.EqualTo("ee") ) channel = "em";
  if ( !flavor.EqualTo("b") && !flavor.EqualTo("c") ) flavor = "udsg";

  const int nBins = 15;
  const double bins[nBins+1] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1000, 1400, 2000};

  map<TString, TString> files;
  files["t#bar{t}"]             = "TTbar";
  files["lowDY"]                = "lowDY";
  files["highDY"]               = "highDY";
  files["STtchannel"]           = "STtchannel";
  files["SaTtchannel"]          = "SaTtchannel";
  files["STschannel"]           = "STschannel";
  files["STtWchannel"]          = "STtWchannel";
  files["SaTtWchannel"]         = "STtWchannel";
  files["Wjet"]                 = "Wjet";
  files["gluonkk-M3000"]        = "gluonkk-M3000";
  files["Z' 3 TeV (width=0.1)"] = "zprime-M3000-W300";

  setStyle();
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  c->SetGrid();
  TH1F* h = new TH1F("h", "h", 200, 0, 2000);

  h->GetXaxis()->SetTitle("Jet p_{T} (GeV)");
  h->GetXaxis()->SetRangeUser(0, 1400);
  h->GetYaxis()->SetTitle("b-tagging #varepsilon");
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetRangeUser(0, 0.5);
  h->Draw();

  //TLegend* leg = new TLegend(.5,.9-2*2*0.04,.7,.9);
  TLegend* leg = new TLegend(.5,.9-2*0.04,.7,.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);

  //hists for plots
  /*map<TString, TH1F*> hists;
  for (map<TString, TString>::iterator it = files.begin(); it != files.end(); it++) {

    TString dataset = it->second;
    int j=-1;

    if (dataset.EqualTo("TTbar")) j = 1;
    else if (dataset.EqualTo("zprime-M3000-W300")) j = 0;
    else continue;

    TFile* file = TFile::Open( dataset + "_" + channel + ".root" );
    TH1F* jetPt = (TH1F*) file->Get( "jetPt_" + flavor );
    jetPt->Sumw2();
    jetPt = (TH1F*) jetPt->Rebin(nBins, "jetPt", bins);

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
  leg->Draw();*/

  TFile* file = TFile::Open( "btag_eff.root" );
  TH1F* h1 = (TH1F*) file->Get( channel + "/" + flavor + "_LWP_" + channel );
  TH1F* h2 = (TH1F*) file->Get( channel + "/" + flavor + "_MWP_" + channel );

  h1->SetMarkerStyle(20);
  h1->SetMarkerColor(1);
  h1->SetLineColor(1);
  h1->Draw("pesame");

  h2->SetMarkerStyle(24);
  h2->SetMarkerColor(1);
  h2->SetLineColor(1);
  h2->Draw("pesame");

  leg->AddEntry(h1, "CSVL", "PL");
  leg->AddEntry(h2, "CSVM", "PL");
  leg->Draw();

  TLatex text;
  text.SetNDC();

  text.SetTextSize(0.05);
  text.SetTextFont(61);
  text.DrawLatex(0.18, 0.96, "CMS");

  text.SetTextSize(0.04);
  text.DrawLatex(0.2, 0.87, flavor + "-Jets");

  //TString channel_text = "e#mu";
  //if (channel.EqualTo("mm")) channel_text = "#mu#mu";
  //else if (channel.EqualTo("ee")) channel_text = "ee";
  //text.DrawLatex(0.2, 0.83, channel_text + " channel");
  TLatex text2;
  text2.SetNDC();
  text2.SetTextSize(0.04);
  text2.DrawLatex(0.2, 0.83, "e#mu channel");

  text.SetTextFont(52);
  text.DrawLatex(0.29, 0.96, "Simulation");

  text.SetTextSize(0.035);
  text.SetTextFont(42);
  text.DrawLatex(0.85, 0.96, "(13 TeV)");

  c->Print("btag_efficiency_" + flavor + "_" + channel + ".pdf");

  if (!createFile) return;

  //hists for root file
  map<TString, TH1F*> effs;

  const int nChannels = 2;
  TString channels[nChannels] = {"mm", "em"};
  for (int i=0; i<nChannels; i++) {

    const int nFlavors = 3;
    TString flavors[nFlavors] = {"b", "c", "udsg"};
    for (int j=0; j<nFlavors; j++) {

      TString LWP_num = flavors[j] + "_LWP_" + channels[i];
      effs[LWP_num] = new TH1F("LWP_num", "LWP_num", nBins, bins);
      TString MWP_num = flavors[j] + "_MWP_" + channels[i];
      effs[MWP_num] = new TH1F("MWP_num", "MWP_num", nBins, bins);
      TString denom = flavors[j] + "_denom_" + channels[i];
      effs[denom] = new TH1F("denom", "denom", nBins, bins);

      for (map<TString, TString>::iterator it = files.begin(); it != files.end(); it++) {

        TFile* file = TFile::Open( it->second + "_" + channels[i] + ".root" );

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

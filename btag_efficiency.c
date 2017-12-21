// Chad Harrington and Bahareh Roozbahani 5/11/2017 - btagging efficiency

void setStyle();

void btag_efficiency(TString dir = "off/") {

  vector<double> bins = {0, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1000, 1400, 2000};

  map<TString, TString> channels = { {"mm","#mu#mu"}, {"ee","ee"}, {"em","e#mu"} }; // channel, label
  map<TString, float>   flavors  = { {"b",1.4}, {"c",1.4}, {"udsg",0.5} };          // flavor, ymax

  vector <TString> files = {
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
  };

  setStyle();
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  c->SetGrid();
  TH1F* h = new TH1F("h", "h", bins.size()-1, &bins[0]);

  h->GetXaxis()->SetTitle("Jet p_{T} (GeV)");
  h->GetXaxis()->SetRangeUser(0, 1400);
  h->GetYaxis()->SetTitle("b-tagging #varepsilon");
  h->GetYaxis()->SetTitleOffset(1.2);

  TLegend* leg = new TLegend(.5,.9-2*0.04,.7,.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);

  TLatex text;
  text.SetNDC();

  map<TString, TH1F*> effs;
  for (auto const& it_chan : channels) {
    TString chan = it_chan.first, clabel = it_chan.second;

    for (auto const& it_flav : flavors) {
      TString flav = it_flav.first;  float ymax = it_flav.second;

      TString LWP = flav + "_LWP_" + chan;
      TString MWP = flav + "_MWP_" + chan;

      effs[LWP] = new TH1F(LWP, LWP, bins.size()-1, &bins[0]);
      effs[LWP]->Sumw2();
      effs[LWP]->SetBit(TH1::kIsAverage);
      effs[MWP] = new TH1F(MWP, MWP, bins.size()-1, &bins[0]);
      effs[MWP]->Sumw2();
      effs[MWP]->SetBit(TH1::kIsAverage);

      for (auto const& fname : files) {
        TFile* file = TFile::Open( dir + chan + "/" + fname + "_" + chan + ".root" );

        TH1F* denom = (TH1F*) file->Get( "jetPt_" + flav );
        denom->Sumw2();
        denom->SetBit(TH1::kIsAverage);
        denom = (TH1F*) denom->Rebin(bins.size()-1, "denom", &bins[0]);

        TH1F* numL = (TH1F*) file->Get( "jetPt_bTagL_" + flav );
        numL->Sumw2();
        numL->SetBit(TH1::kIsAverage);
        numL = (TH1F*) numL->Rebin(bins.size()-1, "numL", &bins[0]);
        numL->Divide(denom);
        effs[LWP]->Add(numL);

        TH1F* numM = (TH1F*) file->Get( "jetPt_bTagM_" + flav );
        numM->Sumw2();
        numM->SetBit(TH1::kIsAverage);
        numM = (TH1F*) numM->Rebin(bins.size()-1, "numM", &bins[0]);
        numM->Divide(denom);
        effs[MWP]->Add(numM);
      }
      h->GetYaxis()->SetRangeUser(0, ymax);
      h->Draw();

      effs[LWP]->SetMarkerStyle(20);
      effs[LWP]->SetMarkerColor(1);
      effs[LWP]->SetLineColor(1);
      effs[LWP]->Draw("pesame");

      effs[MWP]->SetMarkerStyle(24);
      effs[MWP]->SetMarkerColor(1);
      effs[MWP]->SetLineColor(1);
      effs[MWP]->Draw("pesame");

      leg->AddEntry(effs[LWP], "CSVL", "PL");
      leg->AddEntry(effs[MWP], "CSVM", "PL");
      leg->Draw();

      text.SetTextSize(0.05);  text.SetTextFont(61); text.DrawLatex(0.18, 0.96, "CMS");
      text.SetTextSize(0.04);  text.SetTextFont(52); text.DrawLatex(0.29, 0.96, "Simulation");
      text.SetTextSize(0.035); text.SetTextFont(42); text.DrawLatex(0.85, 0.96, "(13 TeV)");

      text.SetTextSize(0.04);  text.SetTextFont(62);
      text.DrawLatex(0.2, 0.87, flav + "-Jets");
      text.DrawLatex(0.2, 0.83, clabel + " channel");

      c->Print("plots/btag_eff_" + flav + "_" + chan + ".pdf");
      leg->Clear();
      c->Clear();
    }
  }

  TFile* outFile = new TFile("btag_eff.root","RECREATE");
  outFile->cd();

  for (auto const& it_chan : channels) outFile->mkdir( it_chan.first + "/" );

  for (auto const& it_eff : effs) {
    outFile->cd();
    TString name = it_eff.first;
    TString postfix = name(name.Last('_')+1, 2);

    outFile->cd("btag_eff.root:/" + postfix);
    it_eff.second->Write(name);
  }

  outFile->Write();
  delete outFile;
  outFile = 0;

  //////////////////////////////////////////////////
  //Plot individual efficiencies on same histogram//
  //////////////////////////////////////////////////

  vector <TString> filesToPlot = { "TTbar0-700", "TTbar700-1000", "TTbar1000-inf", "zprime_M-3000_W-300" };
  leg->SetY1NDC(.9-4*0.04);

  for (auto const& it_chan : channels) {
    TString chan = it_chan.first, clabel = it_chan.second;

    for (auto const& it_flav : flavors) {
      TString flav = it_flav.first;  float ymax = it_flav.second;

      TH1F* h_LWP = new TH1F("LWP", "LWP", bins.size()-1, &bins[0]);
      h_LWP->Sumw2();
      h_LWP->SetBit(TH1::kIsAverage);
      TH1F* h_MWP = new TH1F("MWP", "MWP", bins.size()-1, &bins[0]);
      h_MWP->Sumw2();
      h_MWP->SetBit(TH1::kIsAverage);

      TH1F* h_sigLWP = new TH1F("sigLWP", "sigLWP", bins.size()-1, &bins[0]);
      h_sigLWP->Sumw2();
      h_sigLWP->SetBit(TH1::kIsAverage);
      TH1F* h_sigMWP = new TH1F("sigMWP", "sigMWP", bins.size()-1, &bins[0]);
      h_sigMWP->Sumw2();
      h_sigMWP->SetBit(TH1::kIsAverage);

      for (auto const& fname : filesToPlot) {
        TFile* file = TFile::Open( dir + chan + "/" + fname + "_" + chan + ".root" );

        TH1F* denom = (TH1F*) file->Get( "jetPt_" + flav );
        denom->Sumw2();
        denom->SetBit(TH1::kIsAverage);
        denom = (TH1F*) denom->Rebin(bins.size()-1, "denom", &bins[0]);

        TH1F* numL = (TH1F*) file->Get( "jetPt_bTagL_" + flav );
        numL->Sumw2();
        numL->SetBit(TH1::kIsAverage);
        numL = (TH1F*) numL->Rebin(bins.size()-1, "numL", &bins[0]);
        numL->Divide(denom);

        TH1F* numM = (TH1F*) file->Get( "jetPt_bTagM_" + flav );
        numM->Sumw2();
        numM->SetBit(TH1::kIsAverage);
        numM = (TH1F*) numM->Rebin(bins.size()-1, "numM", &bins[0]);
        numM->Divide(denom);

        if ( fname.Contains("zprime", TString::kIgnoreCase) ) {
          h_sigLWP->Add(numL);
          h_sigMWP->Add(numM);
        }
        else {
          h_LWP->Add(numL);
          h_MWP->Add(numM);
        }
      }
      h->GetYaxis()->SetRangeUser(0, ymax);
      h->Draw();

      h_LWP->SetMarkerStyle(20);
      h_LWP->SetMarkerColor(kBlack);
      h_LWP->SetLineColor(kBlack);
      h_LWP->Draw("pesame");

      h_MWP->SetMarkerStyle(24);
      h_MWP->SetMarkerColor(kBlack);
      h_MWP->SetLineColor(kBlack);
      h_MWP->Draw("pesame");

      h_sigLWP->SetMarkerStyle(21);
      h_sigLWP->SetMarkerColor(kRed);
      h_sigLWP->SetLineColor(kRed);
      h_sigLWP->Draw("pesame");

      h_sigMWP->SetMarkerStyle(25);
      h_sigMWP->SetMarkerColor(kRed);
      h_sigMWP->SetLineColor(kRed);
      h_sigMWP->Draw("pesame");

      leg->AddEntry(h_LWP, "CSVL,  ttbar", "PL");
      leg->AddEntry(h_MWP, "CSVM, ttbar", "PL");
      leg->AddEntry(h_sigLWP, "CSVL,  Z' (10%, 3 Tev)", "PL");
      leg->AddEntry(h_sigMWP, "CSVM, Z' (10%, 3 Tev)", "PL");
      leg->Draw();

      text.SetTextSize(0.05);  text.SetTextFont(61); text.DrawLatex(0.18, 0.96, "CMS");
      text.SetTextSize(0.04);  text.SetTextFont(52); text.DrawLatex(0.29, 0.96, "Simulation");
      text.SetTextSize(0.035); text.SetTextFont(42); text.DrawLatex(0.85, 0.96, "(13 TeV)");

      text.SetTextSize(0.04);  text.SetTextFont(62);
      text.DrawLatex(0.2, 0.87, flav + "-Jets");
      text.DrawLatex(0.2, 0.83, clabel + " channel");

      c->Print("plots/btag_eff_" + flav + "_" + chan + "_individual.pdf");
      leg->Clear();
      c->Clear();
    }
  }
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

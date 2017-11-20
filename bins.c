//Chad Harrington 10/23/2017

bool createBins(TH1D* h, vector<double>& bins, const int& upper);
bool createBins(TH1D* h1, TH1D* h2, vector<double>& bins, const int& upper);
int dwidth = 10, rebin = 4, minwidth = dwidth*rebin;

void bins(TString hname = "7_sT_met", TString dir="off/", TString channel="mm", bool chi2=false) {

  TH1D* h=0, *h_rebin=0, *h_data=0, *h_data_rebin=0;
  vector<TString> files = { "TTbar0-700", "TTbar700-1000", "TTbar1000-inf" };

  if (chi2) { //look at total background and data
    files.insert( files.begin(), {"DYhigh", "DYlow", "STschannel", "STtWchannel", "STtchannel", "SaTtWchannel", "SaTtchannel", "WJets", "WW", "WZ", "ZZ"} );

    TString dataname = channel=="ee" ? "Ele" : "Muon";
    TFile* dataFile = TFile::Open(dir + channel + "/" + dataname + "_" + channel + ".root");
    h_data = (TH1D*) dataFile->FindObjectAny(hname);
    h_data_rebin = (TH1D*) h_data->Rebin(rebin, "h_data_rebin");
  }

  for (int i=0, n=files.size(); i<n; i++) {
    TFile* inFile = TFile::Open(dir + channel + "/" + files[i] + "_" + channel + ".root");

    if (h==0) h = (TH1D*) inFile->FindObjectAny(hname);
    else      h->Add( (TH1D*) inFile->FindObjectAny(hname) );
  }
  int upper = h->GetBinLowEdge(h->GetNbinsX()+1);

  h_rebin = (TH1D*) h->Rebin(rebin, "h_rebin");
  vector<double> bins;

  if (chi2) {
    while (createBins(h_rebin, h_data_rebin, bins, upper)) {
      h_rebin = (TH1D*) h->Rebin(bins.size()-1, "h_rebin", &bins[0]);
      h_data_rebin = (TH1D*) h_data->Rebin(bins.size()-1, "h_data_rebin", &bins[0]);
    }
  }
  else {
    while (createBins(h_rebin, bins, upper)) h_rebin = (TH1D*) h->Rebin(bins.size()-1, "h_rebin", &bins[0]);
  }

  for (int i=0, n=bins.size(); i<n; i++) cout << bins[i] << " ";
  cout << endl;

  int nBins = h_rebin->GetNbinsX();
  for (int i=1; i<=nBins; i++) {
    cout << h_rebin->GetBinLowEdge(i) << "\t" << h_rebin->GetBinLowEdge(i+1) << "\t" << h_rebin->GetBinContent(i) << "\t" << h_rebin->GetBinError(i)/h_rebin->GetBinContent(i)*100 << endl;
    if (chi2) cout << "\t\t" << h_data_rebin->GetBinContent(i) << "\t" << h_data_rebin->GetBinError(i)/h_data_rebin->GetBinContent(i)*100 << " (Data)" << endl;
  }
  cout << nBins << endl;
}

//returns true if we need to create more bins (i.e. error/content > 20% for a bin)
bool createBins(TH1D* h, vector<double>& bins, const int& upper) {
  bins.clear();
  bins.push_back(0);

  int nBins = h->GetNbinsX();
  for (int i=1; i<=nBins; i++) {
    if ( h->GetBinContent(i) == 0 || h->GetBinError(i)/h->GetBinContent(i) > 0.2 ) {
      if (i==nBins) bins.pop_back();
      else {
        bins.push_back( h->GetBinLowEdge(i+1) + dwidth );
        while ( bins.back() + minwidth <= upper - minwidth ) bins.push_back( bins.back() + minwidth );
      }
      if (bins.back() != upper) bins.push_back(upper);
      return true;
    }
    else bins.push_back( h->GetBinLowEdge(i+1) );
  }
  return false;
}

//returns true if we need to create more bins (i.e. error/content > 30% for a bin)
//data and background
bool createBins(TH1D* h1, TH1D* h2, vector<double>& bins, const int& upper) {
  bins.clear();
  bins.push_back(0);

  int nBins = h1->GetNbinsX();
  for (int i=1; i<=nBins; i++) {
    if ( h1->GetBinContent(i) == 0 || h2->GetBinContent(i) == 0 || h1->GetBinError(i)/h1->GetBinContent(i) > 0.3 || h2->GetBinError(i)/h2->GetBinContent(i) > 0.3) {
      if (i==nBins) bins.pop_back();
      else {
        bins.push_back( h1->GetBinLowEdge(i+1) + dwidth );
        while ( bins.back() + minwidth <= upper - minwidth ) bins.push_back( bins.back() + minwidth );
      }
      if (bins.back() != upper) bins.push_back(upper);
      return true;
    }
    else bins.push_back( h1->GetBinLowEdge(i+1) );
  }
  return false;
}

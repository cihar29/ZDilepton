//Chad Harrington 10/23/2017

bool createBins(TH1D* h, vector<double>& bins);

void bins(TString channel="mm") {

  TString files[] = {
    "TTbar0-700",
    "TTbar700-1000",
    "TTbar1000-inf"
  };

  TH1D* h = 0;
  TString hname = "5_sT_met";

  int n = sizeof(files)/sizeof(*files);
  for (int i=0; i<n; i++) {

    TFile* inFile = TFile::Open(channel + "/" + files[i] + "_" + channel + ".root");

    if (h==0) h = (TH1D*) inFile->FindObjectAny(hname);
    else      h->Add( (TH1D*) inFile->FindObjectAny(hname) );
  }
  h->Rebin(2);
  vector<double> bins;

  while (createBins(h, bins)) h = (TH1D*) h->Rebin(bins.size()-1, "h", &bins[0]);

  for (int i=0, n=bins.size(); i<n; i++) cout << bins[i] << " "; cout << endl;

  //int nBins = h->GetNbinsX();
  //for (int i=1; i<=nBins; i++) cout << h->GetBinCenter(i) << "\t" << h->GetBinContent(i) << "\t" << h->GetBinError(i)/h->GetBinContent(i) << endl;
}

//returns true if we need to create more bins (i.e. error/content > 20% for a bin)
bool createBins(TH1D* h, vector<double>& bins) {
  bins.clear();

  int nBins = h->GetNbinsX();
  int upper = h->GetBinLowEdge(nBins+1);

  for (int i=1; i<=nBins; i++) {
    if ( h->GetBinContent(i)==0 || h->GetBinError(i)/h->GetBinContent(i) < 0.2 ) bins.push_back( h->GetBinLowEdge(i) );
    else {
      double new_width = h->GetBinWidth(i) + h->GetBinWidth(i+1);

      if (h->GetBinLowEdge(i) + new_width < upper-new_width) {
        bins.push_back( h->GetBinLowEdge(i) );
        while ( bins.back() + new_width < upper-new_width ) bins.push_back( bins.back() + new_width );
      }
      bins.push_back( upper );
      return true;
    }
  }
  bins.push_back( upper );
  return false;
}

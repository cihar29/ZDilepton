//Chad Harrington 10/23/2017

bool createBins(TH1D* h, vector<double>& bins, const int& upper);
int dwidth = 10, rebin = 4, minwidth = dwidth*rebin;
float threshold = 0.2;

void bins(TString hname = "7_sT_met", TString dir="off/", TString channel="mm") {

  TString files[] = { "TTbar0-700", "TTbar700-1000", "TTbar1000-inf" };

  TH1D* h = 0;

  int n = sizeof(files)/sizeof(*files);
  for (int i=0; i<n; i++) {

    TFile* inFile = TFile::Open(dir + channel + "/" + files[i] + "_" + channel + ".root");

    if (h==0) h = (TH1D*) inFile->FindObjectAny(hname);
    else      h->Add( (TH1D*) inFile->FindObjectAny(hname) );
  }
  int upper = h->GetBinLowEdge(h->GetNbinsX()+1);
  h->Rebin(rebin);
  vector<double> bins;

  while (createBins(h, bins, upper)) h = (TH1D*) h->Rebin(bins.size()-1, "h", &bins[0]);

  for (int i=0, n=bins.size(); i<n; i++) cout << bins[i] << " ";
  cout << endl;

  int nBins = h->GetNbinsX();
  for (int i=1; i<=nBins; i++)
    cout << h->GetBinLowEdge(i) << "\t" << h->GetBinLowEdge(i+1) << "\t" << h->GetBinContent(i) << "\t" << h->GetBinError(i)/h->GetBinContent(i) << endl;
  cout << nBins << endl;
}

//returns true if we need to create more bins (i.e. error/content > 20% for a bin)
bool createBins(TH1D* h, vector<double>& bins, const int& upper) {
  bins.clear();
  bins.push_back(0);

  int j=1;
  while (h->GetBinContent(j)==0) j++;

  int nBins = h->GetNbinsX();
  for (int i=j; i<=nBins; i++) {
    if ( h->GetBinError(i)/h->GetBinContent(i) > threshold ) {
      if (i==nBins) bins.pop_back();
      else {
        bins.push_back( h->GetBinLowEdge(i+1) + minwidth );
        while ( bins.back() + minwidth < upper - minwidth ) bins.push_back( bins.back() + minwidth );
      }
      if (bins.back() != upper) bins.push_back(upper);
      return true;
    }
    else bins.push_back( h->GetBinLowEdge(i+1) );
  }
  return false;
}


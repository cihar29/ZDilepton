void signalRMS() {

  vector<int> masses = { 500, 750, 1000, 1250, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 6500, 7000 };
  vector<TString> signals = { "gkk", "zp1", "zp10", "zp30" };

  ofstream ofile( "signalRMS.txt" );
  ofile << "#mass";
  for ( const TString& sig : signals ) ofile << "\t" << sig;
  ofile << endl;

  for ( const int& mass : masses ) {
    ofile << mass;

    for ( const TString& sig : signals ) {
      float width;
      TString fname;
      if ( sig == "gkk" ) {
        width = 0.14*mass;
        fname = Form("gluon_M-%i.root", mass);
      }
      else {
        if      ( sig == "zp1" )  width = 0.01*mass;
        else if ( sig == "zp10" ) width = 0.1*mass;
        else if ( sig == "zp30" ) width = 0.3*mass;

        TString width_str = Form("%.1f", width);
        width_str.ReplaceAll( ".5", "p5" );
        width_str.ReplaceAll( ".0", "" );
        fname = Form("zprime_M-%i_W-%s.root", mass, width_str.Data());
      }
      TFile* file = TFile::Open( "/uscms_data/d3/cihar29/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/root_trees/" + fname );
      if (file == 0) {
        ofile << "\tnan";
        continue;
      }
      TTree* tree = (TTree*) file->Get("T");

      TH1F* h_mass = new TH1F("mass", "mass", 200, mass-3*width, mass+3*width);
      tree->Draw("gen_mass>>mass", "gen_PID>1000000 && gen_status==22", "histsame");

      TF1* f = new TF1( "f", "gaus" );
      h_mass->Fit("f", "Q");
      cout << mass << "\t" << sig << "\t" << f->GetParameter(2) << "\t" << f->GetChisquare() / f->GetNDF() << endl;

      ofile << "\t" << Form("%.0f", f->GetParameter(2));
      delete h_mass;
      delete f;
    }
    ofile << endl;
  }

  ofile.close();
}
//use this command  make cleaner:
//cat signalRMS.txt | column -t

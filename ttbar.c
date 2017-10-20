//Chad Harrington - 10/19/2017

void ttbar() {
  
  TFile* f_incl = TFile::Open("root_trees/TTbar_incl.root");

  TTree* tree = (TTree*) f_incl->Get("T");
  tree->SetBranchStatus("*",1);
  TH2F* mu_mass = (TH2F*) f_incl->Get("fullMu_mass");

  TFile *f1 = new TFile("TTbar0-700.root","recreate");
  TTree *t1 = tree->CloneTree(0);
  TH1D  *htotal1 = new TH1D("mass_totalEvts","mass_totalEvts",500,0,5000);
  TH1D  *hmet1 = new TH1D("mass_met","mass_met",500,0,5000);
  TH1D  *mu1_double = mu_mass->ProjectionY("fullMu", 1, 70);  //Want mu to be TH1F* and match other mu distributions
  TH1F  mu1;  mu1_double->SetTitle("fullMu");  mu1_double->Copy(mu1);  delete mu1_double;

  TFile *f2 = new TFile("TTbar700-1000_incl.root","recreate");
  TTree *t2 = tree->CloneTree(0);
  TH1D  *htotal2 = new TH1D("mass_totalEvts","mass_totalEvts",500,0,5000);
  TH1D  *hmet2 = new TH1D("mass_met","mass_met",500,0,5000);
  TH1D  *mu2_double = mu_mass->ProjectionY("fullMu", 71, 100);
  TH1F  mu2;  mu2_double->SetTitle("fullMu");  mu2_double->Copy(mu2);  delete mu2_double;

  TFile *f3 = new TFile("TTbar1000-inf_incl.root","recreate");
  TTree *t3 = tree->CloneTree(0);
  TH1D  *htotal3 = new TH1D("mass_totalEvts","mass_totalEvts",500,0,5000);
  TH1D  *hmet3 = new TH1D("mass_met","mass_met",500,0,5000);
  TH1D  *mu3_double = mu_mass->ProjectionY("fullMu", 101, 500);
  TH1F  mu3;  mu3_double->SetTitle("fullMu");  mu3_double->Copy(mu3);  delete mu3_double;

  delete mu_mass;

  TString intnames[] = { "totalEvts", "dilep_cut", "leppt_cut", "dilepmass_cut", "jetpteta_cut", "met_cut" };
  TString doublenames[] = { "topPtWeightNOM", "topPtWeightDN", "pdfUP", "pdfDN", "q2UP", "q2DN" };

  vector< vector<int> > ints1 = { {0}, {0}, {0}, {0}, {0}, {0} };
  vector< vector<int> > ints2 = { {0}, {0}, {0}, {0}, {0}, {0} };
  vector< vector<int> > ints3 = { {0}, {0}, {0}, {0}, {0}, {0} };
  vector< vector<double> > doubles1 = { {0.}, {0.}, {0.}, {0.}, {0.}, {0.} };
  vector< vector<double> > doubles2 = { {0.}, {0.}, {0.}, {0.}, {0.}, {0.} };
  vector< vector<double> > doubles3 = { {0.}, {0.}, {0.}, {0.}, {0.}, {0.} };

  int xs = 831760;
  for (int i=0; i<6; i++) {

    TString hiname = "mass_" + intnames[i];
    hiname.ReplaceAll("_cut", "");
    TString hdname = "mass_" + doublenames[i];

    TH1D* hi = (TH1D*) f_incl->Get(hiname);
    TH1D* hd = (TH1D*) f_incl->Get(hdname);

    ints1[i][0] = int( hi->Integral(1,70) );   //bin 70 is 690-700
    ints2[i][0] = int( hi->Integral(71,100) ); //bin 100 is 990-1000
    ints3[i][0] = int( hi->Integral(101,500) );
    double itotal = hi->Integral(1,500);

    doubles1[i][0] = hd->Integral(1,70);
    doubles2[i][0] = hd->Integral(71,100);
    doubles3[i][0] = hd->Integral(101,500);
    double dtotal = hd->Integral(1,500);

    f1->WriteObject(&ints1[i], intnames[i]);
    f2->WriteObject(&ints2[i], intnames[i]);
    f3->WriteObject(&ints3[i], intnames[i]);

    f1->WriteObject(&doubles1[i], doublenames[i]);
    f2->WriteObject(&doubles2[i], doublenames[i]);
    f3->WriteObject(&doubles3[i], doublenames[i]);

    if (intnames[i] == "totalEvts") {
      cout << "TTbar0-700\t" << ints1[i][0]/itotal*xs << endl;
      cout << "TTbar700-1000\t" << ints2[i][0]/itotal*xs << endl;
      cout << "TTbar1000-inf\t" << ints3[i][0]/itotal*xs << endl;

      for (int i=1; i<=500; i++) {
        double content = hi->GetBinContent(i);

        if      (i<=70)  htotal1->SetBinContent(i, content);
        else if (i<=100) htotal2->SetBinContent(i, content);
        else             htotal3->SetBinContent(i, content);
      }
      htotal1->SetEntries(ints1[i][0]);
      htotal2->SetEntries(ints2[i][0]);
      htotal3->SetEntries(ints3[i][0]);
    }
    cout << "TTbar0-700_" + doublenames[i] << "\t" << doubles1[i][0]/dtotal*xs << endl;
    cout << "TTbar700-1000_" + doublenames[i] << "\t" << doubles2[i][0]/dtotal*xs << endl;
    cout << "TTbar1000-inf_" + doublenames[i] << "\t" << doubles3[i][0]/dtotal*xs << endl;
  }

  int nGen=20, gen_status[nGen], gen_PID[nGen];
  float gen_pt[nGen], gen_mass[nGen], gen_eta[nGen], gen_phi[nGen];

  tree->SetBranchAddress("nGen", &nGen);
  tree->SetBranchAddress("gen_status", gen_status);
  tree->SetBranchAddress("gen_PID", gen_PID);
  tree->SetBranchAddress("gen_pt", gen_pt);
  tree->SetBranchAddress("gen_mass", gen_mass);
  tree->SetBranchAddress("gen_eta", gen_eta);
  tree->SetBranchAddress("gen_phi", gen_phi);

  Long64_t nEntries = tree->GetEntries();
  for (Long64_t n=0; n<nEntries; n++) {
    tree->GetEntry(n);

    //first t's
    TLorentzVector t, tbar;
    for (int i=0; i<2; i++) {
      if      (gen_PID[i]== 6 && gen_status[i]<30) t.SetPtEtaPhiM(gen_pt[i], gen_eta[i], gen_phi[i], gen_mass[i]);
      else if (gen_PID[i]==-6 && gen_status[i]<30) tbar.SetPtEtaPhiM(gen_pt[i], gen_eta[i], gen_phi[i], gen_mass[i]);
    }
    double mass = (t+tbar).M();

    if      (mass < 700)  { t1->Fill(); hmet1->Fill(mass); }
    else if (mass < 1000) { t2->Fill(); hmet2->Fill(mass); }
    else                  { t3->Fill(); hmet3->Fill(mass); }
  }

  f1->Write();
  f2->Write();
  f3->Write();
  delete f_incl;
  delete f1;
  delete f2;
  delete f3;
  cout << "done" << endl;
}

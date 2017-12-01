double delta(double delta1, double delta2, string plus_minus) {
// pos_neg should be either "+" or "-"

  double del = 0. ;
  if (plus_minus=="+") {
    double dd = (delta1 > delta2) ? delta1 : delta2 ;
    if (dd > 0.) del = dd ;
  }
  else if (plus_minus=="-") {
    double dd = (delta1 < delta2) ? delta1 : delta2 ;
    if (dd < 0.) del = dd ;
  }
  return del;
}

void post_table(TString folder = "all_st/gkk/", TString region = "rev") {

  TString dir = "/uscms_data/d3/cihar29/Analysis/CMSSW_8_1_0/src/theta/utils2/2017/";
  TFile* file = TFile::Open(dir + folder + "histo_fits.root");

  vector<TString> channels = { "mm", "ee", "em" }, sets = { "ttbar", "dy", "st", "vv", "wjet", "data" }, systematics;
  map< TString, map<TString, double> > m_NOM, m_errUP, m_errDN;         //channel, set
  map< TString, map<TString, map<TString, double> > > m_sysUP, m_sysDN; //channel, set, sys

  m_NOM["ll"]["bkg"] = 0;
  for (auto const& chan : channels) m_NOM[chan]["bkg"] = 0;
  for (auto const& set : sets)      m_NOM["ll"][set] = 0;

  //Get Systematics
  TIter nextkey(file->GetListOfKeys());
  TKey* key;
  while ( (key = (TKey*)nextkey()) ) {
    TString keyname = key->GetName();

    if (keyname.Contains("UP")) {
      int index = keyname.Index("__", 2, keyname.Index("__")+2, TString::kExact);  //index of second "__"
      TString sys = keyname(index+2, keyname.Index("UP")-index-2);

      unsigned int idx;
      for (idx=0; idx<systematics.size(); idx++) { if (systematics[idx] == sys) break; }
      if (idx == systematics.size()) {
        systematics.push_back(sys);

                                            m_sysUP["ll"]["bkg"][sys] = 0;  m_sysDN["ll"]["bkg"][sys] = 0;
        for (auto const& chan : channels) { m_sysUP[chan]["bkg"][sys] = 0;  m_sysDN[chan]["bkg"][sys] = 0; }
        for (auto const& set : sets)      { m_sysUP["ll"][set][sys] = 0;    m_sysDN["ll"][set][sys] = 0;   }
      }
    }
  }

  for (auto const& chan : channels) {
    for (auto const& set : sets) {
      TH1D* h = (TH1D*) file->FindObjectAny(chan + region + "__" + set);
      m_NOM[chan][set] = h->Integral();
      m_NOM["ll"][set] += m_NOM[chan][set];

      if (set == "data") continue;
      m_NOM[chan]["bkg"] += m_NOM[chan][set];
      m_NOM["ll"]["bkg"] += m_NOM[chan][set];

      for (auto const& sys : systematics) {
        TH1D* hUP = (TH1D*) file->FindObjectAny(chan + region + "__" + set + "__" + sys + "UP");
        TH1D* hDN = (TH1D*) file->FindObjectAny(chan + region + "__" + set + "__" + sys + "DN");

        m_sysUP[chan][set][sys] = hUP->Integral();             m_sysDN[chan][set][sys] = hDN->Integral();
        m_sysUP["ll"][set][sys]   += m_sysUP[chan][set][sys];  m_sysDN["ll"][set][sys]   += m_sysDN[chan][set][sys];
        m_sysUP[chan]["bkg"][sys] += m_sysUP[chan][set][sys];  m_sysDN[chan]["bkg"][sys] += m_sysDN[chan][set][sys];
        m_sysUP["ll"]["bkg"][sys] += m_sysUP[chan][set][sys];  m_sysDN["ll"]["bkg"][sys] += m_sysDN[chan][set][sys];
      }
    }
  }

  for (auto const& i_chan : m_NOM) {
    TString chan = i_chan.first;

    for (auto const& i_set : i_chan.second) {
      TString set = i_set.first;
      double nom = i_set.second, errorUP=0, errorDN=0;

      for (auto const& sys : systematics) {
        double deltaUP = m_sysUP[chan][set][sys] - nom;
        double deltaDN = m_sysDN[chan][set][sys] - nom;

        double d1 = delta(deltaUP, deltaDN, "+");
        double d2 = delta(deltaUP, deltaDN, "-");

        errorUP += d1*d1;
        errorDN += d2*d2;
      }
      m_errUP[chan][set] = sqrt(errorUP);
      m_errDN[chan][set] = sqrt(errorDN);
    }
  }

  ofstream outfile("postfit_" + region + ".txt");

  outfile << "\\begin{tabular}{ |c|c|c|c|c| }\n";
  outfile << "\\multicolumn{5}{c}{" << region << "} \\\\\n";
  outfile << "\\hline\n";
  outfile << "Sample & $\\mu\\mu$ Channel & ee Channel & e$\\mu$ Channel & Combined \\\\\n";
  outfile << "\\hline\n";

  outfile << "t$\\bar{\\textrm{t}}$ &" << Form("$%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ \\\\\n",
    m_NOM["mm"]["ttbar"], m_errUP["mm"]["ttbar"], m_errDN["mm"]["ttbar"], m_NOM["ee"]["ttbar"], m_errUP["ee"]["ttbar"], m_errDN["ee"]["ttbar"],
    m_NOM["em"]["ttbar"], m_errUP["em"]["ttbar"], m_errDN["em"]["ttbar"], m_NOM["ll"]["ttbar"], m_errUP["ll"]["ttbar"], m_errDN["ll"]["ttbar"]
  ) << endl;
  outfile << "Z/$\\gamma^{*}\\rightarrow l^{+}l^{-}$ &" << Form("$%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ \\\\\n",
    m_NOM["mm"]["dy"], m_errUP["mm"]["dy"], m_errDN["mm"]["dy"], m_NOM["ee"]["dy"], m_errUP["ee"]["dy"], m_errDN["ee"]["dy"],
    m_NOM["em"]["dy"], m_errUP["em"]["dy"], m_errDN["em"]["dy"], m_NOM["ll"]["dy"], m_errUP["ll"]["dy"], m_errDN["ll"]["dy"]
  ) << endl;
  outfile << "Single-Top &" << Form("$%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ \\\\\n",
    m_NOM["mm"]["st"], m_errUP["mm"]["st"], m_errDN["mm"]["st"], m_NOM["ee"]["st"], m_errUP["ee"]["st"], m_errDN["ee"]["st"],
    m_NOM["em"]["st"], m_errUP["em"]["st"], m_errDN["em"]["st"], m_NOM["ll"]["st"], m_errUP["ll"]["st"], m_errDN["ll"]["st"]
  ) << endl;
  outfile << "VV &" << Form("$%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ \\\\\n",
    m_NOM["mm"]["vv"], m_errUP["mm"]["vv"], m_errDN["mm"]["vv"], m_NOM["ee"]["vv"], m_errUP["ee"]["vv"], m_errDN["ee"]["vv"],
    m_NOM["em"]["vv"], m_errUP["em"]["vv"], m_errDN["em"]["vv"], m_NOM["ll"]["vv"], m_errUP["ll"]["vv"], m_errDN["ll"]["vv"]
  ) << endl;
  outfile << "W+Jets &" << Form("$%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ \\\\\n",
    m_NOM["mm"]["wjet"], m_errUP["mm"]["wjet"], m_errDN["mm"]["wjet"], m_NOM["ee"]["wjet"], m_errUP["ee"]["wjet"], m_errDN["ee"]["wjet"],
    m_NOM["em"]["wjet"], m_errUP["em"]["wjet"], m_errDN["em"]["wjet"], m_NOM["ll"]["wjet"], m_errUP["ll"]["wjet"], m_errDN["ll"]["wjet"]
  ) << endl;
  outfile << "\\hline\n";

  outfile << "Total Bkg &" << Form("$%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ \\\\\n",
    m_NOM["mm"]["bkg"], m_errUP["mm"]["bkg"], m_errDN["mm"]["bkg"], m_NOM["ee"]["bkg"], m_errUP["ee"]["bkg"], m_errDN["ee"]["bkg"],
    m_NOM["em"]["bkg"], m_errUP["em"]["bkg"], m_errDN["em"]["bkg"], m_NOM["ll"]["bkg"], m_errUP["ll"]["bkg"], m_errDN["ll"]["bkg"]
  ) << endl;
  outfile << "Data &" << Form("$%.0f \\pm %.1f$ & $%.0f \\pm %.1f$ & $%.0f \\pm %.1f$ & $%.0f \\pm %.1f$ \\\\\n",
    m_NOM["mm"]["data"], sqrt(m_NOM["mm"]["data"]), m_NOM["ee"]["data"], sqrt(m_NOM["ee"]["data"]),
    m_NOM["em"]["data"], sqrt(m_NOM["em"]["data"]), m_NOM["ll"]["data"], sqrt(m_NOM["ll"]["data"])
  ) << endl;
  outfile << "Data/Bkg &" << Form("$%.2f^{+%.2f}_{-%.2f}$ & $%.2f^{+%.2f}_{-%.2f}$ & $%.2f^{+%.2f}_{-%.2f}$ & $%.2f^{+%.2f}_{-%.2f}$ \\\\\n",
    m_NOM["mm"]["data"]/m_NOM["mm"]["bkg"],
    m_NOM["mm"]["data"]/m_NOM["mm"]["bkg"] * sqrt( 1/m_NOM["mm"]["data"] + m_errUP["mm"]["bkg"]*m_errUP["mm"]["bkg"]/(m_NOM["mm"]["bkg"]*m_NOM["mm"]["bkg"]) ),
    m_NOM["mm"]["data"]/m_NOM["mm"]["bkg"] * sqrt( 1/m_NOM["mm"]["data"] + m_errDN["mm"]["bkg"]*m_errDN["mm"]["bkg"]/(m_NOM["mm"]["bkg"]*m_NOM["mm"]["bkg"]) ),
    m_NOM["ee"]["data"]/m_NOM["ee"]["bkg"],
    m_NOM["ee"]["data"]/m_NOM["ee"]["bkg"] * sqrt( 1/m_NOM["ee"]["data"] + m_errUP["ee"]["bkg"]*m_errUP["ee"]["bkg"]/(m_NOM["ee"]["bkg"]*m_NOM["ee"]["bkg"]) ),
    m_NOM["ee"]["data"]/m_NOM["ee"]["bkg"] * sqrt( 1/m_NOM["ee"]["data"] + m_errDN["ee"]["bkg"]*m_errDN["ee"]["bkg"]/(m_NOM["ee"]["bkg"]*m_NOM["ee"]["bkg"]) ),
    m_NOM["em"]["data"]/m_NOM["em"]["bkg"],
    m_NOM["em"]["data"]/m_NOM["em"]["bkg"] * sqrt( 1/m_NOM["em"]["data"] + m_errUP["em"]["bkg"]*m_errUP["em"]["bkg"]/(m_NOM["em"]["bkg"]*m_NOM["em"]["bkg"]) ),
    m_NOM["em"]["data"]/m_NOM["em"]["bkg"] * sqrt( 1/m_NOM["em"]["data"] + m_errDN["em"]["bkg"]*m_errDN["em"]["bkg"]/(m_NOM["em"]["bkg"]*m_NOM["em"]["bkg"]) ),
    m_NOM["ll"]["data"]/m_NOM["ll"]["bkg"],
    m_NOM["ll"]["data"]/m_NOM["ll"]["bkg"] * sqrt( 1/m_NOM["ll"]["data"] + m_errUP["ll"]["bkg"]*m_errUP["ll"]["bkg"]/(m_NOM["ll"]["bkg"]*m_NOM["ll"]["bkg"]) ),
    m_NOM["ll"]["data"]/m_NOM["ll"]["bkg"] * sqrt( 1/m_NOM["ll"]["data"] + m_errDN["ll"]["bkg"]*m_errDN["ll"]["bkg"]/(m_NOM["ll"]["bkg"]*m_NOM["ll"]["bkg"]) )
  ) << endl;
  outfile << "\\hline\n";
  outfile << "\\end{tabular}\n" << endl;
}

//Chad Harrington - 11/15/2017

string readFile(istream& file, const string& start, vector< vector<string> >& vec);
double delta(double delta1, double delta2, string plus_minus);

string boldline = "====================================================================================================================================================================";
string divider = "--------------------------------------------------------------------------------------------------------------------------------------------------------------------";

void combine_logs(string dir="off/") {

  //channel, sample, cut, value
  map<string, map<string, vector<double> > > m_nom, m_stat, m_sysUP, m_sysDN;

  vector<double> cuts; cuts.assign(24, 0.);
  map<string, vector<double> > samples = { {"data", cuts}, {"bkgd", cuts}, {"ttbar", cuts}, {"dy", cuts}, {"st", cuts},
                                           {"vv", cuts}, {"wjet", cuts}, {"zprime", cuts}, {"gluon", cuts} };
  m_nom   = { {"ll", samples}, {"mm", samples}, {"ee", samples}, {"em", samples} };
  m_stat  = { {"ll", samples}, {"mm", samples}, {"ee", samples}, {"em", samples} };
  m_sysUP = { {"ll", samples}, {"mm", samples}, {"ee", samples}, {"em", samples} };
  m_sysDN = { {"ll", samples}, {"mm", samples}, {"ee", samples}, {"em", samples} };

  //sys, sample, cut, value
  map<string, map<string, vector<double> > > ll_sysUP, ll_sysDN;  //individual systematics, as opposed to total systematics
  vector<string> cutnames;

  ofstream outfile(dir + "logs/ll_cutflow.txt");

  ifstream mmfile(dir + "logs/mm_cutflow.txt");
  ifstream eefile(dir + "logs/ee_cutflow.txt");
  ifstream emfile(dir + "logs/em_cutflow.txt");

  string headers[] = {"Cut Flow Table: Summary", "Cut Flow Table: Background", "Cut Flow Table: Statistical Error", "Systematics: luminosity",
    "Systematics: \u03C3(single-top)", "Systematics: \u03C3(diboson)", "Systematics: \u03BC trigger", "Systematics: \u03BC ID", "Systematics: \u03BC Isolation",
    "Systematics: e trigger", "Systematics: e ID", "Systematics: e Isolation", "Systematics: JEC", "Systematics: JER", "Systematics: b-tagging",
    "Systematics: mis-tagging", "Systematics: pileup", "Systematics: top pT modeling", "Systematics: PDF", "Systematics: Q2 ttbar", "Systematics: Q2 DY",
    "Systematics: Q2 Single-top", "Systematics: Q2 signal"};

  int nheaders = sizeof(headers)/sizeof(headers[0]);
  for (int ih=0; ih<nheaders; ih++) {
    string header = headers[ih];

    vector< vector<string> > v_mm, v_ee, v_em;
    vector< vector<string> > v_driver;
    string label = "";
    if      (header.find("e ") != string::npos)     { label = readFile(eefile, header, v_ee); v_driver = v_ee; }
    else if (header.find("\u03BC") != string::npos) { label = readFile(mmfile, header, v_mm); v_driver = v_mm; }
    else                                            { label = readFile(mmfile, header, v_mm); readFile(eefile, header, v_ee); v_driver = v_mm; }
    if (header != "Systematics: e trigger")         readFile(emfile, header, v_em);

    header = "                                              Combined " + header;
    outfile << boldline << endl << header << endl << boldline << endl;
    outfile << label << endl;

    int start_cut, signal_cut; //Correct Channel, = 1 Jet, = 0 btags, metR
    string vert_line;
    if (header.find("Cut Flow") != string::npos) {
      start_cut = 6; signal_cut = 16; vert_line = "|||";
    }
    else {
      start_cut = 12; signal_cut = 32; vert_line = "||";
    }

    if (header.find("Systematics:") != string::npos) { ll_sysUP[header] = samples; ll_sysDN[header] = samples; }
    for (int i=0, n=v_driver.size(); i<n; i++) {
      if (i==start_cut || i==signal_cut) outfile << divider << endl;

      if (header.find("Summary") != string::npos) cutnames.push_back(v_driver[i][0]);
      outfile << v_driver[i][0];
      int nCol = v_driver[i].size();
      for (int j=1; j<nCol; j++) {

        string sample = "";
        if (header.find("Summary") != string::npos) {
          if      (j==1) sample = "data";
          else if (j==2) sample = "bkgd";
          else if (j==3) { outfile << Form("|||       %1.3f       ", m_nom["ll"]["data"][i]/m_nom["ll"]["bkgd"][i]); continue; }
          else if (j==4) sample = "zprime";
          else if (j==5) sample = "gluon";
        }
        else {
          if      (j==1) sample = "ttbar";
          else if (j==2) sample = "dy";
          else if (j==3) sample = "st";
          else if (j==4) sample = "vv";
          else if (j==5) sample = "wjet";
          else if (j==6) sample = "bkgd";
          else if (j==7) sample = "zprime";
          else if (j==8) sample = "gluon";
        }

        if (header.find("Summary") != string::npos || header.find("Background") != string::npos) {
          m_nom["mm"][sample][i] = stod( v_mm[i][j].erase( v_mm[i][j].find('(') ) ); m_nom["ee"][sample][i] = stod( v_ee[i][j].erase( v_ee[i][j].find('(') ) );
          m_nom["em"][sample][i] = stod( v_em[i][j].erase( v_em[i][j].find('(') ) );

          if (i < start_cut) {
            if (sample == "data") m_nom["ll"][sample][i] = m_nom["mm"][sample][i] + m_nom["ee"][sample][i];
            else                  m_nom["ll"][sample][i] = m_nom["mm"][sample][i];
          }
          else                    m_nom["ll"][sample][i] = m_nom["mm"][sample][i] + m_nom["ee"][sample][i] + m_nom["em"][sample][i];

          if (sample == "data") outfile << Form("|||  %12.0f (%1.6f)  ", m_nom["ll"][sample][i], m_nom["ll"][sample][i]/m_nom["ll"][sample][0]);
          else                  outfile << Form("|||  %12.1f (%1.6f)  ", m_nom["ll"][sample][i], m_nom["ll"][sample][i]/m_nom["ll"][sample][0]);
        }

        else if (header.find("Statistical Error") != string::npos) {
          m_stat["mm"][sample][i] = stod( v_mm[i][j] ); m_stat["ee"][sample][i] = stod( v_ee[i][j] ); m_stat["em"][sample][i] = stod( v_em[i][j] );

          if (i < start_cut) m_stat["ll"][sample][i] = m_stat["mm"][sample][i];
          else               m_stat["ll"][sample][i] = sqrt( m_stat["mm"][sample][i]*m_stat["mm"][sample][i] + 
                               m_stat["ee"][sample][i]*m_stat["ee"][sample][i] + m_stat["em"][sample][i]*m_stat["em"][sample][i] );

          outfile << Form("|||  %9.2f  ", m_stat["ll"][sample][i]);
        }

        else if (header.find("Systematics:") != string::npos) {
          if (i%2 == 0) {
            double mm_up=0, mm_dn=0, ee_up=0, ee_dn=0, em_up=0, em_dn=0;
            if (header.find("e ") == string::npos) {
              while (v_mm[i][j].length() > 0 && v_mm[i][j].at(0) == ' ') v_mm[i][j].erase(0, 1);
              mm_up = stod( v_mm[i][j].substr(0, v_mm[i][j].find(' ')) );
              mm_dn = stod( v_mm[i][j].erase(0, v_mm[i][j].find(' ')) );
            }
            if (header.find("\u03BC") == string::npos) {
              while (v_ee[i][j].length() > 0 && v_ee[i][j].at(0) == ' ') v_ee[i][j].erase(0, 1);
              ee_up = stod( v_ee[i][j].substr(0, v_ee[i][j].find(' ')) );
              ee_dn = stod( v_ee[i][j].erase(0, v_ee[i][j].find(' ')) );
            }
            if (header.find("e trigger") == string::npos) {
              while (v_em[i][j].length() > 0 && v_em[i][j].at(0) == ' ') v_em[i][j].erase(0, 1);
              em_up = stod( v_em[i][j].substr(0, v_em[i][j].find(' ')) );
              em_dn = stod( v_em[i][j].erase(0, v_em[i][j].find(' ')) );
            }
            int idx = i/2;
            if (i < start_cut) {
              if (header.find("e ") != string::npos) { ll_sysUP[header][sample][idx] = ee_up; ll_sysDN[header][sample][idx] = ee_dn; }
              else                                   { ll_sysUP[header][sample][idx] = mm_up; ll_sysDN[header][sample][idx] = mm_dn; }
            }
            else {
              ll_sysUP[header][sample][idx] = mm_up + ee_up + em_up;
              ll_sysDN[header][sample][idx] = mm_dn + ee_dn + em_dn;
            }
            outfile << Form("|| %9.2f %9.2f ", ll_sysUP[header][sample][idx], ll_sysDN[header][sample][idx]);

            double ll_d1 = delta(ll_sysUP[header][sample][idx], ll_sysDN[header][sample][idx], "+"), mm_d1 = delta(mm_up, mm_dn, "+"), ee_d1 = delta(ee_up, ee_dn, "+"), em_d1 = delta(em_up, em_dn, "+");
            double ll_d2 = delta(ll_sysUP[header][sample][idx], ll_sysDN[header][sample][idx], "-"), mm_d2 = delta(mm_up, mm_dn, "-"), ee_d2 = delta(ee_up, ee_dn, "-"), em_d2 = delta(em_up, em_dn, "-");

            m_sysUP["ll"][sample][idx] += ll_d1*ll_d1; m_sysUP["mm"][sample][idx] += mm_d1*mm_d1; m_sysUP["ee"][sample][idx] += ee_d1*ee_d1; m_sysUP["em"][sample][idx] += em_d1*em_d1;
            m_sysDN["ll"][sample][idx] += ll_d2*ll_d2; m_sysDN["mm"][sample][idx] += mm_d2*mm_d2; m_sysDN["ee"][sample][idx] += ee_d2*ee_d2; m_sysDN["em"][sample][idx] += em_d2*em_d2;
          }
          else {
            int idx = (i-1)/2;
            if (m_nom["ll"][sample][idx] == 0) outfile << Form("|| %9.2f %9.2f ", 0., 0.);
            else                               outfile << Form("|| %9.2f %9.2f ", ll_sysUP[header][sample][idx]/m_nom["ll"][sample][idx]*100, ll_sysDN[header][sample][idx]/m_nom["ll"][sample][idx]*100);
          }
        }
      }
      outfile << endl;
    }
  outfile << endl;
  }
/*
  for (auto const& it_chan : m_nom) {
    string chan = it_chan.first;
    map<string, vector<double> > mm_nom = it_chan.second;

    for (auto const& it_sample : mm_nom) {
      string sample = it_sample.first;
      int n = it_sample.second.size();

      cout << chan << " " << sample << endl;
      for (int i=0; i<n; i++)        
        cout << cutnames[i] << m_nom[chan][sample][i] << endl;
      cout << endl;
    }
  }
*/
  /////////////////////////////////////////
  ///////////TOTAL SYSTEMATICS/////////////
  /////////////////////////////////////////

  outfile << boldline << "\n                                    Combined Total Systematics:\n" << boldline << endl;
  outfile << "                           ||        ttbar        ||      Drell-Yan      ||     Single-Top      ||       Diboson       ||       W+Jets        ||      background     || zprime_M-3000_W-300 || gluon_M-3000" << endl;

  for (int i=0, n=cutnames.size(); i<n; i++) {
    if (i==6 || i==16) outfile << divider << endl;
    cutnames[i].pop_back();

    outfile << cutnames[i] << Form("|| %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f ",
      sqrt(m_sysUP["ll"]["ttbar"][i]),  sqrt(m_sysDN["ll"]["ttbar"][i]),  sqrt(m_sysUP["ll"]["dy"][i]),    sqrt(m_sysDN["ll"]["dy"][i]),   sqrt(m_sysUP["ll"]["st"][i]), sqrt(m_sysDN["ll"]["st"][i]),
      sqrt(m_sysUP["ll"]["vv"][i]),     sqrt(m_sysDN["ll"]["vv"][i]),     sqrt(m_sysUP["ll"]["wjet"][i]),  sqrt(m_sysDN["ll"]["wjet"][i]), sqrt(m_sysUP["ll"]["bkgd"][i]), sqrt(m_sysDN["ll"]["bkgd"][i]),
      sqrt(m_sysUP["ll"]["zprime"][i]), sqrt(m_sysDN["ll"]["zprime"][i]), sqrt(m_sysUP["ll"]["gluon"][i]), sqrt(m_sysDN["ll"]["gluon"][i])
    ) << endl;

    outfile << Form("                           || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f ",
      sqrt(m_sysUP["ll"]["ttbar"][i])/m_nom["ll"]["ttbar"][i]*100,   sqrt(m_sysDN["ll"]["ttbar"][i])/m_nom["ll"]["ttbar"][i]*100,
      sqrt(m_sysUP["ll"]["dy"][i])/m_nom["ll"]["dy"][i]*100,         sqrt(m_sysDN["ll"]["dy"][i])/m_nom["ll"]["dy"][i]*100,
      sqrt(m_sysUP["ll"]["st"][i])/m_nom["ll"]["st"][i]*100,         sqrt(m_sysDN["ll"]["st"][i])/m_nom["ll"]["st"][i]*100,
      sqrt(m_sysUP["ll"]["vv"][i])/m_nom["ll"]["vv"][i]*100,         sqrt(m_sysDN["ll"]["vv"][i])/m_nom["ll"]["vv"][i]*100,
      sqrt(m_sysUP["ll"]["wjet"][i])/m_nom["ll"]["wjet"][i]*100,     sqrt(m_sysDN["ll"]["wjet"][i])/m_nom["ll"]["wjet"][i]*100,
      sqrt(m_sysUP["ll"]["bkgd"][i])/m_nom["ll"]["bkgd"][i]*100,     sqrt(m_sysDN["ll"]["bkgd"][i])/m_nom["ll"]["bkgd"][i]*100,
      sqrt(m_sysUP["ll"]["zprime"][i])/m_nom["ll"]["zprime"][i]*100, sqrt(m_sysDN["ll"]["zprime"][i])/m_nom["ll"]["zprime"][i]*100,
      sqrt(m_sysUP["ll"]["gluon"][i])/m_nom["ll"]["gluon"][i]*100,   sqrt(m_sysDN["ll"]["gluon"][i])/m_nom["ll"]["gluon"][i]*100
    ) << endl;
  }
  outfile << endl;

  int icut = 23;
  outfile << boldline << "\n           Combined Total Systematics: " << cutnames[icut] << endl << boldline << endl;
  outfile << "                           ||        ttbar        ||      Drell-Yan      ||     Single-Top      ||       Diboson       ||       W+Jets        ||      background     || zprime_M-3000_W-300 || gluon_M-3000" << endl;

  for (auto const& it_sys : ll_sysUP) {
    map<string, vector<double> > m_UP = it_sys.second, m_DN = ll_sysDN[it_sys.first];
    TString sysname = it_sys.first;
    sysname.ReplaceAll("                                              Combined Systematics: ", "");

    outfile << Form("%-27s|| %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f ", sysname.Data(),
      m_UP["ttbar"][icut], m_DN["ttbar"][icut], m_UP["dy"][icut],   m_DN["dy"][icut],   m_UP["st"][icut],     m_DN["st"][icut],     m_UP["vv"][icut],    m_DN["vv"][icut],
      m_UP["wjet"][icut],  m_DN["wjet"][icut],  m_UP["bkgd"][icut], m_DN["bkgd"][icut], m_UP["zprime"][icut], m_DN["zprime"][icut], m_UP["gluon"][icut], m_DN["gluon"][icut]
    ) << endl;

    outfile << Form("                           || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f ",
      m_UP["ttbar"][icut]/m_nom["ll"]["ttbar"][icut]*100,   m_DN["ttbar"][icut]/m_nom["ll"]["ttbar"][icut]*100,
      m_UP["dy"][icut]/m_nom["ll"]["dy"][icut]*100,         m_DN["dy"][icut]/m_nom["ll"]["dy"][icut]*100,
      m_UP["st"][icut]/m_nom["ll"]["st"][icut]*100,         m_DN["st"][icut]/m_nom["ll"]["st"][icut]*100,
      m_UP["vv"][icut]/m_nom["ll"]["vv"][icut]*100,         m_DN["vv"][icut]/m_nom["ll"]["vv"][icut]*100,
      m_UP["wjet"][icut]/m_nom["ll"]["wjet"][icut]*100,     m_DN["wjet"][icut]/m_nom["ll"]["wjet"][icut]*100,
      m_UP["bkgd"][icut]/m_nom["ll"]["bkgd"][icut]*100,     m_DN["bkgd"][icut]/m_nom["ll"]["bkgd"][icut]*100,
      m_UP["zprime"][icut]/m_nom["ll"]["zprime"][icut]*100, m_DN["zprime"][icut]/m_nom["ll"]["zprime"][icut]*100,
      m_UP["gluon"][icut]/m_nom["ll"]["gluon"][icut]*100,   m_DN["gluon"][icut]/m_nom["ll"]["gluon"][icut]*100
    ) << endl;
  }
  outfile << Form("%-27s|| %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f ", "Total Sys",
    sqrt(m_sysUP["ll"]["ttbar"][icut]),  sqrt(m_sysDN["ll"]["ttbar"][icut]),  sqrt(m_sysUP["ll"]["dy"][icut]),    sqrt(m_sysDN["ll"]["dy"][icut]),   sqrt(m_sysUP["ll"]["st"][icut]), sqrt(m_sysDN["ll"]["st"][icut]),
    sqrt(m_sysUP["ll"]["vv"][icut]),     sqrt(m_sysDN["ll"]["vv"][icut]),     sqrt(m_sysUP["ll"]["wjet"][icut]),  sqrt(m_sysDN["ll"]["wjet"][icut]), sqrt(m_sysUP["ll"]["bkgd"][icut]), sqrt(m_sysDN["ll"]["bkgd"][icut]),
    sqrt(m_sysUP["ll"]["zprime"][icut]), sqrt(m_sysDN["ll"]["zprime"][icut]), sqrt(m_sysUP["ll"]["gluon"][icut]), sqrt(m_sysDN["ll"]["gluon"][icut])
  ) << endl;
  outfile << Form("                           || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f ",
    sqrt(m_sysUP["ll"]["ttbar"][icut])/m_nom["ll"]["ttbar"][icut]*100,   sqrt(m_sysDN["ll"]["ttbar"][icut])/m_nom["ll"]["ttbar"][icut]*100,
    sqrt(m_sysUP["ll"]["dy"][icut])/m_nom["ll"]["dy"][icut]*100,         sqrt(m_sysDN["ll"]["dy"][icut])/m_nom["ll"]["dy"][icut]*100,
    sqrt(m_sysUP["ll"]["st"][icut])/m_nom["ll"]["st"][icut]*100,         sqrt(m_sysDN["ll"]["st"][icut])/m_nom["ll"]["st"][icut]*100,
    sqrt(m_sysUP["ll"]["vv"][icut])/m_nom["ll"]["vv"][icut]*100,         sqrt(m_sysDN["ll"]["vv"][icut])/m_nom["ll"]["vv"][icut]*100,
    sqrt(m_sysUP["ll"]["wjet"][icut])/m_nom["ll"]["wjet"][icut]*100,     sqrt(m_sysDN["ll"]["wjet"][icut])/m_nom["ll"]["wjet"][icut]*100,
    sqrt(m_sysUP["ll"]["bkgd"][icut])/m_nom["ll"]["bkgd"][icut]*100,     sqrt(m_sysDN["ll"]["bkgd"][icut])/m_nom["ll"]["bkgd"][icut]*100,
    sqrt(m_sysUP["ll"]["zprime"][icut])/m_nom["ll"]["zprime"][icut]*100, sqrt(m_sysDN["ll"]["zprime"][icut])/m_nom["ll"]["zprime"][icut]*100,
    sqrt(m_sysUP["ll"]["gluon"][icut])/m_nom["ll"]["gluon"][icut]*100,   sqrt(m_sysDN["ll"]["gluon"][icut])/m_nom["ll"]["gluon"][icut]*100
  ) << endl;
  outfile << endl;

  ///////////////////////////////////////// 
  ///////////      LATEX      /////////////
  /////////////////////////////////////////

  double sig_zprime = 0.272788, sig_gluon = 0.16757;
  for (auto const& it_chan : m_nom) {
    string chan = it_chan.first;
    for (int i=0; i<24; i++) {
      m_nom[chan]["zprime"][i] *= sig_zprime; m_stat[chan]["zprime"][i] *= sig_zprime; m_sysUP[chan]["zprime"][i] *= sig_zprime; m_sysDN[chan]["zprime"][i] *= sig_zprime;
      m_nom[chan]["gluon"][i] *= sig_gluon; m_stat[chan]["gluon"][i] *= sig_gluon; m_sysUP[chan]["gluon"][i] *= sig_gluon; m_sysDN[chan]["gluon"][i] *= sig_gluon;
    }
  }

  for (int i=16; i<24; i++) {
    TString cut = cutnames[i];
    cut.ReplaceAll(">=", "$\\geq$");

  outfile << "\\begin{tabular}{ |c|c|c|c|c| }\n";
  outfile << "\\multicolumn{5}{c}{" << cut << "} \\\\\n";
  outfile << "\\hline\n";
  outfile << "Sample & $\\mu\\mu$ Channel & ee Channel & e$\\mu$ Channel & Combined \\\\\n";
  outfile << "\\hline\n";

  outfile << "Z$'$ (10\\%, 3 TeV) &" << Form("$%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ \\\\\n",
    m_nom["mm"]["zprime"][i], sqrt(m_sysUP["mm"]["zprime"][i] + m_stat["mm"]["zprime"][i]*m_stat["mm"]["zprime"][i]), sqrt(m_sysDN["mm"]["zprime"][i] + m_stat["mm"]["zprime"][i]*m_stat["mm"]["zprime"][i]),
    m_nom["ee"]["zprime"][i], sqrt(m_sysUP["ee"]["zprime"][i] + m_stat["ee"]["zprime"][i]*m_stat["ee"]["zprime"][i]), sqrt(m_sysDN["ee"]["zprime"][i] + m_stat["ee"]["zprime"][i]*m_stat["ee"]["zprime"][i]),
    m_nom["em"]["zprime"][i], sqrt(m_sysUP["em"]["zprime"][i] + m_stat["em"]["zprime"][i]*m_stat["em"]["zprime"][i]), sqrt(m_sysDN["em"]["zprime"][i] + m_stat["em"]["zprime"][i]*m_stat["em"]["zprime"][i]),
    m_nom["ll"]["zprime"][i], sqrt(m_sysUP["ll"]["zprime"][i] + m_stat["ll"]["zprime"][i]*m_stat["ll"]["zprime"][i]), sqrt(m_sysDN["ll"]["zprime"][i] + m_stat["ll"]["zprime"][i]*m_stat["ll"]["zprime"][i])
  ) << endl;
  outfile << "$\\textrm{g}_{\\textrm{kk}}$ (3 TeV) &" << Form("$%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ \\\\\n",
    m_nom["mm"]["gluon"][i], sqrt(m_sysUP["mm"]["gluon"][i] + m_stat["mm"]["gluon"][i]*m_stat["mm"]["gluon"][i]), sqrt(m_sysDN["mm"]["gluon"][i] + m_stat["mm"]["gluon"][i]*m_stat["mm"]["gluon"][i]),
    m_nom["ee"]["gluon"][i], sqrt(m_sysUP["ee"]["gluon"][i] + m_stat["ee"]["gluon"][i]*m_stat["ee"]["gluon"][i]), sqrt(m_sysDN["ee"]["gluon"][i] + m_stat["ee"]["gluon"][i]*m_stat["ee"]["gluon"][i]),
    m_nom["em"]["gluon"][i], sqrt(m_sysUP["em"]["gluon"][i] + m_stat["em"]["gluon"][i]*m_stat["em"]["gluon"][i]), sqrt(m_sysDN["em"]["gluon"][i] + m_stat["em"]["gluon"][i]*m_stat["em"]["gluon"][i]),
    m_nom["ll"]["gluon"][i], sqrt(m_sysUP["ll"]["gluon"][i] + m_stat["ll"]["gluon"][i]*m_stat["ll"]["gluon"][i]), sqrt(m_sysDN["ll"]["gluon"][i] + m_stat["ll"]["gluon"][i]*m_stat["ll"]["gluon"][i])
  ) << endl;
  outfile << "\\hline\n";

  outfile << "t$\\bar{\\textrm{t}}$ &" << Form("$%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ \\\\\n",
    m_nom["mm"]["ttbar"][i], sqrt(m_sysUP["mm"]["ttbar"][i] + m_stat["mm"]["ttbar"][i]*m_stat["mm"]["ttbar"][i]), sqrt(m_sysDN["mm"]["ttbar"][i] + m_stat["mm"]["ttbar"][i]*m_stat["mm"]["ttbar"][i]),
    m_nom["ee"]["ttbar"][i], sqrt(m_sysUP["ee"]["ttbar"][i] + m_stat["ee"]["ttbar"][i]*m_stat["ee"]["ttbar"][i]), sqrt(m_sysDN["ee"]["ttbar"][i] + m_stat["ee"]["ttbar"][i]*m_stat["ee"]["ttbar"][i]),
    m_nom["em"]["ttbar"][i], sqrt(m_sysUP["em"]["ttbar"][i] + m_stat["em"]["ttbar"][i]*m_stat["em"]["ttbar"][i]), sqrt(m_sysDN["em"]["ttbar"][i] + m_stat["em"]["ttbar"][i]*m_stat["em"]["ttbar"][i]),
    m_nom["ll"]["ttbar"][i], sqrt(m_sysUP["ll"]["ttbar"][i] + m_stat["ll"]["ttbar"][i]*m_stat["ll"]["ttbar"][i]), sqrt(m_sysDN["ll"]["ttbar"][i] + m_stat["ll"]["ttbar"][i]*m_stat["ll"]["ttbar"][i])
  ) << endl;
  outfile << "Z/$\\gamma^{*}\\rightarrow l^{+}l^{-}$ &" << Form("$%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ \\\\\n",
    m_nom["mm"]["dy"][i], sqrt(m_sysUP["mm"]["dy"][i] + m_stat["mm"]["dy"][i]*m_stat["mm"]["dy"][i]), sqrt(m_sysDN["mm"]["dy"][i] + m_stat["mm"]["dy"][i]*m_stat["mm"]["dy"][i]),
    m_nom["ee"]["dy"][i], sqrt(m_sysUP["ee"]["dy"][i] + m_stat["ee"]["dy"][i]*m_stat["ee"]["dy"][i]), sqrt(m_sysDN["ee"]["dy"][i] + m_stat["ee"]["dy"][i]*m_stat["ee"]["dy"][i]),
    m_nom["em"]["dy"][i], sqrt(m_sysUP["em"]["dy"][i] + m_stat["em"]["dy"][i]*m_stat["em"]["dy"][i]), sqrt(m_sysDN["em"]["dy"][i] + m_stat["em"]["dy"][i]*m_stat["em"]["dy"][i]),
    m_nom["ll"]["dy"][i], sqrt(m_sysUP["ll"]["dy"][i] + m_stat["ll"]["dy"][i]*m_stat["ll"]["dy"][i]), sqrt(m_sysDN["ll"]["dy"][i] + m_stat["ll"]["dy"][i]*m_stat["ll"]["dy"][i])
  ) << endl;
  outfile << "Single-Top &" << Form("$%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ \\\\\n",
    m_nom["mm"]["st"][i], sqrt(m_sysUP["mm"]["st"][i] + m_stat["mm"]["st"][i]*m_stat["mm"]["st"][i]), sqrt(m_sysDN["mm"]["st"][i] + m_stat["mm"]["st"][i]*m_stat["mm"]["st"][i]),
    m_nom["ee"]["st"][i], sqrt(m_sysUP["ee"]["st"][i] + m_stat["ee"]["st"][i]*m_stat["ee"]["st"][i]), sqrt(m_sysDN["ee"]["st"][i] + m_stat["ee"]["st"][i]*m_stat["ee"]["st"][i]),
    m_nom["em"]["st"][i], sqrt(m_sysUP["em"]["st"][i] + m_stat["em"]["st"][i]*m_stat["em"]["st"][i]), sqrt(m_sysDN["em"]["st"][i] + m_stat["em"]["st"][i]*m_stat["em"]["st"][i]),
    m_nom["ll"]["st"][i], sqrt(m_sysUP["ll"]["st"][i] + m_stat["ll"]["st"][i]*m_stat["ll"]["st"][i]), sqrt(m_sysDN["ll"]["st"][i] + m_stat["ll"]["st"][i]*m_stat["ll"]["st"][i])
  ) << endl;
  outfile << "VV &" << Form("$%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ \\\\\n",
    m_nom["mm"]["vv"][i], sqrt(m_sysUP["mm"]["vv"][i] + m_stat["mm"]["vv"][i]*m_stat["mm"]["vv"][i]), sqrt(m_sysDN["mm"]["vv"][i] + m_stat["mm"]["vv"][i]*m_stat["mm"]["vv"][i]),
    m_nom["ee"]["vv"][i], sqrt(m_sysUP["ee"]["vv"][i] + m_stat["ee"]["vv"][i]*m_stat["ee"]["vv"][i]), sqrt(m_sysDN["ee"]["vv"][i] + m_stat["ee"]["vv"][i]*m_stat["ee"]["vv"][i]),
    m_nom["em"]["vv"][i], sqrt(m_sysUP["em"]["vv"][i] + m_stat["em"]["vv"][i]*m_stat["em"]["vv"][i]), sqrt(m_sysDN["em"]["vv"][i] + m_stat["em"]["vv"][i]*m_stat["em"]["vv"][i]),
    m_nom["ll"]["vv"][i], sqrt(m_sysUP["ll"]["vv"][i] + m_stat["ll"]["vv"][i]*m_stat["ll"]["vv"][i]), sqrt(m_sysDN["ll"]["vv"][i] + m_stat["ll"]["vv"][i]*m_stat["ll"]["vv"][i])
  ) << endl;
  outfile << "W+Jets &" << Form("$%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ \\\\\n",
    m_nom["mm"]["wjet"][i], sqrt(m_sysUP["mm"]["wjet"][i] + m_stat["mm"]["wjet"][i]*m_stat["mm"]["wjet"][i]), sqrt(m_sysDN["mm"]["wjet"][i] + m_stat["mm"]["wjet"][i]*m_stat["mm"]["wjet"][i]),
    m_nom["ee"]["wjet"][i], sqrt(m_sysUP["ee"]["wjet"][i] + m_stat["ee"]["wjet"][i]*m_stat["ee"]["wjet"][i]), sqrt(m_sysDN["ee"]["wjet"][i] + m_stat["ee"]["wjet"][i]*m_stat["ee"]["wjet"][i]),
    m_nom["em"]["wjet"][i], sqrt(m_sysUP["em"]["wjet"][i] + m_stat["em"]["wjet"][i]*m_stat["em"]["wjet"][i]), sqrt(m_sysDN["em"]["wjet"][i] + m_stat["em"]["wjet"][i]*m_stat["em"]["wjet"][i]),
    m_nom["ll"]["wjet"][i], sqrt(m_sysUP["ll"]["wjet"][i] + m_stat["ll"]["wjet"][i]*m_stat["ll"]["wjet"][i]), sqrt(m_sysDN["ll"]["wjet"][i] + m_stat["ll"]["wjet"][i]*m_stat["ll"]["wjet"][i])
  ) << endl;
  outfile << "\\hline\n";

  double mm_bkgd_errorUP = sqrt(m_sysUP["mm"]["bkgd"][i] + m_stat["mm"]["bkgd"][i]*m_stat["mm"]["bkgd"][i]), mm_bkgd_errorDN = sqrt(m_sysDN["mm"]["bkgd"][i] + m_stat["mm"]["bkgd"][i]*m_stat["mm"]["bkgd"][i]);
  double ee_bkgd_errorUP = sqrt(m_sysUP["ee"]["bkgd"][i] + m_stat["ee"]["bkgd"][i]*m_stat["ee"]["bkgd"][i]), ee_bkgd_errorDN = sqrt(m_sysDN["ee"]["bkgd"][i] + m_stat["ee"]["bkgd"][i]*m_stat["ee"]["bkgd"][i]);
  double em_bkgd_errorUP = sqrt(m_sysUP["em"]["bkgd"][i] + m_stat["em"]["bkgd"][i]*m_stat["em"]["bkgd"][i]), em_bkgd_errorDN = sqrt(m_sysDN["em"]["bkgd"][i] + m_stat["em"]["bkgd"][i]*m_stat["em"]["bkgd"][i]);
  double ll_bkgd_errorUP = sqrt(m_sysUP["ll"]["bkgd"][i] + m_stat["ll"]["bkgd"][i]*m_stat["ll"]["bkgd"][i]), ll_bkgd_errorDN = sqrt(m_sysDN["ll"]["bkgd"][i] + m_stat["ll"]["bkgd"][i]*m_stat["ll"]["bkgd"][i]);

  outfile << "Total Bkg &" << Form("$%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ & $%.1f^{+%.1f}_{-%.1f}$ \\\\\n",
    m_nom["mm"]["bkgd"][i], mm_bkgd_errorUP, mm_bkgd_errorDN,
    m_nom["ee"]["bkgd"][i], ee_bkgd_errorUP, ee_bkgd_errorDN,
    m_nom["em"]["bkgd"][i], em_bkgd_errorUP, em_bkgd_errorDN,
    m_nom["ll"]["bkgd"][i], ll_bkgd_errorUP, ll_bkgd_errorDN
  ) << endl;
  outfile << "Data &" << Form("$%.0f \\pm %.1f$ & $%.0f \\pm %.1f$ & $%.0f \\pm %.1f$ & $%.0f \\pm %.1f$ \\\\\n",
    m_nom["mm"]["data"][i], sqrt(m_nom["mm"]["data"][i]),
    m_nom["ee"]["data"][i], sqrt(m_nom["ee"]["data"][i]),
    m_nom["em"]["data"][i], sqrt(m_nom["em"]["data"][i]),
    m_nom["ll"]["data"][i], sqrt(m_nom["ll"]["data"][i])
  ) << endl;
  outfile << "Data/Bkg &" << Form("$%.2f^{+%.2f}_{-%.2f}$ & $%.2f^{+%.2f}_{-%.2f}$ & $%.2f^{+%.2f}_{-%.2f}$ & $%.2f^{+%.2f}_{-%.2f}$ \\\\\n",
    m_nom["mm"]["data"][i]/m_nom["mm"]["bkgd"][i],
    m_nom["mm"]["data"][i]/m_nom["mm"]["bkgd"][i] * sqrt( 1/m_nom["mm"]["data"][i] + mm_bkgd_errorUP*mm_bkgd_errorUP/(m_nom["mm"]["bkgd"][i]*m_nom["mm"]["bkgd"][i]) ),
    m_nom["mm"]["data"][i]/m_nom["mm"]["bkgd"][i] * sqrt( 1/m_nom["mm"]["data"][i] + mm_bkgd_errorDN*mm_bkgd_errorDN/(m_nom["mm"]["bkgd"][i]*m_nom["mm"]["bkgd"][i]) ),
    m_nom["ee"]["data"][i]/m_nom["ee"]["bkgd"][i],
    m_nom["ee"]["data"][i]/m_nom["ee"]["bkgd"][i] * sqrt( 1/m_nom["ee"]["data"][i] + ee_bkgd_errorUP*ee_bkgd_errorUP/(m_nom["ee"]["bkgd"][i]*m_nom["ee"]["bkgd"][i]) ),
    m_nom["ee"]["data"][i]/m_nom["ee"]["bkgd"][i] * sqrt( 1/m_nom["ee"]["data"][i] + ee_bkgd_errorDN*ee_bkgd_errorDN/(m_nom["ee"]["bkgd"][i]*m_nom["ee"]["bkgd"][i]) ),
    m_nom["em"]["data"][i]/m_nom["em"]["bkgd"][i],
    m_nom["em"]["data"][i]/m_nom["em"]["bkgd"][i] * sqrt( 1/m_nom["em"]["data"][i] + em_bkgd_errorUP*em_bkgd_errorUP/(m_nom["em"]["bkgd"][i]*m_nom["em"]["bkgd"][i]) ),
    m_nom["em"]["data"][i]/m_nom["em"]["bkgd"][i] * sqrt( 1/m_nom["em"]["data"][i] + em_bkgd_errorDN*em_bkgd_errorDN/(m_nom["em"]["bkgd"][i]*m_nom["em"]["bkgd"][i]) ),
    m_nom["ll"]["data"][i]/m_nom["ll"]["bkgd"][i],
    m_nom["ll"]["data"][i]/m_nom["ll"]["bkgd"][i] * sqrt( 1/m_nom["ll"]["data"][i] + ll_bkgd_errorUP*ll_bkgd_errorUP/(m_nom["ll"]["bkgd"][i]*m_nom["ll"]["bkgd"][i]) ),
    m_nom["ll"]["data"][i]/m_nom["ll"]["bkgd"][i] * sqrt( 1/m_nom["ll"]["data"][i] + ll_bkgd_errorDN*ll_bkgd_errorDN/(m_nom["ll"]["bkgd"][i]*m_nom["ll"]["bkgd"][i]) )
  ) << endl;
  outfile << "\\hline\n";

  outfile << "S(Z$'$)/Bkg &" << Form("%.5f & %.5f & %.5f & %.5f \\\\\n",
    m_nom["mm"]["zprime"][i]/m_nom["mm"]["bkgd"][i], m_nom["ee"]["zprime"][i]/m_nom["ee"]["bkgd"][i], m_nom["em"]["zprime"][i]/m_nom["em"]["bkgd"][i], m_nom["ll"]["zprime"][i]/m_nom["ll"]["bkgd"][i]
  ) << endl;
  outfile << "S($\\textrm{g}_{\\textrm{kk}}$)/Bkg &" << Form("%.5f & %.5f & %.5f & %.5f \\\\\n",
    m_nom["mm"]["gluon"][i]/m_nom["mm"]["bkgd"][i], m_nom["ee"]["gluon"][i]/m_nom["ee"]["bkgd"][i], m_nom["em"]["gluon"][i]/m_nom["em"]["bkgd"][i], m_nom["ll"]["gluon"][i]/m_nom["ll"]["bkgd"][i]
  ) << endl;
  outfile << "\\hline\n";

  outfile << "S(Z$'$)/$\\sqrt{\\textrm{S(Z$'$) + Bkg}}$ &" << Form("%.3f & %.3f & %.3f & %.3f \\\\\n",
    m_nom["mm"]["zprime"][i]/sqrt(m_nom["mm"]["zprime"][i] + m_nom["mm"]["bkgd"][i]), m_nom["ee"]["zprime"][i]/sqrt(m_nom["ee"]["zprime"][i] + m_nom["ee"]["bkgd"][i]),
    m_nom["em"]["zprime"][i]/sqrt(m_nom["em"]["zprime"][i] + m_nom["em"]["bkgd"][i]), m_nom["ll"]["zprime"][i]/sqrt(m_nom["ll"]["zprime"][i] + m_nom["ll"]["bkgd"][i])
  ) << endl;
  outfile << "S($\\textrm{g}_{\\textrm{kk}})/\\sqrt{\\textrm{S(g}_{\\textrm{kk}}\\textrm{) + Bkg}}$ &" << Form("%.3f & %.3f & %.3f & %.3f \\\\\n",
    m_nom["mm"]["gluon"][i]/sqrt(m_nom["mm"]["gluon"][i] + m_nom["mm"]["bkgd"][i]), m_nom["ee"]["gluon"][i]/sqrt(m_nom["ee"]["gluon"][i] + m_nom["ee"]["bkgd"][i]),
    m_nom["em"]["gluon"][i]/sqrt(m_nom["em"]["gluon"][i] + m_nom["em"]["bkgd"][i]), m_nom["ll"]["gluon"][i]/sqrt(m_nom["ll"]["gluon"][i] + m_nom["ll"]["bkgd"][i])
  ) << endl;
  outfile << "\\hline\n";

  outfile << "\\end{tabular}\n" << endl;
  }

  outfile.close();
  mmfile.close(); eefile.close(); emfile.close();
}

//read file from start until no line. Store contents in vec
//returns column labels
string readFile(istream& file, const string& start, vector< vector<string> >& vec) {
  vec.clear();
  string line, label;

  while (getline(file, line)) { if (line.find(start) != string::npos) break; }
  if (file.eof()) { cout << "END OF FILE: " << start << endl; file.clear(); file.seekg(0, ios::beg); return ""; }
  getline(file, line); //boldline
  getline(file, line); label = line;

  while (getline(file, line)) {
    if (line.length() == 0) return label;
    if (line.at(0) == '-') continue;

    vector<string> subvec;
    int delim_pos;
    while ( (delim_pos = line.find('|')) != -1) {
      subvec.push_back( line.substr(0, delim_pos) );

      line.erase(0, delim_pos+1);
      while (line.length() > 0 && line.at(0) == '|') line.erase(0, 1);
    }
    if (line.length() > 0) subvec.push_back( line.substr(0, line.length()) );
    vec.push_back(subvec);
  }
  return label;
}

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

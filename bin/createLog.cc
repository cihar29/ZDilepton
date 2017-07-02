#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <iomanip>
#include <string>

#include "TString.h"

using namespace std;

void readFile(const string& fileName, vector<pair<string, map<TString, pair<double, double> > > >& cuts, string& channel);

double delta(double delta1, double delta2, string pos_neg);

int main(int argc, char* argv[]) {

  if (argc == 1 || argc == 2) { cout << "Please provide two log files and systematics (e.g. createLog logData.txt logMC.txt jec)." << endl; return -1; }

  string dataFile = argv[1];
  string mcFile = argv[2];
  string channel = "";

  //use vector to preserve order of cuts
  //vector(cut_name, map(dataset, (N, error) ) )
  vector<pair<string, map<TString, pair<double, double> > > > cuts;

  //read data first to initialize mc
  readFile(dataFile, cuts, channel);
  readFile(mcFile, cuts, channel);

  TString zprime="", gluon="";
  map<TString, pair<double, double> >& m_total = cuts[0].second;
  for (map<TString, pair<double, double> >::iterator i_set = m_total.begin(); i_set != m_total.end(); ++i_set) {
    TString dataset = i_set->first;
    if (dataset.Contains("zprime", TString::kIgnoreCase)) zprime = dataset;
    else if (dataset.Contains("gluon", TString::kIgnoreCase)) gluon = dataset;
  }

  ofstream file( channel + "_cutflow.txt" );

  file<<"====================================================================================================================="<< "\n" ;
  file<<"                                              Cut Flow Table: Summary\n" ;
  file<<"====================================================================================================================="<< "\n" ;
  file<< Form("                          |||               Data                |||             Background            |||        %-20s       |||        %-20s", zprime.Data(), gluon.Data() ) << endl;

  for (vector<pair<string, map<TString, pair<double, double> > > >::iterator i_cut = cuts.begin(); i_cut != cuts.end(); ++i_cut) {
    string cutname = i_cut->first;
    file<< Form("%-25s |||      %12.0f (%1.6f)      |||      %12.1f (%1.6f)      |||      %12.1f (%1.6f)      |||      %12.1f (%1.6f)",
                cutname.data(), i_cut->second["data"].first, i_cut->second["data"].first/m_total["data"].first,
                i_cut->second["background"].first, i_cut->second["background"].first/m_total["background"].first,
                i_cut->second[zprime].first, i_cut->second[zprime].first/m_total[zprime].first,
                i_cut->second[gluon].first, i_cut->second[gluon].first/m_total[gluon].first ) << endl;

    if (cutname == "MET Filters" || cutname == ">= 1 jet")
      file << "---------------------------------------------------------------------------------------------------------------------" << endl;
  }

  file<<"\n====================================================================================================================="<< "\n" ;
  file<<"                                              Cut Flow Table: Background\n" ;
  file<<"====================================================================================================================="<< "\n" ;
  file<<"                          |||              ttbar                |||             Drell-Yan             |||           Single-Top              |||           W+Jets" << endl;

  for (vector<pair<string, map<TString, pair<double, double> > > >::iterator i_cut = cuts.begin(); i_cut != cuts.end(); ++i_cut) {
    string cutname = i_cut->first;
    file<< Form("%-25s |||      %12.1f (%1.6f)      |||      %12.1f (%1.6f)      |||      %12.1f (%1.6f)      |||      %12.1f (%1.6f)",
                cutname.data(), i_cut->second["ttbar"].first, i_cut->second["ttbar"].first/m_total["ttbar"].first,
                i_cut->second["Drell-Yan"].first, i_cut->second["Drell-Yan"].first/m_total["Drell-Yan"].first,
                i_cut->second["Single-Top"].first, i_cut->second["Single-Top"].first/m_total["Single-Top"].first,
                i_cut->second["W+Jets"].first, i_cut->second["W+Jets"].first/m_total["W+Jets"].first ) << endl;

    if (cutname == "MET Filters" || cutname == ">= 1 jet")
      file << "---------------------------------------------------------------------------------------------------------------------" << endl;
  }

  //SYSTEMATICS//

  //which systematics sources to consider: read from command line
  vector<string> systematics ;
  for (int i=3; i < argc; i++){
    string sys = argv[i] ;
    systematics.push_back(sys);
  }

  //use maps that holds vector for each systematics sources. Two maps for UP and DOWN variation for a given sys. source.
  //map(systematics_source, vector(cut_name, map(dataset, (N, error) ) )
  map <string, vector<pair<string, map<TString, pair<double, double> > > > >  m_cutsUP;
  map <string, vector<pair<string, map<TString, pair<double, double> > > > >  m_cutsDOWN;
  for (unsigned int i_sys = 0; i_sys != systematics.size(); ++i_sys) {
    string sys = systematics[i_sys];

    //must have the form log_channel_sysUP.txt
    string mcFile_sysUP = mcFile.substr(0, mcFile.find(".txt"))+"_"+sys+"UP.txt" ;
    vector<pair<string, map<TString, pair<double, double> > > > cutsUP;

    readFile(dataFile, cutsUP, channel);
    readFile(mcFile_sysUP, cutsUP, channel);
    m_cutsUP[sys] = cutsUP ;

    //must have the form log_channel_sysDOWN.txt
    string mcFile_sysDOWN = mcFile.substr(0, mcFile.find(".txt"))+"_"+sys+"DOWN.txt" ;
    vector<pair<string, map<TString, pair<double, double> > > > cutsDOWN;

    readFile(dataFile, cutsDOWN, channel);
    readFile(mcFile_sysDOWN, cutsDOWN, channel);
    m_cutsDOWN[sys] = cutsDOWN ;
  }

  //use vector that will hold total systematics (from all sources added in quadrature)
  //vector(cut_name, map(dataset, (+delta_total, -delta_total ) )
  vector<pair<string, map<TString, pair<double, double> > > > cuts_TotalSysPlusMinus = cuts;
  for (vector<pair<string, map<TString, pair<double, double> > > >::iterator i_cut = cuts_TotalSysPlusMinus.begin(); i_cut != cuts_TotalSysPlusMinus.end(); ++i_cut) {

    map<TString, pair<double, double> > &m = i_cut->second;
    for (map<TString, pair<double, double> >::iterator i_set = m.begin(); i_set != m.end(); ++i_set) {

      TString dataset = i_set->first;
      m[dataset].first = 0. ;  // initialize to zero
      m[dataset].second = 0. ; // initialize to zero
    }
  }

  for (unsigned int i_sys = 0; i_sys != systematics.size(); ++i_sys) { // Loop over systematics sources
    string sys = systematics[i_sys];

    file<<"\n===================================================================================================================================="<< "\n" ;
    file<<"                                    Systematics: " << sys <<  "       [absolute UP DOWN on first line and relative UP DOWN in \% on second line] \n" ;
    file<<"===================================================================================================================================="<< "\n" ;
    file<<Form("                         ||        ttbar        ||      Drell-Yan      ||     Single-Top      ||        W+Jets       ||      background     || %-20s|| %-20s",
    zprime.Data(), gluon.Data() ) << endl;

    vector<pair<string, map<TString, pair<double, double> > > > &cutsUP   = m_cutsUP[sys];
    vector<pair<string, map<TString, pair<double, double> > > > &cutsDOWN = m_cutsDOWN[sys];

    for (unsigned int i_cut = 0; i_cut != cuts.size(); ++i_cut) {
      string cutname = cuts[i_cut].first;

      map<TString, pair<double, double> > &NM = cuts[i_cut].second ;      // nominal settings
      map<TString, pair<double, double> > &UP = cutsUP[i_cut].second ;    // UP variation for a given sys source
      map<TString, pair<double, double> > &DN = cutsDOWN[i_cut].second ;  // DOWN variation for a given sys source

      if (sys != "topPt_weight") {
        file<< Form("%-25s|| %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f ",
        cutname.data(),
        UP["ttbar"].first     -NM["ttbar"].first     ,             DN["ttbar"].first     -NM["ttbar"].first,
        UP["Drell-Yan"].first -NM["Drell-Yan"].first ,             DN["Drell-Yan"].first -NM["Drell-Yan"].first,
        UP["Single-Top"].first-NM["Single-Top"].first,             DN["Single-Top"].first-NM["Single-Top"].first,
        UP["W+Jets"].first    -NM["W+Jets"].first    ,             DN["W+Jets"].first    -NM["W+Jets"].first,
        UP["background"].first-NM["background"].first,             DN["background"].first-NM["background"].first,
        UP[zprime].first      -NM[zprime].first      ,             DN[zprime].first      -NM[zprime].first,
        UP[gluon].first       -NM[gluon].first       ,             DN[gluon].first       -NM[gluon].first
        ) << endl;

        file<< Form("                         || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f",
        100.*(UP["ttbar"].first     /NM["ttbar"].first      - 1.)    ,      100.*(DN["ttbar"].first     /NM["ttbar"].first      - 1.),
        100.*(UP["Drell-Yan"].first /NM["Drell-Yan"].first  - 1.)    ,      100.*(DN["Drell-Yan"].first /NM["Drell-Yan"].first  - 1.),
        100.*(UP["Single-Top"].first/NM["Single-Top"].first - 1.)    ,      100.*(DN["Single-Top"].first/NM["Single-Top"].first - 1.),
        100.*(UP["W+Jets"].first    /NM["W+Jets"].first     - 1.)    ,      100.*(DN["W+Jets"].first    /NM["W+Jets"].first     - 1.),
        100.*(UP["background"].first/NM["background"].first - 1.)    ,      100.*(DN["background"].first/NM["background"].first - 1.),
        100.*(UP[zprime].first      /NM[zprime].first       - 1.)    ,      100.*(DN[zprime].first      /NM[zprime].first       - 1.),
        100.*(UP[gluon].first       /NM[gluon].first        - 1.)    ,      100.*(DN[gluon].first       /NM[gluon].first        - 1.)
        ) << endl;
      }
      else { // top pT reweighting systematic only affects ttbar
        file<< Form("%-25s|| %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f ",
        cutname.data(),
        UP["ttbar"].first-NM["ttbar"].first,               DN["ttbar"].first-NM["ttbar"].first,
        0., 0., 0., 0., 0., 0.,
        UP["ttbar"].first-NM["ttbar"].first,               DN["ttbar"].first-NM["ttbar"].first, // total bkgd
        0., 0., 0., 0.
        ) << endl;

        file<< Form("                         || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f",
        100.*(UP["ttbar"].first/NM["ttbar"].first - 1.),                      100.*(DN["ttbar"].first/NM["ttbar"].first - 1.),
        0., 0., 0., 0., 0., 0.,
        100.*((UP["ttbar"].first-NM["ttbar"].first)/NM["background"].first),  100.*((DN["ttbar"].first-NM["ttbar"].first)/NM["background"].first),
        0., 0., 0., 0.
        ) << endl;
      }

      if (cutname == "MET Filters" || cutname == ">= 1 jet")
      file << "-----------------------------------------------------------------------------------------------------------------------------------------" << endl;

      // ... Now add up deltas in quadrature for total systematics
      // ... delta function correctly returns + and - variations
      map<TString, pair<double, double> > &m = cuts_TotalSysPlusMinus[i_cut].second;
      for (map<TString, pair<double, double> >::iterator i_set = m.begin(); i_set != m.end(); ++i_set) {

        TString dataset = i_set->first;
        double deltaUP, deltaDN;

        if (dataset == "ttbar" || sys != "topPt_weight") {
          deltaUP = UP[dataset].first - NM[dataset].first ;
          deltaDN = DN[dataset].first - NM[dataset].first ;
        }
        else {
          deltaUP = (dataset == "background") ? (UP["ttbar"].first - NM["ttbar"].first) : 0.;
          deltaDN = (dataset == "background") ? (DN["ttbar"].first - NM["ttbar"].first) : 0.;
        }

        double d1 = delta(deltaUP, deltaDN, "+") ;
        double d2 = delta(deltaUP, deltaDN, "-") ;

        m[dataset].first  += d1*d1 ;
        m[dataset].second += d2*d2 ;

      } //end dataset loop
    } //end cut loop
  } //end systematics loop

  if (systematics.size() > 0) {

    file<<"\n===================================================================================================================================="<< "\n" ;
    file<<"                                    Total systematics        [absolute +- on first line and relative +- in \% on second line] \n"    ;
    file<<"===================================================================================================================================="<< "\n" ;
    file<<Form("                         ||        ttbar        ||      Drell-Yan      ||     Single-Top      ||        W+Jets       ||      background     || %-20s|| %-20s",
    zprime.Data(), gluon.Data() ) << endl;

    for (unsigned int i_cut = 0; i_cut != cuts.size(); ++i_cut) {
      string cutname = cuts[i_cut].first;
      map<TString, pair<double, double> > &NM = cuts[i_cut].second ;
      map<TString, pair<double, double> > &TS = cuts_TotalSysPlusMinus[i_cut].second ;

      file<< Form("%-25s|| %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f ",
             cutname.data(),
             sqrt(TS["ttbar"].first)     ,     sqrt(TS["ttbar"].second),
             sqrt(TS["Drell-Yan"].first) ,     sqrt(TS["Drell-Yan"].second),
             sqrt(TS["Single-Top"].first),     sqrt(TS["Single-Top"].second),
             sqrt(TS["W+Jets"].first)    ,     sqrt(TS["W+Jets"].second ),
             sqrt(TS["background"].first),     sqrt(TS["background"].second),
             sqrt(TS[zprime].first)      ,     sqrt(TS[zprime].second),
             sqrt(TS[gluon].first)       ,     sqrt(TS[gluon].second)
             ) << endl;

      file<< Form("                         || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f || %9.2f %9.2f",
             100.*sqrt(TS["ttbar"].first)     /NM["ttbar"].first     ,    100.*sqrt(TS["ttbar"].second)     /NM["ttbar"].first,
             100.*sqrt(TS["Drell-Yan"].first) /NM["Drell-Yan"].first ,    100.*sqrt(TS["Drell-Yan"].second) /NM["Drell-Yan"].first,
             100.*sqrt(TS["Single-Top"].first)/NM["Single-Top"].first,    100.*sqrt(TS["Single-Top"].second)/NM["Single-Top"].first,
             100.*sqrt(TS["W+Jets"].first)    /NM["W+Jets"].first    ,    100.*sqrt(TS["W+Jets"].second)    /NM["W+Jets"].first,
             100.*sqrt(TS["background"].first)/NM["background"].first,    100.*sqrt(TS["background"].second)/NM["background"].first,
             100.*sqrt(TS[zprime].first)      /NM[zprime].first      ,    100.*sqrt(TS[zprime].second)      /NM[zprime].first,
             100.*sqrt(TS[gluon].first)       /NM[gluon].first       ,    100.*sqrt(TS[gluon].second)       /NM[gluon].first
             ) << endl;

      if (cutname == "MET Filters" || cutname == ">= 1 jet")
        file << "-----------------------------------------------------------------------------------------------------------------------------------------" << endl;
    }
  }

  file << "\n\\begin{center}\n";
  file << "  \\begin{tabular}{ |c||c|c|c|c|c|c|c|c| }\n";
  file << "  \\hline\n";
  file << Form("  Cut & Data & ttbar & Drell-Yan & Single-Top & W+Jets & Background & %s & %s \\\\", zprime.Data(), gluon.Data() ) << endl;
  file << "  \\hline\\hline\n";

  for (vector<pair<string, map<TString, pair<double, double> > > >::iterator i_cut = cuts.begin(); i_cut != cuts.end(); ++i_cut) {
    string cutname = i_cut->first;
    file << Form("  %s & %.0f $\\pm$ %.1f & %.1f $\\pm$ %.1f & %.1f $\\pm$ %.1f & %.1f $\\pm$ %.1f & %.1f $\\pm$ %.1f & %.1f $\\pm$ %.1f & %.1f $\\pm$ %.1f & %.1f $\\pm$ %.1f \\\\\n  \\hline",
                cutname.data(), i_cut->second["data"].first, i_cut->second["data"].second, i_cut->second["ttbar"].first, i_cut->second["ttbar"].second,
                i_cut->second["Drell-Yan"].first, i_cut->second["Drell-Yan"].second, i_cut->second["Single-Top"].first, i_cut->second["Single-Top"].second,
                i_cut->second["W+Jets"].first, i_cut->second["W+Jets"].second, i_cut->second["background"].first, i_cut->second["background"].second,
                i_cut->second[zprime].first, i_cut->second[zprime].second, i_cut->second[gluon].first, i_cut->second[gluon].second) << endl;

    if (cutname == "MET Filters" || cutname == ">= 1 jet")
      file << "  \\hline\n";
  }

  file << "  \\end{tabular}\n";
  file << "\\end{center}" << endl;

  file.close();

  ofstream e_outfile( channel + "_efficiencies.txt" );

  ifstream data_infile(dataFile);
  ifstream mc_infile(mcFile);

  e_outfile << data_infile.rdbuf();
  data_infile.close();
  e_outfile << mc_infile.rdbuf();
  mc_infile.close();

  e_outfile << "\n====================================================================================================================="<< "\n" ;
  e_outfile << "                                              Cut Flow Table: Summary\n" ;
  e_outfile << "====================================================================================================================="<< "\n" ;
  e_outfile << Form("                          |||    ttbar   |||  Drell-Yan ||| Single-Top |||   W+Jets   ||| %-20s ||| %-20s", zprime.Data(), gluon.Data() ) << endl;

  for (vector<pair<string, map<TString, pair<double, double> > > >::iterator i_cut = cuts.begin(); i_cut != cuts.end(); ++i_cut) {
    string cutname = i_cut->first;
    e_outfile << Form("%-25s |||  %1.6f  |||  %1.6f  |||  %1.6f  |||  %1.6f  |||  %1.6f  |||  %1.6f",
                cutname.data(), i_cut->second["ttbar"].first/m_total["ttbar"].first, i_cut->second["Drell-Yan"].first/m_total["Drell-Yan"].first,
                i_cut->second["Single-Top"].first/m_total["Single-Top"].first, i_cut->second["W+Jets"].first/m_total["W+Jets"].first,
                i_cut->second[zprime].first/m_total[zprime].first, i_cut->second[gluon].first/m_total[gluon].first ) << endl;

    if (cutname == "MET Filters" || cutname == ">= 1 jet")
      e_outfile << "---------------------------------------------------------------------------------------------------------------------" << endl;
  }

  e_outfile << "\n\\begin{center}\n";
  e_outfile << "  \\begin{tabular}{ |c|c| }\n";
  e_outfile << "  \\hline\n";
  e_outfile << "  Process & Efficiency \\\\\n";
  e_outfile << "  \\hline\\hline\n";

  map<TString, pair<double, double> >& m_last = cuts.back().second;

  e_outfile << Form("  ttbar & %1.6f \\\\\n  \\hline", m_last["ttbar"].first/m_total["ttbar"].first) << endl;
  e_outfile << Form("  Drell-Yan & %1.6f \\\\\n  \\hline", m_last["Drell-Yan"].first/m_total["Drell-Yan"].first) << endl;
  e_outfile << Form("  Single-Top & %1.6f \\\\\n  \\hline", m_last["Single-Top"].first/m_total["Single-Top"].first) << endl;
  e_outfile << Form("  W+Jets & %1.6f \\\\\n  \\hline", m_last["W+Jets"].first/m_total["W+Jets"].first) << endl;
  e_outfile << Form("  %-20s & %1.6f \\\\\n  \\hline", zprime.Data(), m_last[zprime].first/m_total[zprime].first) << endl;
  e_outfile << Form("  %-20s & %1.6f \\\\\n  \\hline", gluon.Data(), m_last[gluon].first/m_total[gluon].first) << endl;

  e_outfile << "  \\end{tabular}\n";
  e_outfile << "\\end{center}" << endl;

  e_outfile.close();
}

//vector(cut_name, map(dataset, (N, error) ) )
void readFile(const string& fileName, vector<pair<string, map<TString, pair<double, double> > > >& cuts, string& channel) {

  ifstream file(fileName);
  string line;
  TString dataset;
  double weight=-1;
  vector<pair<string, map<TString, pair<double, double> > > >::iterator i_cut;

  cout << " =====> Reading file:  " << fileName << endl ;

  while (getline(file, line)){

    if (line.length() > 0){
      while (line.at(0) == ' ') line.erase(0, 1);
    }

    int delim_pos = line.find(' ');
    if (delim_pos == -1) continue;
    string str = line.substr(0, delim_pos);

    if (str == "Weight") {
      delim_pos = line.length()-1;
      while (line.at(delim_pos) != ' ') delim_pos--;

      line.erase(0, delim_pos+1);
      weight = stod(line);

      //begin cut interator
      i_cut = cuts.begin();
    }
    else if (str == "Channel:") {
      delim_pos = line.length()-1;
      while (line.at(delim_pos) != ' ') delim_pos--;

      line.erase(0, delim_pos+1);
      channel = line;
    }
    else if (str == "Cut") {
      delim_pos = line.length()-1;
      while (line.at(delim_pos) != ' ') delim_pos--;

      line.erase(0, delim_pos+1);
      dataset = line.data();
      cout << "Filling " << dataset << " using weight " << weight << endl;
    }
    else {
      delim_pos = line.find("|||");
      if (delim_pos == -1) continue;

      str = line.substr(0, delim_pos);
      if (str == "") continue;

      while (str.at(str.length()-1) == ' ') str.erase(str.length()-1, str.length());
      string cut_name = str;

      line.erase(0, delim_pos + 3);
      while (line.at(0) == ' ') line.erase(0, 1);

      str = line.substr(0, line.find("|||"));
      while (str.at(str.length()-1) == ' ') str.erase(str.length()-1, str.length());

      double N = stod(str);

      //data
      if ( dataset.Contains("Muon", TString::kIgnoreCase) || dataset.Contains("Electron", TString::kIgnoreCase) ) {
        cuts.push_back( make_pair(cut_name, map<TString, pair<double, double> >()) );
        map<TString, pair<double, double> >& m = cuts.back().second;
        m["data"] = make_pair( N, sqrt(N) );
      }

      //signal or background
      else {
        TString key;
        map<TString, pair<double, double> >& m = i_cut->second;

        if ( dataset.Contains("ttbar", TString::kIgnoreCase) )     key = "ttbar";
        else if ( dataset.Contains("dy", TString::kIgnoreCase) )   key = "Drell-Yan";
        else if ( dataset.Contains("wjet", TString::kIgnoreCase) ) key = "W+Jets";
        else if ( dataset.Contains("st", TString::kIgnoreCase)
               || dataset.Contains("sat", TString::kIgnoreCase) )  key = "Single-Top";
        else                                                       key = dataset;

        if ( m.find(key) == m.end() ) m[key] = make_pair( N, weight * sqrt(N/weight) );
        else {
          m[key].first += N;
          m[key].second = sqrt( m[key].second*m[key].second + weight*N );
        }

        //total background
        if ( !dataset.Contains("zprime", TString::kIgnoreCase) && !dataset.Contains("gluon", TString::kIgnoreCase) ) {
          key = "background";
          if ( m.find(key) == m.end() ) m[key] = make_pair( N, weight * sqrt(N/weight) );
          else {
            m[key].first += N;
            m[key].second = sqrt( m[key].second*m[key].second + weight*N );
          }
        }

        ++i_cut;
      }
    }
  }
  file.close();
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

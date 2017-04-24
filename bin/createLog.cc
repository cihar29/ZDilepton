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

int main(int argc, char* argv[]) {

  if (argc == 1 || argc == 2) { cout << "Please provide two log files." << endl; return -1; }

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

  ofstream file( channel + "_finalLog.txt" );

  file<<"====================================================================================================================="<< "\n" ;
  file<<"                                              Cut Flow Table: Summary\n" ;
  file<<"====================================================================================================================="<< "\n" ;
  file<< Form("                               |||               Data                |||             Background            |||        %-20s       |||        %-20s", zprime.Data(), gluon.Data() ) << endl;

  for (vector<pair<string, map<TString, pair<double, double> > > >::iterator i_cut = cuts.begin(); i_cut != cuts.end(); ++i_cut)
    file<< Form("%-30s |||      %12.0f (%1.6f)      |||      %12.1f (%1.6f)      |||      %12.1f (%1.6f)      |||      %12.1f (%1.6f)",
                i_cut->first.data(), i_cut->second["data"].first, i_cut->second["data"].first/m_total["data"].first,
                i_cut->second["background"].first, i_cut->second["background"].first/m_total["background"].first,
                i_cut->second[zprime].first, i_cut->second[zprime].first/m_total[zprime].first,
                i_cut->second[gluon].first, i_cut->second[gluon].first/m_total[gluon].first ) << endl;

  file<<"\n====================================================================================================================="<< "\n" ;
  file<<"                                              Cut Flow Table: Background\n" ;
  file<<"====================================================================================================================="<< "\n" ;
  file<<"                               |||              ttbar                |||             Drell-Yan             |||           Single-Top              |||           W+Jets" << endl;

  for (vector<pair<string, map<TString, pair<double, double> > > >::iterator i_cut = cuts.begin(); i_cut != cuts.end(); ++i_cut)
    file<< Form("%-30s |||      %12.1f (%1.6f)      |||      %12.1f (%1.6f)      |||      %12.1f (%1.6f)      |||      %12.1f (%1.6f)",
                i_cut->first.data(), i_cut->second["ttbar"].first, i_cut->second["ttbar"].first/m_total["ttbar"].first,
                i_cut->second["Drell-Yan"].first, i_cut->second["Drell-Yan"].first/m_total["Drell-Yan"].first,
                i_cut->second["Single-Top"].first, i_cut->second["Single-Top"].first/m_total["Single-Top"].first,
                i_cut->second["W+Jets"].first, i_cut->second["W+Jets"].first/m_total["W+Jets"].first ) << endl;


  file << "\n\\begin{center}\n";
  file << "  \\begin{tabular}{ |c||c|c|c|c|c|c|c|c| }\n";
  file << "  \\hline\n";
  file << Form("  Cut & Data & ttbar & Drell-Yan & Single-Top & W+Jets & Background & %s & %s \\\\", zprime.Data(), gluon.Data() ) << endl;
  file << "  \\hline\\hline\n";

  for (vector<pair<string, map<TString, pair<double, double> > > >::iterator i_cut = cuts.begin(); i_cut != cuts.end(); ++i_cut)
    file << Form("  %s & %.0f $\\pm$ %.1f & %.1f $\\pm$ %.1f & %.1f $\\pm$ %.1f & %.1f $\\pm$ %.1f & %.1f $\\pm$ %.1f & %.1f $\\pm$ %.1f & %.1f $\\pm$ %.1f & %.1f $\\pm$ %.1f \\\\\n  \\hline",
                i_cut->first.data(), i_cut->second["data"].first, i_cut->second["data"].second, i_cut->second["ttbar"].first, i_cut->second["ttbar"].second,
                i_cut->second["Drell-Yan"].first, i_cut->second["Drell-Yan"].second, i_cut->second["Single-Top"].first, i_cut->second["Single-Top"].second,
                i_cut->second["W+Jets"].first, i_cut->second["W+Jets"].second, i_cut->second["background"].first, i_cut->second["background"].second,
                i_cut->second[zprime].first, i_cut->second[zprime].second, i_cut->second[gluon].first, i_cut->second[gluon].second) << endl;

  file << "  \\end{tabular}\n";
  file << "\\end{center}" << endl;

  file.close();

  ifstream mc_infile(mcFile);
  ofstream mc_outfile("final_" + mcFile);

  mc_outfile << mc_infile.rdbuf();
  mc_infile.close();

  mc_outfile << "\n====================================================================================================================="<< "\n" ;
  mc_outfile << "                                              Cut Flow Table: Summary\n" ;
  mc_outfile << "====================================================================================================================="<< "\n" ;
  mc_outfile << Form("                               |||    ttbar   |||  Drell-Yan ||| Single-Top |||   W+Jets   ||| %-20s ||| %-20s", zprime.Data(), gluon.Data() ) << endl;

  for (vector<pair<string, map<TString, pair<double, double> > > >::iterator i_cut = cuts.begin(); i_cut != cuts.end(); ++i_cut)
    mc_outfile << Form("%-30s |||  %1.6f  |||  %1.6f  |||  %1.6f  |||  %1.6f  |||  %1.6f  |||  %1.6f",
                i_cut->first.data(), i_cut->second["ttbar"].first/m_total["ttbar"].first, i_cut->second["Drell-Yan"].first/m_total["Drell-Yan"].first,
                i_cut->second["Single-Top"].first/m_total["Single-Top"].first, i_cut->second["W+Jets"].first/m_total["W+Jets"].first,
                i_cut->second[zprime].first/m_total[zprime].first, i_cut->second[gluon].first/m_total[gluon].first ) << endl;

  mc_outfile << "\n\\begin{center}\n";
  mc_outfile << "  \\begin{tabular}{ |c|c| }\n";
  mc_outfile << "  \\hline\n";
  mc_outfile << "  Process & Efficiency \\\\\n";
  mc_outfile << "  \\hline\\hline\n";

  map<TString, pair<double, double> >& m_last = cuts.back().second;

  mc_outfile << Form("  ttbar & %1.6f \\\\\n  \\hline", m_last["ttbar"].first/m_total["ttbar"].first) << endl;
  mc_outfile << Form("  Drell-Yan & %1.6f \\\\\n  \\hline", m_last["Drell-Yan"].first/m_total["Drell-Yan"].first) << endl;
  mc_outfile << Form("  Single-Top & %1.6f \\\\\n  \\hline", m_last["Single-Top"].first/m_total["Single-Top"].first) << endl;
  mc_outfile << Form("  W+Jets & %1.6f \\\\\n  \\hline", m_last["W+Jets"].first/m_total["W+Jets"].first) << endl;
  mc_outfile << Form("  %-20s & %1.6f \\\\\n  \\hline", zprime.Data(), m_last[zprime].first/m_total[zprime].first) << endl;
  mc_outfile << Form("  %-20s & %1.6f \\\\\n  \\hline", gluon.Data(), m_last[gluon].first/m_total[gluon].first) << endl;

  mc_outfile << "  \\end{tabular}\n";
  mc_outfile << "\\end{center}" << endl;

  mc_outfile.close();
}

//vector(cut_name, map(dataset, (N, error) ) )
void readFile(const string& fileName, vector<pair<string, map<TString, pair<double, double> > > >& cuts, string& channel) {

  ifstream file(fileName);
  string line;
  TString dataset;
  double weight=-1;
  vector<pair<string, map<TString, pair<double, double> > > >::iterator i_cut;

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
      dataset = line.substr(0, line.find(".root")).data();
      cout << "Filling " << dataset << " using weight " << weight << endl;
    }
    else if (str == "Initial" || str == "Passed") {
      delim_pos = line.find("|||");
      str = line.substr(0, delim_pos);
      while (str.at(str.length()-1) == ' ') str.erase(str.length()-1, str.length());

      string cut_name = str;

      line.erase(0, delim_pos + 3);
      while (line.at(0) == ' ') line.erase(0, 1);

      str = line.substr(0, line.find("|||"));
      while (str.at(str.length()-1) == ' ') str.erase(str.length()-1, str.length());

      double unwght_N = stod(str);

      //data
      if ( dataset.Contains("Muon", TString::kIgnoreCase) || dataset.Contains("Electron", TString::kIgnoreCase) ) {
        cuts.push_back( make_pair(cut_name, map<TString, pair<double, double> >()) );
        map<TString, pair<double, double> >& m = cuts.back().second;
        m["data"] = make_pair( unwght_N, sqrt(unwght_N) );
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

        if ( m.find(key) == m.end() ) m[key] = make_pair( weight * unwght_N, weight * sqrt(unwght_N) );
        else {
          m[key].first += weight * unwght_N;
          m[key].second = sqrt( m[key].second*m[key].second + weight*weight*unwght_N );
        }

        //total background
        if ( !dataset.Contains("zprime", TString::kIgnoreCase) && !dataset.Contains("gluon", TString::kIgnoreCase) ) {
          key = "background";
          if ( m.find(key) == m.end() ) m[key] = make_pair( weight * unwght_N, weight * sqrt(unwght_N) );
          else {
            m[key].first += weight * unwght_N;
            m[key].second = sqrt( m[key].second*m[key].second + weight*weight*unwght_N );
          }
        }

        ++i_cut;
      }
    }
  }
  file.close();
}

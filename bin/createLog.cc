#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <iomanip>
#include <string>

#include "TString.h"

using namespace std;

void readFile(const string& fileName, vector<tuple<string, double, double, double, double> >& cuts);
enum Tuple{
  cut = 0, data, background, zprime, gluon, numSets
};

int main(int argc, char* argv[]) {

  if (argc == 1 || argc == 2) { cout << "Please provide two log files." << endl; return -1; }

  string dataFile = argv[1];
  string mcFile = argv[2];

  //tuple(cut_name, data, background, zprime, gluon)
  //use vector to preserve order of cuts
  vector<tuple<string, double, double, double, double> > cuts;

  //read data first to initialize mc
  readFile(dataFile, cuts);
  readFile(mcFile, cuts);

  ofstream file( "finalLog.txt" );

  file<<"====================================================================================================================="<< "\n" ;
  file<<"                                              Cut Flow Table: Summary\n" ;
  file<<"====================================================================================================================="<< "\n" ;
  file<<"                               |||               Data                |||             Background            |||             Zprime                |||             Gluon" << endl;

  for (vector<tuple<string, double, double, double, double> >::iterator i_cut = cuts.begin(); i_cut != cuts.end(); ++i_cut)
    file<< Form("%-30s |||      %12.1f (%1.6f)      |||      %12.1f (%1.6f)      |||      %12.1f (%1.6f)      |||      %12.1f (%1.6f)",
                get<cut>(*i_cut).data(), get<data>(*i_cut), get<data>(*i_cut)/get<data>(cuts[0]), get<background>(*i_cut), get<background>(*i_cut)/get<background>(cuts[0]),
                get<zprime>(*i_cut), get<zprime>(*i_cut)/get<zprime>(cuts[0]), get<gluon>(*i_cut), get<gluon>(*i_cut)/get<gluon>(cuts[0]) ) << endl;
  file.close();
}

void readFile(const string& fileName, vector<tuple<string, double, double, double, double> >& cuts) {

  ifstream file(fileName);
  string line;
  TString dataset;
  double weight=-1;
  vector<tuple<string, double, double, double, double> >::iterator i_cut;

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

      //background
      if ( !dataset.Contains("Muon", TString::kIgnoreCase) && !dataset.Contains("gluon", TString::kIgnoreCase) && !dataset.Contains("zprime", TString::kIgnoreCase) )
        get<background>(*i_cut) += weight*stod(str);

      //data
      else if (dataset.Contains("Muon", TString::kIgnoreCase))
        cuts.push_back( make_tuple(cut_name, weight*stod(str), 0, 0, 0) );

      //signals
      else if (dataset.Contains("zprime", TString::kIgnoreCase))
        get<zprime>(*i_cut) = weight*stod(str);
      else if (dataset.Contains("gluon", TString::kIgnoreCase))
        get<gluon>(*i_cut) = weight*stod(str);

      ++i_cut;
    }
  }
  file.close();
}

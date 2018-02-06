//Bahareh Roozbahani Dec 2017
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include <vector>
#include <map>

void setStyle();
void plot_trigger(string testLeg = "muLeg",TString type=""){

  setStyle();
  TCanvas* c = new TCanvas("c", "c", 600, 600);

  TLegend* leg = new TLegend(.4,.65,.8,.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);

  TLegend* leg_eff = new TLegend(.4,.85,.8,.9);
  leg_eff->SetBorderSize(0);
  leg_eff->SetFillColor(0);
  leg_eff->SetFillStyle(0);
  leg_eff->SetTextSize(0.035);
  leg_eff->SetTextFont(42);

  TLatex text;
  text.SetNDC();
  text.SetTextSize(0.045);
  text.SetTextFont(42);

  vector<double> rebin_loose;


  string probtrig, probelep, probetrig, channel;
  if (testLeg=="muLeg") { channel = "em"; probelep = "subMu"; probetrig = "HLT_Mu50_";
    rebin_loose = {0,10,15,20,25,30,40,45,48,50,52,55,60,80,120,200,300,400,800};
 }
  else if (testLeg=="eleLeg"){ channel = "ee"; probelep = "subEle"; probetrig = "HLT_DoubleEle33_CaloIdL_"; 
    rebin_loose = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,36,40,50,60,70,80,90,100,150,200,250,300,500};
 }

  TFile *file_Data, *file_MC, *file;

  file_Data = TFile::Open( Form("root_files/%s/Ele_%s.root", channel.data(), channel.data()) );
  file_MC   = TFile::Open( Form("root_files/%s/TTBar_%s.root", channel.data(), channel.data()) ); 

  TString sample;
  cout<<"processing data or ttbar? " << endl;
  cin >> sample;
  if (sample=="data")  file = file_Data;
  else 		       file = file_MC;
  TDirectory* dir = (TDirectory*) file;
  vector<TString> name_denom_bs , name_denom_as;
  vector<TString> name_nom_bs , name_nom_as;
  
  
  if (type!="") type = type + "_"; //for Onnb and Onbt
  TIter HistList_tag(file->GetDirectory(type+"HLT_Ele115_CaloIdVT_")->GetListOfKeys());
  TIter HistList_probe(file->GetDirectory(type+Form("%s",probetrig.data()))->GetListOfKeys());

  while (TKey* histKey = (TKey*)HistList_tag() ){
    TString keyname = histKey->GetName();
    cout << keyname << endl;
    if (keyname.Contains("_BS",TString::kIgnoreCase)) name_denom_bs.push_back(keyname);
    else name_denom_as.push_back(keyname); 
  }
  while (TKey* histKey = (TKey*)HistList_probe() ){
    TString keyname = histKey->GetName();
    cout << keyname << endl;
    if (keyname.Contains("_BS",TString::kIgnoreCase)) name_nom_bs.push_back(keyname);
    else name_nom_as.push_back(keyname); 
  }
  map<TString, TH1F*> h_pT, h_eta , h_rmin, MC_h_rmin;
  map<TString, TH2F*> h_eta_pT, MC_h_eta_pT;

  for(int i=0; i< int(name_denom_bs.size()); i++){
    if (name_denom_bs[i].Contains(Form("%s_eta_pT_",probelep.data()),TString::kIgnoreCase)){
      h_eta_pT["denom_BS"] = (TH2F*) file->FindKeyAny(name_denom_bs[i])->ReadObj();
      MC_h_eta_pT["denom_BS"] = (TH2F*) file_MC->FindKeyAny(name_denom_bs[i])->ReadObj();
    }
    else if (name_denom_bs[i].Contains(Form("%s_eta",probelep.data()),TString::kIgnoreCase) )
      h_eta["denom_BS"] = (TH1F*) file->FindKeyAny(name_denom_bs[i])->ReadObj();
    else if (name_denom_bs[i].Contains(Form("%s_pT",probelep.data()),TString::kIgnoreCase) )  
      h_pT["denom_BS"] = (TH1F*) file->FindKeyAny(name_denom_bs[i])->ReadObj();
    else if (name_denom_bs[i].Contains(Form("%s_rmin",probelep.data()),TString::kIgnoreCase) ){
      h_rmin["denom_BS"]    = (TH1F*) file->FindKeyAny(name_denom_bs[i])->ReadObj();
      MC_h_rmin["denom_BS"] = (TH1F*) file_MC->FindKeyAny(name_denom_bs[i])->ReadObj();
    }
  }
  for(int i=0; i< int(name_denom_as.size()); i++){
    if (name_denom_as[i].Contains(Form("%s_eta_pT_",probelep.data()),TString::kIgnoreCase)){
      h_eta_pT["denom_AS"] = (TH2F*) file->FindKeyAny(name_denom_as[i])->ReadObj();
      MC_h_eta_pT["denom_AS"] = (TH2F*) file_MC->FindKeyAny(name_denom_as[i])->ReadObj();
    }
    else if (name_denom_as[i].Contains(Form("%s_eta",probelep.data()),TString::kIgnoreCase) ) 
      h_eta["denom_AS"] = (TH1F*) file->FindKeyAny(name_denom_as[i])->ReadObj();
    else if (name_denom_as[i].Contains(Form("%s_pT",probelep.data()),TString::kIgnoreCase) ) 
      h_pT["denom_AS"] = (TH1F*) file->FindKeyAny(name_denom_as[i])->ReadObj();
    else if (name_denom_as[i].Contains(Form("%s_rmin",probelep.data()),TString::kIgnoreCase) ){
      h_rmin["denom_AS"]   = (TH1F*) file->FindKeyAny(name_denom_as[i])->ReadObj();
      MC_h_rmin["denom_AS"]= (TH1F*) file_MC->FindKeyAny(name_denom_as[i])->ReadObj();
    }
  }

  for(int i=0; i< int(name_nom_bs.size()); i++){
    if (name_nom_bs[i].Contains(Form("%s_eta_pT_",probelep.data()),TString::kIgnoreCase)){
      h_eta_pT["nom_BS"] = (TH2F*) file->FindKeyAny(name_nom_bs[i])->ReadObj();
      MC_h_eta_pT["nom_BS"] = (TH2F*) file_MC->FindKeyAny(name_nom_bs[i])->ReadObj();
    }
    else if (name_nom_bs[i].Contains(Form("%s_eta",probelep.data()),TString::kIgnoreCase) ) 
      h_eta["nom_BS"] = (TH1F*) file->FindKeyAny(name_nom_bs[i])->ReadObj();
    else if (name_nom_bs[i].Contains(Form("%s_pT",probelep.data()),TString::kIgnoreCase) )  
      h_pT["nom_BS"] = (TH1F*) file->FindKeyAny(name_nom_bs[i])->ReadObj();
    else if (name_nom_bs[i].Contains(Form("%s_rmin",probelep.data()),TString::kIgnoreCase) ){
      h_rmin["nom_BS"]   = (TH1F*) file->FindKeyAny(name_nom_bs[i])->ReadObj();
      MC_h_rmin["nom_BS"]= (TH1F*) file_MC->FindKeyAny(name_nom_bs[i])->ReadObj();
    }
  }
  for(int i=0; i< int(name_nom_as.size()); i++){
    if (name_nom_as[i].Contains(Form("%s_eta_pT_",probelep.data()),TString::kIgnoreCase)){
      h_eta_pT["nom_AS"] = (TH2F*) file->FindKeyAny(name_nom_as[i])->ReadObj();
      MC_h_eta_pT["nom_AS"] = (TH2F*) file_MC->FindKeyAny(name_nom_as[i])->ReadObj();
    }
    else if (name_nom_as[i].Contains(Form("%s_eta",probelep.data()),TString::kIgnoreCase) ) 
      h_eta["nom_AS"] = (TH1F*) file->FindKeyAny(name_nom_as[i])->ReadObj();
    else if (name_nom_as[i].Contains(Form("%s_pT",probelep.data()),TString::kIgnoreCase) )  
      h_pT["nom_AS"] = (TH1F*) file->FindKeyAny(name_nom_as[i])->ReadObj();
    else if (name_nom_as[i].Contains(Form("%s_rmin",probelep.data()),TString::kIgnoreCase) ){
      h_rmin["nom_AS"] = (TH1F*) file->FindKeyAny(name_nom_as[i])->ReadObj();
      MC_h_rmin["nom_AS"] = (TH1F*) file_MC->FindKeyAny(name_nom_as[i])->ReadObj();
    }
  }

  TH1F* pTh = 0; TH1F* etah = 0; TH1F* rminh = 0;
  int pTBins = h_pT["denom_BS"]->GetNbinsX(); 
  int etaBins = h_eta["denom_BS"]->GetNbinsX();
  int rminBins = h_rmin["denom_BS"]->GetNbinsX();
  pTh = new TH1F("pThist", "pThist", pTBins , h_pT["denom_BS"]->GetBinLowEdge(1), h_pT["denom_BS"]->GetBinLowEdge(pTBins+1));
  etah = new TH1F("etahist", "etahist", etaBins , h_eta["denom_BS"]->GetBinLowEdge(1), h_eta["denom_BS"]->GetBinLowEdge(etaBins+1));
  rminh = new TH1F("rminhist", "rminhist", rminBins , h_rmin["denom_BS"]->GetBinLowEdge(1), h_rmin["denom_BS"]->GetBinLowEdge(rminBins+1));
  pTh->GetYaxis()->SetRangeUser(0.1, int(h_pT["denom_BS"]->GetMaximum()*1.5) );
  etah->GetYaxis()->SetRangeUser(0, int(h_eta["denom_BS"]->GetMaximum()*1.5) );
  rminh->GetYaxis()->SetRangeUser(0, int(h_rmin["denom_BS"]->GetMaximum()*1.5) );

  pTh->GetXaxis()->SetLabelSize(0.03);
  pTh->GetYaxis()->SetLabelSize(0.03);

  etah->GetXaxis()->SetRangeUser(-5, 5);
  etah->GetXaxis()->SetLabelSize(0.03);
  etah->GetYaxis()->SetLabelSize(0.03);

  rminh->GetXaxis()->SetRangeUser(0, 5);
  rminh->GetXaxis()->SetLabelSize(0.03);
  rminh->GetYaxis()->SetLabelSize(0.03);

  if (testLeg=="eleLeg") { 
    pTh->GetXaxis()->SetTitle("p_{T}^{subleading e} (GeV)");
    pTh->GetXaxis()->SetRangeUser(0, 200);
    etah->GetXaxis()->SetTitle("supercluster |#eta^{subleading e}|"); 
    rminh->GetXaxis()->SetTitle("#DeltaR(subleading e , closest jet)");
  }
  else if (testLeg=="muLeg") { 
    pTh->GetXaxis()->SetTitle("p_{T}^{#mu} (GeV)");
    pTh->GetXaxis()->SetRangeUser(0, 800);
    etah->GetXaxis()->SetTitle("|#eta^{#mu}|"); 
    rminh->GetXaxis()->SetTitle("#DeltaR(#mu , closest jet)");
  }

  for (map<TString, TH1F*>::iterator it = h_pT.begin(); it != h_pT.end(); it++) {
    TString name = it->first ;
    h_pT[name]->SetMarkerStyle(22);
    h_eta[name]->SetMarkerStyle(22);
    h_rmin[name]->SetMarkerStyle(22);
    if (name.Contains("denom",TString::kIgnoreCase)){
      h_pT[name]->SetMarkerColor(kRed);   h_pT[name]->SetLineColor(kRed);
      h_eta[name]->SetMarkerColor(kRed);  h_eta[name]->SetLineColor(kRed);
      h_rmin[name]->SetMarkerColor(kRed); h_rmin[name]->SetLineColor(kRed);
    }
    else{
      h_pT[name]->SetMarkerColor(kBlue);  h_pT[name]->SetLineColor(kBlue);
      h_eta[name]->SetMarkerColor(kBlue); h_eta[name]->SetLineColor(kBlue);
      h_rmin[name]->SetMarkerColor(kBlue);h_rmin[name]->SetLineColor(kBlue);
    }
  }
  c->cd();

  TH1F* clone_denom_BS_pT = (TH1F*) h_pT["denom_BS"]->Clone("clone_denom_BS_pT");
  TH1F* clone_nom_BS_pT = (TH1F*) h_pT["nom_BS"]->Clone("clone_nom_BS_pT");

  clone_denom_BS_pT->Rebin(10);
  clone_nom_BS_pT->Rebin(10); 

  pTh->GetYaxis()->SetRangeUser(0.1, int(clone_denom_BS_pT->GetMaximum()*1.5) );
  pTh->Draw();
 
  clone_denom_BS_pT->Draw("samepE");
  clone_nom_BS_pT->Draw("samepE");
  leg->SetTextSize(0.036);
  leg->AddEntry(clone_denom_BS_pT, "Before Trigger" , "pE"); 
  leg->AddEntry(clone_nom_BS_pT, "After Trigger" , "pE");
  text.SetTextSize(0.03);
  if (testLeg=="eleLeg"){
    text.DrawLatex(0.6, 0.96, "supercluster |#eta^{subleading e}| < 2.5");
    text.DrawLatex(0.22, 0.86,  "HLT_DoubleEle33");
  }
  else if (testLeg=="muLeg"){
    text.DrawLatex(0.6, 0.96, "|#eta^{#mu}| < 2.4");
    text.DrawLatex(0.22, 0.86,  "HLT_Mu50");
  }
  leg->Draw();
  text.SetTextSize(0.03);
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e}>0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e}<0.4");
    }
    else if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu}>0.4");
      else 	         text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu}<0.4");
    }
  }
  if (sample=="ttbar")  text.DrawLatex(0.2, 0.96, "MC");
  else       text.DrawLatex(0.2, 0.96, "Data");
  c->Print( Form("%s_pT_BS.pdf",probelep.data()) );
  c->Clear();
  leg->Clear();

  c->cd();

  TH1F* clone_denom_AS_pT = (TH1F*) h_pT["denom_AS"]->Clone("clone_denom_AS_pT");
  TH1F* clone_nom_AS_pT = (TH1F*) h_pT["nom_AS"]->Clone("clone_nom_AS_pT");

  clone_denom_AS_pT->Rebin(10);
  clone_nom_AS_pT->Rebin(10); 


  pTh->GetYaxis()->SetRangeUser(0.1, int(clone_denom_AS_pT->GetMaximum()*1.5) );
  pTh->Draw();
  clone_denom_AS_pT->Draw("samepE");
  clone_nom_AS_pT->Draw("samepE");
  leg->SetTextSize(0.036);
  leg->AddEntry(h_pT["denom_AS"], "Before Trigger" , "pE"); 
  leg->AddEntry(h_pT["nom_AS"], "After Trigger" , "pE");
  text.SetTextSize(0.03);
  if (testLeg=="eleLeg"){
    text.DrawLatex(0.6, 0.96, "supercluster |#eta^{subleading e}| < 2.5");
    text.DrawLatex(0.22, 0.86,  "HLT_DoubleEle33");
  }
  else if (testLeg=="muLeg"){
    text.DrawLatex(0.6, 0.96, "supercluster |#eta^{#mu}| < 2.4");
    text.DrawLatex(0.22, 0.86,  "HLT_Mu50");
  }
  //text.DrawLatex(0.5, 0.46, "applied selections");
  leg->Draw();
  text.SetTextSize(0.03);
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	         text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    else if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	         text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
  }
  if (sample=="ttbar")  text.DrawLatex(0.2, 0.96, "MC");
  else       text.DrawLatex(0.2, 0.96, "Data");
  c->Print( Form("%s_pT_AS.pdf",probelep.data()) );
  c->Clear();
  leg->Clear();

  c->cd();
  etah->Draw();
  h_eta["denom_BS"]->Draw("samepE");
  h_eta["nom_BS"]->Draw("samepE");
  leg->SetTextSize(0.036);
  leg->AddEntry(h_eta["denom_BS"], "Before Trigger" , "pE"); 
  leg->AddEntry(h_eta["nom_BS"], "After Trigger" , "pE");
  text.SetTextSize(0.03);
  if (testLeg=="eleLeg"){
    text.DrawLatex(0.63, 0.96, "p_{T}^{subleading e} > 36 GeV");
    text.DrawLatex(0.22, 0.86,  "HLT_DoubleEle33");
  }
  else if (testLeg=="muLeg"){
    text.DrawLatex(0.63, 0.96, "p_{T}^{#mu} > 53 GeV");
    text.DrawLatex(0.22, 0.86,  "HLT_Mu50");
  }
  leg->Draw();
  text.SetTextSize(0.03);
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} < 0.4");
    }
  }
  if (sample=="ttbar")  text.DrawLatex(0.2, 0.96, "MC");
  else       text.DrawLatex(0.2, 0.96, "Data");
  c->Print( Form("%s_eta_BS.pdf",probelep.data()) );
  c->Clear();
  leg->Clear();

  c->cd();
  etah->GetYaxis()->SetRangeUser(0.1, int(h_eta["denom_AS"]->GetMaximum()*1.5) );
  etah->Draw();
  h_eta["denom_AS"]->Draw("samepE");
  h_eta["nom_AS"]->Draw("samepE");
  leg->SetTextSize(0.036);
  leg->AddEntry(h_eta["denom_AS"], "Before Trigger" , "pE"); 
  leg->AddEntry(h_eta["nom_AS"], "After Trigger" , "pE");
  text.SetTextSize(0.03);
  if (testLeg=="eleLeg"){
    text.DrawLatex(0.63, 0.96, "p_{T}^{subleading e} > 36 GeV");
    text.DrawLatex(0.22, 0.86,  "HLT_DoubleEle33");
  }
  else if (testLeg=="muLeg"){
    text.DrawLatex(0.63, 0.96, "p_{T}^{#mu} > 53 GeV");
    text.DrawLatex(0.22, 0.86,  "HLT_Mu50");
  }
  //text.DrawLatex(0.5, 0.46, "applied selections");
  leg->Draw();
  text.SetTextSize(0.03);
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} < 0.4");
    }
  }
  if (sample=="ttbar")  text.DrawLatex(0.2, 0.96, "MC");
  else       text.DrawLatex(0.2, 0.96, "Data");
  c->Print( Form("%s_eta_AS.pdf",probelep.data()) );
  c->Clear();
  leg->Clear();

  c->cd();
  rminh->Draw();
  h_rmin["denom_BS"]->Draw("samepE");
  h_rmin["nom_BS"]->Draw("samepE");
  leg->SetTextSize(0.036);
  leg->AddEntry(h_rmin["denom_BS"], "Before Trigger" , "pE"); 
  leg->AddEntry(h_rmin["nom_BS"], "After Trigger" , "pE");
  text.SetTextSize(0.03);
  if (testLeg=="eleLeg") {
    text.DrawLatex(0.30, 0.96, "p_{T}^{subleading e} > 36 GeV , supercluster |#eta^{subleading e}| < 2.5");
    text.DrawLatex(0.22, 0.86,  "HLT_DoubleEle33");
  }
  if (testLeg=="muLeg") {
    text.DrawLatex(0.58, 0.96, "p_{T}^{#mu} > 53 GeV , |#eta^{#mu}| < 2.4");
    text.DrawLatex(0.22, 0.86,  "HLT_Mu50");
  }
  leg->Draw();
  text.SetTextSize(0.03);
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    else if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} < 0.4");
    }
  }
  if (sample=="ttbar")  text.DrawLatex(0.2, 0.96, "MC");
  else       text.DrawLatex(0.2, 0.96, "Data");
  c->Print( Form("%s_rmin_BS.pdf",probelep.data()) );
  c->Clear();
  leg->Clear();

  c->cd();
  rminh->GetYaxis()->SetRangeUser(0.1, int(h_rmin["denom_AS"]->GetMaximum()*1.5) );
  rminh->Draw();
  h_rmin["denom_AS"]->Draw("samepE");
  h_rmin["nom_AS"]->Draw("samepE");
  leg->SetTextSize(0.036);
  leg->AddEntry(h_rmin["denom_AS"], "Before Trigger" , "pE"); 
  leg->AddEntry(h_rmin["nom_AS"], "After Trigger" , "pE");
  text.SetTextSize(0.03);
  if (testLeg=="eleLeg"){
    text.DrawLatex(0.30, 0.96, "p_{T}^{subleading e} > 36 GeV , supercluster |#eta^{subleading e}| < 2.5");
    text.DrawLatex(0.22, 0.86,  "HLT_DoubleEle33");
  }
  if (testLeg=="muLeg"){
    text.DrawLatex(0.58, 0.96, "p_{T}^{#mu} > 53 GeV , |#eta^{#mu}| < 2.4");
    text.DrawLatex(0.22, 0.86,  "HLT_Mu50");
  }
  //text.DrawLatex(0.5, 0.46, "applied selections");
  leg->Draw();
  text.SetTextSize(0.03);
  if (sample=="ttbar")  text.DrawLatex(0.2, 0.96, "MC");
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    else if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} < 0.4");
    }
  }
  if (sample=="ttbar")  text.DrawLatex(0.2, 0.96, "MC");
  else       text.DrawLatex(0.2, 0.96, "Data");
  c->Print( Form("%s_rmin_AS.pdf",probelep.data()) );
  c->Clear();
  leg->Clear();

//Efficiency plots//

  TH1F *h_pt_eff_BS, *h_pt_eff_AS;
  TH1F *h_eta_eff_BS, *h_eta_eff_AS;
  TH1F *h_rmin_eff_BS, *h_rmin_eff_AS;
  TGraphAsymmErrors* err_pt_eff_BS= new TGraphAsymmErrors();
  TGraphAsymmErrors* err_pt_eff_AS= new TGraphAsymmErrors();
  TGraphAsymmErrors* err_eta_eff_BS= new TGraphAsymmErrors();
  TGraphAsymmErrors* err_eta_eff_AS= new TGraphAsymmErrors();
  TGraphAsymmErrors* err_rmin_eff_BS= new TGraphAsymmErrors();
  TGraphAsymmErrors* err_rmin_eff_AS= new TGraphAsymmErrors();

  pTh->GetYaxis()->SetRangeUser(0, 1.1);
  pTh->GetYaxis()->SetTitle("efficiency");
  //pTh->GetXaxis()->SetRangeUser(0, 200);
  pTh->Draw();

  leg->SetTextSize(0.04);

  h_pT["denom_BS"] = (TH1F*) h_pT["denom_BS"]->Rebin(rebin_loose.size()-1, "denom_BS", &rebin_loose[0]);
  h_pT["nom_BS"] = (TH1F*) h_pT["nom_BS"]->Rebin(rebin_loose.size()-1, "nom_BS", &rebin_loose[0]);

  err_pt_eff_BS->BayesDivide(h_pT["nom_BS"],h_pT["denom_BS"]);
  text.SetTextSize(0.03);
  if (testLeg=="eleLeg"){
    text.DrawLatex(0.55, 0.96, "supercluster |#eta^{subleading e}| < 2.5");
    leg->AddEntry(err_pt_eff_BS, "HLT_DoubleEle33" , "pE");
    leg->Draw();
  }
  if (testLeg=="muLeg"){
    text.DrawLatex(0.65, 0.96, "|#eta^{#mu}| < 2.4");
    leg_eff->AddEntry(err_pt_eff_BS, "HLT_Mu50" , "pE");
    leg_eff->Draw();
  }
  err_pt_eff_BS->SetMarkerStyle(22);
  err_pt_eff_BS->Draw("samepe");
  text.SetTextSize(0.03);
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    else if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} < 0.4");
    }
  }
  if (sample=="ttbar")  text.DrawLatex(0.2, 0.96, "MC");
  else       text.DrawLatex(0.2, 0.96, "Data");
  c->Print( Form("%s_pt_eff_BS.pdf",probelep.data()) );
  c->Clear();
  leg->Clear();
  leg_eff->Clear();

  pTh->Draw();
  h_pT["denom_AS"] = (TH1F*) h_pT["denom_AS"]->Rebin(rebin_loose.size()-1, "denom_AS", &rebin_loose[0]);
  h_pT["nom_AS"] = (TH1F*) h_pT["nom_AS"]->Rebin(rebin_loose.size()-1, "nom_AS", &rebin_loose[0]);

  err_pt_eff_AS->BayesDivide(h_pT["nom_AS"],h_pT["denom_AS"]);
  text.SetTextSize(0.03);
  if (testLeg=="eleLeg"){
    leg->AddEntry(err_pt_eff_AS, "HLT_DoubleEle33" , "pE");
    text.DrawLatex(0.55, 0.96, "supercluster |#eta^{subleading e}| < 2.5");
    leg->Draw();
  }
  if (testLeg=="muLeg"){
    leg_eff->AddEntry(err_pt_eff_AS, "HLT_Mu50" , "pE");
    text.DrawLatex(0.65, 0.96, "|#eta^{#mu}| < 2.4");
    leg_eff->Draw();
  }
  err_pt_eff_AS->SetMarkerStyle(22);
  err_pt_eff_AS->Draw("samepe");
  text.SetTextSize(0.03);
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    else if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} < 0.4");
    }
  }
  if (sample=="ttbar")  text.DrawLatex(0.2, 0.96, "MC");
  else       text.DrawLatex(0.2, 0.96, "Data");
  //text.DrawLatex(0.5, 0.46, "applied selections");
  c->Print( Form("%s_pt_eff_AS.pdf",probelep.data()) );
  c->Clear();
  leg->Clear();
  leg_eff->Clear();

  h_eta["denom_BS"]->Rebin(); h_eta["nom_BS"]->Rebin();
  etah->GetYaxis()->SetRangeUser(0, 1.1);
  etah->GetYaxis()->SetTitle("efficiency");
  etah->GetXaxis()->SetRangeUser(-5, 5);
  etah->Draw();
  err_eta_eff_BS->BayesDivide(h_eta["nom_BS"],h_eta["denom_BS"]);
  err_eta_eff_BS->SetMarkerStyle(22);
  err_eta_eff_BS->Draw("samepe");
  if (testLeg=="eleLeg"){
    leg->AddEntry(err_eta_eff_BS, "HLT_DoubleEle33" , "pE");
    text.DrawLatex(0.63, 0.96, "p_{T}^{subleading e} > 36 (GeV)");
    leg->Draw();
  }
  if (testLeg=="muLeg"){
    leg_eff->AddEntry(err_eta_eff_BS, "HLT_Mu50" , "pE");
    text.DrawLatex(0.63, 0.96, "p_{T}^{#mu} > 53 (GeV)");
    leg_eff->Draw();
  }
  text.SetTextSize(0.03);
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	         text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    else if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	         text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} < 0.4");
    }
  }
  if (sample=="ttbar")  text.DrawLatex(0.2, 0.96, "MC");
  else       text.DrawLatex(0.2, 0.96, "Data");
  c->Print( Form("%s_eta_eff_BS.pdf",probelep.data()) );
  c->Clear();
  leg->Clear();
  leg_eff->Clear();


  h_eta["denom_AS"]->Rebin(); h_eta["nom_AS"]->Rebin();
  etah->Draw();
  err_eta_eff_AS->BayesDivide(h_eta["nom_AS"],h_eta["denom_AS"]);
  err_eta_eff_AS->SetMarkerStyle(22);
  err_eta_eff_AS->Draw("samepe");
  if (testLeg=="eleLeg"){
    leg->AddEntry(err_eta_eff_AS, "HLT_DoubleEle33" , "pE");
    text.DrawLatex(0.63, 0.96, "p_{T}^{subleading e} > 36 (GeV)");
    leg->Draw();
  }
  if (testLeg=="muLeg"){
    leg_eff->AddEntry(err_eta_eff_AS, "HLT_Mu50" , "pE");
    text.DrawLatex(0.63, 0.96, "p_{T}^{#mu} > 53 (GeV)");
    leg_eff->Draw();
  }
  text.SetTextSize(0.03);
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	         text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    else if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	         text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} < 0.4");
    }
  }
  if (sample=="ttbar")  text.DrawLatex(0.2, 0.96, "MC");
  else       text.DrawLatex(0.2, 0.96, "Data");
  //text.DrawLatex(0.5, 0.46, "applied selections");
  c->Print( Form("%s_eta_eff_AS.pdf",probelep.data()) );
  c->Clear();
  leg->Clear();
  leg_eff->Clear();


  vector<double> remin_loose = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,2.0,2.5,3.0,5.0};

  h_rmin["denom_BS"] = (TH1F*) h_rmin["denom_BS"]->Rebin(remin_loose.size()-1, "rmin_eff_denom_BS", &remin_loose[0]);
  h_rmin["nom_BS"]   = (TH1F*) h_rmin["nom_BS"]->Rebin(remin_loose.size()-1, "rmin_eff_nom_BS", &remin_loose[0]);
  rminh->GetYaxis()->SetRangeUser(0, 1.1);
  rminh->GetYaxis()->SetTitle("efficiency");
  rminh->GetXaxis()->SetRangeUser(-5, 5);
  rminh->Draw();
  err_rmin_eff_BS->BayesDivide(h_rmin["nom_BS"],h_rmin["denom_BS"]);
  err_rmin_eff_BS->SetMarkerStyle(22);
  err_rmin_eff_BS->Draw("samepe");
  text.SetTextSize(0.03);
  if (testLeg=="eleLeg"){
    leg->AddEntry(err_rmin_eff_BS, "HLT_DoubleEle33" , "pE");
    text.DrawLatex(0.3, 0.96, "p_{T}^{subleading e} > 36 (GeV) , supercluster |#eta^{subleading e}|<2.5");
    leg->Draw();
  }
  else if (testLeg=="muLeg"){
    leg_eff->AddEntry(err_rmin_eff_BS, "HLT_Mu50" , "pE");
    text.DrawLatex(0.58, 0.96, "p_{T}^{#mu} > 53 (GeV) , #eta|^{#mu}<2.4");
    leg_eff->Draw();
  }
  text.SetTextSize(0.03);
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    else if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} < 0.4");
    }
  }
  if (sample=="ttbar")  text.DrawLatex(0.2, 0.96, "MC");
  else       text.DrawLatex(0.2, 0.96, "Data");
  c->Print( Form("%s_rmin_eff_BS.pdf",probelep.data()) );
  c->Clear();
  leg->Clear();
  leg_eff->Clear();

  h_rmin["denom_AS"] = (TH1F*) h_rmin["denom_AS"]->Rebin(remin_loose.size()-1, "rmin_eff_denom_AS", &remin_loose[0]);
  h_rmin["nom_AS"]   = (TH1F*) h_rmin["nom_AS"]->Rebin(remin_loose.size()-1, "rmin_eff_nom_AS", &remin_loose[0]);
  rminh->Draw();
  err_rmin_eff_AS->BayesDivide(h_rmin["nom_AS"],h_rmin["denom_AS"]);
  err_rmin_eff_AS->SetMarkerStyle(22);
  err_rmin_eff_AS->Draw("samepe");
  text.SetTextSize(0.03);
  if (testLeg=="eleLeg"){
    leg->AddEntry(err_rmin_eff_AS, "HLT_DoubleEle33" , "pE");
    text.DrawLatex(0.3, 0.96, "p_{T}^{subleading e} > 36 (GeV) , supercluster |#eta^{subleading e}|<2.5");
    leg->Draw();
  }
  else if (testLeg=="muLeg"){
    leg_eff->AddEntry(err_rmin_eff_AS, "HLT_Mu50" , "pE");
    text.DrawLatex(0.58, 0.96, "p_{T}^{#mu} > 53 (GeV) , |#eta^{#mu}|<2.4");
    leg_eff->Draw();
  }
  leg->Draw();
  text.SetTextSize(0.03);
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    else if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} < 0.4");
    }
  }
  if (sample=="ttbar")  text.DrawLatex(0.2, 0.96, "MC");
  else       text.DrawLatex(0.2, 0.96, "Data");
  //text.DrawLatex(0.5, 0.46, "applied selections");
  c->Print( Form("%s_rmin_eff_AS.pdf",probelep.data()) );
  c->Clear();
  leg->Clear();
  leg_eff->Clear();

  h_eta_pT["denom_BS"]->Sumw2();
  h_eta_pT["nom_BS"]->Sumw2();

  MC_h_eta_pT["denom_BS"]->Sumw2();
  MC_h_eta_pT["nom_BS"]->Sumw2();

  TH2F* ratio_2h_loose = (TH2F*) h_eta_pT["denom_BS"]->Clone("loose");
  ratio_2h_loose->Divide(h_eta_pT["nom_BS"],h_eta_pT["denom_BS"],1,1,"b");
  TH2F* back_ratio_2h_loose = (TH2F*)  MC_h_eta_pT["denom_BS"]->Clone("loose");
  back_ratio_2h_loose->Divide(MC_h_eta_pT["nom_BS"],MC_h_eta_pT["denom_BS"],1,1,"b");


  if (testLeg=="eleLeg"){
    ratio_2h_loose->GetXaxis()->SetTitle("supercluster |#eta^{subleading e}|");
    ratio_2h_loose->GetYaxis()->SetTitle("p_{T}^{subleading e} (GeV)");
    ratio_2h_loose->GetXaxis()->SetRangeUser(0,2.5);
  }
  if (testLeg=="muLeg"){
    ratio_2h_loose->GetXaxis()->SetTitle("|#eta^{#mu}|");
    ratio_2h_loose->GetYaxis()->SetTitle("p_{T}^{#mu} (GeV)");
    ratio_2h_loose->GetXaxis()->SetRangeUser(0,2.4);
    ratio_2h_loose->GetYaxis()->SetRangeUser(50,500);
  }
  ratio_2h_loose->GetZaxis()->SetRangeUser(0,1);

  ratio_2h_loose->GetXaxis()->SetTitleOffset(1.0);
  ratio_2h_loose->GetYaxis()->SetTitleOffset(1.0);
  ratio_2h_loose->GetXaxis()->SetLabelSize(0.024);
  ratio_2h_loose->GetYaxis()->SetLabelSize(0.024);
  ratio_2h_loose->SetMarkerSize(0.65);
  gStyle->SetPalette(kBird);
  gStyle->SetPaintTextFormat("6.4f");
  ratio_2h_loose->Draw("COLZ textE");
  text.SetTextSize(0.03);
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    else if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} < 0.4");
    }
  }
  if (sample=="ttbar")  text.DrawLatex(0.2, 0.96, "MC");
  else                  text.DrawLatex(0.2, 0.96, "Data");
  c->Print( Form("%s_2dalleff.pdf",probelep.data()) );
  c->Clear();


  h_eta_pT["denom_AS"]->Sumw2();
  h_eta_pT["nom_AS"]->Sumw2();

  MC_h_eta_pT["denom_AS"]->Sumw2();
  MC_h_eta_pT["nom_AS"]->Sumw2();

  TH2F* ratio_2h_tight = (TH2F*) h_eta_pT["denom_AS"]->Clone("tight");
  ratio_2h_tight->Divide(h_eta_pT["nom_AS"],h_eta_pT["denom_AS"],1,1,"b");
  TH2F* back_ratio_2h_tight = (TH2F*)  MC_h_eta_pT["denom_AS"]->Clone("tight");
  back_ratio_2h_tight->Divide(MC_h_eta_pT["nom_AS"],MC_h_eta_pT["denom_AS"],1,1,"b");


  if (testLeg=="eleLeg"){
    ratio_2h_tight->GetXaxis()->SetTitle("supercluster |#eta^{subleading e}|");
    ratio_2h_tight->GetYaxis()->SetTitle("p_{T}^{subleading e} (GeV)");
    ratio_2h_tight->GetXaxis()->SetRangeUser(0,2.5);
  }
  else if (testLeg=="muLeg"){
    ratio_2h_tight->GetXaxis()->SetTitle("|#eta^{#mu}|");
    ratio_2h_tight->GetYaxis()->SetTitle("p_{T}^{#mu} (GeV)");
    ratio_2h_tight->GetXaxis()->SetRangeUser(0,2.4);
    ratio_2h_tight->GetYaxis()->SetRangeUser(50,500);
  }

  ratio_2h_tight->GetZaxis()->SetRangeUser(0,1);

  ratio_2h_tight->GetXaxis()->SetTitleOffset(1.0);
  ratio_2h_tight->GetYaxis()->SetTitleOffset(1.0);
  ratio_2h_tight->GetXaxis()->SetLabelSize(0.024);
  ratio_2h_tight->GetYaxis()->SetLabelSize(0.024);
  ratio_2h_tight->SetMarkerSize(0.75);
  gStyle->SetPalette(kBird);
  ratio_2h_tight->Draw("COLZ textE");
  text.SetTextSize(0.03);
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    else if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} < 0.4");
    }
  }
  if (sample=="ttbar")  text.DrawLatex(0.2, 0.96, "MC");
  else                  text.DrawLatex(0.2, 0.96, "Data");
  //text.DrawLatex(0.6, 0.96, "applied selections");
  c->Print( Form("%s_2dtighteff.pdf",probelep.data()) );
  c->Clear();


  TH2F* sf_2h_loose = (TH2F*)ratio_2h_loose->Clone("loose");
  sf_2h_loose->Divide(ratio_2h_loose,back_ratio_2h_loose,1,1,"");

  if (testLeg=="eleLeg"){
    sf_2h_loose->GetXaxis()->SetTitle("supercluster |#eta^{subleading e}|");
    sf_2h_loose->GetYaxis()->SetTitle("p_{T}^{subleading e} (GeV)");
    sf_2h_loose->GetXaxis()->SetRangeUser(0,2.5);
  }
  else if (testLeg=="muLeg"){
    sf_2h_loose->GetXaxis()->SetTitle("|#eta^{#mu}|");
    sf_2h_loose->GetYaxis()->SetTitle("p_{T}^{#mu} (GeV)");
    sf_2h_loose->GetXaxis()->SetRangeUser(0,2.4);
  }
  sf_2h_loose->GetZaxis()->SetRangeUser(0,1);

  sf_2h_loose->GetXaxis()->SetTitleOffset(1.0);
  sf_2h_loose->GetYaxis()->SetTitleOffset(1.2);
  sf_2h_loose->SetMarkerSize(0.75);

  sf_2h_loose->GetXaxis()->SetLabelSize(0.024);
  sf_2h_loose->GetYaxis()->SetLabelSize(0.024);
  sf_2h_loose->GetZaxis()->SetLabelSize(0.024);
  gStyle->SetPalette(kBird);
  sf_2h_loose->Draw("COLZ textE");
  text.SetTextSize(0.03);
  text.DrawLatex(0.2, 0.96, "Data/MC Scale Factor");
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    else if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} < 0.4");
    }
  }
  c->Print( Form("%s_2d_BS_SF.pdf",probelep.data()) );
  c->Clear();


  TH2F* sf_2h_tight = (TH2F*)ratio_2h_tight->Clone("tight");
  sf_2h_tight->Divide(ratio_2h_tight,back_ratio_2h_tight,1,1,"");

  if (testLeg=="eleLeg"){
    sf_2h_tight->GetXaxis()->SetTitle("supercluster |#eta^{subleading e}|");
    sf_2h_tight->GetYaxis()->SetTitle("p_{T}^{subleading e} (GeV)");
    sf_2h_tight->GetXaxis()->SetRangeUser(0,2.5);
  }
  if (testLeg=="muLeg"){
    sf_2h_tight->GetXaxis()->SetTitle("|#eta^{#mu}|");
    sf_2h_tight->GetYaxis()->SetTitle("p_{T}^{#mu} (GeV)");
    sf_2h_tight->GetXaxis()->SetRangeUser(0,2.4);
    sf_2h_tight->GetYaxis()->SetRangeUser(25,800);
  }
  sf_2h_tight->GetZaxis()->SetRangeUser(0,1);

  sf_2h_tight->GetXaxis()->SetTitleOffset(1.0);
  sf_2h_tight->GetYaxis()->SetTitleOffset(1.2);
  sf_2h_tight->SetMarkerSize(0.65);

  sf_2h_tight->GetXaxis()->SetLabelSize(0.024);
  sf_2h_tight->GetYaxis()->SetLabelSize(0.024);
  sf_2h_tight->GetZaxis()->SetLabelSize(0.024);
  gStyle->SetPalette(kBird);
  sf_2h_tight->Draw("COLZ textE");
  text.SetTextSize(0.03);
  text.DrawLatex(0.2, 0.96, "Data/MC Scale Factor");
  //text.DrawLatex(0.6, 0.96, "applied selections");
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    else if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} < 0.4");
    }
  }
  c->Print( Form("%s_2d_cleanSF.pdf",probelep.data()) );
  TFile* outFile = new TFile("electronTrigSF.root","RECREATE");
  if (sample=="data"){

    outFile->cd();
    sf_2h_loose->SetDrawOption("textE");
    sf_2h_loose->Write("loose_ScaleFactor");
    sf_2h_tight->SetDrawOption("textE");
    sf_2h_tight->Write("tight_ScaleFactor");
  }
  c->Clear();

  rminh->Draw();
  rminh->GetYaxis()->SetTitle("Data/MC SF");

  //h_rmin["denom_BS"] = (TH1F*) h_rmin["denom_BS"]->Rebin(remin_loose.size()-1, "rmin_clean_before", &remin_loose[0]);
  //h_rmin["nom_BS"]  = (TH1F*)  h_rmin["nom_BS"]->Rebin(remin_loose.size()-1, "rmin_clean_after", &remin_loose[0]);
  MC_h_rmin["denom_BS"] = (TH1F*) MC_h_rmin["denom_BS"]->Rebin(remin_loose.size()-1, "MC_rmin_clean_before", &remin_loose[0]);
  MC_h_rmin["nom_BS"]  = (TH1F*)  MC_h_rmin["nom_BS"]->Rebin(remin_loose.size()-1, "MC_rmin_clean_after", &remin_loose[0]);

  h_rmin["denom_BS"]->Sumw2();
  h_rmin["nom_BS"]->Sumw2();

  MC_h_rmin["denom_BS"]->Sumw2();
  MC_h_rmin["nom_BS"]->Sumw2();

  TH1F* ratio_rmin_BS = (TH1F*)h_rmin["nom_BS"]->Clone("nom_BS");
  ratio_rmin_BS->Divide(h_rmin["nom_BS"],h_rmin["denom_BS"],1,1,"b");

  TH1F* MC_ratio_rmin_BS = (TH1F*)MC_h_rmin["nom_BS"]->Clone("MCnom_BS");
  MC_ratio_rmin_BS->Divide(MC_h_rmin["nom_BS"],MC_h_rmin["denom_BS"],1,1,"b");

  TH1F* clone_ratio_rmin_BS = (TH1F*)ratio_rmin_BS->Clone("ratioco");
  clone_ratio_rmin_BS->Divide(ratio_rmin_BS, MC_ratio_rmin_BS,1,1,"");


  clone_ratio_rmin_BS->SetMarkerStyle(22);
  clone_ratio_rmin_BS->SetMarkerSize(0.5);
  clone_ratio_rmin_BS->SetMarkerColor(kBlack);
  clone_ratio_rmin_BS->SetLineColor(kBlack);

  text.SetTextSize(0.03);
  if (testLeg=="eleLeg") {
    text.DrawLatex(0.3, 0.96, "p_{T}^{subleading e} > 36 (GeV) , supercluster |#eta^{subleading e}|<2.5");
    leg->AddEntry(clone_ratio_rmin_BS, "HLT_DoubleEle33" , "pE");
  }
  else if (testLeg=="muLeg") {
    text.DrawLatex(0.5, 0.96, "p_{T}^{#mu} > 53 (GeV) , |#eta^{#mu}|<2.4");
    leg->AddEntry(clone_ratio_rmin_BS, "HLT_Mu50" , "pE");
  }

  clone_ratio_rmin_BS->Draw("samePE");
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    else if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} < 0.4");
    }
  }
  c->Print( Form("%s_rmin_BS_SF.pdf",probelep.data()) );
  c->Clear();


  rminh->Draw();
  rminh->GetYaxis()->SetTitle("Data/MC SF");

  //h_rmin["denom_AS"] = (TH1F*) h_rmin["denom_AS"]->Rebin(remin_loose.size()-1, "rmin_clean_before", &remin_loose[0]);
  //h_rmin["nom_AS"]  = (TH1F*)  h_rmin["nom_AS"]->Rebin(remin_loose.size()-1, "rmin_clean_after", &remin_loose[0]);
  MC_h_rmin["denom_AS"] = (TH1F*) MC_h_rmin["denom_AS"]->Rebin(remin_loose.size()-1, "MC_rmin_clean_before", &remin_loose[0]);
  MC_h_rmin["nom_AS"]  = (TH1F*)  MC_h_rmin["nom_AS"]->Rebin(remin_loose.size()-1, "MC_rmin_clean_after", &remin_loose[0]);

  h_rmin["denom_AS"]->Sumw2();
  h_rmin["nom_AS"]->Sumw2();

  MC_h_rmin["denom_AS"]->Sumw2();
  MC_h_rmin["nom_AS"]->Sumw2();

  TH1F* ratio_rmin_AS = (TH1F*)h_rmin["nom_AS"]->Clone("nom_AS");
  ratio_rmin_AS->Divide(h_rmin["nom_AS"],h_rmin["denom_AS"],1,1,"b");

  //cout << "ratio_data: " << "\t" << ratio_rmin_AS->GetBinContent(19) << "+_" << "\t" << ratio_rmin_AS->GetBinError(19) << endl;
  TH1F* MC_ratio_rmin_AS = (TH1F*)MC_h_rmin["nom_AS"]->Clone("MCnom_AS");
  MC_ratio_rmin_AS->Divide(MC_h_rmin["nom_AS"],MC_h_rmin["denom_AS"],1,1,"b");
  //cout << "ratio_MC: " << "\t" << MC_ratio_rmin_AS->GetBinContent(19) << "+_" << "\t" << MC_ratio_rmin_AS->GetBinError(19) << endl;

  TH1F* clone_ratio_rmin_AS = (TH1F*)ratio_rmin_AS->Clone("ratioco_AS");
  clone_ratio_rmin_AS->Divide(ratio_rmin_AS, MC_ratio_rmin_AS,1,1,"");

  TF1 *fit_SF_rmin;
  if (type==""){ fit_SF_rmin = new TF1("fit_SF_rmin","pol0",0.,5.); }
  else{
    if (type=="Onbt_")  fit_SF_rmin = new TF1("fit_SF_rmin","pol0",0.0,0.4); 
    else if (type=="Onnb_")  fit_SF_rmin = new TF1("fit_SF_rmin","pol0",0.4,5.0); 
  }
  cout << "ratio_ratio: " << "\t" << ratio_rmin_AS->GetBinContent(19) << "+_" << "\t" << ratio_rmin_AS->GetBinError(19) << endl;

  text.SetTextSize(0.03);
  if (testLeg=="eleLeg"){
    text.DrawLatex(0.3, 0.96, "p_{T}^{subleading e} > 36 (GeV) , supercluster |#eta^{subleading e}|<2.5");
    //text.DrawLatex(0.5, 0.46, "applied selections");
  }
  if (testLeg=="muLeg"){
    text.DrawLatex(0.5, 0.96, "p_{T}^{#mu} > 53 (GeV) , |#eta^{#mu}|<2.4");
    //text.DrawLatex(0.5, 0.46, "applied selections");
  }
  clone_ratio_rmin_AS->SetMarkerStyle(22);
  clone_ratio_rmin_AS->SetMarkerSize(0.5);
  clone_ratio_rmin_AS->SetMarkerColor(kBlack);
  clone_ratio_rmin_AS->SetLineColor(kBlack);
  if (testLeg=="eleLeg") leg->AddEntry(ratio_rmin_AS, "HLT_DoubleEle33" , "pE");
  else if (testLeg=="muLeg") leg->AddEntry(ratio_rmin_AS, "HLT_Mu50" , "pE");
  //cout << "Number of bins  " << ratio_rmin_AS->GetNbinsX() << endl;
  //cout << "content " << ratio_rmin_AS->GetBinContent(19)<< "  error  " << ratio_rmin_AS->GetBinError(19) << endl;
  rminh->GetYaxis()->SetRangeUser(0.0,1.2);
  clone_ratio_rmin_AS->Fit(fit_SF_rmin, "N");
  text.DrawLatex(0.2, 0.38, Form("#chi^{2}/ndof = %4.2f/%i", fit_SF_rmin->GetChisquare(), fit_SF_rmin->GetNDF() ) );
  text.DrawLatex(0.2, 0.32, Form("Prob = %4.3f", fit_SF_rmin->GetProb() ) );
  text.DrawLatex(0.2, 0.26, Form("p0 = %4.3f #pm %4.3f", fit_SF_rmin->GetParameter(0), fit_SF_rmin->GetParError(0) ) );
  clone_ratio_rmin_AS->Draw("samePE");
  if (type!=""){
    if (testLeg=="eleLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{subleading e} < 0.4");
    }
    else if (testLeg=="muLeg"){
      if (type=="Onnb_") text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} > 0.4");
      else 	        text.DrawLatex(0.2, 0.75, "#DeltaR_{min}^{#mu} < 0.4");
    }
  }
  c->Print( Form("%s_rmin_AS_SF.pdf",probelep.data()) );
  c->Clear();

  
}




void setStyle(){

//Style//

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.16);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.04, "XYZ");
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.1);

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  //tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  tdrStyle->SetPaperSize(20.,20.);

  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->SetOptStat(0);

  tdrStyle->cd();

//End Style//
}

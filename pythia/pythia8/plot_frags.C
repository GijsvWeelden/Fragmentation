
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TLegend.h"

#include "./plotUtils.C"

//-------------------------------------------------------------
//
// Plot pythia fragmentation
//
//-------------------------------------------------------------

string format_string(string hadron, double min_pt, double max_pt);
void plot_FFs(TH1F* hIncl, TH1F* hg, TH1F* hq, string hadron, double min_pt, double max_pt);

void plot_frags(void){
  double time = clock();
  gROOT->SetBatch();
  string inName = "PythiaResult";
  TFile *inFile = TFile::Open(TString::Format("./%s.root", inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  std::vector<double> jetPtBins = {20.,40.,60.,80.,100.,120.,160.,200.};
  std::vector<string> Hadrons = {"pi", "K", "p", "pi0", "K0L", "K0S", "K0", "Lambda0"};
	std::vector<string> Partons = {"g", "q"};
  TH2F *hInclJet, *hgJet, *hqJet;

  TH2F* hJetEtaPt = (TH2F*) inFile->Get("hJetEtaPt");
  int nJets = hJetEtaPt->GetEntries();

  for (int iHad = 0; iHad < Hadrons.size(); iHad++){
    std::cout << "Plotting " << Hadrons[iHad] << std::endl;
    hInclJet = (TH2F*) inFile->Get(TString::Format("hJetFrag_%s", Hadrons[iHad].c_str()).Data());
    hgJet = (TH2F*) inFile->Get(TString::Format("hgFrags_%s", Hadrons[iHad].c_str()).Data());
    hqJet = (TH2F*) inFile->Get(TString::Format("hqFrags_%s", Hadrons[iHad].c_str()).Data());

    // Should be normalised by nJets. Should maybe happen in filling macro.
    hInclJet->Scale(1./nJets);
    hgJet->Scale(1./nJets);
    hqJet->Scale(1./nJets);

    double min_pt = 80, max_pt = 100;
    hInclJet->GetYaxis()->SetRangeUser(min_pt, max_pt);
    hgJet->GetYaxis()->SetRangeUser(min_pt, max_pt);
    hqJet->GetYaxis()->SetRangeUser(min_pt, max_pt);

    TH1F* hIncl = (TH1F*) hInclJet->ProjectionX();
    TH1F* hg = (TH1F*) hgJet->ProjectionX();
    TH1F* hq = (TH1F*) hqJet->ProjectionX();
  }

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
  // Do not close the file, it will delete the ratio histogram
  inFile->Close();
}

void plot_FFs(TH1F* hIncl, TH1F* hg, TH1F* hq, string hadron, double min_pt, double max_pt){
  hIncl->SetStats(0);
  hIncl->SetLineColor(GetColor(0));
  hIncl->SetMarkerColor(GetColor(0));
  hIncl->SetMarkerStyle(GetMarker(0));

  hg->SetStats(1);
  hg->SetLineColor(GetColor(0));
  hg->SetMarkerColor(GetColor(1));
  hg->SetMarkerStyle(GetMarker(1));

  hq->SetStats(2);
  hq->SetLineColor(GetColor(0));
  hq->SetMarkerColor(GetColor(2));
  hq->SetMarkerStyle(GetMarker(2));

  auto legend = CreateLegend(0.25, 0.45, 0.6, 0.9, "", 0.05);
  legend->AddEntry(hIncl, "incl.");
  legend->AddEntry(hg, "g");
  legend->AddEntry(hq, "q");

  string yTitle = format_string(hadron, min_pt, max_pt);
  TH1F* frame1 = DrawFrame(0., 1, 0, 3, "", yTitle);
  frame1->Draw();
  hIncl->Draw("same");
  hg->Draw("same");
  hq->Draw("same");
  legend->Draw("same");
  return;
}

string format_string(string hadron, double min_pt, double max_pt){
  string had, s;
    if (hadron == "pi"){
    had = "#pi^{#pm}";
  }
  else if (hadron == "pi0"){
    had = "#pi^{0}";
  }
  else if (hadron == "K0L"){
    had = "K^{0}_{L}";
  }
  else if (hadron == "K0S"){
    had = "K^{0}_{S}";
  }
  else if (hadron == "K0"){
    had = "K^{0}";
  }
  else if (hadron == "Lambda0"){
    had = "#Lambda^{0}";
  }
  // TODO: Does this work?
  s = TString::Format("#it{D}^{%s}(#it{z}) | #it{p}_{T}^{jet} #in (%.0f, %.0f) GeV", had.c_str(), min_pt, max_pt).Data();
  return s;
}
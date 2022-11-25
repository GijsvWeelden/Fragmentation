
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TError.h"
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

string format_hadron_name(string hadron);
void prep_single_hadron(TH2F* hInclJet, TH2F* hgJet, TH2F* hqJet, TH2F* hNJetTypes,
                        TH1F* &hIncl, TH1F* &hg, TH1F* &hq, TH1F* &hNgq,
                        double min_pt, double max_pt
                        );
void load_hists(TFile* inFile, TH2F* &h0, TH2F* &h1, TH2F* &h2, TH2F* &h3,
               string hadron0, string hadron1, string hadron2, string hadron3,
               string IGQ
               );
void prep_mult_hadron(TH2F* h0_2D, TH2F* h1_2D, TH2F* h2_2D, TH2F* h3_2D, TH2F* hNJetTypes,
                      TH1F* &h0_1D, TH1F* &h1_1D, TH1F* &h2_1D, TH1F* &h3_1D,
                      string hadron0, string hadron1, string hadron2, string hadron3,
                      double min_pt, double max_pt, string IGQ
                      );
void plot_FFs_single(TH1F* hIncl, TH1F* hg, TH1F* hq, string hadron, double min_pt, double max_pt);
void plot_FFs_mult(TH1F* h0, TH1F* h1, TH1F* h2, TH1F* h3,
                   string hadron0, string hadron1, string hadron2, string hadron3,
                   double min_pt, double max_pt, string IGQ
                   );
void plot_mult_hadrons(TFile* inFile, TH2F* h0_2D, TH2F* h1_2D, TH2F* h2_2D, TH2F* h3_2D, TH2F* hNJetTypes,
                       TH1F* &h0_1D, TH1F* &h1_1D, TH1F* &h2_1D, TH1F* &h3_1D,
                       string hadron0, string hadron1, string hadron2, string hadron3,
                       std::vector<double> jetPtBins, string IGQ
                       );
int hadrons_to_kkp(string hadron);
void set_kkp_params(TF1* &kkp, double E, int flavour, int hadron);
int prep_kkp(TF1* kkp, string IGQ, string hadron, double E);
void plot_matchDist(TH1F* dist, std::vector<double> matchDist);

void plot_frags(void)
{
  double time = clock();
  gROOT->SetBatch();
  // gErrorIgnoreLevel = kWarning; // Suppress info messages
  bool plotSingleHadron = false, plotMultHadron = true, plotMatchDist = false;
  string inName = "PythiaResultJob12128342_66_pthat80_200";
  TFile *inFile = TFile::Open(TString::Format("../data/%s.root", inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  std::vector<double> jetPtBins = {40.,60.,80.,100.,120.,160.,200.};
  std::vector<string> Hadrons = {"pi", "K", "p", "pi0", "K0L", "K0S", "K0", "Lambda0"};
	std::vector<string> Partons = {"g", "q"};

  TH2F* hNJetTypes = (TH2F*) inFile->Get("hNJetTypes");
  TH1F* hNgqJets = (TH1F*) hNJetTypes->ProjectionX();

  if (plotSingleHadron){
    TH2F *hInclJet, *hgJet, *hqJet;
    TH1F *hIncl, *hg, *hq, *hNgq;

    for (int iHad = 0; iHad < Hadrons.size(); iHad++){
      std::cout << "Plotting " << Hadrons[iHad] << std::endl;
      hInclJet = (TH2F*) inFile->Get(TString::Format("hJetFrag_%s", Hadrons[iHad].c_str()).Data());
      hgJet = (TH2F*) inFile->Get(TString::Format("hgFrags_%s", Hadrons[iHad].c_str()).Data());
      hqJet = (TH2F*) inFile->Get(TString::Format("hqFrags_%s", Hadrons[iHad].c_str()).Data());

      prep_single_hadron(hInclJet, hgJet, hqJet, hNJetTypes, hIncl, hg, hq, hNgq, jetPtBins.front(), jetPtBins.back());
      plot_FFs_single(hIncl, hg, hq, Hadrons[iHad], jetPtBins.front(), jetPtBins.back());

      for (int ipt = 1; ipt < jetPtBins.size(); ipt++){
        double min_pt = jetPtBins[ipt - 1], max_pt = jetPtBins[ipt];
        prep_single_hadron(hInclJet, hgJet, hqJet, hNJetTypes, hIncl, hg, hq, hNgq, min_pt, max_pt);
        plot_FFs_single(hIncl, hg, hq, Hadrons[iHad], min_pt, max_pt);
      }
    }
  } // Single hadron plotting
  else cout << "Skipping single hadron fragmentation." << endl;
  if (plotMultHadron){
    TH2F *h0_2D, *h1_2D, *h2_2D, *h3_2D;
    TH1F *h0_1D, *h1_1D, *h2_1D, *h3_1D;
    string hadron0, hadron1, hadron2, hadron3;

    if (true){
      hadron0 = Hadrons[0];
      hadron1 = Hadrons[1];
      hadron2 = Hadrons[2];
      hadron3 = "";
      plot_mult_hadrons(inFile, h0_2D, h1_2D, h2_2D, h3_2D, hNJetTypes,
                        h0_1D, h1_1D, h2_1D, h3_1D,
                        hadron0, hadron1, hadron2, hadron3,
                        jetPtBins, "g");
      plot_mult_hadrons(inFile, h0_2D, h1_2D, h2_2D, h3_2D, hNJetTypes,
                        h0_1D, h1_1D, h2_1D, h3_1D,
                        hadron0, hadron1, hadron2, hadron3,
                        jetPtBins, "q");
      // plot_mult_hadrons(inFile, h0_2D, h1_2D, h2_2D, h3_2D, hNJetTypes,
      //                   h0_1D, h1_1D, h2_1D, h3_1D,
      //                   hadron0, hadron1, hadron2, hadron3,
      //                   jetPtBins, "i");
    }

    if (false){
      hadron0 = Hadrons[0];
      hadron1 = Hadrons[1];
      hadron2 = Hadrons[2];
      hadron3 = Hadrons[7];

      plot_mult_hadrons(inFile, h0_2D, h1_2D, h2_2D, h3_2D, hNJetTypes,
                        h0_1D, h1_1D, h2_1D, h3_1D,
                        hadron0, hadron1, hadron2, hadron3,
                        jetPtBins, "g");
      plot_mult_hadrons(inFile, h0_2D, h1_2D, h2_2D, h3_2D, hNJetTypes,
                        h0_1D, h1_1D, h2_1D, h3_1D,
                        hadron0, hadron1, hadron2, hadron3,
                        jetPtBins, "q");
      plot_mult_hadrons(inFile, h0_2D, h1_2D, h2_2D, h3_2D, hNJetTypes,
                        h0_1D, h1_1D, h2_1D, h3_1D,
                        hadron0, hadron1, hadron2, hadron3,
                        jetPtBins, "i");
    }
  } // Multiple hadron plotting
  else cout << "Skipping multiple hadron fragmentation comparison." << endl;
  if (plotMatchDist){
    TH1F* hDeltaPartonJet = (TH1F*) inFile->Get("hDeltaPartonJet");
    std::vector<double> matchDist = {0.3, 0.5, 0.7};
    plot_matchDist(hDeltaPartonJet, matchDist);
  }

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
  // Do not close the file, it will delete the ratio histogram
  // inFile->Close();
}

string format_hadron_name(string hadron)
{
  string had = hadron;
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
  else if (hadron == "K"){
    had = "K^{#pm}";
  }
  else if (hadron == "Lambda0"){
    had = "#Lambda^{0}";
  }
  return had;
}

void plot_matchDist(TH1F* dist, std::vector<double> matchDist)
{
  TCanvas *c_matchDist = new TCanvas("matchDist", "Distance between parton and jet", 1600, 1000);
  dist->Scale(1./dist->Integral());
  dist->SetLineColor(GetColor(2));
  dist->SetMarkerColor(GetColor(2));
  TH1F* distClone = (TH1F*) dist->Clone();
  double ymax = dist->GetBinContent(dist->GetMaximumBin());
  ymax *= 1.1;
  auto legend = CreateLegend(0.5, 0.95, 0.6, 0.9, "", 0.05);
  TH1F* frame = DrawFrame(0, 1., 0, ymax, "#Delta", "");

  distClone->GetXaxis()->SetRangeUser(0, matchDist[2]);
  double frac3 = distClone->Integral();
  TLine* L3 = new TLine(matchDist[2], 0, matchDist[2], distClone->GetBinContent(distClone->FindBin(matchDist[2])));
  L3->SetLineColor(GetColor(3));
  legend->AddEntry(L3, TString::Format("#Delta: %.1f, frac: %.1f", matchDist[2], frac3).Data());

  distClone->GetXaxis()->SetRangeUser(0, matchDist[1]);
  double frac2 = distClone->Integral();
  TLine* L2 = new TLine(matchDist[1], 0, matchDist[1], distClone->GetBinContent(distClone->FindBin(matchDist[1])));
  L2->SetLineColor(GetColor(1));
  legend->AddEntry(L2, TString::Format("#Delta: %.1f, frac: %.1f", matchDist[1], frac2).Data());

  distClone->GetXaxis()->SetRangeUser(0, matchDist[0]);
  double frac1 = distClone->Integral();
  TLine* L1 = new TLine(matchDist[0], 0, matchDist[0], distClone->GetBinContent(distClone->FindBin(matchDist[0])));
  L1->SetLineColor(GetColor(0));
  legend->AddEntry(L1, TString::Format("#Delta: %.1f, frac: %.1f", matchDist[0], frac1).Data());

  frame->Draw();
  dist->Draw("same hist");
  L1->Draw("same");
  L2->Draw("same");
  L3->Draw("same");
  legend->Draw("same");
  c_matchDist->SaveAs(TString::Format("../plots/matchDist.pdf"));
  return;
}

void prep_single_hadron(TH2F* hInclJet, TH2F* hgJet, TH2F* hqJet, TH2F* hNJetTypes,
                        TH1F* &hIncl, TH1F* &hg, TH1F* &hq, TH1F* &hNgq,
                        double min_pt, double max_pt)
{
  hInclJet->GetYaxis()->SetRangeUser(min_pt, max_pt);
  hgJet->GetYaxis()->SetRangeUser(min_pt, max_pt);
  hqJet->GetYaxis()->SetRangeUser(min_pt, max_pt);
  hNJetTypes->GetYaxis()->SetRangeUser(min_pt, max_pt);

  hIncl = (TH1F*) hInclJet->ProjectionX();
  hg = (TH1F*) hgJet->ProjectionX();
  hq = (TH1F*) hqJet->ProjectionX();
  hNgq = (TH1F*) hNJetTypes->ProjectionX();

  int nGluons_Pt = (int) hNgq->GetBinContent(1);
  int nQuarks_Pt = (int) hNgq->GetBinContent(2);
  int nJets_Pt = nGluons_Pt + nQuarks_Pt;

  hIncl->Scale(1./nJets_Pt, "width");
  hg->Scale(1./nGluons_Pt, "width");
  hq->Scale(1./nQuarks_Pt, "width");
  return;
}

void prep_mult_hadron(TH2F* h0_2D, TH2F* h1_2D, TH2F* h2_2D, TH2F* h3_2D, TH2F* hNJetTypes,
                      TH1F* &h0_1D, TH1F* &h1_1D, TH1F* &h2_1D, TH1F* &h3_1D,
                      string hadron0, string hadron1, string hadron2, string hadron3,
                      double min_pt, double max_pt, string IGQ)
{
  hNJetTypes->GetYaxis()->SetRangeUser(min_pt, max_pt);
  TH1F* hNgq = (TH1F*) hNJetTypes->ProjectionX();
  int nGluons = (int) hNgq->GetBinContent(1);
  int nQuarks = (int) hNgq->GetBinContent(2);
  int nJets = nGluons + nQuarks;
  int norm = -1;

  if (IGQ == "g"){
    norm = nGluons;
  }
  else if (IGQ == "q"){
    norm = nQuarks;
  }
  else if (IGQ == "i"){
    norm = nJets;
  }
  else{
    cout << "Error: invalid value for IGQ. Aborting histogram prep." << endl;
    return;
  }

  if (hadron0 != ""){
    h0_2D->GetYaxis()->SetRangeUser(min_pt, max_pt);
    h0_1D = (TH1F*) h0_2D->ProjectionX();
    h0_1D->Scale(1./norm, "width");
  }
  if (hadron1 != ""){
    h1_2D->GetYaxis()->SetRangeUser(min_pt, max_pt);
    h1_1D = (TH1F*) h1_2D->ProjectionX();
    h1_1D->Scale(1./norm, "width");
  }
  if (hadron2 != ""){
    h2_2D->GetYaxis()->SetRangeUser(min_pt, max_pt);
    h2_1D = (TH1F*) h2_2D->ProjectionX();
    h2_1D->Scale(1./norm, "width");
  }
  if (hadron3 != ""){
    h3_2D->GetYaxis()->SetRangeUser(min_pt, max_pt);
    h3_1D = (TH1F*) h3_2D->ProjectionX();
    h3_1D->Scale(1./norm, "width");
  }
  return;
}

void load_hists(TFile* inFile, TH2F* &h0, TH2F* &h1, TH2F* &h2, TH2F* &h3,
                string hadron0, string hadron1, string hadron2, string hadron3,
                string IGQ)
{
  string histName;
  if (IGQ == "i"){
    histName = "JetFrag";
  }
  else if (IGQ == "g"){
    histName = "gFrags";
  }
  else if (IGQ == "q"){
    histName = "gFrags";
  }
  if (hadron0 != "") h0 = (TH2F*) inFile->Get(TString::Format("h%s_%s", histName.c_str(), hadron0.c_str()).Data());
  if (hadron1 != "") h1 = (TH2F*) inFile->Get(TString::Format("h%s_%s", histName.c_str(), hadron1.c_str()).Data());
  if (hadron2 != "") h2 = (TH2F*) inFile->Get(TString::Format("h%s_%s", histName.c_str(), hadron2.c_str()).Data());
  if (hadron3 != "") h3 = (TH2F*) inFile->Get(TString::Format("h%s_%s", histName.c_str(), hadron3.c_str()).Data());
  return;
}

// Plot Fragmentation Function for single hadron species. Plot inclusive, gluon, and quark jets
void plot_FFs_single(TH1F* hIncl, TH1F* hg, TH1F* hq,
                     string hadron, double min_pt, double max_pt)
{
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

  auto legend = CreateLegend(0.75, 0.95, 0.6, 0.9, "", 0.05);
  legend->AddEntry(hIncl, "incl.");
  legend->AddEntry(hg, "g");
  legend->AddEntry(hq, "q");

  TCanvas *c_had = new TCanvas(TString::Format("c_%s_pt_%.0f_%.0f",
                                               hadron.c_str(), min_pt, max_pt).Data(),
                               TString::Format("%s_pt_%.0f_%.0f",
                                               hadron.c_str(), min_pt, max_pt).Data(),
                               1600, 1000);
  c_had->SetLogy();

  // cout << "Plot_FFs_single after setup" << endl;

  double xmin = 0., xmax = 1., ymin = 1e-5, ymax;
  ymax = max({hIncl->GetMaximum(), hg->GetMaximum(), hq->GetMaximum()});
  ymax *= 1.1;
  string s = format_hadron_name(hadron);
  string yTitle = TString::Format("#frac{1}{#it{N}_{jets}} #frac{d#it{N}_{%s}}{d#it{z}}", s.c_str()).Data();
  TH1F* frame1 = DrawFrame(xmin, xmax, ymin, ymax, "#it{z}", yTitle);
  frame1->SetTitle(TString::Format("%s in jets (%.0f, %.0f) GeV", s.c_str(), min_pt, max_pt).Data());

  // cout << "Plot_FFs_single constructed frame" << endl;

  frame1->Draw();
  hIncl->Draw("same");
  hg->Draw("same");
  hq->Draw("same");
  legend->Draw("same");
  // cout << "Plot_FFs_single after draw" << endl;
  c_had->SaveAs(TString::Format("../plots/SingleHadrons/%s_pt_%.0f_%.0f.pdf", hadron.c_str(), min_pt, max_pt).Data());
  return;
}

void plot_FFs_mult(TH1F* h0, TH1F* h1, TH1F* h2, TH1F* h3,
                   string hadron0, string hadron1, string hadron2, string hadron3,
                   double min_pt, double max_pt, string IGQ)
{
  auto legend = CreateLegend(0.75, 0.95, 0.6, 0.9, "", 0.05);
  double xmin = 0., xmax = 1., ymin = 1e-5, ymax;
  int kkp_int0, kkp_int1, kkp_int2, kkp_int3;
  const float E = (min_pt + max_pt) / 2; //50;

  TF1* kkp0 = new TF1("kkp0", kkp_func, 0, 1, 3);
  TF1* kkp1 = new TF1("kkp1", kkp_func, 0, 1, 3);
  TF1* kkp2 = new TF1("kkp2", kkp_func, 0, 1, 3);
  TF1* kkp3 = new TF1("kkp3", kkp_func, 0, 1, 3);

  if (hadron0 != ""){
    h0->SetStats(0);
    h0->SetLineColor(GetColor(0));
    h0->SetMarkerColor(GetColor(0));
    h0->SetMarkerStyle(GetMarker(0));
    string had0 = format_hadron_name(hadron0);
    legend->AddEntry(h0, had0.c_str());
    ymax = h0->GetMaximum();
    kkp_int0 = prep_kkp(kkp0, IGQ, hadron0, E);
    kkp0->SetLineColor(GetColor(0));
    legend->AddEntry(kkp0, TString::Format("KKP %s #rightarrow %s (%.0f GeV)", IGQ.c_str(), hadron0.c_str(), E).Data());
  }
  if (hadron1 != ""){
    h1->SetStats(1);
    h1->SetLineColor(GetColor(1));
    h1->SetMarkerColor(GetColor(1));
    h1->SetMarkerStyle(GetMarker(1));
    string had1 = format_hadron_name(hadron1);
    legend->AddEntry(h1, had1.c_str());
    ymax = max({h0->GetMaximum(), h1->GetMaximum()});
    kkp_int1 = prep_kkp(kkp1, IGQ, hadron1, E);
    kkp1->SetLineColor(GetColor(1));
    legend->AddEntry(kkp1, TString::Format("KKP %s #rightarrow %s (%.0f GeV)", IGQ.c_str(), hadron1.c_str(), E).Data());
  }
  if (hadron2 != ""){
    h2->SetStats(2);
    h2->SetLineColor(GetColor(2));
    h2->SetMarkerColor(GetColor(2));
    h2->SetMarkerStyle(GetMarker(2));
    string had2 = format_hadron_name(hadron2);
    legend->AddEntry(h2, had2.c_str());
    ymax = max({h0->GetMaximum(), h1->GetMaximum(), h2->GetMaximum()});
    kkp_int2 = prep_kkp(kkp2, IGQ, hadron2, E);
    kkp2->SetLineColor(GetColor(2));
    legend->AddEntry(kkp2, TString::Format("KKP %s #rightarrow %s (%.0f GeV)", IGQ.c_str(), hadron2.c_str(), E).Data());
  }
  if (hadron3 != ""){
    h3->SetStats(3);
    h3->SetLineColor(GetColor(3));
    h3->SetMarkerColor(GetColor(3));
    h3->SetMarkerStyle(GetMarker(3));
    string had3 = format_hadron_name(hadron3);
    legend->AddEntry(h3, had3.c_str());
    ymax = max({h0->GetMaximum(), h1->GetMaximum(), h2->GetMaximum(), h3->GetMaximum()});
    kkp_int3 = prep_kkp(kkp3, IGQ, hadron3, E);
    kkp3->SetLineColor(GetColor(3));
    legend->AddEntry(kkp3, TString::Format("KKP %s #rightarrow %s (%.0f GeV)", IGQ.c_str(), hadron3.c_str(), E).Data());
  }

  TCanvas *c_mult = new TCanvas(TString::Format("c_mult_%s_pt_%.0f_%.0f",
                                                IGQ.c_str(), min_pt, max_pt).Data(),
                                TString::Format("mult_%s_pt_%.0f_%.0f",
                                                IGQ.c_str(), min_pt, max_pt).Data(),
                                1600, 1000);
  c_mult->SetLogy();

  ymax *= 1.1;
  string yTitle = TString::Format("#frac{1}{#it{N}_{jets}^{%s}} #frac{d#it{N}_{h}}{d#it{z}}", IGQ.c_str()).Data();
  TH1F* frame1 = DrawFrame(xmin, xmax, ymin, ymax, "#it{z}", yTitle);

  string str;
  if (IGQ == "i"){
    str = "inclusive";
  }
  else if (IGQ == "g"){
    str = "gluon";
  }
  else if (IGQ == "q"){
    str = "quark";
  }
  else{
    cout << "Error: invalid value for IGQ. Aborting histogram plot." << endl;
    return;
  }
  frame1->SetTitle(TString::Format("Hadrons in %s jets (%.0f, %.0f) GeV", str.c_str(), min_pt, max_pt).Data());

  frame1->Draw();
  if (hadron0 != ""){
    h0->Draw("same");
    if (kkp_int0 == 0) kkp0->Draw("same");
  }
  if (hadron1 != ""){
    h1->Draw("same");
    if (kkp_int1 == 0) kkp1->Draw("same");
  }
  if (hadron2 != ""){
    h2->Draw("same");
    if (kkp_int2 == 0) kkp2->Draw("same");
  }
  if (hadron3 != ""){
    h3->Draw("same");
    if (kkp_int3 == 0) kkp3->Draw("same");
  }
  legend->Draw("same");

  if (hadron0 != "") hadron0 = hadron0.append("_");
  if (hadron1 != "") hadron1 = hadron1.append("_");
  if (hadron2 != "") hadron2 = hadron2.append("_");
  if (hadron3 != "") hadron3 = hadron3.append("_");

  c_mult->SaveAs(TString::Format("../plots/MultipleHadrons/%s_%s%s%s%spt_%.0f_%.0f.pdf",
                                 str.c_str(), hadron0.c_str(), hadron1.c_str(), hadron2.c_str(), hadron3.c_str(),
                                 min_pt, max_pt).Data()
                                 );
  return;
}

int hadrons_to_kkp(string hadron){
  int kkp = -1;
  if (hadron == "pi"){
    kkp = 1;
  }
  else if (hadron == "K"){
    kkp = 2;
  }
  else if (hadron == "K0"){
    kkp = 3;
  }
  else if (hadron == "p"){
    kkp = 4;
  }
  else if (hadron == "pi0"){
    kkp = 5;
  }
  else if (hadron == "n"){
    kkp = 6;
  }
  else if (hadron == "h"){
    kkp = 7;
  }
  return kkp;
}

void set_kkp_params(TF1* &kkp, double E, int flavour, int hadron){
  kkp->SetParameter(0, E); // Q for fragmentation
  kkp->SetParameter(1, flavour); // parton flavour
  kkp->SetParameter(2, hadron); // hadron type
  return;
}

int prep_kkp(TF1* kkp, string IGQ, string hadron, double E){
  int flavour, had;
  if (IGQ == "g"){
    flavour = 0;
  }
  else if (IGQ == "q"){
    flavour = 1;
  }
  else{
    cout << "Not quark or gluon KKP requested. Skipping." << endl;
    return -1;
  }
  had = hadrons_to_kkp(hadron);
  if (had < 0){
    cout << "KKP requested for non-KKP hadron. Skipping." << endl;
    return -1;
  }

  set_kkp_params(kkp, E, flavour, had);
  return 0;
}

void plot_mult_hadrons(TFile* inFile, TH2F* h0_2D, TH2F* h1_2D, TH2F* h2_2D, TH2F* h3_2D, TH2F* hNJetTypes,
                       TH1F* &h0_1D, TH1F* &h1_1D, TH1F* &h2_1D, TH1F* &h3_1D,
                       string hadron0, string hadron1, string hadron2, string hadron3,
                       std::vector<double> jetPtBins, string IGQ
                       )
{
    cout << "Loading hists" << endl;
    load_hists(inFile, h0_2D, h1_2D, h2_2D, h3_2D,
               hadron0, hadron1, hadron2, hadron3, IGQ);
    cout << "Prepping hists" << endl;
    prep_mult_hadron(h0_2D, h1_2D, h2_2D, h3_2D, hNJetTypes,
                     h0_1D, h1_1D, h2_1D, h3_1D,
                     hadron0, hadron1, hadron2, hadron3,
                     jetPtBins.front(), jetPtBins.back(), IGQ);
    cout << "Plotting FFs" << endl;
    plot_FFs_mult(h0_1D, h1_1D, h2_1D, h3_1D,
                  hadron0, hadron1, hadron2, hadron3,
                  jetPtBins.front(), jetPtBins.back(), IGQ);

    for (int ipt = 1; ipt < jetPtBins.size(); ipt++){
      double min_pt = jetPtBins[ipt - 1], max_pt = jetPtBins[ipt];
      prep_mult_hadron(h0_2D, h1_2D, h2_2D, h3_2D, hNJetTypes,
                       h0_1D, h1_1D, h2_1D, h3_1D,
                       hadron0, hadron1, hadron2, hadron3,
                       min_pt, max_pt, IGQ);
      plot_FFs_mult(h0_1D, h1_1D, h2_1D, h3_1D,
                    hadron0, hadron1, hadron2, hadron3,
                    min_pt, max_pt, IGQ);
    }
}

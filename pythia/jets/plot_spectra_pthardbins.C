#ifdef __CLING__
#include "utils.C"
#endif

void plot_spectra_pthardbins(void) {
#ifndef __CLING__
  gROOT->LoadMacro("~/macros/utils.C");
#endif
  gStyle->SetLineWidth(2);
  gStyle->SetOptStat(0);
  set_legend_font(62,0.04);

  //Char_t *fname_base="output/pythia_full_jets_13000GeV_%.0f_%.0f_all.root";
  //Char_t *fname_base="output/pythia_jets_13000GeV_%.0f_%.0f_all.root";
  const Char_t *fname_base="output/pythia_charged_jet_spectra_5020GeV_%.0f_%.0f_all.root";
  //Char_t *fname_base="output/pythia_charged_jet_spectra_13000GeV_%.0f_%.0f_all.root";
  //const Char_t *fname_base="output/pythia_charged_jet_spectra_P11_altbins_7000GeV_%.0f_%.0f_all.root";
  //Char_t *hname = "hJetPtEtaPhi";
  const Char_t *hname = "hJetPtEtaPhi_R2";
  //Char_t *fname_base="output/pythia_full_jet_spectra_5020GeV_%.0f_%.0f_all.root";
  //Char_t *hname = "hJetPtEtaPhi_R2";
  //Float_t pthard[]={5,15,25,40,100,250};
  Float_t pthard[]={5,11,21,36,57,84,117,152,191,234,300};
  //Float_t pthard[]={5,11,25,60,120,160,200,500};
  //Float_t pthard[]={5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70,84,117,152,191,234,300};

  const Int_t do_scaling = 1;
  const Float_t max_eta = 0.5;
  const Float_t ymin = 1e-10;

  Int_t nbins = sizeof(pthard)/sizeof(Float_t)-1;
  cout << nbins << endl;
 
  TCanvas *c1= new TCanvas("c1","c1: pt spectra",1200,800);
  c1->Divide(2,1);

  TLatex *ltx = new TLatex;
  ltx->SetNDC();

  TH1D *pt_jet_sum = 0;
  TH1F *pt_pion_sum = 0;

  for (Int_t ibin = 0; ibin < nbins; ibin++) {
    TFile *fin = new TFile(Form(fname_base,pthard[ibin],pthard[ibin+1]));
    auto hXSec = dynamic_cast<TProfile*>(fin->Get("hXSec"));
    auto hNEvent = dynamic_cast<TH1F*>(fin->Get("hNEvent"));
    cout << "pthard " << pthard[ibin] << " - " << pthard[ibin+1]<< " xsec " << hXSec->GetBinContent(1) << " events " << hNEvent->GetBinContent(1) << endl;
    Float_t scale = hXSec->GetBinContent(1)/hNEvent->GetBinContent(1);

    c1->cd(1);
    gPad->SetLogy();
    auto hPtPion = dynamic_cast<TH1F*>(fin->Get("hPtPion"));
    
    if (do_scaling) {
      hPtPion->Sumw2();
      hPtPion->Scale(scale);
    }
    if (ibin == 0) {      
      pt_pion_sum = new TH1F(*hPtPion);
      hPtPion->Draw();
      hPtPion->SetMinimum(ymin);
      hPtPion->GetXaxis()->SetRange(1,125);
      gPad->Update();
      ltx->DrawLatex(ndc_x(0.6),ndc_y(0.94),"p_{T,hard}");
    }
    else {
      pt_pion_sum->Add(hPtPion);
      hPtPion->SetLineColor(1+ibin);
      hPtPion->SetMarkerColor(1+ibin);
      hPtPion->Draw("same");
    }

    auto hptp = hPtPion;

    gPad->Update();
    draw_legend_l(0.6,0.88,hptp,Form("%.0f - %.0f",pthard[ibin],pthard[ibin+1]),ibin);

    c1->cd(2);
    gPad->SetLogy();
    auto hJetPtEtaPhi = (TH3F*)  fin->Get(hname);
    //cout << "File " << fin->GetName() << " hJetPtEtaPhi " << hJetPtEtaPhi << endl;
    Int_t minetabin = hJetPtEtaPhi->GetYaxis()->FindBin(-max_eta+0.001);
    Int_t maxetabin = hJetPtEtaPhi->GetYaxis()->FindBin(max_eta-0.001);
    hJetPtEtaPhi->GetYaxis()->SetRange(minetabin, maxetabin);
    auto px = dynamic_cast<TH1D*>(hJetPtEtaPhi->Project3D(Form("x_%d",ibin)));
    px->SetName(Form("px_%d",ibin));
    if (do_scaling) {
      px->Sumw2();
      px->Scale(scale);
    }
    if (ibin == 0) {
      auto pxc = px->DrawCopy();
      pxc->SetMinimum(ymin);
      pxc->GetXaxis()->SetRange(1,175);
      pt_jet_sum = new TH1D(*px);
      pt_jet_sum->SetName("pt_jet_sum");
    }
    else {
      px->SetLineColor(1+ibin);
      px->SetMarkerColor(1+ibin);
      px->DrawCopy("same");
      pt_jet_sum->Add(px);
    }
  }
  c1->cd(1);
  pt_pion_sum->SetLineWidth(2);
  pt_pion_sum->SetMarkerStyle(20);
  pt_pion_sum->SetMarkerSize(0.5);
  pt_pion_sum->Draw("same");

  c1->cd(2);
  pt_jet_sum->SetLineWidth(2);
  pt_jet_sum->SetMarkerStyle(20);
  pt_jet_sum->SetMarkerSize(0.5);
  pt_jet_sum->Draw("same");

  Int_t max_ptbin = 150;
  Int_t min_ptbin[] = {41, 61, 101};
  Int_t n_ptint = sizeof(min_ptbin)/sizeof(min_ptbin[0]);
  for (Int_t ipt = 0; ipt < n_ptint; ipt++) {
    cout << "Integral " << pt_jet_sum->GetBinLowEdge(min_ptbin[ipt]) << " " << pt_jet_sum->GetXaxis()->GetBinUpEdge(max_ptbin) << " " << pt_jet_sum->Integral(min_ptbin[ipt],max_ptbin) << endl;
  }
}

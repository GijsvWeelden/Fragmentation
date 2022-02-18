{
  Float_t min_pt = 60; //100; //60; //100;
  Float_t max_pt = 80; //120; //80; //120;
  Float_t max_y = 10; //9;
  //const Char_t *declust_label="anti-k_{T}";
  //const Char_t *hname="hJetPtEtaZgAKtDeclust_R04";
  const Char_t *declust_label="k_{T}";
  const Char_t *hname="hJetPtEtaZgKtDeclust_R04";
  //const Char_t *declust_label="CA";
  //const Char_t *hname="hJetPtEtaZg_R04";
  gStyle->SetOptStat(0);
  gROOT->LoadMacro("jewel_plot_utils.C");

  TCanvas *c1 = new TCanvas("c1","c1: Canvas 1", 600,600);

  //TFile *fin_vac = TFile::Open("analysis/results_2017/jet_shapes_pp_2tev76_nobkgsub_full.root");
  //TFile *fin_med = TFile::Open("analysis/results_2017/jet_shapes_AA_Tdef_0cent10_recoil_constsub_eventwise_full.root");
  // These files had the argument for softdrop swapped: beta = 0.1, zcut = 0.
  //TFile *fin_vac = TFile::Open("analysis/results_2017/jet_shapes_pp_2tev76_eventwise_nobkgsub_charged.root");
  //TFile *fin_med = TFile::Open("analysis/results_2017/jet_shapes_AA_Tdef_0cent10_recoil_constsub_eventwise_charged.root");
  TFile *fin_vac = TFile::Open("analysis/results_2017/jet_shapes_constsub_eventwise_charged_nobkg_pp_2tev76.root");
  TFile *fin_med = TFile::Open("analysis/results_2017/jet_shapes_constsub_eventwise_charged_AA_Tdef_0cent10_recoil.root");
  
  hJetPtEtaZgVac = (TH3F*) fin_vac->Get(hname);
  hZgVac = project_shape(hJetPtEtaZgVac, min_pt, max_pt, "hZgVac",1);
  hZgVac->GetXaxis()->SetRangeUser(0.1,0.5);
  hZgVac->Scale(1./hZgVac->Integral(hZgVac->FindBin(0.101),hZgVac->FindBin(0.499))/hZgVac->GetBinWidth(2));
  hZgVac->SetMarkerColor(4);
  hZgVac->SetMinimum(0);
  hZgVac->SetMaximum(max_y);
  hZgVac->Draw();

  hJetPtEtaZgMed = (TH3F*) fin_med->Get(hname);
  hZgMed = project_shape(hJetPtEtaZgMed, min_pt, max_pt, "hZgMed",1);
  hZgMed->GetXaxis()->SetRangeUser(0.1,0.5);
  hZgMed->Scale(1./hZgMed->Integral(hZgMed->FindBin(0.101),hZgMed->FindBin(0.499))/hZgMed->GetBinWidth(2));
  hZgMed->SetMarkerColor(2);
  hZgMed->SetLineColor(2);
  hZgMed->Draw("same");

  TLegend *leg = new TLegend(0.58,0.4,0.85,0.6);
  leg->AddEntry(hZgVac,"pp","pl");
  leg->AddEntry(hZgMed,"Pb-Pb","pl");
  leg->Draw();

  TLatex *ltx = new TLatex;
  ltx->SetNDC();
  ltx->SetTextFont(42);
  ltx->DrawLatex(0.25,0.83,"JEWEL");
  ltx->DrawLatex(0.25,0.77,Form("%.0f < p_{T,jet} < %.0f GeV/c",min_pt,max_pt));
  ltx->DrawLatex(0.25,0.71,declust_label);
  ltx->DrawLatex(0.25,0.65,"charged jets");

  TCanvas *c2 = new TCanvas("c2","c2: Canvas 2", 600,600);
  hRatio = new TH1D(*((TH1D*)hZgMed));
  hRatio->Divide(hZgVac);
  hRatio->Draw();

}


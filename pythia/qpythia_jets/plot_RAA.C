{
  Float_t min_pt = 20;
  Float_t max_pt = 40;
  Float_t max_deta = 0.2; // also used for other projection, max_dphi
  Float_t Rval = 0.2;

  TFile *fin_qpythia_vac = new TFile("pythia_charged_jet_spectra_2760GeV_0_merged.root");
  TFile *fin_qpythia_med = new TFile("pythia_charged_jet_spectra_2760GeV_50_merged.root");

  TCanvas *c1 = new TCanvas("c1","c1: spectra",500,800);
  c1->cd();
  gPad->SetLogy();

  hJetPtEtaPhi_vac = (TH3F*) fin_qpythia_vac->Get(Form("hJetPtEtaPhi_R%d",Int_t(Rval*10)));
  hJetPt_vac = hJetPtEtaPhi_vac->Project3D("x_vac");
  hJetPt_vac->Draw();
 

  hJetPtEtaPhi_med = (TH3F*) fin_qpythia_med->Get(Form("hJetPtEtaPhi_R%d",Int_t(Rval*10)));
  hJetPt_med = hJetPtEtaPhi_med->Project3D("x_med");
  hJetPt_med->SetLineColor(2);
  hJetPt_med->Draw("same");


  TCanvas *c2 = new TCanvas("c2","c2: RAA",600,400);
  c2->cd();

  raa = new TH1D(*((TH1D*)hJetPt_med));
  raa->SetName("raa_jet_h");
  raa->Divide(hJetPt_vac);
  raa->SetYTitle("R_{AA}");
  raa->SetMinimum(0);
  raa->SetMaximum(1);
  raa->Draw();
}

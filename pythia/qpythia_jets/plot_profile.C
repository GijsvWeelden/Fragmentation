{
  Float_t min_pt = 100;
  Float_t max_pt = 120;
  Float_t max_deta = 0.3; // also used for other projection, max_dphi
  Int_t Rint = 3;
  Float_t Rval = Rint/10.;

  TFile *fin_qpythia_vac = new TFile("pythia_charged_jet_spectra_2760GeV_0_merged.root");
  TFile *fin_qpythia_med = new TFile("pythia_charged_jet_spectra_2760GeV_50_merged.root");
  //TFile *fin_qpythia = new TFile("output/pythia_charged_jet_spectra_2760GeV_40_100_0_all.root");
  //TFile *fin_qpythia = new TFile("output/pythia_charged_jet_spectra_2760GeV_40_100_0_9145417.root");
  //TFile *fin_qpythia = new TFile("qPythiaJetSpectra.root");
  
  TCanvas *c1 = new TCanvas("c1","c1: spectra");
  h3DSpecPythia_vac = (TH3F*) fin_qpythia_vac->Get(Form("hLeadJetPtEtaPhi_R%d",Rint));
  hJetPtPythia_vac = h3DSpecPythia_vac->Project3D("x_vac");
  hJetPtPythia_vac->Draw();
  gPad->SetLogy();

  h3DSpecPythia_med = (TH3F*) fin_qpythia_med->Get(Form("hLeadJetPtEtaPhi_R%d",Rint));
  hJetPtPythia_med = h3DSpecPythia_med->Project3D("x_med");
  hJetPtPythia_med->SetLineColor(2);
  hJetPtPythia_med->Draw("same");


  TCanvas *c2 = new TCanvas("c2","c2: deta, dphi",1000,500);
  c2->Divide(2,1);
  TCanvas *c3 = new TCanvas("c3","c3: ratios",1000,500);
  c3->Divide(2,1);

  c2->cd(1);
  //h3DPythia = (TH3F*) fin_qpythia->Get("hPtJetDeltaEtaDeltaPhiN_R2");


  h3DPythia_med = (TH3F*) fin_qpythia_med->Get(Form("hPtJetDeltaEtaDeltaPhiPt_R%d",Rint));

  h3DPythia_med->GetXaxis()->SetRange( h3DPythia_med->GetXaxis()->FindBin(min_pt+0.001),  h3DPythia_med->GetXaxis()->FindBin(max_pt-0.001)); 
  h3DPythia_med->GetYaxis()->SetRange( h3DPythia_med->GetYaxis()->FindBin(-max_deta+0.001),  h3DPythia_med->GetYaxis()->FindBin(max_deta-0.001)); 
  h2=h3DPythia_med->Project3D("z_med");
  Float_t njet_pythia_med = hJetPtPythia_med->Integral(hJetPtPythia_med->FindBin(min_pt+0.001),hJetPtPythia_med->FindBin(max_pt-0.001)); 
  cout << "njet_pythia medium " << njet_pythia_med << endl;
  h2->Scale(1./ njet_pythia_med / h2->GetBinWidth(1));
  h2->GetXaxis()->SetRangeUser(-0.5,0.5);
  h2->Draw();
  h2->SetLineColor(2);
  cout << "PYTHIA dphi integral in R " << h2->Integral( h2->FindBin(-Rval+0.001),h2->FindBin(Rval-0.001))*h2->GetBinWidth(1) << endl;

  dphi_ratio = new TH1D(*((TH1D*)h2));
  dphi_ratio->SetName("dphi_ratio");

  h3DPythia_vac = (TH3F*) fin_qpythia_vac->Get(Form("hPtJetDeltaEtaDeltaPhiPt_R%d",Rint));
  h3DPythia_vac->GetXaxis()->SetRange( h3DPythia_vac->GetXaxis()->FindBin(min_pt+0.001),  h3DPythia_vac->GetXaxis()->FindBin(max_pt-0.001)); 
  h3DPythia_vac->GetYaxis()->SetRange( h3DPythia_vac->GetYaxis()->FindBin(-max_deta+0.001),  h3DPythia_vac->GetYaxis()->FindBin(max_deta-0.001)); 
  h2=h3DPythia_vac->Project3D("z_vac");
  Float_t njet_pythia_vac = hJetPtPythia_vac->Integral(hJetPtPythia_vac->FindBin(min_pt+0.001),hJetPtPythia_vac->FindBin(max_pt-0.001)); 
  cout << "njet_pythia vacuum " << njet_pythia_vac << endl;
  h2->Scale(1./ njet_pythia_vac / h2->GetBinWidth(1));
  h2->Draw("same");
  cout << "PYTHIA dphi integral in R " << h2->Integral( h2->FindBin(-Rval+0.001),h2->FindBin(Rval-0.001))*h2->GetBinWidth(1) << endl;

  c3->cd(1);
  dphi_ratio->Divide(h2);
  dphi_ratio->Draw();
  
  c2->cd(2);

  h3DPythia_med->GetYaxis()->SetRange( 1,h3DPythia_med->GetYaxis()->GetNbins());
  h3DPythia_med->GetZaxis()->SetRange(h3DPythia_med->GetZaxis()->FindBin(-max_deta+0.001),  h3DPythia_med->GetZaxis()->FindBin(max_deta-0.001)); 
  h2=h3DPythia_med->Project3D("y_med");
  h2->Scale(1./ njet_pythia_med / h2->GetBinWidth(1));
  h2->SetLineColor(2);
  h2->Draw("same");
  cout << "PYTHIA deta integral in R " << h2->Integral(47,54)*h2->GetBinWidth(1) << endl;
  h2->GetXaxis()->SetRangeUser(-0.5,0.5);

  deta_ratio = new TH1D(*((TH1D*)h2));
  deta_ratio->SetName("deta_ratio");

  h3DPythia_vac->GetYaxis()->SetRange( 1,h3DPythia_vac->GetYaxis()->GetNbins());
  h3DPythia_vac->GetZaxis()->SetRange(h3DPythia_vac->GetZaxis()->FindBin(-max_deta+0.001),  h3DPythia_vac->GetZaxis()->FindBin(max_deta-0.001)); 
  h2=h3DPythia_vac->Project3D("y_vac");
  h2->Scale(1./ njet_pythia_vac / h2->GetBinWidth(1));
  h2->Draw("same");
  cout << "PYTHIA deta integral in R " << h2->Integral(47,54)*h2->GetBinWidth(1) << endl;

  c3->cd(2);
  deta_ratio->Divide(h2);
  deta_ratio->Draw();
}

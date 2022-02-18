Float_t etamin = -1;
Float_t etamax = 1;

TH1* project_shape(TH3F* pt_eta_shape, Float_t ptmin, Float_t ptmax, const Char_t *hname, Int_t nrebin = 2, Int_t normalise = 1) {
  Int_t ptminbin = pt_eta_shape->GetXaxis()->FindBin(ptmin + 0.001);
  Int_t ptmaxbin = pt_eta_shape->GetXaxis()->FindBin(ptmax - 0.001);
  
  Int_t etaminbin = pt_eta_shape->GetYaxis()->FindBin(etamin + 0.001);
  Int_t etamaxbin = pt_eta_shape->GetYaxis()->FindBin(etamax - 0.001);
  pt_eta_shape->GetXaxis()->SetRange(ptminbin, ptmaxbin);
  pt_eta_shape->GetYaxis()->SetRange(etaminbin, etamaxbin);
  TH1 *hshape = pt_eta_shape->Project3D("z");
  hshape->UseCurrentStyle();
  hshape->SetName(hname);
  hshape->Rebin(nrebin);
  hshape->SetLineWidth(2);
  hshape->SetMarkerStyle(20);
  if (normalise && hshape->Integral() != 0)
    hshape->Scale(1./hshape->Integral(),"width");

  return hshape;
}

TH1* project_shape_2D(TH2F* pt_shape, Float_t ptmin, Float_t ptmax, const Char_t *tag, Int_t nrebin = 2) {
  Int_t ptminbin = pt_shape->GetXaxis()->FindBin(ptmin + 0.001);
  Int_t ptmaxbin = pt_shape->GetXaxis()->FindBin(ptmax - 0.001);
  TH1 *hshape = pt_shape->ProjectionY(Form("mass_%s",tag),ptminbin,ptmaxbin);
  hshape->Rebin(nrebin);
  hshape->SetLineWidth(2);
  hshape->SetMarkerStyle(20);
  if (hshape->Integral() != 0)
    hshape->Scale(1./hshape->Integral(),"width");
  return hshape;
}

TH2* project_mass_vs_moments(TH3F* pt_shape, Float_t ptmin, Float_t ptmax, const char *tag="", Int_t nrebin = 2) {
  //
  // projects yz, with a selection on pt (x axis)
  //  main use case: histo with pt, mass, angularity/moment
  //  output: 2D histo
  //
  Int_t ptminbin = pt_shape->GetXaxis()->FindBin(ptmin + 0.001);
  Int_t ptmaxbin = pt_shape->GetXaxis()->FindBin(ptmax - 0.001);
  cout << "ptminbin " << ptminbin << " max " << ptmaxbin << endl;
  pt_shape->GetXaxis()->SetRange(ptminbin, ptmaxbin);
  cout << "Range set, projecting " << endl;
  TH2* hshape = (TH2*) pt_shape->Project3D("yz");
  hshape->SetName(Form("%s_%s",hshape->GetName(),tag));
  cout << "projection done " << endl;
  return hshape;
}

TH1* project_shape_from_3D(TH3F* pt_shape, Float_t ptmin, Float_t ptmax, const char *tag="", Int_t nrebin = 2) {
  //
  // projects z axis, with a selection on pt (x axis)
  //  main use case: histo with pt, mass, angularity/moment
  //
  Int_t ptminbin = pt_shape->GetXaxis()->FindBin(ptmin + 0.001);
  Int_t ptmaxbin = pt_shape->GetXaxis()->FindBin(ptmax - 0.001);
  pt_shape->GetXaxis()->SetRange(ptminbin, ptmaxbin);
  TH1* hshape = pt_shape->Project3D("z");
  hshape->SetName(Form("%s_%s",hshape->GetName(),tag));
  hshape->Rebin(nrebin);
  return hshape;
}

TH1* get_jet_mass(TH3F* pt_eta_mass, Float_t ptmin, Float_t ptmax, const Char_t *tag, Int_t nrebin = 2) {
  Int_t ptminbin = pt_eta_mass->GetXaxis()->FindBin(ptmin + 0.001);
  Int_t ptmaxbin = pt_eta_mass->GetXaxis()->FindBin(ptmax - 0.001);
  
  Int_t etaminbin = pt_eta_mass->GetYaxis()->FindBin(etamin + 0.001);
  Int_t etamaxbin = pt_eta_mass->GetYaxis()->FindBin(etamax - 0.001);
  pt_eta_mass->GetXaxis()->SetRange(ptminbin, ptmaxbin);
  pt_eta_mass->GetYaxis()->SetRange(etaminbin, etamaxbin);
  TH1 * hmass = pt_eta_mass->Project3D("z");
  hmass->SetName(Form("mass_%s",tag));
  hmass->Rebin(nrebin);
  hmass->SetLineWidth(2);
  hmass->SetMarkerStyle(20);
  hmass->Scale(1./hmass->Integral(),"width");

  return hmass;
}

TH1* get_jet_mass_2D(TH2F* pt_mass, Float_t ptmin, Float_t ptmax, const Char_t *tag, Int_t nrebin = 2) {
  Int_t ptminbin = pt_mass->GetXaxis()->FindBin(ptmin + 0.001);
  Int_t ptmaxbin = pt_mass->GetXaxis()->FindBin(ptmax - 0.001);
  TH1* hmass = pt_mass->ProjectionY(Form("mass_%s",tag),ptminbin,ptmaxbin);
  hmass->Rebin(nrebin);
  hmass->SetLineWidth(2);
  hmass->SetMarkerStyle(20);
  hmass->Scale(1./hmass->Integral(),"width");
  return hmass;
}

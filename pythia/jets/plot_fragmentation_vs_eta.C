{
  gROOT->LoadMacro("~/macros/utils.C");
  set_legend_font(62,0.05);

  //Char_t *fname_base="output/pythia_jets_%.0f_%.0f_all.root";
  Char_t *fname_base="output/pythia_full_jet_ff_5020GeV_%.0f_%.0f_all.root";
  Float_t pthard[]={5,15,25,40,100};

  const Int_t do_scaling = 1;

  Int_t npthardbins = sizeof(pthard)/sizeof(Float_t)-1;
  cout << npthardbins << endl;
 
  TLatex *ltx = new TLatex;
  ltx->SetNDC();

  Float_t ptjetmin = 15;
  Float_t ptjetmax = 20;
  Float_t etajetmax[] = {1,2,3,4,5};
  Int_t netabins = sizeof(etajetmax)/sizeof(etajetmax[0])-1;


  TH3F *hsum_pion = 0;
  TH3F *hsum_gamma = 0;
  TH3F *hsum_jet = 0;

  for (Int_t ipthbin = 0; ipthbin < npthardbins; ipthbin++) {
    TFile *fin = new TFile(Form(fname_base,pthard[ipthbin],pthard[ipthbin+1]));
    cout << "pthard " << pthard[ipthbin] << " - " << pthard[ipthbin+1]<< " xsec " << hXSec->GetBinContent(1) << " events " << hNEvent->GetBinContent(1) << endl;
    Float_t scale = hXSec->GetBinContent(1)/hNEvent->GetBinContent(1);

    hEtaPtJetPtPion = (TH3F*) fin->Get("hEtaPtJetPtPion");
    hEtaPtJetPtPion->Sumw2();
    hEtaPtJetPtPion->Scale(scale);

    hEtaPtJetPtGamma = (TH3F*) fin->Get("hEtaPtJetPtGamma");
    hEtaPtJetPtGamma->Sumw2();
    hEtaPtJetPtGamma->Scale(scale);

    hJetPtEtaPhi = (TH3F*) fin->Get("hJetPtEtaPhi");
    hJetPtEtaPhi->Sumw2();
    hJetPtEtaPhi->Scale(scale);

  pt_jet_alleta = hJetPtEtaPhi->Project3D("x_all");
  cout << "Jet yld 20 GeV (all eta)" << pt_jet_alleta->GetBinContent(20) << endl;
  pt_pion_alleta = hEtaPtJetPtPion->Project3D("z_all");
  cout << "Pion yld 20 GeV (all eta)" << pt_pion_alleta->GetBinContent(20) << endl;

    if (ipthbin == 0) {
      hsum_jet = new TH3F(*hJetPtEtaPhi);
      hsum_jet->SetName("hsum_jet");

      hsum_pion = new TH3F(*hEtaPtJetPtPion);
      hsum_pion->SetName("hsum_pion");
      hsum_gamma = new TH3F(*hEtaPtJetPtGamma);
      hsum_gamma->SetName("hsum_gamma");
    }
    else {
      hsum_jet->Add(hJetPtEtaPhi);
      hsum_pion->Add(hEtaPtJetPtPion);
      hsum_gamma->Add(hEtaPtJetPtGamma);
    }
  }


  TLegend *leg = new TLegend(0.15,0.15,0.3,0.4,"#eta #in","brNDC");
  leg->SetTextSize(0.04);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);

  THStack *hjets_all = new THStack("hjets_all","jets");
  THStack *hpions_all = new THStack("hpions_all","pions");

  pt_jet_alleta = hsum_jet->Project3D("x_all");
  cout << "Jet yld 20 GeV (all eta)" << pt_jet_alleta->GetBinContent(20) << endl;
  pt_pion_alleta = hsum_pion->Project3D("z_all");
  cout << "Pion yld 20 GeV (all eta)" << pt_pion_alleta->GetBinContent(20) << endl;

  for (Int_t ietabin = 0; ietabin < netabins; ietabin++) {
    Int_t minetabin = hsum_jet->GetYaxis()->FindBin(etajetmax[ietabin]+0.001);
    Int_t maxetabin = hsum_jet->GetYaxis()->FindBin(etajetmax[ietabin+1]-0.001);
    cout << "minetabin " << minetabin << " max " << maxetabin << endl;
    hsum_jet->GetYaxis()->SetRange(minetabin, maxetabin);

    px = hsum_jet->Project3D(Form("x_%d",ietabin));
    px->SetLineColor(ietabin+1);
    px->SetMarkerColor(ietabin+1);
    hjets_all->Add(px);
    cout << "Jet yld 20 GeV " << px->GetBinContent(20) << endl;
    
    leg->AddEntry(px,Form("%.1f - %.2f",etajetmax[ietabin],etajetmax[ietabin+1]),"l");

    hsum_pion->GetXaxis()->SetRange(minetabin, maxetabin);
    cout << "check eta range " << hsum_pion->GetXaxis()->GetBinLowEdge(minetabin) << " " << hsum_pion->GetXaxis()->GetBinUpEdge(maxetabin) << endl;

    px = hsum_pion->Project3D(Form("z_all_%d",ietabin));
    cout << "pion yld 20 GeV " << px->GetBinContent(20) << endl;
    px->SetLineColor(ietabin+1);
    px->SetMarkerColor(ietabin+1);
    hpions_all->Add(px);
  }

  TCanvas *c4 = new TCanvas("c4","c4: jets and pions",800,600);
  c4->Divide(2,1);
  c4->cd(1);
  gPad->SetLogy();
  hjets_all->Draw("nostack");
  hjets_all->GetYaxis()->SetTitle("d#sigma/dp_{T} (mb c/GeV)");
  hjets_all->GetXaxis()->SetTitle("p_{T,jet} (GeV/c}");

  leg->Draw();

  c4->cd(2);
  gPad->SetLogy();
  hpions_all->Draw("nostack");
  hpions_all->GetYaxis()->SetTitle("d#sigma/dp_{T} (mb c/GeV)");
  hpions_all->GetXaxis()->SetTitle("p_{T} (GeV/c}");
  
  Int_t netabinshist = hsum_jet->GetYaxis()->GetNbins();
  hsum_jet->GetYaxis()->SetRange(1, netabinshist);

  Int_t minptbin = hsum_pion->GetYaxis()->FindBin(ptjetmin+0.0001);
  Int_t maxptbin = hsum_pion->GetYaxis()->FindBin(ptjetmax-0.0001);

  hsum_pion->GetYaxis()->SetRange(minptbin, maxptbin);
  hsum_gamma->GetYaxis()->SetRange(minptbin, maxptbin);
  hsum_jet->GetXaxis()->SetRange(minptbin, maxptbin);

  hjet_eta = hsum_jet->Project3D("y_jets");

  THStack *hpions = new THStack("hpions","pions in jets");
  THStack *hgamma = new THStack("hgamma","gammas in jets");

  // fragmentation functions
  for (Int_t ietabin = 0; ietabin < netabins; ietabin++) {
    Int_t minetabin = hsum_pion->GetXaxis()->FindBin(etajetmax[ietabin]+0.001);
    Int_t maxetabin = hsum_pion->GetXaxis()->FindBin(etajetmax[ietabin+1]-0.001);
    cout << "minetabin " << minetabin << " max " << maxetabin << endl;
    hsum_pion->GetXaxis()->SetRange(minetabin, maxetabin);

    px = hsum_pion->Project3D(Form("z_%d",ietabin));
    px->Scale(1./hjet_eta->Integral(minetabin, maxetabin));
    px->SetLineColor(ietabin+1);
    px->SetMarkerColor(ietabin+1);
    hpions->Add(px);

    hsum_gamma->GetXaxis()->SetRange(minetabin, maxetabin);

    px = hsum_gamma->Project3D(Form("z_%d",ietabin));
    px->Scale(1./hjet_eta->Integral(minetabin, maxetabin));
    px->SetLineColor(ietabin+1);
    px->SetMarkerColor(ietabin+1);
    hgamma->Add(px);

  }
  TCanvas *c1 = new TCanvas("c1","pt",800,600);
  c1->Divide(2,1);
    
  c1->cd(1);
  gPad->SetLogy();
  hpions->Draw("nostack");
  hpions->GetXaxis()->SetRangeUser(0,ptjetmax);
  hpions->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hpions->GetYaxis()->SetTitle("1/N_{jet} dN/dp_{T} (c/GeV)");
  ltx->DrawLatex(0.5,0.8,Form("%.0f < p_{T,jet} < %.f GeV",ptjetmin,ptjetmax));
  ltx->DrawLatex(0.7,0.7,"#pi^{#pm}");
  leg->Draw();

  c1->cd(2);
  gPad->SetLogy();
  hgamma->Draw("nostack");
  hgamma->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hgamma->GetYaxis()->SetTitle("1/N_{jet} dN/dp_{T} (c/GeV)");
  hgamma->GetXaxis()->SetRangeUser(0,ptjetmax);
  ltx->DrawLatex(0.5,0.8,Form("%.0f < p_{T,jet} < %.f GeV",ptjetmin,ptjetmax));
  ltx->DrawLatex(0.7,0.7,"#gamma");


  // inverted

  Int_t nptbin = hsum_pion->GetYaxis()->GetNbins();

  hsum_pion->GetYaxis()->SetRange(1, nptbin);
  hsum_gamma->GetYaxis()->SetRange(1, nptbin);
  hsum_jet->GetXaxis()->SetRange(1, nptbin);

  hjet_eta = hsum_jet->Project3D("y_jets");

  THStack *hpions_jets = new THStack("hpions_jets","jets from pions");
  THStack *hgamma_jets = new THStack("hgamma_jets","jets from gammas");

  Float_t min_ptpi = 5;
  ptpibin = hsum_pion->GetZaxis()->FindBin(min_ptpi+0.001);
  nptbin = hsum_pion->GetZaxis()->GetNbins();
  hsum_pion->GetZaxis()->SetRange(ptpibin, nptbin);
  hsum_gamma->GetZaxis()->SetRange(ptpibin, nptbin);
  // jet distributions
  for (Int_t ietabin = 0; ietabin < netabins; ietabin++) {
    Int_t minetabin = hsum_pion->GetXaxis()->FindBin(etajetmax[ietabin]+0.001);
    Int_t maxetabin = hsum_pion->GetXaxis()->FindBin(etajetmax[ietabin+1]-0.001);
    cout << "minetabin " << minetabin << " max " << maxetabin << endl;
    hsum_pion->GetXaxis()->SetRange(minetabin, maxetabin);

    px = hsum_pion->Project3D(Form("y_%d",ietabin));
    px->SetLineColor(ietabin+1);
    px->SetMarkerColor(ietabin+1);
    hpions_jets->Add(px);

    hsum_gamma->GetXaxis()->SetRange(minetabin, maxetabin);

    px = hsum_gamma->Project3D(Form("y_%d",ietabin));
    px->SetLineColor(ietabin+1);
    px->SetMarkerColor(ietabin+1);
    hgamma_jets->Add(px);

  }

  TCanvas *c2 = new TCanvas("c2","c2: jet pt",800,600);
  c2->Divide(2,1);
    
  c2->cd(1);
  gPad->SetLogy();
  hpions_jets->Draw("nostack");
  hpions_jets->GetYaxis()->SetTitle("d#sigma/dp_{T,jet} (mb c/GeV)");
  hpions_jets->GetXaxis()->SetTitle("p_{T,jet} (GeV/c}");
  ltx->DrawLatex(0.5,0.8,Form("p_{T,#pi} > %.f GeV",min_ptpi));
  ltx->DrawLatex(0.7,0.7,"#pi^{#pm}");

  leg->Draw();
  c2->cd(2);
  gPad->SetLogy();
  hgamma_jets->Draw("nostack");
  hgamma_jets->GetYaxis()->SetTitle("d#sigma/dp_{T,jet} (mb c/GeV)");
  hgamma_jets->GetXaxis()->SetTitle("p_{T,jet} (GeV/c}");
  ltx->DrawLatex(0.5,0.8,Form("p_{T,#gamma} > %.f GeV",min_ptpi));
  ltx->DrawLatex(0.7,0.7,"#gamma");


  TCanvas *c3 = new TCanvas("c3","c3: eta",600,600); 
  hjet_eta->Draw();
}

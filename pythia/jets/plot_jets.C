{
  gROOT->LoadMacro("~/macros/utils.C");
  set_legend_font(62,0.05);

  Char_t *fname_base="output/pythia_jets_%.0f_%.0f_all.root";
  Float_t pthard[]={5,15,25,40,100,250};

  const Int_t do_scaling = 1;

  Int_t nbins = sizeof(pthard)/sizeof(Float_t)-1;
  cout << nbins << endl;
 
  TCanvas *c1= new TCanvas("c1","c1: pt spectra",1000,600);
  c1->Divide(2,1);

  TCanvas *c2 = new TCanvas("c2","c2: frag funcs",800,800);
  c2->Divide(2,2);

  TCanvas *c3 = new TCanvas("c3","c3: frag func in z",500,500);
  TLatex *ltx = new TLatex;
  ltx->SetNDC();

  Float_t ptjetmin = 15;
  Float_t ptjetmax = 20;

  THStack *hpions = new THStack("hpions","pions in jets");
  THStack *hkaons = new THStack("hkaons","kaons in jets");
  THStack *hprotons = new THStack("hprotons","protons in jets");
  THStack *hdzero = new THStack("hdzero","DZero in jets");

  THStack *hdzero_z = new THStack("hdzero_z","DZero in jets");

  for (Int_t ibin = 0; ibin < nbins; ibin++) {
    TFile *fin = new TFile(Form(fname_base,pthard[ibin],pthard[ibin+1]));
    cout << "pthard " << pthard[ibin] << " - " << pthard[ibin+1]<< " xsec " << hXSec->GetBinContent(1) << " events " << hNEvent->GetBinContent(1) << endl;
    Float_t scale = hXSec->GetBinContent(1)/hNEvent->GetBinContent(1);

    c1->cd(1);
    gPad->SetLogy();
    if (do_scaling) {
      hPtPion->Sumw2();
      hPtPion->Scale(scale);
    }
    if (ibin == 0) {      
      hPtPion->Draw();
      gPad->Update();
      ltx->DrawLatex(ndc_x(0.6),ndc_y(0.86),"p_{T,hard}");
    }
    else {
      hPtPion->SetLineColor(1+ibin);
      hPtPion->SetMarkerColor(1+ibin);
      hPtPion->Draw("same");
    }

    hptp = hPtPion;

    gPad->Update();
    draw_legend_l(0.6,0.8,hptp,Form("%.0f - %.0f",pthard[ibin],pthard[ibin+1]),ibin);

    c1->cd(2);
    gPad->SetLogy();
    px = hJetPtEtaPhi->Project3D(Form("x_%d"));
    px = new TH1D(*((TH1D*)px));
    px->SetName(Form("px_%d"));
    if (do_scaling) {
      px->Sumw2();
      px->Scale(scale);
    }
    if (ibin == 0)
      px->DrawCopy();
    else {
      px->SetLineColor(1+ibin);
      px->SetMarkerColor(1+ibin);
      px->DrawCopy("same");
    }


    Int_t minbin = hPtJetPtPion->GetXaxis()->FindBin(ptjetmin+0.0001);
    Int_t maxbin = hPtJetPtPion->GetXaxis()->FindBin(ptjetmax-0.0001);

    hpy = hPtJetPtPion->ProjectionY(Form("pion_%d",ibin),minbin,maxbin);
    if (do_scaling) {
      hpy->Sumw2();
      hpy->Scale(scale);
    }
    hpy->SetLineColor(1+ibin);
    hpy->SetMarkerColor(1+ibin);
    hpions->Add(hpy);

    hpy = hPtJetPtKaon->ProjectionY(Form("kaon_%d",ibin),minbin,maxbin);
    if (do_scaling) {
      hpy->Sumw2();
      hpy->Scale(scale);
    }
    hpy->SetLineColor(1+ibin);
    hpy->SetMarkerColor(1+ibin);
    hkaons->Add(hpy);

    hpy = hPtJetPtProton->ProjectionY(Form("proton_%d",ibin),minbin,maxbin);
    if (do_scaling) {
      hpy->Sumw2();
      hpy->Scale(scale);
    }
    hpy->SetLineColor(1+ibin);
    hpy->SetMarkerColor(1+ibin);
    hprotons->Add(hpy);

    hpy = hPtJetPtD->ProjectionY(Form("dzero_%d",ibin),minbin,maxbin);
    if (do_scaling) {
      hpy->Sumw2();
      hpy->Scale(scale);
    }
    hpy->SetLineColor(1+ibin);
    hpy->SetMarkerColor(1+ibin);
    hdzero->Add(hpy);

    hpy = hPtJetZD->ProjectionY(Form("dzero_z_%d",ibin),minbin,maxbin);
    if (do_scaling) {
      hpy->Sumw2();
      hpy->Scale(scale);
    }
    hpy->SetLineColor(1+ibin);
    hpy->SetMarkerColor(1+ibin);
    hdzero_z->Add(hpy);
  }

  c2->cd(1);
  gPad->SetLogy();
  hpions->Draw("nostack");
  ltx->DrawLatex(0.5,0.8,Form("%.0f < p_{T,jet} < %.f GeV",ptjetmin,ptjetmax));
  ltx->DrawLatex(0.7,0.7,"#pi^{#pm}");

  c2->cd(2);
  gPad->SetLogy();

  hkaons->Draw("nostack");
  ltx->DrawLatex(0.7,0.8,"K^{#pm}");

  c2->cd(3);
  gPad->SetLogy();
  hprotons->Draw("nostack");
  ltx->DrawLatex(0.7,0.8,"p+#bar{p}");

  c2->cd(4);
  gPad->SetLogy();
  hdzero->Draw("nostack");
  ltx->DrawLatex(0.7,0.8,"D^{0}+#bar{D}^{0}");
  ltx->DrawLatex(0.5,0.75,Form("p_{T,jet} %.0f - %.0f GeV",ptjetmin,ptjetmax));

  c3->cd();
  gPad->SetLogy();
  hdzero_z->Draw("nostack");
  ltx->DrawLatex(0.7,0.8,"D^{0}+#bar{D}^{0}");
  ltx->DrawLatex(0.55,0.75,Form("p_{T,jet} %.0f - %.0f GeV",ptjetmin,ptjetmax));
}

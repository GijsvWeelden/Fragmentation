{
  gROOT->LoadMacro("~/macros/utils.C");
  set_legend_font(62,0.05);

  //Char_t *fname_base="output/pythia_full_jet_spectra_8000GeV_%.0f_%.0f_all.root";
  //Char_t *fname_base="output/pythia_full_jet_spectra_5020GeV_%.0f_%.0f_all.root";
  Char_t *fname_base="output/pythia_charged_jet_spectra_13000GeV_%.0f_%.0f_all.root";
  Char_t *hname = "hJetPtTrigPt_R4";
  Float_t pthard[]={5,15,25,40,100};

  const Int_t do_scaling = 1;

  Int_t nbins = sizeof(pthard)/sizeof(Float_t)-1;
  cout << nbins << endl;
 
  TCanvas *c1= new TCanvas("c1","c1: pt spectra",1000,600);
  c1->Divide(2,1);
  
  Float_t pt_trig = 10;

  TLatex *ltx = new TLatex;
  ltx->SetNDC();

  TH2F *pt_jet_trig_sum = 0;
  TH1D *pt_sum = 0;

  for (Int_t ibin = 0; ibin < nbins; ibin++) {
    TFile *fin = new TFile(Form(fname_base,pthard[ibin],pthard[ibin+1]));
    cout << "pthard " << pthard[ibin] << " - " << pthard[ibin+1]<< " xsec " << hXSec->GetBinContent(1) << " events " << hNEvent->GetBinContent(1) << endl;
    Float_t scale = hXSec->GetBinContent(1)/hNEvent->GetBinContent(1);

    c1->cd(1);
    gPad->SetLogy();
    cout << "make pion projection " << endl;
    hPtPion = hEtaPtPion->ProjectionY(Form("hPtPion_%d",ibin));
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
      hPtPion->DrawCopy("same");
    }

    hptp = hPtPion;

    gPad->Update();
    draw_legend_l(0.6,0.8,hptp,Form("%.0f - %.0f",pthard[ibin],pthard[ibin+1]),ibin);
    cout << "ptpion drawn" << endl;

    c1->cd(2);
    gPad->SetLogy();
    hJetPtTrigPt = (TH2F*)  fin->Get(hname);
    if (do_scaling) {
      hJetPtTrigPt->Sumw2(scale);
      hJetPtTrigPt->Scale(scale);
    }

    px = hJetPtTrigPt->ProjectionX();
    px = new TH1D(*((TH1D*)px));
    px->SetName(Form("px_%d"));
    //px->Sumw2();

    if (ibin == 0) {
      px->DrawCopy();
      pt_sum = new TH1D(*((TH1D*)px));
      pt_sum->SetName("pt_sum");

      pt_jet_trig_sum = new TH2F(*hJetPtTrigPt);
      pt_jet_trig_sum->SetName("pt_jet_trig_sum");
    }
    else {
      px->SetLineColor(1+ibin);
      px->SetMarkerColor(1+ibin);
      px->DrawCopy("same");
      pt_sum->Add(px);
      pt_jet_trig_sum->Add(hJetPtTrigPt);
    }
    cout << "bin done" << endl;
  }

  pt_sum->SetMarkerStyle(25);
  pt_sum->Draw("same");

  TCanvas *c2 = new TCanvas("c2","c2: trigger jets", 500,500);
  gPad->SetLogz();
  pt_jet_trig_sum->Draw("colz");

  TCanvas *c3 = new TCanvas("c3","c3: spectra", 800,500);
  c3->Divide(2,1);
  c3->cd(1);
  gPad->SetLogy();
  pt = pt_jet_trig_sum->ProjectionX("PtAll",0,pt_jet_trig_sum->GetYaxis()->GetNbins()+1);
  pt->Draw();

  Int_t minbin = pt_jet_trig_sum->GetYaxis()->FindBin(pt_trig+0.001);
  trigH = pt_jet_trig_sum->ProjectionX("trigH",minbin, pt_jet_trig_sum->GetYaxis()->GetNbins()+1);
  trigH->SetLineColor(2);
  trigH->Draw("same");

  c3->cd(2);
  
  rat = new TH1D(*trigH);
  rat->SetName("ratio");
  rat->Divide(pt);
  rat->Draw();
}

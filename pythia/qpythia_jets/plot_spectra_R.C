{
  gROOT->LoadMacro("~/macros/utils.C");
  set_legend_font(62,0.05);

  //Char_t *fname_base="output/pythia_full_jets_13000GeV_%.0f_%.0f_all.root";
  //Char_t *fname_base="output/pythia_full_jet_spectra_5020GeV_%.0f_%.0f_all.root";
  Char_t *fname_base="output/pythia_charged_jet_spectra_2760GeV_%.0f_%.0f_all.root";
  Float_t pthard[]={5,15,25,40,100,250};

  TFile *fout = new TFile("qpythia_sum_2TeV76.root","RECREATE");

  const Int_t nR = 3;
  const Float_t Rvals[nR]={0.2,0.3,0.4};

  Int_t nbins = sizeof(pthard)/sizeof(Float_t)-1;
  cout << nbins << endl;
 
  TCanvas *c1= new TCanvas("c1","c1: pt spectra",1000,600);
  c1->Divide(2,1);

  Float_t ptjetmin = 15;
  Float_t ptjetmax = 20;

  TLatex *ltx = new TLatex;
  ltx->SetNDC();

  TH1D *pt_jet_sum[nR] = {0};

  for (Int_t ibin = 0; ibin < nbins; ibin++) {
    TFile *fin = new TFile(Form(fname_base,pthard[ibin],pthard[ibin+1]));
    cout << "pthard " << pthard[ibin] << " - " << pthard[ibin+1]<< " xsec " << hXSec->GetBinContent(1) << " events " << hNEvent->GetBinContent(1) << endl;
    Float_t scale = hXSec->GetBinContent(1)/hNEvent->GetBinContent(1);

    c1->cd(1);
    gPad->SetLogy();
    hPtPion->Sumw2();
    hPtPion->Scale(scale);

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

    for (Int_t iR = 0; iR < nR; iR++) {
      hJetPtEtaPhi = (TH3F*)  fin->Get(Form("hJetPtEtaPhi_R%.0f",10*Rvals[iR]));
      px = hJetPtEtaPhi->Project3D(Form("x_%d"));
      //px = new TH1D(*((TH1D*)px));
      //px->SetName(Form("px_%d"));
      px->Sumw2();
      px->Scale(scale);
      if (ibin == 0) {
	//px->DrawCopy();
	cout << "Clone projection iR " << iR << endl;
	pt_jet_sum[iR] = new TH1D(*((TH1D*)px));
	pt_jet_sum[iR]->SetName(Form("pt_jet_sum_R%.0f",10*Rvals[iR]));
      }
      else {
	//px->SetLineColor(1+ibin);
	//px->SetMarkerColor(1+ibin);
	//px->DrawCopy("same");
	cout << "Add iR " << iR << endl;
	pt_jet_sum[iR]->Add(px);
      }
    }
  }

  c1->cd(2);
  for (Int_t iR = 0; iR < nR; iR++) {
    pt_jet_sum[iR]->SetLineWidth(2);
    pt_jet_sum[iR]->SetMarkerColor(4);
    pt_jet_sum[iR]->SetMarkerStyle(20+iR);
    if (iR == 0) 
      pt_jet_sum[iR]->Draw();
    else
      pt_jet_sum[iR]->Draw("same");
    if (fout) {
      fout->cd();
      pt_jet_sum[iR]->Write();
    }
  }
  
  TCanvas *c3 = new TCanvas("c3","c3: ratios",500,500);
  ptrat_02_04 = new TH1D(*(pt_jet_sum[0]));
  ptrat_02_04->SetTitle("ptrat_02_04");
  ptrat_02_04->Divide(pt_jet_sum[2]);
  ptrat_02_04->Draw();

}

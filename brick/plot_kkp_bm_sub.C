{
  // NB this does not run standalone, called by plot_kkp_bm.C
  gROOT->SetStyle("Plain");
  gStyle->SetNdivisions(505);
  gStyle->SetNdivisions(505,"Y");
  gStyle->SetMarkerSize(1.2);

  const Float_t E = 50;

  TCanvas *c2 = new TCanvas("c2","c2: FF",500,500);
  gPad->SetLogy();
  h1 = gPad->DrawFrame(0,0.01,1,1000);
  h1->SetXTitle("z");
  h1->SetYTitle("D(z)");

  TF1 *frag_f = new TF1("frag_f",kkp_func,0,1,3);
  frag_f->SetParameter(0,E); // Q for fragmentation
  frag_f->SetParameter(1,1); // parton flavour 1=u
  frag_f->SetParameter(2,1); // hadron type: 1=pi+pi-; 4 = p+pbar

  frag_f->Draw("same");

  TF1 *frag_f_prot = new TF1("frag_f_prot",kkp_func,0,1,3);
  frag_f_prot->SetParameter(0,E); // Q for fragmentation
  frag_f_prot->SetParameter(1,1); // parton flavour 1=u
  frag_f_prot->SetParameter(2,4); // hadron type: 1=pi+pi-; 4 = p+pbar
  frag_f_prot->SetLineStyle(2);
  frag_f_prot->Draw("same");

  TF1 *frag_f_glue = new TF1("frag_f_glue",kkp_func,0,1,3);
  frag_f_glue->SetParameter(0,E); // Q for fragmentation
  frag_f_glue->SetParameter(1,0); // parton flavour 0=g
  frag_f_glue->SetParameter(2,1); // hadron type: 1=pi+pi-; 4 = p+pbar
  frag_f_glue->SetLineColor(2);
  frag_f_glue->Draw("same");

  TF1 *frag_f_glue_prot = new TF1("frag_f_glue_prot",kkp_func,0,1,3);
  frag_f_glue_prot->SetParameter(0,E); // Q for fragmentation
  frag_f_glue_prot->SetParameter(1,0); // parton flavour 0=g
  frag_f_glue_prot->SetParameter(2,4); // hadron type: 1=pi+pi-; 4 = p+pbar
  frag_f_glue_prot->SetLineColor(kGreen+1);
  frag_f_glue_prot->SetLineStyle(2);
  frag_f_glue_prot->Draw("same");

  TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->SetBorderSize(0);
  leg->AddEntry(frag_f,"u #rightarrow #pi","l");
  leg->AddEntry(frag_f_prot,"u #rightarrow p","l");
  leg->AddEntry(frag_f_glue,"g #rightarrow #pi","l");
  leg->AddEntry(frag_f_glue_prot,"g #rightarrow p","l");
  leg->Draw();
}

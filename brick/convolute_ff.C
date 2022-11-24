{
  gROOT->SetStyle("Plain");
  gROOT->LoadMacro("brick_utils.C");

  TGraph *whdg_rad_gr = 0;
  TGraph *whdg_edge_gr = 0;
  TGraph *whdg_tot_gr = 0;
  read_whdg(whdg_rad_gr, whdg_edge_gr, whdg_tot_gr);

  TF1 *frag_f = new TF1("frag_f","exp([0]*x)",0,1);
  frag_f->SetParameter(0,-8);

  normalize(whdg_rad_gr, whdg_edge_gr);
  whdg_ff_gr = convolute_ff(10, whdg_rad_gr, whdg_edge_gr, frag_f);

  whdg_ff_gr->Draw("la");
  gPad->SetLogy();

  h1 = whdg_ff_gr->GetHistogram();
  h1->SetXTitle("z");
  h1->SetYTitle("dN/dz");

  frag_f->SetLineWidth(1);
  frag_f->SetLineStyle(2);
  frag_f->Draw("same");

  TGraph *asw_rad_gr = 0;
  TGraph *asw_edge_gr = 0;
  read_asw(asw_rad_gr, asw_edge_gr);
  normalize(asw_rad_gr, asw_edge_gr);
  asw_ff_gr = convolute_ff(10, asw_rad_gr, asw_edge_gr, frag_f);
  asw_ff_gr->SetLineColor(2);
  asw_ff_gr->Draw("l");  

  /*
  asw_ff_h = convolute_ff_mc(10, asw_rad_gr, asw_edge_gr, frag_f);
  asw_ff_h->SetLineColor(2);
  asw_ff_h->Draw("same");  
  */
}

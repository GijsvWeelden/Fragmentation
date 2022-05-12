{
  gROOT->SetStyle("Plain");
  gStyle->SetNdivisions(505);
  gStyle->SetNdivisions(505,"Y");
  gStyle->SetMarkerSize(1.2);
  gROOT->LoadMacro("brick_utils.C");
  gSystem->Load("kkp/libkkp.so");
  gROOT->LoadMacro("kkp_func.C");
  gROOT->LoadMacro("~/macros/utils.C");

  const Float_t E = 10;

  TCanvas *c1 = new TCanvas("c1","c1: P(DE)",500,500);
  TCanvas *c2 = new TCanvas("c2","c2: FF ratio",500,500);
  h1 = gPad->DrawFrame(0,0,1,1);
  h1->SetXTitle("z");
  h1->SetYTitle("D(z)/pp");
  
  TF1 *frag_f = new TF1("frag_f",kkp_func,0,1,2);
  frag_f->SetParameter(0,E); // Q for fragmentation
  frag_f->SetParameter(1,1); // parton flavour 1=u

  frag_f->SetLineWidth(1);
  frag_f->SetLineStyle(2);
  //frag_f->Draw("same");

  TF1 *frag_f_glue = new TF1("frag_f_glue",kkp_func,0,1,2);
  frag_f_glue->SetParameter(0,E); // Q for fragmentation
  frag_f_glue->SetParameter(1,0); // parton flavour 0=g

  c1->cd();
  h1 = gPad->DrawFrame(0,0,1,1);
  h1->SetXTitle("#DeltaE/E");
  h1->SetYTitle("P(#DeltaE/E)");

  TLatex *ltx = new TLatex();
  ltx->SetNDC();
  ltx->DrawLatex(ndc_x(0.4),ndc_y(0.85),"R_{8} #approx 0.2");
  bw_cont_gr = new TGraph();
  bw_edge_gr = new TGraph();
  bw_edge_gr->SetPoint(0,0,0.2);
  bw_edge_gr->SetPoint(1,1,0.8);
  for (Int_t i = 0; i < 10; i++) {
    bw_cont_gr->SetPoint(i,(i+0.5)*0.1,0);
  }
  bw_edge_gr->SetMarkerStyle(20);
  bw_edge_gr->Draw("p");
  bw_cont_gr->SetLineWidth(2);
  bw_cont_gr->Draw("l");
  Float_t bw_R8 = calc_R(bw_cont_gr,bw_edge_gr);
  cout << "R8 B&W " << bw_R8 << endl;
  draw_legend_l(0.4,0.8,bw_cont_gr,"Black&White",0,2,0.04,62);

  bw2_cont_gr = new TGraph();
  bw2_edge_gr = new TGraph();
  bw2_edge_gr->SetPoint(0,0,0.19);
  bw2_edge_gr->SetPoint(1,1,0.61);
  for (Int_t i = 0; i < 10; i++) {
    bw2_cont_gr->SetPoint(i,(i+0.5)*0.1,0.2);
  }
  bw2_edge_gr->SetMarkerStyle(20);
  bw2_edge_gr->SetMarkerColor(2);
  bw2_edge_gr->Draw("p");
  bw2_cont_gr->SetLineWidth(2);
  bw2_cont_gr->SetLineColor(2);
  bw2_cont_gr->Draw("l");
  Float_t bw2_R8 = calc_R(bw2_cont_gr,bw2_edge_gr);
  cout << "R8 B&W+cont " << bw2_R8 << endl;
  draw_legend_l(0.4,0.8,bw2_cont_gr,"B&W+constant",1,2,0.04,62);

  g_cont_gr = new TGraph();
  g_edge_gr = new TGraph();
  //Float_t p0 = 0.17;
  Float_t p0 = 0.0;
  g_edge_gr->SetPoint(0,0,p0);
  g_edge_gr->SetPoint(1,1,0);
  Float_t g_mean = 0.08;
  Float_t g_sig = 0.03;
  Float_t A = (1.-p0)/sqrt(2*TMath::Pi())/g_sig;
  for (Int_t i = 0; i < 50; i++) {
    Float_t x = (i+0.5)*0.02;
    g_cont_gr->SetPoint(i,x,A*exp(-0.5*(x-g_mean)*(x-g_mean)/g_sig/g_sig));
  }
  g_edge_gr->SetMarkerStyle(20);
  g_edge_gr->SetMarkerColor(4);
  g_edge_gr->Draw("p");
  g_cont_gr->SetLineWidth(2);
  g_cont_gr->SetLineColor(4);
  g_cont_gr->Draw("l");
  Float_t g_R8 = calc_R(g_cont_gr,g_edge_gr);
  cout << "R8 gaus " << g_R8 << endl;
  draw_legend_l(0.4,0.8,g_cont_gr,"Gaus",2,2,0.04,62);

  c2->cd();
  bw_ff_gr = convolute_ff(E, bw_cont_gr, bw_edge_gr, frag_f,49,0.02,0.98);
  bw_rat_gr = divide_graph(bw_ff_gr,frag_f);
  bw_rat_gr->SetLineWidth(2);
  bw_rat_gr->Draw("l");

  bw2_ff_gr = convolute_ff(E, bw2_cont_gr, bw2_edge_gr, frag_f,49,0.02,0.98);
  bw2_rat_gr = divide_graph(bw2_ff_gr,frag_f);
  bw2_rat_gr->SetLineWidth(2);
  bw2_rat_gr->SetLineColor(2);
  bw2_rat_gr->Draw("l");

  g_ff_gr = convolute_ff(E, g_cont_gr, g_edge_gr, frag_f,49,0.02,0.98);
  g_rat_gr = divide_graph(g_ff_gr,frag_f);
  g_rat_gr->SetLineWidth(2);
  g_rat_gr->SetLineColor(4);
  g_rat_gr->Draw("l");

}


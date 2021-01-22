#include "style.C"
#include "legend_utils.C"
#include "jewel_plot_utils.C"
#include "TLegend.h"

void plot_jet_variables(void) {
  // pt range
  Float_t min_pt = 80;
  Float_t max_pt = 120;
  Float_t jetR = 0.4;//2;//5;

  // Which observables to plot?
  bool plot_mass = true;
  bool plot_Rg = false;
  bool plot_rz = false; // Not yet implemented
  bool plot_zg = false;

  // TChain AA
  int numFiles = 20;
  TChain *chain = new TChain("jetprops");
  for (int fileNum = 1; fileNum <= numFiles; fileNum++) {
    chain->AddFile(Form("../run_AA_2tev76_norecoil/jet_shapes_constsub_eventwise_tree_%d_full_nobkg.root", fileNum));
  }

  // TODO: how many bins?
  // Initialise desired histograms
  TH1F *hmass_AA = new TH1F("hmass_AA", Form("PbPb2tev76 mass, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,40.0);
  TH1F *hmass_pp = new TH1F("hmass_pp", Form("pp2tev76 mass, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,40.0);
  TH1F *hRg_AA = new TH1F("hRg_AA", Form("PbPb2tev76 Rg, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,0.5);
  TH1F *hRg_pp = new TH1F("hRg_pp", Form("pp2tev76 Rg, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,0.5);
  TH1F *hrz_AA = new TH1F("hrz_AA", Form("PbPb2tev76 rz, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,30.0);
  TH1F *hrz_pp = new TH1F("hrz_pp", Form("pp2tev76 rz, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,30.0);
  TH1F *hzg_AA = new TH1F("hzg_AA", Form("PbPb2tev76 zg, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,0.5);
  TH1F *hzg_pp = new TH1F("hzg_pp", Form("pp2tev76 zg, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,0.5);

  /* Declare leaves, branches and addresses */
  // Declaration of leaf types
  Int_t           ievt;
  Int_t           ijet;
  Float_t         evwt;
  Float_t         pt;
  Float_t         eta;
  Float_t         phi;
  Float_t         dphi;
  Int_t           nconst;
  Float_t         zg;
  Float_t         Rg;
  Int_t           nSD;
  Float_t         mass;
  Float_t         mz2;
  Float_t         mr;
  Float_t         mr2;
  Float_t         rz;
  Float_t         r2z;
  // List of branches
  TBranch        *b_ievt;   //!
  TBranch        *b_ijet;   //!
  TBranch        *b_evwt;   //!
  TBranch        *b_pt;   //!
  TBranch        *b_eta;   //!
  TBranch        *b_phi;   //!
  TBranch        *b_dphi;   //!
  TBranch        *b_nconst;   //!
  TBranch        *b_zg;   //!
  TBranch        *b_Rg;   //!
  TBranch        *b_nSD;   //!
  TBranch        *b_mass;   //!
  TBranch        *b_mz2;   //!
  TBranch        *b_mr;   //!
  TBranch        *b_mr2;   //!
  TBranch        *b_rz;   //!
  TBranch        *b_r2z;   //!
  // Set branch addresses
  chain->SetBranchAddress("ievt", &ievt, &b_ievt);
  chain->SetBranchAddress("ijet", &ijet, &b_ijet);
  chain->SetBranchAddress("evwt", &evwt, &b_evwt);
  chain->SetBranchAddress("pt", &pt, &b_pt);
  chain->SetBranchAddress("eta", &eta, &b_eta);
  chain->SetBranchAddress("phi", &phi, &b_phi);
  chain->SetBranchAddress("dphi", &dphi, &b_dphi);
  chain->SetBranchAddress("nconst", &nconst, &b_nconst);
  chain->SetBranchAddress("zg", &zg, &b_zg);
  chain->SetBranchAddress("Rg", &Rg, &b_Rg);
  chain->SetBranchAddress("nSD", &nSD, &b_nSD);
  chain->SetBranchAddress("mass", &mass, &b_mass);
  chain->SetBranchAddress("mz2", &mz2, &b_mz2);
  chain->SetBranchAddress("mr", &mr, &b_mr);
  chain->SetBranchAddress("mr2", &mr2, &b_mr2);
  chain->SetBranchAddress("rz", &rz, &b_rz);
  chain->SetBranchAddress("r2z", &r2z, &b_r2z);

  Long64_t nentries = chain->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  /* Analysis here */
  for (Long64_t entry = 0; entry < nentries; entry++){
    Long64_t local_entry = chain->LoadTree(entry);
    if (local_entry < 0) break;
    nb = chain->GetEntry(entry);
    nbytes+=nb;
    if (pt >= min_pt && pt < max_pt){
      if (plot_mass) hmass_AA->Fill(mass,evwt);
      if (plot_Rg) hRg_AA->Fill(Rg,evwt);
      if (plot_zg) hzg_AA->Fill(zg,evwt);
    }
  }

  // TChain pp
  int numFiles_pp = 1;
  chain->Reset();
  chain->AddFile("../run_pp_2tev76/jet_tree_pp2tev76_full_nobkg.root");
  /*
     for (int fileNum = 1; fileNum <= numFiles_pp; fileNum++) {
     chain->AddFile(Form("../run_pp_2tev76/jet_shapes_constsub_eventwise_tree_%d_full_nobkg.root", fileNum));
     }
   */

  // Reset branch addresses
  chain->SetBranchAddress("ievt", &ievt, &b_ievt);
  chain->SetBranchAddress("ijet", &ijet, &b_ijet);
  chain->SetBranchAddress("evwt", &evwt, &b_evwt);
  chain->SetBranchAddress("pt", &pt, &b_pt);
  chain->SetBranchAddress("eta", &eta, &b_eta);
  chain->SetBranchAddress("phi", &phi, &b_phi);
  chain->SetBranchAddress("dphi", &dphi, &b_dphi);
  chain->SetBranchAddress("nconst", &nconst, &b_nconst);
  chain->SetBranchAddress("zg", &zg, &b_zg);
  chain->SetBranchAddress("Rg", &Rg, &b_Rg);
  chain->SetBranchAddress("nSD", &nSD, &b_nSD);
  chain->SetBranchAddress("mass", &mass, &b_mass);
  chain->SetBranchAddress("mz2", &mz2, &b_mz2);
  chain->SetBranchAddress("mr", &mr, &b_mr);
  chain->SetBranchAddress("mr2", &mr2, &b_mr2);
  chain->SetBranchAddress("rz", &rz, &b_rz);
  chain->SetBranchAddress("r2z", &r2z, &b_r2z);

  nentries = chain->GetEntries();

  for (Long64_t entry = 0; entry < nentries; entry++){
    Long64_t local_entry = chain->LoadTree(entry);
    if (local_entry < 0) break;
    nb = chain->GetEntry(entry);
    nbytes+=nb;
    if (pt >= min_pt && pt < max_pt){
      if (plot_mass) hmass_pp->Fill(mass,evwt);
      if (plot_Rg) hRg_pp->Fill(Rg,evwt);
      if (plot_zg) hzg_pp->Fill(zg,evwt);
    }
  }

  if (plot_mass){
    cout << "Plotting mass" << endl;
    // Normalise histograms
    hmass_pp->Scale(1./hmass_pp->Integral());
    hmass_AA->Scale(1./hmass_AA->Integral());
    // Make canvases and draw histograms
    TCanvas *c1 = new TCanvas("c1","c1: mass",600,600);
    //c1->DrawFrame(-10,-10,10,10);
    //c1->SetLeftMargin(0.5);
    //gStyle->SetPadLeftMargin(0.05);
    //gStyle->SetOptStat(0);
    hmass_pp->SetMarkerStyle(kFullDotMedium);
    hmass_pp->SetMarkerColor(kRed);
    hmass_pp->SetLineColor(kRed);
    hmass_pp->SetTitle("");
    //hmass_AA->SetTitle("");
    hmass_AA->GetXaxis()->SetTitle("mass");
    hmass_AA->GetYaxis()->SetTitle("fraction");
    hmass_AA->Draw("E1");
    hmass_pp->Draw("E1,SAMES");
    //c1->BuildLegend(0.5,0.9,0.5,0.9);
    // Make and draw legend
    //raw_legend_m(0.4,0.86,hmass_AA,"AA",1);
    /*TLegend leg(.1,.7,.3,.9,"Lab. Lesson 1");
      leg.SetFillColor(0);
      leg.AddEntry(hmass_AA,"AA");
      leg.AddEntry(hmass_pp,"pp");
      leg.DrawClone("Same");*/
    //gStyle->SetOptStat();
    //c1->BuildLegend();
    //legend->Draw("SAMES");
    //c1->Update();
    // Compute ratio of histograms and draw ratio plot
    auto massratio = new TRatioPlot(hmass_AA, hmass_pp);
    massratio->Draw();
    massratio->GetLowerRefYaxis()->SetRangeUser(-0.2,6.);
    massratio->GetLowerRefYaxis()->SetTitle("PbPb/pp");
    //c1->BuildLegend();
    //c1->Update();
  }
  if (plot_Rg){
    cout << "Plotting Rg" << endl;
    // Normalise histograms
    hRg_pp->Scale(1./hRg_pp->Integral());
    hRg_AA->Scale(1./hRg_AA->Integral());
    // Make canvases and draw histograms
    TCanvas *c2 = new TCanvas("c2","c2: Rg",600,600);
    hRg_pp->SetMarkerStyle(kFullDotMedium);
    hRg_pp->SetMarkerColor(kRed);
    hRg_pp->SetLineColor(kRed);
    hRg_AA->GetXaxis()->SetTitle("Rg");
    hRg_AA->GetYaxis()->SetTitle("fraction");
    hRg_AA->GetXaxis()->SetRangeUser(0.1,0.6);
    hRg_AA->GetYaxis()->SetRangeUser(0.,0.025);
    hRg_AA->Draw("E1");
    hRg_pp->Draw("E1,SAMES");
    // Compute ratio of histograms and draw ratio plot
    auto Rg_ratio = new TRatioPlot(hRg_AA, hRg_pp);
    Rg_ratio->Draw();
    Rg_ratio->GetLowerRefYaxis()->SetRangeUser(-0.2,2.5);
    Rg_ratio->GetLowerRefYaxis()->SetTitle("PbPb/pp");
    c2->Update();
  }
  if (plot_zg){
    cout << "Plotting zg" << endl;
    // Normalise histograms
    hzg_pp->Scale(1./hzg_pp->Integral());
    hzg_AA->Scale(1./hzg_AA->Integral());
    // Make canvases and draw histograms
    TCanvas *c3 = new TCanvas("c3","c3: zg",600,600);
    hzg_pp->SetMarkerStyle(kFullDotMedium);
    hzg_pp->SetMarkerColor(kRed);
    hzg_pp->SetLineColor(kRed);
    hzg_AA->GetXaxis()->SetTitle("zg");
    hzg_AA->GetYaxis()->SetTitle("fraction");
    hzg_AA->GetXaxis()->SetRangeUser(0.1,0.6);
    hzg_AA->Draw("E1");
    hzg_pp->Draw("E1,SAMES");
    // Compute ratio of histograms and draw ratio plot
    auto zg_ratio = new TRatioPlot(hzg_AA, hzg_pp);
    zg_ratio->Draw();
    zg_ratio->GetLowerRefYaxis()->SetRangeUser(-0.2,2.);
    zg_ratio->GetLowerRefYaxis()->SetTitle("PbPb/pp");
    c3->Update();
  }

}

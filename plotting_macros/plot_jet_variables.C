#include "style.C"
#include "legend_utils.C"
#include "jewel_plot_utils.C"
#include "TLegend.h"
#include <typeinfo>

void plot_jet_variables(void) {
  //gInterpreter->LoadMacro("./plot.C"<PathToYourMacro>)

  // pt range
  // pt_bins = [0,20,40,60,80,100,120,160,200];
  Float_t min_pt = 60;
  Float_t max_pt = 80;
  Float_t jetR = 0.4;

  // Which observables to plot?
  bool plot_nconst = true;
  bool plot_zg = true;
  bool plot_Rg = true;
  bool plot_nSD = true;
  bool plot_mass = true;
  bool plot_mz2 = true;
  bool plot_mr = true;
  bool plot_mr2 = true;
  bool plot_rz = true;
  bool plot_r2z = true;

  // TODO: Can we auto-select the bins?
  // TODO: Set right limits
  // Initialise desired histograms
  TH1F *hnconst_AA = new TH1F("hnconst_AA", Form("PbPb2tev76 nconst, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,30.0);
  TH1F *hnconst_pp = new TH1F("hnconst_pp", Form("pp2tev76 nconst, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,30.0);
  TH1F *hzg_AA = new TH1F("hzg_AA", Form("PbPb2tev76 zg, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,0.5);
  TH1F *hzg_pp = new TH1F("hzg_pp", Form("pp2tev76 zg, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,0.5);
  TH1F *hRg_AA = new TH1F("hRg_AA", Form("PbPb2tev76 Rg, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,0.5);
  TH1F *hRg_pp = new TH1F("hRg_pp", Form("pp2tev76 Rg, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,0.5);
  TH1F *hnSD_AA = new TH1F("hnSD_AA", Form("PbPb2tev76 nSD, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,30.0);
  TH1F *hnSD_pp = new TH1F("hnSD_pp", Form("pp2tev76 nSD, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,30.0);
  TH1F *hmass_AA = new TH1F("hmass_AA", Form("PbPb2tev76 mass, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,40.);
  TH1F *hmass_pp = new TH1F("hmass_pp", Form("pp2tev76 mass, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,40.);
  TH1F *hmz2_AA = new TH1F("hmz2_AA", Form("PbPb2tev76 mz2, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,0.5);
  TH1F *hmz2_pp = new TH1F("hmz2_pp", Form("pp2tev76 mz2, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,0.5);
  TH1F *hmr_AA = new TH1F("hmr_AA", Form("PbPb2tev76 mr, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,0.5);
  TH1F *hmr_pp = new TH1F("hmr_pp", Form("pp2tev76 mr, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,0.5);
  TH1F *hmr2_AA = new TH1F("hmr2_AA", Form("PbPb2tev76 mr2, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,0.5);
  TH1F *hmr2_pp = new TH1F("hmr2_pp", Form("pp2tev76 mr2, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,0.5);
  TH1F *hrz_AA = new TH1F("hrz_AA", Form("PbPb2tev76 rz, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,30.0);
  TH1F *hrz_pp = new TH1F("hrz_pp", Form("pp2tev76 rz, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,30.0);
  TH1F *hr2z_AA = new TH1F("hr2z_AA", Form("PbPb2tev76 r2z, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,0.5);
  TH1F *hr2z_pp = new TH1F("hr2z_pp", Form("pp2tev76 r2z, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,0.5);

  //-------------------------------------------------------------
  //
  // AA analysis
  //
  //-------------------------------------------------------------

  // TChain AA
  int numFiles = 20;
  TChain *chain = new TChain("jetprops");
  for (int fileNum = 1; fileNum <= numFiles; fileNum++) {
    chain->AddFile(Form("../run_AA_2tev76_norecoil/jet_shapes_constsub_eventwise_tree_%d_full_nobkg.root", fileNum));
  }

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

  // Analysis here
  for (Long64_t entry = 0; entry < nentries; entry++){
    Long64_t local_entry = chain->LoadTree(entry);
    if (local_entry < 0) break;
    nb = chain->GetEntry(entry);
    nbytes+=nb;
    if (pt >= min_pt && pt < max_pt){
      // TODO: Add histogram filling
      if (plot_nconst) hnconst_AA->Fill(nconst,evwt);
      if (plot_zg) hzg_AA->Fill(zg,evwt);
      if (plot_Rg) hRg_AA->Fill(Rg,evwt);
      if (plot_nSD) hnSD_AA->Fill(nSD,evwt);
      if (plot_mass) hmass_AA->Fill(mass,evwt);
      if (plot_mz2) hmz2_AA->Fill(mz2,evwt);
      if (plot_mr) hmr_AA->Fill(mr,evwt);
      if (plot_mr2) hmr2_AA->Fill(mr2,evwt);
      if (plot_rz) hrz_AA->Fill(rz,evwt);
      if (plot_r2z) hr2z_AA->Fill(r2z,evwt);
    }
  }

  //-------------------------------------------------------------
  //
  // pp analysis
  //
  //-------------------------------------------------------------

  // TChain pp
  int numFiles_pp = 1;
  chain->Reset();
  chain->AddFile("../run_pp_2tev76/jet_tree_pp2tev76_full_nobkg.root");
  /* for (int fileNum = 1; fileNum <= numFiles_pp; fileNum++) {
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
      // TODO: add more histograms
      if (plot_nconst) hnconst_pp->Fill(nconst,evwt);
      if (plot_zg) hzg_pp->Fill(zg,evwt);
      if (plot_Rg) hRg_pp->Fill(Rg,evwt);
      if (plot_nSD) hnSD_pp->Fill(nSD,evwt);
      if (plot_mass) hmass_pp->Fill(mass,evwt);
      if (plot_mz2) hmz2_pp->Fill(mz2,evwt);
      if (plot_mr) hmr_pp->Fill(mr,evwt);
      if (plot_mr2) hmr2_pp->Fill(mr2,evwt);
      if (plot_rz) hrz_pp->Fill(rz,evwt);
      if (plot_r2z) hr2z_pp->Fill(r2z,evwt);
      /* 
         if (plot_mass) hmass_pp->Fill(mass,evwt);
         if (plot_Rg) hRg_pp->Fill(Rg,evwt);
         if (plot_zg) hzg_pp->Fill(zg,evwt);
       */  }
  }

  //-------------------------------------------------------------
  //
  // Plotting variables
  //
  //-------------------------------------------------------------

  if (plot_nconst){
    cout << "Plotting nconst" << endl;
    hnconst_pp->Scale(1./hnconst_pp->Integral());
    hnconst_AA->Scale(1./hnconst_AA->Integral());
    TH1F *ratio = (TH1F*)hnconst_AA->Clone("Ratio");

    // Check normalisation within precision
    if (hnconst_pp->Integral() -  1.0 > 0.01 || hnconst_AA->Integral() - 1.0 > 0.01){
      cout << "WARNING: NORMALISATION PROBLEM!" << endl << "pp = " << hnconst_pp->Integral() << " ( " << hnconst_pp->Integral()-1.0  << ")" << endl 
        << "AA = " << hnconst_AA->Integral() << " (type: " << hnconst_AA->Integral()-1.0 << ")" << endl;
    }

    // Plot settings
    hnconst_AA->SetStats(0);
    hnconst_AA->SetMarkerStyle(kPlus);
    hnconst_AA->SetTitle(Form("Jet nconst, pt #in [%.0f,%.0f), R = %.1f", min_pt, max_pt, jetR));
    // y axis
    hnconst_AA->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN}{dnconst}");
    hnconst_AA->GetYaxis()->SetTitleFont(43);
    hnconst_AA->GetYaxis()->SetTitleSize(20);//0.05);
    hnconst_AA->GetYaxis()->SetTitleOffset(2.5);
    hnconst_AA->GetYaxis()->SetLabelFont(43);
    hnconst_AA->GetYaxis()->SetLabelSize(14);

    hnconst_pp->SetStats(0);
    hnconst_pp->SetMarkerStyle(kMultiply);
    hnconst_pp->SetLineColor(kRed);
    hnconst_pp->SetMarkerColor(kRed);

    ratio->Divide(hnconst_pp);
    ratio->SetStats(0);
    ratio->SetMarkerStyle(kDot);
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerColor(kBlack);
    ratio->SetTitle("");
    // x axis
    ratio->GetXaxis()->SetTitle("nconst");
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleOffset(2.5);
    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelSize(14);
    // y axis
    ratio->GetYaxis()->SetTitle("PbPb/pp");
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleOffset(2.5);
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(14);

    // Draw histograms
    TCanvas *c_nconst = new TCanvas("c_nconst","nconst",600,900);
    TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
    pad1->SetBottomMargin(0);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();
    hnconst_AA->Draw();
    hnconst_pp->Draw("same");
    auto legend = new TLegend(0.8,0.8,0.9,0.9);
    legend->AddEntry(hnconst_AA,"AA");
    legend->AddEntry(hnconst_pp,"pp");
    legend->Draw();
    c_nconst->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();
    ratio->Draw();
    c_nconst->cd();
  } //-------------------------------------------------------------
  if (plot_zg){
    cout << "Plotting zg" << endl;
    hzg_pp->Scale(1./hzg_pp->Integral());
    hzg_AA->Scale(1./hzg_AA->Integral());
    TH1F *ratio = (TH1F*)hzg_AA->Clone("Ratio");

    // Check normalisation within precision
    if (hzg_pp->Integral() -  1.0 > 0.01 || hzg_AA->Integral() - 1.0 > 0.01){
      cout << "WARNING: NORMALISATION PROBLEM!" << endl << "pp = " << hzg_pp->Integral() << " ( " << hzg_pp->Integral()-1.0  << ")" << endl 
        << "AA = " << hzg_AA->Integral() << " (type: " << hzg_AA->Integral()-1.0 << ")" << endl;
    }

    // Plot settings
    hzg_AA->SetStats(0);
    hzg_AA->SetMarkerStyle(kPlus);
    hzg_AA->SetTitle(Form("Jet zg, pt #in [%.0f,%.0f), R = %.1f", min_pt, max_pt, jetR));
    // y axis
    hzg_AA->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN}{dzg}");
    hzg_AA->GetYaxis()->SetTitleFont(43);
    hzg_AA->GetYaxis()->SetTitleSize(20);//0.05);
    hzg_AA->GetYaxis()->SetTitleOffset(2.5);
    hzg_AA->GetYaxis()->SetLabelFont(43);
    hzg_AA->GetYaxis()->SetLabelSize(14);

    hzg_pp->SetStats(0);
    hzg_pp->SetMarkerStyle(kMultiply);
    hzg_pp->SetLineColor(kRed);
    hzg_pp->SetMarkerColor(kRed);

    ratio->Divide(hzg_pp);
    ratio->SetStats(0);
    ratio->SetMarkerStyle(kDot);
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerColor(kBlack);
    ratio->SetTitle("");
    // x axis
    ratio->GetXaxis()->SetTitle("zg");
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleOffset(2.5);
    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelSize(14);
    // y axis
    ratio->GetYaxis()->SetTitle("PbPb/pp");
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleOffset(2.5);
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(14);

    // Draw histograms
    TCanvas *c_zg = new TCanvas("c_zg","zg",600,900);
    TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
    pad1->SetBottomMargin(0);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();
    hzg_AA->Draw();
    hzg_pp->Draw("same");
    auto legend = new TLegend(0.8,0.8,0.9,0.9);
    legend->AddEntry(hzg_AA,"AA");
    legend->AddEntry(hzg_pp,"pp");
    legend->Draw();
    c_zg->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();
    ratio->Draw();
    c_zg->cd();
  } //-------------------------------------------------------------
  if (plot_Rg){
    cout << "Plotting Rg" << endl;
    hRg_pp->Scale(1./hRg_pp->Integral());
    hRg_AA->Scale(1./hRg_AA->Integral());
    TH1F *ratio = (TH1F*)hRg_AA->Clone("Ratio");

    // Check normalisation within precision
    if (hRg_pp->Integral() -  1.0 > 0.01 || hRg_AA->Integral() - 1.0 > 0.01){
      cout << "WARNING: NORMALISATION PROBLEM!" << endl << "pp = " << hRg_pp->Integral() << " ( " << hRg_pp->Integral()-1.0  << ")" << endl 
        << "AA = " << hRg_AA->Integral() << " (type: " << hRg_AA->Integral()-1.0 << ")" << endl;
    }

    // Plot settings
    hRg_AA->SetStats(0);
    hRg_AA->SetMarkerStyle(kPlus);
    hRg_AA->SetTitle(Form("Jet Rg, pt #in [%.0f,%.0f), R = %.1f", min_pt, max_pt, jetR));
    // y axis
    hRg_AA->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN}{dRg}");
    hRg_AA->GetYaxis()->SetTitleFont(43);
    hRg_AA->GetYaxis()->SetTitleSize(20);//0.05);
    hRg_AA->GetYaxis()->SetTitleOffset(2.5);
    hRg_AA->GetYaxis()->SetLabelFont(43);
    hRg_AA->GetYaxis()->SetLabelSize(14);

    hRg_pp->SetStats(0);
    hRg_pp->SetMarkerStyle(kMultiply);
    hRg_pp->SetLineColor(kRed);
    hRg_pp->SetMarkerColor(kRed);

    ratio->Divide(hRg_pp);
    ratio->SetStats(0);
    ratio->SetMarkerStyle(kDot);
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerColor(kBlack);
    ratio->SetTitle("");
    // x axis
    ratio->GetXaxis()->SetTitle("Rg");
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleOffset(2.5);
    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelSize(14);
    // y axis
    ratio->GetYaxis()->SetTitle("PbPb/pp");
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleOffset(2.5);
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(14);

    // Draw histograms
    TCanvas *c_Rg = new TCanvas("c_Rg","Rg",600,900);
    TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
    pad1->SetBottomMargin(0);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();
    hRg_AA->Draw();
    hRg_pp->Draw("same");
    auto legend = new TLegend(0.8,0.8,0.9,0.9);
    legend->AddEntry(hRg_AA,"AA");
    legend->AddEntry(hRg_pp,"pp");
    legend->Draw();
    c_Rg->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();
    ratio->Draw();
    c_Rg->cd();
  } //-------------------------------------------------------------
  if (plot_nSD){
    cout << "Plotting nSD" << endl;
    hnSD_pp->Scale(1./hnSD_pp->Integral());
    hnSD_AA->Scale(1./hnSD_AA->Integral());
    TH1F *ratio = (TH1F*)hnSD_AA->Clone("Ratio");

    // Check normalisation within precision
    if (hnSD_pp->Integral() -  1.0 > 0.01 || hnSD_AA->Integral() - 1.0 > 0.01){
      cout << "WARNING: NORMALISATION PROBLEM!" << endl << "pp = " << hnSD_pp->Integral() << " ( " << hnSD_pp->Integral()-1.0  << ")" << endl 
        << "AA = " << hnSD_AA->Integral() << " (type: " << hnSD_AA->Integral()-1.0 << ")" << endl;
    }

    // Plot settings
    hnSD_AA->SetStats(0);
    hnSD_AA->SetMarkerStyle(kPlus);
    hnSD_AA->SetTitle(Form("Jet nSD, pt #in [%.0f,%.0f), R = %.1f", min_pt, max_pt, jetR));
    // y axis
    hnSD_AA->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN}{dnSD}");
    hnSD_AA->GetYaxis()->SetTitleFont(43);
    hnSD_AA->GetYaxis()->SetTitleSize(20);//0.05);
    hnSD_AA->GetYaxis()->SetTitleOffset(2.5);
    hnSD_AA->GetYaxis()->SetLabelFont(43);
    hnSD_AA->GetYaxis()->SetLabelSize(14);

    hnSD_pp->SetStats(0);
    hnSD_pp->SetMarkerStyle(kMultiply);
    hnSD_pp->SetLineColor(kRed);
    hnSD_pp->SetMarkerColor(kRed);

    ratio->Divide(hnSD_pp);
    ratio->SetStats(0);
    ratio->SetMarkerStyle(kDot);
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerColor(kBlack);
    ratio->SetTitle("");
    // x axis
    ratio->GetXaxis()->SetTitle("nSD");
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleOffset(2.5);
    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelSize(14);
    // y axis
    ratio->GetYaxis()->SetTitle("PbPb/pp");
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleOffset(2.5);
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(14);

    // Draw histograms
    TCanvas *c_nSD = new TCanvas("c_nSD","nSD",600,900);
    TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
    pad1->SetBottomMargin(0);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();
    hnSD_AA->Draw();
    hnSD_pp->Draw("same");
    auto legend = new TLegend(0.8,0.8,0.9,0.9);
    legend->AddEntry(hnSD_AA,"AA");
    legend->AddEntry(hnSD_pp,"pp");
    legend->Draw();
    c_nSD->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();
    ratio->Draw();
    c_nSD->cd();
  } //-------------------------------------------------------------
  if (plot_mass){
    cout << "Plotting mass" << endl;
    hmass_pp->Scale(1./hmass_pp->Integral());
    hmass_AA->Scale(1./hmass_AA->Integral());
    TH1F *ratio = (TH1F*)hmass_AA->Clone("Ratio");

    // Check normalisation within precision
    if (hmass_pp->Integral() -  1.0 > 0.01 || hmass_AA->Integral() - 1.0 > 0.01){
      cout << "WARNING: NORMALISATION PROBLEM!" << endl << "pp = " << hmass_pp->Integral() << " ( " << hmass_pp->Integral()-1.0  << ")" << endl 
        << "AA = " << hmass_AA->Integral() << " (type: " << hmass_AA->Integral()-1.0 << ")" << endl;
    }

    // Plot settings
    hmass_AA->SetStats(0);
    hmass_AA->SetMarkerStyle(kPlus);
    hmass_AA->SetTitle(Form("Jet mass, pt #in [%.0f,%.0f), R = %.1f", min_pt, max_pt, jetR));
    // y axis
    hmass_AA->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN}{dm}");
    hmass_AA->GetYaxis()->SetTitleFont(43);
    hmass_AA->GetYaxis()->SetTitleSize(20);//0.05);
    hmass_AA->GetYaxis()->SetTitleOffset(2.5);
    hmass_AA->GetYaxis()->SetLabelFont(43);
    hmass_AA->GetYaxis()->SetLabelSize(14);

    hmass_pp->SetStats(0);
    hmass_pp->SetMarkerStyle(kMultiply);
    hmass_pp->SetLineColor(kRed);
    hmass_pp->SetMarkerColor(kRed);

    ratio->Divide(hmass_pp);
    ratio->SetStats(0);
    ratio->SetMarkerStyle(kDot);
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerColor(kBlack);
    ratio->SetTitle("");
    // x axis
    ratio->GetXaxis()->SetTitle("mass(MeV/c)");
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleOffset(2.5);
    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelSize(14);
    // y axis
    ratio->GetYaxis()->SetTitle("PbPb/pp");
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleOffset(2.5);
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(14);

    // Draw histograms
    TCanvas *c_mass = new TCanvas("c_mass","mass",600,900);
    TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
    pad1->SetBottomMargin(0);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();
    hmass_AA->Draw();
    hmass_pp->Draw("same");
    auto legend = new TLegend(0.8,0.8,0.9,0.9);
    legend->AddEntry(hmass_AA,"AA");
    legend->AddEntry(hmass_pp,"pp");
    legend->Draw();
    c_mass->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();
    ratio->Draw();
    c_mass->cd();
  } //-------------------------------------------------------------
  if (plot_mz2){
    cout << "Plotting mz2" << endl;
    hmz2_pp->Scale(1./hmz2_pp->Integral());
    hmz2_AA->Scale(1./hmz2_AA->Integral());
    TH1F *ratio = (TH1F*)hmz2_AA->Clone("Ratio");

    // Check normalisation within precision
    if (hmz2_pp->Integral() -  1.0 > 0.01 || hmz2_AA->Integral() - 1.0 > 0.01){
      cout << "WARNING: NORMALISATION PROBLEM!" << endl << "pp = " << hmz2_pp->Integral() << " ( " << hmz2_pp->Integral()-1.0  << ")" << endl 
        << "AA = " << hmz2_AA->Integral() << " (type: " << hmz2_AA->Integral()-1.0 << ")" << endl;
    }

    // Plot settings
    hmz2_AA->SetStats(0);
    hmz2_AA->SetMarkerStyle(kPlus);
    hmz2_AA->SetTitle(Form("Jet mz2, pt #in [%.0f,%.0f), R = %.1f", min_pt, max_pt, jetR));
    // y axis
    hmz2_AA->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN}{dmz2}");
    hmz2_AA->GetYaxis()->SetTitleFont(43);
    hmz2_AA->GetYaxis()->SetTitleSize(20);//0.05);
    hmz2_AA->GetYaxis()->SetTitleOffset(2.5);
    hmz2_AA->GetYaxis()->SetLabelFont(43);
    hmz2_AA->GetYaxis()->SetLabelSize(14);

    hmz2_pp->SetStats(0);
    hmz2_pp->SetMarkerStyle(kMultiply);
    hmz2_pp->SetLineColor(kRed);
    hmz2_pp->SetMarkerColor(kRed);

    ratio->Divide(hmz2_pp);
    ratio->SetStats(0);
    ratio->SetMarkerStyle(kDot);
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerColor(kBlack);
    ratio->SetTitle("");
    // x axis
    ratio->GetXaxis()->SetTitle("mz2");
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleOffset(2.5);
    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelSize(14);
    // y axis
    ratio->GetYaxis()->SetTitle("PbPb/pp");
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleOffset(2.5);
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(14);

    // Draw histograms
    TCanvas *c_mz2 = new TCanvas("c_mz2","mz2",600,900);
    TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
    pad1->SetBottomMargin(0);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();
    hmz2_AA->Draw();
    hmz2_pp->Draw("same");
    auto legend = new TLegend(0.8,0.8,0.9,0.9);
    legend->AddEntry(hmz2_AA,"AA");
    legend->AddEntry(hmz2_pp,"pp");
    legend->Draw();
    c_mz2->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();
    ratio->Draw();
    c_mz2->cd();
  } //-------------------------------------------------------------
  if (plot_mr){
    cout << "Plotting mr" << endl;
    hmr_pp->Scale(1./hmr_pp->Integral());
    hmr_AA->Scale(1./hmr_AA->Integral());
    TH1F *ratio = (TH1F*)hmr_AA->Clone("Ratio");

    // Check normalisation within precision
    if (hmr_pp->Integral() -  1.0 > 0.01 || hmr_AA->Integral() - 1.0 > 0.01){
      cout << "WARNING: NORMALISATION PROBLEM!" << endl << "pp = " << hmr_pp->Integral() << " ( " << hmr_pp->Integral()-1.0  << ")" << endl 
        << "AA = " << hmr_AA->Integral() << " (type: " << hmr_AA->Integral()-1.0 << ")" << endl;
    }

    // Plot settings
    hmr_AA->SetStats(0);
    hmr_AA->SetMarkerStyle(kPlus);
    hmr_AA->SetTitle(Form("Jet mr, pt #in [%.0f,%.0f), R = %.1f", min_pt, max_pt, jetR));
    // y axis
    hmr_AA->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN}{dmr}");
    hmr_AA->GetYaxis()->SetTitleFont(43);
    hmr_AA->GetYaxis()->SetTitleSize(20);//0.05);
    hmr_AA->GetYaxis()->SetTitleOffset(2.5);
    hmr_AA->GetYaxis()->SetLabelFont(43);
    hmr_AA->GetYaxis()->SetLabelSize(14);

    hmr_pp->SetStats(0);
    hmr_pp->SetMarkerStyle(kMultiply);
    hmr_pp->SetLineColor(kRed);
    hmr_pp->SetMarkerColor(kRed);

    ratio->Divide(hmr_pp);
    ratio->SetStats(0);
    ratio->SetMarkerStyle(kDot);
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerColor(kBlack);
    ratio->SetTitle("");
    // x axis
    ratio->GetXaxis()->SetTitle("mr");
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleOffset(2.5);
    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelSize(14);
    // y axis
    ratio->GetYaxis()->SetTitle("PbPb/pp");
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleOffset(2.5);
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(14);

    // Draw histograms
    TCanvas *c_mr = new TCanvas("c_mr","mr",600,900);
    TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
    pad1->SetBottomMargin(0);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();
    hmr_AA->Draw();
    hmr_pp->Draw("same");
    auto legend = new TLegend(0.8,0.8,0.9,0.9);
    legend->AddEntry(hmr_AA,"AA");
    legend->AddEntry(hmr_pp,"pp");
    legend->Draw();
    c_mr->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();
    ratio->Draw();
    c_mr->cd();
  } //-------------------------------------------------------------
  if (plot_mr2){
    cout << "Plotting mr2" << endl;
    hmr2_pp->Scale(1./hmr2_pp->Integral());
    hmr2_AA->Scale(1./hmr2_AA->Integral());
    TH1F *ratio = (TH1F*)hmr2_AA->Clone("Ratio");

    // Check normalisation within precision
    if (hmr2_pp->Integral() -  1.0 > 0.01 || hmr2_AA->Integral() - 1.0 > 0.01){
      cout << "WARNING: NORMALISATION PROBLEM!" << endl << "pp = " << hmr2_pp->Integral() << " ( " << hmr2_pp->Integral()-1.0  << ")" << endl 
        << "AA = " << hmr2_AA->Integral() << " (type: " << hmr2_AA->Integral()-1.0 << ")" << endl;
    }

    // Plot settings
    hmr2_AA->SetStats(0);
    hmr2_AA->SetMarkerStyle(kPlus);
    hmr2_AA->SetTitle(Form("Jet mr2, pt #in [%.0f,%.0f), R = %.1f", min_pt, max_pt, jetR));
    // y axis
    hmr2_AA->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN}{dmr2}");
    hmr2_AA->GetYaxis()->SetTitleFont(43);
    hmr2_AA->GetYaxis()->SetTitleSize(20);//0.05);
    hmr2_AA->GetYaxis()->SetTitleOffset(2.5);
    hmr2_AA->GetYaxis()->SetLabelFont(43);
    hmr2_AA->GetYaxis()->SetLabelSize(14);

    hmr2_pp->SetStats(0);
    hmr2_pp->SetMarkerStyle(kMultiply);
    hmr2_pp->SetLineColor(kRed);
    hmr2_pp->SetMarkerColor(kRed);

    ratio->Divide(hmr2_pp);
    ratio->SetStats(0);
    ratio->SetMarkerStyle(kDot);
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerColor(kBlack);
    ratio->SetTitle("");
    // x axis
    ratio->GetXaxis()->SetTitle("mr2");
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleOffset(2.5);
    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelSize(14);
    // y axis
    ratio->GetYaxis()->SetTitle("PbPb/pp");
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleOffset(2.5);
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(14);

    // Draw histograms
    TCanvas *c_mr2 = new TCanvas("c_mr2","mr2",600,900);
    TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
    pad1->SetBottomMargin(0);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();
    hmr2_AA->Draw();
    hmr2_pp->Draw("same");
    auto legend = new TLegend(0.8,0.8,0.9,0.9);
    legend->AddEntry(hmr2_AA,"AA");
    legend->AddEntry(hmr2_pp,"pp");
    legend->Draw();
    c_mr2->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();
    ratio->Draw();
    c_mr2->cd();
  } //-------------------------------------------------------------
  if (plot_rz){
    cout << "Plotting rz" << endl;
    hrz_pp->Scale(1./hrz_pp->Integral());
    hrz_AA->Scale(1./hrz_AA->Integral());
    TH1F *ratio = (TH1F*)hrz_AA->Clone("Ratio");

    // Check normalisation within precision
    if (hrz_pp->Integral() -  1.0 > 0.01 || hrz_AA->Integral() - 1.0 > 0.01){
      cout << "WARNING: NORMALISATION PROBLEM!" << endl << "pp = " << hrz_pp->Integral() << " ( " << hrz_pp->Integral()-1.0  << ")" << endl 
        << "AA = " << hrz_AA->Integral() << " (type: " << hrz_AA->Integral()-1.0 << ")" << endl;
    }

    // Plot settings
    hrz_AA->SetStats(0);
    hrz_AA->SetMarkerStyle(kPlus);
    hrz_AA->SetTitle(Form("Jet rz, pt #in [%.0f,%.0f), R = %.1f", min_pt, max_pt, jetR));
    // y axis
    hrz_AA->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN}{drz}");
    hrz_AA->GetYaxis()->SetTitleFont(43);
    hrz_AA->GetYaxis()->SetTitleSize(20);//0.05);
    hrz_AA->GetYaxis()->SetTitleOffset(2.5);
    hrz_AA->GetYaxis()->SetLabelFont(43);
    hrz_AA->GetYaxis()->SetLabelSize(14);

    hrz_pp->SetStats(0);
    hrz_pp->SetMarkerStyle(kMultiply);
    hrz_pp->SetLineColor(kRed);
    hrz_pp->SetMarkerColor(kRed);

    ratio->Divide(hrz_pp);
    ratio->SetStats(0);
    ratio->SetMarkerStyle(kDot);
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerColor(kBlack);
    ratio->SetTitle("");
    // x axis
    ratio->GetXaxis()->SetTitle("rz");
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleOffset(2.5);
    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelSize(14);
    // y axis
    ratio->GetYaxis()->SetTitle("PbPb/pp");
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleOffset(2.5);
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(14);

    // Draw histograms
    TCanvas *c_rz = new TCanvas("c_rz","rz",600,900);
    TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
    pad1->SetBottomMargin(0);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();
    hrz_AA->Draw();
    hrz_pp->Draw("same");
    auto legend = new TLegend(0.8,0.8,0.9,0.9);
    legend->AddEntry(hrz_AA,"AA");
    legend->AddEntry(hrz_pp,"pp");
    legend->Draw();
    c_rz->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();
    ratio->Draw();
    c_rz->cd();
  } //-------------------------------------------------------------
  if (plot_r2z){
    cout << "Plotting r2z" << endl;
    hr2z_pp->Scale(1./hr2z_pp->Integral());
    hr2z_AA->Scale(1./hr2z_AA->Integral());
    TH1F *ratio = (TH1F*)hr2z_AA->Clone("Ratio");

    // Check normalisation within precision
    if (hr2z_pp->Integral() -  1.0 > 0.01 || hr2z_AA->Integral() - 1.0 > 0.01){
      cout << "WARNING: NORMALISATION PROBLEM!" << endl << "pp = " << hr2z_pp->Integral() << " ( " << hr2z_pp->Integral()-1.0  << ")" << endl 
        << "AA = " << hr2z_AA->Integral() << " (type: " << hr2z_AA->Integral()-1.0 << ")" << endl;
    }

    // Plot settings
    hr2z_AA->SetStats(0);
    hr2z_AA->SetMarkerStyle(kPlus);
    hr2z_AA->SetTitle(Form("Jet r2z, pt #in [%.0f,%.0f), R = %.1f", min_pt, max_pt, jetR));
    // y axis
    hr2z_AA->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN}{dr2z}");
    hr2z_AA->GetYaxis()->SetTitleFont(43);
    hr2z_AA->GetYaxis()->SetTitleSize(20);//0.05);
    hr2z_AA->GetYaxis()->SetTitleOffset(2.5);
    hr2z_AA->GetYaxis()->SetLabelFont(43);
    hr2z_AA->GetYaxis()->SetLabelSize(14);

    hr2z_pp->SetStats(0);
    hr2z_pp->SetMarkerStyle(kMultiply);
    hr2z_pp->SetLineColor(kRed);
    hr2z_pp->SetMarkerColor(kRed);

    ratio->Divide(hr2z_pp);
    ratio->SetStats(0);
    ratio->SetMarkerStyle(kDot);
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerColor(kBlack);
    ratio->SetTitle("");
    // x axis
    ratio->GetXaxis()->SetTitle("r2z");
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleOffset(2.5);
    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelSize(14);
    // y axis
    ratio->GetYaxis()->SetTitle("PbPb/pp");
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleOffset(2.5);
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(14);

    // Draw histograms
    TCanvas *c_r2z = new TCanvas("c_r2z","r2z",600,900);
    TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
    pad1->SetBottomMargin(0);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();
    hr2z_AA->Draw();
    hr2z_pp->Draw("same");
    auto legend = new TLegend(0.8,0.8,0.9,0.9);
    legend->AddEntry(hr2z_AA,"AA");
    legend->AddEntry(hr2z_pp,"pp");
    legend->Draw();
    c_r2z->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();
    ratio->Draw();
    c_r2z->cd();
  } //-------------------------------------------------------------
}

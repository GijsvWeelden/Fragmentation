#include "style.C"
#include "legend_utils.C"
#include "jewel_plot_utils.C"
void plot_jet_variables(void) {
  set_legend_font(62,0.04);
  // pt range
  Float_t min_pt = 40;
  Float_t max_pt = 60;
  Float_t jetR = 0.4;//2;//5;

  TLatex *ltx = new TLatex;
  ltx->SetNDC();
  ltx->SetTextSize(0.045);

  // TChain AA
  int numFiles = 3;
  TChain *chain = new TChain("jetprops");
  for (int fileNum = 1; fileNum <= numFiles; fileNum++) {
    chain->AddFile(Form("../run_AA_2tev76_norecoil/jet_shapes_constsub_eventwise_tree_%d_full_nobkg.root", fileNum));
  }

  // Initialise desired histograms
  TH1F *hmass_AA = new TH1F("hmass_AA", Form("PbPb2tev76 mass, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,30.0);
  TH1F *hmass_pp = new TH1F("hmass_pp", Form("pp2tev76 mass, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,30.0);

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

  for (Long64_t entry = 0; entry < nentries; entry++){
    Long64_t local_entry = chain->LoadTree(entry);
    if (local_entry < 0) break;
    nb = chain->GetEntry(entry);
    nbytes+=nb;
    if (pt >= min_pt && pt < max_pt){
      hmass_AA->Fill(mass,evwt);
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
      hmass_pp->Fill(mass,evwt);
    }
  }

  // Draw histograms
  // Make canvases
  TCanvas *c1 = new TCanvas("c1","c1: mass",600,600);
  hmass_pp->SetMarkerStyle(kFullDotLarge);
  hmass_pp->SetMarkerColor(kRed);
  hmass_pp->SetLineColor(kRed);
  hmass_AA->DrawNormalized("E1");
  hmass_pp->DrawNormalized("E1,SAMES");
  //hmass_AA->SetAxisRange(0, 3.5e-6, "Y");
  //hmass_pp->SetAxisRange(0, 3.5e-6, "Y");
  //TFile hmass_save("hmass.root","RECREATE");
  //hmass->Write();

}

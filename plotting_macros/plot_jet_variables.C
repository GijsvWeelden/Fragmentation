#include "style.C"
#include "legend_utils.C"
#include "jewel_plot_utils.C"
void plot_jet_variables(void) {
  set_legend_font(62,0.04);
  // pt range
  Float_t min_pt = 40;
  Float_t max_pt = 60;
  Float_t jetR = 0.4;//2;//5;

  // Make canvases
  TCanvas *c1 = new TCanvas("c1","c1: mass",600,600);
  //TCanvas *c2 = new TCanvas("c2","c2: angularity",600,600);
  //TCanvas *c3 = new TCanvas("c3","c3: PtD",600,600);

  TLatex *ltx = new TLatex;
  ltx->SetNDC();
  ltx->SetTextSize(0.045);

  // Make TChain
  int numFiles = 3;
  TChain *chain = new TChain("jetprops");
  for (int fileNum = 1; fileNum <= numFiles; fileNum++) {
    chain->AddFile(Form("../run_AA_2tev76_norecoil/jet_shapes_constsub_eventwise_tree_%d_full_nobkg.root", fileNum));
  }

  // Initialise desired histograms
  TH1F *hmass = new TH1F("hmass", Form("PbPb2tev76 mass, pt #in [%.0f,%.0f)", min_pt, max_pt),100,0.0,30.0);
  //hmass->SetName("hmass");
  //hmass->SetTitle(Form("PbPb2tev76 mass, pt #in [%.0f,%.0f)", min_pt, max_pt));

  Long64_t nentries = chain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  //cout << "nentries: " << nentries;

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

  int mass_counter = 0;
  float_t mass_avg = 0;

  for (Long64_t entry = 0; entry < nentries; entry++){
    Long64_t local_entry = chain->LoadTree(entry);
    if (local_entry < 0) break;
    nb = chain->GetEntry(entry);
    nbytes+=nb;
    // float pt = chain->GetEntry(entry);
    /*
       if (local_entry < 10){
       cout << "\npt: " << pt;
       cout << "\nentry: " << entry << ", local_entry: " << local_entry;// << ", pt: " << pt;
       }
     */
    //chain->Draw("mass");
    if (pt >= min_pt && pt < max_pt){
      hmass->Fill(mass);
    }
  }

  // TODO
  // Change to fit needs
  hmass->Draw();
  //TFile hmass_save("hmass.root","RECREATE");
  //hmass->Write();

  /*
  // TODO: how to make histograms?
  for (Int_t iFile = 0; iFile < nFiles; iFile++) {
  c1->cd();
  TFile *fin = new TFile(fnames[iFile]);
  if (!fin || !fin->IsOpen())
  cout << "Error opening file " << fnames[iFile] << endl;
  //pt_eta_mass = (TH3F*) fin->Get("hLeadJetPtEtaMass_R04");
  auto pt_eta_mass = (TH3F*) fin->Get(Form("hJetPtEtaMass_R%02d",Int_t(10*jetR)));
  TH1 *hmass = 0;
  TH2F *pt_mass = 0;
  if (pt_eta_mass != 0) {
  hmass = project_shape( pt_eta_mass,min_pt,max_pt,Form("mass_%d",iFile));
  }
  else {
  //pt_mass = (TH2F*) fin->Get("hPtLeadJetMass_R04");
  pt_mass = (TH2F*) fin->Get(Form("hPtJetMass_R%02d",Int_t(10*jetR)));
  hmass = project_shape_2D(pt_mass, min_pt, max_pt, Form("mass_%d",iFile));
  }
  hmass->UseCurrentStyle();
  hmass->GetYaxis()->SetNdivisions(503);
  hmass->SetMarkerColor(colors[iFile]);
  hmass->SetLineColor(colors[iFile]);
  hmass->SetMarkerStyle(markers[iFile]);
  if (iFile == 0) {
  hmass->Draw();
  hmass->SetYTitle("1/#it{N}_{jet} d#it{N}/d#it{M}_{jet} (#it{c}^{2}/GeV)");
  hmass->SetXTitle("#it{M}_{jet} (GeV/#it{c}^{2})");
  hmass->SetMaximum(0.249);
  }
  else
  hmass->Draw("same");

  gPad->Update();
  if (iFile == 0) {
  ltx->DrawLatex(ndc_x(0.05),ndc_y(0.93),"JEWEL");
  ltx->DrawLatex(ndc_x(0.4),ndc_y(0.93),Form("#it{R}=%.1f, %.0f < #it{p}_{T, jet}^{ch} < %.0f",jetR,min_pt,max_pt));
  }
  draw_legend_m(0.4,0.86,hmass,labels[iFile],iFile);

  //
  // Angularity
  //
  c2->cd();
  //pt_eta_ang = (TH3F*) fin->Get("hLeadJetPtEtaAng_R04");
  auto pt_eta_ang = (TH3F*) fin->Get(Form("hJetPtEtaAng_R%02d",Int_t(10*jetR)));
  TH1 *hang = 0;
  TH2F *pt_ang = 0;
  if (pt_eta_ang != 0) {
  hang = project_shape( pt_eta_ang,min_pt,max_pt,Form("ang_%d",iFile));
  }
  else {
  //pt_ang = (TH2F*) fin->Get("hPtLeadJetAng_R04");
  pt_ang = (TH2F*) fin->Get(Form("hPtJetAng_R%02d",Int_t(10*jetR)));
  if (pt_ang == 0)
  cout << "No angularity hist found!!!" << endl;
  hang = project_shape_2D(pt_ang, min_pt, max_pt, Form("ang_%d",iFile));
  }
  hang->UseCurrentStyle();
  hang->SetMarkerColor(colors[iFile]);
  hang->SetLineColor(colors[iFile]);
  hang->SetMarkerStyle(markers[iFile]);
  if (iFile == 0) {
  hang->Draw();
  hang->SetYTitle("1/#it{N}_{jet} d#it{N}/d#it{g}");
  hang->SetXTitle("#it{g}");
  hang->GetXaxis()->SetRangeUser(0,0.3);
  hang->SetMaximum(25);
  }
  else
  hang->Draw("same");

  gPad->Update();
  if (iFile == 0) {
    ltx->DrawLatex(ndc_x(0.05),ndc_y(0.93),"JEWEL");
    ltx->DrawLatex(ndc_x(0.4),ndc_y(0.93),Form("#it{R}=%.1f, %.0f < #it{p}_{T, jet}^{ch} < %.0f",jetR,min_pt,max_pt));
  }
  draw_legend_m(0.4,0.86,hang,labels[iFile],iFile);

  //
  //   p_T,D
  // 

  //cout << "File " << fin->GetName() << " integral of ptd hist" << 
  c3->cd();
  //pt_eta_PtD = (TH3F*) fin->Get("hLeadJetPtEtaPtD_R04");
  auto pt_eta_PtD = (TH3F*) fin->Get(Form("hJetPtEtaPtD_R%02d",Int_t(10*jetR)));
  TH1 *hPtD = 0;
  TH2F *pt_PtD = 0;
  if (pt_eta_PtD != 0) {
    hPtD = project_shape( pt_eta_PtD,min_pt,max_pt,Form("PtD_%d",iFile));
  }
  else {
    pt_PtD = (TH2F*) fin->Get(Form("hPtLeadJetpTD_R%02d",Int_t(10*jetR)));
    if (pt_PtD == 0)
      cout << "No PtD hist found!!!" << endl;
    hPtD = project_shape_2D(pt_PtD, min_pt, max_pt, Form("PtD_%d",iFile));
  }
  hPtD->SetMarkerColor(colors[iFile]);
  hPtD->SetLineColor(colors[iFile]);
  hPtD->SetMarkerStyle(markers[iFile]);
  if (iFile == 0) {
    auto h1 = gPad->DrawFrame(0,0,1.3,6);
    h1->SetXTitle("#it{p}_{T,D}");
    h1->SetYTitle("1/#it{N}_{jet} d#it{N}/d#it{p}_{T,D}");
    hPtD->Draw("same");
    //hPtD->SetMaximum(0.24);
  }
  else //if (iFile < 3)
    hPtD->Draw("same");

  gPad->Update();
  if (iFile == 0) {
    ltx->DrawLatex(ndc_x(0.05),ndc_y(0.93),"JEWEL");
    ltx->DrawLatex(ndc_x(0.4),ndc_y(0.93),Form("#it{R}=%.1f, %.0f < #it{p}_{T, jet}^{ch} < %.0f",jetR,min_pt,max_pt));
  }
  //if (iFile < 3)
  draw_legend_m(0.4,0.86,hPtD,labels[iFile],iFile);

}*/
}

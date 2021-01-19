#include "style.C"
#include "legend_utils.C"
#include "jewel_plot_utils.C"
void plot_shapes(void) {
  //gROOT->Macro("~/macros/style.C");
  /*
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLineWidth(2);
  gStyle->SetOptTitle(0);
  gStyle->SetNdivisions(505);
  gStyle->SetNdivisions(505,"Y");
  gStyle->SetTitleOffset(1.2,"X");
  gStyle->SetTitleOffset(1.1,"Y");
  gStyle->SetTitleSize(0.04,"X");
  gStyle->SetTitleSize(0.04,"Y");
  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadTopMargin(0.05);
  */

  //gROOT->LoadMacro("~/macros/legend_utils.C");
  set_legend_font(62,0.04);
  //gROOT->LoadMacro("jewel_plot_utils.C");

  //TFile *fin_recoil0 = new TFile("output_0_10_recoil0/jet_shapes_out_0_10_recoil0_all.root");
  //TFile *fin_med = new TFile("output_0_10_recoil0/jet_shapes_out_0_10_medium_recoil0_all_Escheme_mass.root");
  //TString medLabel("Pb+Pb, no rc");
  //Char_t *outLabel = "charged_noRec_R04_40pt60";

  /*
  Char_t *fnames[] = {"output_0_10_vacuum/jet_shapes_out_all_vacuum.hepmc.root",
		      "output_0_10_recoil0/jet_shapes_out_all_medium_recoil0.hepmc.root",
		      "output_0_10_recoil1/jet_shapes_bkgsub_out_all_medium_recoil1.hepmc.root"};
  */
  //"output_0_10_recoil1/jet_shapes_constsub_out_all_medium_recoil1.hepmc.root"};

  // Used for Hard Probes results: 
  /*
  Char_t *fnames[] = {"analysis/results_May2015/jet_shapes_out_all_vacuum.hepmc.root",
		      "analysis/results_May2015/jet_shapes_out_all_medium_recoil0.hepmc.root",
		      "analysis/results_May2015/jet_shapes_bkgsub_out_all_medium_recoil1.hepmc.root",
                      "analysis/results_May2015/jet_shapes_constsub_out_all_medium_recoil1.hepmc.root"};
  //"output_0_10_recoil1/jet_shapes_constsub_out_all_medium_recoil1.hepmc.root"}
		     //"output_0_10_recoil1/jet_shapes_bkgsub_all_rhom_medium_recoil1.hepmc.root"};
  Char_t *labels[] = {"pp","PbPb no recoil","PbPb recoil (der sub)","PbPb recoil (const sub)"};
 */
 
  Char_t *fnames[] = {"analysis/results_2017/jet_shapes_pp_2tev76_nobkgsub_charged.root",
		      //"analysis/results_2017/jet_shapes_AA_Tdef_0cent10_norecoil_nobkgsub_charged.root",
		      "analysis/results_2017/jet_shapes_AA_Tdef_0cent10_recoil_bkgsub_charged.root",
		      //"analysis/results_2017/jet_shapes_AA_Tdef_0cent10_recoil_constsub_eventwise_full_deltaRcut.root"};
		      //"analysis/results_2017/jet_shapes_AA_Tdef_0cent10_recoil_constsub_eventwise_charged_deltaRcut.root"};
		      "analysis/results_2017/jet_shapes_AA_Tdef_0cent10_recoil_constsub_eventwise_charged.root",
		      //"analysis/results_2017/jet_shapes_AA_Tdef_0cent10_recoil_constsub_eventwise_maxR1_charged.root",
                      "analysis/results_2017/jet_shapes_AA_Tdef_0cent10_recoil_constsub_charged.root"};
  
  Char_t *labels[] = {"pp",/*"PbPb no recoil",*/"PbPb recoil (der sub)","PbPb recoil (const sub evt)","PbPb recoil (const sub)"};

  Int_t colors[] = {4,2,kRed,kOrange};
  Int_t markers[] = {20,24,21,22};

  Int_t nFiles = sizeof(fnames)/sizeof(fnames[0]);
  cout << "nFiles " << nFiles << endl;

  Float_t min_pt = 40;
  Float_t max_pt = 60;
  Float_t jetR = 0.5;//2;//5;

  TCanvas *c1 = new TCanvas("c1","c1: mass",600,600);
  TCanvas *c2 = new TCanvas("c2","c2: angularity",600,600);
  TCanvas *c3 = new TCanvas("c3","c3: PtD",600,600);
  
  TLatex *ltx = new TLatex;
  ltx->SetNDC();
  ltx->SetTextSize(0.045);

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

  }
}

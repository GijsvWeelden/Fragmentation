
#include <vector> 
#include <iostream> 
#include <time.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"

//-------------------------------------------------------------
//
// jetprops_hists.C extracts data from trees (chains),
// and stores it in histograms, for several pt bins
//
//-------------------------------------------------------------

void make_hists(TChain *chain, string setting, vector<string> obs, vector<double> ptBins){
  // Setting = pp, AA_norecoil, AA_recoil

  for (auto iobs=0; iobs<obs.size(); iobs++){
    // Declare lists
    TList* list = new TList();
    list->SetName(TString::Format("%s_%s",
          setting.c_str(), obs[iobs].c_str()).Data());
    list->SetOwner(true);

    // Loop over pt bins
    for (auto ipt=0; ipt<ptBins.size()-1; ipt++){
      TH1F *hist = new TH1F(TString::Format("%s_h%s_pt%.f_%.f",
            setting.c_str(), obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data(),
          TString::Format("Jet %s (%s), pt #in [%.f_%.f), R = 0.4",
            obs[iobs].c_str(), setting.c_str() ,ptBins[ipt], ptBins[ipt+1]).Data(),
          100,0,0.5);
      hist->Sumw2();
      if (obs[iobs] == "dphi"){
        hist->SetBins(100,0,3.2);
      }
      else if (obs[iobs] == "nconst" || obs[iobs] == "nSD" || obs[iobs] == "mass"){
        hist->SetBins(100,0,100);
      }
      chain->Draw(TString::Format("%s>>%s_h%s_pt%.f_%.f",
            obs[iobs].c_str(), setting.c_str(), obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data(),
          TString::Format("evwt*(pt>%.f && pt<%.f)",
            ptBins[ipt], ptBins[ipt+1]).Data());
      list->Add(hist);
    }
    list->Write(TString::Format("%s_%s",
          setting.c_str(), obs[iobs].c_str()).Data(),1);
  }
}

//-------------------------------------------------------------
//
// Main Function
//
//-------------------------------------------------------------
void jetprops_hists(void){
  double time = clock();

  // pt bins and observables
  std::vector<double> ptBins = {0,20,40,60,80,100,120,160,200};
  std::vector<string> obs = {"dphi","nconst","zg","Rg","nSD","mass","mz2","mr","mr2","rz","r2z"};

  // TODO: make these input from shell
  int numFiles_AAnr = 20;
  int numFiles_AAr = 0;
  int numFiles_pp = 10;

  // Declare output file
  // TODO: make this input from shell
  string outName = "compare_2tev76_ppAAnr";
  TFile *outFile = new TFile(Form("%s.root", outName.c_str()),"RECREATE");
  outFile->cd();

  //-------------------------------------------------------------
  //
  // AA (no recoil)
  //
  //-------------------------------------------------------------

  // TChain AA
  TChain *chain = new TChain("jetprops");
  for (int fileNum = 1; fileNum <= numFiles_AAnr; fileNum++) {
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

  if (numFiles_AAnr != 0) make_hists(chain, "AA_norecoil", obs, ptBins);

  //-------------------------------------------------------------
  //
  // AA (recoil)
  //
  //-------------------------------------------------------------

  // TChain AA recoil
  // int numFiles_pp = 1;
  chain->Reset();
  for (int fileNum = 1; fileNum <= numFiles_AAr; fileNum++) {
    chain->AddFile(Form("../run_AA_2tev76_recoil/jet_shapes_constsub_eventwise_tree_%d_full.root", fileNum));
  }

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

  make_hists(chain, "AA_recoil", obs, ptBins);

  //-------------------------------------------------------------
  //
  // pp
  //
  //-------------------------------------------------------------

  // TChain pp
  chain->Reset();
  //chain->AddFile("../run_pp_2tev76/jet_tree_pp2tev76_full_nobkg.root");
  for (int fileNum = 1; fileNum <= numFiles_pp; fileNum++) {
    chain->AddFile(Form("../run_pp_2tev76/jet_shapes_constsub_eventwise_tree_%d_full_nobkg.root", fileNum));
  }

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

  make_hists(chain, "pp", obs, ptBins);

  //-------------------------------------------------------------
  //
  // End of file
  //
  //-------------------------------------------------------------

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << endl;
  outFile->Close();
}

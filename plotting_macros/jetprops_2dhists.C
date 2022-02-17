
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
// and stores it in 2d histograms (pt - obs).
// This code distinguishes between leading and away-side jets (dphi > pi/2)
//
// If you want to look at a single variable, using the histograms
// is much faster than reading in the whole tree
//
//-------------------------------------------------------------

void make_hists(TChain *chain, string setting, vector<string> obs, TList *list){
  // Setting = pp, AA_norecoil, AA_recoil

  double pi = 3.14159265358979323846;

  for (auto iobs=0; iobs<obs.size(); iobs++){
    TH2F *hL = new TH2F(TString::Format("%s_%s_Leading",
          setting.c_str(), obs[iobs].c_str()).Data(),
        TString::Format("%s leading jet %s", setting.c_str(), obs[iobs].c_str()).Data(),
        100, 0, 0.5, 200, 0, 200);
    hL->Sumw2();
    TH2F *hS = new TH2F(TString::Format("%s_%s_Subleading",
          setting.c_str(), obs[iobs].c_str()).Data(),
        TString::Format("%s wayside jet(s) %s (d#phi>#pi/2)", setting.c_str(), obs[iobs].c_str()).Data(),
        100, 0, 0.5, 200, 0, 200);
    hS->Sumw2();

    hL->GetYaxis()->SetTitle("pt");
    hL->GetXaxis()->SetTitle(TString::Format("%s", obs[iobs].c_str()).Data());

    hS->GetYaxis()->SetTitle("pt");
    hS->GetXaxis()->SetTitle(TString::Format("%s", obs[iobs].c_str()).Data());

    if (obs[iobs] == "dphi"){
      hL->SetBins(100, 0, 3.2, 200, 0, 200);
      hS->SetBins(100, 0, 3.2, 200, 0, 200);
    }
    else if (obs[iobs] == "nconst" || obs[iobs] == "nSD" || obs[iobs] == "mass"){
      hL->SetBins(100, 0, 100, 200, 0, 200);
      hS->SetBins(100, 0, 100, 200, 0, 200);
    }

    chain->Draw(TString::Format("pt:%s>>%s_%s_Leading",
          obs[iobs].c_str(), setting.c_str(), obs[iobs].c_str()).Data(),
        "evwt*(ijet==0)");
    chain->Draw(TString::Format("pt:%s>>%s_%s_Subleading",
          obs[iobs].c_str(), setting.c_str(), obs[iobs].c_str()).Data(),
        TString::Format("evwt*(ijet>0 && dphi>%f)",pi/2).Data());

    list->Add(hL);
    list->Add(hS);
  }
  cout << TString::Format("Name: %s", setting.c_str()).Data() << endl;
}

//-------------------------------------------------------------
//
// Main Function
//
//-------------------------------------------------------------
void jetprops_2dhists(void){
  double time = clock();

  // pt bins and observables
  std::vector<string> obs = {"dphi","nconst","zg","Rg","nSD","mass","mz2","mr","mr2","rz","r2z","t2t1","t3t2","t2dist","t3dist"};

  // TODO: make these input from shell
  int numFiles_AAnr = 20;
  int numFiles_AAr = 20;
  int numFiles_pp = 10;

  // Declare output file
  // TODO: make this input from shell
  string outName = "2dhists_2tev76_ppAAnrAAr";
  TFile *outFile = new TFile(Form("%s.root", outName.c_str()),"RECREATE");
  outFile->cd();

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
  Float_t         t2t1;
  Float_t         t3t2;
  Float_t         t2dist;
  Float_t         t3dist[3];
  // List of branches
  TBranch        *b_ievt;   //!
  TBranch        *b_ijet;   //!
  TBranch        *b_evwt;   //!
  TBranch        *b_pt;     //!
  TBranch        *b_eta;    //!
  TBranch        *b_phi;    //!
  TBranch        *b_dphi;   //!
  TBranch        *b_nconst; //!
  TBranch        *b_zg;     //!
  TBranch        *b_Rg;     //!
  TBranch        *b_nSD;    //!
  TBranch        *b_mass;   //!
  TBranch        *b_mz2;    //!
  TBranch        *b_mr;     //!
  TBranch        *b_mr2;    //!
  TBranch        *b_rz;     //!
  TBranch        *b_r2z;    //!
  TBranch        *b_t2t1;   //!
  TBranch        *b_t3t2;   //!
  TBranch        *b_t2dist; //!
  TBranch        *b_t3dist; //!

  //-------------------------------------------------------------
  //
  // AA (no recoil)
  //
  //-------------------------------------------------------------

  // TChain AA_norecoil
  TChain *chain = new TChain("jetprops");
  for (int fileNum = 1; fileNum <= numFiles_AAnr; fileNum++) {
    chain->AddFile(Form("../run_AA_2tev76_norecoil/jet_shapes_constsub_eventwise_tree_%d_full_nobkg.root", fileNum));
  }
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
  chain->SetBranchAddress("t2t1", &t2t1, &b_t2t1);
  chain->SetBranchAddress("t3t2", &t3t2, &b_t3t2);
  chain->SetBranchAddress("t2dist", &t2dist, &b_t2dist);
  chain->SetBranchAddress("t3dist", &t3dist, &b_t3dist);

  if (numFiles_AAnr != 0){
    TList *list = new TList();
    list->SetName("AA_norecoil");
    list->SetOwner(true);
    make_hists(chain, "AA_norecoil", obs, list);
    list->Write("AA_norecoil",1);
    delete list;
  }

  //-------------------------------------------------------------
  //
  // AA (recoil)
  //
  //-------------------------------------------------------------
  chain->Reset();
  for (int fileNum = 1; fileNum <= numFiles_AAr; fileNum++) {
    chain->AddFile(Form("../run_AA_2tev76_recoil/jet_shapes_constsub_eventwise_tree_%d_full_nobkg.root", fileNum));
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
  chain->SetBranchAddress("t2t1", &t2t1, &b_t2t1);
  chain->SetBranchAddress("t3t2", &t3t2, &b_t3t2);
  chain->SetBranchAddress("t2dist", &t2dist, &b_t2dist);
  chain->SetBranchAddress("t3dist", &t3dist, &b_t3dist);

  if (numFiles_AAr != 0){
    TList *list = new TList();
    list->SetName("AA_recoil");
    list->SetOwner(true);
    make_hists(chain, "AA_recoil", obs, list);
    list->Write("AA_recoil",1);
    delete list;
  }

  //-------------------------------------------------------------
  //
  // pp
  //
  //-------------------------------------------------------------
  chain->Reset();
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
  chain->SetBranchAddress("t2t1", &t2t1, &b_t2t1);
  chain->SetBranchAddress("t3t2", &t3t2, &b_t3t2);
  chain->SetBranchAddress("t2dist", &t2dist, &b_t2dist);
  chain->SetBranchAddress("t3dist", &t3dist, &b_t3dist);

  if (numFiles_pp != 0){
    TList *list = new TList();
    list->SetName("pp");
    list->SetOwner(true);
    make_hists(chain, "pp", obs, list);
    list->Write("pp",1);
    delete list;
  }

  //-------------------------------------------------------------
  //
  // End of file
  //
  //-------------------------------------------------------------

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds" << endl;
  outFile->Close();
}

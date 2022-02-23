
#include <vector>
#include <iostream>
#include <time.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"

//-------------------------------------------------------------
//
// Reads in tree with jet data. Saves the variables in 2D hists of type pt:observable
// Distinguishes between leading and away-side jets and also considers the complete jet sample
// Saves hists in 2dhists_2tev76_ppAAnrAAr
//
//-------------------------------------------------------------

void read_chain(TChain *chain, string setting, vector<string> obs);
void save_hists(TChain *chain, string setting, vector<string> obs);
TH2F *make_hists(string name, string title, string obs);

void trees_to_hists(void){
 double time = clock();
 std::vector<string> observables = {"dphi","nconst","zg","Rg","nSD","mass","mz2","mr","mr2","rz","r2z","ptD","t2t1","t3t2","t2dist","t3dist"};
  // TODO: make these input from shell
  int numFiles_AAnr = 20;
  int numFiles_AAr = 20;
  int numFiles_pp = 10;
  string outName = "2dhists_2tev76_ppAAnrAAr";
  TFile *outFile = new TFile(Form("%s.root", outName.c_str()),"RECREATE");
  outFile->cd();

  TChain *chain = new TChain("jetprops");
  if (numFiles_AAnr > 0){
    for (int fileNum = 1; fileNum <= numFiles_AAnr; fileNum++) {
      chain->AddFile(Form("../run_AA_2tev76_norecoil/jet_shapes_constsub_eventwise_tree_%d_full_nobkg.root", fileNum));
    }
    read_chain(chain, "AA_norecoil", observables);
  }

  chain->Reset();
  if (numFiles_AAr > 0){
    for (int fileNum = 1; fileNum <= numFiles_AAr; fileNum++) {
      chain->AddFile(Form("../run_AA_2tev76_recoil/jet_shapes_constsub_eventwise_tree_%d_full_nobkg.root", fileNum));
    }
    read_chain(chain, "AA_recoil", observables);
  }

  chain->Reset();
  if (numFiles_pp > 0){
    for (int fileNum = 1; fileNum <= numFiles_pp; fileNum++) {
      chain->AddFile(Form("../run_pp_2tev76/jet_shapes_constsub_eventwise_tree_%d_full_nobkg.root", fileNum));
    }
    read_chain(chain, "pp", observables);
  }

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << endl;
  // outFile->Close();
}

void read_chain(TChain *chain, string setting, vector<string> obs){
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
  Float_t         ptD;
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
  TBranch        *b_ptD;    //!
  TBranch        *b_t2t1;   //!
  TBranch        *b_t3t2;   //!
  TBranch        *b_t2dist; //!
  TBranch        *b_t3dist; //!
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
  chain->SetBranchAddress("ptD", &ptD, &b_ptD);
  chain->SetBranchAddress("t2t1", &t2t1, &b_t2t1);
  chain->SetBranchAddress("t3t2", &t3t2, &b_t3t2);
  chain->SetBranchAddress("t2dist", &t2dist, &b_t2dist);
  chain->SetBranchAddress("t3dist", &t3dist, &b_t3dist);

  save_hists(chain, setting, obs);
}

void save_hists(TChain *chain, string setting, vector<string> obs){
  TList *list = new TList();
  list->SetName(setting.c_str());
  list->SetOwner(true);

  for (auto iobs=0; iobs<obs.size(); iobs++){
    if (obs[iobs] == "t3dist") continue; // Is an array. Handled separately
    TH2F *hComp = make_hists(TString::Format("%s_%s",
                                       setting.c_str(), obs[iobs].c_str()).Data(),
                       TString::Format("Jet %s (%s)",
                                       obs[iobs].c_str(), setting.c_str()).Data(),
                       obs[iobs].c_str()
                       );
    TH2F *hLead = make_hists(TString::Format("%s_%s_leading",
                                       setting.c_str(), obs[iobs].c_str()).Data(),
                       TString::Format("Leading jet %s (%s)",
                                       obs[iobs].c_str(), setting.c_str()).Data(),
                       obs[iobs].c_str()
                       );
    TH2F *hAway = make_hists(TString::Format("%s_%s_awayside",
                                       setting.c_str(), obs[iobs].c_str()).Data(),
                       TString::Format("Awayside jet %s (%s)",
                                       obs[iobs].c_str(), setting.c_str()).Data(),
                       obs[iobs].c_str()
                       );
    // Draw hComp
    chain->Draw(TString::Format("pt:%s>>%s_%s",
                                obs[iobs].c_str(),
                                setting.c_str(),
                                obs[iobs].c_str()
                                ).Data()
                );
    // Draw hLead
    chain->Draw(TString::Format("pt:%s>>%s_%s_leading",
                                obs[iobs].c_str(),
                                setting.c_str(),
                                obs[iobs].c_str()
                                ).Data(),
                                "(evwt*(ijet==0))"
                );
    // Draw hAway
    chain->Draw(TString::Format("pt:%s>>%s_%s_awayside",
                                obs[iobs].c_str(),
                                setting.c_str(),
                                obs[iobs].c_str()
                                ).Data(),
                                TString::Format("(evwt*(ijet>0 && dphi>%f))", // TODO: Should be abs(dphi)
                                                TMath::PiOver2()).Data()
    );
    list->Add(hComp);
    list->Add(hLead);
    list->Add(hAway);
  }
  for (int i=0; i < 3; ++i){
    string obs = "t3dist[" + std::to_string(i) + "]";
    TH2F *hComp = make_hists(TString::Format("%s_%s",
                                             setting.c_str(), obs.c_str()).Data(),
                             TString::Format("Jet %s (%s)",
                                             obs.c_str(), setting.c_str()).Data(),
                             obs.c_str()
                             );
    TH2F *hLead = make_hists(TString::Format("%s_%s_leading",
                                             setting.c_str(), obs.c_str()).Data(),
                             TString::Format("Leading jet %s (%s)",
                                             obs.c_str(), setting.c_str()).Data(),
                             obs.c_str()
                             );
    TH2F *hAway = make_hists(TString::Format("%s_%s_awayside",
                                             setting.c_str(), obs.c_str()).Data(),
                             TString::Format("Awayside jet %s (%s)",
                                             obs.c_str(), setting.c_str()).Data(),
                             obs.c_str()
                             );
    // Draw hComp
    chain->Draw(TString::Format("pt:%s>>%s_%s",
                                obs.c_str(),
                                setting.c_str(),
                                obs.c_str()
                                ).Data()
                );
    // Draw hLead
    chain->Draw(TString::Format("pt:%s>>%s_%s_leading",
                                obs.c_str(),
                                setting.c_str(),
                                obs.c_str()
                                ).Data(),
                                "(evwt*(ijet==0))"
                );
    // Draw hAway
    chain->Draw(TString::Format("pt:%s>>%s_%s_awayside",
                                obs.c_str(),
                                setting.c_str(),
                                obs.c_str()
                                ).Data(),
                                TString::Format("(evwt*(ijet>0 && dphi>%f))",
                                                TMath::PiOver2()
                                                ).Data()
                );
    list->Add(hComp);
    list->Add(hLead);
    list->Add(hAway);
  }
  cout << TString::Format("Name: %s", setting.c_str()).Data() << endl;
  list->Write(setting.c_str(),1);
  delete list;
}

TH2F *make_hists(string name, string title, string obs){
  TH2F *hist = new TH2F(name.c_str(), title.c_str(),
                        100, 0, 0.5, 200, 0, 200);
  hist->Sumw2();
  hist->GetYaxis()->SetTitle("pt");
  hist->GetXaxis()->SetTitle(TString::Format("%s", obs.c_str()).Data());

  if (obs == "dphi"){
    hist->SetBins(100, 0, 3.2, 200, 0, 200);
  }
  else if (obs == "nconst" || obs == "nSD" || obs == "mass"){
    hist->SetBins(100, 0, 100, 200, 0, 200);
  }
  else if (obs == "t2dist" || obs == "t3dist"){
    hist->SetBins(100, 0, 2, 200, 0, 200);
  }
  else if (obs == "ptD" || obs == "t2dist" || obs == "t3dist[0]" || obs == "t3dist[1]" || obs == "t3dist[2]"){
    hist->SetBins(100, 0, 1, 200, 0, 200);
  }
  return hist;
}

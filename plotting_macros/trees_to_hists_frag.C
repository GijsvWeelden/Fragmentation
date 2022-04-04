
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
  std::vector<string> observables = {"frag", "orth"};
  string suffix = "charged_nobkg"; // charged/full, nobkg/<nothing>
  int numFiles_AAnr = 20;
  int numFiles_AAr = 20;
  int numFiles_pp = 10;
  string outName = "2dhists_5tev02_ppAAnrAAr";
  TFile *outFile = new TFile(Form("%s_%s.root", outName.c_str(), suffix.c_str()),"RECREATE");
  outFile->cd();

  TChain *chain = new TChain("jetprops");
  if (numFiles_AAnr > 0){
    for (int fileNum = 1; fileNum <= numFiles_AAnr; fileNum++) {
      chain->AddFile(Form("../run_AA_2tev76_norecoil/jet_frag_%d_%s.root", fileNum, suffix.c_str()));
    }
    read_chain(chain, "AA_norecoil", observables);
  }

  chain->Reset();
  if (numFiles_AAr > 0){
    for (int fileNum = 1; fileNum <= numFiles_AAr; fileNum++) {
      chain->AddFile(Form("../run_AA_2tev76_recoil/jet_frag_%d_%s.root", fileNum, suffix.c_str()));
    }
    read_chain(chain, "AA_recoil", observables);
  }

  chain->Reset();
  if (numFiles_pp > 0){
    for (int fileNum = 1; fileNum <= numFiles_pp; fileNum++) {
      chain->AddFile(Form("../run_pp_2tev76/jet_frag_%d_%s.root", fileNum, suffix.c_str()));
    }
    read_chain(chain, "pp", observables);
  }

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << endl;
  outFile->Close();
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
  Float_t         frag;
  Float_t         orth;
  Int_t           pdg;
  // List of branches
  TBranch        *b_ievt;   //!
  TBranch        *b_ijet;   //!
  TBranch        *b_evwt;   //!
  TBranch        *b_pt;     //!
  TBranch        *b_eta;    //!
  TBranch        *b_phi;    //!
  TBranch        *b_dphi;   //!
  TBranch        *b_nconst; //!
  TBranch        *b_frag;   //!
  TBranch        *b_orth;   //!
  TBranch        *b_pdg;    //!
  // Set branch addresses
  chain->SetBranchAddress("ievt", &ievt, &b_ievt);
  chain->SetBranchAddress("ijet", &ijet, &b_ijet);
  chain->SetBranchAddress("evwt", &evwt, &b_evwt);
  chain->SetBranchAddress("pt", &pt, &b_pt);
  chain->SetBranchAddress("eta", &eta, &b_eta);
  chain->SetBranchAddress("phi", &phi, &b_phi);
  chain->SetBranchAddress("dphi", &dphi, &b_dphi);
  chain->SetBranchAddress("nconst", &nconst, &b_nconst);
  chain->SetBranchAddress("frag", &frag, &b_frag);
  chain->SetBranchAddress("orth", &orth, &b_orth);
  chain->SetBranchAddress("pdg", &pdg, &b_pdg);

  save_hists(chain, setting, obs);
}

void save_hists(TChain *chain, string setting, vector<string> obs){
  TList *list = new TList();
  list->SetName(setting.c_str());
  list->SetOwner(true);

  for (auto iobs=0; iobs<obs.size(); iobs++){
    TH2F *hIncl = make_hists(TString::Format("%s_%s",
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
    // Draw hIncl
    chain->Draw(TString::Format("pt:%s>>%s_%s",
                                obs[iobs].c_str(),
                                setting.c_str(),
                                obs[iobs].c_str()
                                ).Data(),
                                "evwt"
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
                                TString::Format("(evwt*(ijet>0 && abs(dphi)>%f))",
                                                TMath::PiOver2()).Data()
    );
    list->Add(hIncl);
    list->Add(hLead);
    list->Add(hAway);
  }

  for (auto iobs=0; iobs<obs.size(); iobs++){
    TH2F *hIncl = make_hists(TString::Format("%s_%s",
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
    // Draw hIncl
    chain->Draw(TString::Format("nconst:%s>>%s_%s",
                                obs[iobs].c_str(),
                                setting.c_str(),
                                obs[iobs].c_str()
                                ).Data(),
                                "evwt"
                );
    // Draw hLead
    chain->Draw(TString::Format("nconst:%s>>%s_%s_leading",
                                obs[iobs].c_str(),
                                setting.c_str(),
                                obs[iobs].c_str()
                                ).Data(),
                                "(evwt*(ijet==0))"
                );
    // Draw hAway
    chain->Draw(TString::Format("nconst:%s>>%s_%s_awayside",
                                obs[iobs].c_str(),
                                setting.c_str(),
                                obs[iobs].c_str()
                                ).Data(),
                                TString::Format("(evwt*(ijet>0 && abs(dphi)>%f))",
                                                TMath::PiOver2()).Data()
    );
    list->Add(hIncl);
    list->Add(hLead);
    list->Add(hAway);
  }

  cout << TString::Format("Name: %s", setting.c_str()).Data() << endl;
  list->Write(setting.c_str(),1);
  delete list;
}

TH2F *make_hists(string name, string title, string obs, string secondary){
  TH2F *hist = new TH2F(name.c_str(), title.c_str(),
                        100, 0, 1.0, 200, 0, 200); // Assuming frag, pt
  hist->Sumw2();
  hist->GetYaxis()->SetTitle(TString::Format("%s", secondary.c_str()).Data()); //"pt");
  hist->GetXaxis()->SetTitle(TString::Format("%s", obs.c_str()).Data());

  if (obs == "orth"){
    if (secondary == "nconst")
      hist->SetBins(100, 0, 50, 100, 0, 100);
    else
      hist->SetBins(100, 0, 50, 200, 0, 200);
  }
  else if (secondary == "nconst"){
    hist->SetBins(100, 0, 1.0, 100, 0, 100);
  }
  return hist;
}

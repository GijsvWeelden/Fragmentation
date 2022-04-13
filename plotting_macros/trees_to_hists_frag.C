
#include <vector>
#include <iostream>
#include <time.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "ROOT/RDataFrame.hxx"

//-------------------------------------------------------------
//
// Reads in tree with jet data. Saves the variables in 2D hists of type pt:observable
// Distinguishes between leading and away-side jets and also considers the complete jet sample
// Note that this requires all data be combined into one file with the hadd command
// See the notes for more details
//
//-------------------------------------------------------------

void read_chain(TChain *chain, string setting, vector<string> obs);
void save_hists(TTree *tree, string setting, vector<string> obs, TList* L);//TFile *outFile);
void hadron_frag(TTree *tree, string setting, string jetType, TList* L); //TFile *outFile);
TH2F *make_hists(string name, string title, string obs, string secondary);

void trees_to_hists_frag(void){
  gROOT->SetBatch();
  double time = clock();
  double t2, dt;
  string jetType = "full"; //"charged";
  string suffix = jetType + "_nobkg"; // charged/full, nobkg/<nothing>
  string sNN = "5tev02"; // e.g. 5tev02
  string settings = "ppAAnrAAr";
  std::vector<string> observables = {"frag"}; //{"frag", "orth"};

  string outName = "2dhist_test"; //"2dhists_frag_" + sNN + "_" + jetType; //suffix;
  TFile *outFile = new TFile(Form("%s.root", outName.c_str()),"RECREATE");
  if (outFile) std::cout << "Output file: " << outName << ".root" << std::endl;
  else{
    std::cout << "ERROR: Output file not created. Aborting" << std::endl;
    return;
  }

  TFile* ppFile = TFile::Open(TString::Format("../run_pp_%s/jet_frag_%s_nobkg.root", sNN.c_str(), jetType.c_str()).Data(), "read");
  if (ppFile){
    TTree *ppT = (TTree*)ppFile->Get("jetprops");
    outFile->cd();
    TList *L = new TList();
    L->SetName("pp");//setting.c_str());
    L->SetOwner(true);
    //hadron_frag(ppT, "pp", jetType, L);
    //save_hists(ppT, "pp", observables, L);
    L->Write("pp",1);
    delete L;
    delete ppT;
    ppFile->Close();
  }
  dt = (clock() - time)/CLOCKS_PER_SEC;
  std::cout << "Time taken for pp: " << dt << std::endl;
  t2 = clock();

  //TFile nrFile(TString::Format("../run_AA_%s_norecoil/jet_frag_%s_nobkg.root", sNN.c_str(), jetType.c_str()).Data(),"read");
  TFile* nrFile = TFile::Open(TString::Format("../run_AA_%s_norecoil/jet_frag_%s_nobkg.root", sNN.c_str(), jetType.c_str()).Data(), "read");
  if (nrFile){
    TTree *nrT = (TTree*)nrFile->Get("jetprops");
    outFile->cd();
    TList* L = new TList();
    L->SetName("AAnr");
    L->SetOwner(true);
    //hadron_frag(nrT, "AAnr", jetType, L);
    //save_hists(nrT, "AAnr", observables, L);
    delete L;
    delete nrT;
    nrFile->Close();
  }
  dt = (clock() - t2)/CLOCKS_PER_SEC;
  std::cout << "Time taken for NR: " << dt << std::endl;
  t2 = clock();

  //TFile rNoBkgFile(TString::Format("../run_AA_%s_recoil/jet_frag_%s_nobkg.root", sNN.c_str(), jetType.c_str()).Data(),"read");
  TFile* rNoBkgFile = TFile::Open(TString::Format("../run_AA_%srecoil/jet_frag_%s_nobkg.root", sNN.c_str(), jetType.c_str()).Data(), "read");
  if (rNoBkgfile){
    TTree *rNoBkgT = (TTree*)rNoBkgFile.Get("jetprops");
    outFile->cd();
    TList* L = new TList();
    L->SetName("AAr_nobkg");
    L->SetOwner(true);
    //hadron_frag(rNoBkgT, "AAr_nobkg", jetType, L);
    //save_hists(rNoBkgT, "AAr_nobkg", observables, L);
    delete L;
    delete rNoBkgT;
    rNoBkgFile->Close();
  }
  dt = (clock() - t2)/CLOCKS_PER_SEC;
  std::cout << "Time taken for R_nobkg: " << dt << std::endl;
  t2 = clock();

  //TFile rBkgFile(TString::Format("../run_AA_%s_recoil/jet_frag_%s.root", sNN.c_str(), jetType.c_str()).Data(),"read");
  TFile* rBkgFile = TFile::Open(TString::Format("../run_AA_%srecoil/jet_frag_%s.root", sNN.c_str(), jetType.c_str()).Data(), "read");
  if (rBkgFile){
    TTree *rBkgT = (TTree*)rBkgFile.Get("jetprops");
    outFile->cd();
    TList* L = new TList();
    L->SetName("AAr_bkg");
    L->SetOwner(true);
    //hadron_frag(rBkgT, "AAr_bkg", jetType, L);
    //save_hists(rBkgT, "AAr_bkg", observables, L);
    delete L;
    delete rBkgT;
    rBkgFile->Close();
  }
  dt = (clock() - t2)/CLOCKS_PER_SEC;
  std::cout << "Time taken for R_bkg: " << dt << std::endl;

  time = (clock() - time)/CLOCKS_PER_SEC;
  std::cout << "Total time taken: " << time << "seconds" << std::endl;
  //outFile->Close();
}

void save_hists(TTree *tree, string setting, vector<string> obs, TList* L){//TFile *outFile){
  std::cout << "Save hists " << setting << std::endl;
  for (auto iobs=0; iobs<obs.size(); iobs++){
    TH2F *hIncl = make_hists(TString::Format("%s_%s",
          setting.c_str(), obs[iobs].c_str()).Data(),
        TString::Format("Jet %s (%s)",
          obs[iobs].c_str(), setting.c_str()).Data(),
        obs[iobs].c_str(),
        "pt"
        );
    TH2F *hLead = make_hists(TString::Format("%s_%s_leading",
          setting.c_str(), obs[iobs].c_str()).Data(),
        TString::Format("Leading jet %s (%s)",
          obs[iobs].c_str(), setting.c_str()).Data(),
        obs[iobs].c_str(),
        "pt"
        );
    TH2F *hAway = make_hists(TString::Format("%s_%s_awayside",
          setting.c_str(), obs[iobs].c_str()).Data(),
        TString::Format("Awayside jet %s (%s)",
          obs[iobs].c_str(), setting.c_str()).Data(),
        obs[iobs].c_str(),
        "pt"
        );
    // Draw hIncl
    tree->Draw(TString::Format("pt:%s>>%s_%s",
          obs[iobs].c_str(),
          setting.c_str(),
          obs[iobs].c_str()
          ).Data(),
        "evwt"
        );
    // Draw hLead
    tree->Draw(TString::Format("pt:%s>>%s_%s_leading",
          obs[iobs].c_str(),
          setting.c_str(),
          obs[iobs].c_str()
          ).Data(),
        "(evwt*(ijet==0))"
        );
    // Draw hAway
    tree->Draw(TString::Format("pt:%s>>%s_%s_awayside",
          obs[iobs].c_str(),
          setting.c_str(),
          obs[iobs].c_str()
          ).Data(),
        TString::Format("(evwt*(ijet>0 && abs(dphi)>%f))",
          TMath::PiOver2()).Data()
        );
    L->Add(hIncl);
    L->Add(hLead);
    L->Add(hAway);
  }
  for (auto iobs=0; iobs<obs.size(); iobs++){
    TH2F *hIncl = make_hists(TString::Format("%s_%s",
          setting.c_str(), obs[iobs].c_str()).Data(),
        TString::Format("Jet %s (%s)",
          obs[iobs].c_str(), setting.c_str()).Data(),
        obs[iobs].c_str(),
        "nconst"
        );
    TH2F *hLead = make_hists(TString::Format("%s_%s_leading",
          setting.c_str(), obs[iobs].c_str()).Data(),
        TString::Format("Leading jet %s (%s)",
          obs[iobs].c_str(), setting.c_str()).Data(),
        obs[iobs].c_str(),
        "nconst"
        );
    TH2F *hAway = make_hists(TString::Format("%s_%s_awayside",
          setting.c_str(), obs[iobs].c_str()).Data(),
        TString::Format("Awayside jet %s (%s)",
          obs[iobs].c_str(), setting.c_str()).Data(),
        obs[iobs].c_str(),
        "nconst"
        );
    // Draw hIncl
    tree->Draw(TString::Format("nconst:%s>>%s_%s",
          obs[iobs].c_str(),
          setting.c_str(),
          obs[iobs].c_str()
          ).Data(),
        "evwt"
        );
    // Draw hLead
    tree->Draw(TString::Format("nconst:%s>>%s_%s_leading",
          obs[iobs].c_str(),
          setting.c_str(),
          obs[iobs].c_str()
          ).Data(),
        "(evwt*(ijet==0))"
        );
    // Draw hAway
    tree->Draw(TString::Format("nconst:%s>>%s_%s_awayside",
          obs[iobs].c_str(),
          setting.c_str(),
          obs[iobs].c_str()
          ).Data(),
        TString::Format("(evwt*(ijet>0 && abs(dphi)>%f))",
          TMath::PiOver2()).Data()
        );
    L->Add(hIncl);
    L->Add(hLead);
    L->Add(hAway);
  }
  std::cout << TString::Format("Name: %s", setting.c_str()).Data() << std::endl;
}

TH2F *make_hists(string name, string title, string obs, string secondary){
  TH2F *hist = new TH2F(name.c_str(), title.c_str(), 100, 0, 1.0, 200, 0, 200); // Assuming frag, pt
  hist->Sumw2();
  hist->GetYaxis()->SetTitle(TString::Format("%s", secondary.c_str()).Data());
  hist->GetXaxis()->SetTitle(TString::Format("%s", obs.c_str()).Data());
  std::cout << "After making the hist." << std::endl;

  const double x[48] = {1e-5,
    2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4, 1e-3,
    2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3, 1e-2,
    2e-2, 3e-2, 4e-2, 5e-2, 6e-2, 7e-2, 8e-2, 9e-2, 1e-1,
    2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1, 1,
    2, 3, 4, 5, 6, 7, 8, 9, 10,
    20, 30};
  double y1[101];
  double y2[201];
  for (int i = 0; i <= 101; ++i){
    y1[i] = i;
    y2[i] = i;
    if (i!=0) y2[i+100] = i+100;
  }

  if (obs == "orth"){
    if (secondary == "nconst")
      //hist->SetBins(26, *x, 100, 0, 100);
      hist->SetBins(47, x, 100, y1);
    else
      //hist->SetBins(26, *x, 200, 0, 200);
      hist->SetBins(47, x, 200, y2);
  }
  else if (secondary == "nconst"){
    hist->SetBins(100, 0, 1.0, 100, 0, 100);
  }
  std::cout << "After rebinning." << std::endl;
  return hist;
}

void hadron_frag(TTree *tree, string setting, string jetType, TList* L){//TFile *outFile){
  std::vector<int> chPDG = {211, 321, 2212};
  std::vector<string> chHadrons = {"pi", "K", "p"};
  std::vector<int> nPDG = {111, 130, 310, 311, 3122};
  std::vector<string> nHadrons = {"pi0", "K0L", "K0S", "K0", "Lambda0"};

  std::cout << "Hadron spectra: " << setting << std::endl;
  for (int i = 0; i < chHadrons.size(); i++){
    string hadron = chHadrons[i];
    int code = chPDG[i];
    std::cout << hadron << " = " << code << std::endl;
    TH2F *hFrag = make_hists(TString::Format("%s_%s_frag", setting.c_str(), hadron.c_str()).Data(),
        TString::Format("%s fragmentation (%s)", hadron.c_str(), setting.c_str()).Data(),
        "frag", "pt");
    //TH2F *orth = make_hists(TString::Format("%s_%s_orth", setting.c_str(), hadron.c_str()).Data(),
    //                        TString::Format("%s #it{j}_{T} (%s)", hadron.c_str(), setting.c_str()).Data(),
    //                        "orth", "pt");
    tree->Draw(TString::Format("pt:frag>>%s_%s_frag", setting.c_str(), hadron.c_str()).Data(),
        TString::Format("evwt*(abs(pdg) == %d)", code).Data());
    L->Add(hFrag);
  }
  if (jetType == "charged") return;
  for (int i = 0; i < nHadrons.size(); i++){
    string hadron = nHadrons[i];
    int code = nPDG[i];
    std::cout << hadron << " = " << code << std::endl;
    TH2F *hFrag = make_hists(TString::Format("%s_%s_frag", setting.c_str(), hadron.c_str()).Data(),
        TString::Format("%s fragmentation (%s)", hadron.c_str(), setting.c_str()).Data(),
        "frag", "pt");
    //TH2F *orth = make_hists(TString::Format("%s_%s_orth", setting.c_str(), hadron.c_str()).Data(),
    //                        TString::Format("%s #it{j}_{T} (%s)", hadron.c_str(), setting.c_str()).Data(),
    //                        "orth", "pt");
    tree->Draw(TString::Format("pt:frag>>%s_%s_frag", setting.c_str(), hadron.c_str()).Data(),
        TString::Format("evwt*(abs(pdg) == %d)", code).Data());
    L->Add(hFrag);
  }
}


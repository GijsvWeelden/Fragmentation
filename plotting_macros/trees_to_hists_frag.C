
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
// Note that this requires all data be combined into one file with the hadd command
// See the notes for more details
//
//-------------------------------------------------------------

void read_chain(TChain *chain, string setting, vector<string> obs);
void save_hists(TTree *tree, string setting, vector<string> obs, TFile *outFile);
TH2F *make_hists(string name, string title, string obs, string secondary);

void trees_to_hists_frag(void){
  gROOT->SetBatch();
  double time = clock();
  double t2, dt;
  string jetType = "charged";
  string suffix = jetType + "_nobkg"; // charged/full, nobkg/<nothing>
  string sNN = "5tev02"; // e.g. 5tev02
  string settings = "ppAAnrAAr";
  std::vector<string> observables = {"frag", "orth"};

  string outName = "2dhist_test"; //"2dhists_frag_" + sNN + "_" + jetType; //suffix;
  TFile *outFile = new TFile(Form("%s.root", outName.c_str()),"RECREATE");
  if (outFile) std::cout << "Output file: " << outName << ".root" << std::endl;
  else{
    std::cout << "ERROR: Output file not created. Aborting" << std::endl;
    return;
  }

  TFile ppFile(TString::Format("../run_pp_%s/jet_frag_%s_nobkg.root", sNN.c_str(), jetType.c_str()).Data(),"read");
  TTree *ppT = (TTree*)ppFile.Get("jetprops");
  save_hists(ppT, "pp", observables, outFile);
  delete ppT;
  ppFile.Close();
  dt = (clock() - time)/CLOCKS_PER_SEC;
  std::cout << "Time taken for pp: " << dt << std::endl;
  t2 = clock();

  TFile nrFile(TString::Format("../run_AA_%s_norecoil/jet_frag_%s_nobkg.root", sNN.c_str(), jetType.c_str()).Data(),"read");
  TTree *nrT = (TTree*)nrFile.Get("jetprops");
  save_hists(nrT, "AAnr", observables, outFile);
  delete nrT;
  nrFile.Close();
  dt = (clock() - t2)/CLOCKS_PER_SEC;
  std::cout << "Time taken for NR: " << dt << std::endl;
  t2 = clock();

  TFile rNoBkgFile(TString::Format("../run_AA_%s_recoil/jet_frag_%s_nobkg.root", sNN.c_str(), jetType.c_str()).Data(),"read");
  TTree *rNoBkgT = (TTree*)rNoBkgFile.Get("jetprops");
  save_hists(rNoBkgT, "AAr_nobkg", observables, outFile);
  delete rNoBkgT;
  rNoBkgFile.Close();
  dt = (clock() - t2)/CLOCKS_PER_SEC;
  std::cout << "Time taken for R_nobkg: " << dt << std::endl;
  t2 = clock();

  TFile rBkgFile(TString::Format("../run_AA_%s_recoil/jet_frag_%s.root", sNN.c_str(), jetType.c_str()).Data(),"read");
  TTree *rBkgT = (TTree*)rBkgFile.Get("jetprops");
  save_hists(rBkgT, "AAr_bkg", observables, outFile);
  delete rBkgT;
  rBkgFile.Close();
  dt = (clock() - t2)/CLOCKS_PER_SEC;
  std::cout << "Time taken for R_bkg: " << dt << std::endl;

  time = (clock() - time)/CLOCKS_PER_SEC;
  std::cout << "Total time taken: " << time << "seconds" << std::endl;
  //outFile->Close();
}

void save_hists(TTree *tree, string setting, vector<string> obs, TFile *outFile){
  outFile->cd();
  TList *list = new TList();
  list->SetName(setting.c_str());
  list->SetOwner(true);
  std::cout << "made list" << std::endl;
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
    std::cout << "before draw" << std::endl;

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
    std::cout << "before list add" << std::endl;
    list->Add(hIncl);
    list->Add(hLead);
    list->Add(hAway);
  }
  std::cout << "before nconst loop" << std::endl;
  for (auto iobs=0; iobs<obs.size(); iobs++){
    std::cout << "more hists" << std::endl;
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
    list->Add(hIncl);
    list->Add(hLead);
    list->Add(hAway);
  }

  std::cout << TString::Format("Name: %s", setting.c_str()).Data() << std::endl;
  list->Write(setting.c_str(),1);
  delete list;
}

TH2F *make_hists(string name, string title, string obs, string secondary){
  TH2F *hist = new TH2F(name.c_str(), title.c_str(), 100, 0, 1.0, 200, 0, 200); // Assuming frag, pt
  hist->Sumw2();
  hist->GetYaxis()->SetTitle(TString::Format("%s", secondary.c_str()).Data());
  hist->GetXaxis()->SetTitle(TString::Format("%s", obs.c_str()).Data());

  if (obs == "orth"){
    if (secondary == "nconst")
      hist->SetBins(100, 0, 25, 100, 0, 100);
    else
      hist->SetBins(100, 0, 25, 200, 0, 200);
  }
  else if (secondary == "nconst"){
    hist->SetBins(100, 0, 1.0, 100, 0, 100);
  }
  return hist;
}

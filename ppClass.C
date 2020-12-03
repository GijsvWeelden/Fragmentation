#define ppClass_cxx
#include "ppClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void ppClass::Loop()
{
  //   In a ROOT session, you can do:
  //      root> .L ppClass.C
  //      root> ppClass t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  float_t ptL = 0;
  // ptCuts=1: high pt, ptL>120 GeV, ptS>30 GeV
  // ptCuts=0: low pt, ptL>35 GeV, ptS>10 GeV
  bool ptCuts = 0;

  if(ptCuts){
    TH1F *hdphi = new TH1F("Hist_dphi_highpt","Acoplanarity dphi: ptL> 120 GeV, ptS> 30 GeV, |eta_jet|<2",100,0,3.2);
    TH1F *hdphiw = new TH1F("Hist_dphi_highpt","Acoplanarity dphi (weighted): ptL> 120 GeV, ptS> 30 GeV, |eta_jet|<2",100,0,3.2);
  }
  else{
    TH1F *hdphi = new TH1F("Hist_dphi_lowpt","Acoplanarity dphi: ptL> 35 GeV, ptS> 10 GeV, |eta_jet|<2",100,0,3.2);
    TH1F *hdphiw = new TH1F("Hist_dphi_lowpt","Acoplanarity dphi (weighted): ptL> 35 GeV, ptS> 10 GeV, |eta_jet|<2",100,0,3.2);
  }

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    /* If this block outputs all zeroes, restart the ssh connection  */
    /*
       if (ientry < 10){
       cout << "jentry:" << jentry << "\nientry:" << ientry << "\nievt/ijet:" << ievt << "/" << ijet;
       cout << "\npt/phi:" << pt << "/" << phi <<  "\ndphi:" << dphi;
       cout << "\n";
       }
     */
    // Identify leading jet
    if(ijet==0){
      ptL = pt;
      //cout << "Leading jet ievt(pt): " << ievt << "(" << pt << ")\n";
      continue;
    }
    // Kinematic cuts on subleading jets
    /*
       else if(ptL>120 && pt>30){
    //cout << "Acoplanarity: " << abs(dphi) << "\n";
    hdphi->Fill(abs(dphi));
    //hdphiw->Fill(abs(dphi),evwt);
    }
    else if(ptL>35 && pt>10){
    //hdphi->Fill(abs(dphi));
    hdphiw->Fill(abs(dphi),evwt);
    }

     */
    else{
      if(ptCuts && ptL>120 && pt>30){
        hdphi->Fill(abs(dphi));
        hdphiw->Fill(abs(dphi),evwt);
      }
      else if(ptL>35 && pt>10){
        hdphi->Fill(abs(dphi));
        hdphiw->Fill(abs(dphi),evwt);
      }
    }
  }
  hdphi->Draw();
  // TODO: write the histogram to its own .root file
  string houtname = "pp2tev76_hdphi";
  string hwoutname = houtname + "w";

  if(ptCuts){
    houtname += "high_pt";
    hwoutname += "high_pt";
    else{
      houtname += "low_pt";
      hwoutname+= "low_pt";
      outputname += ".root";

      hdphi->Write(houtname,TObject::kOverwrite);
      hdphiw->Write(hwoutname,TObject::kOverwrite);
    }

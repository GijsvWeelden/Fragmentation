#define ppClass_cxx
#include "ppClass.h"
#include <TH1.h>
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
  Long64_t nbytes = 0, nb = 0;
  float_t ptL = 0;
  // ptCuts=1: high pt, ptL>120 GeV, ptS>30 GeV
  // ptCuts=0: low pt, ptL>35 GeV, ptS>10 GeV
  bool ptCuts = 0;
  /*
     string htitle = "Acoplanarity dphi";
     string hwtitle = "Acoplanarity dphi (weighted):";

     if(ptCuts){
     htitle += "ptL>120 GeV, ptS>30 GeV";
     hwtitle += "ptL>120 GeV, ptS>30 GeV";
     }
     else{
     htitle += "ptL>35 GeV, ptS>10 GeV";
     hwtitle += "ptL>35 GeV, ptS>10 GeV"; 
     }   

     htitle += ", |eta_jet|<2";
     hwtitle += ", |eta_jet|<2";
     string title1 = &htitle;
     string title2 = &hwtitle;
   */
  TH1F *hdphi = new TH1F("Hist_dphi","a",100,0,3.2);
  TH1F *hwdphi = new TH1F("Hist_w_dphi","b",100,0,3.2);
  if(ptCuts){
    hdphi->SetTitle("Acoplanarity dphi: ptL> 120 GeV, ptS> 30 GeV, |eta_jet|<2");
    hwdphi->SetTitle("Acoplanarity dphi (weighted): ptL> 120 GeV, ptS> 30 GeV, |eta_jet|<2");
  }
  else{
    hdphi->SetTitle("Acoplanarity dphi: ptL> 35 GeV, ptS> 10 GeV, |eta_jet|<2");
    hwdphi->SetTitle("Acoplanarity dphi (weighted): ptL> 35 GeV, ptS> 10 GeV, |eta_jet|<2");
  }
  //TH1F *hdphi = new TH1F("Hist_dphi_lowpt","Acoplanarity dphi: ptL> 35 GeV, ptS> 10 GeV, |eta_jet|<2",100,0,3.2);
  //TH1F *hwdphi = new TH1F("Hist_dphi_lowpt","Acoplanarity dphi (weighted): ptL> 35 GeV, ptS> 10 GeV, |eta_jet|<2",100,0,3.2);


  /*
     if(ptCuts){
     TH1F *hdphi = new TH1F("Hist_dphi_highpt","Acoplanarity dphi: ptL> 120 GeV, ptS> 30 GeV, |eta_jet|<2",100,0,3.2);
     TH1F *hwdphi = new TH1F("Hist_dphi_highpt","Acoplanarity dphi (weighted): ptL> 120 GeV, ptS> 30 GeV, |eta_jet|<2",100,0,3.2);
     }
     else{
     TH1F *hdphi = new TH1F("Hist_dphi_lowpt","Acoplanarity dphi: ptL> 35 GeV, ptS> 10 GeV, |eta_jet|<2",100,0,3.2);
     TH1F *hwdphi = new TH1F("Hist_dphi_lowpt","Acoplanarity dphi (weighted): ptL> 35 GeV, ptS> 10 GeV, |eta_jet|<2",100,0,3.2);
     }
   */

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
    //hwdphi->Fill(abs(dphi),evwt);
    }
    else if(ptL>35 && pt>10){
    //hdphi->Fill(abs(dphi));
    hwdphi->Fill(abs(dphi),evwt);
    }

     */
    else{
      if(ptCuts && ptL>120 && pt>30){
        hdphi->Fill(abs(dphi));
        hwdphi->Fill(abs(dphi),evwt);
      }
      else if(ptL>35 && pt>10){
        hdphi->Fill(abs(dphi));
        hwdphi->Fill(abs(dphi),evwt);
      }
    }
  }
  // TODO: add canvas
  TCanvas *c1 = new TCanvas();
  hdphi->Draw();
  TCanvas *c2 = new TCanvas();
  hwdphi->Draw();
  // TODO: write the histogram to its own .root file
  /*
     string houtname = "pp2tev76_hdphi";
     string hwoutname = houtname + "w";

     if(ptCuts){
     houtname += "high_pt";
     hwoutname += "high_pt";
     }
     else{
     houtname += "low_pt";
     hwoutname+= "low_pt";
     }
     houtname += ".root";
     hwoutname += ".root";
   */

  if(ptCuts){
    TFile hdphi_save("pp2tev76_hdphi_highpt.root","RECREATE");
    hdphi->Write();//hdphi_save,TObject::kOverwrite); 
    TFile hwdphi_save("pp2tev76_hwdphi_highpt.root","RECREATE");
    hwdphi->Write();//hwdphi,TObject::kOverwrite);
  }
  else{
    TFile hdphi_save("pp2tev76_hdphi_lowpt.root","RECREATE");
    hdphi->Write();//hdphi_save,TObject::kOverwrite);
    TFile hwdphi_save("pp2tev76_hwdphi_lowpt.root","RECREATE");
    hwdphi->Write();//hwdphi_save,TObject::kOverwrite);
  }
}

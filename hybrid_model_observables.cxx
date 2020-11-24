/*
 *  Program to analyse root files created by the "analyze..." scripts.
 *
 *  Observables:
 *  - Acoplanarity (Delta_phi)
 *  - Differential jet shape (rho)
 *  
 *  Author: Gijs van Weelden, Nikhef
 *
 */

#include "TPDGCode.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "THn.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TCanvas.h"

#include <iostream>
#include <string>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
#include "getopt.h"

// Acoplanarity
float Delta_phi (float jet1_phi, float jet2_phi){
  float Delta = jet1_phi - jet2_phi;
  Delta = abs(Delta);
  float pi = 3.14159;
  if(Delta > pi){
    Delta -= pi;
  }
  return Delta;
}

/*
* float rho (){
*   for (ijet=0; ijet<Njets; ijet++;){
*     // Sort tracks
*     for (itrack=0; itrack<Ntrack; itrack++){
*       sumtracks += pttrack
*     }
*   }
*   rho = / Njets / delta_r
* }
*/

int main (int argc, char** argv){// string &ppFile, string &AAFile){ //, string &outputFile){
  // Read in ppfile
  auto ppLoad = TFile::Open("jet_tree_pp2tev76_full_nobkg.root");
  if (!ppLoad || ppLoad->IsZombie()){
    return 1;
  }
  
  auto ppHist = new TH1I();

  // Read in the tree
  TTreeReader ppReader("jetprops",ppLoad);
  
  TTreeReaderValue<int> pp_evt(ppReader, "ievt");
  TTreeReaderValue<int> pp_jet(ppReader, "ijet");
  TTreeReaderValue<float> pp_phi(ppReader, "phi");
  
  // Check if there are 2 jets in a single event
  while(ppReader.Next()){
    ppHist->Fill(*pp_evt + *pp_phi);
  }

  // Save histogram as .root and .pdf file
  TFile* ppSave = new TFile("ppHist.root","RECREATE");
  ppHist->Write();
  ppSave->Close();
  delete ppSave;
  TCanvas *ppCanvas = new TCanvas();
  ppHist->Draw();
  ppCanvas->Print("ppHist.pdf","PDF");

  // ROOT->GetListOfCanvases()->Draw();
  // TFile *ppFile = new TFile("jet_tree_pp_2tev76_full_nobkg.root");
  // TTree *ppTree = new (TTree), ppFile.Get("jetprops");
  // Read in AAfile
  // TFile *AAFile = new TFile("jet_tree_AA_2tev76_norecoil_full_nobkg.root");
  // TTree *AATree = new (TTree), AAFile.Get("jetprops");

  cout << "Loaded files and trees\n";

  // ppFile->Close();
  // AAFile->Close();

  cout << "Closed files\n";
// calculate observables
// write/plot to output file
// Save output
// Close output
  return 0;
}


//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 30 16:14:26 2020 by ROOT version 6.16/00
// from TTree jetprops/Jet properties
// found on file: jet_tree_pp2tev76_full_nobkg.root
//////////////////////////////////////////////////////////

#ifndef ppClass_h
#define ppClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ppClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

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

   ppClass(TTree *tree=0);
   virtual ~ppClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ppClass_cxx
ppClass::ppClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("jet_tree_pp2tev76_full_nobkg.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("jet_tree_pp2tev76_full_nobkg.root");
      }
      f->GetObject("jetprops",tree);

   }
   Init(tree);
}

ppClass::~ppClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ppClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ppClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ppClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ievt", &ievt, &b_ievt);
   fChain->SetBranchAddress("ijet", &ijet, &b_ijet);
   fChain->SetBranchAddress("evwt", &evwt, &b_evwt);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("dphi", &dphi, &b_dphi);
   fChain->SetBranchAddress("nconst", &nconst, &b_nconst);
   fChain->SetBranchAddress("zg", &zg, &b_zg);
   fChain->SetBranchAddress("Rg", &Rg, &b_Rg);
   fChain->SetBranchAddress("nSD", &nSD, &b_nSD);
   fChain->SetBranchAddress("mass", &mass, &b_mass);
   fChain->SetBranchAddress("mz2", &mz2, &b_mz2);
   fChain->SetBranchAddress("mr", &mr, &b_mr);
   fChain->SetBranchAddress("mr2", &mr2, &b_mr2);
   fChain->SetBranchAddress("rz", &rz, &b_rz);
   fChain->SetBranchAddress("r2z", &r2z, &b_r2z);
   Notify();
}

Bool_t ppClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ppClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ppClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ppClass_cxx

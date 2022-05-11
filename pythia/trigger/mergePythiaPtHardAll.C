
void mergePythiaPtHardAll(){
  // This is used to merge the analysis-output from different 
  // data samples/pt_hard bins


  // Based on: $ALICE_ROOT/PWGJE/runJetSpectrumUnfolding

  // Updated by: Marco van Leeuwen, Nikhef

  Int_t debug = 1;

  //Char_t *tag = "charged_jet_spectra_13000GeV";
  Char_t *tag = "full_jet_spectra_5020GeV";

  const char *cFile[] = {
    Form("output/pythia_%s_5_15_all.root",tag),
    Form("output/pythia_%s_15_25_all.root",tag),
    Form("output/pythia_%s_25_40_all.root",tag),
    //Form("output/pythia_%s_40_100_all.root",tag)};
    Form("output/pythia_%s_40_80_all.root",tag),
    Form("output/pythia_%s_80_150_all.root",tag),
    Form("output/pythia_%s_150_250_all.root",tag),
    Form("output/pythia_%s_250_400_all.root",tag)};

  TFile *fOut = new TFile(Form("pythia_%s_merged.root",tag),"RECREATE");

  const Int_t nBins = sizeof(cFile)/sizeof(cFile[0]);

  TObjArray *outList = new TObjArray();


  for(int ib = 0;ib < nBins;++ib){
    TFile *fIn = TFile::Open(cFile[ib]);
    //TH1* hTrials = (TH1F*)fIn->FindObject("fh1Trials");
    //TProfile* hXsec = (TProfile*)fIn->FindObject("fh1Xsec");
    TH1* hTrials = (TH1F*)fIn->Get("hNEvent");
    TProfile* hXsec = (TProfile*)fIn->Get("hXSec");
    if (hTrials == 0 || hXsec == 0) {
      cout << "WARNING: no nTrials or Xsec histo; skipping" << endl;
      break;
    }
  
    Double_t xsection = hXsec->GetBinContent(1);
    Double_t nTrials = hTrials->Integral();
    Double_t scaleFactor = xsection/nTrials;

    keyList = fIn->GetListOfKeys();
    TIter keyIter(keyList); 
	
    // only merge top level dir, no subdir structure traversal
    while (key = (TKey*) keyIter()) {
	
      
      obj = key->ReadObj();
	if (debug >= 2)
	  cout << "Adding object (ib " << ib << " ) " << obj->GetName() << " type " << obj->Class()->GetName() << endl;

	if(obj->InheritsFrom("TH1")){
	  if (debug >= 3)
	    cout << "Case TH1" << endl;
	  TH1 *h1 = (TH1*)obj;
	  if(ib==0){
	    h1Add = (TH1*)h1->Clone(h1->GetName());
	    h1Add->Sumw2();
	    h1Add->Scale(scaleFactor);
	    h1Add->SetDirectory(fOut);
	    outList->Add(h1Add);
	  }
	  else{
	    h1Add = (TH1*) outList->FindObject(h1->GetName());
	    if (h1Add)
	      h1Add->Add(h1,scaleFactor);
	    else 
	      cout << "Could not find histo: " << h1->GetName() << endl;
	  }
	}
	else if(obj->InheritsFrom("THnSparse")){
	  if (debug >= 3)
	    cout << "Case THnSParse" << endl;
	  THnSparse *hn = (THnSparse*)obj;
	  if(ib==0){
	    hnAdd = (THnSparse*)hn->Clone(hn->GetName());
	    hnAdd->Sumw2();
	    hnAdd->Scale(scaleFactor);
	    outList->Add(hnAdd);
	  }
	  else{
	    hnAdd = (THnSparse*) outList->FindObject(hn->GetName());
	    if (hnAdd)
	      hnAdd->Add(hn,scaleFactor);
	    else 
	      cout << "Could not find (sparse) histo: " << hn->GetName() << endl;
	  }
	}
	else{
	  if (ib==0) 
	    cout << "Skipping " << obj->GetName() << " -- not a TH1 or THnSparse" << endl;
	}
      
    }
  }
  fOut->cd();
  
  if (debug >= 1)
    cout << "outList has " << outList->GetEntries() << " entries" << endl;
  if (debug >= 2)
    outList->Print();
  outList->Write();
  outList->Delete();
  delete outList;

  if (debug >= 1) 
    cout << "Closing output file" << endl;
  fOut->cd();
  //lOut->Write(lOut->GetName(),TObject::kSingleKey);
  fOut->Close();
}


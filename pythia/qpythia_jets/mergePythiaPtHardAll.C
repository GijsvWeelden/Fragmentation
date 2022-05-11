void removeOutliers(TH1* h1, Int_t minBinValue = 2) {
  // set bins with less than minBinValue entries to 0
  Int_t nbins = h1->GetNbinsX()+2;
  TH2* h2 = dynamic_cast<TH2*>(h1);
  if (h2 && h2->GetNbinsY() != 0) 
    nbins *= h2->GetNbinsY()+2;
  TH3* h3 = dynamic_cast<TH3*>(h1);
  if (h3 && h3->GetNbinsZ() != 0) 
    nbins *= h3->GetNbinsZ()+2;
  cout << "hist " << h1->GetName() << " nbins " << nbins << endl;
  for (Int_t ibin = 0; ibin < nbins; ibin++) {
    if (h1->GetBinContent(ibin) < minBinValue && h1->GetBinContent(ibin) > 0) {
      //cout << "Histo: " << h1->GetName() << " setting bin " << ibin << " to zero"<< endl;
      h1->SetBinContent(ibin,0);
    }
  }
}

void mergePythiaPtHardAll(){
  // This is used to merge the analysis-output from different 
  // data samples/pt_hard bins


  // Based on: $ALICE_ROOT/PWGJE/runJetSpectrumUnfolding

  // Updated by: Marco van Leeuwen, Nikhef

  Int_t debug = 1;

  Float_t pthard[] = {5, 15, 25, 40, 80, 120, 250};
  Float_t qhat = 0; //50;

  //Char_t *fin_name_format = "output/pythia_jet_shapes_sd_%.0f_%.0f_%.0f.root";
  //Char_t *fout_name = Form("pythia_jet_shapes_sd_%.0f_merged.root",qhat);
  //Char_t *fin_name_format = "output/pythia_jets_%.0f_%.0f_%.0f_full.root";
  //Char_t *fout_name = Form("qpythia_jets_%.0f_full_merged.root",qhat);
  //Char_t *fin_name_format = "output/pythia_jets_%.0f_%.0f_%.0f_charged.root";
  //Char_t *fout_name = Form("qpythia_jets_%.0f_charged_merged.root",qhat);

  // Jet shapes
  //Char_t *fin_name_format = "output_v11/pythia_jet_shapes_%.0f_%.0f_%.0f.root";
  //Char_t *fout_name = Form("qpythia_v11_jet_shapes_%.0f_merged.root",qhat);
  //Char_t *fin_name_format = "output/pythia_softdrop_shapes_fulljets_%.0f_%.0f_%.0f.root";
  //Char_t *fout_name = Form("qpythia_softdrop_shapes_fulljets_%.0f_merged.root",qhat);
  //Char_t *fin_name_format = "output/pythia_softdrop_shapes_chargedjets_%.0f_%.0f_%.0f.root";
  //Char_t *fout_name = Form("qpythia_softdrop_shapes_chargedjets_%.0f_merged.root",qhat);

  // Recoil shapes
  //Char_t *fin_name_format = "output_v11/pythia_recoil_shapes_%.0f_%.0f_%.0f.root";
  //Char_t *fout_name = Form("qpythia_v11_recoil_shapes_%.0f_merged.root",qhat);
  //Char_t *fin_name_format = "output/pythia_5tev02_recoil_shapes_fulljets_%.0f_%.0f_%.0f.root";
  //Char_t *fout_name = Form("qpythia_5tev02_recoil_shapes_fulljets_%.0f_merged.root",qhat);
 
  // Strangeness
  Char_t *fin_name_format = "output/pythia_5tev02_fulljets_strangeness_%.0f_%.0f_%.0f.root";
  Char_t *fout_name = Form("qpythia_5tev02_fulljets_strangeness_%.0f_merged.root",qhat);

  Int_t minBinValue = 0; // Binning too narrow for this procedure?

  //Char_t *fin_name_format = "output/pythia_jet_shapes_%.0f_%.0f_%.0f.root";
  //Char_t *fout_name = Form("pythia_jet_shapes_%.0f_merged.root",qhat);
  //Int_t minBinValue = 3;

  TFile *fOut = new TFile(fout_name,"RECREATE");

  const Int_t nBins = sizeof(pthard)/sizeof(pthard[0])-1;

  TObjArray *outList = new TObjArray();

  for(int ib = 0;ib < nBins;++ib){
    TFile *fIn = TFile::Open(Form(fin_name_format, pthard[ib],pthard[ib+1],qhat));
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
          if (minBinValue > 0)
            removeOutliers(h1, minBinValue);
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


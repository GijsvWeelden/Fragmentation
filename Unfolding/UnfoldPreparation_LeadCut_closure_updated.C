void UnfoldPreparation_EF_embed_LeadCut_updated(TString strIn="AnalysisResults.root",Double_t Rjet =0.1,Double_t PT_MIN = 0.0, Double_t DPT_MIN = -20.5,Double_t DPT_MAX=30.5,TString strL = "AliAnalysisTaskEmcalJetEnergyFlow_tracks_caloClusters_emcalCells_Embed_embed_Const_histos",Int_t binWidthPt=5, Int_t skipBins = 0){

//Open the file and get the list where the response matrix is stored in THnSparse format
TFile *f = TFile::Open(strIn.Data());
AliEmcalList *ali;
f->GetObject(strL.Data(),ali);

//Get the response matrix in THn format

THnSparseD* RM_thn = static_cast<THnSparseD*>(ali->FindObject(Form("ResponseMatrix_R%03d_0",int(Rjet*100))));

//THnSparseD* RM_miss = static_cast<THnSparseD*>(ali->FindObject(Form("MismatchResponseMatrix_R%03d_0",int(Rjet*100))));

//RM_thn->Add(RM_miss,-1.);

// Two dimensional unfolding (Dpt vs Pt on Gen and Det level)
const Int_t ndim =4;

Int_t nDim = RM_thn->GetNdimensions();
Int_t iPt_gen   = 0;
Int_t iPt_det  = 1;
Int_t iDpt_gen  = 2;
Int_t iDpt_det = 3;

//These are the original(before rebinning) projections of the THnSparse which encapsulates the respose. Maybe I don't need
//to project but instead I can just take the calculated 2D histograms
TH2D *Dpt_pt_true = dynamic_cast<TH2D*>(RM_thn->Projection(iDpt_gen,iPt_gen,"E"));
Dpt_pt_true->SetName(Form("DeltaPt_Pt_gen_R%03d",int(Rjet*100)));
TH2D *Dpt_pt_meas = dynamic_cast<TH2D*>(RM_thn->Projection(iDpt_det,iPt_det,"E"));
Dpt_pt_meas->SetName(Form("DeltaPt_Pt_det_R%03d",int(Rjet*100)));

//Setting bin configuration for Measured & True

Double_t pt_min = PT_MIN; // Lower Pt limit for detector level, Perhaps we can change this to evaluate systematics
Double_t pt_max = 100.0;  //Upper Pt limit for detector level (Gen level is higher by 10 GeV (feed-in effect)).
Double_t pt_min_gen = PT_MIN-5.;// Lower Pt limit for generator level
Double_t Dpt_min = DPT_MIN; //Lower Dpt limit for detector level
Double_t Dpt_max = DPT_MAX; //Upper Dpt limit for detector level(Gen level is higher by 5 GeV (feed-in effect))

      //Detector, Particle level
Int_t nPtBinWidth[2] = {5,5};
Double_t ptmin[2] = {pt_min,pt_min_gen};
Double_t ptmax[2] = {pt_max, pt_max+5};
Int_t nBinPt[2] = {int(ptmax[0]-ptmin[0])/nPtBinWidth[0],int(ptmax[1]-ptmin[1])/nPtBinWidth[1]};

     //Detector, Particle level
Int_t nDptBinWidth[2] = {1,1};
//Double_t Dptmin[2] = {Dpt_min ,Dpt_min};
Double_t Dptmin[2] = {Dpt_min ,-.5};
Double_t Dptmax[2] = {Dpt_max,Dpt_max+2};
Int_t nBinDpt[2] = {int(Dptmax[0]-Dptmin[0])/nDptBinWidth[0],int(Dptmax[1]-Dptmin[1])/nDptBinWidth[1]};

//Rebining the Detector level axes
TH2D *h2Dpt_Rebined_Meas = new TH2D(Form("h2Dpt_Rebined_Meas_R%03d",int(Rjet*100)),Form("#DeltaP_{t} distribution R =%.2f (Detector level -Rebinned);P_{t} (GeV/c);#DeltaP_{t} (GeV/c)",Rjet),nBinPt[0],ptmin[0],ptmax[0],nBinDpt[0],Dptmin[0],Dptmax[0]);
//Rebining the Particle level axes
TH2D *h2Dpt_Rebined_True = new TH2D(Form("h2Dpt_Rebined_True_R%03d",int(Rjet*100)),Form("#DeltaP_{t} distribution R =%.2f (Generator level -Rebinned);P_{t} (GeV/c);#DeltaP_{t} (GeV/c)",Rjet),nBinPt[1],ptmin[1],ptmax[1],nBinDpt[1],Dptmin[1],Dptmax[1]);
TH2D *h2Meas_m = new TH2D(Form("h2Meas_m_R%03d",int(Rjet*100)),Form("h2Meas_m_R =%.2f",Rjet),nBinPt[0],ptmin[0],ptmax[0],nBinDpt[0],Dptmin[0],Dptmax[0]);
TH2D *h2True_m = new TH2D(Form("h2True_m_R%03d",int(Rjet*100)),Form("h2True (measured) R =%.2f",Rjet),nBinPt[1],ptmin[1],ptmax[1],nBinDpt[1],Dptmin[1],Dptmax[1]);
//Feed-out Response
TH2D *h2Dpt_miss = new TH2D(Form("h2Dpt_miss_R%03d",int(Rjet*100)),Form("h2Dpt_miss_R =%.2f",Rjet),nBinPt[1],ptmin[1],ptmax[1],nBinDpt[1],Dptmin[1],Dptmax[1]);
TH2D *h2Dpt_miss_Det = new TH2D(Form("h2Dpt_miss_Det_R%03d",int(Rjet*100)),Form("h2Dpt_miss_Det_R =%.2f",Rjet),nBinPt[0],ptmin[0],ptmax[0],nBinDpt[0],Dptmin[0],Dptmax[0]);

//Now that the dimensions and binning configurations of the new 2D distributions at the two levels are defined, it's time to fill the histograms
//First the detector level histogram
for(Int_t ix=1;ix<=h2Dpt_Rebined_Meas->GetNbinsX();ix++){

      Double_t xlow = h2Dpt_Rebined_Meas->GetXaxis()->GetBinLowEdge(ix);
      Double_t xup = h2Dpt_Rebined_Meas->GetXaxis()->GetBinUpEdge(ix);
      Int_t jxlow = Dpt_pt_meas->GetXaxis()->FindBin(xlow+0.000001);
      Int_t jxup = Dpt_pt_meas->GetXaxis()->FindBin(xup-0.000001);
      for(Int_t iy = 1; iy<=h2Dpt_Rebined_Meas->GetNbinsY(); iy++) {
        Double_t ylow = h2Dpt_Rebined_Meas->GetYaxis()->GetBinLowEdge(iy);
        Double_t yup =h2Dpt_Rebined_Meas->GetYaxis()->GetBinUpEdge(iy);
        Int_t jylow = Dpt_pt_meas->GetYaxis()->FindBin(ylow+0.000001);
        Int_t jyup = Dpt_pt_meas->GetYaxis()->FindBin(yup-0.000001);
  
        Double_t err = 0.;
        Double_t con = Dpt_pt_meas->IntegralAndError(jxlow,jxup,jylow,jyup,err);
        h2Dpt_Rebined_Meas->SetBinContent(ix,iy,con);
        h2Dpt_Rebined_Meas->SetBinError(ix,iy,err);
         } //End of ybin loop
         } //End of xbin loop
                           

// Repeating the process for the particle level
  for(Int_t ix = 1; ix<=h2Dpt_Rebined_True->GetNbinsX(); ix++) {
      Double_t xlow = h2Dpt_Rebined_True->GetXaxis()->GetBinLowEdge(ix);
      Double_t xup = h2Dpt_Rebined_True->GetXaxis()->GetBinUpEdge(ix);
      Int_t jxlow = Dpt_pt_true->GetXaxis()->FindBin(xlow+0.000001);
      Int_t jxup = Dpt_pt_true->GetXaxis()->FindBin(xup-0.000001);
      for(Int_t iy = 1; iy<=h2Dpt_Rebined_True->GetNbinsY(); iy++) {
        Double_t ylow = h2Dpt_Rebined_True->GetYaxis()->GetBinLowEdge(iy);
        Double_t yup = h2Dpt_Rebined_True->GetYaxis()->GetBinUpEdge(iy);
        Int_t jylow = Dpt_pt_true->GetYaxis()->FindBin(ylow+0.000001);
        Int_t jyup = Dpt_pt_true->GetYaxis()->FindBin(yup-0.000001);
  
        Double_t err = 0.;
        Double_t con = Dpt_pt_true->IntegralAndError(jxlow,jxup,jylow,jyup,err);
        h2Dpt_Rebined_True->SetBinContent(ix,iy,con);
        h2Dpt_Rebined_True->SetBinError(ix,iy,err);
      }
    }


//And now create the RooUnfoldResponse object and fill it
    //response object for RooUnfold
    RooUnfoldResponse *fResponse = new RooUnfoldResponse(Form("Response_R%03d",int(Rjet*100)),"RM");
    fResponse->Setup(h2Dpt_Rebined_Meas,h2Dpt_Rebined_True);
  
    //Fill RooUnfoldResponse object
    Int_t* coord = new Int_t[nDim]; //Carries the bin coordinates
    Int_t nbin = RM_thn->GetNbins();
    for(Int_t bin=0; bin<nbin; bin++) {
      Double_t w = RM_thn->GetBinContent(bin,coord);
      Double_t Dpt_meas = RM_thn->GetAxis(3)->GetBinCenter(coord[3]);
      Double_t Dpt_true = RM_thn->GetAxis(2)->GetBinCenter(coord[2]);
      Double_t pt_meas = RM_thn->GetAxis(1)->GetBinCenter(coord[1]);
      Double_t pt_true = RM_thn->GetAxis(0)->GetBinCenter(coord[0]);
      if(Dpt_meas>=Dptmin[0] && Dpt_meas<Dptmax[0]
         && Dpt_true>=Dptmin[1] && Dpt_true<Dptmax[1]
         && pt_meas>=ptmin[0] && pt_meas<ptmax[0]
         && pt_true>=ptmin[1] && pt_true<ptmax[1]
        ){
        fResponse->Fill(pt_meas,Dpt_meas,pt_true,Dpt_true,w);
        h2Meas_m->Fill(pt_meas,Dpt_meas,w);
        h2True_m->Fill(pt_true,Dpt_true,w);}
      else {
        fResponse->Miss(pt_true,Dpt_true,w);
//        std::cout<<"Miss values:"<<pt_true<<" "<<pt_meas<<" " << Dpt_true<<" "<<Dpt_meas<<"\n";
        h2Dpt_miss->Fill(pt_true,Dpt_true,w);
        h2Dpt_miss_Det->Fill(pt_meas,Dpt_meas,w);
      }
    }
          TH2D* h2True = (TH2D*) h2True_m->Clone();
        h2True->Add(h2Dpt_miss,1);
        h2True->SetName(Form("h2True_R%03d",int(Rjet*100)));
        TH2D* h2Meas = (TH2D*) h2Meas_m->Clone();
        h2Meas->Add(h2Dpt_miss_Det,1);
        h2Meas->SetName(Form("h2Meas_R%03d",int(Rjet*100)));

    delete [] coord;

    TFile *fout = new TFile(Form("ResponseMatrix_PbPb_ConstSub_LeadCut_Updated_ptmin%d_Dptmin%d_%d.root",int(PT_MIN),int(DPT_MIN),int(DPT_MAX)),"UPDATE");
    RM_thn->Write();
    fResponse->Write(Form("Response_R%03d",int(Rjet*100)));
    Dpt_pt_meas->Write(Dpt_pt_meas->GetName(),TObject::kOverwrite);
    Dpt_pt_true->Write(Dpt_pt_true->GetName(),TObject::kOverwrite);
    h2True_m->Write(h2True_m->GetName(),TObject::kOverwrite);
    h2Meas_m->Write(h2Meas_m->GetName(),TObject::kOverwrite);
    h2Meas->Write(h2Meas->GetName(),TObject::kOverwrite);
    h2True->Write(h2True->GetName(),TObject::kOverwrite);
    h2Dpt_Rebined_Meas->Write(h2Dpt_Rebined_Meas->GetName(),TObject::kOverwrite);
    h2Dpt_Rebined_True->Write(h2Dpt_Rebined_True->GetName(),TObject::kOverwrite);
    h2Dpt_miss->Write(h2Dpt_miss->GetName(),TObject::kOverwrite);
    h2Dpt_miss_Det->Write(h2Dpt_miss_Det->GetName(),TObject::kOverwrite);
    fout->Write();
    fout->Close();


}

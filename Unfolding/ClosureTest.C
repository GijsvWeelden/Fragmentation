void SetStyle(TH1* h1);
TH2D* RM_normalization(TH2D* input_RM);
void normaliseHistRowByRow(TH2F* hist);

void ClosureTest(string testFileName = "AnalysisResults.root", string responseFileName = "responseMatrix_ppMC_ptmin0.root", Double_t RJET = 0.1, Int_t PT_LOW = 0 , Int_t PT_HIGH = 100, Int_t N_ITER = 1, Int_t BINWIDTH = 1)
{
  //Response Matrix visualisation
  TFile* responseFile = TFile::Open(responseFileName.Data());
  TH2D* hDetector = static_cast<TH2D*>(responseFile->Get("hDetector"));
  TH2D* hTruth    = static_cast<TH2D*>(responseFile->Get("hTruth"));

  Int_t nBinsPt      = hDetector->GetNbinsX();
  Double_t ptmin     = hDetector->GetXaxis()->GetXmin();
  Double_t ptmax     = hDetector->GetXaxis()->GetXmax();

  Int_t nBinsZ       = hDetector->GetNbinsY();
  Double_t zmin      = hDetector->GetYaxis()->GetXmin();
  Double_t zmax      = hDetector->GetYaxis()->GetXmax();

  Int_t nBinsPt_gen  = hTruth->GetNbinsX();
  Double_t ptmin_gen = hTruth->GetXaxis()->GetXmin();
  Double_t ptmax_gen = hTruth->GetXaxis()->GetXmax();

  Int_t nBinsZ_gen   = hTruth->GetNbinsY();
  Double_t zmin_gen  = hTruth->GetYaxis()->GetXmin();
  Double_t zmax_gen  = hTruth->GetYaxis()->GetXmax();

  Int_t ptWidth = int(hDetector->GetXaxis()->GetBinWidth(1));
  Int_t zWidth  = int(hDetector->GetYaxis()->GetBinWidth(1));

  Double_t zLow = zmin;
  Double_t zHigh = zmax;

  Double_t zLow_gen = zmin_gen;
  Double_t zHigh_gen = zmax_gen;

  //Plotting the RM
  //Set drawing style
  gStyle->SetOptStat(0);
  gStyle->SetTextFont(42);
  gStyle->SetTitleFont(42);
  gStyle->SetTextSize(18);
  gStyle->SetOptTitle(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLabelOffset(.001);
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetTitleSize(0.05);
  gStyle->SetTextSize(16);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLegendTextSize(0.025);
  gStyle->SetTitleXOffset(.7);
  gStyle->SetTitleYOffset(.7);
  gStyle->SetTitleOffset(.7);

  const Int_t nDim  = 4;
  Int_t ptDetAxis   = 0;
  Int_t zDetAxis    = 1;
  Int_t ptTruthAxis = 2;
  Int_t zTruthAxis  = 3;

  // Visualisation of response matrix
  THnSparse* responseMatrix = static_cast<THnSparse*>(responseFile->Get("responseMatrix"));
  responseMatrix ->UseCurrentStyle();
  Int_t ptLowBin_gen  = responseMatrix->GetAxis(ptTruthAxis)->FindBin(PT_LOW);
  Int_t ptHighBin_gen = responseMatrix->GetAxis(ptTruthAxis)->FindBin(PT_HIGH + 10) - 1;
  Int_t ptLowBin      = responseMatrix->GetAxis(ptDetAxis)->FindBin(PT_LOW);
  Int_t ptHighBin     = responseMatrix->GetAxis(ptDetAxis)->FindBin(PT_HIGH) - 1;
  responseMatrix->GetAxis(ptTruthAxis)->SetRange(ptLowBin_gen, ptHighBin_gen);
  responseMatrix->GetAxis(ptDetAxis)->SetRange(Ptlow_bin, Pthigh_bin);
  TH2D *ptRM = dynamic_cast<TH2D*>(responseMatrix->Projection(ptTruthAxis, ptDetAxis, "E"));
  ptRM->SetName("ptRM");
  ptRM->SetTitle("Response matrix projected onto jet pt");
  // RM_normalization(ptRM);
  normaliseHistRowByRow(ptRM);
  ptRM->GetZaxis()->SetMaximum(1);
  ptRM->GetZaxis()->SetMinimum(1e-5);

  Int_t zLow_bin_gen  = responseMatrix->GetAxis(zTruthAxis)->FindBin(zLow_gen);
  Int_t zHigh_bin_gen = responseMatrix->GetAxis(zTruthAxis)->FindBin(zHigh_gen)-1;
  Int_t zLow_bin      = responseMatrix->GetAxis(zDetAxis)->FindBin(zLow);
  Int_t zHigh_bin     = responseMatrix->GetAxis(zDetAxis)->FindBin(zHigh)-1;
  responseMatrix->GetAxis(zTruthAxis)->SetRange(zLow_bin_gen, zHigh_bin_gen);
  responseMatrix->GetAxis(zDetAxis)->SetRange(zLow_bin,zHigh_bin);
  TH2D *zRM = dynamic_cast<TH2D*>(responseMatrix->Projection(zTruthAxis, zDetAxis, "E"));
  zRM->SetName("zRM");
  zRM->SetTitle("Response matrix projected onto z");
  // RM_normalization(zRM);
  normaliseHistRowByRow(zRM);
  zRM->GetZaxis()->SetMaximum(1);
  zRM->GetZaxis()->SetMinimum(1e-5);

  // Introducing the test or pseudodata distributions
  // If the closure is trivial then the input here is the same as it was used in the UnfoldPreparation macro
  // Creating the histograms needed for the closure test
  string detectorTitle = "Detector level #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}";
  string truthTitle = "Particle level #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}";
  // Detector level - In Detector and in Truth acceptance (Response->Fill)
  TH2D *hDetectorMeasured = new TH2D("hDetectorMeasured", TString::Format("Measured %s", detectorTitle.c_str()), nBinsPt, ptmin, ptmax, nBinsZ, zmin, zmax);
  // Detector level - In Truth, but out of Detector acceptance (Response->Miss) (should be empty)
  TH2D *hDetectorMissed = new TH2D("hDetectorMissed", TString::Format("Missed %s", detectorTitle.c_str()), nBinsPt, ptmin, ptmax, nBinsZ, zmin, zmax);
  // Detector level - In Detector, but out of Truth acceptance (Response->Fake)
  TH2D *hDetectorFake = new TH2D("hDetectorFake", TString::Format("Fake %s", detectorTitle.c_str()),nBinsPt,ptmin,ptmax,nBinsZ,zmin,zmax);
  // Truth level - In Detector and in Truth acceptance (Response->Fill)
  TH2D *hTruthMeasured = new TH2D("hTruthMeasured", TString::Format("Measured %s", truthTitle.c_str()), nBinsPt_gen, ptmin_gen, ptmax_gen, nBinsZ_gen, zmin_gen, zmax_gen);
  //Truth level distribution within true but outside detector range (Response->Miss)
  TH2D *hTruthMissed = new TH2D("hTruthMissed", TString::Format("Missed %s", truthTitle.c_str()), nBinsPt_gen, ptmin_gen, ptmax_gen, nBinsZ_gen, zmin_gen, zmax_gen);
  // Truth level - In Detector, but out of Truth acceptance (Response->Fake) (should be empty)
  TH2D *hTruthFake = new TH2D("hTruthFake", TString::Format("Fake %s", truthTitle.c_str()), nBinsPt_gen, ptmin_gen, ptmax_gen, nBinsZ_gen, zmin_gen, zmax_gen);

  //Open the test/pseudodata file to retrieve the responseMatrix, its projections and to fill the above histograms
  TFile* testFile = TFile::Open(testFileName.Data());
  string testDirName = "jet-fragmentation/matching/jets";
  TDirectory* testDir = (TDirectory*)testFile->Get(TString::Format("%s", testDirName).Data());
  string testResponseMatrixName = "matchDetJetPtTrackProjPartJetPtTrackProj";
  THnSparseF* testResponseMatrix = (THnSparseF*)matchJetsDir->Get(TString::Format("%s", testResponseMatrixName.c_str()).Data());

  //Response matrix projections from the test file before rebinning
  TH2D *testResponseMatrixTruthProjection = dynamic_cast<TH2D*>(testResponseMatrix->Projection(zTruthAxis, ptTruthAxis, "E"));
  testResponseMatrixTruthProjection->SetName("testResponseMatrixTruthProjection");
  TH2D *testResponseMatrixDetectorProjection = dynamic_cast<TH2D*>(testResponseMatrix->Projection(zDetAxis, ptDetAxis, "E"));
  testResponseMatrixDetectorProjection->SetName("testResponseMatrixDetectorProjection");

  TH2D *testResponseMatrixDetectorProjection_rebinned = new TH2D("testResponseMatrixDetectorProjection_rebinned", "Detector level projection of test RM", nBinsPt, ptmin, ptmax, nBinsZ, zmin, zmax);
  TH2D *testResponseMatrixTruthProjection_rebinned = new TH2D("testResponseMatrixTruthProjection_rebinned", "Truth level projection of test RM", nBinsPt_gen, ptmin_gen, ptmax_gen, nBinsZ_gen, zmin_gen, zmax_gen);

  //Rebinning the detector level distribution from the test file
  for (Int_t ix = 1; ix <= testResponseMatrixDetectorProjection_rebinned->GetNbinsX(); ix++) {
    Double_t xMin = testResponseMatrixDetectorProjection_rebinned->GetXaxis()->GetBinLowEdge(ix);
    Double_t xMax = testResponseMatrixDetectorProjection_rebinned->GetXaxis()->GetBinUpEdge(ix);
    Int_t xBinLow = testResponseMatrixDetectorProjection->GetXaxis()->FindBin(xMin + 0.000001);
    Int_t xBinUp  = testResponseMatrixDetectorProjection->GetXaxis()->FindBin(xup - 0.000001);
    for (Int_t iy = 1; iy <= testResponseMatrixDetectorProjection_rebinned->GetNbinsY(); iy++) {
      Double_t yMin = testResponseMatrixDetectorProjection_rebinned->GetYaxis()->GetBinLowEdge(iy);
      Double_t yMax = testResponseMatrixDetectorProjection_rebinned->GetYaxis()->GetBinUpEdge(iy);
      Int_t yBinLow = testResponseMatrixDetectorProjection->GetYaxis()->FindBin(yMin + 0.000001);
      Int_t yBinUp  = testResponseMatrixDetectorProjection->GetYaxis()->FindBin(yMax - 0.000001);

      Double_t binError = 0.;
      Double_t binContent = testResponseMatrixDetectorProjection->IntegralAndError(xBinLow, xBinUp, yBinLow, yBinUp, binError);
      testResponseMatrixDetectorProjection_rebinned->SetBinContent(ix, iy, binContent);
      testResponseMatrixDetectorProjection_rebinned->SetBinError(ix, iy, binError);
    } // for iy
  }   // for ix

  //Rebinning the truth level distribution from the test file
  for(Int_t ix = 1; ix <= testResponseMatrixTruthProjection_rebinned->GetNbinsX(); ix++) {
    Double_t xMin = testResponseMatrixTruthProjection_rebinned->GetXaxis()->GetBinLowEdge(ix);
    Double_t xMax = testResponseMatrixTruthProjection_rebinned->GetXaxis()->GetBinUpEdge(ix);
    Int_t xBinLow = testResponseMatrixTruthProjection->GetXaxis()->FindBin(xMin + 0.000001);
    Int_t xBinUp  = testResponseMatrixTruthProjection->GetXaxis()->FindBin(xMax - 0.000001);
    for(Int_t iy = 1; iy <= testResponseMatrixTruthProjection_rebinned->GetNbinsY(); iy++) {
      Double_t yMin = testResponseMatrixTruthProjection_rebinned->GetYaxis()->GetBinLowEdge(iy);
      Double_t yMax = testResponseMatrixTruthProjection_rebinned->GetYaxis()->GetBinUpEdge(iy);
      Int_t yBinLow = testResponseMatrixTruthProjection->GetYaxis()->FindBin(yMin + 0.000001);
      Int_t yBinUp  = testResponseMatrixTruthProjection->GetYaxis()->FindBin(yMax - 0.000001);

      Double_t binError = 0.;
      Double_t binContent = testResponseMatrixTruthProjection->IntegralAndError(xBinLow, xBinUp, yBinLow, yBinUp, binError);
      testResponseMatrixTruthProjection_rebinned->SetBinContent(ix, iy, binContent);
      testResponseMatrixTruthProjection_rebinned->SetBinError(ix, iy, binError);
    } // for iy
  }   // for ix


  //Fill the measured, missed and Fake distribution at each level
  Int_t* coord = new Int_t[nDim]; //Carries the bin coordinates
  Int_t  nBins = testResponseMatrix->GetNbins();
  bool inAcceptanceDetector = false;
  bool inAcceptanceTruth    = false;

  for(Int_t bin = 0; bin < nBins; bin++) {
    inAcceptanceDetector = false;
    inAcceptanceTruth    = false;
    Double_t w          = testResponseMatrix->GetBinContent(bin, coord);
    Double_t zMeasured  = testResponseMatrix->GetAxis(zDetAxis)->GetBinCenter(coord[zDetAxis]);
    Double_t zTruth     = testResponseMatrix->GetAxis(zTruthAxis)->GetBinCenter(coord[zTruthAxis]);
    Double_t ptMeasured = testResponseMatrix->GetAxis(ptDetAxis)->GetBinCenter(coord[ptDetAxis]);
    Double_t ptTruth    = testResponseMatrix->GetAxis(ptTruthAxis)->GetBinCenter(coord[ptTruthAxis]);

    inAcceptanceDetector = (zMeasured  >= zmin[0]) * (zMeasured  < zmax[0]) * (ptMeasured >= ptmin[0]) * (ptMeasured < ptmax[0]);
    inAcceptanceTruth = (zTruth >= zmin[1]) * (zTruth < zmax[1]) * (ptTruth >= ptmin[1]) * (ptTruth < ptmax[1]);

    if(inAcceptanceTruth && inAcceptanceDetector) {
      hDetectorMeasured->Fill(ptMeasured, zMeasured, w);
      hTruthMeasured->Fill(ptTruth, zTruth, w);
    }
    else if(inAcceptanceTruth) {
      hTruthMissed->Fill(ptTruth, zTruth, w);
      hDetectorMissed->Fill(ptMeasured, zMeasured, w);
    }
    else if (inAcceptanceDetector) {
      hTruthFake->Fill(ptTruth, zTruth, w);
      hDetectorFake->Fill(ptMeasured, zMeasured, w);
    }
  } // for each bin
  delete [] coord;

  // Complete truth distribution
  TH2D* h2TrueFull = (TH2D*) hTruthMeasured->Clone("h2TrueFull");
  h2TrueFull->Add(hTruthMissed, 1);
  // Complete measured distribution
  TH2D* h2MeasuredFull = (TH2D*) hDetectorMeasured->Clone("h2MeasuredFull");
  h2MeasuredFull->Add(hDetectorFake, 1);

  //This is the section where the actual unfolding starts and RooUnfold is involved
  string ruResponseName = "Response";
  string unfoldName = TString::Format("Bayesian_Unfolding_ptmin%d", int(ptmin)).Data();
  string unfoldTitle = unfoldName;
  int doSmoothing = 0;
  //Load the response object from the Response file
  RooUnfoldResponse* ruResponse = static_cast<RooUnfoldResponse*>(responseFile->Get(TString::Format("%s", ruResponseName.c_str()).Data()));
  //Create the Bayesian unfolding object
  RooUnfoldBayes Unfolding_bayes(ruResponse, testResponseMatrixDetectorProjection_rebinned, N_ITER, doSmoothing, Unfoldname, unfoldTitle);
  //Extracting the unfolded distribution as a histogram
  TH2D* unfoldedTestZ = (TH2D*)Unfolding_bayes.Hreco(RooUnfold::kCovToy);
  unfoldedTestZ->SetName("unfoldedTestZ_R%03d");
  auto Unf_Bayes = &Unfolding_bayes;
  //Applying the Response matrix to the unfolded result
  TH2D* refoldedTestZ = (TH2D*)ruResponse->ApplyToTruth(unfoldedTestZ, "refoldedTestZ");
  //Include the fakes in order to compare with the input distribution
  refoldedTestZ->Add(hDetectorFake);

  // The unfolding is complete, the following sections refer to projections in one dimension and calculating the ratio between distributions
  Int_t ptBinLow      = hDetector->GetXaxis()->FindBin(PT_LOW);
  Int_t ptBinHigh     = hDetector->GetXaxis()->FindBin(PT_HIGH) - 1;
  Int_t ptBinLow_gen  = hTruth->GetXaxis()->FindBin(PT_LOW);
  Int_t ptBinHigh_gen = hTruth->GetXaxis()->FindBin(PT_HIGH) - 1;
  string projectionName = TString::Format("d#it{N}/d#it{z}, #it{p}_{T}^{jet} #in [%d, %d] GeV/#it{c}", PT_LOW, PT_HIGH).Data();
  string ratioName      = TString::Format("d#it{N}/d#it{z} ratio, #it{p}_{T}^{jet} #in [%d, %d] GeV/#it{c}", PT_LOW, PT_HIGH).Data();

  //--------------Training sample (Sample creating the response)--------------------//
    //Truth level
  TH1D* trainingTruth = hTruth->ProjectionY("trainingTruth", ptBinLow_gen, ptBinHigh_gen);
  trainingTruth->SetTitle(TString::Format("%s (truth level, training sample); #it{z}", projectionName.c_str()).Data());
  trainingTruth->Scale(1./trainingTruth->Integral(), "width");
    //Detector level
  TH1D*trainingDetector = hDetector->ProjectionY("trainingDetector", ptBinLow, ptBinHigh);
  trainingDetector->SetTitle(TString::Format("%s (MC Detector level, training sample); #it{z}", projectionName.c_str()).Data());
  trainingDetector->Scale(1./trainingDetector->Integral(), "width");

  //--------------Test sample (MC sample treated as pseudodata)--------------------//
    //Truth level
  TH1D* testTruth = testResponseMatrixTruthProjection_rebinned->ProjectionY("testTruth", ptBinLow_gen, ptBinHigh_gen);
  testTruth->SetTitle(TString::Format("%s (truth level, test sample); #it{z}", projectionName.c_str()).Data());
  testTruth->Scale(1./testTruth->Integral(), "width");
    //Detector level
  TH1D* testDetector = testResponseMatrixDetectorProjection_rebinned->ProjectionY("testDetector", ptBinLow, ptBinHigh);
  testDetector->SetTitle(TString::Format("%s (detector level, test sample); #it{z}", projectionName.c_str()).Data());
  testDetector->Scale(1./testDetector->Integral(), "width");
    //Unfolded
  TH1D* testUnfolded = unfoldedTestZ->ProjectionY("testUnfolded");
  testUnfolded->SetTitle(TString::Format("%s (unfolded, test sample); #it{z}", projectionName.c_str()));
  testUnfolded->Scale(1./testUnfolded->Integral(), "width");
    //Refolded
  TH1D* testRefolded = refoldedTestZ->ProjectionY("testRefolded");
  testRefolded->SetTitle(TString::Format("%s (refolded, test sample); #it{z}", projectionName.c_str()));
  testRefolded->Scale(1./testRefolded->Integral(), "width");

  // Ratios of the projections
  // Training distributions - Detector response corrections
  TH1D* trainingTruthOverDetector = (TH1D*)trainingTruth->Clone("trainingTruthOverDetector");
  trainingTruthOverDetector->Reset();
  for (int iBin = 0; iBin < trainingDetector->GetNbinsX(); iBin++) {
    trainingTruthOverDetector->SetBinContent(iBin, trainingDetector->GetBinContent(iBin));
    trainingTruthOverDetector->SetBinError(iBin, trainingDetector->GetBinError(iBin));
  }
  trainingTruthOverDetector->Divide(trainingTruth, trainingTruthOverDetector);
  trainingTruthOverDetector->SetTitle(TString::Format("Truth over Detector %s (training sample); #it{z}; #frac{Truth}{Detector}", ratioName.c_str()));

  //Test distributions - Detector response corrections
  TH1D* testUnfoldedOverDetector = (TH1D*)testUnfolded->Clone("testUnfoldedOverDetector");
  testUnfoldedOverDetector->Reset();
  for (int iBin = 0; iBin < testDetector->GetNbinsX(); iBin++) {
    testUnfoldedOverDetector->SetBinContent(iBin, testDetector->GetBinContent(iBin));
    testUnfoldedOverDetector->SetBinError(iBin, testDetector->GetBinError(iBin));
  }
  testUnfoldedOverDetector->Divide(testUnfolded, testUnfoldedOverDetector);
  testUnfoldedOverDetector->SetTitle(TString::Format("Unfolded over Detector %s (test sample); #it{z}; #frac{Unfolded}{Detector}", ratioName.c_str()));
  // testUnfoldedOverDetector->SetName(TString::Format("Ratio_Unftest_R%03d",int(Rjet*100)));

  //Closure on the truth level
  TH1D* testUnfoldedOverTruth = (TH1D*)testUnfolded->Clone("testUnfoldedOverTruth");
  testUnfoldedOverTruth->Reset();
  testUnfoldedOverTruth->Divide(testUnfolded, testTruth);
  testUnfoldedOverTruth->SetTitle(TString::Format("Unfolded over Truth %s (test sample); #it{z}; #frac{Unfolded}{Truth}", ratioName.c_str()));
  // testUnfoldedOverTruth->SetName(TString::Format("Ratio_Unftest_truth_R%03d",int(Rjet*100)));

  //Closure on the detector level
  TH1D* testRefoldedOverDetector = (TH1D*)testRefolded->Clone("testRefoldedOverDetector");
  testRefoldedOverDetector->Reset();
  testRefoldedOverDetector->Divide(testRefolded, testDetector);
  testRefoldedOverDetector->SetTitle(TString::Format("Refolded over Detector %s (test sample); #it{z}; #frac{Refolded}{Detector}",ratioName.c_str()));
  // testRefoldedOverDetector->SetName(TString::Format("Ratio_Reftest_R%03d",int(Rjet*100)));

  //Saving the output in a new file
  string outFileName = TString::Format("closureTest_binWidth%d_pt%d-%d_ptmin%d_z%d-%d_nIter%d.root",
                                       BINWIDTH, PT_LOW, PT_HIGH, ptmin, zmin, zmax, N_ITER).Data();
  // TString Outputname = Form("ClosureTest_ConstSub_wFakes_new_Binw%d_LC8_pt_%d_%d_ptmin%d_zmin%d_max%d_Iter%d.root",BinWidth,PT_LOW,PT_HIGH,int(ptmin),int(zmin),int(zmax),Niter);
  TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");

  // Response objects
  ptRM->Write(ptRM->GetName(), TObject::kOverwrite);
  zRM->Write(zRM->GetName(), TObject::kOverwrite);
  responseMatrix->Write(responseMatrix->GetName(), TObject::kOverwrite);
  hTruthMeasured->Write(hTruthMeasured->GetName(), TObject::kOverwrite);
  hDetectorMeasured->Write(hDetectorMeasured->GetName(), TObject::kOverwrite);
  hDetectorMeasured->Write(hDetectorMeasured->GetName(), TObject::kOverwrite);
  hTruthMeasured->Write(hTruthMeasured->GetName(), TObject::kOverwrite);
  hTruthMissed->Write(hTruthMissed->GetName(), TObject::kOverwrite);
  hDetectorMissed->Write(hDetectorMissed->GetName(), TObject::kOverwrite);
  hDetectorFake->Write(hDetectorFake->GetName(), TObject::kOverwrite);
  hTruthFake->Write(hTruthFake->GetName(), TObject::kOverwrite);
  // Unfolding objects
  testResponseMatrixTruthProjection_rebinned->Write(testResponseMatrixTruthProjection_rebinned->GetName(), TObject::kOverwrite);
  testResponseMatrixDetectorProjection_rebinned->Write(testResponseMatrixDetectorProjection_rebinned->GetName(), TObject::kOverwrite);
  Unf_Bayes->Write(Unf_Bayes->GetName(), TObject::kOverwrite);
  unfoldedTestZ->Write(unfoldedTestZ->GetName(), TObject::kOverwrite);
  refoldedTestZ->Write(refoldedTestZ->GetName(), TObject::kOverwrite);
  trainingDetector->Write(trainingDetector->GetName(), TObject::kOverwrite);
  trainingTruth->Write(trainingTruth->GetName(), TObject::kOverwrite);
  testDetector->Write(testDetector->GetName(), TObject::kOverwrite);
  testTruth->Write(testTruth->GetName(), TObject::kOverwrite);
  testUnfolded->Write(testUnfolded->GetName(), TObject::kOverwrite);
  testRefolded->Write(testRefolded->GetName(), TObject::kOverwrite);

  trainingTruthOverDetector->Write(trainingTruthOverDetector->GetName(), TObject::kOverwrite);
  testUnfoldedOverDetector->Write(testUnfoldedOverDetector->GetName(), TObject::kOverwrite);
  testUnfoldedOverTruth->Write(testUnfoldedOverTruth->GetName(), TObject::kOverwrite);
  testRefoldedOverDetector->Write(testRefoldedOverDetector->GetName(), TObject::kOverwrite);
  testFile->Close();
  responseFile->Close();
  outFile->Close();
}

void SetStyle(TH1* h1){
  gStyle->SetOptStat(0); gStyle->SetTextFont(42); gStyle->SetTitleFont(42);
  gStyle->SetOptTitle(1); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
  gStyle->SetTextSize(18);gStyle->SetLabelOffset(.001); gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"X"); gStyle->SetTitleSize(0.05);
  gStyle->SetLegendTextSize(0.025); gStyle->SetTitleXOffset(.7);
  gStyle->SetTitleYOffset(.7); gStyle->SetTitleOffset(.7);
  h1->SetMarkerStyle(33); h1->SetMarkerSize(2); h1->SetLineWidth(2);
}

TH2D* RM_normalization(TH2D* input_RM){
  Int_t nBinsDpt[2]= {input_RM->GetXaxis()->GetNbins(),input_RM->GetYaxis()->GetNbins()};
  TH2D* RM_norm =(TH2D*)input_RM->Clone();
  RM_norm->Reset();
  for(int iy=1;iy<=nBinsDpt[1];iy++){
    Double_t sum = input_RM->Integral(1,nBinsDpt[0],iy,iy);
    for(int ix=1;ix<=nBinsDpt[0];ix++){
        Double_t binContent = input_RM->GetBinContent(ix,iy);
        RM_norm->SetBinContent(ix,iy,binContent/sum);
      }
    }
  return RM_norm;
}

// Normalise 2D histogram row-by-row
void normaliseHistRowByRow(TH2F* hist)
{
  int firstColBin = 1, lastColBin = hist->GetNbinsX();
  int firstRowBin = 1, lastRowBin = hist->GetNbinsY();
  for (int iRow = 1; iRow <= lastRowBin; iRow++) {
    double integral = hist->Integral(firstColBin, lastColBin, iRow, iRow);
    if (integral < 1) { continue; }
    for (int iCol = 1; iCol <= lastColBin; iCol++) {
      double binContent = hist->GetBinContent(iCol, iRow);
      binContent /= integral;
      hist->SetBinContent(iCol, iRow, binContent);
    }
  }
}

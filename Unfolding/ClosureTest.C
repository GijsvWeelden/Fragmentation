void SetStyle(TH1* h1);
TH2F* RM_normalization(TH2F* input_RM);
void normaliseHistRowByRow(TH2F* hist);

void ClosureTest(string testFileName = "AnalysisResults.root", string responseFileName = "RooUnfoldResponse.root", Int_t PT_LOW = 20 , Int_t PT_HIGH = 80, Int_t N_ITER = 1)
{
  //Response Matrix visualisation
  TFile* responseFile = TFile::Open(TString::Format("%s", responseFileName.c_str()).Data());
  TH2F* hDetector = (TH2F*)responseFile->Get("hDetector");
  TH2F* hTruth    = (TH2F*)responseFile->Get("hTruth");

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
  responseMatrix->SetName("responseMatrix");
  responseMatrix ->UseCurrentStyle();
  Int_t ptLowBin_gen  = responseMatrix->GetAxis(ptTruthAxis)->FindBin(PT_LOW);
  Int_t ptHighBin_gen = responseMatrix->GetAxis(ptTruthAxis)->FindBin(PT_HIGH) - 1;
  Int_t ptLowBin      = responseMatrix->GetAxis(ptDetAxis)->FindBin(PT_LOW);
  Int_t ptHighBin     = responseMatrix->GetAxis(ptDetAxis)->FindBin(PT_HIGH) - 1;
  responseMatrix->GetAxis(ptTruthAxis)->SetRange(ptLowBin_gen, ptHighBin_gen);
  responseMatrix->GetAxis(ptDetAxis)->SetRange(ptLowBin, ptHighBin);
  TH2F *ptRM = (TH2F*)responseMatrix->Projection(ptTruthAxis, ptDetAxis, "E");
  ptRM->SetName("ptRM");
  ptRM->SetTitle("Response matrix projected onto jet pt");
  normaliseHistRowByRow(ptRM);
  ptRM->GetZaxis()->SetRangeUser(1e-5, 1);

  Int_t zLow_bin_gen  = responseMatrix->GetAxis(zTruthAxis)->FindBin(zLow_gen);
  Int_t zHigh_bin_gen = responseMatrix->GetAxis(zTruthAxis)->FindBin(zHigh_gen)-1;
  Int_t zLow_bin      = responseMatrix->GetAxis(zDetAxis)->FindBin(zLow);
  Int_t zHigh_bin     = responseMatrix->GetAxis(zDetAxis)->FindBin(zHigh)-1;
  responseMatrix->GetAxis(zTruthAxis)->SetRange(zLow_bin_gen, zHigh_bin_gen);
  responseMatrix->GetAxis(zDetAxis)->SetRange(zLow_bin,zHigh_bin);
  TH2F *zRM = (TH2F*)responseMatrix->Projection(zTruthAxis, zDetAxis, "E");
  zRM->SetName("zRM");
  zRM->SetTitle("Response matrix projected onto z");
  normaliseHistRowByRow(zRM);
  zRM->GetZaxis()->SetRangeUser(1e-5, 1);

  // Introducing the test or pseudodata distributions
  // If the closure is trivial then the input here is the same as it was used in the UnfoldPreparation macro
  // Creating the histograms needed for the closure test
  string detectorTitle = "Detector level #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}";
  string truthTitle = "Particle level #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}";
  // Detector level - In Detector and in Truth acceptance (Response->Fill)
  TH2F *hDetectorMeasured = new TH2F("hDetectorMeasured", TString::Format("Measured %s", detectorTitle.c_str()), nBinsPt, ptmin, ptmax, nBinsZ, zmin, zmax);
  // Detector level - In Truth, but out of Detector acceptance (Response->Miss) (should be empty)
  TH2F *hDetectorMissed = new TH2F("hDetectorMissed", TString::Format("Missed %s", detectorTitle.c_str()), nBinsPt, ptmin, ptmax, nBinsZ, zmin, zmax);
  // Detector level - In Detector, but out of Truth acceptance (Response->Fake)
  TH2F *hDetectorFake = new TH2F("hDetectorFake", TString::Format("Fake %s", detectorTitle.c_str()),nBinsPt,ptmin,ptmax,nBinsZ,zmin,zmax);
  // Truth level - In Detector and in Truth acceptance (Response->Fill)
  TH2F *hTruthMeasured = new TH2F("hTruthMeasured", TString::Format("Measured %s", truthTitle.c_str()), nBinsPt_gen, ptmin_gen, ptmax_gen, nBinsZ_gen, zmin_gen, zmax_gen);
  //Truth level distribution within true but outside detector range (Response->Miss)
  TH2F *hTruthMissed = new TH2F("hTruthMissed", TString::Format("Missed %s", truthTitle.c_str()), nBinsPt_gen, ptmin_gen, ptmax_gen, nBinsZ_gen, zmin_gen, zmax_gen);
  // Truth level - In Detector, but out of Truth acceptance (Response->Fake) (should be empty)
  TH2F *hTruthFake = new TH2F("hTruthFake", TString::Format("Fake %s", truthTitle.c_str()), nBinsPt_gen, ptmin_gen, ptmax_gen, nBinsZ_gen, zmin_gen, zmax_gen);

  //Open the test/pseudodata file to retrieve the responseMatrix, its projections and to fill the above histograms
  TFile* testFile = TFile::Open(TString::Format("%s", testFileName.c_str()).Data());
  string testDirName = "jet-fragmentation/matching/jets";
  TDirectory* testDir = (TDirectory*)testFile->Get(TString::Format("%s", testDirName.c_str()).Data());
  string testResponseMatrixName = "matchDetJetPtTrackProjPartJetPtTrackProj";
  THnSparseF* testResponseMatrix = (THnSparseF*)testDir->Get(TString::Format("%s", testResponseMatrixName.c_str()).Data());

  //Response matrix projections from the test file before rebinning
  TH2F *testResponseMatrixTruthProjection = (TH2F*)testResponseMatrix->Projection(zTruthAxis, ptTruthAxis, "E");
  testResponseMatrixTruthProjection->SetName("testResponseMatrixTruthProjection");
  TH2F *testResponseMatrixDetectorProjection = (TH2F*)testResponseMatrix->Projection(zDetAxis, ptDetAxis, "E");
  testResponseMatrixDetectorProjection->SetName("testResponseMatrixDetectorProjection");

  TH2F *testResponseMatrixDetectorProjection_rebinned = new TH2F("testResponseMatrixDetectorProjection_rebinned", "Detector level projection of test RM", nBinsPt, ptmin, ptmax, nBinsZ, zmin, zmax);
  TH2F *testResponseMatrixTruthProjection_rebinned = new TH2F("testResponseMatrixTruthProjection_rebinned", "Truth level projection of test RM", nBinsPt_gen, ptmin_gen, ptmax_gen, nBinsZ_gen, zmin_gen, zmax_gen);

  //Rebinning the detector level distribution from the test file
  for (Int_t ix = 1; ix <= testResponseMatrixDetectorProjection_rebinned->GetNbinsX(); ix++) {
    Double_t xMin = testResponseMatrixDetectorProjection_rebinned->GetXaxis()->GetBinLowEdge(ix);
    Double_t xMax = testResponseMatrixDetectorProjection_rebinned->GetXaxis()->GetBinUpEdge(ix);
    Int_t xBinLow = testResponseMatrixDetectorProjection->GetXaxis()->FindBin(xMin + 0.000001);
    Int_t xBinUp  = testResponseMatrixDetectorProjection->GetXaxis()->FindBin(xMax - 0.000001);
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

    inAcceptanceDetector = (zMeasured  >= zmin) * (zMeasured  < zmax) * (ptMeasured >= ptmin) * (ptMeasured < ptmax);
    inAcceptanceTruth = (zTruth >= zmin_gen) * (zTruth < zmax_gen) * (ptTruth >= ptmin_gen) * (ptTruth < ptmax_gen);

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
  TH2F* h2TrueFull = (TH2F*) hTruthMeasured->Clone("h2TrueFull");
  h2TrueFull->Add(hTruthMissed, 1);
  // Complete measured distribution
  TH2F* h2MeasuredFull = (TH2F*) hDetectorMeasured->Clone("h2MeasuredFull");
  h2MeasuredFull->Add(hDetectorFake, 1);

  //This is the section where the actual unfolding starts and RooUnfold is involved
  string ruResponseName = "Response";
  string unfoldName = TString::Format("Bayesian_Unfolding_ptmin%d", int(ptmin)).Data();
  string unfoldTitle = unfoldName;
  int doSmoothing = 0;
  //Load the response object from the Response file
  RooUnfoldResponse* ruResponse = static_cast<RooUnfoldResponse*>(responseFile->Get(TString::Format("%s", ruResponseName.c_str()).Data()));
  //Create the Bayesian unfolding object
  RooUnfoldBayes Unfolding_bayes(ruResponse, testResponseMatrixDetectorProjection_rebinned, N_ITER, doSmoothing, unfoldName.c_str(), unfoldTitle.c_str());
  //Extracting the unfolded distribution as a histogram
  TH2F* unfoldedTest = (TH2F*)Unfolding_bayes.Hreco(RooUnfold::kCovToy);
  unfoldedTest->SetName("unfoldedTest");
  auto Unf_Bayes = &Unfolding_bayes;
  //Applying the Response matrix to the unfolded result
  TH2F* refoldedTest = (TH2F*)ruResponse->ApplyToTruth(unfoldedTest, "refoldedTest");
  //Include the fakes in order to compare with the input distribution
  refoldedTest->Add(hDetectorFake);

  // The unfolding is complete, the following sections refer to projections in one dimension and calculating the ratio between distributions
  Int_t ptBinLow      = hDetector->GetXaxis()->FindBin(PT_LOW);
  Int_t ptBinHigh     = hDetector->GetXaxis()->FindBin(PT_HIGH) - 1;
  Int_t ptBinLow_gen  = hTruth->GetXaxis()->FindBin(PT_LOW);
  Int_t ptBinHigh_gen = hTruth->GetXaxis()->FindBin(PT_HIGH) - 1;
  string projectionNameZ  = TString::Format("d#it{N}/d#it{z}, #it{p}_{T}^{jet} #in [%d, %d] GeV/#it{c}", PT_LOW, PT_HIGH).Data();
  string ratioNameZ       = TString::Format("d#it{N}/d#it{z} ratio, #it{p}_{T}^{jet} #in [%d, %d] GeV/#it{c}", PT_LOW, PT_HIGH).Data();
  string projectionNamePt = TString::Format("d#it{N}/d#it{p}_{T}, #it{p}_{T}^{jet} #in [%d, %d] GeV/#it{c}", PT_LOW, PT_HIGH).Data();
  string ratioNamePt      = TString::Format("d#it{N}/d#it{p}_{T} ratio, #it{p}_{T}^{jet} #in [%d, %d] GeV/#it{c}", PT_LOW, PT_HIGH).Data();

  //--------------Training sample (Sample creating the response)--------------------//
    //Truth level
  TH1D* trainingPtTruth = hTruth->ProjectionX("trainingPtTruth", ptBinLow_gen, ptBinHigh_gen);
  trainingPtTruth->SetTitle(TString::Format("%s (truth level, training sample); #it{p}_{T}", projectionNamePt.c_str()).Data());
  trainingPtTruth->Scale(1./trainingPtTruth->Integral(), "width");

  TH1D* trainingZTruth = hTruth->ProjectionY("trainingZTruth", ptBinLow_gen, ptBinHigh_gen);
  trainingZTruth->SetTitle(TString::Format("%s (truth level, training sample); #it{z}", projectionNameZ.c_str()).Data());
  trainingZTruth->Scale(1./trainingZTruth->Integral(), "width");
    //Detector level
  TH1D* trainingPtDetector = hDetector->ProjectionX("trainingPtDetector", ptBinLow, ptBinHigh);
  trainingPtDetector->SetTitle(TString::Format("%s (MC Detector level, training sample); #it{p}_{T}", projectionNamePt.c_str()).Data());
  trainingPtDetector->Scale(1./trainingPtDetector->Integral(), "width");

  TH1D* trainingZDetector = hDetector->ProjectionY("trainingZDetector", ptBinLow, ptBinHigh);
  trainingZDetector->SetTitle(TString::Format("%s (MC Detector level, training sample); #it{z}", projectionNameZ.c_str()).Data());
  trainingZDetector->Scale(1./trainingZDetector->Integral(), "width");

  //--------------Test sample (MC sample treated as pseudodata)--------------------//
    //Truth level
  TH1D* testPtTruth = testResponseMatrixTruthProjection_rebinned->ProjectionX("testPtTruth", ptBinLow_gen, ptBinHigh_gen);
  testPtTruth->SetTitle(TString::Format("%s (truth level, test sample); #it{p}_{T}", projectionNamePt.c_str()).Data());
  testPtTruth->Scale(1./testPtTruth->Integral(), "width");

  TH1D* testZTruth = testResponseMatrixTruthProjection_rebinned->ProjectionY("testZTruth", ptBinLow_gen, ptBinHigh_gen);
  testZTruth->SetTitle(TString::Format("%s (truth level, test sample); #it{z}", projectionNameZ.c_str()).Data());
  testZTruth->Scale(1./testZTruth->Integral(), "width");
    //Detector level
  TH1D* testPtDetector = testResponseMatrixDetectorProjection_rebinned->ProjectionX("testPtDetector", ptBinLow, ptBinHigh);
  testPtDetector->SetTitle(TString::Format("%s (detector level, test sample); #it{p}_{T}", projectionNamePt.c_str()).Data());
  testPtDetector->Scale(1./testPtDetector->Integral(), "width");

  TH1D* testZDetector = testResponseMatrixDetectorProjection_rebinned->ProjectionY("testZDetector", ptBinLow, ptBinHigh);
  testZDetector->SetTitle(TString::Format("%s (detector level, test sample); #it{z}", projectionNameZ.c_str()).Data());
  testZDetector->Scale(1./testZDetector->Integral(), "width");
    //Unfolded
  TH1D* testPtUnfolded = unfoldedTest->ProjectionX("testPtUnfolded");
  testPtUnfolded->SetTitle(TString::Format("%s (unfolded, test sample); #it{p}_{T}", projectionNamePt.c_str()));
  testPtUnfolded->Scale(1./testPtUnfolded->Integral(), "width");

  TH1D* testZUnfolded = unfoldedTest->ProjectionY("testZUnfolded");
  testZUnfolded->SetTitle(TString::Format("%s (unfolded, test sample); #it{z}", projectionNameZ.c_str()));
  testZUnfolded->Scale(1./testZUnfolded->Integral(), "width");
    //Refolded
  TH1D* testPtRefolded = refoldedTest->ProjectionX("testPtRefolded");
  testPtRefolded->SetTitle(TString::Format("%s (refolded, test sample); #it{p}_{T}", projectionNameZ.c_str()));
  testPtRefolded->Scale(1./testPtRefolded->Integral(), "width");

  TH1D* testZRefolded = refoldedTest->ProjectionY("testZRefolded");
  testZRefolded->SetTitle(TString::Format("%s (refolded, test sample); #it{z}", projectionNameZ.c_str()));
  testZRefolded->Scale(1./testZRefolded->Integral(), "width");

  // Ratios of the projections
  // Training distributions - Detector response corrections
  TH1D* trainingPtTruthOverDetector = (TH1D*)trainingPtTruth->Clone("trainingPtTruthOverDetector");
  trainingPtTruthOverDetector->Reset();
  for (int iBin = 0; iBin < trainingPtDetector->GetNbinsX(); iBin++) {
    trainingPtTruthOverDetector->SetBinContent(iBin, trainingPtDetector->GetBinContent(iBin));
    trainingPtTruthOverDetector->SetBinError(iBin, trainingPtDetector->GetBinError(iBin));
  }
  trainingPtTruthOverDetector->Divide(trainingPtTruth, trainingPtTruthOverDetector);
  trainingPtTruthOverDetector->SetTitle(TString::Format("Truth over Detector %s (training sample); #it{p}_{T}; #frac{Truth}{Detector}", ratioNamePt.c_str()));
  trainingPtTruthOverDetector->Print();

  TH1D* trainingZTruthOverDetector = (TH1D*)trainingZTruth->Clone("trainingZTruthOverDetector");
  trainingZTruthOverDetector->Reset();
  for (int iBin = 0; iBin < trainingZDetector->GetNbinsX(); iBin++) {
    trainingZTruthOverDetector->SetBinContent(iBin, trainingZDetector->GetBinContent(iBin));
    trainingZTruthOverDetector->SetBinError(iBin, trainingZDetector->GetBinError(iBin));
  }
  trainingZTruthOverDetector->Divide(trainingZTruth, trainingZTruthOverDetector);
  trainingZTruthOverDetector->SetTitle(TString::Format("Truth over Detector %s (training sample); #it{z}; #frac{Truth}{Detector}", ratioNameZ.c_str()));
  trainingZTruthOverDetector->Print();

  //Test distributions - Detector response corrections
  TH1D* testPtUnfoldedOverDetector = (TH1D*)testPtUnfolded->Clone("testPtUnfoldedOverDetector");
  testPtUnfoldedOverDetector->Reset();
  for (int iBin = 0; iBin < testPtDetector->GetNbinsX(); iBin++) {
    testPtUnfoldedOverDetector->SetBinContent(iBin, testPtDetector->GetBinContent(iBin));
    testPtUnfoldedOverDetector->SetBinError(iBin, testPtDetector->GetBinError(iBin));
  }
  testPtUnfoldedOverDetector->Divide(testPtUnfolded, testPtUnfoldedOverDetector);
  testPtUnfoldedOverDetector->SetTitle(TString::Format("Unfolded over Detector %s (test sample); #it{z}; #frac{Unfolded}{Detector}", ratioNamePt.c_str()));
  testPtUnfoldedOverDetector->Print();

  TH1D* testZUnfoldedOverDetector = (TH1D*)testZUnfolded->Clone("testZUnfoldedOverDetector");
  testZUnfoldedOverDetector->Reset();
  for (int iBin = 0; iBin < testZDetector->GetNbinsX(); iBin++) {
    testZUnfoldedOverDetector->SetBinContent(iBin, testZDetector->GetBinContent(iBin));
    testZUnfoldedOverDetector->SetBinError(iBin, testZDetector->GetBinError(iBin));
  }
  testZUnfoldedOverDetector->Divide(testZUnfolded, testZUnfoldedOverDetector);
  testZUnfoldedOverDetector->SetTitle(TString::Format("Unfolded over Detector %s (test sample); #it{z}; #frac{Unfolded}{Detector}", ratioNameZ.c_str()));
  testZUnfoldedOverDetector->Print();

  //Closure on the truth level
  TH1D* testPtUnfoldedOverTruth = (TH1D*)testPtUnfolded->Clone("testPtUnfoldedOverTruth");
  testPtUnfoldedOverTruth->Reset();
  // for (int iBin = 0; iBin < testPtTruth->GetNbinsX(); iBin++) {
  //   testZUnfoldedOverDetector->SetBinContent(iBin, testPtTruth->GetBinContent(iBin));
  //   testZUnfoldedOverDetector->SetBinError(iBin, testPtTruth->GetBinError(iBin));
  // }
  testPtUnfoldedOverTruth->Divide(testPtUnfolded, testPtTruth);
  testPtUnfoldedOverTruth->SetTitle(TString::Format("Unfolded over Truth %s (test sample); #it{p}_{T}; #frac{Unfolded}{Truth}", ratioNamePt.c_str()));
  testPtUnfoldedOverTruth->Print();

  TH1D* testZUnfoldedOverTruth = (TH1D*)testZUnfolded->Clone("testZUnfoldedOverTruth");
  testZUnfoldedOverTruth->Reset();
  testZUnfoldedOverTruth->Divide(testZUnfolded, testZTruth);
  testZUnfoldedOverTruth->SetTitle(TString::Format("Unfolded over Truth %s (test sample); #it{z}; #frac{Unfolded}{Truth}", ratioNameZ.c_str()));
  testZUnfoldedOverTruth->Print();

  //Closure on the detector level
  TH1D* testPtRefoldedOverDetector = (TH1D*)testPtRefolded->Clone("testPtRefoldedOverDetector");
  testPtRefoldedOverDetector->Reset();
  testPtRefoldedOverDetector->Divide(testPtRefolded, testPtDetector);
  testPtRefoldedOverDetector->SetTitle(TString::Format("Refolded over Detector %s (test sample); #it{p}_{T}; #frac{Refolded}{Detector}",ratioNamePt.c_str()));
  testPtRefoldedOverDetector->Print();

  TH1D* testZRefoldedOverDetector = (TH1D*)testZRefolded->Clone("testZRefoldedOverDetector");
  testZRefoldedOverDetector->Reset();
  testZRefoldedOverDetector->Divide(testZRefolded, testZDetector);
  testZRefoldedOverDetector->SetTitle(TString::Format("Refolded over Detector %s (test sample); #it{z}; #frac{Refolded}{Detector}",ratioNameZ.c_str()));
  testZRefoldedOverDetector->Print();

  //Saving the output in a new file
  string outFileName = TString::Format("closureTest_pt%d-%d_ptmin%.0f_z%.0f-%.0f_nIter%d.root",
                                       PT_LOW, PT_HIGH, ptmin, zmin, zmax, N_ITER).Data();
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
  unfoldedTest->Write(unfoldedTest->GetName(), TObject::kOverwrite);
  refoldedTest->Write(refoldedTest->GetName(), TObject::kOverwrite);
  trainingPtDetector->Write(trainingPtDetector->GetName(), TObject::kOverwrite);
  trainingZDetector->Write(trainingZDetector->GetName(), TObject::kOverwrite);
  trainingPtTruth->Write(trainingPtTruth->GetName(), TObject::kOverwrite);
  trainingZTruth->Write(trainingZTruth->GetName(), TObject::kOverwrite);
  testPtDetector->Write(testPtDetector->GetName(), TObject::kOverwrite);
  testZDetector->Write(testZDetector->GetName(), TObject::kOverwrite);
  testPtTruth->Write(testPtTruth->GetName(), TObject::kOverwrite);
  testZTruth->Write(testZTruth->GetName(), TObject::kOverwrite);
  testPtUnfolded->Write(testPtUnfolded->GetName(), TObject::kOverwrite);
  testZUnfolded->Write(testZUnfolded->GetName(), TObject::kOverwrite);
  testPtRefolded->Write(testPtRefolded->GetName(), TObject::kOverwrite);
  testZRefolded->Write(testZRefolded->GetName(), TObject::kOverwrite);

  trainingPtTruthOverDetector->Write(trainingPtTruthOverDetector->GetName(), TObject::kOverwrite);
  trainingZTruthOverDetector->Write(trainingZTruthOverDetector->GetName(), TObject::kOverwrite);
  testPtUnfoldedOverDetector->Write(testPtUnfoldedOverDetector->GetName(), TObject::kOverwrite);
  testZUnfoldedOverDetector->Write(testZUnfoldedOverDetector->GetName(), TObject::kOverwrite);
  testPtUnfoldedOverTruth->Write(testPtUnfoldedOverTruth->GetName(), TObject::kOverwrite);
  testZUnfoldedOverTruth->Write(testZUnfoldedOverTruth->GetName(), TObject::kOverwrite);
  testPtRefoldedOverDetector->Write(testPtRefoldedOverDetector->GetName(), TObject::kOverwrite);
  testZRefoldedOverDetector->Write(testZRefoldedOverDetector->GetName(), TObject::kOverwrite);
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

TH2F* RM_normalization(TH2F* input_RM){
  Int_t nBinsDpt[2]= {input_RM->GetXaxis()->GetNbins(),input_RM->GetYaxis()->GetNbins()};
  TH2F* RM_norm =(TH2F*)input_RM->Clone();
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

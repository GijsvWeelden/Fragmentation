
// This task:
// * Loads the response from the task output
// * Creates the RooUnfoldResponse object used in unfolding

Int_t ptDetAxis   = 0;
Int_t zDetAxis    = 1;
Int_t ptTruthAxis = 2;
Int_t zTruthAxis  = 3;

void ProjectionRooUnfoldResponse(string inName = "AnalysisResults.root", Double_t PT_MIN = 20.0, Double_t Z_MIN = 0, Double_t Z_MAX = 1, Int_t BINWIDTHPT = 5)
{
  //Open the file and get the list where the response matrix is stored in THnSparse format
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  TDirectory* matchJetsDir = (TDirectory*)inFile->Get("jet-fragmentation/matching/jets");
  if(!matchJetsDir){
    std::cout << "Directory " << matchJetsDir << " not found. Aborting program." << std::endl;
    return;
  }
  THnSparseF* responseMatrix = (THnSparseF*)matchJetsDir->Get("matchDetJetPtTrackProjPartJetPtTrackProj");
  if(!responseMatrix){
    std::cout << "Response matrix " << responseMatrix << " not found. Aborting program." << std::endl;
    return;
  }
  responseMatrix->Print();

  // True and measured spectra as constructed from the response matrix
  // responseMatrix->GetAxis(ptDetAxis)->SetRange(
  //   3, responseMatrix->GetAxis(ptDetAxis)->GetNbins()
  // );
  // responseMatrix->GetAxis(ptTruthAxis)->SetRange(
  //   3, responseMatrix->GetAxis(ptTruthAxis)->GetNbins()
  // );
  TH2F* responseMatrixTruthProjection = (TH2F*)responseMatrix->Projection(zTruthAxis, ptTruthAxis, "E");
  responseMatrixTruthProjection->SetName("responseMatrixTruthProjection");
  TH2F* responseMatrixDetectorProjection = (TH2F*)responseMatrix->Projection(zDetAxis, ptDetAxis, "E");
  responseMatrixDetectorProjection->SetName("responseMatrixDetectorProjection");

  // Create RooUnfoldResponse and fill it
  responseMatrixDetectorProjection->Print();
  responseMatrixTruthProjection->Print();
  RooUnfoldResponse *ruResponse = new RooUnfoldResponse("Response", "RM");
  ruResponse->Setup(responseMatrixDetectorProjection, responseMatrixTruthProjection);

  Int_t* coord = new Int_t[responseMatrix->GetNdimensions()]; //Carries the bin coordinates
  for (int iBin = 0; iBin < responseMatrix->GetNbins(); iBin++) {
    Double_t w          = responseMatrix->GetBinContent(iBin, coord);
    Double_t ptTruth    = responseMatrix->GetAxis(ptTruthAxis)->GetBinCenter(coord[ptTruthAxis]);
    Double_t zTruth     = responseMatrix->GetAxis(zTruthAxis)->GetBinCenter(coord[zTruthAxis]);
    Double_t ptMeasured = responseMatrix->GetAxis(ptDetAxis)->GetBinCenter(coord[ptDetAxis]);
    Double_t zMeasured  = responseMatrix->GetAxis(zDetAxis)->GetBinCenter(coord[zDetAxis]);

    ruResponse->Fill(ptMeasured, zMeasured, ptTruth, zTruth, w);
  }
  delete [] coord;

  TFile* outFile = new TFile("RooUnfoldResponse.root", "RECREATE");
  responseMatrix->Write("responseMatrix");
  responseMatrixTruthProjection->Write();
  responseMatrixDetectorProjection->Write();
  ruResponse->Write();
  outFile->Write();
  outFile->Close();
}

void ProjectionClosureTest(string testFileName = "AnalysisResults.root", string responseFileName = "RooUnfoldResponse.root", Int_t PT_LOW = 20 , Int_t PT_HIGH = 80, Int_t N_ITER = 1)
{
  // Training RM
  TFile* responseFile = TFile::Open(TString::Format("%s", responseFileName.c_str()).Data());
  THnSparse* responseMatrix = static_cast<THnSparse*>(responseFile->Get("responseMatrix"));
  responseMatrix->SetName("responseMatrix");
  // responseMatrix->GetAxis(ptDetAxis)->SetRange(
  //   3, responseMatrix->GetAxis(ptDetAxis)->GetNbins()
  // );
  // responseMatrix->GetAxis(ptTruthAxis)->SetRange(
  //   3, responseMatrix->GetAxis(ptTruthAxis)->GetNbins()
  // );

  // Test RM
  TFile* testFile = TFile::Open(TString::Format("%s", testFileName.c_str()).Data());
  TDirectory* testDir = (TDirectory*)testFile->Get("jet-fragmentation/matching/jets");
  THnSparseF* testResponseMatrix = (THnSparseF*)testDir->Get("matchDetJetPtTrackProjPartJetPtTrackProj");
  // testResponseMatrix->GetAxis(ptDetAxis)->SetRange(
  //   3, testResponseMatrix->GetAxis(ptDetAxis)->GetNbins()
  // );
  // testResponseMatrix->GetAxis(ptTruthAxis)->SetRange(
  //   3, testResponseMatrix->GetAxis(ptTruthAxis)->GetNbins()
  // );
  TH2F *testResponseMatrixTruthProjection = (TH2F*)testResponseMatrix->Projection(zTruthAxis, ptTruthAxis, "E");
  testResponseMatrixTruthProjection->SetName("testResponseMatrixTruthProjection");
  TH2F *testResponseMatrixDetectorProjection = (TH2F*)testResponseMatrix->Projection(zDetAxis, ptDetAxis, "E");
  testResponseMatrixDetectorProjection->SetName("testResponseMatrixDetectorProjection");

  // This is the section where the actual unfolding starts and RooUnfold is involved
  string ruResponseName = "Response";
  string unfoldName = "Bayesian_Unfolding";
  string unfoldTitle = unfoldName;
  int doSmoothing = 0;
  // Load the response object from the Response file
  // Create the Bayesian unfolding object
  // Extract the unfolded distribution as a histogram
  // Apply the Response matrix to the unfolded result
  RooUnfoldResponse* ruResponse = static_cast<RooUnfoldResponse*>(responseFile->Get("Response"));
  RooUnfoldBayes Unfolding_bayes(ruResponse, testResponseMatrixDetectorProjection, N_ITER, doSmoothing, unfoldName.c_str(), unfoldTitle.c_str());
  TH2F* unfoldedTest = (TH2F*)Unfolding_bayes.Hreco(RooUnfold::kCovToy);
  unfoldedTest->SetName("unfoldedTest");
  auto Unf_Bayes = &Unfolding_bayes;
  TH2F* refoldedTest = (TH2F*)ruResponse->ApplyToTruth(unfoldedTest, "refoldedTest");
  // Include the fakes in order to compare with the input distribution
  // refoldedTest->Add(hDetectorFake);

  // The unfolding is complete, the following sections refer to projections in one dimension and calculating the ratio between distributions

  //--------------Training sample (Sample creating the response)--------------------//
    //Truth level
  TH1D* trainingPtTruth = (TH1D*)responseMatrix->Projection(ptTruthAxis);
  trainingPtTruth->SetName("trainingPtTruth");
  trainingPtTruth->SetTitle("Training truth pt");
  trainingPtTruth->Scale(1./trainingPtTruth->Integral(), "width");

  TH1D* trainingZTruth = (TH1D*)responseMatrix->Projection(zTruthAxis);
  trainingZTruth->SetName("trainingZTruth");
  trainingZTruth->SetTitle("Training truth z");
  trainingZTruth->Scale(1./trainingZTruth->Integral(), "width");
    //Detector level
  TH1D* trainingPtDetector = (TH1D*)responseMatrix->Projection(ptDetAxis);
  trainingPtDetector->SetName("trainingPtDetector");
  trainingPtDetector->SetTitle("Training detector pt");
  trainingPtDetector->Scale(1./trainingPtDetector->Integral(), "width");

  TH1D* trainingZDetector = (TH1D*)responseMatrix->Projection(zDetAxis);
  trainingZDetector->SetName("trainingZDetector");
  trainingZDetector->SetTitle("Training detector z");
  trainingZDetector->Scale(1./trainingZDetector->Integral(), "width");

  //--------------Test sample (MC sample treated as pseudodata)--------------------//
    //Truth level
  TH1D* testPtTruth = (TH1D*)testResponseMatrix->Projection(ptTruthAxis);
  testPtTruth->SetName("testPtTruth");
  testPtTruth->SetTitle("Test truth pt");
  testPtTruth->Scale(1./testPtTruth->Integral(), "width");

  TH1D* testZTruth = (TH1D*)testResponseMatrix->Projection(zTruthAxis);
  testZTruth->SetName("testZTruth");
  testZTruth->SetTitle("Test truth z");
  testZTruth->Scale(1./testZTruth->Integral(), "width");
    //Detector level
  TH1D* testPtDetector = (TH1D*)testResponseMatrix->Projection(ptDetAxis);
  testPtDetector->SetName("testPtDetector");
  testPtDetector->SetTitle("Test detector pt");
  testPtDetector->Scale(1./testPtDetector->Integral(), "width");

  TH1D* testZDetector = (TH1D*)testResponseMatrix->Projection(zDetAxis);
  testZDetector->SetName("testZDetector");
  testZDetector->SetTitle("Test truth z");
  testZDetector->Scale(1./testZDetector->Integral(), "width");
    //Unfolded
  TH1D* testPtUnfolded = unfoldedTest->ProjectionX("testPtUnfolded");
  testPtUnfolded->SetTitle("Unfolded pt");
  testPtUnfolded->Scale(1./testPtUnfolded->Integral(), "width");

  TH1D* testZUnfolded = unfoldedTest->ProjectionY("testZUnfolded");
  testZUnfolded->SetTitle("Unfolded z");
  testZUnfolded->Scale(1./testZUnfolded->Integral(), "width");
    //Refolded
  TH1D* testPtRefolded = refoldedTest->ProjectionX("testPtRefolded");
  testPtRefolded->SetTitle("Refolded pt");
  testPtRefolded->Scale(1./testPtRefolded->Integral(), "width");

  TH1D* testZRefolded = refoldedTest->ProjectionY("testZRefolded");
  testZRefolded->SetTitle("Refolded z");
  testZRefolded->Scale(1./testZRefolded->Integral(), "width");

  // Ratios of the projections
  // Training distributions - Detector response corrections
  TH1D* trainingPtTruthOverDetector = (TH1D*)trainingPtTruth->Clone("trainingPtTruthOverDetector");
  trainingPtTruthOverDetector->Divide(trainingPtDetector);
  trainingPtTruthOverDetector->SetTitle("Training pt truth/detector; #it{p}_{T}; #frac{Truth}{Detector}");
  trainingPtTruthOverDetector->Print();

  TH1D* trainingZTruthOverDetector = (TH1D*)trainingZTruth->Clone("trainingZTruthOverDetector");
  trainingZTruthOverDetector->Divide(trainingZDetector);
  trainingZTruthOverDetector->SetTitle("Training z truth/detector; #it{z}; #frac{Truth}{Detector}");
  trainingZTruthOverDetector->Print();

  //Test distributions - Detector response corrections
  TH1D* testPtUnfoldedOverDetector = (TH1D*)testPtUnfolded->Clone("testPtUnfoldedOverDetector");
  testPtUnfoldedOverDetector->Divide(testPtDetector);
  testPtUnfoldedOverDetector->SetTitle("Test pt unfolded/detector; #it{p}_{T}; #frac{Unfolded}{Detector}");
  testPtUnfoldedOverDetector->Print();

  TH1D* testZUnfoldedOverDetector = (TH1D*)testZUnfolded->Clone("testZUnfoldedOverDetector");
  testZUnfoldedOverDetector->Divide(testZDetector);
  testZUnfoldedOverDetector->SetTitle("Test z unfolded/detector; #it{z}; #frac{Unfolded}{Detector}");
  testZUnfoldedOverDetector->Print();

  //Closure on the truth level
  TH1D* testPtUnfoldedOverTruth = (TH1D*)testPtUnfolded->Clone("testPtUnfoldedOverTruth");
  testPtUnfoldedOverTruth->Divide(testPtTruth);
  testPtUnfoldedOverTruth->SetTitle("Test pt unfolded/truth; #it{p}_{T}; #frac{Unfolded}{Truth}");
  testPtUnfoldedOverTruth->Print();

  TH1D* testZUnfoldedOverTruth = (TH1D*)testZUnfolded->Clone("testZUnfoldedOverTruth");
  testZUnfoldedOverTruth->Divide(testZTruth);
  testZUnfoldedOverTruth->SetTitle("Test z unfolded/truth; #it{z}; #frac{Unfolded}{Truth}");
  testZUnfoldedOverTruth->Print();

  //Closure on the detector level
  TH1D* testPtRefoldedOverDetector = (TH1D*)testPtRefolded->Clone("testPtRefoldedOverDetector");
  testPtRefoldedOverDetector->Divide(testPtDetector);
  testPtRefoldedOverDetector->SetTitle("Test pt refolded/detector; #it{p}_{T}; #frac{Refolded}{Detector}");
  testPtRefoldedOverDetector->Print();

  TH1D* testZRefoldedOverDetector = (TH1D*)testZRefolded->Clone("testZRefoldedOverDetector");
  testZRefoldedOverDetector->Divide(testZDetector);
  testZRefoldedOverDetector->SetTitle("Test z unfolded/detector; #it{z}; #frac{Refolded}{Detector}");
  testZRefoldedOverDetector->Print();

  //Saving the output in a new file
  string outFileName = "closureTest_projection.root";
  TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");

  // Response objects
  responseMatrix->Write(responseMatrix->GetName(), TObject::kOverwrite);
  testResponseMatrix->Write(testResponseMatrix->GetName(), TObject::kOverwrite);
  testResponseMatrixTruthProjection->Write(testResponseMatrixTruthProjection->GetName(), TObject::kOverwrite);
  testResponseMatrixDetectorProjection->Write(testResponseMatrixDetectorProjection->GetName(), TObject::kOverwrite);
  ruResponse->Write(ruResponse->GetName(), TObject::kOverwrite);
  Unf_Bayes->Write(Unf_Bayes->GetName(), TObject::kOverwrite);

  trainingPtDetector->Write(trainingPtDetector->GetName(), TObject::kOverwrite);
  trainingZDetector->Write(trainingZDetector->GetName(), TObject::kOverwrite);
  trainingPtTruth->Write(trainingPtTruth->GetName(), TObject::kOverwrite);
  trainingZTruth->Write(trainingZTruth->GetName(), TObject::kOverwrite);

  testPtDetector->Write(testPtDetector->GetName(), TObject::kOverwrite);
  testZDetector->Write(testZDetector->GetName(), TObject::kOverwrite);
  testPtTruth->Write(testPtTruth->GetName(), TObject::kOverwrite);
  testZTruth->Write(testZTruth->GetName(), TObject::kOverwrite);

  unfoldedTest->Write(unfoldedTest->GetName(), TObject::kOverwrite);
  refoldedTest->Write(refoldedTest->GetName(), TObject::kOverwrite);

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

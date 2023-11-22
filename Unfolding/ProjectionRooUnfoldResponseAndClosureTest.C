
// This task:
// * Loads the response from the task output
// * Creates the RooUnfoldResponse object used in unfolding

using std::cout, std::endl;

Int_t ptDetAxis   = 0;
Int_t zDetAxis    = 1;
Int_t ptTruthAxis = 2;
Int_t zTruthAxis  = 3;

void CheckBinning(TH1* hist)
{
  cout << hist->GetName() << " "
    << hist->GetXaxis()->GetBinLowEdge(1) << " - " << hist->GetXaxis()->GetBinLowEdge(hist->GetNbinsX() + 1)
    << "(" << hist->GetXaxis()->GetBinWidth(1) << ")"
    << endl;
}
void CheckBinning(TH2* hist)
{
  cout << hist->GetName() << " "
    << hist->GetXaxis()->GetBinLowEdge(1) << " - " << hist->GetXaxis()->GetBinLowEdge(hist->GetNbinsX() + 1)
    << "(" << hist->GetXaxis()->GetBinWidth(1) << "), "
    << hist->GetYaxis()->GetBinLowEdge(1) << " - " << hist->GetYaxis()->GetBinLowEdge(hist->GetNbinsY() + 1)
    << "(" << hist->GetYaxis()->GetBinWidth(1) << ")"
    << endl;
}
void CompareHists(TH1* histA, TH1* histB)
{
  cout << histA->GetName() << " - " << histB->GetName() << endl;
  for (int i = 0; i < histA->GetNbinsX(); i++)
  {
    double binDiff = histA->GetBinContent(i) - histB->GetBinContent(i);
    if ( abs(binDiff) > 0.) cout << "Bin " << i << " = " << binDiff << " ";
  }
  cout << endl;
}
void CompareHists(TH2* histA, TH2* histB)
{
  cout << histA->GetName() << " - " << histB->GetName() << endl;
  for (int i = 0; i < histA->GetNbinsX(); i++)
  {
    for (int j = 0; j < histA->GetNbinsY(); j++) {
      double binDiff = histA->GetBinContent(i, j) - histB->GetBinContent(i, j);
      if ( abs(binDiff) > 0.) cout << "Bin (" << i << ", " << j << ") = " << binDiff << " ";
    }
  }
  cout << endl;
}
void NormaliseHist(TH2* hist, bool rowbyrow = false)
{
  int firstColBin = 1, lastColBin = hist->GetNbinsX();
  int firstRowBin = 1, lastRowBin = hist->GetNbinsY();

  if (rowbyrow) {
    for (int iRow = 1; iRow <= lastRowBin; iRow++) {
      double integral = hist->Integral(firstColBin, lastColBin, iRow, iRow);
      if (integral < 1e-9) { continue; } // Protection against divide by 0
      for (int iCol = 1; iCol <= lastColBin; iCol++) {
        double binContent = hist->GetBinContent(iCol, iRow);
        binContent /= integral;
        hist->SetBinContent(iCol, iRow, binContent);
      }
    }
  } else {
    for (int iCol = 1; iCol <= lastColBin; iCol++) {
      double integral = hist->Integral(firstRowBin, lastRowBin, iCol, iCol);
      if (integral < 1e-9) { continue; } // Protection against divide by 0
      for (int iRow = 1; iRow <= lastRowBin; iRow++) {
        double binContent = hist->GetBinContent(iCol, iRow);
        binContent /= integral;
        hist->SetBinContent(iCol, iRow, binContent);
      }
    }
  }
}
void NormaliseHistColByCol(TH2* hist)
{
  NormaliseHist(hist, false);
}
void NormaliseHistRowByRow(TH2* hist)
{
  NormaliseHist(hist, true);
}
void THnMinMax(THnSparse* hist)
{
  int minBin = 0, maxBin = 0;
  double minBinContent = 9e10, maxBinContent = -9e10;
  for (int i = 1; i <= hist->GetNbins(); i++)
  {
    double binContent = hist->GetBinContent(i);
    if (binContent < minBinContent && binContent > 0) {
      minBinContent = binContent;
      minBin = i;
    }
    if (binContent > maxBinContent) {
      maxBinContent = binContent;
      maxBin = i;
    }
  }
  cout << "Min value (!0)" << minBinContent << " in bin " << minBin << endl;
  cout << "Max value " << maxBinContent << " in bin " << maxBin << endl;
}

// Checks whether what RooUnfold considers measured/truth actually corresponds to the measured/truth distributions given as input
void CheckRooUnfoldInput(string inName = "RooUnfoldResponse.root", bool verbose = true)
{
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  // int ptTruthMin = 10, ptTruthMax = 110;
  // int ptDetectorMin = 10, ptDetectorMax = 110;
  TH2D* rmTruth = (TH2D*)inFile->Get("rmTruth");
  TH1F* rmTruthPt = (TH1F*)rmTruth->ProjectionX("rmTruthPt");
  TH1F* rmTruthZ = (TH1F*)rmTruth->ProjectionY("rmTruthZ");
  // TH1F* rmTruthZ = (TH1F*)rmTruth->ProjectionY("rmTruthZ", rmTruth->GetYaxis()->FindBin(ptTruthMin), rmTruth->GetYaxis()->FindBin(ptTruthMax));
  TH2D* rmDet = (TH2D*)inFile->Get("rmDet");
  TH1F* rmDetPt = (TH1F*)rmDet->ProjectionX("rmDetPt");
  TH1F* rmDetZ = (TH1F*)rmDet->ProjectionY("rmDetZ");
  // TH1F* rmMeasuredZ = (TH1F*)rmMeasured->ProjectionY("rmMeasuredZ", rmMeasured->GetYaxis()->FindBin(ptDetectorMin), rmMeasured->GetYaxis()->FindBin(ptDetectorMin));

  RooUnfoldResponse* response = (RooUnfoldResponse*)inFile->Get("Response");
  TH2D* ruMeasured = (TH2D*)response->Hmeasured();
  TH1D* ruMeasuredPt = (TH1D*)ruMeasured->ProjectionX("ruMeasuredPt");
  TH1D* ruMeasuredZ = (TH1D*)ruMeasured->ProjectionY("ruMeasuredZ");
  // TH1D* ruMeasuredZ = (TH1D*)ruMeasured->ProjectionY("ruMeasuredZ", ruMeasured->GetYaxis()->FindBin(ptDetectorMin), ruMeasured->GetYaxis()->FindBin(ptDetectorMax));

  TH2D* ruTruth = (TH2D*)response->Htruth();
  TH1D* ruTruthPt = (TH1D*)ruTruth->ProjectionX("ruTruthPt");
  TH1D* ruTruthZ = (TH1D*)ruTruth->ProjectionY("ruTruthZ");
  // TH1D* ruTruthZ = (TH1D*)ruTruth->ProjectionY("ruTruthZ", ruTruth->GetYaxis()->FindBin(ptTruthMin), ruTruth->GetYaxis()->FindBin(ptTruthMax));

  TH2D* ruApplied = (TH2D*)response->ApplyToTruth();
  TH2D* hFake = (TH2D*)inFile->Get("hFake");
  ruApplied->Add(hFake);
  TH1D* ruAppliedPt = (TH1D*)ruApplied->ProjectionX("ruAppliedPt");
  TH1D* ruAppliedZ = (TH1D*)ruApplied->ProjectionY("ruAppliedZ");
  // TH1D* ruAppliedZ = (TH1D*)ruApplied->ProjectionY("ruAppliedZ", ruApplied->GetYaxis()->FindBin(ptDetectorMin), ruApplied->GetYaxis()->FindBin(ptDetectorMax));

  TH1D* ptMeasuredOverResponse = (TH1D*)rmDetPt->Clone("ptMeasuredOverResponse");
  ptMeasuredOverResponse->Divide(ruMeasuredPt);
  ptMeasuredOverResponse->SetTitle("RM Measured / RU Measured");
  TH1D* zMeasuredOverResponse = (TH1D*)rmDetZ->Clone("zMeasuredOverResponse");
  zMeasuredOverResponse->Divide(ruMeasuredZ);
  zMeasuredOverResponse->SetTitle("RM Measured / RU Measured");
  TH1D* ptTruthOverResponse = (TH1D*)rmTruthPt->Clone("ptTruthOverResponse");
  ptTruthOverResponse->Divide(ruTruthPt);
  ptTruthOverResponse->SetTitle("RM Truth / RU Truth");
  TH1D* zTruthOverResponse = (TH1D*)rmTruthZ->Clone("zTruthOverResponse");
  zTruthOverResponse->Divide(ruTruthZ);
  zTruthOverResponse->SetTitle("RM Truth / RU Truth");
  TH1D* ptMeasuredOverApplied = (TH1D*)rmDetPt->Clone("ptMeasuredOverApplied");
  ptMeasuredOverApplied->Divide(ruAppliedPt);
  ptMeasuredOverApplied->SetTitle("RM Measured / RU Applied");
  TH1D* zMeasuredOverApplied = (TH1D*)rmDetZ->Clone("zMeasuredOverApplied");
  zMeasuredOverApplied->Divide(ruAppliedZ);
  zMeasuredOverApplied->SetTitle("RM Measured / RU Applied");

  TCanvas* myCanvas = new TCanvas("myCanvas", "myCanvas", 800, 1200);
  myCanvas->Divide(2, 3);

  myCanvas->cd(1);
  ptMeasuredOverResponse->Draw();
  myCanvas->cd(2);
  zMeasuredOverResponse->Draw();
  myCanvas->cd(3);
  ptTruthOverResponse->Draw();
  myCanvas->cd(4);
  zTruthOverResponse->Draw();
  myCanvas->cd(5);
  ptMeasuredOverApplied->Draw();
  myCanvas->cd(6);
  zMeasuredOverApplied->Draw();

  myCanvas->SaveAs("./checkRooUnfold.pdf");

  if (verbose) {
    cout << "Comparing histograms bin-by-bin" << endl;
    CompareHists(rmDetPt, ruMeasuredPt);
    CompareHists(rmDetZ, ruMeasuredZ);
    CompareHists(rmTruthPt, ruTruthPt);
    CompareHists(rmTruthZ, ruTruthZ);
    CompareHists(rmDetPt, ruAppliedPt);
    CompareHists(rmDetZ, ruAppliedZ);
  }
}

void CreateRooUnfoldResponse(string inName = "AnalysisResults.root", bool makePlots = false)
//, Double_t PT_MIN = 20.0, Double_t Z_MIN = 0, Double_t Z_MAX = 1, Int_t BINWIDTHPT = 5)
{
  //Open the file and get the response matrix in THnSparse format
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  THnSparseF* responseMatrix = (THnSparseF*)inFile->Get("jet-fragmentation/matching/jets/matchDetJetPtTrackProjPartJetPtTrackProj");
  if(!responseMatrix){
    std::cout << "Response matrix " << responseMatrix << " not found. Aborting program." << std::endl;
    return;
  }
  responseMatrix->Print();

  // !!! Make sure that all entries are >~ 1. Only for tests !!!
  THnMinMax(responseMatrix);
  responseMatrix->Scale(1e10);
  THnMinMax(responseMatrix);

  THnSparseF* tmpRM = (THnSparseF*)responseMatrix->Clone("tmpRM");
  double ptTruthMin = 10, ptTruthMax = 300;
  double ptDetMin = 10, ptDetMax = 300;
  double zTruthMin = 0, zTruthMax = 1;
  double zDetMin = 0, zDetMax = 1;
  // double ptTruthMin = 10, ptTruthMax = 300;
  // double ptDetMin = 10, ptDetMax = 300;
  // int rebinX = 2, rebinY = 2;
  int ptTruthMinBin = tmpRM->GetAxis(ptTruthAxis)->FindBin(ptTruthMin + 0.1);
  int ptTruthMaxBin = tmpRM->GetAxis(ptTruthAxis)->FindBin(ptTruthMax + 0.1);
  int zTruthMinBin  = tmpRM->GetAxis(zTruthAxis)->FindBin(zTruthMin + 1e-3);
  int zTruthMaxBin  = tmpRM->GetAxis(zTruthAxis)->FindBin(zTruthMax + 1e-3);

  int ptDetMinBin   = tmpRM->GetAxis(ptDetAxis)->FindBin(ptDetMin + 0.1);
  int ptDetMaxBin   = tmpRM->GetAxis(ptDetAxis)->FindBin(ptDetMax + 0.1);
  int zDetMinBin    = tmpRM->GetAxis(zDetAxis)->FindBin(zDetMin + 1e-3);
  int zDetMaxBin    = tmpRM->GetAxis(zDetAxis)->FindBin(zDetMax + 1e-3);

  tmpRM->GetAxis(ptTruthAxis)->SetRange(ptTruthMinBin, ptTruthMaxBin - 1);
  tmpRM->GetAxis(ptDetAxis)->SetRange(ptDetMinBin, ptDetMaxBin - 1);
  tmpRM->GetAxis(zTruthAxis)->SetRange(zTruthMinBin, zTruthMaxBin - 1);
  tmpRM->GetAxis(zDetAxis)->SetRange(zDetMinBin, zDetMaxBin - 1);

  TH2F* rmTruth = (TH2F*)tmpRM->Projection(zTruthAxis, ptTruthAxis, "E");
  rmTruth->Reset();
  rmTruth->SetName("rmTruth");
  // rmTruth->Rebin2D(rebinX, rebinY);
  rmTruth->Print();
  TH2F* rmDet = (TH2F*)tmpRM->Projection(zDetAxis, ptDetAxis, "E");
  rmDet->Reset();
  rmDet->SetName("rmDet");
  // rmDet->Rebin2D(rebinX, rebinY);
  rmDet->Print();

  TH2F* hMiss = new TH2F("hMiss", "hMiss",
                         rmTruth->GetNbinsX(),
                         rmTruth->GetXaxis()->GetBinLowEdge(1),
                         rmTruth->GetXaxis()->GetBinLowEdge(rmTruth->GetNbinsX() + 1),
                         rmTruth->GetNbinsY(),
                         rmTruth->GetYaxis()->GetBinLowEdge(1),
                         rmTruth->GetYaxis()->GetBinLowEdge(rmTruth->GetNbinsY() + 1)
                        );
  TH2F* hFake = new TH2F("hFake", "hFake",
                         rmDet->GetNbinsX(),
                         rmDet->GetXaxis()->GetBinLowEdge(1),
                         rmDet->GetXaxis()->GetBinLowEdge(rmDet->GetNbinsX() + 1),
                         rmDet->GetNbinsY(),
                         rmDet->GetYaxis()->GetBinLowEdge(1),
                         rmDet->GetYaxis()->GetBinLowEdge(rmDet->GetNbinsY() + 1)
                        );

  bool checkBinning = true;
  if (checkBinning) {
    CheckBinning(rmTruth);
    CheckBinning(hMiss);
    CheckBinning(rmDet);
    CheckBinning(hFake);
  }

  // Create RooUnfoldResponse and fill it
  RooUnfoldResponse *ruResponse = new RooUnfoldResponse("Response", "RM");
  ruResponse->Setup(rmDet, rmTruth);

  Int_t* coord = new Int_t[responseMatrix->GetNdimensions()]; //Carries the bin coordinates
  for (int iBin = 0; iBin < responseMatrix->GetNbins(); iBin++) {
    Double_t w          = responseMatrix->GetBinContent(iBin, coord);
    Double_t ptTruth    = responseMatrix->GetAxis(ptTruthAxis)->GetBinCenter(coord[ptTruthAxis]);
    Double_t zTruth     = responseMatrix->GetAxis(zTruthAxis)->GetBinCenter(coord[zTruthAxis]);
    Double_t ptMeasured = responseMatrix->GetAxis(ptDetAxis)->GetBinCenter(coord[ptDetAxis]);
    Double_t zMeasured  = responseMatrix->GetAxis(zDetAxis)->GetBinCenter(coord[zDetAxis]);

    bool inAcceptanceTruth = false;
    bool inAcceptanceDetector = false;

    inAcceptanceTruth = (ptTruth >= ptTruthMin) * (ptTruth < ptTruthMax);
    inAcceptanceDetector = (ptMeasured >= ptDetMin) * (ptMeasured < ptDetMax);

    if (inAcceptanceTruth) {
      rmTruth->Fill(ptTruth, zTruth, w);
    }
    if (inAcceptanceDetector) {
      rmDet->Fill(ptMeasured, zMeasured, w);
    }

    if (inAcceptanceTruth && inAcceptanceDetector) {
      ruResponse->Fill(ptMeasured, zMeasured, ptTruth, zTruth, w);
    }
    else if (inAcceptanceTruth && !inAcceptanceDetector) {
      ruResponse->Miss(ptTruth, zTruth, w);
      hMiss->Fill(ptTruth, zTruth, w);
      // cout << "Miss: ptT=" << ptTruth << ", zT=" << zTruth << ", ptM=" << ptMeasured << ", zM=" << zMeasured << endl;
    }
    else if (!inAcceptanceTruth && inAcceptanceDetector) {
      ruResponse->Fake(ptMeasured, zMeasured, w);
      hFake->Fill(ptMeasured, zMeasured, w);
      // cout << "Fake: ptT=" << ptTruth << ", zT=" << zTruth << ", ptM=" << ptMeasured << ", zM=" << zMeasured << endl;
    }
  }
  delete [] coord;

  TFile* outFile = new TFile("RooUnfoldResponse.root", "RECREATE");
  responseMatrix->Write("responseMatrix");
  rmTruth->Write();
  rmDet->Write();
  hMiss->Write();
  hFake->Write();
  ruResponse->Write();
  outFile->Write();
  outFile->Close();

  if (!makePlots) return;

  rmTruth->SetStats(0);
  rmTruth->SetTitle("RM truth");
  rmDet->SetStats(0);
  rmDet->SetTitle("RM detector");

  TH1D* rmTruthPt = rmTruth->ProjectionX("rmTruthPt");
  rmTruthPt->SetTitle("Truth pt");
  rmTruthPt->SetStats(0);
  TH1D* rmTruthZ = rmTruth->ProjectionY("rmTruthZ");
  rmTruthZ->SetTitle("Truth z");
  rmTruthZ->SetStats(0);
  TH1D* rmDetectorPt = rmDet->ProjectionX("rmDetectorPt");
  rmDetectorPt->SetTitle("Detector pt");
  rmDetectorPt->SetStats(0);
  TH1D* rmDetectorZ = rmDet->ProjectionY("rmDetectorZ");
  rmDetectorZ->SetTitle("Detector z");
  rmDetectorZ->SetStats(0);

  // Plot the distributions
  TCanvas* myCanvas = new TCanvas("myCanvas", "myCanvas", 800, 1200);
  myCanvas->Divide(2, 3);
  auto pad1 = myCanvas->cd(1);
  pad1->SetLogz();
  NormaliseHistColByCol(rmTruth);
  rmTruth->SetMinimum(1e-10);
  rmTruth->SetMaximum(1e10);
  rmTruth->Draw("colz");

  auto pad2 = myCanvas->cd(2);
  pad2->SetLogz();
  NormaliseHistColByCol(rmDet);
  rmDet->SetMinimum(1e-10);
  rmDet->SetMaximum(1e10);
  rmDet->Draw("colz");

  auto pad3 = myCanvas->cd(3);
  pad3->SetLogy();
  rmTruthPt->Scale(1./rmTruthPt->Integral());
  rmTruthPt->Draw();

  auto pad4 = myCanvas->cd(4);
  pad4->SetLogy();
  rmTruthZ->Scale(1./rmTruthZ->Integral());
  rmTruthZ->Draw();

  auto pad5 = myCanvas->cd(5);
  pad5->SetLogy();
  rmDetectorPt->Scale(1./rmDetectorPt->Integral());
  rmDetectorPt->Draw();

  auto pad6 = myCanvas->cd(6);
  pad6->SetLogy();
  rmDetectorZ->Scale(1./rmDetectorZ->Integral());
  rmDetectorZ->Draw();

  myCanvas->SaveAs("./RM_pt-z.pdf");
}

void ClosureTest(string testFileName = "AnalysisResults.root", string responseFileName = "RooUnfoldResponse.root", int nIter = 3)//, Int_t PT_LOW = 20 , Int_t PT_HIGH = 80, Int_t N_ITER = 1)
{
  // Training RM
  TFile* responseFile = TFile::Open(TString::Format("%s", responseFileName.c_str()).Data());
  THnSparse* responseMatrix = static_cast<THnSparse*>(responseFile->Get("responseMatrix"));
  responseMatrix->SetName("responseMatrix");
  TH2F* responseMatrixTruthProjection = (TH2F*)responseFile->Get("responseMatrixTruthProjection");
  TH2F* responseMatrixDetectorProjection = (TH2F*)responseFile->Get("responseMatrixDetectorProjection");

  // Test RM
  TFile* testFile = TFile::Open(TString::Format("%s", testFileName.c_str()).Data());
  THnSparseF* testResponseMatrix = (THnSparseF*)testFile->Get("jet-fragmentation/matching/jets/matchDetJetPtTrackProjPartJetPtTrackProj");

  // !!! These ranges have to be same as the ones used in CreateRooUnfoldResponse() !!!
  double ptTruthMin    = responseMatrixTruthProjection->GetXaxis()->GetXmin();
  double ptTruthMax    = responseMatrixTruthProjection->GetXaxis()->GetXmax();
  double ptDetectorMin = responseMatrixDetectorProjection->GetXaxis()->GetXmin();
  double ptDetectorMax = responseMatrixDetectorProjection->GetXaxis()->GetXmax();
  double zTruthMin    = responseMatrixTruthProjection->GetYaxis()->GetXmin();
  double zTruthMax    = responseMatrixTruthProjection->GetYaxis()->GetXmax();
  double zDetectorMin = responseMatrixDetectorProjection->GetYaxis()->GetXmin();
  double zDetectorMax = responseMatrixDetectorProjection->GetYaxis()->GetXmax();
  // double ptTruthMin = 10, ptTruthMax = 300;
  // double ptDetectorMin = 10, ptDetectorMax = 300;
  testResponseMatrix->GetAxis(ptTruthAxis)->SetRange(testResponseMatrix->GetAxis(ptTruthAxis)->FindBin(ptTruthMin), testResponseMatrix->GetAxis(ptTruthAxis)->FindBin(ptTruthMax) - 1); // 50 bins
  testResponseMatrix->GetAxis(ptDetAxis)->SetRange(testResponseMatrix->GetAxis(ptDetAxis)->FindBin(ptDetectorMin), testResponseMatrix->GetAxis(ptDetAxis)->FindBin(ptDetectorMax) - 1); // 42 bins
  testResponseMatrix->GetAxis(zTruthAxis)->SetRange(testResponseMatrix->GetAxis(zTruthAxis)->FindBin(zTruthMin), testResponseMatrix->GetAxis(zTruthAxis)->FindBin(zTruthMax) - 1); // 50 bins
  testResponseMatrix->GetAxis(zDetAxis)->SetRange(testResponseMatrix->GetAxis(zDetAxis)->FindBin(zDetectorMin), testResponseMatrix->GetAxis(zDetAxis)->FindBin(zDetectorMax) - 1); // 42 bins

  TH2F *testResponseMatrixTruthProjection = (TH2F*)testResponseMatrix->Projection(zTruthAxis, ptTruthAxis, "E");
  testResponseMatrixTruthProjection->SetName("testResponseMatrixTruthProjection");
  TH2F *testResponseMatrixDetectorProjection = (TH2F*)testResponseMatrix->Projection(zDetAxis, ptDetAxis, "E");
  testResponseMatrixDetectorProjection->SetName("testResponseMatrixDetectorProjection");

  // This is the section where the actual unfolding starts and RooUnfold is involved
  string ruResponseName = "Response";
  string unfoldName = "Bayesian_Unfolding";
  string unfoldTitle = unfoldName;
  int doSmoothing = 0;
  // Load the response object from the Response file, create the Bayesian unfolding object,
  // Extract the unfolded distribution as a histogram, apply the Response matrix to the unfolded result
  RooUnfoldResponse* ruResponse = static_cast<RooUnfoldResponse*>(responseFile->Get("Response"));
  RooUnfoldBayes Unfolding_bayes(ruResponse, testResponseMatrixDetectorProjection, nIter, doSmoothing, unfoldName.c_str(), unfoldTitle.c_str());
  TH2F* unfoldedTest = (TH2F*)Unfolding_bayes.Hreco(RooUnfold::kCovToy);
  unfoldedTest->SetName("unfoldedTest");
  auto Unf_Bayes = &Unfolding_bayes;
  TH2F* refoldedTest = (TH2F*)ruResponse->ApplyToTruth(unfoldedTest, "refoldedTest");
  // Include the fakes in order to compare with the input distribution
  TH2F* hFake = (TH2F*)responseFile->Get("hFake");
  refoldedTest->Add(hFake);

  // The unfolding is complete, the following sections refer to projections in one dimension and calculating the ratio between distributions
  // !!! These ranges can be whatever the fuck you like !!!
  int PTLOW = 40, PTHIGH = 60;
  // int PTLOW = 10, PTHIGH = 300;
  int ptLowTruth = responseMatrixTruthProjection->GetXaxis()->FindBin(PTLOW);
  int ptHighTruth = responseMatrixTruthProjection->GetXaxis()->FindBin(PTHIGH);
  int ptLowDetector = responseMatrixDetectorProjection->GetXaxis()->FindBin(PTLOW);
  int ptHighDetector = responseMatrixDetectorProjection->GetXaxis()->FindBin(PTHIGH);

  //--------------Training sample (Sample creating the response)--------------------//
    //Truth level
  TH1D* trainingPtTruth = (TH1D*)responseMatrixTruthProjection->ProjectionX("trainingPtTruth");
  // trainingPtTruth->SetName("trainingPtTruth");
  trainingPtTruth->SetTitle("Training truth pt");
  // trainingPtTruth->Scale(1./trainingPtTruth->Integral(), "width");

  TH1D* trainingZTruth = (TH1D*)responseMatrixTruthProjection->ProjectionY("trainingZTruth", ptLowTruth, ptHighTruth);
  // trainingZTruth->SetName("trainingZTruth");
  trainingZTruth->SetTitle("Training truth z");
  // trainingZTruth->Scale(1./trainingZTruth->Integral(), "width");
    //Detector level
  TH1D* trainingPtDetector = (TH1D*)responseMatrixDetectorProjection->ProjectionX("trainingPtDetector");
  // trainingPtDetector->SetName("trainingPtDetector");
  trainingPtDetector->SetTitle("Training detector pt");
  // trainingPtDetector->Scale(1./trainingPtDetector->Integral(), "width");

  TH1D* trainingZDetector = (TH1D*)responseMatrixDetectorProjection->ProjectionY("trainingZDetector", ptLowDetector, ptHighDetector);
  // trainingZDetector->SetName("trainingZDetector");
  trainingZDetector->SetTitle("Training detector z");
  // trainingZDetector->Scale(1./trainingZDetector->Integral(), "width");

  //--------------Test sample (MC sample treated as pseudodata)--------------------//
    //Truth level
  TH1D* testPtTruth = (TH1D*)testResponseMatrixTruthProjection->ProjectionX("testPtTruth");
  // testPtTruth->SetName("testPtTruth");
  testPtTruth->SetTitle("Test truth pt");
  // testPtTruth->Scale(1./testPtTruth->Integral(), "width");

  TH1D* testZTruth = (TH1D*)testResponseMatrixTruthProjection->ProjectionY("testZTruth", ptLowTruth, ptHighTruth);
  // testZTruth->SetName("testZTruth");
  testZTruth->SetTitle("Test truth z");
  // testZTruth->Scale(1./testZTruth->Integral(), "width");
    //Detector level
  TH1D* testPtDetector = (TH1D*)testResponseMatrixDetectorProjection->ProjectionX("testPtDetector");
  // testPtDetector->SetName("testPtDetector");
  testPtDetector->SetTitle("Test detector pt");
  // testPtDetector->Scale(1./testPtDetector->Integral(), "width");

  TH1D* testZDetector = (TH1D*)testResponseMatrixDetectorProjection->ProjectionY("testZDetector", ptLowDetector, ptHighDetector);
  // testZDetector->SetName("testZDetector");
  testZDetector->SetTitle("Test truth z");
  // testZDetector->Scale(1./testZDetector->Integral(), "width");
    //Unfolded
  TH1D* testPtUnfolded = unfoldedTest->ProjectionX("testPtUnfolded");
  testPtUnfolded->SetTitle("Unfolded pt");
  // testPtUnfolded->Scale(1./testPtUnfolded->Integral(), "width");

  TH1D* testZUnfolded = unfoldedTest->ProjectionY("testZUnfolded", ptLowTruth, ptHighTruth);
  testZUnfolded->SetTitle("Unfolded z");
  // testZUnfolded->Scale(1./testZUnfolded->Integral(), "width");
    //Refolded
  TH1D* testPtRefolded = refoldedTest->ProjectionX("testPtRefolded");
  testPtRefolded->SetTitle("Refolded pt");
  // testPtRefolded->Scale(1./testPtRefolded->Integral(), "width");

  TH1D* testZRefolded = refoldedTest->ProjectionY("testZRefolded", ptLowDetector, ptHighDetector);
  testZRefolded->SetTitle("Refolded z");
  // testZRefolded->Scale(1./testZRefolded->Integral(), "width");

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
  // testPtUnfoldedOverDetector->Rebin(2);
  testPtUnfoldedOverDetector->Divide(testPtDetector);
  testPtUnfoldedOverDetector->SetTitle("Test pt unfolded/detector; #it{p}_{T}; #frac{Unfolded}{Detector}");
  testPtUnfoldedOverDetector->Print();

  TH1D* testZUnfoldedOverDetector = (TH1D*)testZUnfolded->Clone("testZUnfoldedOverDetector");
  // testZUnfoldedOverDetector->Rebin(2);
  testZUnfoldedOverDetector->Divide(testZDetector);
  testZUnfoldedOverDetector->SetTitle("Test z unfolded/detector; #it{z}; #frac{Unfolded}{Detector}");
  testZUnfoldedOverDetector->Print();

  //Closure on the truth level
  TH1D* testPtUnfoldedOverTruth = (TH1D*)testPtUnfolded->Clone("testPtUnfoldedOverTruth");
  testPtTruth->Rebin(2);
  testPtUnfoldedOverTruth->Divide(testPtTruth);
  testPtUnfoldedOverTruth->SetTitle("Test pt unfolded/truth; #it{p}_{T}; #frac{Unfolded}{Truth}");
  testPtUnfoldedOverTruth->Print();
  CompareHists(testPtUnfolded, testPtTruth);

  TH1D* testZUnfoldedOverTruth = (TH1D*)testZUnfolded->Clone("testZUnfoldedOverTruth");
  testZTruth->Rebin(2);
  testZUnfoldedOverTruth->Divide(testZTruth);
  testZUnfoldedOverTruth->SetTitle("Test z unfolded/truth; #it{z}; #frac{Unfolded}{Truth}");
  testZUnfoldedOverTruth->Print();
  CompareHists(testZUnfolded, testZTruth);

  //Closure on the detector level
  TH1D* testPtRefoldedOverDetector = (TH1D*)testPtRefolded->Clone("testPtRefoldedOverDetector");
  testPtRefoldedOverDetector->Divide(testPtDetector);
  testPtRefoldedOverDetector->SetTitle("Test pt refolded/detector; #it{p}_{T}; #frac{Refolded}{Detector}");
  testPtRefoldedOverDetector->Print();
  CompareHists(testPtRefolded, testPtDetector);

  TH1D* testZRefoldedOverDetector = (TH1D*)testZRefolded->Clone("testZRefoldedOverDetector");
  testZRefoldedOverDetector->Divide(testZDetector);
  testZRefoldedOverDetector->SetTitle("Test z unfolded/detector; #it{z}; #frac{Refolded}{Detector}");
  testZRefoldedOverDetector->Print();
  CompareHists(testZRefolded, testZDetector);

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

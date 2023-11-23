
// This task:
// * Loads the response from the task output
// * Creates the RooUnfoldResponse object used in unfolding

using std::cout, std::endl;

int ptDetAxis   = 0;
int zDetAxis    = 1;
int ptTruthAxis = 2;
int zTruthAxis  = 3;

void CheckBinning(TH1* hist)
{
  cout << hist->GetName() << " "
    << hist->GetXaxis()->GetBinLowEdge(1) << " - " << hist->GetXaxis()->GetBinLowEdge(hist->GetNbinsX() + 1)
    << " (" << hist->GetNbinsX() << " * " << hist->GetXaxis()->GetBinWidth(1) << ")"
    << endl;
}
void CheckBinning(TH2* hist)
{
  cout << hist->GetName() << " "
    << hist->GetXaxis()->GetBinLowEdge(1) << " - " << hist->GetXaxis()->GetBinLowEdge(hist->GetNbinsX() + 1)
    << " (" << hist->GetNbinsX() << " * " << hist->GetXaxis()->GetBinWidth(1) << "), "
    << hist->GetYaxis()->GetBinLowEdge(1) << " - " << hist->GetYaxis()->GetBinLowEdge(hist->GetNbinsY() + 1)
    << " (" << hist->GetNbinsY() << " * " << hist->GetYaxis()->GetBinWidth(1) << ")"
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
TH1* HistDiff(TH1* histA, TH1* histB, string outName = "outHist", string outTitle = "")
{
  TH1* outHist = (TH1*)histA->Clone(outName.c_str());
  outHist->SetTitle(outTitle.c_str());
  outHist->Add(histB, -1.);
  return outHist;
}
TH1* HistRatio(TH1* histA, TH1* histB, string outName = "outHist", string outTitle = "")
{
  TH1* outHist = (TH1*)histA->Clone(outName.c_str());
  outHist->SetTitle(outTitle.c_str());
  outHist->Divide(histB);
  return outHist;
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
  double minBinContent = 9e12, maxBinContent = -9e12;
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
  cout << "Min value (!0) " << minBinContent << " in bin " << minBin << "\n";
  cout << "Max value " << maxBinContent << " in bin " << maxBin << endl;
}

// Checks whether what RooUnfold considers measured/truth actually corresponds to the measured/truth distributions given as input
void CheckRooUnfoldInput(string inName = "RooUnfoldResponse.root")
{
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  TH2D* hKinEff = (TH2D*)inFile->Get("hKinEff");

  TH2D* rmTruth = (TH2D*)inFile->Get("rmTruth");
  TH1F* rmTruthPt = (TH1F*)rmTruth->ProjectionX("rmTruthPt");
  TH1F* rmTruthZ = (TH1F*)rmTruth->ProjectionY("rmTruthZ");

  TH2D* rmDet = (TH2D*)inFile->Get("rmDet");
  TH1F* rmDetPt = (TH1F*)rmDet->ProjectionX("rmDetPt");
  TH1F* rmDetZ = (TH1F*)rmDet->ProjectionY("rmDetZ");

  RooUnfoldResponse* response = (RooUnfoldResponse*)inFile->Get("Response");
  TH2D* ruMeasured = (TH2D*)response->Hmeasured();
  TH1D* ruMeasuredPt = (TH1D*)ruMeasured->ProjectionX("ruMeasuredPt");
  TH1D* ruMeasuredZ = (TH1D*)ruMeasured->ProjectionY("ruMeasuredZ");

  TH2D* ruTruth = (TH2D*)response->Htruth();
  TH1D* ruTruthPt = (TH1D*)ruTruth->ProjectionX("ruTruthPt");
  TH1D* ruTruthZ = (TH1D*)ruTruth->ProjectionY("ruTruthZ");

  TH2D* ruApplied = (TH2D*)response->ApplyToTruth();
  TH2D* hFake = (TH2D*)inFile->Get("hFake");
  ruApplied->Add(hFake);
  TH1D* ruAppliedPt = (TH1D*)ruApplied->ProjectionX("ruAppliedPt");
  TH1D* ruAppliedZ = (TH1D*)ruApplied->ProjectionY("ruAppliedZ");

  TCanvas* spectraCanvas = new TCanvas("spectraCanvas", "spectraCanvas", 1600, 1200);
  spectraCanvas->Divide(4, 2);
  auto pad1 = spectraCanvas->cd(1);
  pad1->SetLogy();
  rmTruthPt->Draw();
  auto pad2 = spectraCanvas->cd(2);
  pad2->SetLogy();
  rmDetPt->Draw();
  auto pad3 = spectraCanvas->cd(3);
  pad3->SetLogy();
  rmTruthZ->Draw();
  auto pad4 = spectraCanvas->cd(4);
  pad4->SetLogy();
  rmDetZ->Draw();
  auto pad5 = spectraCanvas->cd(5);
  pad5->SetLogy();
  ruTruthPt->Draw();
  auto pad6 = spectraCanvas->cd(6);
  pad6->SetLogy();
  ruAppliedPt->Draw();
  auto pad7 = spectraCanvas->cd(7);
  pad7->SetLogy();
  ruTruthZ->Draw();
  auto pad8 = spectraCanvas->cd(8);
  pad8->SetLogy();
  ruAppliedZ->Draw();
  spectraCanvas->SaveAs("./checkRooUnfold-spectra.pdf");

  TH1D* ptMeasuredOverResponse = (TH1D*)HistRatio(rmDetPt, ruMeasuredPt, "ptMeasuredOverResponse", "RM Measured / RU Measured");
  TH1D* zMeasuredOverResponse  = (TH1D*)HistRatio(rmDetZ, ruMeasuredZ, "zMeasuredOverResponse", "RM Measured / RU Measured");
  TH1D* ptTruthOverResponse    = (TH1D*)HistRatio(rmTruthPt, ruTruthPt, "ptTruthOverResponse", "RM Truth / RU Truth");
  TH1D* zTruthOverResponse     = (TH1D*)HistRatio(rmTruthZ, ruTruthZ, "zTruthOverResponse", "RM Truth / RU Truth");
  TH1D* ptMeasuredOverApplied  = (TH1D*)HistRatio(rmDetPt, ruAppliedPt, "ptMeasuredOverApplied", "RM Measured / RU Applied");
  TH1D* zMeasuredOverApplied   = (TH1D*)HistRatio(rmDetZ, ruAppliedZ, "zMeasuredOverApplied", "RM Measured / RU Applied");

  TCanvas* ratioCanvas = new TCanvas("ratioCanvas", "ratioCanvas", 1200, 1200);
  ratioCanvas->Divide(2, 3);
  ratioCanvas->cd(1);
  ptMeasuredOverResponse->Draw();
  ratioCanvas->cd(2);
  zMeasuredOverResponse->Draw();
  ratioCanvas->cd(3);
  ptTruthOverResponse->Draw();
  ratioCanvas->cd(4);
  zTruthOverResponse->Draw();
  ratioCanvas->cd(5);
  ptMeasuredOverApplied->Draw();
  ratioCanvas->cd(6);
  zMeasuredOverApplied->Draw();
  ratioCanvas->SaveAs("./checkRooUnfold-ratios.pdf");

  TH1D* ptMeasuredMinusResponse = (TH1D*) HistDiff(rmDetPt, ruMeasuredPt, "ptMeasuredMinusResponse", "RM Measured - RU Measured; #it{p}_{T}; Measured - Response");
  TH1D* zMeasuredMinusResponse  = (TH1D*) HistDiff(rmDetZ, ruMeasuredZ, "zMeasuredMinusResponse", "RM Measured - RU Measured; #it{z}; Measured - Response");
  TH1D* ptTruthMinusResponse    = (TH1D*) HistDiff(rmTruthPt, ruTruthPt, "ptTruthMinusResponse", "RM Truth - RU Truth; #it{p}_{T}; Truth - Response");
  TH1D* zTruthMinusResponse     = (TH1D*) HistDiff(rmTruthZ, ruTruthZ, "zTruthMinusResponse", "RM Truth - RU Truth; #it{z}; Truth - Response");
  TH1D* ptMeasuredMinusApplied  = (TH1D*) HistDiff(rmDetPt, ruAppliedPt, "ptMeasuredMinusApplied", "RM Measured - RU Applied; #it{p}_{T}; Measured - Applied");
  TH1D* zMeasuredMinusApplied   = (TH1D*) HistDiff(rmDetZ, ruAppliedZ, "zMeasuredMinusApplied", "RM Measured - RU Applied; #it{z}; Measured - Applied");

  TCanvas* diffCanvas = new TCanvas("diffCanvas", "diffCanvas", 1200, 1200);
  diffCanvas->Divide(2, 3);
  diffCanvas->cd(1);
  ptMeasuredMinusResponse->Draw();
  diffCanvas->cd(2);
  zMeasuredMinusResponse->Draw();
  diffCanvas->cd(3);
  ptTruthMinusResponse->Draw();
  diffCanvas->cd(4);
  zTruthMinusResponse->Draw();
  diffCanvas->cd(5);
  ptMeasuredMinusApplied->Draw();
  diffCanvas->cd(6);
  zMeasuredMinusApplied->Draw();
  diffCanvas->SaveAs("./checkRooUnfold-diff.pdf");
}

void CreateRooUnfoldResponse(string inName = "AnalysisResults.root",
                             const double PT_TRUTH_MIN = 10., const double PT_TRUTH_MAX = 300.,
                             const double PT_DET_MIN = 10.,   const double PT_DET_MAX = 300.,
                             const double BINWIDTH_PT = 5.,
                             const double Z_TRUTH_MIN = 0.,   const double Z_TRUTH_MAX = 1.,
                             const double Z_DET_MIN = 0.,     const double Z_DET_MAX = 1.,
                             const double BINWIDTH_Z = 0.025,
                             bool makePlots = false)
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
  THnMinMax(responseMatrix);
  // !!! Make sure that all entries are >~ 1. Only for tests !!!
  responseMatrix->Scale(1e10);
  THnMinMax(responseMatrix);

  THnSparseF* tmpRM = (THnSparseF*)responseMatrix->Clone("tmpRM");
  int ptTruthMinBin = tmpRM->GetAxis(ptTruthAxis)->FindBin(PT_TRUTH_MIN + 0.1);
  int ptTruthMaxBin = tmpRM->GetAxis(ptTruthAxis)->FindBin(PT_TRUTH_MAX + 0.1);
  int zTruthMinBin  = tmpRM->GetAxis(zTruthAxis)->FindBin(Z_TRUTH_MIN + 1e-3);
  int zTruthMaxBin  = tmpRM->GetAxis(zTruthAxis)->FindBin(Z_TRUTH_MAX + 1e-3);
  int ptDetMinBin   = tmpRM->GetAxis(ptDetAxis)->FindBin(PT_DET_MIN + 0.1);
  int ptDetMaxBin   = tmpRM->GetAxis(ptDetAxis)->FindBin(PT_DET_MAX + 0.1);
  int zDetMinBin    = tmpRM->GetAxis(zDetAxis)->FindBin(Z_DET_MIN + 1e-3);
  int zDetMaxBin    = tmpRM->GetAxis(zDetAxis)->FindBin(Z_DET_MAX + 1e-3);

  // The actual bounds used for filling the distributions, derived from the input bounds
  double ptTruthMin = tmpRM->GetAxis(ptTruthAxis)->GetBinLowEdge(ptTruthMinBin);
  double ptTruthMax = tmpRM->GetAxis(ptTruthAxis)->GetBinLowEdge(ptTruthMaxBin);
  double zTruthMin = tmpRM->GetAxis(zTruthAxis)->GetBinLowEdge(zTruthMinBin);
  double zTruthMax = tmpRM->GetAxis(zTruthAxis)->GetBinLowEdge(zTruthMaxBin);
  double ptDetMin = tmpRM->GetAxis(ptDetAxis)->GetBinLowEdge(ptDetMinBin);
  double ptDetMax = tmpRM->GetAxis(ptDetAxis)->GetBinLowEdge(ptDetMaxBin);
  double zDetMin = tmpRM->GetAxis(zDetAxis)->GetBinLowEdge(zDetMinBin);
  double zDetMax = tmpRM->GetAxis(zDetAxis)->GetBinLowEdge(zDetMaxBin);

  int nbinsPtTruth = (int)((ptTruthMax - ptTruthMin)/(BINWIDTH_PT));
  int nbinsZTruth  = (int)((zTruthMax - zTruthMin)/(BINWIDTH_Z));
  int nbinsPtDet   = (int)((ptDetMax - ptDetMin)/(BINWIDTH_PT));
  int nbinsZDet    = (int)((zDetMax - zDetMin)/(BINWIDTH_Z));
  cout << nbinsPtTruth << endl;

  TH2F* rmTruth = new TH2F("rmTruth", "rmTruth; #it{p}_{T}; #it{z}",
                           nbinsPtTruth, ptTruthMin, ptTruthMax,
                           nbinsZTruth, zTruthMin, zTruthMax
                          );
  TH2F* hMiss   = new TH2F("hMiss", "hMiss; #it{p}_{T}; #it{z}",
                           nbinsPtTruth, ptTruthMin, ptTruthMax,
                           nbinsZTruth, zTruthMin, zTruthMax
                          );
  TH2F* hKinEff = new TH2F("hKinEff", "hKinEff; #it{p}_{T}; #it{z}",
                           nbinsPtTruth, ptTruthMin, ptTruthMax,
                           nbinsZTruth, zTruthMin, zTruthMax
                          );

  TH2F* rmDet = new TH2F("rmDet", "rmDet; #it{p}_{T}; #it{z}",
                         nbinsPtDet, ptDetMin, ptDetMax,
                         nbinsZDet, zDetMin, zDetMax
                         );
  TH2F* hFake = new TH2F("hFake", "hFake; #it{p}_{T}; #it{z}",
                         nbinsPtDet, ptDetMin, ptDetMax,
                         nbinsZDet, zDetMin, zDetMax
                        );


  bool checkBinning = true;
  if (checkBinning) {
    CheckBinning(rmTruth);
    CheckBinning(hMiss);
    CheckBinning(rmDet);
    CheckBinning(hFake);
    CheckBinning(hKinEff);
  }

  // Create RooUnfoldResponse and fill it
  RooUnfoldResponse *ruResponse = new RooUnfoldResponse("Response", "RM");
  ruResponse->Setup(rmDet, rmTruth);

  int* coord = new int[responseMatrix->GetNdimensions()]; //Carries the bin coordinates
  for (int iBin = 1; iBin <= responseMatrix->GetNbins(); iBin++) {
    double w          = responseMatrix->GetBinContent(iBin, coord);
    double ptTruth    = responseMatrix->GetAxis(ptTruthAxis)->GetBinCenter(coord[ptTruthAxis]);
    double zTruth     = responseMatrix->GetAxis(zTruthAxis)->GetBinCenter(coord[zTruthAxis]);
    double ptMeasured = responseMatrix->GetAxis(ptDetAxis)->GetBinCenter(coord[ptDetAxis]);
    double zMeasured  = responseMatrix->GetAxis(zDetAxis)->GetBinCenter(coord[zDetAxis]);

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
      hKinEff->Fill(ptTruth, zTruth, w);
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

  hKinEff->Divide(rmTruth);

  TFile* outFile = new TFile("RooUnfoldResponse.root", "RECREATE");
  responseMatrix->Write("responseMatrix");
  rmTruth->Write();
  rmDet->Write();
  hMiss->Write();
  hFake->Write();
  hKinEff->Write();
  ruResponse->Write();
  outFile->Write();
  // outFile->Close();

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

void ClosureTest(string testFileName = "AnalysisResults.root", string responseFileName = "RooUnfoldResponse.root",
                 double PTLOW = 40., double PTHIGH = 60.,
                 int nIter = 3)
{
  // Training RM
  TFile* responseFile = TFile::Open(TString::Format("%s", responseFileName.c_str()).Data());
  THnSparse* responseMatrix = static_cast<THnSparse*>(responseFile->Get("responseMatrix"));
  responseMatrix->SetName("responseMatrix");
  TH2F* rmTruth = (TH2F*)responseFile->Get("rmTruth");
  TH2F* rmDet = (TH2F*)responseFile->Get("rmDet");

  // Test RM
  TFile* testFile = TFile::Open(TString::Format("%s", testFileName.c_str()).Data());
  THnSparseF* testRM = (THnSparseF*)testFile->Get("jet-fragmentation/matching/jets/matchDetJetPtTrackProjPartJetPtTrackProj");
  // !!! Make sure that all entries are >~ 1. Only for tests !!!
  THnMinMax(responseMatrix);
  responseMatrix->Scale(1e10);
  THnMinMax(responseMatrix);

  // !!! These ranges have to be same as the ones used in CreateRooUnfoldResponse() !!!
  double ptTruthMin = rmTruth->GetXaxis()->GetXmin();
  double ptTruthMax = rmTruth->GetXaxis()->GetXmax();
  double ptDetMin   = rmDet->GetXaxis()->GetXmin();
  double ptDetMax   = rmDet->GetXaxis()->GetXmax();
  double zTruthMin  = rmTruth->GetYaxis()->GetXmin();
  double zTruthMax  = rmTruth->GetYaxis()->GetXmax();
  double zDetMin    = rmDet->GetYaxis()->GetXmin();
  double zDetMax    = rmDet->GetYaxis()->GetXmax();

  TH2F* testTruth  = new TH2F("testTruth", "testTruth",
                              rmTruth->GetNbinsX(),
                              rmTruth->GetXaxis()->GetBinLowEdge(1),
                              rmTruth->GetXaxis()->GetBinLowEdge(rmTruth->GetNbinsX() + 1),
                              rmTruth->GetNbinsY(),
                              rmTruth->GetYaxis()->GetBinLowEdge(1),
                              rmTruth->GetYaxis()->GetBinLowEdge(rmTruth->GetNbinsY() + 1)
                             );
  TH2F* testMiss   = new TH2F("testMiss", "testMiss",
                              rmTruth->GetNbinsX(),
                              rmTruth->GetXaxis()->GetBinLowEdge(1),
                              rmTruth->GetXaxis()->GetBinLowEdge(rmTruth->GetNbinsX() + 1),
                              rmTruth->GetNbinsY(),
                              rmTruth->GetYaxis()->GetBinLowEdge(1),
                              rmTruth->GetYaxis()->GetBinLowEdge(rmTruth->GetNbinsY() + 1)
                             );
  TH2F* testKinEff = new TH2F("testKinEff", "testKinEff",
                              rmTruth->GetNbinsX(),
                              rmTruth->GetXaxis()->GetBinLowEdge(1),
                              rmTruth->GetXaxis()->GetBinLowEdge(rmTruth->GetNbinsX() + 1),
                              rmTruth->GetNbinsY(),
                              rmTruth->GetYaxis()->GetBinLowEdge(1),
                              rmTruth->GetYaxis()->GetBinLowEdge(rmTruth->GetNbinsY() + 1)
                             );

  TH2F* testDet  = new TH2F("testDet", "testDet",
                            rmDet->GetNbinsX(),
                            rmDet->GetXaxis()->GetBinLowEdge(1),
                            rmDet->GetXaxis()->GetBinLowEdge(rmDet->GetNbinsX() + 1),
                            rmDet->GetNbinsY(),
                            rmDet->GetYaxis()->GetBinLowEdge(1),
                            rmDet->GetYaxis()->GetBinLowEdge(rmDet->GetNbinsY() + 1)
                           );
  TH2F* testFake = new TH2F("testFake", "testFake",
                            rmDet->GetNbinsX(),
                            rmDet->GetXaxis()->GetBinLowEdge(1),
                            rmDet->GetXaxis()->GetBinLowEdge(rmDet->GetNbinsX() + 1),
                            rmDet->GetNbinsY(),
                            rmDet->GetYaxis()->GetBinLowEdge(1),
                            rmDet->GetYaxis()->GetBinLowEdge(rmDet->GetNbinsY() + 1)
                           );

  int* coord = new int[testRM->GetNdimensions()]; //Carries the bin coordinates
  for (int iBin = 1; iBin <= testRM->GetNbins(); iBin++) {
    double w          = testRM->GetBinContent(iBin, coord);
    double ptTruth    = testRM->GetAxis(ptTruthAxis)->GetBinCenter(coord[ptTruthAxis]);
    double zTruth     = testRM->GetAxis(zTruthAxis)->GetBinCenter(coord[zTruthAxis]);
    double ptMeasured = testRM->GetAxis(ptDetAxis)->GetBinCenter(coord[ptDetAxis]);
    double zMeasured  = testRM->GetAxis(zDetAxis)->GetBinCenter(coord[zDetAxis]);

    bool inAcceptanceTruth = false;
    bool inAcceptanceDetector = false;

    inAcceptanceTruth = (ptTruth >= ptTruthMin) * (ptTruth < ptTruthMax);
    inAcceptanceDetector = (ptMeasured >= ptDetMin) * (ptMeasured < ptDetMax);

    if (inAcceptanceTruth) {
      testTruth->Fill(ptTruth, zTruth, w);
    }
    if (inAcceptanceDetector) {
      testDet->Fill(ptMeasured, zMeasured, w);
    }

    if (inAcceptanceTruth && inAcceptanceDetector) {
      testKinEff->Fill(ptTruth, zTruth, w);
    }
    else if (inAcceptanceTruth && !inAcceptanceDetector) {
      testMiss->Fill(ptTruth, zTruth, w);
      // cout << "Miss: ptT=" << ptTruth << ", zT=" << zTruth << ", ptM=" << ptMeasured << ", zM=" << zMeasured << endl;
    }
    else if (!inAcceptanceTruth && inAcceptanceDetector) {
      testFake->Fill(ptMeasured, zMeasured, w);
      // cout << "Fake: ptT=" << ptTruth << ", zT=" << zTruth << ", ptM=" << ptMeasured << ", zM=" << zMeasured << endl;
    }
  }
  delete [] coord;
  testKinEff->Divide(testDet);

  // This is the section where the actual unfolding starts and RooUnfold is involved
  string ruResponseName = "Response";
  string unfoldName = "Bayesian_Unfolding";
  string unfoldTitle = unfoldName;
  int doSmoothing = 0;
  // Load the response object from the Response file, create the Bayesian unfolding object,
  // Extract the unfolded distribution as a histogram, apply the Response matrix to the unfolded result
  RooUnfoldResponse* ruResponse = static_cast<RooUnfoldResponse*>(responseFile->Get("Response"));
  RooUnfoldBayes Unfolding_bayes(ruResponse, testDet, nIter, doSmoothing, unfoldName.c_str(), unfoldTitle.c_str());
  TH2F* unfoldedTest = (TH2F*)Unfolding_bayes.Hreco(RooUnfold::kCovToy);
  unfoldedTest->SetName("unfoldedTest");
  auto Unf_Bayes = &Unfolding_bayes;
  TH2F* refoldedTest = (TH2F*)ruResponse->ApplyToTruth(unfoldedTest, "refoldedTest");
  // Include the fakes in order to compare with the input distribution
  // TH2F* hFake = (TH2F*)responseFile->Get("hFake");
  // refoldedTest->Add(hFake);

  // The unfolding is complete, the following sections refer to projections in one dimension and calculating the ratio between distributions
  int ptLowTruth = rmTruth->GetXaxis()->FindBin(PTLOW);
  int ptHighTruth = rmTruth->GetXaxis()->FindBin(PTHIGH);
  int ptLowDet = rmDet->GetXaxis()->FindBin(PTLOW);
  int ptHighDet = rmDet->GetXaxis()->FindBin(PTHIGH);

  //--------------Training sample (Sample creating the response)--------------------//
    //Truth level
  TH1D* trainingPtTruth = (TH1D*)rmTruth->ProjectionX("trainingPtTruth");
  trainingPtTruth->SetTitle("Training truth pt");
  // trainingPtTruth->Scale(1./trainingPtTruth->Integral(), "width");

  TH1D* trainingZTruth = (TH1D*)rmTruth->ProjectionY("trainingZTruth", ptLowTruth, ptHighTruth - 1);
  trainingZTruth->SetTitle("Training truth z");
  // trainingZTruth->Scale(1./trainingZTruth->Integral(), "width");

    //Detector level
  TH1D* trainingPtDetector = (TH1D*)rmDet->ProjectionX("trainingPtDetector");
  trainingPtDetector->SetTitle("Training detector pt");
  // trainingPtDetector->Scale(1./trainingPtDetector->Integral(), "width");

  TH1D* trainingZDetector = (TH1D*)rmDet->ProjectionY("trainingZDetector", ptLowDet, ptHighDet - 1);
  trainingZDetector->SetTitle("Training detector z");
  // trainingZDetector->Scale(1./trainingZDetector->Integral(), "width");

  //--------------Test sample (MC sample treated as pseudodata)--------------------//
    //Truth level
  TH1D* testPtTruth = (TH1D*)testTruth->ProjectionX("testPtTruth");
  testPtTruth->SetTitle("Test truth pt");
  // testPtTruth->Scale(1./testPtTruth->Integral(), "width");

  TH1D* testZTruth = (TH1D*)testTruth->ProjectionY("testZTruth", ptLowTruth, ptHighTruth - 1);
  testZTruth->SetTitle("Test truth z");
  // testZTruth->Scale(1./testZTruth->Integral(), "width");

    //Detector level
  TH1D* testPtDetector = (TH1D*)testDet->ProjectionX("testPtDetector");
  testPtDetector->SetTitle("Test detector pt");
  // testPtDetector->Scale(1./testPtDetector->Integral(), "width");

  TH1D* testZDetector = (TH1D*)testDet->ProjectionY("testZDetector", ptLowDet, ptHighDet - 1);
  testZDetector->SetTitle("Test truth z");
  // testZDetector->Scale(1./testZDetector->Integral(), "width");

    //Unfolded
  TH1D* testPtUnfolded = unfoldedTest->ProjectionX("testPtUnfolded");
  testPtUnfolded->SetTitle("Unfolded pt");
  // testPtUnfolded->Scale(1./testPtUnfolded->Integral(), "width");

  TH1D* testZUnfolded = unfoldedTest->ProjectionY("testZUnfolded", ptLowTruth, ptHighTruth - 1);
  testZUnfolded->SetTitle("Unfolded z");
  // testZUnfolded->Scale(1./testZUnfolded->Integral(), "width");

    //Refolded
  TH1D* testPtRefolded = refoldedTest->ProjectionX("testPtRefolded");
  testPtRefolded->SetTitle("Refolded pt");
  // testPtRefolded->Scale(1./testPtRefolded->Integral(), "width");

  TH1D* testZRefolded = refoldedTest->ProjectionY("testZRefolded", ptLowDet, ptHighDet - 1);
  testZRefolded->SetTitle("Refolded z");
  // testZRefolded->Scale(1./testZRefolded->Integral(), "width");

  /*
  // Ratios of the projections
  // Training distributions - Detector response corrections
  TH1D* trainingPtTruthOverDetector = (TH1D*)trainingPtTruth->Clone("trainingPtTruthOverDetector");
  trainingPtTruthOverDetector->Divide(trainingPtDetector);
  trainingPtTruthOverDetector->SetTitle("Training pt truth/detector; #it{p}_{T}; #frac{Truth}{Detector}");

  TH1D* trainingZTruthOverDetector = (TH1D*)trainingZTruth->Clone("trainingZTruthOverDetector");
  trainingZTruthOverDetector->Divide(trainingZDetector);
  trainingZTruthOverDetector->SetTitle("Training z truth/detector; #it{z}; #frac{Truth}{Detector}");

  //Test distributions - Detector response corrections
  TH1D* testPtUnfoldedOverDetector = (TH1D*)testPtUnfolded->Clone("testPtUnfoldedOverDetector");
  // testPtUnfoldedOverDetector->Rebin(2);
  testPtUnfoldedOverDetector->Divide(testPtDetector);
  testPtUnfoldedOverDetector->SetTitle("Test pt unfolded/detector; #it{p}_{T}; #frac{Unfolded}{Detector}");

  TH1D* testZUnfoldedOverDetector = (TH1D*)testZUnfolded->Clone("testZUnfoldedOverDetector");
  // testZUnfoldedOverDetector->Rebin(2);
  testZUnfoldedOverDetector->Divide(testZDetector);
  testZUnfoldedOverDetector->SetTitle("Test z unfolded/detector; #it{z}; #frac{Unfolded}{Detector}");
  */

  TCanvas* spectraCanvas = new TCanvas("spectraCanvas", "spectraCanvas", 2000, 1000);
  spectraCanvas->Divide(4, 2);
  auto pad1 = spectraCanvas->cd(1);
  pad1->SetLogy();
  testPtTruth->Draw();
  auto pad2 = spectraCanvas->cd(2);
  pad2->SetLogy();
  testPtDetector->Draw();
  auto pad3 = spectraCanvas->cd(3);
  pad3->SetLogy();
  testZTruth->Draw();
  auto pad4 = spectraCanvas->cd(4);
  pad4->SetLogy();
  testZDetector->Draw();
  auto pad5 = spectraCanvas->cd(5);
  pad5->SetLogy();
  testPtUnfolded->Draw();
  auto pad6 = spectraCanvas->cd(6);
  pad6->SetLogy();
  testPtRefolded->Draw();
  auto pad7 = spectraCanvas->cd(7);
  pad7->SetLogy();
  testZUnfolded->Draw();
  auto pad8 = spectraCanvas->cd(8);
  pad8->SetLogy();
  testZRefolded->Draw();
  spectraCanvas->SaveAs("./closureTest-spectra.pdf");

  TH1D* testPtUnfoldedOverTruth = (TH1D*)HistRatio(testPtUnfolded, testPtTruth, "testPtUnfoldedOverTruth", "Test pt unfolded/truth; #it{p}_{T}; #frac{Unfolded}{Truth}");
  TH1D* testZUnfoldedOverTruth = (TH1D*)HistRatio(testZUnfolded, testZTruth, "testZUnfoldedOverTruth", "Test z unfolded/truth; #it{z}; #frac{Unfolded}{Truth}");
  TH1D* testPtRefoldedOverDetector = (TH1D*)HistRatio(testPtRefolded, testPtDetector, "testPtRefoldedOverDetector", "Test pt refolded/detector; #it{p}_{T}; #frac{Refolded}{Detector}");
  TH1D* testZRefoldedOverDetector = (TH1D*)HistRatio(testZRefolded, testZDetector, "testZRefoldedOverDetector", "Test z refolded/detector; #it{z}; #frac{Refolded}{Detector}");

  TCanvas* ratioCanvas = new TCanvas("ratioCanvas", "ratioCanvas", 1000, 1000);
  ratioCanvas->Divide(2, 2);
  ratioCanvas->cd(1);
  testPtUnfoldedOverTruth->Draw();
  ratioCanvas->cd(2);
  testZUnfoldedOverTruth->Draw();
  ratioCanvas->cd(3);
  testPtRefoldedOverDetector->Draw();
  ratioCanvas->cd(4);
  testZRefoldedOverDetector->Draw();
  ratioCanvas->SaveAs("./closureTest-ratio.pdf");

  TH1D* testPtUnfoldedMinusTruth = (TH1D*)HistDiff(testPtUnfolded, testPtTruth, "testPtUnfoldedMinusTruth", "Test pt unfolded - truth; #it{p}_{T}; Unfolded - Truth");
  TH1D* testZUnfoldedMinusTruth = (TH1D*)HistDiff(testZUnfolded, testZTruth, "testZUnfoldedMinusTruth", "Test z unfolded - truth; #it{z}; Unfolded - Truth");
  TH1D* testPtRefoldedMinusDetector = (TH1D*)HistDiff(testPtRefolded, testPtDetector, "testPtRefoldedMinusDetector", "Test pt refolded - detector; #it{p}_{T}; Refolded - Detector");
  TH1D* testZRefoldedMinusDetector = (TH1D*)HistDiff(testZRefolded, testZDetector, "testZRefoldedMinusDetector", "Test z refolded - detector; #it{z}; Refolded - Detector");

  TCanvas* diffCanvas = new TCanvas("diffCanvas", "diffCanvas", 1000, 1000);
  diffCanvas->Divide(2, 2);
  diffCanvas->cd(1);
  testPtUnfoldedMinusTruth->Draw();
  diffCanvas->cd(2);
  testZUnfoldedMinusTruth->Draw();
  diffCanvas->cd(3);
  testPtRefoldedMinusDetector->Draw();
  diffCanvas->cd(4);
  testZRefoldedMinusDetector->Draw();
  diffCanvas->SaveAs("./closureTest-diff.pdf");

  //Saving the output in a new file
  string outFileName = "closureTest_projection.root";
  TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");

  // Response objects
  responseMatrix->Write(responseMatrix->GetName(), TObject::kOverwrite);
  rmTruth->Write(rmTruth->GetName(), TObject::kOverwrite);
  rmDet->Write(rmDet->GetName(), TObject::kOverwrite);
  testRM->Write(testRM->GetName(), TObject::kOverwrite);
  testTruth->Write(testTruth->GetName(), TObject::kOverwrite);
  testDet->Write(testDet->GetName(), TObject::kOverwrite);
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

  // trainingPtTruthOverDetector->Write(trainingPtTruthOverDetector->GetName(), TObject::kOverwrite);
  // trainingZTruthOverDetector->Write(trainingZTruthOverDetector->GetName(), TObject::kOverwrite);
  // testPtUnfoldedOverDetector->Write(testPtUnfoldedOverDetector->GetName(), TObject::kOverwrite);
  // testZUnfoldedOverDetector->Write(testZUnfoldedOverDetector->GetName(), TObject::kOverwrite);
  testPtUnfoldedOverTruth->Write(testPtUnfoldedOverTruth->GetName(), TObject::kOverwrite);
  testZUnfoldedOverTruth->Write(testZUnfoldedOverTruth->GetName(), TObject::kOverwrite);
  testPtRefoldedOverDetector->Write(testPtRefoldedOverDetector->GetName(), TObject::kOverwrite);
  testZRefoldedOverDetector->Write(testZRefoldedOverDetector->GetName(), TObject::kOverwrite);

  testPtUnfoldedMinusTruth->Write(testPtUnfoldedMinusTruth->GetName(), TObject::kOverwrite);
  testZUnfoldedMinusTruth->Write(testZUnfoldedMinusTruth->GetName(), TObject::kOverwrite);
  testPtRefoldedMinusDetector->Write(testPtRefoldedMinusDetector->GetName(), TObject::kOverwrite);
  testZRefoldedMinusDetector->Write(testZRefoldedMinusDetector->GetName(), TObject::kOverwrite);

  // testFile->Close();
  // responseFile->Close();
  // outFile->Close();
}

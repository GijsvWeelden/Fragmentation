
// This task:
// * Loads the response from the task output
// * Creates the RooUnfoldResponse object used in unfolding

using std::cout, std::endl, std::string;

int ptDetAxis   = 0;
int zDetAxis    = 1;
int ptTruthAxis = 2;
int zTruthAxis  = 3;

void CheckBinning(const TH1* hist)
{
  cout << hist->GetName() << " "
    << hist->GetXaxis()->GetBinLowEdge(1) << " - " << hist->GetXaxis()->GetBinLowEdge(hist->GetNbinsX() + 1)
    << " (" << hist->GetNbinsX() << " * " << hist->GetXaxis()->GetBinWidth(1) << ")"
    << endl;
}
void CheckBinning(const TH2* hist)
{
  cout << hist->GetName() << " "
    << hist->GetXaxis()->GetBinLowEdge(1) << " - " << hist->GetXaxis()->GetBinLowEdge(hist->GetNbinsX() + 1)
    << " (" << hist->GetNbinsX() << " * " << hist->GetXaxis()->GetBinWidth(1) << "), "
    << hist->GetYaxis()->GetBinLowEdge(1) << " - " << hist->GetYaxis()->GetBinLowEdge(hist->GetNbinsY() + 1)
    << " (" << hist->GetNbinsY() << " * " << hist->GetYaxis()->GetBinWidth(1) << ")"
    << endl;
}
void CompareHists(const TH1* histA, const TH1* histB)
{
  cout << histA->GetName() << " - " << histB->GetName() << endl;
  for (int i = 0; i < histA->GetNbinsX(); i++)
  {
    double binDiff = histA->GetBinContent(i) - histB->GetBinContent(i);
    if ( abs(binDiff) > 0.) cout << "Bin " << i << " = " << binDiff << " ";
  }
  cout << endl;
}
void CompareHists(const TH2* histA, const TH2* histB)
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
TH1* HistDiff(const TH1* histA, const TH1* histB, string outName = "outHist", string outTitle = "")
{
  TH1* outHist = (TH1*)histA->Clone(outName.c_str());
  outHist->SetTitle(outTitle.c_str());
  outHist->Add(histB, -1.);
  return outHist;
}
TH1* HistRatio(const TH1* histA, const TH1* histB, string outName = "outHist", string outTitle = "")
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
void THnMinMax(const THnSparse* hist)
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
  TH1F* rmTruthZ = (TH1F*)rmTruth->ProjectionX("rmTruthZ");
  TH1F* rmTruthPt = (TH1F*)rmTruth->ProjectionY("rmTruthPt");

  TH2D* rmDet = (TH2D*)inFile->Get("rmDet");
  TH1F* rmDetZ = (TH1F*)rmDet->ProjectionX("rmDetZ");
  TH1F* rmDetPt = (TH1F*)rmDet->ProjectionY("rmDetPt");

  RooUnfoldResponse* response = (RooUnfoldResponse*)inFile->Get("Response");
  TH2D* ruMeasured = (TH2D*)response->Hmeasured();
  TH1D* ruMeasuredZ = (TH1D*)ruMeasured->ProjectionX("ruMeasuredZ");
  TH1D* ruMeasuredPt = (TH1D*)ruMeasured->ProjectionY("ruMeasuredPt");

  TH2D* ruTruth = (TH2D*)response->Htruth();
  TH1D* ruTruthZ = (TH1D*)ruTruth->ProjectionX("ruTruthZ");
  TH1D* ruTruthPt = (TH1D*)ruTruth->ProjectionY("ruTruthPt");

  TH2D* ruApplied = (TH2D*)response->ApplyToTruth();
  TH2D* hFake = (TH2D*)inFile->Get("hFake");
  ruApplied->Add(hFake);
  TH1D* ruAppliedZ = (TH1D*)ruApplied->ProjectionX("ruAppliedZ");
  TH1D* ruAppliedPt = (TH1D*)ruApplied->ProjectionY("ruAppliedPt");

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

  TH1D* ptMeasuredMinusResponse = (TH1D*) HistDiff(rmDetPt, ruMeasuredPt, "ptMeasuredMinusResponse", "RM Measured - RU Measured; #it{p}_{T, ch. jet}; Measured - Response");
  TH1D* zMeasuredMinusResponse  = (TH1D*) HistDiff(rmDetZ, ruMeasuredZ, "zMeasuredMinusResponse", "RM Measured - RU Measured; #it{z}; Measured - Response");
  TH1D* ptTruthMinusResponse    = (TH1D*) HistDiff(rmTruthPt, ruTruthPt, "ptTruthMinusResponse", "RM Truth - RU Truth; #it{p}_{T, ch. jet}; Truth - Response");
  TH1D* zTruthMinusResponse     = (TH1D*) HistDiff(rmTruthZ, ruTruthZ, "zTruthMinusResponse", "RM Truth - RU Truth; #it{z}; Truth - Response");
  TH1D* ptMeasuredMinusApplied  = (TH1D*) HistDiff(rmDetPt, ruAppliedPt, "ptMeasuredMinusApplied", "RM Measured - RU Applied; #it{p}_{T, ch. jet}; Measured - Applied");
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
  string saveName = "RooUnfoldResponse.root";
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

  TH2F* rmTruth = new TH2F("rmTruth", "rmTruth; #it{z}; #it{p}_{T, ch. jet}",
                           nbinsZTruth, zTruthMin, zTruthMax,
                           nbinsPtTruth, ptTruthMin, ptTruthMax
                          );
  TH2F* hMiss   = new TH2F("hMiss", "hMiss; #it{z}; #it{p}_{T, ch. jet}",
                           nbinsZTruth, zTruthMin, zTruthMax,
                           nbinsPtTruth, ptTruthMin, ptTruthMax
                          );
  TH2F* hKinEff = new TH2F("hKinEff", "hKinEff; #it{z}; #it{p}_{T, ch. jet}",
                           nbinsZTruth, zTruthMin, zTruthMax,
                           nbinsPtTruth, ptTruthMin, ptTruthMax
                          );

  TH2F* rmDet   = new TH2F("rmDet", "rmDet; #it{z}; #it{p}_{T, ch. jet}",
                           nbinsZDet, zDetMin, zDetMax,
                           nbinsPtDet, ptDetMin, ptDetMax
                           );
  TH2F* hFake   = new TH2F("hFake", "hFake; #it{z}; #it{p}_{T, ch. jet}",
                           nbinsZDet, zDetMin, zDetMax,
                           nbinsPtDet, ptDetMin, ptDetMax
                          );


  bool checkBinning = false;
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
      rmTruth->Fill(zTruth, ptTruth, w);
    }
    if (inAcceptanceDetector) {
      rmDet->Fill(zMeasured, ptMeasured, w);
    }

    if (inAcceptanceTruth && inAcceptanceDetector) {
      ruResponse->Fill(zMeasured, ptMeasured, zTruth, ptTruth, w);
      hKinEff->Fill(zTruth, ptTruth, w);
    }
    else if (inAcceptanceTruth && !inAcceptanceDetector) {
      ruResponse->Miss(zTruth, ptTruth, w);
      hMiss->Fill(zTruth, ptTruth, w);
      // cout << "Miss: ptT=" << ptTruth << ", zT=" << zTruth << ", ptM=" << ptMeasured << ", zM=" << zMeasured << endl;
    }
    else if (!inAcceptanceTruth && inAcceptanceDetector) {
      ruResponse->Fake(zMeasured, ptMeasured, w);
      hFake->Fill(zMeasured, ptMeasured, w);
      // cout << "Fake: ptT=" << ptTruth << ", zT=" << zTruth << ", ptM=" << ptMeasured << ", zM=" << zMeasured << endl;
    }
  }
  delete [] coord;

  hKinEff->Divide(rmTruth);

  TFile* outFile = new TFile(TString::Format("./%s", saveName.c_str()).Data(), "RECREATE");
  responseMatrix->Write("responseMatrix");
  rmTruth->Write();
  rmDet->Write();
  hMiss->Write();
  hFake->Write();
  hKinEff->Write();
  ruResponse->Write();
  outFile->Write();

  if (!makePlots) return;
  /*
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
  */
} // CreateRooUnfoldResponse

void ClosureTest(string testFileName = "AnalysisResults.root", string responseFileName = "RooUnfoldResponse.root",
                 double PTLOW = 40., double PTHIGH = 60.,
                 int itmin = 2, int itmax = 3,
                 bool makePlots = false, string drawOption = "text90 hist")
{
  double time = clock();
  // Training RM
  TFile* responseFile = TFile::Open(TString::Format("%s", responseFileName.c_str()).Data(), "READ");
  THnSparse* responseMatrix = static_cast<THnSparse*>(responseFile->Get("responseMatrix"));
  responseMatrix->SetName("responseMatrix");
  TH2F* rmTruth = (TH2F*)responseFile->Get("rmTruth");
  TH2F* rmDet = (TH2F*)responseFile->Get("rmDet");

  // Test RM
  TFile* testFile = TFile::Open(TString::Format("%s", testFileName.c_str()).Data(), "READ");
  THnSparseF* testRM = (THnSparseF*)testFile->Get("jet-fragmentation/matching/jets/matchDetJetPtTrackProjPartJetPtTrackProj");

  // !!! These ranges have to be same as the ones used in CreateRooUnfoldResponse() !!!
  double zTruthMin  = rmTruth->GetXaxis()->GetXmin();
  double zTruthMax  = rmTruth->GetXaxis()->GetXmax();
  double zDetMin    = rmDet->GetXaxis()->GetXmin();
  double zDetMax    = rmDet->GetXaxis()->GetXmax();
  double ptTruthMin = rmTruth->GetYaxis()->GetXmin();
  double ptTruthMax = rmTruth->GetYaxis()->GetXmax();
  double ptDetMin   = rmDet->GetYaxis()->GetXmin();
  double ptDetMax   = rmDet->GetYaxis()->GetXmax();

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

  TH2F* testDet    = new TH2F("testDet", "testDet",
                              rmDet->GetNbinsX(),
                              rmDet->GetXaxis()->GetBinLowEdge(1),
                              rmDet->GetXaxis()->GetBinLowEdge(rmDet->GetNbinsX() + 1),
                              rmDet->GetNbinsY(),
                              rmDet->GetYaxis()->GetBinLowEdge(1),
                              rmDet->GetYaxis()->GetBinLowEdge(rmDet->GetNbinsY() + 1)
                             );
  TH2F* testFake   = new TH2F("testFake", "testFake",
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
      testTruth->Fill(zTruth, ptTruth, w);
    }
    if (inAcceptanceDetector) {
      testDet->Fill(zMeasured, ptMeasured, w);
    }

    if (inAcceptanceTruth && inAcceptanceDetector) {
      testKinEff->Fill(zTruth, ptTruth, w);
    }
    else if (inAcceptanceTruth && !inAcceptanceDetector) {
      testMiss->Fill(zTruth, ptTruth, w);
      // cout << "Miss: ptT=" << ptTruth << ", zT=" << zTruth << ", ptM=" << ptMeasured << ", zM=" << zMeasured << endl;
    }
    else if (!inAcceptanceTruth && inAcceptanceDetector) {
      testFake->Fill(zMeasured, ptMeasured, w);
      // cout << "Fake: ptT=" << ptTruth << ", zT=" << zTruth << ", ptM=" << ptMeasured << ", zM=" << zMeasured << endl;
    }
  }
  delete [] coord;
  testKinEff->Divide(testDet);

  // Output file
  string outFileName = "./closureTest.root";
  TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");
  responseMatrix->Write("trainingRM", TObject::kOverwrite);
  rmTruth->Write(rmTruth->GetName(), TObject::kOverwrite);
  rmDet->Write(rmDet->GetName(), TObject::kOverwrite);
  testRM->Write("testRM", TObject::kOverwrite);
  testTruth->Write(testTruth->GetName(), TObject::kOverwrite);
  testDet->Write(testDet->GetName(), TObject::kOverwrite);

  // Make projections of truth, detector distributions to compare to unfolding results
  int ptTruthLowBin  = rmTruth->GetYaxis()->FindBin(PTLOW);
  int ptTruthHighBin = rmTruth->GetYaxis()->FindBin(PTHIGH);
  int ptDetLowBin    = rmDet->GetYaxis()->FindBin(PTLOW);
  int ptDetHighBin   = rmDet->GetYaxis()->FindBin(PTHIGH);

  double ptTruthLow  = rmTruth->GetYaxis()->GetBinLowEdge(ptTruthLowBin);
  double ptTruthHigh = rmTruth->GetYaxis()->GetBinLowEdge(ptTruthHighBin);
  double ptDetLow    = rmDet->GetYaxis()->GetBinLowEdge(ptDetLowBin);
  double ptDetHigh   = rmDet->GetYaxis()->GetBinLowEdge(ptDetHighBin);

  // Training sample (Sample creating the response)
  TH1D* trainingZTruth = (TH1D*)rmTruth->ProjectionX("trainingZTruth", ptTruthLowBin, ptTruthHighBin - 1);
  trainingZTruth->SetTitle("Training truth z; #it{z}; 1/N dN/d#it{z}");
  TH1D* trainingPtTruth = (TH1D*)rmTruth->ProjectionY("trainingPtTruth");
  trainingPtTruth->SetTitle("Training truth pt; #it{p}_{T}; 1/N dN/d#it{p}_{T}");
  TH1D* trainingZDetector = (TH1D*)rmDet->ProjectionX("trainingZDetector", ptDetLowBin, ptDetHighBin - 1);
  trainingZDetector->SetTitle("Training detector z; #it{z}; 1/N dN/d#it{z}");
  TH1D* trainingPtDetector = (TH1D*)rmDet->ProjectionY("trainingPtDetector");
  trainingPtDetector->SetTitle("Training detector pt; #it{p}_{T}; 1/N dN/d#it{p}_{T}");

  trainingPtTruth->Write(trainingPtTruth->GetName(), TObject::kOverwrite);
  trainingZTruth->Write(trainingZTruth->GetName(), TObject::kOverwrite);
  trainingPtDetector->Write(trainingPtDetector->GetName(), TObject::kOverwrite);
  trainingZDetector->Write(trainingZDetector->GetName(), TObject::kOverwrite);

  // Test sample (MC sample treated as pseudodata)
  TH1D* testZTruth = (TH1D*)testTruth->ProjectionX("testZTruth", ptTruthLowBin, ptTruthHighBin - 1);
  testZTruth->SetTitle("Test truth z");
  testZTruth->Write(testZTruth->GetName(), TObject::kOverwrite);
  TH1D* testPtTruth = (TH1D*)testTruth->ProjectionY("testPtTruth");
  testPtTruth->SetTitle("Test truth pt");
  testPtTruth->Write(testPtTruth->GetName(), TObject::kOverwrite);

  TH1D* testZDetector = (TH1D*)testDet->ProjectionX("testZDetector", ptDetLowBin, ptDetHighBin - 1);
  testZDetector->SetTitle("Test truth z");
  testZDetector->Write(testZDetector->GetName(), TObject::kOverwrite);
  TH1D* testPtDetector = (TH1D*)testDet->ProjectionY("testPtDetector");
  testPtDetector->SetTitle("Test detector pt");
  testPtDetector->Write(testPtDetector->GetName(), TObject::kOverwrite);

  // This is the section where the actual unfolding starts and RooUnfold is involved
  string ruResponseName = "Response";
  string unfoldName = "Bayesian_Unfolding";
  string unfoldTitle = unfoldName;
  int doSmoothing = 0;
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovToy;

  int nIter = itmax + 1 - itmin;
  RooUnfoldBayes ruBayes[nIter];
  TH2F* unfolded[nIter];
  TH2F* refolded[nIter];
  TH2D* covmat[nIter];
  TH2D* pearson[nIter];

  TH1D* unfoldedPt[nIter];
  TH1D* unfoldedZ[nIter];
  TH1D* refoldedPt[nIter];
  TH1D* refoldedZ[nIter];
  TH1D* ptUnfoldedOverTruth[nIter];
  TH1D* zUnfoldedOverTruth[nIter];
  TH1D* ptRefoldedOverDetector[nIter];
  TH1D* zRefoldedOverDetector[nIter];
  TH1D* ptUnfoldedMinusTruth[nIter];
  TH1D* zUnfoldedMinusTruth[nIter];
  TH1D* ptRefoldedMinusDetector[nIter];
  TH1D* zRefoldedMinusDetector[nIter];
  TH1D* ptTruthDiffRatio[nIter];
  TH1D* zTruthDiffRatio[nIter];
  TH1D* ptDetectorDiffRatio[nIter];
  TH1D* zDetectorDiffRatio[nIter];

  // Load the response object from the Response file
  RooUnfoldResponse* ruResponse = static_cast<RooUnfoldResponse*>(responseFile->Get("Response"));
  ruResponse->Write(ruResponse->GetName(), TObject::kOverwrite);
  for (int iter = 0; iter < nIter; iter++)
  {
    // Create the Bayesian unfolding object,
    // Extract the unfolded distribution as a histogram
    // Apply the Response matrix to the unfolded result
    // Retrieve the covariance matrix and Pearson coefficients

    // int iteration = iter + 1; // Used for names and titles
    int iteration = iter + itmin; // Used for names and titles
    cout << "Unfolding iteration " << iteration << endl;
    ruBayes[iter] = RooUnfoldBayes(ruResponse, testDet, iteration, doSmoothing,
                                   TString::Format("%s_iter%d", unfoldName.c_str(), iteration).Data(),
                                   TString::Format("%s_iter%d", unfoldTitle.c_str(), iteration).Data()
                                  );
    unfolded[iter] = (TH2F*)ruBayes[iter].Hreco(errorTreatment);
    unfolded[iter]->SetName(TString::Format("unfolded_iter%d", iteration).Data());
    refolded[iter] = (TH2F*)ruResponse->ApplyToTruth(unfolded[iter], TString::Format("refolded_iter%d", iteration).Data());
    // Include the fakes in order to compare with the input distribution
    refolded[iter]->Add(testFake);
    TH2D htmp(ruBayes[iter].Ereco(errorTreatment));
    covmat[iter] = (TH2D*)htmp.Clone(TString::Format("covmat_iter%d", iteration).Data());

    // Calculate the Pearson coefficients
    // One block per pt bin
    pearson[iter] = (TH2D*)covmat[iter]->Clone("pearson");
    pearson[iter]->Reset();
    int nCellsPerBlock = rmTruth->GetNbinsX();
    int nBlocks = rmTruth->GetNbinsY();
    for (int iPtx = 1; iPtx <= nBlocks; iPtx++) {
      for (int iPty = 1; iPty <= nBlocks; iPty++) {
        for (int iZx = 1; iZx <= nCellsPerBlock; iZx++) {
          for (int iZy = 1; iZy <= nCellsPerBlock; iZy++) {
            int xBin = (iPtx - 1) * nCellsPerBlock + iZx;
            int yBin = (iPty - 1) * nCellsPerBlock + iZy;
            double cxy = covmat[iter]->GetBinContent(xBin, yBin);
            double cxx = covmat[iter]->GetBinContent(xBin, xBin);
            double cyy = covmat[iter]->GetBinContent(yBin, yBin);
            double p = cxy/TMath::Sqrt(cxx * cyy);
            // double p = covmat[iter]->GetBinContent(iZx, iZy)/TMath::Sqrt(covmat[iter]->GetBinContent(iZx, iZx) * covmat[iter]->GetBinContent(iZy, iZy));
            pearson[iter]->SetBinContent(xBin, yBin, p);
          }
        }
      }
    }

    ruBayes[iter].Write(ruBayes[iter].GetName(), TObject::kOverwrite);
    unfolded[iter]->Write(unfolded[iter]->GetName(), TObject::kOverwrite);
    refolded[iter]->Write(refolded[iter]->GetName(), TObject::kOverwrite);
    covmat[iter]->Write(covmat[iter]->GetName(), TObject::kOverwrite);
    pearson[iter]->Write(pearson[iter]->GetName(), TObject::kOverwrite);

    // ------------------ Compare 1D spectra for each iteration of the unfolding ------------------
    unfoldedZ[iter] = (TH1D*)unfolded[iter]->ProjectionX(TString::Format("unfoldedZ_iter%d", iteration).Data(), ptTruthLowBin, ptTruthHighBin - 1);
    unfoldedZ[iter]->SetTitle(TString::Format("Unfolded z iteration %d, #it{p}_{T, ch. jet}: %.0f - %.0f", iteration, ptTruthLow, ptTruthHigh).Data());
    unfoldedZ[iter]->Write(unfoldedZ[iter]->GetName(), TObject::kOverwrite);
    unfoldedPt[iter] = (TH1D*)unfolded[iter]->ProjectionY(TString::Format("unfoldedPt_iter%d", iteration).Data());
    unfoldedPt[iter]->SetTitle(TString::Format("Unfolded pt iteration %d; #it{p}_{T, ch. jet}", iteration).Data());
    unfoldedPt[iter]->Write(unfoldedPt[iter]->GetName(), TObject::kOverwrite);

    refoldedZ[iter] = (TH1D*)refolded[iter]->ProjectionX(TString::Format("refoldedZ_iter%d", iteration).Data(), ptTruthLowBin, ptTruthHighBin - 1);
    refoldedZ[iter]->SetTitle(TString::Format("Refolded z iteration %d, #it{p}_{T, ch. jet}: %.0f - %.0f", iteration, ptDetLow, ptDetHigh).Data());
    refoldedZ[iter]->Write(refoldedZ[iter]->GetName(), TObject::kOverwrite);
    refoldedPt[iter] = (TH1D*)refolded[iter]->ProjectionY(TString::Format("refoldedPt_iter%d", iteration).Data());
    refoldedPt[iter]->SetTitle(TString::Format("Refolded pt iteration %d; #it{p}_{T, ch. jet}", iteration).Data());
    refoldedPt[iter]->Write(refoldedPt[iter]->GetName(), TObject::kOverwrite);

    // Ratio histograms
    ptUnfoldedOverTruth[iter] =
      (TH1D*)HistRatio(unfoldedPt[iter], testPtTruth,
                       TString::Format("ptUnfoldedOverTruth_iter%d", iteration).Data(),
                       TString::Format("Test pt unfolded/truth, iter %d; #it{p}_{T, ch. jet}; #frac{Unfolded}{Truth}", iteration).Data()
                      );
    zUnfoldedOverTruth[iter] =
      (TH1D*)HistRatio(unfoldedZ[iter], testZTruth,
                       TString::Format("zUnfoldedOverTruth_iter%d", iteration).Data(),
                       TString::Format("Test z unfolded/truth, iter %d, #it{p}_{T, ch. jet}: %.0f - %.0f; #it{z}; #frac{Unfolded}{Truth}", iteration, ptTruthLow, ptTruthHigh).Data()
                      );
    ptRefoldedOverDetector[iter] =
      (TH1D*)HistRatio(refoldedPt[iter], testPtDetector,
                       TString::Format("ptRefoldedOverDetector_iter%d", iteration).Data(),
                       TString::Format("Test pt refolded/detector, iter %d; #it{p}_{T, ch. jet}; #frac{Refolded}{Detector}", iteration).Data()
                      );
    zRefoldedOverDetector[iter] =
      (TH1D*)HistRatio(refoldedZ[iter], testZDetector,
                       TString::Format("zRefoldedOverDetector_iter%d", iteration).Data(),
                       TString::Format("Test z refolded/detector, iter %d, #it{p}_{T, ch. jet}: %.0f - %.0f; #it{z}; #frac{Refolded}{Detector}", iteration, ptDetLow, ptDetHigh).Data()
                      );

    ptUnfoldedOverTruth[iter]->Write(ptUnfoldedOverTruth[iter]->GetName(), TObject::kOverwrite);
    zUnfoldedOverTruth[iter]->Write(zUnfoldedOverTruth[iter]->GetName(), TObject::kOverwrite);
    ptRefoldedOverDetector[iter]->Write(ptRefoldedOverDetector[iter]->GetName(), TObject::kOverwrite);
    zRefoldedOverDetector[iter]->Write(zRefoldedOverDetector[iter]->GetName(), TObject::kOverwrite);

    // Diff histograms
    ptUnfoldedMinusTruth[iter] =
      (TH1D*)HistDiff(unfoldedPt[iter], testPtTruth,
                       TString::Format("ptUnfoldedMinusTruth_iter%d", iteration).Data(),
                       TString::Format("Test pt unfolded - truth, iter %d; #it{p}_{T, ch. jet}; Unfolded - Truth", iteration).Data()
                      );
    zUnfoldedMinusTruth[iter] =
      (TH1D*)HistDiff(unfoldedZ[iter], testZTruth,
                       TString::Format("zUnfoldedMinusTruth_iter%d", iteration).Data(),
                       TString::Format("Test z unfolded - truth, iter %d, #it{p}_{T, ch. jet}: %.0f - %.0f; #it{z}; Unfolded - Truth", iteration, ptTruthLow, ptTruthHigh).Data()
                      );
    ptRefoldedMinusDetector[iter] =
      (TH1D*)HistDiff(refoldedPt[iter], testPtDetector,
                       TString::Format("ptRefoldedMinusDetector_iter%d", iteration).Data(),
                       TString::Format("Test pt refolded - detector, iter %d; #it{p}_{T, ch. jet}; Refolded - Detector", iteration).Data()
                      );
    zRefoldedMinusDetector[iter] =
      (TH1D*)HistDiff(refoldedZ[iter], testZDetector,
                       TString::Format("zRefoldedMinusDetector_iter%d", iteration).Data(),
                       TString::Format("Test z refolded - detector, iter %d, #it{p}_{T, ch. jet}: %.0f - %.0f; #it{z}; Refolded - Detector", iteration, ptDetLow, ptDetHigh).Data()
                      );

    ptUnfoldedMinusTruth[iter]->Write(ptUnfoldedMinusTruth[iter]->GetName(), TObject::kOverwrite);
    zUnfoldedMinusTruth[iter]->Write(zUnfoldedMinusTruth[iter]->GetName(), TObject::kOverwrite);
    ptRefoldedMinusDetector[iter]->Write(ptRefoldedMinusDetector[iter]->GetName(), TObject::kOverwrite);
    zRefoldedMinusDetector[iter]->Write(zRefoldedMinusDetector[iter]->GetName(), TObject::kOverwrite);

    // Diff ratio histograms
    ptTruthDiffRatio[iter]
      = (TH1D*)HistRatio(ptUnfoldedMinusTruth[iter], testPtTruth,
                         TString::Format("ptTruthDiffRatio_iter%d", iteration).Data(),
                         TString::Format("Test pt (unfolded - truth)/truth, iter %d; #it{p}_{T, ch. jet}; #frac{Unfolded - Truth}{Truth}", iteration).Data()
                        );
    zTruthDiffRatio[iter]
      = (TH1D*)HistRatio(zUnfoldedMinusTruth[iter], testZTruth,
                         TString::Format("zTruthDiffRatio_iter%d", iteration).Data(),
                         TString::Format("Test z (unfolded - truth)/truth, iter %d, #it{p}_{T, ch. jet}: %.0f - %.0f; #it{z}; #frac{Unfolded - Truth}{Truth}", iteration, ptTruthLow, ptTruthHigh).Data()
                        );
    ptDetectorDiffRatio[iter]
      = (TH1D*)HistRatio(ptRefoldedMinusDetector[iter], testPtDetector,
                         TString::Format("ptDetectorDiffRatio_iter%d", iteration).Data(),
                         TString::Format("Test pt (refolded - detector)/detector, iter %d; #it{z}; #frac{Refolded - Detector}{Detector}", iteration).Data()
                        );
    zDetectorDiffRatio[iter]
      = (TH1D*)HistRatio(zRefoldedMinusDetector[iter], testZDetector,
                         TString::Format("zDetectorDiffRatio_iter%d", iteration).Data(),
                         TString::Format("Test z (refolded - detector)/detector, iter %d, #it{p}_{T, ch. jet}: %.0f - %.0f; #it{z}; #frac{Refolded - Detector}{Detector}", iteration, ptDetLow, ptDetHigh).Data()
                        );

    ptTruthDiffRatio[iter]->Write(ptTruthDiffRatio[iter]->GetName(), TObject::kOverwrite);
    zTruthDiffRatio[iter]->Write(zTruthDiffRatio[iter]->GetName(), TObject::kOverwrite);
    ptDetectorDiffRatio[iter]->Write(ptDetectorDiffRatio[iter]->GetName(), TObject::kOverwrite);
    zDetectorDiffRatio[iter]->Write(zDetectorDiffRatio[iter]->GetName(), TObject::kOverwrite);

    if (!makePlots) { continue; }
    // ----------------------- Plotting -----------------------
    TCanvas* spectraCanvas = new TCanvas(TString::Format("spectraCanvas_iter%d", iteration).Data(),
                                         TString::Format("spectraCanvas_iter%d", iteration).Data(),
                                         2000, 1000);
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
    unfoldedPt[iter]->Draw();
    auto pad6 = spectraCanvas->cd(6);
    pad6->SetLogy();
    refoldedPt[iter]->Draw();
    auto pad7 = spectraCanvas->cd(7);
    pad7->SetLogy();
    unfoldedZ[iter]->Draw();
    auto pad8 = spectraCanvas->cd(8);
    pad8->SetLogy();
    refoldedZ[iter]->Draw();
    spectraCanvas->SaveAs(TString::Format("closureTest-spectra-iter%d.pdf", iteration).Data());

    TCanvas* covCanvas = new TCanvas(TString::Format("covCanvas_iter%d", iteration).Data(),
                                     TString::Format("covCanvas_iter%d", iteration).Data(),
                                     1000, 1000);
    covCanvas->Divide(2);
    covCanvas->cd(1);
    covmat[iter]->Draw("colz");
    covCanvas->cd(2);
    pearson[iter]->Draw("colz");
    covCanvas->SaveAs(TString::Format("cov-iter%d.pdf", iteration).Data());

    TCanvas* ratioCanvas = new TCanvas(TString::Format("ratioCanvas_iter%d", iteration).Data(),
                                       TString::Format("ratioCanvas_iter%d", iteration).Data(),
                                       1000, 1000);
    ratioCanvas->Divide(2, 2);
    ratioCanvas->cd(1);
    ptUnfoldedOverTruth[iter]->Draw(drawOption.c_str());
    ratioCanvas->cd(2);
    zUnfoldedOverTruth[iter]->Draw(drawOption.c_str());
    ratioCanvas->cd(3);
    ptRefoldedOverDetector[iter]->Draw(drawOption.c_str());
    ratioCanvas->cd(4);
    zRefoldedOverDetector[iter]->Draw(drawOption.c_str());
    ratioCanvas->SaveAs(TString::Format("closureTest-ratio-iter%d.pdf", iteration).Data());

    TCanvas* diffCanvas = new TCanvas(TString::Format("diffCanvas_iter%d", iteration).Data(),
                                      TString::Format("diffCanvas_iter%d", iteration).Data(),
                                      1000, 1000);
    diffCanvas->Divide(2, 2);
    diffCanvas->cd(1);
    ptUnfoldedMinusTruth[iter]->Draw(drawOption.c_str());
    diffCanvas->cd(2);
    zUnfoldedMinusTruth[iter]->Draw(drawOption.c_str());
    diffCanvas->cd(3);
    ptRefoldedMinusDetector[iter]->Draw(drawOption.c_str());
    diffCanvas->cd(4);
    zRefoldedMinusDetector[iter]->Draw(drawOption.c_str());
    diffCanvas->SaveAs(TString::Format("closureTest-diff-iter%d.pdf", iteration).Data());

    TCanvas* diffratioCanvas = new TCanvas(TString::Format("diffratioCanvas_iter%d", iteration).Data(),
                                           TString::Format("diffratioCanvas_iter%d", iteration).Data(),
                                           1000, 1000);
    diffratioCanvas->Divide(2, 2);
    diffratioCanvas->cd(1);
    ptTruthDiffRatio[iter]->Draw(drawOption.c_str());
    diffratioCanvas->cd(2);
    zTruthDiffRatio[iter]->Draw(drawOption.c_str());
    diffratioCanvas->cd(3);
    ptDetectorDiffRatio[iter]->Draw(drawOption.c_str());
    diffratioCanvas->cd(4);
    zDetectorDiffRatio[iter]->Draw(drawOption.c_str());
    diffratioCanvas->SaveAs(TString::Format("closureTest-diffratio-iter%d.pdf", iteration).Data());
  } // for iter
  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
} // ClosureTest


// This task:
// * Loads the response from the task output
// * Creates the RooUnfoldResponse object used in unfolding

void CreateRooUnfoldResponse(string inName = "AnalysisResults.root", Double_t PT_MIN = 20.0,
                             Double_t Z_MIN = 0, Double_t Z_MAX = 1, Int_t BINWIDTHPT = 5){

  //Open the file and get the list where the response matrix is stored in THnSparse format
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  TDirectory* fragDir = (TDirectory*)inFile->Get("jet-fragmentation");
  if(!fragDir){
    std::cout << "Directory " << fragDir << " not found. Aborting program." << std::endl;
    return;
  }
  TDirectory* matchDir = (TDirectory*)fragDir->Get("matching");
  if(!matchDir){
    std::cout << "Directory " << matchDir << " not found. Aborting program." << std::endl;
    return;
  }
  TDirectory* matchJetsDir = (TDirectory*)matchDir->Get("jets");
  if(!matchJetsDir){
    std::cout << "Directory " << matchJetsDir << " not found. Aborting program." << std::endl;
    return;
  }

  string responseMatrixName = "matchDetJetPtTrackProjPartJetPtTrackProj";
  THnSparseF* responseMatrix = (THnSparseF*)matchJetsDir->Get(TString::Format("%s", responseMatrixName.c_str()).Data());
  if(!responseMatrix){
    std::cout << "Response matrix " << responseMatrix << " not found. Aborting program." << std::endl;
    return;
  }
  responseMatrix->Print();

  Int_t ptDetAxis   = 0;
  Int_t zDetAxis    = 1;
  Int_t ptTruthAxis = 2;
  Int_t zTruthAxis  = 3;

  // True and measured spectra as constructed from the response matrix
  TH2F* responseMatrixTruthProjection = (TH2F*)responseMatrix->Projection(zTruthAxis, ptTruthAxis, "E");
  responseMatrixTruthProjection->SetName("responseMatrixTruthProjection");
  TH2F* responseMatrixDetectorProjection = (TH2F*)responseMatrix->Projection(zDetAxis, ptDetAxis, "E");
  responseMatrixDetectorProjection->SetName("responseMatrixDetectorProjection");
  // Kinematic range
  Double_t ptMinDet    = PT_MIN;         // Lower Pt limit for detector level, Perhaps we can change this to evaluate systematics
  Double_t ptMaxDet    = 100.0;          // Upper Pt limit for detector level
  Double_t ptMinTruth  = 5.;             // Lower Pt limit for truth level, TODO: is this a suitable value?
  Double_t ptMaxTruth  = ptMaxDet + 50;  // Upper Pt limit for truth level, higher than detector level to allow for feed-in of jets
  Double_t zMinDet     = Z_MIN;          // Lower z limit for detector level
  Double_t zMaxDet     = Z_MAX;          // Upper z limit for detector level
  Double_t zMinTruth   = Z_MIN;          // Lower z limit for truth level
  Double_t zMaxTruth   = Z_MAX;          // Upper z limit for truth level

  // Rebinning settings: (detector level, particle level)
  Int_t    nPtBinWidth[2] = { 5, 5 };
  Double_t ptmin[2]       = { ptMinDet, ptMinTruth };
  Double_t ptmax[2]       = { ptMaxDet, ptMaxTruth };
  Int_t    nBinsPt[2]     = { int( (ptmax[0]-ptmin[0]) / nPtBinWidth[0]),
                              int( (ptmax[1]-ptmin[1]) / nPtBinWidth[1]) };

  Double_t nZBinWidth[2] = { 5e-2, 5e-2 };
  Double_t zmin[2]       = { zMinDet, zMinTruth };
  Double_t zmax[2]       = { zMaxDet, zMaxTruth };
  Int_t    nBinsZ[2]     = { int( (zmax[0] - zmin[0]) / nZBinWidth[0]),
                             int( (zmax[1] - zmin[1]) / nZBinWidth[1]) };

  // Create histograms used to construct RooUnfoldResponse and check what it actually does
  string detectorTitle = "Detector level #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}";
  string truthTitle = "Particle level #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}";

  // Detector level - Rebinned projection of response matrix
  TH2F* hDetector         = new TH2F("hDetector", TString::Format("%s", detectorTitle.c_str()).Data(),
                                     nBinsPt[0], ptmin[0], ptmax[0], nBinsZ[0], zmin[0], zmax[0]);
  // Detector level - In Detector and in Truth acceptance (response->Fill)
  TH2F* hDetectorMeasured = new TH2F("hDetectorMeasured", TString::Format("Measured %s", detectorTitle.c_str()).Data(),
                                     nBinsPt[0], ptmin[0], ptmax[0], nBinsZ[0], zmin[0], zmax[0]);
  // Detector level - In Truth, but out of Detector acceptance (response->Miss) (should be empty)
  TH2F* hDetectorMissed   = new TH2F("hDetectorMissed", TString::Format("Missed %s", detectorTitle.c_str()).Data(),
                                     nBinsPt[0], ptmin[0], ptmax[0], nBinsZ[0], zmin[0], zmax[0]);
  // Detector level - In Detector, but out of Truth acceptance (response->Fake)
  TH2F* hDetectorFake     = new TH2F("hDetectorFake", TString::Format("Fake %s", detectorTitle.c_str()).Data(),
                                     nBinsPt[0], ptmin[0], ptmax[0], nBinsZ[0], zmin[0], zmax[0]);

  // Truth level - Rebinned projection of response matrix
  TH2F* hTruth            = new TH2F("hTruth", TString::Format("%s", truthTitle.c_str()).Data(),
                                     nBinsPt[1], ptmin[1], ptmax[1], nBinsZ[1], zmin[1], zmax[1]);
  // Truth level - In Detector and in Truth acceptance (response->Fill)
  TH2F* hTruthMeasured    = new TH2F("hTruthMeasured", TString::Format("Measured %s", truthTitle.c_str()).Data(),
                                     nBinsPt[1], ptmin[1], ptmax[1], nBinsZ[1], zmin[1], zmax[1]);
  // Truth level - In Truth, but out of Detector acceptance (response->Miss)
  TH2F* hTruthMissed      = new TH2F("hTruthMissed", TString::Format("Missed %s", truthTitle.c_str()).Data(),
                                     nBinsPt[1], ptmin[1], ptmax[1], nBinsZ[1], zmin[1], zmax[1]);
  // Truth level - In Detector, but out of Truth acceptance (response->Fake) (should be empty)
  TH2F* hTruthFake        = new TH2F("hTruthFake", TString::Format("Fake %s", truthTitle.c_str()).Data(),
                                     nBinsPt[1], ptmin[1], ptmax[1], nBinsZ[1], zmin[1], zmax[1]);

  // Rebinning detector and truth
  // These must be done separately, as they have different bin numbers and sizes
  for (Int_t ix = 1; ix <= hDetector->GetNbinsX(); ix++) {
    Double_t xMin    = hDetector->GetXaxis()->GetBinLowEdge(ix);
    Double_t xMax    = hDetector->GetXaxis()->GetBinUpEdge(ix);
    Int_t    xBinLow = responseMatrixDetectorProjection->GetXaxis()->FindBin(xMin + 0.000001);
    Int_t    xBinUp  = responseMatrixDetectorProjection->GetXaxis()->FindBin(xMax - 0.000001);
    for (Int_t iy = 1; iy <= hDetector->GetNbinsY(); iy++) {
      Double_t yMin    = hDetector->GetYaxis()->GetBinLowEdge(iy);
      Double_t yMax    = hDetector->GetYaxis()->GetBinUpEdge(iy);
      Int_t    yBinLow = responseMatrixDetectorProjection->GetYaxis()->FindBin(yMin + 0.000001);
      Int_t    yBinUp  = responseMatrixDetectorProjection->GetYaxis()->FindBin(yMax - 0.000001);

      Double_t binError = 0.;
      Double_t binContent = responseMatrixDetectorProjection->IntegralAndError(xBinLow, xBinUp, yBinLow, yBinUp, binError);
      hDetector->SetBinContent(ix, iy, binContent);
      hDetector->SetBinError(ix, iy, binError);
    } // for iy
  }   // for ix
  for (Int_t ix = 1; ix <= hTruth->GetNbinsX(); ix++) {
    Double_t xMin    = hTruth->GetXaxis()->GetBinLowEdge(ix);
    Double_t xMax    = hTruth->GetXaxis()->GetBinUpEdge(ix);
    Int_t    xBinLow = responseMatrixTruthProjection->GetXaxis()->FindBin(xMin + 0.000001);
    Int_t    xBinUp  = responseMatrixTruthProjection->GetXaxis()->FindBin(xMax - 0.000001);
    for (Int_t iy = 1; iy <= hTruth->GetNbinsY(); iy++) {
      Double_t yMin    = hTruth->GetYaxis()->GetBinLowEdge(iy);
      Double_t yMax    = hTruth->GetYaxis()->GetBinUpEdge(iy);
      Int_t    yBinLow = responseMatrixTruthProjection->GetYaxis()->FindBin(yMin + 0.000001);
      Int_t    yBinUp  = responseMatrixTruthProjection->GetYaxis()->FindBin(yMax - 0.000001);

      Double_t binError = 0.;
      Double_t binContent = responseMatrixTruthProjection->IntegralAndError(xBinLow, xBinUp, yBinLow, yBinUp, binError);
      hTruth->SetBinContent(ix, iy, binContent);
      hTruth->SetBinError(ix, iy, binError);
    } // for iy
  }   // for ix

  // Create RooUnfoldResponse and fill it
  hDetector->Print();
  hTruth->Print();
  RooUnfoldResponse *ruResponse = new RooUnfoldResponse("Response", "RM");
  ruResponse->Setup(hDetector, hTruth);

  Int_t* coord = new Int_t[responseMatrix->GetNdimensions()]; //Carries the bin coordinates
  Int_t  nBins  = responseMatrix->GetNbins();
  for(Int_t bin = 0; bin < nBins; bin++) {
    Double_t w          = responseMatrix->GetBinContent(bin, coord);
    Double_t ptTruth    = responseMatrix->GetAxis(ptTruthAxis)->GetBinCenter(coord[ptTruthAxis]);
    Double_t zTruth     = responseMatrix->GetAxis(zTruthAxis)->GetBinCenter(coord[zTruthAxis]);
    Double_t ptMeasured = responseMatrix->GetAxis(ptDetAxis)->GetBinCenter(coord[ptDetAxis]);
    Double_t zMeasured  = responseMatrix->GetAxis(zDetAxis)->GetBinCenter(coord[zDetAxis]);

    bool inAcceptanceDetector = false;
    bool inAcceptanceTruth    = false;
    inAcceptanceDetector = (zMeasured  >= zmin[0]) * (zMeasured  < zmax[0]) * (ptMeasured >= ptmin[0]) * (ptMeasured < ptmax[0]);
    inAcceptanceTruth = (zTruth >= zmin[1]) * (zTruth < zmax[1]) * (ptTruth >= ptmin[1]) * (ptTruth < ptmax[1]);

    if(inAcceptanceDetector && inAcceptanceTruth) {
      ruResponse->Fill(ptMeasured, zMeasured, ptTruth, zTruth, w);
      hDetectorMeasured->Fill(ptMeasured, zMeasured, w);
      hTruthMeasured->Fill(ptTruth, zTruth, w);
    }
    else if (inAcceptanceTruth) {
      ruResponse->Miss(ptTruth, zTruth, w);
      hDetectorMissed->Fill(ptMeasured, zMeasured, w);
      hTruthMissed->Fill(ptTruth, zTruth, w);
    }
    else if (inAcceptanceDetector) {
      ruResponse->Fake(ptMeasured, zMeasured, w);
      hDetectorFake->Fill(ptMeasured, zMeasured, w);
      hTruthFake->Fill(ptTruth, zTruth, w);
    }
  } // for each bin
  TH2F* h2TrueFull = (TH2F*)hTruthMeasured->Clone("h2TrueFull");
  h2TrueFull->Add(hTruthMissed, 1);
  TH2F* h2MeasFull = (TH2F*)hDetectorMeasured->Clone("h2MeasFull");
  h2MeasFull->Add(hDetectorMissed, 1);

  delete [] coord;

  TFile* outFile = new TFile("RooUnfoldResponse.root", "RECREATE");
  responseMatrix->Write("responseMatrix");
  responseMatrixTruthProjection->Write();
  responseMatrixDetectorProjection->Write();
  ruResponse->Write();

  hDetector->Write();
  hDetectorMeasured->Write();
  hDetectorMissed->Write();
  hDetectorFake->Write();
  h2MeasFull->Write();

  hTruth->Write();
  hTruthMeasured->Write();
  hTruthMissed->Write();
  hTruthFake->Write();
  h2TrueFull->Write();

  outFile->Write();
  outFile->Close();
}

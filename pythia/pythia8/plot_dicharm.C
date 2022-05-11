{						
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetLineWidth(2);
  gStyle->SetFrameLineWidth(2);

  const Float_t ptMin = 5.0;
  //TFile *fin = TFile::Open("/dcache/alice/marcovl/pythia8charm/PythiaDiCharm_5457514.burrell.nikhef.nl.root");
  TFile *fin = TFile::Open("PythiaDiCharmAllNew.root");
  const Char_t *label = "Pythia8 hardQCD D-#bar{D}";
  TString ptlabel; ptlabel.Form("p_{T, D} > %.1f GeV/c",ptMin);

  auto hPt1Pt2EtaCharm = dynamic_cast<TH3F*>(fin->Get("hPt1Pt2EtaCharm"));
  Int_t ptBin = hPt1Pt2EtaCharm->GetXaxis()->FindBin(ptMin+0.001);
  auto hPtPt = hPt1Pt2EtaCharm->Project3D("xy");
  hPt1Pt2EtaCharm->GetXaxis()->SetRange(ptBin,hPt1Pt2EtaCharm->GetXaxis()->GetNbins()+1);
  auto hEta = hPt1Pt2EtaCharm->Project3D("z_1");
  hEta->SetMinimum(0);
  // TODO add case with D0bar > pt; D0 < pt

  hPt1Pt2EtaCharm->GetYaxis()->SetRange(ptBin,hPt1Pt2EtaCharm->GetYaxis()->GetNbins()+1);
  auto hEta2 = hPt1Pt2EtaCharm->Project3D("z_2");
  hEta2->SetLineColor(2);

  Int_t nBinsEta = hEta->GetNbinsX()/2;
  TH1F *hYieldEtaAcc = new TH1F("hYieldEtaAcc","Yield vs eta acc;#eta_{max,acc};dN_{D#bar{D}}/d#eta",nBinsEta,0,hEta->GetXaxis()->GetXmax());
  for (Int_t iEta = 0; iEta < nBinsEta; iEta++) {
    Int_t binMin = nBinsEta-iEta;
    Int_t binMax = nBinsEta+iEta+1;
    Float_t yld = hEta->Integral(binMin, binMax)/(hEta->GetXaxis()->GetBinUpEdge(binMax)-hEta->GetBinLowEdge(binMin));
    cout << "eta " << hEta->GetBinLowEdge(binMin) << " - " << hEta->GetXaxis()->GetBinUpEdge(binMax) << " yield " << yld << endl;
    hYieldEtaAcc->Fill(hEta->GetBinCenter(binMax), yld);
  }

  TCanvas *c1 = new TCanvas("c1","c1: eta max");
  hEta->Draw();
  hEta2->Draw("same");
  TLatex *ltx = new TLatex;
  ltx->SetNDC();
  ltx->SetTextFont(42);
  ltx->DrawLatex(0.14,0.82,label);

  TCanvas *c1b = new TCanvas("c1b","c1b: eta acc");
  hYieldEtaAcc->SetMinimum(0);
  hYieldEtaAcc->SetMaximum(800);
  hYieldEtaAcc->Draw();
  ltx->DrawLatex(0.14,0.82,label);
  ltx->DrawLatex(0.14,0.75,ptlabel);

  auto hPt1Pt2DEtaCharm = dynamic_cast<TH3F*>(fin->Get("hPt1Pt2DEtaCharm"));
  hPt1Pt2DEtaCharm->GetXaxis()->SetRange(ptBin,hPt1Pt2DEtaCharm->GetXaxis()->GetNbins()+1);
  auto hDEta = hPt1Pt2DEtaCharm->Project3D("z_1");
  // TODO add case with D0bar > pt; D0 < pt
  
  TCanvas *c2 = new TCanvas("c2","c2: delta eta");
  hDEta->Draw();
  ltx->DrawLatex(0.14,0.82,label);
  ltx->DrawLatex(0.14,0.75,ptlabel);

  auto hPt1Pt2DPhiCharm = dynamic_cast<TH3F*>(fin->Get("hPt1Pt2DPhiCharm"));
  hPt1Pt2DPhiCharm->GetXaxis()->SetRange(ptBin,hPt1Pt2DPhiCharm->GetXaxis()->GetNbins()+1);
  auto hDPhi = hPt1Pt2DPhiCharm->Project3D("z_1");
  // TODO add case with D0bar > pt; D0 < pt
  
  TCanvas *c3 = new TCanvas("c3","c3: delta phi");
  hDPhi->SetMinimum(0);
  hDPhi->Draw();
  ltx->DrawLatex(0.14,0.82,label);
  ltx->DrawLatex(0.14,0.75,ptlabel);

  TCanvas *c4 = new TCanvas("c4","c4: pt vs pt");
  c4->SetLogz();
  hPtPt->Draw("colz");
  ltx->DrawLatex(0.5,0.85,label);
}

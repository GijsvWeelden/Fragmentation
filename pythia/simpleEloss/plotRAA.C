{
// Basic approach: fixed fraction eloss for each jet; divided over particles
//TFile *fin = TFile::Open("Output_jetbyjet_relDE015.root");
// Jet Eloss drawn from exponential distribution; particle Eloss derived, no fluct
TFile *fin = TFile::Open("Output_jetbyjet_relDE015_expo.root");
// Jet Eloss drawn from exponential distribution; particle Eloss also 
//TFile *fin = TFile::Open("Output_jetbyjet_relDE015_expo_expo.root");
// Energy per particle: 3 or 6 GeV fixed
//TFile *fin = TFile::Open("Output_byparticle_DE3.root");
// Energy per particle: 4 GeV with exp fluctuations
//TFile *fin = TFile::Open("Output_byparticle_DE4_exp.root");
// Energy per particle: relative: 0.3 with exp fluctuations
//TFile *fin = TFile::Open("Output_byparticle_relDE03_exp.root");
// Energy per particle: relative: 0.15 with exp fluctuations
//TFile *fin = TFile::Open("Output_byparticle_relDE015_exp.root");
// Energy per particle: relative: 0.25 with max DE=6 and exp fluctuations
TFile *fin = TFile::Open("Output_byparticle_relDE025_maxDE6_exp_all.root");

cout << "Using file " << fin->GetName() << endl;
TCanvas *c1 = new TCanvas("c1","c1: spectra",500,700);
TCanvas *c1 = new TCanvas("c1","c1: spectra",500,700);
c1->SetLogy();

auto hJetEtaPt = dynamic_cast<TH2D*>(fin->Get("hJetEtaPt"));
auto hJetPt = hJetEtaPt->ProjectionY();
hJetPt->Draw("");
auto hJetEtaPtLoss = dynamic_cast<TH2D*>(fin->Get("hJetEtaPtLoss"));
hJetEtaPtLoss->SetLineColor(2);
hJetEtaPtLoss->SetMarkerColor(2);
auto hJetPtLoss = hJetEtaPtLoss->ProjectionY();
hJetPtLoss->Draw("same");

auto hJetRAA = new TH1D(*hJetPtLoss);
hJetRAA->SetName("hJetRAA");
hJetRAA->Divide(hJetPt);

auto hTrackEtaPt = dynamic_cast<TH2D*>(fin->Get("hTrackEtaPt"));
auto hTrackPt = hTrackEtaPt->ProjectionY();
hTrackPt->Draw("same");

auto hTrackEtaPtLoss = dynamic_cast<TH2D*>(fin->Get("hTrackEtaPtLoss"));
hTrackEtaPtLoss->SetLineColor(2);
hTrackEtaPtLoss->SetMarkerColor(2);
auto hTrackPtLoss = hTrackEtaPtLoss->ProjectionY();
hTrackPtLoss->Draw("same");

auto hTrackRAA = new TH1D(*hTrackPtLoss);
hTrackRAA->SetName("hTrackRAA");
hTrackRAA->Divide(hTrackPt);

TCanvas *c2 = new TCanvas("c2","c2: RAA",700,500);
gPad->DrawFrame(0,0,200,2);
hJetRAA->Draw("same");
hTrackRAA->SetLineColor(4);
hTrackRAA->SetMarkerColor(4);
hTrackRAA->Draw("same");

TLegend *leg = new TLegend(0.5,0.6,0.85,0.8);
leg->SetBorderSize(0);
leg->AddEntry(hTrackRAA,"particles","l");
leg->AddEntry(hJetRAA,"jets","l");
leg->Draw();
}

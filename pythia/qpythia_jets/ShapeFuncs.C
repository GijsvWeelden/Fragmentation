Float_t angularity(AliEmcalJet *jet, TClonesArray *particles) {
  if (!jet->GetNumberOfTracks())
    return 0;
  Double_t den=0.;
  Double_t num = 0.;
  AliVParticle *vp1 = 0x0;
  for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
    vp1 = static_cast<TParticle*>(jet->TrackAt(i, particles));
    Double_t dphi = vp1->Phi()-jet->Phi();
    if(dphi<-1.*TMath::Pi()) dphi+=TMath::TwoPi();
    if(dphi>TMath::Pi()) dphi-=TMath::TwoPi();
    Double_t dr2 = (vp1->Eta()-jet->Eta())*(vp1->Eta()-jet->Eta()) + dphi*dphi;
    Double_t dr = TMath::Sqrt(dr2);
    num=num+vp1->Pt()*dr;
    den=den+vp1->Pt();
  }
  return num/den;
}


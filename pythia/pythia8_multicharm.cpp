#include <iostream>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"

#define nEvents 3 

using namespace Pythia8;

int main(int /*argc*/, char** /*argv*/)
{


//PYTHIA SETTINGS

	TString name;

	int mecorr=1;
	
	// Generator. Process selection. LHC initialization. Histogram.
	Pythia pythia;

	pythia.readString("Random:setSeed = on");
	pythia.readString("Random:seed = 0");

	pythia.readString("ProcessLevel::all = off"); // No event generation

	// Multicharm modes from Stefano
	pythia.readString("4444:oneChannel = 1 1 0 4432 211");
	pythia.readString("4444:onMode = off");
	pythia.readString("4444:onIfMatch = 4432 211");

	pythia.readString("4432:oneChannel = 1 1 0 4332 211");
	pythia.readString("4432:onMode = off");
	pythia.readString("4432:onIfMatch = 4332 211");
	
	pythia.readString("4332:oneChannel = 1 1 0 3334 211");
	pythia.readString("4332:onMode = off");
	pythia.readString("4332:onIfMatch = 3334 211");

	/*
	  // Omega and Xi decay
	pythia.readString("3334:onMode = off");
	pythia.readString("3334:onIfAll = 3112 321");

	pythia.readString("3112:onMode = off");
	pythia.readString("3112:onIfAll = 2212 311");
	*/
	pythia.readString("310:mayDecay  = off"); //K0s
	pythia.readString("3122:mayDecay = off"); //labda0
	pythia.readString("3112:mayDecay = off"); //sigma-
	pythia.readString("3212:mayDecay = off"); //sigma0
	pythia.readString("3222:mayDecay = off"); //sigma+
	pythia.readString("3312:mayDecay = off"); //xi-
	pythia.readString("3322:mayDecay = off"); //xi+
	pythia.readString("3334:mayDecay = off"); //omega-

	pythia.init();

        // Output histograms
	TFile* outFile = new TFile("PythiaMultiCharm.root","RECREATE");
        TH2F *hEtaPt = new TH2F("hEtaPt","Pt vs Eta for all particles;#eta;p_{T} (GeV/c)",40,-2,2,50,0,10);


	double mOmega = pythia.particleData.m0(4444);
        double tauOmega = pythia.particleData.tau0(4444);  // Pythia default is 0.1 mm = 100 micron
        cout << " Omega mass " << mOmega << ", tau0 " << tauOmega << endl;
        cout << " tau0 for Omegacc " << pythia.particleData.tau0(4432) << endl;
	//Begin event loop
	for(int iEvent = 0; iEvent < nEvents; iEvent++)
	{
		double fourvec[4];
		double pvec[3] = {1.0, 0, 0};
		double E = sqrt(mOmega*mOmega + pvec[0]*pvec[0] + pvec[1]*pvec[1] + pvec[2]*pvec[2]);
		pythia.event.reset();
		// first decay happens at vertex; how can we propogate?
		pythia.event.append(4444, 11, 0, 0, pvec[0], pvec[1], pvec[2], E, mOmega);

                pythia.moreDecays(); // do decay
		//if(!pythia.next()) continue;
		
                Double_t ptSumPythia = 0;
                Int_t nPartPythia = 0;
                int nPart = pythia.event.size();
		cout << "Found " << nPart << " particles " << endl;
		for(int iPart = 0; iPart < nPart; iPart++) 
		{
                    const Particle &part = pythia.event[iPart];
		    cout << " particle " << part.id() << " prod vtx " << part.xProd() << " " << part.yProd() << " " << part.zProd() << endl;
			double rprod = TMath::Sqrt(part.xProd()*part.xProd() + part.yProd()*part.yProd());
			double x0 = part.xProd();
			double y0 = part.yProd();
			double x1 = part.xProd()+1;
			double y1 = part.yProd()+1*part.py()/part.px();
			double d0 = fabs((x1-x0)*y0-x0*(y1-y0))/sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)); // Could define x,y1 such that they the distance is 1 mm?
                        cout << " rprod " << rprod*1000. << " d0 " << d0*1000. << " mum" << endl;
			if(part.isFinal())
			{

				hEtaPt->Fill(part.eta(),part.pT());
                                nPartPythia++;
			}
		}
                
                if ((iEvent%1000)==0) 
                   cout << "Pythia event: " << nPartPythia << " particles" << endl;
	}
//End event loop

   
    	outFile->Write();
        cout << "Histos written to file " << outFile->GetName() << endl;
        outFile->Close();
}

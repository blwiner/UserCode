#include "L1Ntuple.h"
#include "L1AnalysisDataFormat.h"
#include "hist.C"
#include "Style.C"

#include "TLegend.h"
#include "TMath.h"
#include "TText.h"
#include "TH2.h"
#include "TAxis.h"
#include "TString.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <set>

// -- Huge prescale value for seeds "for lower PU"
#define INFTY 10000



/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Notes:
  
    -> This needs to be run within the UserCode/L1TriggerDPG package
	 
	 -> In also requires:
	           L1AnalaysisDataFormat.h
				  getLumi_out_pixelCorrLumi_*PU_stdCorr.txt (for HPF Data)
				  
	 -> General running format (from UserCode/L1TriggerDPG/macro/ area)
	 
	    linux>  root initL1Analysis.C      (loading libraries etc) 
		 root>   .L L1Menu2015.C++          (compiles the script) 			  
       root>   RunL1_HFW(Bool_t calcThreshold=true,Bool_t useL1Extra=true,Int_t usedL1Menu=0,Float_t targetlumi=200,Int_t whichFileAndLumiToUse=1,int which_jet_seed_to_use=0, Int_t pMenu2015 = 0, Int_t pNevts = -1)
 
       Definition of input quantities:
                  calcThreshold = flag for whether rate vs threshold plots are made (uses a lot more CPU)
						useL1Extra    = flag for whether to use l1extra informaiton or quantities from GT
						targetlumi    = units of E32
						whichFileAndLumiToUse = Specifies input file (set with Switch Statement)  NOTE: This is ugly and should be improved.
						which_jet_seed_to_use = Also a file specifier for selecting different versions of the same sample  NOTE:  ditto
						pMenu2015     = Selects different trigger thresholds to use for the full menu emulation
						pNevts        = number of events to run over (-1 ==> run over all events in file)
						
   There are a lot of details that are not yet summarized here....
	
					

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
//HFW
//Cross
TH1F *h_SingleMu_ETM_byThreshold; 
TH2F *h2_SingleMu_ETM_byThreshold;
TH1F *h_SingleMu_CJet_byThreshold;
TH2F *h2_SingleMu_CJet_byThreshold;
TH1F *h_SingleEG_ETM_byThreshold; 
TH2F *h2_SingleEG_ETM_byThreshold;
TH1F *h_SingleEG_CJet_byThreshold;
TH2F *h2_SingleEG_CJet_byThreshold;

TH1F *h_Mu_EG_byThreshold;
TH1F *h_EG_Mu_byThreshold;
TH2F *h2_Mu_EG_byThreshold;
//Jets
TH1F *h_SingleJet_byThreshold;
TH1F *h_DoubleJet_byThreshold;
TH1F *h_QuadJetCentral_byThreshold;
TH1F *h_SingleTau_byThreshold;
TH1F *h_DoubleTau_byThreshold;
//Sums
TH1F *h_HTT_byThreshold;
TH1F *h_ETM_byThreshold;
//EGamma
TH1F *h_SingleEG_byThreshold;
TH1F *h_SingleIsoEG_byThreshold;
TH1F *h_DoubleEG_byThreshold;
TH2F *h2_DoubleEG_byThreshold;
//Muons
TH1F *h_SingleMu_byThreshold;
TH1F *h_DoubleMu_byThreshold;
TH2F *h2_DoubleMu_byThreshold;

// Plots for Trigger Quantities
TH1F *h_Mu_Nmu,      *h_Mu_Et,     *h_Mu_Eta,     *h_Mu_Phi;
TH1F *h_isoEG_Nele,  *h_isoEG_Et,  *h_isoEG_Eta,  *h_isoEG_Phi;
TH1F *h_nIsoEG_Nele, *h_nIsoEG_Et, *h_nIsoEG_Eta, *h_nIsoEG_Phi;
TH1F *h_CJet_Njet, *h_CJet_Et, *h_CJet_Eta, *h_CJet_Phi;
TH1F *h_FJet_Njet, *h_FJet_Et, *h_FJet_Eta, *h_FJet_Phi;
TH1F *h_TJet_Njet, *h_TJet_Et, *h_TJet_Eta, *h_TJet_Phi;
TH1F *h_Sum_ETT,   *h_Sum_ETM, *h_Sum_PhiETM;
TH1F *h_Sum_HTT,   *h_Sum_HTM, *h_Sum_PhiHTM;


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Int_t NPAGS = 6;
TH2F *cor_PAGS;
TH1F *h_PAGS_pure;
TH1F *h_PAGS_shared;

const Int_t N128 = 128;			// could be > 128 for "test seeds"
Int_t kOFFSET = 0;
Bool_t TheTriggerBits[N128];	// contains the emulated triggers for each event
TH1F *h_All, *h_Trig;		// one bin for each trigger. Fill bin i if event fires trigger i.
TH1F *h_Pure;		// one bin for each trigger. Fill bin i if event fires trigger i and NO OTHER TRIGGER.
TH2F *h_Corr;
TH1F *h_Cumm;

Int_t Menu2015 = 0;

// Methods to scale L1 jets for new HCAL LUTs and estimate the rate changes 

// correction by 5% overall (from HCAL January 2012)
Double_t CorrectedL1JetPtByFactor(Double_t JetPt, Bool_t theL1JetCorrection=false) {

	Double_t JetPtcorr = JetPt;

	if (theL1JetCorrection) {
		JetPtcorr = JetPt*1.05;
	}
	return JetPtcorr;
}

// correction by 8% for forward jets (from HCAL January 2012)
Double_t CorrectedL1FwdJetPtByFactor(Bool_t isFwdJet, Double_t JetPt, Bool_t theL1JetCorrection=false) {

	Double_t JetPtcorr = JetPt;

	if (theL1JetCorrection) {
		if (isFwdJet) { JetPtcorr = JetPt*1.08; }
	}
	return JetPtcorr;
}

// correction for HF bins (from HCAL January 2012)
Size_t   JetHFiEtabins   = 13;
Int_t    JetHFiEtabin[]  = {29,30,31,32,33,34,35,36,37,38,39,40,41};
Double_t JetHFiEtacorr[] = {0.982,0.962,0.952, 0.943,0.947,0.939, 0.938,0.935,0.934, 0.935,0.942,0.923,0.914};

Double_t CorrectedL1JetPtByHFtowers(Double_t JetiEta,Double_t JetPt, Bool_t theL1JetCorrection=false) {

	Double_t JetPtcorr   = JetPt;

	if (theL1JetCorrection) {
		Int_t    iJetiEtabin = 0;
		for (iJetiEtabin=0; iJetiEtabin<JetHFiEtabins; iJetiEtabin++) {
			if (JetHFiEtabin[iJetiEtabin]==JetiEta) {
				JetPtcorr = JetPt * (1+(1-JetHFiEtacorr[iJetiEtabin]));
			}
		}
	}
	return JetPtcorr;
}

// correction for RCT->GCT bins (from HCAL January 2012)
// HF from 29-41, first 3 HF trigger towers 3 iEtas, last highest eta HF trigger tower 4 iEtas; each trigger tower is 0.5 eta, RCT iEta from 0->21 (left->right)
Double_t JetRCTHFiEtacorr[]  = {0.965,0.943,0.936,0.929}; // from HF iEta=29 to 41 (smaller->higher HF iEta)

Double_t CorrectedL1JetPtByGCTregions(Double_t JetiEta,Double_t JetPt, Bool_t theL1JetCorrection=false) {

	Double_t JetPtcorr   = JetPt;

	if (theL1JetCorrection) {

		if ((JetiEta>=7 && JetiEta<=14)) {
			JetPtcorr = JetPt * 1.05;
		}

		if ((JetiEta>=4 && JetiEta<=6) || (JetiEta>=15 && JetiEta<=17)) {
			JetPtcorr = JetPt * 0.95;
		}

		if (JetiEta==0 || JetiEta==21) {
			JetPtcorr = JetPt * (1+(1-JetRCTHFiEtacorr[3]));
		}
		else if (JetiEta==1 || JetiEta==20) {
			JetPtcorr = JetPt * (1+(1-JetRCTHFiEtacorr[2]));
		}
		else if (JetiEta==2 || JetiEta==19) {
			JetPtcorr = JetPt * (1+(1-JetRCTHFiEtacorr[1]));
		}
		else if (JetiEta==3 || JetiEta==18) {
			JetPtcorr = JetPt * (1+(1-JetRCTHFiEtacorr[0]));
		}
	}

	return JetPtcorr;
}

// methods for the correlation conditions

size_t PHIBINS = 18;
Double_t PHIBIN[] = {10,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350};

size_t ETABINS = 23;
Double_t ETABIN[] = {-5.,-4.5,-4.,-3.5,
	-3.,-2.172,-1.74,-1.392,-1.044,-0.696,-0.348,
	0,
	0.348,0.696,1.044,1.392,1.74,2.172,3.,
	3.5,4.,4.5,5.};

size_t ETAMUBINS = 65;
Double_t ETAMU[] = { -2.45,-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45 };

Int_t etaMuIdx(Double_t eta) {
	size_t etaIdx = 0.;
	for (size_t idx=0; idx<ETAMUBINS; idx++) {
		if (eta>=ETAMU[idx] and eta<ETAMU[idx+1])
			etaIdx = idx;
	}
	return int(etaIdx);
}

Int_t etaINjetCoord(Double_t eta){
	size_t etaIdx = 0.;
	for (size_t idx=0; idx<ETABINS; idx++) {
		if (eta>=ETABIN[idx] and eta<ETABIN[idx+1])
			etaIdx = idx;
	}
	return int(etaIdx);
}

Double_t degree(Double_t radian) {
	if (radian<0)
		return 360.+(radian/TMath::Pi()*180.);
	else
		return radian/TMath::Pi()*180.;
}

Int_t phiINjetCoord(Double_t phi) {
	size_t phiIdx = 0;
	Double_t phidegree = degree(phi);
	for (size_t idx=0; idx<PHIBINS; idx++) {
		if (phidegree>=PHIBIN[idx] and phidegree<PHIBIN[idx+1])
			phiIdx = idx;
		else if (phidegree>=PHIBIN[PHIBINS-1] || phidegree<=PHIBIN[0])
			phiIdx = idx;
	}
	phiIdx = phiIdx + 1;
	if (phiIdx == 18)  phiIdx = 0;
	return int(phiIdx);
}

Bool_t correlateInPhi(Int_t jetphi, Int_t muphi, Int_t delta=1) {

	Bool_t correlationINphi = fabs(muphi-jetphi)<fabs(2 +delta-1) || fabs(muphi-jetphi)>fabs(PHIBINS-2 - (delta-1) );
	return correlationINphi;

}

Bool_t correlateInEta(Int_t mueta, Int_t jeteta, Int_t delta=1) {
	Bool_t correlationINeta = fabs(mueta-jeteta)<2 + delta-1;
	return correlationINeta;
}

// set the errors properly
void CorrectScale(TH1F* h, Float_t scal) {

	Int_t nbins = h -> GetNbinsX();

	for (Int_t i=1; i<= nbins; i++)  {
		Float_t val = h -> GetBinContent(i);
		Float_t er = sqrt(val);
		val = val * scal;
		er = er * scal;
		h -> SetBinContent(i,val);
		h -> SetBinError(i,er);
	}
}

class L1Menu2015 : public L1Ntuple {
	public :

	L1Menu2015(Int_t aL1Menu,Float_t aTargetLumi, Float_t aNumberOfUserdLumiSections, Float_t aLumiForThisSetOfLumiSections, std::string aL1NtupleFileName,Float_t aAveragePU, Float_t aZeroBiasPrescale,Bool_t aL1JetCorrection) : 
     	theL1Menu(aL1Menu), 
		theTargetLumi(aTargetLumi), 
		theNumberOfUserdLumiSections(aNumberOfUserdLumiSections),
		theLumiForThisSetOfLumiSections(aLumiForThisSetOfLumiSections),
		theL1NtupleFileName(aL1NtupleFileName),
		theAveragePU(aAveragePU),
		theZeroBiasPrescale(aZeroBiasPrescale),
		theL1JetCorrection(aL1JetCorrection)
		{}

	~L1Menu2015() {}

	Int_t theL1Menu;

	// The luminosity for which we want the rates, in units 1e32 (this sets also the right prescales for the 2 menus and some descriptions correctly).
	// 70. for 7e33, 50 for 5e33, etc. Use 70.001 for the "emergency columns" to the corresponding target luminosity.
	// For the moment we have pre-scales for 5e33,6e33,7e33 plus emergency pre-scales for 5e33 and 7e33
	Float_t theTargetLumi;

	// the setting below are/will be specific for each L1Ntuple file used
	Float_t theNumberOfUserdLumiSections;
	Float_t theLumiForThisSetOfLumiSections;
	std::string theL1NtupleFileName;
	Float_t theAveragePU;
	Float_t theZeroBiasPrescale;
	Bool_t theL1JetCorrection;

	void MyInit();
	void FilL1Bits();
	
	L1Analysis::L1AnalysisDataFormat myEvt_;

	std::map<std::string, int> Counts;
	std::map<std::string, int> Prescales;
	std::map<std::string, bool> Biased;

        std::map<std::string, int> BitMapping;
	
	typedef struct {
	   float primTh ;
	   float secTh;
	   float triTh;
	   float quadTh;
	   float etaCut;
	   int minQual;
	} trigPar;
	
	std::map<std::string, trigPar> trigParList;

	std::map<std::string, float> WeightsPAGs;


	void InsertInMenu(std::string L1name, Bool_t value);

	Int_t L1BitNumber(std::string l1name);

	Bool_t EvalMenu(double lumiWeight);
	void EvalThresh(double lumiWeight);

// -- Cross
	Bool_t Mu_EG(Float_t mucut, Float_t EGcut, Int_t minMuQual = 4 );
	Bool_t MuOpen_EG(Float_t mucut, Float_t EGcut );
	Bool_t Mu_JetCentral(Float_t mucut, Float_t jetcut );
	Bool_t Mu_DoubleJetCentral(Float_t mucut, Float_t jetcut );
	Bool_t Mu_JetCentral_LowerTauTh(Float_t mucut, Float_t jetcut, Float_t taucut );
	Bool_t Muer_JetCentral(Float_t mucut, Float_t jetcut, Float_t etacut = 2.1, Int_t minMuQual=4 );
	Bool_t Muer_JetCentral_LowerTauTh(Float_t mucut, Float_t jetcut, Float_t taucut );
	Bool_t Mu_HTT(Float_t mucut, Float_t HTcut );
	Bool_t Muer_ETM(Float_t mucut, Float_t ETMcut, Float_t etacut = 2.1, Int_t minMuQual=4 );
	Bool_t EG_FwdJet(Float_t EGcut, Float_t FWcut ) ;
	Bool_t EG_JetCentral(Float_t EGcut, Float_t jetcut );
	Bool_t EG_HT(Float_t EGcut, Float_t HTcut );
	Bool_t EG_DoubleJetCentral(Float_t EGcut, Float_t jetcut );
	Bool_t EG_ETM(Float_t EGcut, Float_t ETMcut );
	Bool_t DoubleEG_HT(Float_t EGcut, Float_t HTcut );
	Bool_t EGEta2p1_JetCentral(Float_t EGcut, Float_t jetcut);		// delta
	Bool_t EGEta2p1_JetCentral_LowTauTh(Float_t EGcut, Float_t jetcut, Float_t taucut);          // delta
	Bool_t IsoEGEta2p1_JetCentral_LowTauTh(Float_t EGcut, Float_t jetcut, Float_t taucut);          // delta
	Bool_t EGEta2p1_DoubleJetCentral(Float_t EGcut, Float_t jetcut);	// delta
	Bool_t EGEta2p1_DoubleJetCentral_TripleJetCentral(Float_t EGcut, Float_t jetcut2, Float_t jetcut3);   

	Bool_t HTT_HTM(Float_t HTTcut, Float_t HTMcut);
	Bool_t JetCentral_ETM(Float_t jetcut, Float_t ETMcut);
	Bool_t DoubleJetCentral_ETM(Float_t jetcut1, Float_t jetcut2, Float_t ETMcut);
	Bool_t DoubleMu_EG(Float_t mucut, Float_t EGcut );
	Bool_t Mu_DoubleEG(Float_t mucut, Float_t EGcut);

	Bool_t Muer_TripleJetCentral(Float_t mucut, Float_t jet1, Float_t jet2, Float_t jet3);
	Bool_t Mia(Float_t mucut, Float_t jet1, Float_t jet2);	// delta
	Bool_t Mu_JetCentral_delta(Float_t mucut, Float_t ptcut);	// delta
	Bool_t Mu_JetCentral_deltaOut(Float_t mucut, Float_t ptcut); // delta


// -- Jets 
	Bool_t SingleJet(Float_t cut);
	Bool_t SingleTauJet(Float_t cut, Float_t etaCut=4.5);
	Bool_t SingleJetCentral(Float_t cut);
	Bool_t DoubleJetCentral(Float_t cut1, Float_t cut2);
	Bool_t DoubleJet_Eta1p7_deltaEta4(Float_t cut1, Float_t cut2);
	Bool_t TripleJetCentral(Float_t cut1, Float_t cut2, Float_t cut3);
	Bool_t TripleJet_VBF(Float_t cut1, Float_t cut2, Float_t cut3);

	Bool_t QuadJetCentral(Float_t cut1, Float_t cut2, Float_t cut3, Float_t cut4);
	Bool_t DoubleTauJetEta(Float_t cut1, Float_t cut2, Float_t etaCut=4.5);

// -- Sums
	Bool_t ETT(Float_t ETTcut);
	Bool_t HTT(Float_t HTTcut);
	Bool_t ETM(Float_t ETMcut);

// -- Egamma
	Bool_t SingleEG(Float_t cut);
	Bool_t SingleEG_Eta(Float_t cut, Float_t etaCut=4.5);
	Bool_t SingleIsoEG_Eta(Float_t cut, Float_t etaCut=4.5);

	Bool_t DoubleEG(Float_t cut1, Float_t cut2);
	Bool_t TripleEG(Float_t cut1, Float_t cut2, Float_t cut3);

// -- Muons 
	Bool_t SingleMu(Float_t ptcut, Int_t qualmin=4);
	Bool_t SingleMuEta(Float_t ptcut, Float_t etaCut=2.1, Int_t qualmin=4);
	Bool_t DoubleMu(Float_t cut1, Float_t cut2, Int_t qualmin=4);	// on top of DoubleMu3
	Bool_t DoubleMuHighQEtaCut(Float_t ptcut, Float_t etacut);
	Bool_t TripleMu(Float_t cut1, Float_t cut2, Float_t cut3, Int_t qualmin);	// on top of DoubleMu3
	Bool_t DoubleMuXOpen(Float_t ptcut);	// on top of SingleMu7
	Bool_t Onia(Float_t ptcut1, Float_t ptcut2, Float_t etacut, Int_t delta);   

	void Loop(Bool_t calcThreshold, Bool_t useL1Extra, TString lsRunFile, const int n_events_=-1);
	void fillDataStructure(Bool_t UseL1Extra=true);

	private :

	Bool_t PhysicsBits[128];
	Bool_t first;

	Int_t insert_ibin;
	Bool_t insert_val[100];
	std::string insert_names[100];

	Int_t NBITS_TRIGS;


};


void L1Menu2015::fillDataStructure(Bool_t useL1Extra) {
   
	 // printf("Entering fillDataStucture \n");
	 myEvt_.Reset();

    if(useL1Extra) {

// Grab the iso first
		 for(unsigned int i=0; i<l1extra_->nIsoEm; i++) {


      	 myEvt_.Bxel.push_back(l1extra_->isoEmBx.at(i));
      	 myEvt_.Etel.push_back(l1extra_->isoEmEt.at(i));
      	 myEvt_.Phiel.push_back(phiINjetCoord(l1extra_->isoEmPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etael.push_back(etaINjetCoord(l1extra_->isoEmEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
       	 myEvt_.Isoel.push_back(true);

// Histogram Quantities          		 
          if(l1extra_->isoEmBx.at(i)==0) {
			    h_isoEG_Et->Fill(myEvt_.Etel.at(myEvt_.Nele));
				 h_isoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele));
				 h_isoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele));
			 }	 
			 myEvt_.Nele++;

		 }
		 h_isoEG_Nele->Fill(l1extra_->nIsoEm);	
		 
		 
		 for(unsigned int i=0; i<l1extra_->nNonIsoEm; i++) {

      	 myEvt_.Bxel.push_back(l1extra_->nonIsoEmBx.at(i));
      	 myEvt_.Etel.push_back(l1extra_->nonIsoEmEt.at(i));
      	 myEvt_.Phiel.push_back(phiINjetCoord(l1extra_->nonIsoEmPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etael.push_back(etaINjetCoord(l1extra_->nonIsoEmEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
			 myEvt_.Isoel.push_back(false);

// Histogram Quantities          			 
          if(l1extra_->nonIsoEmBx.at(i)==0) {
			    h_nIsoEG_Et-> Fill(myEvt_.Etel.at(myEvt_.Nele));
				 h_nIsoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele));
				 h_nIsoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele));
			 }	 
			 myEvt_.Nele++;
		 }	
		 h_nIsoEG_Nele->Fill(l1extra_->nNonIsoEm);	 
		 // printf("Number of Electrons in myEvt %i \n",myEvt_.Nele);

		 for(unsigned int i=0; i< l1extra_->nCenJets; i++) {

      	 
      	 myEvt_.Bxjet.push_back(l1extra_->cenJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1extra_->cenJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1extra_->cenJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1extra_->cenJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);
			 
// Histogram Quantities          			 
          if(l1extra_->cenJetBx.at(i)==0) {
			    h_CJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet));
				 h_CJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet));
				 h_CJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet));
			 }
			 myEvt_.Njet++;	 
		 }
		 h_CJet_Njet->Fill(l1extra_->nCenJets);
		 
		 
		 for(unsigned int i=0; i< l1extra_->nFwdJets; i++) {

      	 
      	 myEvt_.Bxjet.push_back(l1extra_->fwdJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1extra_->fwdJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1extra_->fwdJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1extra_->fwdJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(true);

// Histogram Quantities          			 
          if(l1extra_->fwdJetBx.at(i)==0) {
			    h_FJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet));
				 h_FJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet));
				 h_FJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet));
			 }	 
			 myEvt_.Njet++;
		 }
		 h_FJet_Njet->Fill(l1extra_->nFwdJets);
		 
		 for(unsigned int i=0; i< l1extra_->nTauJets; i++) {

      	 myEvt_.Bxjet.push_back(l1extra_->tauJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1extra_->tauJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1extra_->tauJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1extra_->tauJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(true);
      	 myEvt_.Fwdjet.push_back(false);

// Histogram Quantities          			 
          if(l1extra_->tauJetBx.at(i)==0) {
			    h_TJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet));
				 h_TJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet));
				 h_TJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet));
			 }	 
			 myEvt_.Njet++;

		 }		 		 
		 h_TJet_Njet->Fill(l1extra_->nTauJets);
		// printf("Number of Jets in myEvt %i \n",myEvt_.Njet);

	   // Fill energy sums  (Are overflow flags accessible in l1extra?)
   	 for(unsigned int i=0; i< l1extra_->nMet; i++) {
  		     if(l1extra_->metBx.at(i)==0) {
			    myEvt_.ETT     = l1extra_->et.at(i) ; 
      	    myEvt_.ETM     = l1extra_->met.at(i) ; 
   	       myEvt_.PhiETM  = l1extra_->metPhi.at(i)  ;
				 
// Histogram Quantities
             h_Sum_ETT->Fill(myEvt_.ETT);
				 h_Sum_ETM->Fill(myEvt_.ETM);
				 h_Sum_PhiETM->Fill(myEvt_.PhiETM);				 
			  }
		 }	  	  		 	
   	 myEvt_.OvETT   = gt_ -> OvETT	;
		 myEvt_.OvETM   = gt_ -> OvETM	;  
   	  
   	 for(unsigned int i=0; i< l1extra_->nMht; i++) {
  		     if(l1extra_->mhtBx.at(i)==0) {   	 
   	       myEvt_.HTT     = l1extra_->ht.at(i) ; 
   	       myEvt_.HTM     = l1extra_->mht.at(i) ; 
   	       myEvt_.PhiHTM  = l1extra_->mhtPhi.at(i) ;
				 
// Histogram Quantities
             h_Sum_HTT->Fill(myEvt_.HTT);
				 h_Sum_HTM->Fill(myEvt_.HTM);
				 h_Sum_PhiHTM->Fill(myEvt_.PhiHTM);		
		     }
		 }	  
   	 myEvt_.OvHTM   = gt_ -> OvHTM	; 
	    myEvt_.OvHTT   = gt_ -> OvHTT	;
	
	
// Get the muon information  (FOR NOW TAKE THIS FROM THE GMT TO BE CONSISTENT WITH BRISTOL)
		 for(unsigned int i=0; i<gmt_->N; i++) {

      	 myEvt_.Bxmu.push_back(gmt_->CandBx[i]);
      	 myEvt_.Ptmu.push_back(gmt_->Pt[i]);
      	 myEvt_.Phimu.push_back(gmt_->Phi[i]); 
      	 myEvt_.Etamu.push_back(gmt_->Eta[i]); 
			 myEvt_.Qualmu.push_back(gmt_->Qual[i]);

// Histogram Quantities
			 
          if(myEvt_.Bxmu.at(myEvt_.Nmu)==0) {
			    h_Mu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu));
				 h_Mu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu));
				 h_Mu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu));
			 }	 
			 myEvt_.Nmu++; 
		 }	      
       h_Mu_Nmu->Fill(myEvt_.Nmu);	
		 
		 
// Extract the quantities from the GT
// ==================================		 
	 } else {	 
		 for(int i=0; i< gt_->Nele; i++) {

      	 myEvt_.Nele++;
      	 myEvt_.Bxel.push_back(gt_->Bxel[i]);
      	 myEvt_.Etel.push_back(gt_->Rankel[i]);
      	 myEvt_.Phiel.push_back(gt_->Phiel[i]);
      	 myEvt_.Etael.push_back(gt_->Etael[i]);
      	 myEvt_.Isoel.push_back(gt_->Isoel[i]);

		 }
		 // printf("Number of Electrons in myEvt %i \n",myEvt_.Nele);

		 for(int i=0; i< gt_->Njet; i++) {

      	 myEvt_.Njet++;
      	 myEvt_.Bxjet.push_back(gt_->Bxjet[i]);
      	 myEvt_.Etjet.push_back((gt_->Rankjet[i])*4.);
      	 myEvt_.Phijet.push_back(gt_->Phijet[i]);
      	 myEvt_.Etajet.push_back(gt_->Etajet[i]);
      	 myEvt_.Taujet.push_back(gt_->Taujet[i]);
      	 myEvt_.Fwdjet.push_back(gt_->Fwdjet[i]);

		 }
		// printf("Number of Jets in myEvt %i \n",myEvt_.Njet);

	   // Fill energy sums
   	 myEvt_.ETT     = (gt_ -> RankETT)/2. ; 
   	 myEvt_.OvETT   = gt_ -> OvETT	; 
   	 myEvt_.HTT     = (gt_ -> RankHTT)/2. ; 
   	 myEvt_.OvHTT   = gt_ -> OvHTT	; 
   	 myEvt_.ETM     = (gt_ -> RankETM)/2. ; 
   	 myEvt_.PhiETM  = gt_ -> PhiETM  ; 
   	 myEvt_.OvETM   = gt_ -> OvETM	; 
   	 myEvt_.HTM     = (gt_ -> RankHTM)*2. ; 
   	 myEvt_.PhiHTM  = gt_ -> PhiHTM  ; 
   	 myEvt_.OvHTM   = gt_ -> OvHTM	; 
   }
	 
	 return;
}

void L1Menu2015::InsertInMenu(std::string L1name, Bool_t value) {

	Bool_t post_prescale = false;

	Int_t prescale = 1;

	std::map<std::string, int>::const_iterator it = Prescales.find(L1name);
	if (it == Prescales.end() ) {
		std::cout << " --- NO PRESCALE DEFINED FOR " << L1name << " ---  SET P = 1 " << std::endl;
	}
	else {
		prescale = Prescales[L1name];
	}

	if (prescale >0) {
		Counts[L1name] ++;
		Int_t n = Counts[L1name];
		if ( n % prescale == 0) post_prescale = value; 
	}

	insert_names[insert_ibin] = L1name;
	insert_val[insert_ibin] = post_prescale ;

	insert_ibin ++;

}

Int_t L1Menu2015::L1BitNumber(std::string l1name) {

	std::map<std::string, int>::const_iterator it = BitMapping.find(l1name);
	if (it == BitMapping.end() ) {
		std::cout << " Wrong L1 name, not in BitMapping " << l1name << std::endl;
		return -1;
	}

	return BitMapping[l1name];
}

void L1Menu2015::FilL1Bits() {
        //printf("Run %i Event %i Alg Fired:",event_->run,event_->event);
	for (Int_t ibit=0; ibit < 128; ibit++) {
		PhysicsBits[ibit] = 0;
		if (ibit<64) {
			PhysicsBits[ibit] = (gt_->tw1[2]>>ibit)&1;
		}
		else {
			PhysicsBits[ibit] = (gt_->tw2[2]>>(ibit-64))&1;
		}
		//if(PhysicsBits[ibit]) printf(" %i ",ibit);
	}
	//printf("\n");
}       

void L1Menu2015::MyInit() {


	NBITS_TRIGS=0;

// ---- The bit std::mapping
//  Note: These are just for labeling purposes...they are not used for any array indexing in the code.
	
	BitMapping["L1_ZeroBias"] = 0 ;
	
	BitMapping["L1_SingleEG"] =  10;
	BitMapping["L1_SingleIsoEG"] =  11;
	BitMapping["L1_DoubleEG"] = 12 ;

	BitMapping["L1_SingleMu"] =  20;
	BitMapping["L1_DoubleMu"] =  21;
	
	BitMapping["L1_EG_Mu"] =  60;
	BitMapping["L1_Mu_EG"] =  61;

	BitMapping["L1_SingleTau"] = 30 ;
	BitMapping["L1_DoubleTau"] = 31 ;
	
	BitMapping["L1_SingleJet"] = 40 ;
	BitMapping["L1_DoubleJet"] = 41 ;
	BitMapping["L1_QuadJetC"] =   42;
	
	BitMapping["L1_ETM"] =  50 ;
	
	BitMapping["L1_HTT"] =  51 ;

	BitMapping["L1_SingleMu_ETM"] =  70 ;
	BitMapping["L1_SingleEG_ETM"] =  71 ;
	BitMapping["L1_SingleMu_CJet"] = 72 ;
	BitMapping["L1_SingleEG_CJet"] = 73 ;
	

//  DEFINE THE DEFAULT PARAMETERS (Can Override below)

// Define Trigger Parameters (Default 7E33 Menu)
	trigParList["L1_SingleEG"].primTh   = 22.;
	trigParList["L1_SingleEG"].secTh    = -1.;
	trigParList["L1_SingleEG"].triTh    = -1.;
	trigParList["L1_SingleEG"].quadTh   = -1.;
	trigParList["L1_SingleEG"].etaCut   = -1.;
	trigParList["L1_SingleEG"].minQual  = -1.;

	trigParList["L1_SingleIsoEG"].primTh   = 18.;
	trigParList["L1_SingleIsoEG"].secTh    = -1.;
	trigParList["L1_SingleIsoEG"].triTh    = -1.;
	trigParList["L1_SingleIsoEG"].quadTh   = -1.;
	trigParList["L1_SingleIsoEG"].etaCut   =  4.5; //corresponds to eta<2.17  
	trigParList["L1_SingleIsoEG"].minQual  = -1.;

	trigParList["L1_DoubleEG"].primTh   = 13.;
	trigParList["L1_DoubleEG"].secTh    =  7.;
	trigParList["L1_DoubleEG"].triTh    = -1.;
	trigParList["L1_DoubleEG"].quadTh   = -1.;
	trigParList["L1_DoubleEG"].etaCut   = -1.;
	trigParList["L1_DoubleEG"].minQual  = -1.;

	trigParList["L1_SingleEG_ETM"].primTh   = 20.;
	trigParList["L1_SingleEG_ETM"].secTh    = 20.;
	trigParList["L1_SingleEG_ETM"].triTh    = -1.;
	trigParList["L1_SingleEG_ETM"].quadTh   = -1.;
	trigParList["L1_SingleEG_ETM"].etaCut   = -1.;
	trigParList["L1_SingleEG_ETM"].minQual  = -1.;	  

	trigParList["L1_SingleEG_CJet"].primTh   = 20.;
	trigParList["L1_SingleEG_CJet"].secTh    = 32.;
	trigParList["L1_SingleEG_CJet"].triTh    = -1.;
	trigParList["L1_SingleEG_CJet"].quadTh   = -1.;
	trigParList["L1_SingleEG_CJet"].etaCut   = -1.;
	trigParList["L1_SingleEG_CJet"].minQual  = -1.;	  	 

	trigParList["L1_SingleMu"].primTh   = 16.;
	trigParList["L1_SingleMu"].secTh    = -1.;
	trigParList["L1_SingleMu"].triTh    = -1.;
	trigParList["L1_SingleMu"].quadTh   = -1.;
	trigParList["L1_SingleMu"].etaCut   =  2.1;
	trigParList["L1_SingleMu"].minQual  =  5;

	trigParList["L1_DoubleMu"].primTh   = 10.;
	trigParList["L1_DoubleMu"].secTh    =  0.;
	trigParList["L1_DoubleMu"].triTh    = -1.;
	trigParList["L1_DoubleMu"].quadTh   = -1.;
	trigParList["L1_DoubleMu"].etaCut   =  5.0;
	trigParList["L1_DoubleMu"].minQual  =  5;

	trigParList["L1_SingleMu_ETM"].primTh   = 20.;
	trigParList["L1_SingleMu_ETM"].secTh    = 20.;
	trigParList["L1_SingleMu_ETM"].triTh    = -1.;
	trigParList["L1_SingleMu_ETM"].quadTh   = -1.;
	trigParList["L1_SingleMu_ETM"].etaCut   =  5.0;
	trigParList["L1_SingleMu_ETM"].minQual  =  5;	  

	trigParList["L1_SingleMu_CJet"].primTh   = 20.;
	trigParList["L1_SingleMu_CJet"].secTh    = 32.;
	trigParList["L1_SingleMu_CJet"].triTh    = -1.;
	trigParList["L1_SingleMu_CJet"].quadTh   = -1.;
	trigParList["L1_SingleMu_CJet"].etaCut   =  5.0;
	trigParList["L1_SingleMu_CJet"].minQual  =  5;	

	trigParList["L1_EG_Mu"].primTh   = 3.5; //First threshold is on muon
	trigParList["L1_EG_Mu"].secTh    = 12.; //Second threshol is on EG
	trigParList["L1_EG_Mu"].triTh    = -1.;
	trigParList["L1_EG_Mu"].quadTh   = -1.;
	trigParList["L1_EG_Mu"].etaCut   =  5.0; //No meaning currently
	trigParList["L1_EG_Mu"].minQual  =  5;

	trigParList["L1_Mu_EG"].primTh   = 12.; //First threshold is on muon
	trigParList["L1_Mu_EG"].secTh    =  7.; //Second threshol is on EG
	trigParList["L1_Mu_EG"].triTh    = -1.;
	trigParList["L1_Mu_EG"].quadTh   = -1.;
	trigParList["L1_Mu_EG"].etaCut   =  5.0;
	trigParList["L1_Mu_EG"].minQual  =  5;

	trigParList["L1_SingleTau"].primTh   = 30.;
	trigParList["L1_SingleTau"].secTh    = -1.;
	trigParList["L1_SingleTau"].triTh    = -1.;
	trigParList["L1_SingleTau"].quadTh   = -1.;
	trigParList["L1_SingleTau"].etaCut   = -1.;
	trigParList["L1_SingleTau"].minQual  = -1.;

	trigParList["L1_DoubleTau"].primTh   = 44.;
	trigParList["L1_DoubleTau"].secTh    = 44.;
	trigParList["L1_DoubleTau"].triTh    = -1.;
	trigParList["L1_DoubleTau"].quadTh   = -1.;
	trigParList["L1_DoubleTau"].etaCut   =  4.5; //corresponds to eta<2.17 
	trigParList["L1_DoubleTau"].minQual  = -1.;

	trigParList["L1_SingleJet"].primTh   = 128.;
	trigParList["L1_SingleJet"].secTh    = -1.;
	trigParList["L1_SingleJet"].triTh    = -1.;
	trigParList["L1_SingleJet"].quadTh   = -1.;
	trigParList["L1_SingleJet"].etaCut   = -1.;
	trigParList["L1_SingleJet"].minQual  = -1.;

	trigParList["L1_DoubleJet"].primTh   = 56.;
	trigParList["L1_DoubleJet"].secTh    = 56.;
	trigParList["L1_DoubleJet"].triTh    = -1.;
	trigParList["L1_DoubleJet"].quadTh   = -1.;
	trigParList["L1_DoubleJet"].etaCut   = -1.;  //note this is implemented with Central Jets only
	trigParList["L1_DoubleJet"].minQual  = -1.;


	trigParList["L1_QuadJetC"].primTh   = 36.;
	trigParList["L1_QuadJetC"].secTh    = 36.;
	trigParList["L1_QuadJetC"].triTh    = 36.;
	trigParList["L1_QuadJetC"].quadTh   = 36.;
	trigParList["L1_QuadJetC"].etaCut   = -1.;
	trigParList["L1_QuadJetC"].minQual  = -1.;

	trigParList["L1_ETM"].primTh   = 36.;
	trigParList["L1_ETM"].secTh    = -1.;
	trigParList["L1_ETM"].triTh    = -1.;
	trigParList["L1_ETM"].quadTh   = -1.;
	trigParList["L1_ETM"].etaCut   = -1.;
	trigParList["L1_ETM"].minQual  = -1.;	

	trigParList["L1_HTT"].primTh   =150.;
	trigParList["L1_HTT"].secTh    = -1.;
	trigParList["L1_HTT"].triTh    = -1.;
	trigParList["L1_HTT"].quadTh   = -1.;
	trigParList["L1_HTT"].etaCut   = -1.;
	trigParList["L1_HTT"].minQual  = -1.;


// Set Prescales
	Prescales["L1_ZeroBias"] = 1 ;

	Prescales["L1_SingleEG"] =  1;
	Prescales["L1_SingleIsoEG"] =  1;
	Prescales["L1_DoubleEG"] = 1 ;

	Prescales["L1_SingleMu"] =  1;
	Prescales["L1_DoubleMu"] =  1;

	Prescales["L1_EG_Mu"] =  1;
	Prescales["L1_Mu_EG"] =  1;

	Prescales["L1_SingleTau"] = 1 ;
	Prescales["L1_DoubleTau"] = 1 ;

	Prescales["L1_SingleJet"] = 1 ;
	Prescales["L1_DoubleJet"] = 1 ;
	Prescales["L1_QuadJetC"] =   1;

	Prescales["L1_ETM"] =  1 ;

	Prescales["L1_HTT"] =  1 ;

	Prescales["L1_SingleMu_ETM"] =  1 ;
	Prescales["L1_SingleEG_ETM"] =  1 ;
	Prescales["L1_SingleMu_CJet"] = 1 ;
	Prescales["L1_SingleEG_CJet"] = 1 ;

	if (theL1Menu == 1) {  //Override default choices

// Define Trigger Parameters (Menu v1)
	  trigParList["L1_SingleEG"].primTh      = 39.;

	  trigParList["L1_SingleIsoEG"].primTh   = 29.;

	  trigParList["L1_DoubleEG"].primTh      = 20.;
	  trigParList["L1_DoubleEG"].secTh       =(7./13.)*trigParList["L1_DoubleEG"].primTh;
	  
	  trigParList["L1_SingleEG_ETM"].primTh  = 20.;
	  trigParList["L1_SingleEG_ETM"].secTh   = 20.;

	  trigParList["L1_SingleEG_CJet"].primTh = 20.;
	  trigParList["L1_SingleEG_CJet"].secTh  = 32.;

	  trigParList["L1_SingleMu"].primTh      = 35.;

	  trigParList["L1_DoubleMu"].primTh      = 26.;
	  trigParList["L1_DoubleMu"].secTh       =(0./10.)*trigParList["L1_DoubleMu"].primTh;
	  
	  trigParList["L1_SingleMu_ETM"].primTh  = 20.;
	  trigParList["L1_SingleMu_ETM"].secTh   = 20.;

	  trigParList["L1_SingleMu_CJet"].primTh = 20.;
	  trigParList["L1_SingleMu_CJet"].secTh  = 32.;

	  trigParList["L1_EG_Mu"].secTh          = 20.; //Second threshol is on EG
	  trigParList["L1_EG_Mu"].primTh         =(3.5/12.)*trigParList["L1_EG_Mu"].secTh; //First threshold is on muon

	  trigParList["L1_Mu_EG"].primTh         = 21.; //First threshold is on muon
	  trigParList["L1_Mu_EG"].secTh          =(7./12.)*trigParList["L1_Mu_EG"].primTh; //Second threshol is on EG

	  trigParList["L1_SingleTau"].primTh     = 30.;

	  trigParList["L1_DoubleTau"].primTh     = 56.;
	  trigParList["L1_DoubleTau"].secTh      =(44./44.)*trigParList["L1_DoubleTau"].primTh;
	  
	  trigParList["L1_SingleJet"].primTh     = 184.;

	  trigParList["L1_DoubleJet"].primTh     = 124.;
	  trigParList["L1_DoubleJet"].secTh      =(56./56.)*trigParList["L1_DoubleJet"].primTh;

	  trigParList["L1_QuadJetC"].primTh      = 96.;
	  trigParList["L1_QuadJetC"].secTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].triTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].quadTh      =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  
	  trigParList["L1_ETM"].primTh           = 81.;
	  trigParList["L1_HTT"].primTh           =436.;
	  
	} else if (theL1Menu == 2) {  //Override default choices// Thresholds defined by 8 TeV HPF 66 PU: Initial Choice (Scale 1.)

// Define Trigger Parameters (Menu v1)
	  trigParList["L1_SingleEG"].primTh      = 36.;

	  trigParList["L1_SingleIsoEG"].primTh   = 28.;

	  trigParList["L1_DoubleEG"].primTh      = 19.;
	  trigParList["L1_DoubleEG"].secTh       =(7./13.)*trigParList["L1_DoubleEG"].primTh;
	  
	  trigParList["L1_SingleEG_ETM"].primTh  = 20.;
	  trigParList["L1_SingleEG_ETM"].secTh   = 20.;

	  trigParList["L1_SingleEG_CJet"].primTh = 20.;
	  trigParList["L1_SingleEG_CJet"].secTh  = 32.;

	  trigParList["L1_SingleMu"].primTh      = 20.;

	  trigParList["L1_DoubleMu"].primTh      = 13.;
	  trigParList["L1_DoubleMu"].secTh       =(0./10.)*trigParList["L1_DoubleMu"].primTh;
	  
	  trigParList["L1_SingleMu_ETM"].primTh  = 20.;
	  trigParList["L1_SingleMu_ETM"].secTh   = 20.;

	  trigParList["L1_SingleMu_CJet"].primTh = 20.;
	  trigParList["L1_SingleMu_CJet"].secTh  = 32.;

	  trigParList["L1_EG_Mu"].secTh          = 16.; //Second threshol is on EG
	  trigParList["L1_EG_Mu"].primTh         =(3.5/12.)*trigParList["L1_EG_Mu"].secTh; //First threshold is on muon

	  trigParList["L1_Mu_EG"].primTh         = 15.; //First threshold is on muon
	  trigParList["L1_Mu_EG"].secTh          =(7./12.)*trigParList["L1_Mu_EG"].primTh; //Second threshol is on EG

	  trigParList["L1_SingleTau"].primTh     = 30.;

	  trigParList["L1_DoubleTau"].primTh     = 48.;
	  trigParList["L1_DoubleTau"].secTh      =(44./44.)*trigParList["L1_DoubleTau"].primTh;
	  
	  trigParList["L1_SingleJet"].primTh     = 192.;

	  trigParList["L1_DoubleJet"].primTh     = 88.;
	  trigParList["L1_DoubleJet"].secTh      =(56./56.)*trigParList["L1_DoubleJet"].primTh;

	  trigParList["L1_QuadJetC"].primTh      = 61.;
	  trigParList["L1_QuadJetC"].secTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].triTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].quadTh      =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  
	  trigParList["L1_ETM"].primTh           = 58.;
	  trigParList["L1_HTT"].primTh           =555.;

	} else if (theL1Menu == 3) {  //Override default choices// Thresholds defined by 8 TeV HPF 66 PU: Scaled Choice (Scale 1.75 )

// Define Trigger Parameters (Menu v1)
	  trigParList["L1_SingleEG"].primTh      = 29.;

	  trigParList["L1_SingleIsoEG"].primTh   = 23.;

	  trigParList["L1_DoubleEG"].primTh      = 17.;
	  trigParList["L1_DoubleEG"].secTh       =(7./13.)*trigParList["L1_DoubleEG"].primTh;

	  trigParList["L1_SingleMu"].primTh      = 15.;

	  trigParList["L1_DoubleMu"].primTh      = 10.;
	  trigParList["L1_DoubleMu"].secTh       =(0./10.)*trigParList["L1_DoubleMu"].primTh;

	  trigParList["L1_EG_Mu"].secTh          = 15.; //Second threshol is on EG
	  trigParList["L1_EG_Mu"].primTh         =(3.5/12.)*trigParList["L1_EG_Mu"].secTh; //First threshold is on muon

	  trigParList["L1_Mu_EG"].primTh         = 13.; //First threshold is on muon
	  trigParList["L1_Mu_EG"].secTh          =(7./12.)*trigParList["L1_Mu_EG"].primTh; //Second threshol is on EG

	  trigParList["L1_DoubleTau"].primTh     = 48.;
	  trigParList["L1_DoubleTau"].secTh      =(44./44.)*trigParList["L1_DoubleTau"].primTh;
	  
	  trigParList["L1_SingleJet"].primTh     = 164.;

	  trigParList["L1_DoubleJet"].primTh     = 80.;
	  trigParList["L1_DoubleJet"].secTh      =(56./56.)*trigParList["L1_DoubleJet"].primTh;

	  trigParList["L1_QuadJetC"].primTh      = 60.;
	  trigParList["L1_QuadJetC"].secTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].triTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].quadTh      =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  
	  trigParList["L1_ETM"].primTh           = 54.;
	  trigParList["L1_HTT"].primTh           =515.;


	} else if (theL1Menu == 4) {  //Override default choices// Thresholds defined by 8 TeV HPF 45 PU: Initial Choice (Scale 1.0)

// Define Trigger Parameters (Menu v1)
	  trigParList["L1_SingleEG"].primTh      = 34.;

	  trigParList["L1_SingleIsoEG"].primTh   = 27.;

	  trigParList["L1_DoubleEG"].primTh      = 19.;
	  trigParList["L1_DoubleEG"].secTh       =(7./13.)*trigParList["L1_DoubleEG"].primTh;
	  
	  trigParList["L1_SingleEG_ETM"].primTh  = 20.;
	  trigParList["L1_SingleEG_ETM"].secTh   = 20.;

	  trigParList["L1_SingleEG_CJet"].primTh = 20.;
	  trigParList["L1_SingleEG_CJet"].secTh  = 32.;

	  trigParList["L1_SingleMu"].primTh      = 21.;

	  trigParList["L1_DoubleMu"].primTh      = 13.;
	  trigParList["L1_DoubleMu"].secTh       =(0./10.)*trigParList["L1_DoubleMu"].primTh;
	  
	  trigParList["L1_SingleMu_ETM"].primTh  = 20.;
	  trigParList["L1_SingleMu_ETM"].secTh   = 20.;

	  trigParList["L1_SingleMu_CJet"].primTh = 20.;
	  trigParList["L1_SingleMu_CJet"].secTh  = 32.;

	  trigParList["L1_EG_Mu"].secTh          = 17.; //Second threshol is on EG
	  trigParList["L1_EG_Mu"].primTh         =(3.5/12.)*trigParList["L1_EG_Mu"].secTh; //First threshold is on muon

	  trigParList["L1_Mu_EG"].primTh         = 14.; //First threshold is on muon
	  trigParList["L1_Mu_EG"].secTh          =(7./12.)*trigParList["L1_Mu_EG"].primTh; //Second threshol is on EG

	  trigParList["L1_SingleTau"].primTh     = 30.;

	  trigParList["L1_DoubleTau"].primTh     = 52.;
	  trigParList["L1_DoubleTau"].secTh      =(44./44.)*trigParList["L1_DoubleTau"].primTh;
	  
	  trigParList["L1_SingleJet"].primTh     = 164.;

	  trigParList["L1_DoubleJet"].primTh     = 84.;
	  trigParList["L1_DoubleJet"].secTh      =(56./56.)*trigParList["L1_DoubleJet"].primTh;

	  trigParList["L1_QuadJetC"].primTh      = 48.;
	  trigParList["L1_QuadJetC"].secTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].triTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].quadTh      =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  
	  trigParList["L1_ETM"].primTh           = 51.;
	  trigParList["L1_HTT"].primTh           =335.;

	} else if (theL1Menu == 5) {  //Override default choices// Thresholds defined by 8 TeV HPF 45 PU: Scaled Choice (Scale 1.7) 

// Define Trigger Parameters (Menu v1)
	  trigParList["L1_SingleEG"].primTh      = 29.;

	  trigParList["L1_SingleIsoEG"].primTh   = 23.;

	  trigParList["L1_DoubleEG"].primTh      = 17.;
	  trigParList["L1_DoubleEG"].secTh       =(7./13.)*trigParList["L1_DoubleEG"].primTh;
	  
	  trigParList["L1_SingleMu"].primTh      = 17.;

	  trigParList["L1_DoubleMu"].primTh      = 10.;
	  trigParList["L1_DoubleMu"].secTh       =(0./10.)*trigParList["L1_DoubleMu"].primTh;

	  trigParList["L1_EG_Mu"].secTh          = 15.; //Second threshol is on EG
	  trigParList["L1_EG_Mu"].primTh         =(3.5/12.)*trigParList["L1_EG_Mu"].secTh; //First threshold is on muon

	  trigParList["L1_Mu_EG"].primTh         = 13.; //First threshold is on muon
	  trigParList["L1_Mu_EG"].secTh          =(7./12.)*trigParList["L1_Mu_EG"].primTh; //Second threshol is on EG

	  trigParList["L1_DoubleTau"].primTh     = 48.;
	  trigParList["L1_DoubleTau"].secTh      =(44./44.)*trigParList["L1_DoubleTau"].primTh;
	  
	  trigParList["L1_SingleJet"].primTh     = 152.;

	  trigParList["L1_DoubleJet"].primTh     = 76.;
	  trigParList["L1_DoubleJet"].secTh      =(56./56.)*trigParList["L1_DoubleJet"].primTh;

	  trigParList["L1_QuadJetC"].primTh      = 48.;
	  trigParList["L1_QuadJetC"].secTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].triTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].quadTh      =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  
	  trigParList["L1_ETM"].primTh           = 47.;
	  trigParList["L1_HTT"].primTh           =306.;	  	  

	} else if (theL1Menu == 6) {  //Menu v2 (Keep single lepton fixed) 66 PU Scale 0.75

// Define Trigger Parameters (Menu v2)
//	  trigParList["L1_SingleEG"].primTh      = 22.;

//	  trigParList["L1_SingleIsoEG"].primTh   = 18.;

	  trigParList["L1_DoubleEG"].primTh      = 20.;
	  trigParList["L1_DoubleEG"].secTh       =(7./13.)*trigParList["L1_DoubleEG"].primTh;

//	  trigParList["L1_SingleMu"].primTh      = 16.;

	  trigParList["L1_DoubleMu"].primTh      = 14.;
	  trigParList["L1_DoubleMu"].secTh       =(0./10.)*trigParList["L1_DoubleMu"].primTh;

	  trigParList["L1_EG_Mu"].secTh          = 17.; //Second threshol is on EG
	  trigParList["L1_EG_Mu"].primTh         =(3.5/12.)*trigParList["L1_EG_Mu"].secTh; //First threshold is on muon

	  trigParList["L1_Mu_EG"].primTh         = 16.; //First threshold is on muon
	  trigParList["L1_Mu_EG"].secTh          =(7./12.)*trigParList["L1_Mu_EG"].primTh; //Second threshol is on EG

	  trigParList["L1_DoubleTau"].primTh     = 52.;
	  trigParList["L1_DoubleTau"].secTh      =(44./44.)*trigParList["L1_DoubleTau"].primTh;
	  
	  trigParList["L1_SingleJet"].primTh     = 200.;

	  trigParList["L1_DoubleJet"].primTh     = 96.;
	  trigParList["L1_DoubleJet"].secTh      =(56./56.)*trigParList["L1_DoubleJet"].primTh;

	  trigParList["L1_QuadJetC"].primTh      = 64.;
	  trigParList["L1_QuadJetC"].secTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].triTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].quadTh      =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  
	  trigParList["L1_ETM"].primTh           = 61.;
	  trigParList["L1_HTT"].primTh           =567.;	  	  

	} else if (theL1Menu == 7) {  //Menu v2 (Keep single lepton fixed) 45 PU  Scale 0.85

// Define Trigger Parameters (Menu v2)
//	  trigParList["L1_SingleEG"].primTh      = 22.;

//	  trigParList["L1_SingleIsoEG"].primTh   = 18.;

	  trigParList["L1_DoubleEG"].primTh      = 20.;
	  trigParList["L1_DoubleEG"].secTh       =(7./13.)*trigParList["L1_DoubleEG"].primTh;

//	  trigParList["L1_SingleMu"].primTh      = 16.;

	  trigParList["L1_DoubleMu"].primTh      = 14.;
	  trigParList["L1_DoubleMu"].secTh       =(0./10.)*trigParList["L1_DoubleMu"].primTh;

	  trigParList["L1_EG_Mu"].secTh          = 17.; //Second threshol is on EG
	  trigParList["L1_EG_Mu"].primTh         =(3.5/12.)*trigParList["L1_EG_Mu"].secTh; //First threshold is on muon

	  trigParList["L1_Mu_EG"].primTh         = 15.; //First threshold is on muon
	  trigParList["L1_Mu_EG"].secTh          =(7./12.)*trigParList["L1_Mu_EG"].primTh; //Second threshol is on EG

	  trigParList["L1_DoubleTau"].primTh     = 56.;
	  trigParList["L1_DoubleTau"].secTh      =(44./44.)*trigParList["L1_DoubleTau"].primTh;
	  
	  trigParList["L1_SingleJet"].primTh     = 172.;

	  trigParList["L1_DoubleJet"].primTh     = 88.;
	  trigParList["L1_DoubleJet"].secTh      =(56./56.)*trigParList["L1_DoubleJet"].primTh;

	  trigParList["L1_QuadJetC"].primTh      = 48.;
	  trigParList["L1_QuadJetC"].secTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].triTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].quadTh      =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  
	  trigParList["L1_ETM"].primTh           = 53.;
	  trigParList["L1_HTT"].primTh           =343.;	  	  

	} else if (theL1Menu == 8) {  //Menu v3 (Keep double lepton fixed) 66 PU Scale 0.5

// Define Trigger Parameters (Menu v3)
	  trigParList["L1_SingleEG"].primTh      = 47.;

	  trigParList["L1_SingleIsoEG"].primTh   = 35.;

//	  trigParList["L1_DoubleEG"].primTh      = 13.;
//	  trigParList["L1_DoubleEG"].secTh       =(7./13.)*trigParList["L1_DoubleEG"].primTh;

	  trigParList["L1_SingleMu"].primTh      = 41.;

//	  trigParList["L1_DoubleMu"].primTh      = 10.;
//	  trigParList["L1_DoubleMu"].secTh       =(0./10.)*trigParList["L1_DoubleMu"].primTh;

//	  trigParList["L1_EG_Mu"].secTh          = 12.; //Second threshol is on EG
//	  trigParList["L1_EG_Mu"].primTh         =(3.5/12.)*trigParList["L1_EG_Mu"].secTh; //First threshold is on muon

//	  trigParList["L1_Mu_EG"].primTh         = 12.; //First threshold is on muon
//	  trigParList["L1_Mu_EG"].secTh          =(7./12.)*trigParList["L1_Mu_EG"].primTh; //Second threshol is on EG

	  trigParList["L1_DoubleTau"].primTh     = 52.;
	  trigParList["L1_DoubleTau"].secTh      =(44./44.)*trigParList["L1_DoubleTau"].primTh;
	  
	  trigParList["L1_SingleJet"].primTh     = 224.;

	  trigParList["L1_DoubleJet"].primTh     = 104.;
	  trigParList["L1_DoubleJet"].secTh      =(56./56.)*trigParList["L1_DoubleJet"].primTh;

	  trigParList["L1_QuadJetC"].primTh      = 64.;
	  trigParList["L1_QuadJetC"].secTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].triTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].quadTh      =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  
	  trigParList["L1_ETM"].primTh           = 64.;
	  trigParList["L1_HTT"].primTh           =586.;	  	  

	} else if (theL1Menu == 9) {  //Menu v3 (Keep double lepton fixed) 45 PU Scale 0.5

// Define Trigger Parameters (Menu v3)
	  trigParList["L1_SingleEG"].primTh      = 40.;

	  trigParList["L1_SingleIsoEG"].primTh   = 31.;

//	  trigParList["L1_DoubleEG"].primTh      = 13.;
//	  trigParList["L1_DoubleEG"].secTh       =(7./13.)*trigParList["L1_DoubleEG"].primTh;

	  trigParList["L1_SingleMu"].primTh      = 41.;

//	  trigParList["L1_DoubleMu"].primTh      = 10.;
//	  trigParList["L1_DoubleMu"].secTh       =(0./10.)*trigParList["L1_DoubleMu"].primTh;

//	  trigParList["L1_EG_Mu"].secTh          = 12.; //Second threshol is on EG
//	  trigParList["L1_EG_Mu"].primTh         =(3.5/12.)*trigParList["L1_EG_Mu"].secTh; //First threshold is on muon

//	  trigParList["L1_Mu_EG"].primTh         = 12.; //First threshold is on muon
//	  trigParList["L1_Mu_EG"].secTh          =(7./12.)*trigParList["L1_Mu_EG"].primTh; //Second threshol is on EG

	  trigParList["L1_DoubleTau"].primTh     = 56.;
	  trigParList["L1_DoubleTau"].secTh      =(44./44.)*trigParList["L1_DoubleTau"].primTh;
	  
	  trigParList["L1_SingleJet"].primTh     = 184.;

	  trigParList["L1_DoubleJet"].primTh     = 96.;
	  trigParList["L1_DoubleJet"].secTh      =(56./56.)*trigParList["L1_DoubleJet"].primTh;

	  trigParList["L1_QuadJetC"].primTh      = 52.;
	  trigParList["L1_QuadJetC"].secTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].triTh       =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  trigParList["L1_QuadJetC"].quadTh      =(36./36.)*trigParList["L1_QuadJetC"].primTh;
	  
	  trigParList["L1_ETM"].primTh           = 56.;
	  trigParList["L1_HTT"].primTh           =352.;	  	  
	}





/*
// -- test  to see where we stand if we would get rid of all p'ed seeds
// -- (in 2011 we spent ~ 20% of the rate in monitoring / control p'ed seeds..)

for (std::map<std::string, int>::iterator it=Prescales.begin(); it != Prescales.end(); it++) {
std::string name = it -> first;
Int_t p = it -> second;
if (p > 1 ) Prescales[name] = 0;
}
*/


for (std::map<std::string, int>::iterator it=Prescales.begin(); it != Prescales.end(); it++) {
	std::string name = it -> first;
	Counts[name] = 0;
	Biased[name] = false; 
}


// -- The "Biased" table is only used for the final print-out
// -- set true for seeds for which the rate estimation is biased by
// -- the sample (because of the seeds enabled in the high PU run)

// Biased["L1_TripleMu0"] = true;
// Biased["L1_DoubleMu_10_Open"] = true;
// Biased["L1_SingleEG5"] = true;
// Biased["L1_TripleEG7"] = true;
// Biased["L1_TripleEG_10_7_5"] = true;
// Biased["L1_SingleJet36"] = true;
// Biased["L1_DoubleJetC36"] = true;
// Biased["L1_DoubleJetC44_Eta1p74_WdEta4"] = true;
// Biased["L1_QuadJetC36"] = true;
// Biased["L1_QuadJetC40"] = true;
// Biased["L1_DoubleTauJet44er"] = true;
// Biased["L1_Mu3p5_DoubleEG5"] = true;
// Biased["L1_QuadJetC32"] = true;
// Biased["L1_TripleJet28_Central"] = true;

}

Bool_t L1Menu2015::SingleMuEta(Float_t ptcut, Float_t etaCut, Int_t qualmin) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) { 
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < qualmin) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        
		if (fabs(eta) > etaCut) continue;
		if (pt >= ptcut) muon = true;
	}

	Bool_t ok = muon;
	return ok;

}

Bool_t L1Menu2015::SingleMu(Float_t ptcut, Int_t qualmin) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;

	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		//    BX = 0, +/- 1 or +/- 2
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);		
		if ( qual < qualmin) continue;
		if (pt >= ptcut) muon = true;
	}

	Bool_t ok = muon;
	return ok;

}

Bool_t L1Menu2015::DoubleMuHighQEtaCut(Float_t ptcut, Float_t etacut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Int_t nmu=0;
	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		//    BX = 0, +/- 1 or +/- 2
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);		
		if ( qual < 4) continue;
		Float_t eta = myEvt_.Etamu.at(imu);		
		if (fabs(eta) > etacut) continue;
		if (pt >= ptcut) nmu ++;
	}

	Bool_t ok = (nmu >= 2 ) ;
	return ok;

}

Bool_t L1Menu2015::Onia(Float_t ptcut1, Float_t ptcut2, Float_t etacut, Int_t delta) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Int_t Nmu = myEvt_.Nmu;
	Int_t n1=0;
	Int_t n2=0;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);		
		if ( qual < 4) continue;
		Float_t eta = myEvt_.Etamu.at(imu);		
		if (fabs(eta) > etacut) continue;
		if (pt >= ptcut1) n1 ++;
		if (pt >= ptcut2) n2++;
	}

	Bool_t ok = (n1 >=1 && n2 >= 2 ) ;
	if (! ok) return false;

	// -- now the CORRELATION condition
	Bool_t CORREL = false;

	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        
		if (fabs(eta) > etacut) continue;
		if (pt < ptcut1) continue;
		Int_t ieta1 = etaMuIdx(eta);

		for (Int_t imu2=0; imu2 < Nmu; imu2++) {
			if (imu2 == imu) continue;
			Int_t bx2 = myEvt_.Bxmu.at(imu2);		
			if (bx2 != 0) continue;
			Float_t pt2 = myEvt_.Ptmu.at(imu2);			
			Int_t qual2 = myEvt_.Qualmu.at(imu2);        
			if ( qual2 < 4) continue;
			Float_t eta2 = myEvt_.Etamu.at(imu2);        
			if (fabs(eta2) > etacut) continue;
			if (pt2 < ptcut2) continue;
			Int_t ieta2 = etaMuIdx(eta2);

			Float_t deta = ieta1 - ieta2; 
		// std::cout << "eta 1 2 delta " << ieta1 << " " << ieta2 << " " << deta << std::endl;
			if ( fabs(deta) <= delta)  CORREL = true;
		// if (fabs ( eta - eta2) <=  1.7) CORREL = true; 
		}

	}

	return CORREL;

}

Bool_t L1Menu2015::DoubleMu(Float_t cut1, Float_t cut2, Int_t qualmin) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias
	if (! raw) return false;  

	Int_t n1=0;
	Int_t n2=0;
	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
	        if ( qual < qualmin) continue;
		if (qual < 4  && qual != 3 ) continue;
		if (pt >= cut1) n1 ++;
		if (pt >= cut2) n2 ++;
	}

	Bool_t ok = (n1 >= 1 && n2 >= 2 );
	return ok;

}

Bool_t L1Menu2015::TripleMu(Float_t cut1, Float_t cut2, Float_t cut3, Int_t qualmin) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t n3=0;
	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < qualmin) continue;
		if (pt >= cut1) n1 ++;
		if (pt >= cut2) n2 ++;
		if (pt >= cut3) n3 ++;
	}

	Bool_t ok = ( n1 >= 1 && n2 >= 2 && n3 >= 3 );
	return ok;

}

Bool_t L1Menu2015::DoubleMuXOpen(Float_t cut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( (qual >= 5 || qual == 3 ) && pt >= cut ) n1 ++;
		if ( pt >= 0 ) n2 ++;
	}

	Bool_t ok = ( n1 >= 1 && n2 >= 2 );
	return ok;
}


Bool_t L1Menu2015::EvalMenu(double lumiWeight) {

	insert_ibin = 0;
	
	if(Menu2015==0) {  //Menu v1
	

	   InsertInMenu("L1_SingleEG",      SingleEG_Eta(trigParList["L1_SingleEG"].primTh,trigParList["L1_SingleEG"].etaCut) );
	   InsertInMenu("L1_SingleIsoEG",   SingleIsoEG_Eta(trigParList["L1_SingleIsoEG"].primTh,trigParList["L1_SingleIsoEG"].etaCut) );
	   InsertInMenu("L1_DoubleEG",      DoubleEG(trigParList["L1_DoubleEG"].primTh,trigParList["L1_DoubleEG"].secTh) );
//	   InsertInMenu("L1_SingleEG_ETM",  EG_ETM(trigParList["L1_SingleEG_ETM"].primTh,trigParList["L1_SingleEG_ETM"].secTh) );
//	   InsertInMenu("L1_SingleEG_CJet", EG_JetCentral(trigParList["L1_SingleEG_CJet"].primTh,trigParList["L1_SingleEG_CJet"].secTh) );	   


	   InsertInMenu("L1_SingleMu",      SingleMuEta(trigParList["L1_SingleMu"].primTh,trigParList["L1_SingleMu"].etaCut,trigParList["L1_SingleMu"].minQual));
	   InsertInMenu("L1_DoubleMu",      DoubleMu(trigParList["L1_DoubleMu"].primTh,trigParList["L1_DoubleMu"].secTh,trigParList["L1_DoubleMu"].minQual));
//	   InsertInMenu("L1_SingleMu_ETM",  Muer_ETM(trigParList["L1_SingleMu_ETM"].primTh,trigParList["L1_SingleMu_ETM"].secTh,trigParList["L1_SingleMu_ETM"].etaCut, trigParList["L1_SingleMu_ETM"].minQual) );
//	   InsertInMenu("L1_SingleMu_CJet", Muer_JetCentral(trigParList["L1_SingleMu_CJet"].primTh,trigParList["L1_SingleMu_CJet"].secTh,trigParList["L1_SingleMu_CJet"].etaCut, trigParList["L1_SingleMu_CJet"].minQual) );

	   
	   InsertInMenu("L1_EG_Mu",         Mu_EG(trigParList["L1_EG_Mu"].primTh,trigParList["L1_EG_Mu"].secTh, trigParList["L1_EG_Mu"].minQual ) );
	   InsertInMenu("L1_Mu_EG",         Mu_EG(trigParList["L1_Mu_EG"].primTh,trigParList["L1_Mu_EG"].secTh, trigParList["L1_Mu_EG"].minQual) );
	   
	   InsertInMenu("L1_SingleJet",     SingleJet(trigParList["L1_SingleJet"].primTh) );
      InsertInMenu("L1_DoubleJet",     DoubleJetCentral(trigParList["L1_DoubleJet"].primTh,trigParList["L1_DoubleJet"].secTh) );
	   InsertInMenu("L1_QuadJetC",      QuadJetCentral(trigParList["L1_QuadJetC"].primTh,trigParList["L1_QuadJetC"].secTh,trigParList["L1_QuadJetC"].triTh,trigParList["L1_QuadJetC"].quadTh) );
//	   InsertInMenu("L1_SingleTau",     SingleTauJet(trigParList["L1_SingleTau"].primTh) );	   	   
	   InsertInMenu("L1_DoubleTau",     DoubleTauJetEta(trigParList["L1_DoubleTau"].primTh,trigParList["L1_DoubleTau"].secTh,trigParList["L1_DoubleTau"].etaCut) );

	   InsertInMenu("L1_ETM",           ETM(trigParList["L1_ETM"].primTh) );           
	   InsertInMenu("L1_HTT",           HTT(trigParList["L1_HTT"].primTh) );

	} else if(Menu2015==1){  //Menu v2 (Only single lepton triggers)


	   InsertInMenu("L1_SingleEG",      SingleEG_Eta(trigParList["L1_SingleEG"].primTh,trigParList["L1_SingleEG"].etaCut) );
	   InsertInMenu("L1_SingleIsoEG",   SingleIsoEG_Eta(trigParList["L1_SingleIsoEG"].primTh,trigParList["L1_SingleIsoEG"].etaCut) );
//	   InsertInMenu("L1_DoubleEG",      DoubleEG(trigParList["L1_DoubleEG"].primTh,trigParList["L1_DoubleEG"].secTh) );
//	   InsertInMenu("L1_SingleEG_ETM",  EG_ETM(trigParList["L1_SingleEG_ETM"].primTh,trigParList["L1_SingleEG_ETM"].secTh) );
//	   InsertInMenu("L1_SingleEG_CJet", EG_JetCentral(trigParList["L1_SingleEG_CJet"].primTh,trigParList["L1_SingleEG_CJet"].secTh) );	   


	   InsertInMenu("L1_SingleMu",      SingleMuEta(trigParList["L1_SingleMu"].primTh,trigParList["L1_SingleMu"].etaCut,trigParList["L1_SingleMu"].minQual));
//	   InsertInMenu("L1_DoubleMu",      DoubleMu(trigParList["L1_DoubleMu"].primTh,trigParList["L1_DoubleMu"].secTh,trigParList["L1_DoubleMu"].minQual));
//	   InsertInMenu("L1_SingleMu_ETM",  Muer_ETM(trigParList["L1_SingleMu_ETM"].primTh,trigParList["L1_SingleMu_ETM"].secTh,trigParList["L1_SingleMu_ETM"].etaCut, trigParList["L1_SingleMu_ETM"].minQual) );
//	   InsertInMenu("L1_SingleMu_CJet", Muer_JetCentral(trigParList["L1_SingleMu_CJet"].primTh,trigParList["L1_SingleMu_CJet"].secTh,trigParList["L1_SingleMu_CJet"].etaCut, trigParList["L1_SingleMu_CJet"].minQual) );

	   
//	   InsertInMenu("L1_EG_Mu",         Mu_EG(trigParList["L1_EG_Mu"].primTh,trigParList["L1_EG_Mu"].secTh, trigParList["L1_EG_Mu"].minQual ) );
//	   InsertInMenu("L1_Mu_EG",         Mu_EG(trigParList["L1_Mu_EG"].primTh,trigParList["L1_Mu_EG"].secTh, trigParList["L1_Mu_EG"].minQual) );
	   
//	   InsertInMenu("L1_SingleJet",     SingleJet(trigParList["L1_SingleJet"].primTh) );
//      InsertInMenu("L1_DoubleJet",     DoubleJetCentral(trigParList["L1_DoubleJet"].primTh,trigParList["L1_DoubleJet"].secTh) );
//	   InsertInMenu("L1_QuadJetC",      QuadJetCentral(trigParList["L1_QuadJetC"].primTh,trigParList["L1_QuadJetC"].secTh,trigParList["L1_QuadJetC"].triTh,trigParList["L1_QuadJetC"].quadTh) );
//	   InsertInMenu("L1_SingleTau",     SingleTauJet(trigParList["L1_SingleTau"].primTh) );	   	   
//	   InsertInMenu("L1_DoubleTau",     DoubleTauJetEta(trigParList["L1_DoubleTau"].primTh,trigParList["L1_DoubleTau"].secTh,trigParList["L1_DoubleTau"].etaCut) );

//	   InsertInMenu("L1_ETM",           ETM(trigParList["L1_ETM"].primTh) );           
//	   InsertInMenu("L1_HTT",           HTT(trigParList["L1_HTT"].primTh) );

	} else if(Menu2015==2){  //Menu v3 (Only dilepton lepton triggers)


//	   InsertInMenu("L1_SingleEG",      SingleEG_Eta(trigParList["L1_SingleEG"].primTh,trigParList["L1_SingleEG"].etaCut) );
//	   InsertInMenu("L1_SingleIsoEG",   SingleIsoEG_Eta(trigParList["L1_SingleIsoEG"].primTh,trigParList["L1_SingleIsoEG"].etaCut) );
	   InsertInMenu("L1_DoubleEG",      DoubleEG(trigParList["L1_DoubleEG"].primTh,trigParList["L1_DoubleEG"].secTh) );
//	   InsertInMenu("L1_SingleEG_ETM",  EG_ETM(trigParList["L1_SingleEG_ETM"].primTh,trigParList["L1_SingleEG_ETM"].secTh) );
//	   InsertInMenu("L1_SingleEG_CJet", EG_JetCentral(trigParList["L1_SingleEG_CJet"].primTh,trigParList["L1_SingleEG_CJet"].secTh) );	   


//	   InsertInMenu("L1_SingleMu",      SingleMuEta(trigParList["L1_SingleMu"].primTh,trigParList["L1_SingleMu"].etaCut,trigParList["L1_SingleMu"].minQual));
	   InsertInMenu("L1_DoubleMu",      DoubleMu(trigParList["L1_DoubleMu"].primTh,trigParList["L1_DoubleMu"].secTh,trigParList["L1_DoubleMu"].minQual));
//	   InsertInMenu("L1_SingleMu_ETM",  Muer_ETM(trigParList["L1_SingleMu_ETM"].primTh,trigParList["L1_SingleMu_ETM"].secTh,trigParList["L1_SingleMu_ETM"].etaCut, trigParList["L1_SingleMu_ETM"].minQual) );
//	   InsertInMenu("L1_SingleMu_CJet", Muer_JetCentral(trigParList["L1_SingleMu_CJet"].primTh,trigParList["L1_SingleMu_CJet"].secTh,trigParList["L1_SingleMu_CJet"].etaCut, trigParList["L1_SingleMu_CJet"].minQual) );

	   
	   InsertInMenu("L1_EG_Mu",         Mu_EG(trigParList["L1_EG_Mu"].primTh,trigParList["L1_EG_Mu"].secTh, trigParList["L1_EG_Mu"].minQual ) );
	   InsertInMenu("L1_Mu_EG",         Mu_EG(trigParList["L1_Mu_EG"].primTh,trigParList["L1_Mu_EG"].secTh, trigParList["L1_Mu_EG"].minQual) );
	   
//	   InsertInMenu("L1_SingleJet",     SingleJet(trigParList["L1_SingleJet"].primTh) );
//      InsertInMenu("L1_DoubleJet",     DoubleJetCentral(trigParList["L1_DoubleJet"].primTh,trigParList["L1_DoubleJet"].secTh) );
//	   InsertInMenu("L1_QuadJetC",      QuadJetCentral(trigParList["L1_QuadJetC"].primTh,trigParList["L1_QuadJetC"].secTh,trigParList["L1_QuadJetC"].triTh,trigParList["L1_QuadJetC"].quadTh) );
//	   InsertInMenu("L1_SingleTau",     SingleTauJet(trigParList["L1_SingleTau"].primTh) );	   	   
//	   InsertInMenu("L1_DoubleTau",     DoubleTauJetEta(trigParList["L1_DoubleTau"].primTh,trigParList["L1_DoubleTau"].secTh,trigParList["L1_DoubleTau"].etaCut) );

//	   InsertInMenu("L1_ETM",           ETM(trigParList["L1_ETM"].primTh) );           
//	   InsertInMenu("L1_HTT",           HTT(trigParList["L1_HTT"].primTh) );

   } else {
	
	   printf("No Menu2015 defined\n");
   }

	Int_t NN = insert_ibin;

	Int_t kOFFSET_old = kOFFSET;
	for (Int_t k=0; k < NN; k++) {
		TheTriggerBits[k + kOFFSET_old] = insert_val[k];
	}
	kOFFSET += insert_ibin;

	if (first) {

		NBITS_TRIGS = NN;

		for (Int_t ibin=0; ibin < insert_ibin; ibin++) {
			TString l1name = (TString)insert_names[ibin];
			h_Trig -> GetXaxis() -> SetBinLabel(ibin+1, l1name );
		}
		h_Trig -> GetXaxis() -> SetBinLabel(NN+1, "Triggered") ;

		for (Int_t k=1; k <= kOFFSET - kOFFSET_old ; k++) {
			h_All -> GetXaxis() -> SetBinLabel(k +kOFFSET_old , h_Trig -> GetXaxis() -> GetBinLabel(k) );
			h_Cumm -> GetXaxis() -> SetBinLabel(k +kOFFSET_old , h_Trig -> GetXaxis() -> GetBinLabel(k) );
		   h_Corr -> GetXaxis() -> SetBinLabel(k +kOFFSET_old , h_Trig -> GetXaxis() -> GetBinLabel(k) );
			h_Corr -> GetYaxis() -> SetBinLabel(k +kOFFSET_old , h_Trig -> GetXaxis() -> GetBinLabel(k) );
		}
	}

	Bool_t res = false;
	for (Int_t i=0; i < NN; i++) {
		res = res || insert_val[i] ;
		if (insert_val[i]) h_Trig -> Fill(i,lumiWeight);
	}
	if (res) h_Trig -> Fill(NN,lumiWeight);

	return res;
}

void L1Menu2015::EvalThresh(double lumiWeight) {
 
// Flag for whether do evaluate 2-D trigger space (This can be CPU intensive) 
   bool TwoDimScan = true;  

/* notes:

  1) Histograms are binned so that integer values are at the bin centers.
  2) For EG and Muon, the bin widths are 1 GeV (Not ideal for muons which have a non uniform distribution of values)
  3) For Jets, the bin widths are 4 GeV
  4) For EG and Jets we use the lower edge of the bin to define the threshold and avoid rounding errors
         e.g For EG:  Lower edge is 15.5 GeV so test will be whether EG > 15.5 but then this will be plotted
			             at bin center which is 16 GeV.  This avoids cases where numbers like  15.99999999999999   
							 would fail threshold or be plotted in the incorrect bin.
			For jets it is similar excepts bins step by 4 GeV.
  5) For Muons, we take the threholds from the bin centers and live with small rnding problems since muon pts are nonuniform.
  6) For non-primary thresholds, we scale by ratios. These ratio are multiplied by the bin_center values since they 
     correspond to the true threshold.
	      -> No attempt is currently made to adjust the non-primary thresholds so that they find actual integer values (or factors of 4 in the case of jets) 
  
  
*/
	//------- SingleMu -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleMu = h_SingleMu_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleMu+1; bin++){
		const float bin_low_edge = h_SingleMu_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleMu_byThreshold->GetBinCenter(bin);
		if(SingleMuEta(bin_center,trigParList["L1_SingleMu"].etaCut,trigParList["L1_SingleMu"].minQual)){
			h_SingleMu_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ DoubleMu ------------------------------------------------------------------------------------------------------
	const unsigned n_bins_DoubleMu = h_DoubleMu_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_DoubleMu+1; bin++){
		const float bin_low_edge = h_DoubleMu_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_DoubleMu_byThreshold->GetBinCenter(bin);
		if(DoubleMu(bin_center,bin_center*(trigParList["L1_DoubleMu"].secTh/trigParList["L1_DoubleMu"].primTh,trigParList["L1_DoubleMu"].minQual) )){
			h_DoubleMu_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ DoubleMu ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(TwoDimScan) {
	  const unsigned n_bins_DoubleMu_X = h2_DoubleMu_byThreshold->GetNbinsX();
	  const unsigned n_bins_DoubleMu_Y = h2_DoubleMu_byThreshold->GetNbinsY();
	  for(unsigned bin=1; bin <= n_bins_DoubleMu_X+1; bin++){	
	     const float bin_low_edge_X = h2_DoubleMu_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
		  const float bin_center_X = h2_DoubleMu_byThreshold->GetXaxis()->GetBinCenter(bin);
		  for(unsigned ybin=1; ybin <= bin; ybin++){	  //stop the y-axis scan at the diagonal since this is symmetric
		     const float bin_low_edge_Y = h2_DoubleMu_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		     const float bin_center_Y = h2_DoubleMu_byThreshold->GetYaxis()->GetBinCenter(ybin);
		     if(DoubleMu(bin_center_X, bin_center_Y,trigParList["L1_DoubleMu"].minQual)){
			     h2_DoubleMu_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		     }
		  }	
	  }
	}
	//----------------------------------------------------------------------------------------------------------------------



	//------- SingleMu_ETM -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleMu_ETM = h_SingleMu_ETM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleMu_ETM+1; bin++){
		const float bin_low_edge = h_SingleMu_ETM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleMu_ETM_byThreshold->GetBinCenter(bin);
		if(Muer_ETM(bin_center,bin_center*(trigParList["L1_SingleMu_ETM"].secTh/trigParList["L1_SingleMu_ETM"].primTh),trigParList["L1_SingleMu_ETM"].etaCut,trigParList["L1_SingleMu_ETM"].minQual)){
			h_SingleMu_ETM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ SingleMu_ETM ---- 2-D Evaluation ------------------------------------------------------------------------------
   if(TwoDimScan) {
		const unsigned n_bins_SingleMu_ETM_X = h2_SingleMu_ETM_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleMu_ETM_Y = h2_SingleMu_ETM_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleMu_ETM_X+1; bin++){	
	   	const float bin_low_edge_X = h2_SingleMu_ETM_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleMu_ETM_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleMu_ETM_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_SingleMu_ETM_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleMu_ETM_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(Muer_ETM(bin_center_X, bin_center_Y,trigParList["L1_SingleMu_ETM"].etaCut,trigParList["L1_SingleMu_ETM"].minQual)){
			   	h2_SingleMu_ETM_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------


	//------- SingleMu_CJet -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleMu_CJet = h_SingleMu_CJet_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleMu_CJet+1; bin++){
		const float bin_low_edge = h_SingleMu_CJet_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleMu_CJet_byThreshold->GetBinCenter(bin);
		if(Muer_JetCentral(bin_center,bin_center*(trigParList["L1_SingleMu_CJet"].secTh/trigParList["L1_SingleMu_CJet"].primTh),trigParList["L1_SingleMu_CJet"].etaCut,trigParList["L1_SingleMu_CJet"].minQual)){
			h_SingleMu_CJet_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ SingleMu_CJet ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(TwoDimScan) {
		const unsigned n_bins_SingleMu_CJet_X = h2_SingleMu_CJet_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleMu_CJet_Y = h2_SingleMu_CJet_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleMu_CJet_X+1; bin++){	
	   	const float bin_low_edge_X = h2_SingleMu_CJet_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleMu_CJet_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleMu_CJet_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_SingleMu_CJet_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleMu_CJet_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(Muer_JetCentral(bin_center_X, bin_center_Y,trigParList["L1_SingleMu_CJet"].etaCut,trigParList["L1_SingleMu_CJet"].minQual)){
			   	h2_SingleMu_CJet_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------


	//--------- SingleEG ---------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleEG = h_SingleEG_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleEG+1; bin++){
		const float bin_low_edge = h_SingleEG_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleEG_byThreshold->GetBinCenter(bin);
		if(SingleEG_Eta(bin_low_edge,trigParList["L1_SingleEG"].etaCut)){
			h_SingleEG_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------
	
	//--------- SingleIsoEG ------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleIsoEG = h_SingleIsoEG_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleIsoEG+1; bin++){
		const float bin_low_edge = h_SingleIsoEG_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleIsoEG_byThreshold->GetBinCenter(bin);
		if(SingleIsoEG_Eta(bin_low_edge,trigParList["L1_SingleIsoEG"].etaCut)){
			h_SingleIsoEG_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------	

	//------ DoubleEG ------------------------------------------------------------------------------------------------------
	const unsigned n_bins_DoubleEG = h_DoubleEG_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_DoubleEG+1; bin++){
		const float bin_low_edge = h_DoubleEG_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_DoubleEG_byThreshold->GetBinCenter(bin);
		if(DoubleEG(bin_low_edge, bin_center*(trigParList["L1_DoubleEG"].secTh/trigParList["L1_DoubleEG"].primTh))){
			h_DoubleEG_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ DoubleEG ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(TwoDimScan) {
		const unsigned n_bins_DoubleEG_X = h2_DoubleEG_byThreshold->GetNbinsX();
		const unsigned n_bins_DoubleEG_Y = h2_DoubleEG_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_DoubleEG_X+1; bin++){	
	   	const float bin_low_edge_X = h2_DoubleEG_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_DoubleEG_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= bin; ybin++){	         //stop the y-axis scan at the diagonal since this is symmetric
		   	const float bin_low_edge_Y = h2_DoubleEG_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_DoubleEG_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(DoubleEG(bin_low_edge_X, bin_low_edge_Y)){
			   	h2_DoubleEG_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------


	//------ SingleEG_ETM --------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleEG_ETM = h_SingleEG_ETM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleEG_ETM+1; bin++){
		const float bin_low_edge = h_SingleEG_ETM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleEG_ETM_byThreshold->GetBinCenter(bin);
		if(EG_ETM(bin_low_edge, bin_center*(trigParList["L1_SingleEG_ETM"].secTh/trigParList["L1_SingleEG_ETM"].primTh))){
			h_SingleEG_ETM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ SingleEG_ETM ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(TwoDimScan) {
		const unsigned n_bins_SingleEG_ETM_X = h2_SingleEG_ETM_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleEG_ETM_Y = h2_SingleEG_ETM_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleEG_ETM_X+1; bin++){	
	   	const float bin_low_edge_X = h2_SingleEG_ETM_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleEG_ETM_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleEG_ETM_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_SingleEG_ETM_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleEG_ETM_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(EG_ETM(bin_low_edge_X, bin_low_edge_Y)){
			   	h2_SingleEG_ETM_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------



	//------ SingleEG_CJet --------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleEG_CJet = h_SingleEG_CJet_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleEG_CJet+1; bin++){
		const float bin_low_edge = h_SingleEG_CJet_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleEG_CJet_byThreshold->GetBinCenter(bin);
		if(EG_JetCentral(bin_low_edge, bin_center*(trigParList["L1_SingleEG_CJet"].secTh/trigParList["L1_SingleEG_CJet"].primTh))){
			h_SingleEG_CJet_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------


	//------ SingleEG_CJet ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(TwoDimScan) {	
		const unsigned n_bins_SingleEG_CJet_X = h2_SingleEG_CJet_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleEG_CJet_Y = h2_SingleEG_CJet_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleEG_CJet_X+1; bin++){	
	   	const float bin_low_edge_X = h2_SingleEG_CJet_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleEG_CJet_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleEG_CJet_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_SingleEG_CJet_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleEG_CJet_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(EG_JetCentral(bin_low_edge_X, bin_low_edge_Y)){
			   	h2_SingleEG_CJet_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------


	//-------- Mu_EG -------------------------------------------------------------------------------------------------------
	const unsigned n_bins_Mu_EG = h_Mu_EG_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_Mu_EG+1; bin++){
		const float bin_low_edge = h_Mu_EG_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_Mu_EG_byThreshold->GetBinCenter(bin);
		if(Mu_EG(bin_center, bin_center*(trigParList["L1_Mu_EG"].secTh/trigParList["L1_Mu_EG"].primTh), trigParList["L1_Mu_EG"].minQual)){
			h_Mu_EG_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//------- EG_Mu --------------------------------------------------------------------------------------------------------
	const unsigned n_bins_EG_Mu = h_EG_Mu_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_EG_Mu+1; bin++){
		const float bin_low_edge = h_EG_Mu_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_EG_Mu_byThreshold->GetBinCenter(bin);
		if(Mu_EG( bin_center*(trigParList["L1_EG_Mu"].primTh/trigParList["L1_EG_Mu"].secTh), bin_low_edge, trigParList["L1_EG_Mu"].minQual) ){
			h_EG_Mu_byThreshold->Fill(bin_center,lumiWeight);
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------

	//------ Mu_EG ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(TwoDimScan) {	
		const unsigned n_bins_Mu_EG_X = h2_Mu_EG_byThreshold->GetNbinsX();
		const unsigned n_bins_Mu_EG_Y = h2_Mu_EG_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_Mu_EG_X+1; bin++){	
	   	const float bin_low_edge_X = h2_Mu_EG_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_Mu_EG_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_Mu_EG_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_Mu_EG_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_Mu_EG_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(Mu_EG(bin_center_X, bin_low_edge_Y, trigParList["L1_Mu_EG"].minQual)){
			   	h2_Mu_EG_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------



	//------ SingleJet -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleJet = h_SingleJet_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleJet+1; bin++){
		const float bin_low_edge = h_SingleJet_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleJet_byThreshold->GetBinCenter(bin);
		if(SingleJet(bin_low_edge)){
			h_SingleJet_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ DoubleJet -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_DoubleJet = h_DoubleJet_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_DoubleJet+1; bin++){
		const float bin_low_edge = h_DoubleJet_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_DoubleJet_byThreshold->GetBinCenter(bin);
		if(DoubleJetCentral(bin_low_edge,bin_center*(trigParList["L1_DoubleJet"].secTh/trigParList["L1_DoubleJet"].primTh))){
			h_DoubleJet_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------
	
	
	//------- QuadJet ------------------------------------------------------------------------------------------------------
	const unsigned n_bins_QuadJetCentral = h_QuadJetCentral_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_QuadJetCentral+1; bin++){
		const float bin_low_edge = h_QuadJetCentral_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_QuadJetCentral_byThreshold->GetBinCenter(bin);
		if(QuadJetCentral(bin_low_edge, bin_center*(trigParList["L1_QuadJetC"].secTh/trigParList["L1_QuadJetC"].primTh) , bin_center*(trigParList["L1_QuadJetC"].triTh/trigParList["L1_QuadJetC"].primTh), bin_center*(trigParList["L1_QuadJetC"].quadTh/trigParList["L1_QuadJetC"].primTh))){
			h_QuadJetCentral_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------
	
	//------ SingleTau -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleTau = h_SingleTau_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleTau+1; bin++){
		const float bin_low_edge = h_SingleTau_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleTau_byThreshold->GetBinCenter(bin);
		if(SingleTauJet(bin_low_edge)){
			h_SingleTau_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------			
	
	//------ DoubleTau -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_DoubleTau = h_DoubleTau_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_DoubleTau+1; bin++){
		const float bin_low_edge = h_DoubleTau_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_DoubleTau_byThreshold->GetBinCenter(bin);
		if(DoubleTauJetEta(bin_low_edge,bin_center*(trigParList["L1_DoubleTau"].secTh/trigParList["L1_DoubleTau"].primTh),trigParList["L1_DoubleTau"].etaCut)){
			h_DoubleTau_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------	

	//-------- HTT ---------------------------------------------------------------------------------------------------------
	const unsigned n_bins_HTT = h_HTT_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_HTT+1; bin++){
		const float bin_low_edge = h_HTT_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_HTT_byThreshold->GetBinCenter(bin);
		if(HTT(bin_low_edge)){
			h_HTT_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//------- ETM ----------------------------------------------------------------------------------------------------------
	const unsigned n_bins_ETM = h_ETM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_ETM+1; bin++){
		const float bin_low_edge = h_ETM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_ETM_byThreshold->GetBinCenter(bin);
		if(ETM(bin_low_edge)){
			h_ETM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}

}


Bool_t L1Menu2015::Mu_EG(Float_t mucut, Float_t EGcut , Int_t minMuQual) {

	Bool_t raw = PhysicsBits[0];    // ZeroBias
	if (! raw) return false;


	Bool_t eg =false;
	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {   
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < minMuQual) continue;
		if (pt >= mucut) muon = true;
	}

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Bool_t ok = muon && eg;
	return ok;

}

Bool_t L1Menu2015::DoubleMu_EG(Float_t mucut, Float_t EGcut ) {

	Bool_t raw = PhysicsBits[0]; 	// ZeroBias
	if (! raw) return false;

	Bool_t eg =false;
	Bool_t muon = false;
	Int_t  Nmuons = 0;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		// if ( qual < 4) continue;
		if (qual < 4 && qual !=3 ) continue;
		if (pt >= mucut) Nmuons ++;
	}
	if (Nmuons >= 2) muon = true;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Bool_t ok = muon && eg;
	return ok;

}

Bool_t L1Menu2015::Mu_DoubleEG(Float_t mucut, Float_t EGcut ) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias..
	if (! raw) return false;

	Bool_t eg =false;
	Bool_t muon = false;
	Int_t  Nmuons = 0;
	Int_t Nelectrons = 0;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if (pt >= mucut) Nmuons ++;
	}
	if (Nmuons >= 1) muon = true;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) Nelectrons ++;
	}  // end loop over EM objects
	if (Nelectrons >= 2) eg = true;

	Bool_t ok = muon && eg;
	return ok;

}

Bool_t L1Menu2015::MuOpen_EG(Float_t mucut, Float_t EGcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;


	Bool_t eg =false;
	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		if (pt >= mucut) muon = true;
	}

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Bool_t ok = muon && eg;
	return ok;

}

Bool_t L1Menu2015::Mu_JetCentral(Float_t mucut, Float_t jetcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t muon = false;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) jet = true;
	}

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if (pt >= mucut) muon = true;
	}

	Bool_t ok = muon && jet;
	return ok;

}

Bool_t L1Menu2015::Mu_DoubleJetCentral(Float_t mucut, Float_t jetcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t muon = false;

	Int_t n1 = 0;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) n1 ++;
	}
	jet = ( n1 >= 2 );

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if (pt >= mucut) muon = true;
	}

	Bool_t ok = muon && jet;
	return ok;

}

Bool_t L1Menu2015::Mu_JetCentral_LowerTauTh(Float_t mucut, Float_t jetcut, Float_t taucut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t central = false;
	Bool_t tau = false;
	Bool_t muon = false;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (! isTauJet) {  	// look at CentralJet
			if (pt >= jetcut) central = true;
		}
		else   {		// look at TauJets
			if (pt >= taucut) tau = true;
		}
	}
	jet = central || tau  ;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if (pt >= mucut) muon = true;
	}

	Bool_t ok = muon && jet;
	return ok;

}

Bool_t L1Menu2015::Muer_JetCentral(Float_t mucut, Float_t jetcut, Float_t etacut, Int_t minMuQual ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t muon = false;

   
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) jet = true;
	}

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < minMuQual) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        
		if (fabs(eta) > etacut) continue;

		if (pt >= mucut) muon = true;
	}

	Bool_t ok = muon && jet;
	return ok;

}

Bool_t L1Menu2015::Muer_JetCentral_LowerTauTh(Float_t mucut, Float_t jetcut, Float_t taucut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t central = false;
	Bool_t tau = false;
	Bool_t muon = false;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (! isTauJet) {       // look at CentralJet
			if (pt >= jetcut) central = true;
		}
		else   {                // look at TauJets
			if (pt >= taucut) tau = true;
		}
	}
	jet = central || tau  ;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        
		if (fabs(eta) > 2.1) continue;
		if (pt >= mucut) muon = true;
	}

	Bool_t ok = muon && jet;
	return ok;

}

Bool_t L1Menu2015::Mia(Float_t mucut, Float_t jet1, Float_t jet2) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t muon = false;
	Int_t n1 = 0;
	Int_t n2 = 0;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jet1) n1 ++;
		if (pt >= jet2) n2 ++;
	}       
	jet = (n1 >= 1 && n2 >= 2 );

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        
		if (fabs(eta) > 2.1) continue;        
		if (pt >= mucut) muon = true;
	} 

	Bool_t ok = muon && jet;
	if (! ok) return false;

	// now the CORREL condition


	Bool_t CORREL = false;

	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if (pt < mucut) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        
		if (fabs(eta) > 2.1) continue;

		Float_t phimu = myEvt_.Phimu.at(imu);
		Int_t iphi_mu = phiINjetCoord(phimu);
		Float_t etamu = myEvt_.Etamu.at(imu);
		Int_t ieta_mu = etaINjetCoord(etamu);

		for (Int_t ue=0; ue < Nj; ue++) {
			Int_t bxj = myEvt_.Bxjet[ue];        		
			if (bxj != 0) continue;
			Bool_t isFwdJet = myEvt_.Fwdjet[ue];
			if (isFwdJet) continue;
			Float_t rank = myEvt_.Etjet[ue];
			Float_t ptj = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
			if (ptj < jet2) continue;
			Float_t phijet = myEvt_.Phijet[ue];
			Int_t iphi_jet = (int)phijet;
			Float_t etajet = myEvt_.Etajet[ue];
			Int_t ieta_jet = (int)etajet;

			Bool_t corr_phi = correlateInPhi(iphi_jet, iphi_mu);
			Bool_t corr_eta = correlateInEta(ieta_jet, ieta_mu);
			Bool_t corr = corr_phi && corr_eta;
			if (corr) CORREL = true ;
		}
	}

	return CORREL;

}

Bool_t L1Menu2015::Mu_JetCentral_delta(Float_t mucut, Float_t jetcut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t muon = false;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) jet = true;
	}

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if (pt >= mucut) muon = true;
	}

	Bool_t ok = muon && jet;
	if (! ok) return false;

		//  -- now evaluate the delta condition :

	Bool_t CORREL = false;

	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if (pt < mucut) continue;

		Float_t phimu = myEvt_.Phimu.at(imu);
		Int_t iphi_mu = phiINjetCoord(phimu);
		Float_t etamu = myEvt_.Etamu.at(imu);
		Int_t ieta_mu = etaINjetCoord(etamu);

		for (Int_t ue=0; ue < Nj; ue++) {
			Int_t bxj = myEvt_.Bxjet[ue];        		
			if (bxj != 0) continue;
			Bool_t isFwdJet = myEvt_.Fwdjet[ue];
			if (isFwdJet) continue;
			Float_t rank = myEvt_.Etjet[ue];
			Float_t ptj = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
			if (ptj < jetcut) continue;
			Float_t phijet = myEvt_.Phijet[ue];
			Int_t iphi_jet = (int)phijet;
			Float_t etajet = myEvt_.Etajet[ue];
			Int_t ieta_jet = (int)etajet;

			Bool_t corr = correlateInPhi(iphi_jet, iphi_mu, 2) && correlateInEta(ieta_jet, ieta_mu, 2);
			if (corr) CORREL = true ;
		}
	}

	return CORREL;

}

Bool_t L1Menu2015::Mu_JetCentral_deltaOut(Float_t mucut, Float_t jetcut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t muon = false;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) jet = true;
	}

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if (pt >= mucut) muon = true;
	}

	Bool_t ok = muon && jet;
	if (! ok) return false;

		//  -- now evaluate the delta condition :

	Bool_t CORREL = false;

	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if (pt < mucut) continue;

		Float_t phimu = myEvt_.Phimu.at(imu);
		Int_t iphi_mu = phiINjetCoord(phimu);
//          Float_t etamu = myEvt_.Etamu.at(imu);
//          Int_t ieta_mu = etaINjetCoord(etamu);

// 		Int_t PhiOut[3];
// 		PhiOut[0] = iphi_mu;
// 		if (iphi_mu< 17) PhiOut[1] = iphi_mu+1;
// 		if (iphi_mu == 17) PhiOut[1] = 0;
// 		if (iphi_mu > 0) PhiOut[2] = iphi_mu - 1;
// 		if (iphi_mu == 0) PhiOut[2] = 17;

		for (Int_t ue=0; ue < Nj; ue++) {
			Int_t bxj = myEvt_.Bxjet[ue];        		
			if (bxj != 0) continue;
			Bool_t isFwdJet = myEvt_.Fwdjet[ue];
			if (isFwdJet) continue;
			Float_t rank = myEvt_.Etjet[ue];
			Float_t ptj = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
			if (ptj < jetcut) continue;
			Float_t phijet = myEvt_.Phijet[ue];
			Int_t iphi_jet = (int)phijet;
//                  Float_t etajet = myEvt_.Etajet[ue];
//                  Int_t ieta_jet = (int)etajet;

			if (! correlateInPhi(iphi_jet, iphi_mu, 8)) CORREL = true;


		}
	}

	return CORREL;

}

Bool_t L1Menu2015::Muer_TripleJetCentral(Float_t mucut, Float_t jet1, Float_t jet2, Float_t jet3)  {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t muon = false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t n3=0;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jet1) n1 ++;
		if (pt >= jet2) n2 ++;
		if (pt >= jet3) n3 ++;
	}

	jet = ( n1 >= 1 && n2 >= 2 && n3 >= 3 ) ;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		Float_t eta = myEvt_.Etamu.at(imu) ;
		if (fabs(eta) > 2.1 ) continue;
		if (pt >= mucut) muon = true;
	}

	Bool_t ok = muon && jet;
	return ok;

}

Bool_t L1Menu2015::Mu_HTT(Float_t mucut, Float_t HTcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t ht=false;
	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if (pt >= mucut) muon = true;
	}

	Float_t adc = myEvt_.HTT ;
	Float_t TheHTT = adc; // / 2. ;
	ht = (TheHTT >= HTcut) ;

	Bool_t ok = muon && ht;
	return ok;

}

Bool_t L1Menu2015::Muer_ETM(Float_t mucut, Float_t ETMcut, Float_t etacut, Int_t minMuQual ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t etm = false;
	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < minMuQual) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        
		if (fabs(eta) > etacut) continue;

		if (pt >= mucut) muon = true;
	}

	Float_t adc = myEvt_.ETM ;
	Float_t TheETM = adc; // / 2. ;
	etm = (TheETM >= ETMcut);

	Bool_t ok = muon && etm;
	return ok;

}

Bool_t L1Menu2015::EG_FwdJet(Float_t EGcut, Float_t FWcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t jet = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {        
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (!isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= FWcut) jet = true;
	}

	Bool_t ok = ( eg && jet);
	return ok;

}


Bool_t L1Menu2015::EG_JetCentral(Float_t EGcut, Float_t jetcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t eg = false;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) jet = true;
	}

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Bool_t ok = eg && jet;
	return ok;

}


Bool_t L1Menu2015::EG_DoubleJetCentral(Float_t EGcut, Float_t jetcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t jet = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Int_t Nj = myEvt_.Njet ;
	Int_t njets = 0;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) njets ++;
	}
	jet = ( njets >= 2 );

	Bool_t ok = ( eg && jet);
	return ok;

}

Bool_t L1Menu2015::EG_HT(Float_t EGcut, Float_t HTcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t ht = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Float_t adc = myEvt_.HTT ;
	Float_t TheHTT = adc; // / 2. ;
	ht = (TheHTT >= HTcut) ;

	Bool_t ok = ( eg && ht);
	return ok;

}


Bool_t L1Menu2015::EG_ETM(Float_t EGcut, Float_t ETMcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t etm = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Float_t adc = myEvt_.ETM ;;
	Float_t TheETM = adc; // / 2. ;
	etm = (TheETM >= ETMcut) ;

	Bool_t ok = ( eg && etm);
	
	//printf("EG_ETM:  eg %i  etm %i \n",eg,etm);
	
	return ok;

}


Bool_t L1Menu2015::DoubleEG_HT(Float_t EGcut, Float_t HTcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Int_t n1 = 0;
	Bool_t ht = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) n1 ++;
	}  // end loop over EM objects
	eg = ( n1 >= 2 );

	Float_t adc = myEvt_.HTT ;
	Float_t TheHTT = adc; // / 2. ;
	ht = (TheHTT >= HTcut) ;

	Bool_t ok = ( eg && ht);
	return ok;

}

Bool_t L1Menu2015::EGEta2p1_JetCentral(Float_t EGcut, Float_t jetcut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;
	
	Bool_t eg = false;
	Bool_t jet = false;

	Int_t Nele = myEvt_.Nele; 
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) jet = true;
	}

	Bool_t ok = (eg && jet);
	if (! ok) return false;


	//  -- now evaluate the delta condition :

	Bool_t CORREL = false;
	Int_t PhiOut[3];

	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt < EGcut) continue;

		Float_t phiel = myEvt_.Phiel[ue];
		Int_t iphiel = (int)phiel;

		PhiOut[0]=0; PhiOut[1]=0; PhiOut[2]=0;   

		PhiOut[0] = iphiel;
		if (iphiel< 17) PhiOut[1] = iphiel+1;
		if (iphiel == 17) PhiOut[1] = 0;
		if (iphiel > 0) PhiOut[2] = iphiel - 1;
		if (iphiel == 0) PhiOut[2] = 17;

		for (Int_t uj=0; uj < Nj; uj++) {
			Int_t bxj = myEvt_.Bxjet[uj];        		
			if (bxj != 0) continue;
			Bool_t isFwdJet = myEvt_.Fwdjet[uj];
			if (isFwdJet) continue;
			Float_t rankj = myEvt_.Etjet[uj];
			// Float_t ptj = rankj * 4;
			Float_t ptj = rankj; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[uj],rankj*4.,theL1JetCorrection);
			if (ptj < jetcut) continue;
			Float_t phijet = myEvt_.Phijet[uj];
			Int_t iphijet = (int)phijet; 

			if ( iphijet != PhiOut[0] && 
				iphijet != PhiOut[1] &&
				iphijet != PhiOut[2] ) CORREL = true;
		}  // loop over jets

	}  // end loop over EM objects

	return CORREL;
	
}

Bool_t L1Menu2015::EGEta2p1_JetCentral_LowTauTh(Float_t EGcut, Float_t jetcut, Float_t taucut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t jet = false;
	Bool_t central = false;
	Bool_t tau = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (! isTauJet) {
			if (pt >= jetcut) central = true;
		}
		else {
			if (pt >= taucut) tau = true;
		}
	}
	jet = tau || central;

	Bool_t ok = (eg && jet);
	if (! ok) return false;

	//  -- now evaluate the delta condition :

	Bool_t CORREL_CENTRAL = false;
	Bool_t CORREL_TAU = false;
	Int_t PhiOut[3];

	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt < EGcut) continue;

		Float_t phiel = myEvt_.Phiel[ue];
		Int_t iphiel = (int)phiel;

		PhiOut[0]=0; PhiOut[1]=0; PhiOut[2]=0;   

		PhiOut[0] = iphiel;
		if (iphiel< 17) PhiOut[1] = iphiel+1;
		if (iphiel == 17) PhiOut[1] = 0;
		if (iphiel > 0) PhiOut[2] = iphiel - 1;
		if (iphiel == 0) PhiOut[2] = 17;

		for (Int_t uj=0; uj < Nj; uj++) {
			Int_t bxj = myEvt_.Bxjet[uj];        		
			if (bxj != 0) continue;
			Bool_t isFwdJet = myEvt_.Fwdjet[uj];
			if (isFwdJet) continue;
			Bool_t isTauJet = myEvt_.Taujet[uj];
			Float_t rankj = myEvt_.Etjet[uj];
			Float_t ptj = rankj; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[uj],rankj*4.,theL1JetCorrection);
			Float_t phijet = myEvt_.Phijet[uj];
			Int_t iphijet = (int)phijet;

			if (! isTauJet) {

				if (ptj >= jetcut) { 
					if ( iphijet != PhiOut[0] &&
						iphijet != PhiOut[1] &&
						iphijet != PhiOut[2] ) CORREL_CENTRAL = true;
				}

			}
			else {
				if (ptj >= taucut) {
					if ( iphijet != PhiOut[0] &&
						iphijet != PhiOut[1] &&
						iphijet != PhiOut[2] ) CORREL_TAU = true;
				}

			}


		}  // loop over jets

	}  // end loop over EM objects

	Bool_t CORREL = CORREL_CENTRAL || CORREL_TAU ;
	return CORREL;

}

Bool_t L1Menu2015::IsoEGEta2p1_JetCentral_LowTauTh(Float_t EGcut, Float_t jetcut, Float_t taucut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t jet = false;
	Bool_t central = false;
	Bool_t tau = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Bool_t iso = myEvt_.Isoel[ue];
		if ( ! iso) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (! isTauJet) {
			if (pt >= jetcut) central = true;
		}
		else {
			if (pt >= taucut) tau = true;
		}
	}
	jet = tau || central;

	Bool_t ok = (eg && jet);
	if (! ok) return false;

		//  -- now evaluate the delta condition :

	Bool_t CORREL_CENTRAL = false;
	Bool_t CORREL_TAU = false;
	Int_t PhiOut[3];

	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Bool_t iso = myEvt_.Isoel[ue];
		if ( ! iso) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt < EGcut) continue;

		Float_t phiel = myEvt_.Phiel[ue];
		Int_t iphiel = (int)phiel;

		PhiOut[0]=0; PhiOut[1]=0; PhiOut[2]=0;   

		PhiOut[0] = iphiel; 
		if (iphiel< 17) PhiOut[1] = iphiel+1;
		if (iphiel == 17) PhiOut[1] = 0;
		if (iphiel > 0) PhiOut[2] = iphiel - 1;
		if (iphiel == 0) PhiOut[2] = 17;

		for (Int_t uj=0; uj < Nj; uj++) {
			Int_t bxj = myEvt_.Bxjet[uj];        		
			if (bxj != 0) continue;
			Bool_t isFwdJet = myEvt_.Fwdjet[uj];
			if (isFwdJet) continue;
			Bool_t isTauJet = myEvt_.Taujet[uj];
			Float_t rankj = myEvt_.Etjet[uj];
			Float_t ptj = rankj; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[uj],rankj*4.,theL1JetCorrection);
			Float_t phijet = myEvt_.Phijet[uj];
			Int_t iphijet = (int)phijet;

			if (! isTauJet) {

				if (ptj >= jetcut) {
					if ( iphijet != PhiOut[0] &&
						iphijet != PhiOut[1] &&
						iphijet != PhiOut[2] ) CORREL_CENTRAL = true;
				}

			}
			else {
				if (ptj >= taucut) {
					if ( iphijet != PhiOut[0] &&
						iphijet != PhiOut[1] &&
						iphijet != PhiOut[2] ) CORREL_TAU = true;
				}

			}


		}  // loop over jets

	}  // end loop over EM objects

	Bool_t CORREL = CORREL_CENTRAL || CORREL_TAU ;

	return CORREL;

}

Bool_t L1Menu2015::EGEta2p1_DoubleJetCentral(Float_t EGcut, Float_t jetcut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t jet = false;
	Int_t n2=0;

	Int_t Nele = myEvt_.Nele; 
	for (Int_t ue=0; ue < Nele; ue++) { 
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true; 
	}  // end loop over EM objects

	Int_t Nj = myEvt_.Njet ;               
	for (Int_t ue=0; ue < Nj; ue++) {      
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue; 
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) n2 ++;
	}

	jet = (n2 >= 2);

	Bool_t ok = (eg && jet);
	if (! ok) return false;

		//  -- now evaluate the delta condition :

	Bool_t CORREL = false;
	Int_t PhiOut[3];

	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt < EGcut) continue;

		Float_t phiel = myEvt_.Phiel[ue];
		Int_t iphiel = (int)phiel;

		PhiOut[0]=0; PhiOut[1]=0; PhiOut[2]=0; 

		PhiOut[0] = iphiel;
		if (iphiel< 17) PhiOut[1] = iphiel+1;
		if (iphiel == 17) PhiOut[1] = 0;
		if (iphiel > 0) PhiOut[2] = iphiel - 1;
		if (iphiel == 0) PhiOut[2] = 17;

		Int_t npair = 0;

		for (Int_t uj=0; uj < Nj; uj++) {
			Int_t bxj = myEvt_.Bxjet[uj];        		
			if (bxj != 0) continue; 
			Bool_t isFwdJet = myEvt_.Fwdjet[uj];
			if (isFwdJet) continue;
			Float_t rankj = myEvt_.Etjet[uj];
										// Float_t ptj = rankj * 4;
			Float_t ptj = rankj; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[uj],rankj*4.,theL1JetCorrection);
			if (ptj < jetcut) continue;
			Float_t phijet = myEvt_.Phijet[uj];
			Int_t iphijet = (int)phijet;

			if ( iphijet != PhiOut[0] &&
				iphijet != PhiOut[1] &&
				iphijet != PhiOut[2] ) npair ++;

		}  // loop over jets

		if (npair >= 2 ) CORREL = true ;

	}  // end loop over EM objects

	return CORREL;

}

Bool_t L1Menu2015::EGEta2p1_DoubleJetCentral_TripleJetCentral(Float_t EGcut, Float_t jetcut2, Float_t jetcut3) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t jet = false;
	Int_t n2=0;       
	Int_t n3=0;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue; 
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;  
	}  // end loop over EM objects

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut2) n2 ++;
		if (pt >= jetcut3) n3 ++;
	}

	jet = (n2 >= 2 && n3 >= 3 );

	Bool_t ok = (eg && jet);
	return ok;

}

Bool_t L1Menu2015::HTT_HTM(Float_t HTTcut, Float_t HTMcut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t htt = false;
	Bool_t htm = false;
	Float_t adc = myEvt_.HTT;   
	Float_t TheHTT =  adc; // / 2.   ;          
	htt = ( TheHTT >= HTTcut ) ;

	Int_t adc_HTM  = myEvt_.HTM ; 
	Float_t TheHTM = adc_HTM; // * 2.  ;           
	htm = ( TheHTM >= HTMcut );

	Bool_t ok = (htt && htm);
	return ok;

}

Bool_t L1Menu2015::JetCentral_ETM(Float_t jetcut, Float_t ETMcut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t etm = false;
	Bool_t jet = false;

	Float_t adc = myEvt_.ETM ;
	Float_t TheETM = adc; // / 2. ;
	etm = (TheETM >= ETMcut);

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) jet = true;
	}

	Bool_t ok = ( jet && etm );
	return ok;
	
}

Bool_t L1Menu2015::DoubleJetCentral_ETM(Float_t jetcut1, Float_t jetcut2, Float_t ETMcut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t etm = false; 
	Bool_t jet = false;
	Int_t n1=0;
	Int_t n2=0;

	Float_t adc = myEvt_.ETM ;
	Float_t TheETM = adc; // / 2. ;
	etm = (TheETM >= ETMcut);

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut1) n1 ++;
		if (pt >= jetcut2) n2 ++;
	}       
	jet = (n1 >= 1 && n2 >= 2);

	Bool_t ok = ( jet && etm );
	return ok;

}

Bool_t L1Menu2015::SingleJetCentral(Float_t cut ) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;

	Bool_t ok=false;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue; 
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= cut) ok = true;
	} 

	return ok;

}

Bool_t L1Menu2015::SingleJet(Float_t cut ) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;

	Bool_t ok=false;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= cut) ok = true;
	}

	return ok;

}

Bool_t L1Menu2015::DoubleJetCentral(Float_t cut1, Float_t cut2 ) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;


	Int_t n1=0;
	Int_t n2=0;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= cut1) n1++;
		if (pt >= cut2) n2++;
	}
	Bool_t ok = ( n1 >=1 && n2 >= 2);
	//if(ok) printf("Run %i  Event  %i \n",event_->run,event_->event);
	return ok;

}

Bool_t L1Menu2015::DoubleJet_Eta1p7_deltaEta4(Float_t cut1, Float_t cut2 ) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		Float_t eta = myEvt_.Etajet[ue];
		if (eta < 5.5 || eta > 15.5) continue;  // eta = 6 - 15
		if (pt >= cut1) n1++;
		if (pt >= cut2) n2++;
	}
	Bool_t ok = ( n1 >=1 && n2 >= 2);
	if (! ok) return false;

	// -- now the correlation

	Bool_t CORREL = false;

	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		Float_t eta1 = myEvt_.Etajet[ue];
		if (eta1 < 5.5 || eta1 > 15.5) continue;  // eta = 6 - 15
		if (pt < cut1) continue;

		for (Int_t ve=0; ve < Nj; ve++) {
			if (ve == ue) continue;
			Int_t bx2 = myEvt_.Bxjet[ve];        		
			if (bx2 != 0) continue;
			Bool_t isFwdJet2 = myEvt_.Fwdjet[ve];
			if (isFwdJet2) continue;
			Float_t rank2 = myEvt_.Etjet[ve];
			Float_t pt2 = rank2 * 4;
			Float_t eta2 = myEvt_.Etajet[ve];
			if (eta2 < 5.5 || eta2 > 15.5) continue;  // eta = 6 - 15
			if (pt2 < cut2) continue;

			Bool_t corr = correlateInEta((int)eta1, (int)eta2, 4);
			if (corr) CORREL = true;
		}


	}

	return CORREL ;

}


Bool_t L1Menu2015::SingleTauJet(Float_t cut, Float_t etaCut ) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;

	Int_t n1=0;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue; 
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (! isTauJet) continue;
		Float_t rank = myEvt_.Etjet[ue];    // the rank of the electron
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		Float_t eta = myEvt_.Etajet[ue];
      if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16  // eta = 5 - 16
		if (pt >= cut) n1++;
	}  // end loop over jets

	Bool_t ok = ( n1 >=1 );
	return ok;

}

Bool_t L1Menu2015::DoubleTauJetEta(Float_t cut1, Float_t cut2, Float_t etaCut) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue; 
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (! isTauJet) continue;
		Float_t rank = myEvt_.Etjet[ue];    // the rank of the electron
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		Float_t eta = myEvt_.Etajet[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		if (pt >= cut1) n1++;
		if (pt >= cut2) n2++;
	}  // end loop over jets

	Bool_t ok = ( n1 >=1 && n2 >= 2);
	return ok;

}

Bool_t L1Menu2015::TripleJetCentral(Float_t cut1, Float_t cut2, Float_t cut3 ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t n3=0;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= cut1) n1++;
		if (pt >= cut2) n2++;
		if (pt >= cut3) n3++;
	}

	Bool_t ok = ( n1 >=1 && n2 >= 2 && n3 >= 3 );
	return ok;

}

Bool_t L1Menu2015::TripleJet_VBF(Float_t jet1, Float_t jet2, Float_t jet3 ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;        
	Bool_t jetf1=false;           
	Bool_t jetf2=false;   

	Int_t n1=0;
	Int_t n2=0;
	Int_t n3=0;

	Int_t f1=0;
	Int_t f2=0;
	Int_t f3=0;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);

		if (isFwdJet) {
			if (pt >= jet1) f1 ++;
			if (pt >= jet2) f2 ++;
			if (pt >= jet3) f3 ++;              
		} 
		else {
			if (pt >= jet1) n1 ++;
			if (pt >= jet2) n2 ++;
			if (pt >= jet3) n3 ++;
		}    
	}

	jet   = ( n1 >= 1 && n2 >= 2 && n3 >= 3 ) ;        
	jetf1 = ( f1 >= 1 && n2 >= 1 && n3 >= 2 ) ;  // numbers change ofcourse    
	jetf2 = ( n1 >= 1 && f2 >= 1 && n3 >= 2 ) ;  

	Bool_t ok = false;

	if( jet || jetf1 || jetf2 ) ok =true;

	return ok;
}

Bool_t L1Menu2015::QuadJetCentral(Float_t cut1, Float_t cut2, Float_t cut3, Float_t cut4 ) {

// cut1 >= cut2  >= cut3 >= cut4

	// ZeroBias
// Bool_t raw = PhysicsBits[16];  // SingleJet36
// if (! raw) return false;
	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;


	Int_t n1=0;
	Int_t n2=0;
	Int_t n3=0;
	Int_t n4=0;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= cut1) n1++;
		if (pt >= cut2) n2++;
		if (pt >= cut3) n3++;
		if (pt >= cut4) n4++;
	}

	Bool_t ok = ( n1 >=1 && n2 >= 2 && n3 >= 3 && n4 >= 4);
	return ok;

}


Bool_t L1Menu2015::ETM(Float_t ETMcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Float_t adc = myEvt_.ETM ;
	Float_t TheETM = adc; // / 2. ;

	if (TheETM < ETMcut) return false;
	return true;

}

Bool_t L1Menu2015::HTT(Float_t HTTcut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Float_t adc = myEvt_.HTT ;
	Float_t TheHTT = adc; // / 2. ;

	if (TheHTT < HTTcut) return false;
	return true;

}

Bool_t L1Menu2015::ETT(Float_t ETTcut) {

	Float_t adc = myEvt_.ETT ;
	Float_t TheETT = adc; // / 2. ;

	if (TheETT < ETTcut) return false;

	return true;

}


Bool_t L1Menu2015::SingleEG(Float_t cut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t ok=false; 
	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {               
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ; 
		if (pt >= cut) ok = true;
	}  // end loop over EM objects

	return ok; 

}

Bool_t L1Menu2015::SingleIsoEG_Eta(Float_t cut, Float_t etaCut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t ok=false;
	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Bool_t iso = myEvt_.Isoel[ue];
		if (! iso) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= cut) ok = true;
	}  // end loop over EM objects

	return ok;

}

Bool_t L1Menu2015::SingleEG_Eta(Float_t cut, Float_t etaCut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t ok=false;
	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= cut) ok = true;
	}  // end loop over EM objects

	return ok;

}

Bool_t L1Menu2015::DoubleEG(Float_t cut1, Float_t cut2 ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {               
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= cut1) n1++;
		if (pt >= cut2) n2++;
	}  // end loop over EM objects

	Bool_t ok = ( n1 >= 1 && n2 >= 2) ;
	//if(ok) printf("Found doubleEG event Run %i Event %i \n",event_->run,event_->event);
	return ok;

}

Bool_t L1Menu2015::TripleEG(Float_t cut1, Float_t cut2, Float_t cut3 ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t n3=0;
	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= cut1) n1++;
		if (pt >= cut2) n2++;
		if (pt >= cut3) n3++;
	}  // end loop over EM objects

	Bool_t ok = ( n1 >= 1 && n2 >= 2 && n3 >= 3) ;
	return ok;


}


void L1Menu2015::Loop(Bool_t calcThreshold, Bool_t useL1Extra, TString lsRunFile, const int n_events_) {

	
	const Int_t nevents = (n_events_ < 0) ? GetEntries() : n_events_;

	Int_t NPASS = 0; 

//	Int_t nPAG =0;
	first = true;
	int cnt=0;
   int lastLumi = -1;

// Load luminosity information
   int ind = 0;
   double LS[450], IntL[450], PU[450], InstL[450];
	//TString lsRunFile = "getLumi_out_pixCorrLumi_66PU_stdCorr.txt";
   ifstream ifs( lsRunFile );
	while(ifs){
		ifs >> LS[ind];
		ifs >> IntL[ind];
		ifs >> InstL[ind];
		ifs >> PU[ind];
		printf("LS %f  InstL %6.3e \n",LS[ind],InstL[ind]);
		ind++;
	}




	for (Long64_t i=0; i<nevents; i++)
	{     
	//load the i-th event
		Long64_t ientry = LoadTree(i); if (ientry < 0) break;
		GetEntry(i);

      cnt++;
      if(cnt%(int)pow(10.,(double)((int)log10((double)cnt)))==0) printf("Event Number %i\n",cnt);


// Fill my event data
      fillDataStructure(useL1Extra);


//      Fill the physics bits:
                //printf("Entry %i",i);
		FilL1Bits();

		if (first) MyInit();

		//HFW
		//PhysicsBits[0] = true; //force it!
		Bool_t raw = PhysicsBits[0];  // ZeroBias
		if (! raw) continue;


//  --- Reset the emulated "Trigger Bits"
		kOFFSET = 0;
		for (Int_t k=0; k < N128; k++) {
			TheTriggerBits[k] = false;
		}


// Get the instantaneous luminosity for this lumiSection (defaults at set in RunL1_HFW)
			if (event_->run == 198588) theZeroBiasPrescale = 44.;
			if (event_->run == 198603) theZeroBiasPrescale = 92.;
			if (event_->run == 198609) theZeroBiasPrescale = 92.;
		   int thisLumi = event_->lumi; //get this event's lumi
			if(thisLumi != lastLumi) {
			  // Note: If LS not found in list, value of inst. lumi will not be changed from default set in RunL1_HFW.
 			  for(int k=0; k<ind-1; k++) if( LS[k]==thisLumi ) theLumiForThisSetOfLumiSections = InstL[k]/1.0e32; 		   
			  lastLumi = thisLumi;
			}         

/*  Determine the Weight for this event to turn it into a contribution to a rate

  Data:
        scal = theZeroBiasPrescale *                                  <-- All Prescales (L1 and any HLT if imposed)
		         (theTargetLumi/theLumiForThisSetOfLumiSections) /      <-- Ratio of luminosity
					(23.3570304 * theNumberOfUserdLumiSections) /          <-- total time
					 1000.                                                 <-- turn Hz to kHz 

  MC:
*/
// For test L1 Upgrade Triggers all are collected in one place.
      double lumiWeight = theZeroBiasPrescale *                          
                          (theTargetLumi/theLumiForThisSetOfLumiSections) /  
                          (23.3570304 * theNumberOfUserdLumiSections) / 1000.;
		Bool_t pass = EvalMenu(lumiWeight);
		if(calcThreshold) EvalThresh(lumiWeight);

		if (pass) NPASS ++;

// We have done the first event...used elsewhere
		first = false;

	// -- now the pure rate stuff
	// -- kOFFSET now contains the number of triggers we have calculated
      Bool_t anyTrigger = false;
		Bool_t firstTrig  = true;
		
		for (Int_t k=0; k < kOFFSET; k++) {
			if ( ! TheTriggerBits[k] ) continue;
		
		   anyTrigger = true;
			h_All -> Fill(k,lumiWeight);

		// did the event pass another trigger ?
			Bool_t nonpure = false;
			for (Int_t k2=0; k2 < kOFFSET; k2++) {
				if (k2 == k) {
				  h_Corr->Fill(k,k,lumiWeight);
				  //continue;
				} else if ( TheTriggerBits[k2] ) {
				  nonpure = true;
				  h_Corr->Fill(k,k2,lumiWeight);
				}	
			}
			Bool_t pure = !nonpure ;
			if (pure) {
			  h_Pure -> Fill(k,lumiWeight);
			  h_Pure -> Fill(kOFFSET,lumiWeight);
			} 
			
// Calculate the cummulative rates
         if(firstTrig) {
			   for(Int_t k2=k; k2<kOFFSET; k2++) h_Cumm->Fill(k2,lumiWeight);
				firstTrig = false;
			}
			
			
		}
		if(anyTrigger) h_All->Fill(kOFFSET,lumiWeight);
		h_All -> GetXaxis() -> SetBinLabel(kOFFSET+1,"Triggered Evt");
		h_Pure -> GetXaxis() -> SetBinLabel(kOFFSET+1,"Triggered Evt");

	}  // end evt loop


	std::cout << " Prescales for: " << theTargetLumi << ", LumiForThisSetOfLumiSections = " << theLumiForThisSetOfLumiSections << ", L1NtupleFileName = " << theL1NtupleFileName << std::endl;
	std::cout << std::endl << " --------------------------------------------------------- " << std::endl << std::endl;

	
// Loop over triggers and print out configuration
//=================================================
        for(std::map<std::string, trigPar>::iterator itr = trigParList.begin(); itr != trigParList.end(); itr++) {
	   std::cout  << setw(20) << itr->first << "   Parameters(" << setw(5) << (itr->second).primTh << "," 
	                                                << setw(5) << (itr->second).secTh << ","
						        << setw(5) << (itr->second).triTh << ","
						        << setw(5) << (itr->second).quadTh << ","
						        << setw(5) << (itr->second).etaCut << ","
						        << setw(5) << (itr->second).minQual << ")"  << std::endl;
	} 
	std::cout << std::endl;
        	
}



void RunL1_HFW(Bool_t calcThreshold=true,Bool_t useL1Extra=true,Int_t usedL1Menu=0,Float_t targetlumi=200,Int_t whichFileAndLumiToUse=1,int which_jet_seed_to_use=0, Int_t pMenu2015 = 0, Int_t pNevts = -1) {

// Make sure we get the errors correct on histograms   
	TH1::SetDefaultSumw2();


//which_jet_seed_to_use:

//0: default
//1: reEmul
//2: 5 GeV 


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Using the 10-bunches High PU run, 179828. In this run ETM30 and HTT50 were enabled (they were not enabled in the 1-bunch run).
	Float_t NumberOfUserdLumiSections=0; 
	Float_t LumiForThisSetOfLumiSections=0;
	std::string L1NtupleFileName="";
	std::string jobTag="";
	Float_t AveragePU=0;
	Float_t ZeroBiasPrescale=0;
	Bool_t L1JetCorrection=false;
	Menu2015 = pMenu2015;
	Int_t procNevts = pNevts;
	TString lsFileName = "";
	

	if (whichFileAndLumiToUse==1) {
	// -- Run 179828, LS 374 - 394, PU=28:
		NumberOfUserdLumiSections = 21; 
		LumiForThisSetOfLumiSections = 0.437;
		switch(which_jet_seed_to_use){
			case 0:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1Tree_R179828_LS374-394.root";
				jobTag = "R179828_LS374-394";
			break;
			case 1:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1SKIM_R179828_LS374-394_NoGCT_reEmul.root";
				jobTag = "R179828_LS374-394_NoGCT_reEmul";
			break;
			case 2:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1SKIM_R179828_LS374-394_GCT_jet_seed_5GeV.root";
				jobTag = "R179828_LS374-394_GCT_jet_seed_5GeV";
			break;
			default: cout << __LINE__ << ":" << __FILE__ << endl; 
				
		}
		AveragePU = 28;
		ZeroBiasPrescale = 3 * 148; //L1 Prescale * HLT Prescale
		L1JetCorrection=false;
	}
	else if (whichFileAndLumiToUse==2) {

	// -- Run 179828, LS 300 - 320, PU=29:
		NumberOfUserdLumiSections = 21; 
		LumiForThisSetOfLumiSections = 0.463;
		switch(which_jet_seed_to_use){
			case 0:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1Tree_R179828_LS300-320.root";
				jobTag = "R179828_LS300-320";
			break;
			case 1:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1SKIM_R179828_LS300-320_NoGCT_reEmul.root";
				jobTag = "R179828_LS300-320_NoGCT_reEmul";
			break;
			case 2:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1SKIM_R179828_LS300-320_GCT_jet_seed_5GeV.root";
				jobTag = "R179828_LS300-320_GCT_jet_seed_5GeV";
			break;
			default: cout << __LINE__ << ":" << __FILE__ << endl; 
				
		}
		AveragePU = 29;
		ZeroBiasPrescale = 3 * 148; //L1 Prescale * HLT Prescale
		L1JetCorrection=false;
	}
	else if (whichFileAndLumiToUse==3) {

	// -- Run 179828, LS 270 - 290, PU=30:
		NumberOfUserdLumiSections = 21; 
		LumiForThisSetOfLumiSections = 0.470;
		switch(which_jet_seed_to_use){
			case 0:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1Tree_R179828_LS270-290.root";
				jobTag = "R179828_LS270-290";
			break;
			case 1:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1SKIM_R179828_LS270-290_NoGCT_reEmul.root";
				jobTag = "R179828_LS270-290_NoGCT_reEmul";
			break;
			case 2:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1SKIM_R179828_LS270-290_GCT_jet_seed_5GeV.root";
				jobTag = "R179828_LS270-290_GCT_jet_seed_5GeV";
			break;
			default: cout << __LINE__ << ":" << __FILE__ << endl; 
				
		}
		AveragePU = 29;
		ZeroBiasPrescale = 3 * 148; //L1 Prescale * HLT Prescale
		L1JetCorrection=false;
	}
	else if (whichFileAndLumiToUse==4) {

	// -- Run 179828, LS 140 - 160, PU=33:
		NumberOfUserdLumiSections = 21; 
		LumiForThisSetOfLumiSections = 0.509;
		switch(which_jet_seed_to_use){
			case 0:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1Tree_R179828_LS140-160.root";
				jobTag = "R179828_LS140-160";
			break;
			case 1:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1SKIM_R179828_LS140-160_NoGCT_reEmul.root";
				jobTag = "R179828_LS140-160_NoGCT_reEmul";
			break;
			case 2:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1SKIM_R179828_LS140-160_GCT_jet_seed_5GeV.root";
				jobTag = "R179828_LS140-160_GCT_jet_seed_5GeV";
			break;
			default: cout << __LINE__ << ":" << __FILE__ << endl; 
				
		}
		AveragePU = 32;
		ZeroBiasPrescale = 3 * 148; //L1 Prescale * HLT Prescale
		L1JetCorrection=false;
	}
	else if (whichFileAndLumiToUse==5) {

	// -- Run 179828, LS 50 - 70, PU=34:
		NumberOfUserdLumiSections = 21; 
		LumiForThisSetOfLumiSections = 0.529;
		switch(which_jet_seed_to_use){
			case 0:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1Tree_R179828_LS050-070.root";
				jobTag = "R179828_LS050-070";		
			break;
			case 1:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1SKIM_R179828_LS050-070_NoGCT_reEmul.root";
				jobTag = "R179828_LS050-070_NoGCT_reEmul";		
			break;
			case 2:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1SKIM_R179828_LS050-070_GCT_jet_seed_5GeV.root";
				jobTag = "R179828_LS050-070_GCT_jet_seed_5GeV";		
			break;
			default: cout << __LINE__ << ":" << __FILE__ << endl; 
				
		}
		AveragePU = 33;
		ZeroBiasPrescale = 3 * 148; //L1 Prescale * HLT Prescale
		L1JetCorrection=false;
	}
	else if (whichFileAndLumiToUse==6) {

	// -- Run 178803, LS 400 - 420, PU=18, (with bunch trains i.e. possible OOT PU): 
		NumberOfUserdLumiSections = 21; 
		LumiForThisSetOfLumiSections = 0.131;
		switch(which_jet_seed_to_use){
			case 0:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1Tree_R178803_LS400-420.root";
				jobTag = "R178803_LS400-420";
			break;
			case 1:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1SKIM_R178803_LS400-420._NoGCT_reEmul.root";
				jobTag = "R178803_LS400-420_NoGCT_reEmul";
			break;
			case 2:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1SKIM_R178803_LS400-420_GCT_jet_seed_5GeV.root";
				jobTag = "R178803_LS400-420_GCT_jet_seed_5GeV";
			break;
			default: cout << __LINE__ << ":" << __FILE__ << endl; 
				
		}
		AveragePU = 18;
		ZeroBiasPrescale = 29483;
		L1JetCorrection=false;
	}
	else if (whichFileAndLumiToUse==7) {
	// -- MC Run
	   int nevt = 2497500;
		if(procNevts>0) nevt = procNevts;
		NumberOfUserdLumiSections = (float)nevt/(10.*11246.* 23.3);  //lumi section time/number of mc events in file 23.3 gets removed below.
		LumiForThisSetOfLumiSections = 0.529;
		//LumiForThisSetOfLumiSections = targetlumi; //make scale factor 1.0
		switch(which_jet_seed_to_use){
			case 0:
				L1NtupleFileName = "/obsidian/users/winer/cms/L1Trigger/Simulations/L1SKIM_MC_7TeV_Ave32-v3_AllSet0_NoGCT_reEmul.root";
				jobTag = "MCv2_Ave32_Set0";
			break;
			case 1:
				L1NtupleFileName = "/obsidian/users/winer/cms/L1Trigger/Simulations/L1SKIM_MC_7TeV_Ave32-v3_AllSet0_NoGCT_reEmul.root";
				jobTag = "MCv3_Ave32_AllSet0_NoGCT_reEmul";
			break;
			case 2:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1SKIM_MC_7TeV_Ave32-v3_AllSet0_GCT_jet_seed_5GeV.root";
				jobTag = "MCv3_Ave32_AllSet0_GCT_jet_seed_5GeV";
			break;
			default: cout << __LINE__ << ":" << __FILE__ << endl; 
				
		}
		AveragePU = 32;
		ZeroBiasPrescale = 1;  //MC No prescale
		L1JetCorrection=false;
	}	
	else if (whichFileAndLumiToUse==8) {
	// -- MC Run
	   int nevt = 360000;
		if(procNevts>0) nevt = procNevts;
		NumberOfUserdLumiSections = (float)nevt/(2808.*11246.*23.3);  //lumi section time/number of mc events in file 23.3 gets removed below.
		LumiForThisSetOfLumiSections = 200.;
		//LumiForThisSetOfLumiSections = targetlumi; //make scale factor 1.0
		switch(which_jet_seed_to_use){
			case 0:
				L1NtupleFileName = "";
				jobTag = "MCv2_Ave32_Set0";
			break;
			case 1:
				L1NtupleFileName = "/obsidian/users/winer/cms/L1Trigger/Simulations/L1SKIM_MinBias_14TeV_StdGeom_PU50_All0000_NoGCT_reEmul.root";
				jobTag = "MinBias_14TeV_StdGeom_PU50_All0000_NoGCT_reEmul_Test2";
			break;
			case 2:
				L1NtupleFileName = "/uscms_data/d1/winer/batch_tutorial/L1SKIM_MinBias_14TeV_StdGeom_PU50_All0000_GCT_jet_seed_5GeV.root";
				jobTag = "MinBias_14TeV_StdGeom_PU50_All0000_GCT_jet_seed_5GeV";
			break;
			default: cout << __LINE__ << ":" << __FILE__ << endl; 
				
		}
		AveragePU = 50;
		ZeroBiasPrescale = 1;  //MC No prescale
		L1JetCorrection=false;
	}
	else if (whichFileAndLumiToUse==9) {
	// -- MC Run
		NumberOfUserdLumiSections = 63; 
		LumiForThisSetOfLumiSections = 0.1765; //units of e32

		switch(which_jet_seed_to_use){
			case 0:
				L1NtupleFileName = "/obsidian/users/winer/cms/L1Trigger/Simulations/Bristol/8TeV/ZeroBiasHPF1/2012HPF_66_v1/Test2.root";
				jobTag = "ZeroBiasHPF1_66";
			break;
			default: cout << __LINE__ << ":" << __FILE__ << endl; 
				
		}
		lsFileName = "getLumi_out_pixCorrLumi_66PU_stdCorr.txt";
		AveragePU = 66;
		ZeroBiasPrescale = 92;  //MC No prescale
		L1JetCorrection=false;
	}
	else if (whichFileAndLumiToUse==10) {
	// -- MC Run
		NumberOfUserdLumiSections = 143; 
		LumiForThisSetOfLumiSections = 0.061; // Note: wide range in this file...average is not a good thing //units of e32

		switch(which_jet_seed_to_use){
			case 0:
				L1NtupleFileName = "/obsidian/users/winer/cms/L1Trigger/Simulations/Bristol/8TeV/ZeroBiasHPF1/2012HPF_45_v1/Test.root";
				jobTag = "ZeroBiasHPF1_45";
			break;
			default: cout << __LINE__ << ":" << __FILE__ << endl; 
				
		}
		lsFileName = "getLumi_out_pixCorrLumi_45PU_stdCorr.txt";
		AveragePU = 45;
		ZeroBiasPrescale = 92;  //MC No prescale
		L1JetCorrection=false;
	} 
	else if (whichFileAndLumiToUse==11) {
	// -- MC Run
		NumberOfUserdLumiSections = 143; 
		LumiForThisSetOfLumiSections = 0.061; // Note: wide range in this file...average is not a good thing //units of e32

		switch(which_jet_seed_to_use){
			case 0:
				L1NtupleFileName = "/obsidian/users/winer/cms/L1Trigger/Simulations/Bristol/8TeV/UpgradeAlgos/ZeroBiasHPF1/2012HPF_45_v1/Test.root";
				jobTag = "UpgradeAlgo_ZeroBiasHPF1_45";
			break;
			default: cout << __LINE__ << ":" << __FILE__ << endl; 
				
		}
		lsFileName = "getLumi_out_pixCorrLumi_45PU_stdCorr.txt";
		AveragePU = 45;
		ZeroBiasPrescale = 92;  //MC No prescale
		L1JetCorrection=false;
	} 	
	else if (whichFileAndLumiToUse==12) {
	// -- MC Run
	   int nevt = 3100;
		if(procNevts>0) nevt = procNevts;
		NumberOfUserdLumiSections = (float)nevt/(2808.*11246.*23.3);  //lumi section time/number of mc events in file 23.3 gets removed below.
		LumiForThisSetOfLumiSections = 200.;

		switch(which_jet_seed_to_use){
			case 0:
				L1NtupleFileName = "/obsidian/users/winer/cms/L1Trigger/Simulations/L1NT_myminbias_8TeV_pu66_534_all.root";
				jobTag = "RichardTest";
			break;
			default: cout << __LINE__ << ":" << __FILE__ << endl; 
				
		}
		AveragePU = 66;
		ZeroBiasPrescale = 1;  //MC No prescale
		L1JetCorrection=false;
	} 	 	 		 	 
	else {
		std::cout << std::endl << "ERROR: Please define a ntuple file which is in the allowed range! You did use: whichFileAndLumiToUse = " << whichFileAndLumiToUse << " This is not in the allowed range" << std::endl << std::endl;
	}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	std::stringstream histos;
	std::string MenuPar = "Menu2015";
	histos << "L1RateHist_" << jobTag << "_" << MenuPar << "_" << usedL1Menu << "_" << targetlumi << "_" << AveragePU << "_" << LumiForThisSetOfLumiSections << "_rates.root";
        TString outHistName = histos.str();
	TFile* outHist = new TFile(outHistName,"RECREATE");
	outHist->cd();


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//HFW   For comparison to Chris L.
/*
//Cross
	h_SingleMu_ETM_byThreshold  = new TH1F("h_SingleMu_ETM_byThreshold","h_SingleMu_ETM_byThreshold",140,0.0,140.);
	h_SingleMu_CJet_byThreshold = new TH1F("h_SingleMu_CJet_byThreshold","h_SingleMu_CJet_byThreshold",140,0.0,140.);
	h_SingleEG_ETM_byThreshold  = new TH1F("h_SingleEG_ETM_byThreshold","h_SingleEG_ETM_byThreshold",63,0.0,63.);
	h_SingleEG_CJet_byThreshold = new TH1F("h_SingleEG_CJet_byThreshold","h_SingleEG_CJet_byThreshold",63,0.0,63.);

	h_Mu_EG_byThreshold = new TH1F("h_Mu_EG_byThreshold","h_Mu_EG_byThreshold",140,0.0,140.);
	h_EG_Mu_byThreshold = new TH1F("h_EG_Mu_byThreshold","h_EG_Mu_byThreshold",63,0.0,63.);

//Jets
	h_SingleJet_byThreshold      = new TH1F("h_SingleJet_byThreshold","h_SingleJet_byThreshold",400,0.0,400.);
   h_DoubleJet_byThreshold      = new TH1F("h_DoubleJet_byThreshold","h_DoubleJet_byThreshold",400,0.0,400.);
	h_QuadJetCentral_byThreshold = new TH1F("h_QuadJetCentral_byThreshold","h_QuadJetCentral_byThreshold",400,0.0,400.);
	h_SingleTau_byThreshold      = new TH1F("h_SingleTau_byThreshold","h_SingleTau_byThreshold",400,0.0,400.);
	h_DoubleTau_byThreshold      = new TH1F("h_DoubleTau_byThreshold","h_DoubleTau_byThreshold",400,0.0,400.);

//Sums
	h_HTT_byThreshold = new TH1F("h_HTT_byThreshold","h_HTT_byThreshold",750,0.0,750.);
	h_ETM_byThreshold = new TH1F("h_ETM_byThreshold","h_ETM_byThreshold",750,0.0,750.);
//EGamma
	h_SingleEG_byThreshold    = new TH1F("h_SingleEG_byThreshold","h_SingleEG_byThreshold",63,0.0,63.);
	h_SingleIsoEG_byThreshold = new TH1F("h_SingleIsoEG_byThreshold","h_SingleIsoEG_byThreshold",63,0.0,63.);
	h_DoubleEG_byThreshold    = new TH1F("h_DoubleEG_byThreshold","h_DoubleEG_byThreshold",63,0.0,63.);

//Muons
	h_SingleMu_byThreshold = new TH1F("h_SingleMu_byThreshold","h_SingleMu_byThreshold",140,0.0,140.);
	h_DoubleMu_byThreshold = new TH1F("h_DoubleMu_byThreshold","h_DoubleMu_byThreshold",140,0.0,140.);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//HFW  Take Care with bin edge effects in a better way.

//Cross
	h_SingleMu_ETM_byThreshold  = new TH1F("h_SingleMu_ETM_byThreshold","h_SingleMu_ETM_byThreshold",141,-0.5,140.5);
	h2_SingleMu_ETM_byThreshold = new TH2F("h2_SingleMu_ETM_byThreshold","h2_SingleMu_ETM_byThreshold",141,-0.5,140.5,201,-0.25,100.25);
	h_SingleMu_CJet_byThreshold = new TH1F("h_SingleMu_CJet_byThreshold","h_SingleMu_CJet_byThreshold",141,-0.5,140.5);
	h2_SingleMu_CJet_byThreshold= new TH2F("h2_SingleMu_CJet_byThreshold","h2_SingleMu_CJet_byThreshold",141,-0.5,140.5,51,-2.,202.);
	h_SingleEG_ETM_byThreshold  = new TH1F("h_SingleEG_ETM_byThreshold","h_SingleEG_ETM_byThreshold",64,-0.5,63.5);
	h2_SingleEG_ETM_byThreshold = new TH2F("h2_SingleEG_ETM_byThreshold","h2_SingleEG_ETM_byThreshold",64,-0.5,63.5,201,-0.25,100.25);
	h_SingleEG_CJet_byThreshold = new TH1F("h_SingleEG_CJet_byThreshold","h_SingleEG_CJet_byThreshold",64,-0.5,63.5);
	h2_SingleEG_CJet_byThreshold= new TH2F("h2_SingleEG_CJet_byThreshold","h2_SingleEG_CJet_byThreshold",64,-0.5,63.5,51,-2.,202.);

	h_Mu_EG_byThreshold = new TH1F("h_Mu_EG_byThreshold","h_Mu_EG_byThreshold",141,-0.5,140.5);
	h_EG_Mu_byThreshold = new TH1F("h_EG_Mu_byThreshold","h_EG_Mu_byThreshold",64,-0.5,63.5);
   h2_Mu_EG_byThreshold= new TH2F("h2_Mu_EG_byThreshold","h2_Mu_EG_byThreshold",141,-0.5,140.5,64,-0.5,63.5);
	
//Jets
	h_SingleJet_byThreshold      = new TH1F("h_SingleJet_byThreshold","h_SingleJet_byThreshold",101,-2.,402.);
   h_DoubleJet_byThreshold      = new TH1F("h_DoubleJet_byThreshold","h_DoubleJet_byThreshold",101,-2.,402.);
	h_QuadJetCentral_byThreshold = new TH1F("h_QuadJetCentral_byThreshold","h_QuadJetCentral_byThreshold",101,-2.,402.);
	h_SingleTau_byThreshold      = new TH1F("h_SingleTau_byThreshold","h_SingleTau_byThreshold",101,-2.,402.);
	h_DoubleTau_byThreshold      = new TH1F("h_DoubleTau_byThreshold","h_DoubleTau_byThreshold",101,-2.,402.);

//Sums 
	h_HTT_byThreshold = new TH1F("h_HTT_byThreshold","h_HTT_byThreshold",1601,-0.25,800.25);
	h_ETM_byThreshold = new TH1F("h_ETM_byThreshold","h_ETM_byThreshold",201 ,-0.25,100.25);
//EGamma
	h_SingleEG_byThreshold    = new TH1F("h_SingleEG_byThreshold","h_SingleEG_byThreshold",64,-0.5,63.5);
	h_SingleIsoEG_byThreshold = new TH1F("h_SingleIsoEG_byThreshold","h_SingleIsoEG_byThreshold",64,-0.5,63.5);
	h_DoubleEG_byThreshold    = new TH1F("h_DoubleEG_byThreshold","h_DoubleEG_byThreshold",64,-0.5,63.5);
	h2_DoubleEG_byThreshold   = new TH2F("h2_DoubleEG_byThreshold","h2_DoubleEG_byThreshold",64,-0.5,63.5,64,-0.5,63.5);

//Muons
	h_SingleMu_byThreshold = new TH1F("h_SingleMu_byThreshold","h_SingleMu_byThreshold",141,-0.5,140.5);
	h_DoubleMu_byThreshold = new TH1F("h_DoubleMu_byThreshold","h_DoubleMu_byThreshold",141,-0.5,140.5);
	h2_DoubleMu_byThreshold= new TH2F("h2_DoubleMu_byThreshold","h2_DoubleMu_byThreshold",141,-0.5,140.5,141,-0.5,140.5);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Plot of input quantities
   h_Mu_Nmu = new TH1F("h_Mu_Nmu","Num iso Mu",5,-0.5,4.5);
   h_Mu_Et  = new TH1F("h_Mu_Et","iso Mu Pt",141,-0.5,140.5);
   h_Mu_Eta = new TH1F("h_Mu_Eta","iso Mu Eta",100,-2.5,2.5);
   h_Mu_Phi = new TH1F("h_Mu_Phi","iso Mu Phi",100,0.0,TMath::TwoPi());


   h_isoEG_Nele= new TH1F("h_isoEG_Nele","Num iso EG",5,-0.5,4.5);
   h_isoEG_Et  = new TH1F("h_isoEG_Et","iso EG Et",64,-0.5,63.5);
   h_isoEG_Eta = new TH1F("h_isoEG_Eta","iso EG Eta",22,-0.5,21.5);
   h_isoEG_Phi = new TH1F("h_isoEG_Phi","iso EG Phi",18,-0.5,17.5);

   h_nIsoEG_Nele= new TH1F("h_nIsoEG_Nele","Num nonIso EG",5,-0.5,4.5);
   h_nIsoEG_Et  = new TH1F("h_nIsoEG_Et","nonIso EG Et",64,-0.5,63.5);
   h_nIsoEG_Eta = new TH1F("h_nIsoEG_Eta","nonIso EG Eta",22,-0.5,21.5);
   h_nIsoEG_Phi = new TH1F("h_nIsoEG_Phi","nonIso EG Phi",18,-0.5,17.5);
	
   h_CJet_Njet= new TH1F("h_CJet_Njet","Num Central Jets",5,-0.5,4.5);
   h_CJet_Et  = new TH1F("h_CJet_Et","Central Jet Et",257,-0.5,256.5);
   h_CJet_Eta = new TH1F("h_CJet_Eta","Central Jet Eta",22,-0.5,21.5);
   h_CJet_Phi = new TH1F("h_CJet_Phi","Central Jet Phi",18,-0.5,17.5);

   h_FJet_Njet= new TH1F("h_FJet_Njet","Num Forward Jets",5,-0.5,4.5);
   h_FJet_Et  = new TH1F("h_FJet_Et","Forward Jet Et",257,-0.5,256.5);
   h_FJet_Eta = new TH1F("h_FJet_Eta","Forward Jet Eta",22,-0.5,21.5);
   h_FJet_Phi = new TH1F("h_FJet_Phi","Forward Jet Phi",18,-0.5,17.5);

   h_TJet_Njet= new TH1F("h_TJet_Njet","Num Tau Jets",5,-0.5,4.5);
	h_TJet_Et  = new TH1F("h_TJet_Et","Tau Jet Et",257,-0.5,256.5);
   h_TJet_Eta = new TH1F("h_TJet_Eta","Tau Jet Eta",22,-0.5,21.5);
   h_TJet_Phi = new TH1F("h_TJet_Phi","Tau Jet Phi",18,-0.5,17.5);
	
   h_Sum_ETT   = new TH1F("h_Sum_ETT","ET Total",1601,-0.25,800.25);
	h_Sum_ETM   = new TH1F("h_Sum_ETM","ET Miss",201,-0.25,100.25);
	h_Sum_PhiETM= new TH1F("h_Sum_PhiETM","ET Miss Phi",100,-TMath::Pi(),TMath::Pi());	

   h_Sum_HTT   = new TH1F("h_Sum_HTT","HT Total",1601,-0.25,800.25);
	h_Sum_HTM   = new TH1F("h_Sum_HTM","HT Miss",101,-0.25,200.25);
	h_Sum_PhiHTM= new TH1F("h_Sum_PhiHTM","HT Miss Phi",100,-TMath::Pi(),TMath::Pi());


// Histograms for the Trigger Results.
   h_Trig = new TH1F("h_Trig","h_Trig",N128,-0.5,N128+0.5);
	h_All  = new TH1F("h_All","h_All",N128,-0.5,N128+0.5);
	h_Pure = new TH1F("h_Pure","h_Pure",N128,-0.5,N128+0.5);
	h_Corr = new TH2F("h_Corr","h_Corr",N128,-0.5,N128+0.5,N128,-0.5,N128+0.5);
	h_Cumm = new TH1F("h_Cumm","h_Cumm",N128,-0.5,N128+0.5);

// Time to dump the configuration 
	std::cout << std::endl << "L1 menu used = " << usedL1Menu << std::endl;
	std::cout << std::endl << "Target Luminosity = " << targetlumi << std::endl << std::endl;
	std::cout << std::endl << "Using: whichFileAndLumiToUse = " << whichFileAndLumiToUse << std::endl;
	std::cout << "  NumberOfUserdLumiSections        = " << NumberOfUserdLumiSections << std::endl;
	std::cout << "  LumiForThisSetOfLumiSections     = " << LumiForThisSetOfLumiSections << std::endl;
	std::cout << "  L1NtupleFileName                 = " << L1NtupleFileName << std::endl;
	std::cout << "  lsFileName                       = " << lsFileName.Data() << std::endl;
	std::cout << "  AveragePU                        = " << AveragePU << std::endl;
	std::cout << "  Use L1Extra Quantities           = " << useL1Extra << std::endl;
   std::cout << "  Calculate Threshold Plots        = " << calcThreshold  << std::endl;
   std::cout << "  L1JetCorrections (for 2011 data) = " << L1JetCorrection << std::endl << std::endl;


//  Do the heavy lifting
	L1Menu2015 a(usedL1Menu,targetlumi,NumberOfUserdLumiSections,LumiForThisSetOfLumiSections,
			L1NtupleFileName,AveragePU,ZeroBiasPrescale,L1JetCorrection);
	a.Open(L1NtupleFileName);
	a.Loop(calcThreshold,useL1Extra,lsFileName,procNevts);


// Table header	
	printf("L1Bit      L1SeedName   pre-sc     rate (kHz)        cumulative (kHz)        pure (kHz)  \n");
		
	Float_t totalrate     = 0.;
	Float_t totalrate_err = 0.;	
   Float_t finalL1Rate = h_All -> GetBinContent(kOFFSET+1);
	Float_t finalL1Rate_err = h_All -> GetBinError(kOFFSET+1);
	Float_t pureL1Rate = h_Pure -> GetBinContent(kOFFSET+1);
	Float_t pureL1Rate_err = h_Pure -> GetBinError(kOFFSET+1);

	for (Int_t k=1; k < kOFFSET+1; k++) {  // -- kOFFSET now contains the number of triggers we have calculated
		TString name = h_All -> GetXaxis() -> GetBinLabel(k);

		Float_t rate = h_All -> GetBinContent(k);
		Float_t err_rate  = h_All -> GetBinError(k);

		Float_t cumm_rate = h_Cumm -> GetBinContent(k);
      Float_t err_cumm_rate  = h_Cumm -> GetBinError(k);

		Float_t pure_rate = h_Pure -> GetBinContent(k);
		Float_t err_pure_rate = h_Pure -> GetBinError(k);

		std::string L1namest = (std::string)name;
		std::map<std::string, int>::const_iterator it = a.Prescales.find(L1namest);
		Float_t pre;
		if (it == a.Prescales.end() ) {
			std::cout << " --- SET P = 1 FOR SEED :  " << L1namest << std::endl;
			pre = 1;
		}
		else {
			pre = it -> second;
		}

      totalrate += rate;
		totalrate_err += err_rate*err_rate;
		
//print the results
      printf("%2i %20s   %2i %8.2f  +/- %5.2f  %8.2f  +/- %5.2f  %8.2f  +/- %5.2f \n",a.L1BitNumber(L1namest),name.Data(),(int)pre,rate,err_rate,cumm_rate,err_cumm_rate,pure_rate,err_pure_rate);				
		
	}

   printf("\n Total L1 Rate (with overlaps)    = %8.2f +/- %5.2f  kHz\n",finalL1Rate, finalL1Rate_err); 
	printf(  " Total L1 Rate (without overlaps) = %8.2f +/- %5.2f  kHz\n",totalrate,sqrt(totalrate_err));
	printf(  " Total L1 Rate (pure triggers)    = %8.2f +/- %5.2f  kHz\n",pureL1Rate, pureL1Rate_err);
		
	outHist->Write();
	outHist->Close();
}

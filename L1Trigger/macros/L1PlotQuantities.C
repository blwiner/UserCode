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

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//HFW
//Cross


// Plots for Trigger Quantities
TH1F *h_Mu_Nmu,      *h_Mu_Et,     *h_Mu_Eta,     *h_Mu_Phi;
TH1F *h_isoEG_Nele,  *h_isoEG_Et,  *h_isoEG_Eta,  *h_isoEG_Phi;
TH1F *h_nIsoEG_Nele, *h_nIsoEG_Et, *h_nIsoEG_Eta, *h_nIsoEG_Phi;
TH1F *h_CJet_Njet, *h_CJet_Et, *h_CJet_Eta, *h_CJet_Phi;
TH1F *h_FJet_Njet, *h_FJet_Et, *h_FJet_Eta, *h_FJet_Phi;
TH1F *h_TJet_Njet, *h_TJet_Et, *h_TJet_Eta, *h_TJet_Phi;
TH1F *h_Sum_ETT,   *h_Sum_ETM, *h_Sum_PhiETM;
TH1F *h_Sum_HTT,   *h_Sum_HTM, *h_Sum_PhiHTM;

TH1F *h_deltaPhi_MuJet; 
TH1F *h_deltaEta_MuJet; 
TH1F *h_deltaR_MuJet  ; 
TH1F *h_deltaRMin_MuJet;
TH1F *h_deltaR_MuJet_PtEt  ; 
TH1F *h_deltaRMin_MuJet_PtEt;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

Double_t degree(Double_t radian) {
	if (radian<0)
		return 360.+(radian/TMath::Pi()*180.);
	else
		return radian/TMath::Pi()*180.;
}


Int_t etaINjetCoord(Double_t eta){
	size_t etaIdx = 0.;
	for (size_t idx=0; idx<ETABINS; idx++) {
		if (eta>=ETABIN[idx] and eta<ETABIN[idx+1])
			etaIdx = idx;
	}
	return int(etaIdx);
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

class L1PlotQuantities : public L1Ntuple {
	public :

	L1PlotQuantities()  {}

	~L1PlotQuantities() {}



	void MyInit();
	
	L1Analysis::L1AnalysisDataFormat myEvt_;

	void Loop(Bool_t useL1Extra, const int n_events_=-1);
	void fillDataStructure(Bool_t UseL1Extra=true);

	private :


};


void L1PlotQuantities::fillDataStructure(Bool_t useL1Extra) {
   
	 // printf("Entering fillDataStucture \n");
	 myEvt_.Reset();

    if(useL1Extra) {

// Grab the iso first
		 for(unsigned int i=0; i<l1extra_->nIsoEm; i++) {


      	 myEvt_.Bxel.push_back(l1extra_->isoEmBx.at(i));
      	 myEvt_.Etel.push_back(l1extra_->isoEmEt.at(i));
//      	 myEvt_.Phiel.push_back(phiINjetCoord(l1extra_->isoEmPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
//      	 myEvt_.Etael.push_back(etaINjetCoord(l1extra_->isoEmEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Phiel.push_back(l1extra_->isoEmPhi.at(i)); 
      	 myEvt_.Etael.push_back(l1extra_->isoEmEta.at(i)); 
       	 myEvt_.Isoel.push_back(true);

// Histogram Quantities
          h_isoEG_Nele->Fill(l1extra_->nIsoEm);			 
          if(l1extra_->isoEmBx.at(i)==0) {
			    h_isoEG_Et->Fill(myEvt_.Etel.at(myEvt_.Nele));
				 h_isoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele));
				 h_isoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele));
			 }	 
			 myEvt_.Nele++;

		 }
		 
		 for(unsigned int i=0; i<l1extra_->nNonIsoEm; i++) {

      	 myEvt_.Bxel.push_back(l1extra_->nonIsoEmBx.at(i));
      	 myEvt_.Etel.push_back(l1extra_->nonIsoEmEt.at(i));
//      	 myEvt_.Phiel.push_back(phiINjetCoord(l1extra_->nonIsoEmPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
//      	 myEvt_.Etael.push_back(etaINjetCoord(l1extra_->nonIsoEmEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Phiel.push_back(l1extra_->nonIsoEmPhi.at(i)); 
      	 myEvt_.Etael.push_back(l1extra_->nonIsoEmEta.at(i)); 
			 myEvt_.Isoel.push_back(false);

// Histogram Quantities
          h_nIsoEG_Nele->Fill(l1extra_->nNonIsoEm);			 
          if(l1extra_->nonIsoEmBx.at(i)==0) {
			    h_nIsoEG_Et-> Fill(myEvt_.Etel.at(myEvt_.Nele));
				 h_nIsoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele));
				 h_nIsoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele));
			 }	 
			 myEvt_.Nele++;


		 }		 
		 // printf("Number of Electrons in myEvt %i \n",myEvt_.Nele);

		 for(unsigned int i=0; i< l1extra_->nCenJets; i++) {

      	 
      	 myEvt_.Bxjet.push_back(l1extra_->cenJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1extra_->cenJetEt.at(i));
//      	 myEvt_.Phijet.push_back(phiINjetCoord(l1extra_->cenJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
//      	 myEvt_.Etajet.push_back(etaINjetCoord(l1extra_->cenJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Phijet.push_back(l1extra_->cenJetPhi.at(i)); 
      	 myEvt_.Etajet.push_back(l1extra_->cenJetEta.at(i)); 
      	 myEvt_.Taujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);
			 
// Histogram Quantities
          h_CJet_Njet->Fill(l1extra_->nCenJets);			 
          if(l1extra_->cenJetBx.at(i)==0) {
			    h_CJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet));
				 h_CJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet));
				 h_CJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet));
			 }
			 myEvt_.Njet++;	 
		 }
		 
		 
		 for(unsigned int i=0; i< l1extra_->nFwdJets; i++) {

      	 
      	 myEvt_.Bxjet.push_back(l1extra_->fwdJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1extra_->fwdJetEt.at(i));
//      	 myEvt_.Phijet.push_back(phiINjetCoord(l1extra_->fwdJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
//      	 myEvt_.Etajet.push_back(etaINjetCoord(l1extra_->fwdJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Phijet.push_back(l1extra_->fwdJetPhi.at(i)); 
      	 myEvt_.Etajet.push_back(l1extra_->fwdJetEta.at(i)); 
      	 myEvt_.Taujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(true);

// Histogram Quantities
          h_FJet_Njet->Fill(l1extra_->nFwdJets);			 
          if(l1extra_->fwdJetBx.at(i)==0) {
			    h_FJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet));
				 h_FJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet));
				 h_FJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet));
			 }	 
			 myEvt_.Njet++;
		 }
		 
		 for(unsigned int i=0; i< l1extra_->nTauJets; i++) {

      	 myEvt_.Bxjet.push_back(l1extra_->tauJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1extra_->tauJetEt.at(i));
//      	 myEvt_.Phijet.push_back(phiINjetCoord(l1extra_->tauJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
//      	 myEvt_.Etajet.push_back(etaINjetCoord(l1extra_->tauJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Phijet.push_back(l1extra_->tauJetPhi.at(i)); 
      	 myEvt_.Etajet.push_back(l1extra_->tauJetEta.at(i)); 
      	 myEvt_.Taujet.push_back(true);
      	 myEvt_.Fwdjet.push_back(false);

// Histogram Quantities
          h_TJet_Njet->Fill(l1extra_->nTauJets);			 
          if(l1extra_->tauJetBx.at(i)==0) {
			    h_TJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet));
				 h_TJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet));
				 h_TJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet));
			 }	 
			 myEvt_.Njet++;

		 }		 		 
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


void L1PlotQuantities::MyInit() {


}

void L1PlotQuantities::Loop(Bool_t useL1Extra, const int n_events_) {

	
	const Int_t nevents = (n_events_ < 0) ? GetEntries() : n_events_;
   int cnt = 0;

	for (Long64_t i=0; i<nevents; i++)
	{     
	//load the i-th event
		Long64_t ientry = LoadTree(i); if (ientry < 0) break;
		GetEntry(i);

      cnt++;
      if(cnt%(int)pow(10.,(double)((int)log10((double)cnt)))==0) printf("Event Number %i\n",cnt);


// Fill my event data
      fillDataStructure(useL1Extra);


//  Determine Delta R between muons and nearest jet
	   Int_t Nmu = myEvt_.Nmu;
	   for (Int_t imu=0; imu < Nmu; imu++) { 
		   		
		   if (myEvt_.Bxmu.at(imu) != 0) continue;
		   Float_t ptM  = myEvt_.Ptmu.at(imu);			
		   Int_t   qual = myEvt_.Qualmu.at(imu); 
			Float_t etaM = myEvt_.Etamu.at(imu);
			Float_t phiM = myEvt_.Phimu.at(imu);
			if(phiM>TMath::Pi()) phiM -= TMath::TwoPi(); //make sure -pi to pi      

// Only look at high quality muons
		   if ( qual > 4) {
			    
				 Float_t minDeltaR = 999.;
				 Float_t minDeltaR_PtEt = 999.;
			    for(Int_t nj=0; nj<myEvt_.Njet; nj++) {
				 
				     if( !myEvt_.Taujet.at(nj) && !myEvt_.Fwdjet.at(nj) && myEvt_.Bxjet.at(nj)==0) {
  
                    //printf("Muon %i  Jet  %i\n",imu,nj);
// Calculate the delta phi between the jet and the muon
                    Float_t deltaPhi = fabs(phiM - myEvt_.Phijet.at(nj));
						  if(deltaPhi>TMath::Pi()) deltaPhi = TMath::TwoPi() - deltaPhi;
						  //printf("phiM %f  phiJ %f Delta Phi %f \n",phiM,myEvt_.Phijet.at(nj),deltaPhi);
                    h_deltaPhi_MuJet->Fill(deltaPhi);

// Calculate the delta eta between the jet and the muon
                    Float_t deltaEta = fabs(etaM - myEvt_.Etajet.at(nj));
						  //printf("etaM %f  etaJ %f Delta eta %f \n",etaM,myEvt_.Etajet.at(nj),deltaEta);
						  h_deltaEta_MuJet->Fill(deltaEta);
						  
// Now delta R
                    Float_t deltaR = sqrt(deltaPhi*deltaPhi + deltaEta * deltaEta);					
						 // printf("Delta R %f \n",deltaR);	  
                    h_deltaR_MuJet->Fill(deltaR);

// Keep track if this is the closest Jet					  
                    if(deltaR<minDeltaR) {
						      minDeltaR = deltaR;
						  }
						  
// place some cuts on the muon Pt and Jet Et						  
						  if(ptM>=100. && myEvt_.Etjet.at(nj)>=10.) {
						     h_deltaR_MuJet_PtEt->Fill(deltaR);
                       if(deltaR<minDeltaR_PtEt) {
 						        minDeltaR_PtEt = deltaR;
						     }
						  }
					  
					  } //demand it is a central jet
				 } //loop over jets

// Fill minDeltaR between jet and muon				 
				 h_deltaRMin_MuJet->Fill(minDeltaR);
				 h_deltaRMin_MuJet_PtEt->Fill(minDeltaR_PtEt);
				 
			} // muon is high quality
			
	  }//loop over muons      

//     printf("=====================================\n");
	}  // end evt loop
        	
}



void Run(Bool_t useL1Extra=true,Int_t whichDataSetToUse=1,int whichFilesToUse=0,Int_t pNevts = -1) {

// Make sure we get the errors correct on histograms   
	TH1::SetDefaultSumw2();


//whichFilesToUse:

//0: default
//1: reEmul
//2: 5 GeV 


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Using the 10-bunches High PU run, 179828. In this run ETM30 and HTT50 were enabled (they were not enabled in the 1-bunch run).
	std::string L1NtupleFileName="";
	std::string jobTag="";
	Int_t procNevts = pNevts;
	TString lsFileName = "";
	

	if (whichDataSetToUse==1) {

		switch(whichFilesToUse){
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

	}
	else if (whichDataSetToUse==2) {

	// -- Run 179828, LS 300 - 320, PU=29:

		switch(whichFilesToUse){
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

	}
	else if (whichDataSetToUse==3) {

	// -- Run 179828, LS 270 - 290, PU=30:
		switch(whichFilesToUse){
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
	}
	else if (whichDataSetToUse==4) {

	// -- Run 179828, LS 140 - 160, PU=33:
		switch(whichFilesToUse){
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
	}
	else if (whichDataSetToUse==5) {

	// -- Run 179828, LS 50 - 70, PU=34:
		switch(whichFilesToUse){
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
	}
	else if (whichDataSetToUse==6) {

	// -- Run 178803, LS 400 - 420, PU=18, (with bunch trains i.e. possible OOT PU): 
		switch(whichFilesToUse){
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
	}
	else if (whichDataSetToUse==7) {
	// -- MC Run

		switch(whichFilesToUse){
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
	}	
	else if (whichDataSetToUse==8) {
	// -- MC Run
		switch(whichFilesToUse){
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
	}
	else if (whichDataSetToUse==9) {
	// -- MC Run

		switch(whichFilesToUse){
			case 0:
				L1NtupleFileName = "/obsidian/users/winer/cms/L1Trigger/Simulations/Bristol/8TeV/ZeroBiasHPF1/2012HPF_66_v1/Test2.root";
				jobTag = "ZeroBiasHPF1_66";
			break;
			default: cout << __LINE__ << ":" << __FILE__ << endl; 
				
		}
	}
	else if (whichDataSetToUse==10) {
	// -- MC Run

		switch(whichFilesToUse){
			case 0:
				L1NtupleFileName = "/obsidian/users/winer/cms/L1Trigger/Simulations/Bristol/8TeV/ZeroBiasHPF1/2012HPF_45_v1/Test.root";
				jobTag = "ZeroBiasHPF1_45";
			break;
			default: cout << __LINE__ << ":" << __FILE__ << endl; 
				
		}
	} 	 	 
	else {
		std::cout << std::endl << "ERROR: Please define a ntuple file which is in the allowed range! You did use: whichDataSetToUse = " << whichDataSetToUse << " This is not in the allowed range" << std::endl << std::endl;
	}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	std::stringstream histos;
	histos << "L1QntHist_" << jobTag << ".root";
        TString outHistName = histos.str();
	TFile* outHist = new TFile(outHistName,"RECREATE");
	outHist->cd();

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Plot of input quantities
   h_Mu_Nmu = new TH1F("h_Mu_Nmu","Num iso Mu",5,-0.5,4.5);
   h_Mu_Et  = new TH1F("h_Mu_Et","iso Mu Pt",141,-0.5,140.5);
   h_Mu_Eta = new TH1F("h_Mu_Eta","iso Mu Eta",100,-2.5,2.5);
   h_Mu_Phi = new TH1F("h_Mu_Phi","iso Mu Phi",100,0.0,TMath::TwoPi());

   h_isoEG_Nele= new TH1F("h_isoEG_Nele","Num iso EG",5,0.5,4.5);
   h_isoEG_Et  = new TH1F("h_isoEG_Et","iso EG Et",64,-0.5,63.5);
   h_isoEG_Eta = new TH1F("h_isoEG_Eta","iso EG Eta",22,-0.5,21.5);
   h_isoEG_Phi = new TH1F("h_isoEG_Phi","iso EG Phi",18,-0.5,17.5);

   h_nIsoEG_Nele= new TH1F("h_nIsoEG_Nele","Num nonIso EG",5,0.5,4.5);
   h_nIsoEG_Et  = new TH1F("h_nIsoEG_Et","nonIso EG Et",64,-0.5,63.5);
   h_nIsoEG_Eta = new TH1F("h_nIsoEG_Eta","nonIso EG Eta",22,-0.5,21.5);
   h_nIsoEG_Phi = new TH1F("h_nIsoEG_Phi","nonIso EG Phi",18,-0.5,17.5);
	
   h_CJet_Njet= new TH1F("h_CJet_Njet","Num Central Jets",5,0.5,4.5);
   h_CJet_Et  = new TH1F("h_CJet_Et","Central Jet Et",257,-0.5,256.5);
   h_CJet_Eta = new TH1F("h_CJet_Eta","Central Jet Eta",22,-0.5,21.5);
   h_CJet_Phi = new TH1F("h_CJet_Phi","Central Jet Phi",18,-0.5,17.5);

   h_FJet_Njet= new TH1F("h_FJet_Njet","Num Forward Jets",5,0.5,4.5);
   h_FJet_Et  = new TH1F("h_FJet_Et","Forward Jet Et",257,-0.5,256.5);
   h_FJet_Eta = new TH1F("h_FJet_Eta","Forward Jet Eta",22,-0.5,21.5);
   h_FJet_Phi = new TH1F("h_FJet_Phi","Forward Jet Phi",18,-0.5,17.5);

   h_TJet_Njet= new TH1F("h_TJet_Njet","Num Tau Jets",5,0.5,4.5);
	h_TJet_Et  = new TH1F("h_TJet_Et","Tau Jet Et",257,-0.5,256.5);
   h_TJet_Eta = new TH1F("h_TJet_Eta","Tau Jet Eta",22,-0.5,21.5);
   h_TJet_Phi = new TH1F("h_TJet_Phi","Tau Jet Phi",18,-0.5,17.5);
	
   h_Sum_ETT   = new TH1F("h_Sum_ETT","ET Total",1601,-0.25,800.25);
	h_Sum_ETM   = new TH1F("h_Sum_ETM","ET Miss",201,-0.25,100.25);
	h_Sum_PhiETM= new TH1F("h_Sum_PhiETM","ET Miss Phi",100,-TMath::Pi(),TMath::Pi());	

   h_Sum_HTT   = new TH1F("h_Sum_HTT","HT Total",1601,-0.25,800.25);
	h_Sum_HTM   = new TH1F("h_Sum_HTM","HT Miss",101,-0.25,200.25);
	h_Sum_PhiHTM= new TH1F("h_Sum_PhiHTM","HT Miss Phi",100,-TMath::Pi(),TMath::Pi());


// Separation between muons and jets
   h_deltaPhi_MuJet  = new TH1F("h_deltaPhi_MuJet","delta phi Separation mu-Jet",100,0.,TMath::Pi());
	h_deltaEta_MuJet  = new TH1F("h_deltaEta_MuJet","delta eta Separation mu-Jet",100,0.,6.);
	h_deltaR_MuJet    = new TH1F("h_deltaR_MuJet","delta R Separation mu-Jet",100,0.,6.5);
	h_deltaRMin_MuJet = new TH1F("h_deltaRMin_MuJet","Min delta R Separation mu-Jet",100,0.,6.5);

	h_deltaR_MuJet_PtEt    = new TH1F("h_deltaR_MuJet_PtEt","delta R Separation mu-Jet",100,0.,6.5);
	h_deltaRMin_MuJet_PtEt = new TH1F("h_deltaRMin_MuJet_PtEt","Min delta R Separation mu-Jet",100,0.,6.5);
	
//  Do the heavy lifting
	L1PlotQuantities a;
	a.Open(L1NtupleFileName);
	a.Loop(useL1Extra,procNevts);
		
	outHist->Write();
	outHist->Close();
}

#include "TObject.h"
#include "TMath.h"
#include "TFile.h"
#include "TText.h"
#include "TH2.h"
#include "TAxis.h"
#include "TString.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>


class EvaluateL1Menu {
	public :

	EvaluateL1Menu(){}

	~EvaluateL1Menu() {}


   void  LoadL1Menu(TString fileName);
   void  LoadThresholdPlots(TString fileName);	
	void  WriteL1Menu(TString fileName);
   void  DetermineThresholds();
	void  ScaleBandwidth(double scaleFactor);
   double FindThreshold(TH1F* byThreshold, double Threshold);
	
	typedef struct {
	   float primTh ;
	   float secTh;
	   float triTh;
	   float quadTh;
	   float etaCut;
	   int minQual;
		float bandwidth;
		int   scalable;
		bool  locked;
	} trigPar;
	
	std::map<std::string, trigPar> trigParList;
	std::map<std::string, int> Prescales;  //keep definitions used in L1Menu2012.
   std::map<std::string, int> BitMapping;

	private :

   TFile* histFile;
};

void EvaluateL1Menu::LoadL1Menu(TString fileName) {

// Open File
   printf("\n Reading L1 Menu File %s \n",fileName.Data());
   ifstream ifs( fileName );  
	
// Read through Menu
   while(ifs) {
	
	  TString algName;
	  ifs >> algName;
	  
	  ifs >> BitMapping[algName.Data()];	  
	  ifs >> Prescales[algName.Data()];
	  
	  ifs >> trigParList[algName.Data()].primTh;
	  ifs >> trigParList[algName.Data()].secTh;
	  ifs >> trigParList[algName.Data()].triTh;
	  ifs >> trigParList[algName.Data()].quadTh;
	  ifs >> trigParList[algName.Data()].etaCut;
	  ifs >> trigParList[algName.Data()].minQual;
	  ifs >> trigParList[algName.Data()].bandwidth;
	  ifs >> trigParList[algName.Data()].scalable;
	  ifs >> trigParList[algName.Data()].locked;	  
	  
//	  printf("Alg: %20s %2i %2i %6.2f %6.2f %6.2f %6.2f %6.2f %2i %6.2f %2i %2i\n",algName.Data(),BitMapping[algName.Data()],Prescales[algName.Data()],
//	         trigParList[algName.Data()].primTh,trigParList[algName.Data()].secTh,trigParList[algName.Data()].triTh,trigParList[algName.Data()].quadTh,
//				trigParList[algName.Data()].etaCut,trigParList[algName.Data()].minQual,trigParList[algName.Data()].bandwidth,trigParList[algName.Data()].scalable,trigParList[algName.Data()].locked);
	}

  return;
}

void EvaluateL1Menu::LoadThresholdPlots(TString fileName) {

  histFile = new TFile(fileName);
  printf("\n EvaluateL1Menu:: Loading rate vs Threshold plots from %s\n",fileName.Data());
  
  return;
}

void EvaluateL1Menu::WriteL1Menu(TString fileName) {

  FILE *outMenu = fopen(fileName,"w");
  
  printf("\n EvaluateL1Menu:: Writing New L1 Menu to %s\n \n",fileName.Data());
  for (std::map<std::string, trigPar>::iterator it=trigParList.begin(); it != trigParList.end(); it++) {
      TString algName = it->first;
		
      if(algName.Contains("L1")) fprintf(outMenu,"%20s %2i %2i %6.2f %6.2f %6.2f %6.2f %6.2f %2i %6.2f %2i %2i\n",algName.Data(),BitMapping[algName.Data()],Prescales[algName.Data()],
	         trigParList[algName.Data()].primTh,trigParList[algName.Data()].secTh,trigParList[algName.Data()].triTh,trigParList[algName.Data()].quadTh,
				trigParList[algName.Data()].etaCut,trigParList[algName.Data()].minQual,trigParList[algName.Data()].bandwidth,trigParList[algName.Data()].scalable,trigParList[algName.Data()].locked); 
  }

  fclose(outMenu);
  
  return;
}

void EvaluateL1Menu::DetermineThresholds() {

 
  printf("\n============= EvaluateL1Menu:: Finding new Thresholds ===============================================================\n");
  for (std::map<std::string, trigPar>::iterator it=trigParList.begin(); it != trigParList.end(); it++) {

      TString algName = it->first;		
      if(algName.Contains("L1")) {
        if(!trigParList[algName.Data()].locked  && Prescales[algName.Data()]>0) {
// Get name of threshold plot for this algorithm
          TString hname = "h_";
			 hname += algName(3,algName.Length());
			 hname += "_byThreshold";
          TH1F* tmp = (TH1F*)histFile->Get(hname)->Clone();		
			 
// Get dedicated threshold for this trigger
          double threshold = FindThreshold(tmp,trigParList[algName.Data()].bandwidth);			 
			 printf("Changing prim. threshold for %20s from %7.3f ---> %7.3f  to achieve target rate of %6.2f kHz\n",algName.Data(),trigParList[algName.Data()].primTh,threshold,trigParList[algName.Data()].bandwidth);
			 
			 float binWidth = tmp->GetBinWidth(1);
			 // use "floor" function to put non-primary thresholds on correct boundaries.
			 if(trigParList[algName.Data()].secTh>0.)  trigParList[algName.Data()].secTh  = binWidth*floor((threshold*(trigParList[algName.Data()].secTh/trigParList[algName.Data()].primTh))/binWidth+0.5);
			 if(trigParList[algName.Data()].triTh>0.)  trigParList[algName.Data()].triTh  = binWidth*floor((threshold*(trigParList[algName.Data()].triTh/trigParList[algName.Data()].primTh))/binWidth+0.5);
			 if(trigParList[algName.Data()].quadTh>0.) trigParList[algName.Data()].quadTh = binWidth*floor((threshold*(trigParList[algName.Data()].quadTh/trigParList[algName.Data()].primTh))/binWidth+0.5);
			 trigParList[algName.Data()].primTh = threshold;
		  } else if(trigParList[algName.Data()].locked) {
		    printf("Locked:: prim. threshold for %20s at   %7.3f \n",algName.Data(),trigParList[algName.Data()].primTh);
		  }  
		}  
  }
  printf("=====================================================================================================================\n");
  
  return;
}


void EvaluateL1Menu::ScaleBandwidth(double scaleFactor) {

 
  printf("\n============== EvaluateL1Menu:: Scaling Bandwidths ===================================\n");
  for (std::map<std::string, trigPar>::iterator it=trigParList.begin(); it != trigParList.end(); it++) {

      TString algName = it->first;		
      if(algName.Contains("L1")) {
		   if(trigParList[algName.Data()].scalable == 1 && !trigParList[algName.Data()].locked && Prescales[algName.Data()]>0) {
		     printf("Changing allocated bandwidth for %20s from %7.3f ---> %7.3f \n",algName.Data(),trigParList[algName.Data()].bandwidth,scaleFactor*trigParList[algName.Data()].bandwidth);
           trigParList[algName.Data()].bandwidth *= scaleFactor; 			  
		   } else if( (trigParList[algName.Data()].scalable == 0 || trigParList[algName.Data()].locked) && Prescales[algName.Data()]>0){
			  printf("Algorithm Rates are locked for   %20s at %7.3f \n",algName.Data(),trigParList[algName.Data()].bandwidth);
			}
		}  
  }
  printf("========================================================================================\n");
  
  return;
}


double EvaluateL1Menu::FindThreshold(TH1F* byThreshold, double targetRate) {

	double threshold = -1.;
   int i=1;
	bool fnd = false;
   while(i<byThreshold->GetNbinsX() && !fnd) {
	   if(byThreshold->GetBinContent(i)<targetRate) {
          threshold = byThreshold->GetBinCenter(i);
			 fnd = true;
			 
// Check whether previous bin was closer to target, if so, use it.
//          if((byThreshold->GetBinContent(i-1)-targetRate)
//			    (targetRate-byThreshold->GetBinContent(i))) {
//				    threshold = byThreshold->GetBinCenter(i-1);
//  			 }		 			  			 
		} else {
		    i++;
		}		   
	}
	if(!fnd) printf("WARNING: did not find an acceptable threshold\n");

  return threshold;
}

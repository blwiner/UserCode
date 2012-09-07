#ifndef __L1Analysis_L1AnalysisDataFormat_H__
#define __L1Analysis_L1AnalysisDataFormat_H__

#include <vector>
// #include <inttypes.h>
#include <TROOT.h>

namespace L1Analysis
{
  struct L1AnalysisDataFormat
  {
    L1AnalysisDataFormat(){Reset();};
    ~L1AnalysisDataFormat(){};
    
    void Reset()
    { 


    	//PSB info
       Nele = 0;
       Bxel.clear();
       Etel.clear();
       Phiel.clear();
       Etael.clear();
       Isoel.clear();
       
       Njet = 0;
       Bxjet.clear();
       Etjet.clear();
       Phijet.clear();
       Etajet.clear();
       Taujet.clear();
       Fwdjet.clear();  
		 
		 Nmu = 0;
		 Bxmu.clear();
		 Ptmu.clear();
		 Phimu.clear();
		 Etamu.clear();
		 Qualmu.clear();
		 

// ------ ETT, ETM, HTT and HTM from PSB14:

	ETT = -1.;
	OvETT = false;
	HTT = -1.;
	OvHTT = false;
	ETM = -1.;
	PhiETM = -1;
	OvETM = false;
	HTM = -1.;
	PhiHTM = -1;
	OvHTM = false;
    }
   
    // ---- L1AnalysisDataFormat information.
   
    
    //PSB info
    int            Nele;
    std::vector<int>    Bxel;
    std::vector<float>  Etel;
    std::vector<float>  Phiel;
    std::vector<float>  Etael;
    std::vector<bool>   Isoel;
    
    int            Njet;
    std::vector<int>    Bxjet;
    std::vector<float>  Etjet;
    std::vector<float>  Phijet;
    std::vector<float>  Etajet;
    std::vector<bool>   Taujet;
    std::vector<bool>   Fwdjet;

    int            Nmu;
    std::vector<int>    Bxmu;
    std::vector<float>  Ptmu;
    std::vector<float>  Phimu;
    std::vector<float>  Etamu;
    std::vector<int>   Qualmu;

// ------ ETT, ETM, HTT and HTM :

    float ETT;
    bool OvETT;

    float HTT;
    bool OvHTT;

    float ETM;
    int PhiETM;
    bool OvETM;

    float HTM;
    int PhiHTM;
    bool OvHTM;

  }; 
} 
#endif



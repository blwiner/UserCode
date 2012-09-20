runL1Menu2015_45(){
  
   
   gROOT->ProcessLine(".X initL1Analysis.C");
	gROOT->ProcessLine(".L menu/L1Menu2015.C+");

   gSystem->Exec("date"); //beginning time stamp
	
// Define the input parameters
   bool makeThrPlots = false;
	bool useL1Extra   = true;
	int  useMenu      = 10;  //L1 Menu to use
	int  menuThresh   = 0; // Which thresholds to use for the menu
	int  inputSample  = 10;
	int  inputFiles   = 0;	
	float targetLumi  = 200.; //units of e32  
	int  numEvts      =-1; //Number of events to process (-1 == all in file)
	
// Run the job
   RunL1_HFW(makeThrPlots, useL1Extra, useMenu, menuThresh, inputSample, inputFiles, targetLumi, numEvts);		
	
   gSystem->Exec("date"); //ending time stamp
}

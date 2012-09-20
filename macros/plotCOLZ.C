void plotCOLZ(TH2F* hist, Double_t maxVal, Double_t minVal=0., int nLvl=10, float breakTh=-1., int maxBinX=-1, int maxBinY=-1) {

   
   Double_t delta = (maxVal-minVal)/(float)nLvl;
   hist->SetContour(nLvl+1);
	
	Double_t val = minVal;
   for(int ii=0; ii<nLvl+1 ; ii++) {

	  //printf("Setting Level %i to %f\n",ii,val);
     hist->SetContourLevel(ii,val);
	  val += delta;
   }
	
	if(maxBinX>0) hist->GetXaxis()->SetRange(0,maxBinX);
	if(maxBinY>0) hist->GetYaxis()->SetRange(0,maxBinY);
	hist->SetMaximum(maxVal);
	hist->Draw("COLZ");

	TLatex *myTex = new TLatex();
	myTex->SetNDC();
	myTex->SetTextSize(0.03);
	TString txt = "Red: Rate > 100 kHz";
	myTex->SetTextColor(kRed);
	myTex->DrawLatex(0.15,0.73,txt);
	myTex->SetTextColor(kBlack);

	if(breakTh>0.) {	
	   plotBreakLine(hist,breakTh);
	   plotBreakLine(hist,1.0,2);
		

		txt = "Solid line ";
		txt += breakTh;
		txt += " kHz boundary";
		myTex->DrawLatex(0.15,0.7,txt); 
		txt = "Dashed 1 kHz boundary";
		myTex->DrawLatex(0.15,0.67,txt);
		
	}
	
	return;
}	
void plotBreakLine(TH2F* hist, double thres, int lineStyle=1){
   
	TLine* myL = new TLine();
	myL->SetLineStyle(lineStyle);
	myL->SetLineWidth(3);
	
	for(int bx=1; bx < (hist->GetXaxis()->GetLast()+1); bx++){
	   int by=hist->GetYaxis()->GetLast();
		while(by>=1 && hist->GetBinContent(bx,by)<thres) {
		 by--;
		 //printf("Checking bin bx %i by %i th %f\n",bx,by,hist->GetBinContent(bx,by));
		}
		//printf("Making xline at bx %i by %i \n",bx,by);
		if(by>0 && (hist->GetBinContent(bx,by+1)>0 || hist->GetBinContent(bx+1,by+1)<thres) && by< hist->GetYaxis()->GetLast())
        myL->DrawLine(hist->GetXaxis()->GetBinCenter(bx)-hist->GetXaxis()->GetBinWidth(bx)/2.,hist->GetYaxis()->GetBinCenter(by)+hist->GetYaxis()->GetBinWidth(by)/2.,
	 	                hist->GetXaxis()->GetBinCenter(bx)+hist->GetXaxis()->GetBinWidth(bx)/2.,hist->GetYaxis()->GetBinCenter(by)+hist->GetYaxis()->GetBinWidth(by)/2.);  		
	}

	for(int by=1; by < (hist->GetYaxis()->GetLast()+1); by++){
	   int bx=hist->GetXaxis()->GetLast(); //GetNbinsX();
		while(bx>=1 && hist->GetBinContent(bx,by)<thres) {
		 bx--;
		 //printf("Checking bin bx %i by %i th %f\n",bx,by,hist->GetBinContent(bx,by));
		}
		//printf("Making yline at bx %i by %i \n",bx,by);
		if(bx>0 && bx<hist->GetXaxis()->GetLast())
        myL->DrawLine(hist->GetXaxis()->GetBinCenter(bx)+hist->GetXaxis()->GetBinWidth(bx)/2.,hist->GetYaxis()->GetBinCenter(by)-hist->GetYaxis()->GetBinWidth(by)/2.,
	  	                hist->GetXaxis()->GetBinCenter(bx)+hist->GetXaxis()->GetBinWidth(bx)/2.,hist->GetYaxis()->GetBinCenter(by)+hist->GetYaxis()->GetBinWidth(by)/2.);  		
	}

	
	return;
}
void plotAllTH2D(bool save=false, TString code="_Test") {

   double maxVal = 100.;
	double minVal = 0.;
	int    nCont  = 10;
	double breakTh= 10.;
	   
   nslide(1,1);
	gStyle->SetOptStat(0);

// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2_DoubleEG_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
	tmp->SetTitle("Double EG");
	tmp->GetXaxis()->SetTitle("Leading EG E_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("2nd EG E_{T}^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,-1,25);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());

// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2_DoubleMu_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
	tmp->SetTitle("Double Mu");
	tmp->GetXaxis()->SetTitle("Leading Mu P_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("2nd Mu P_{T}^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);	
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,30,20);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());

// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2_Mu_EG_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
   tmp->SetTitle("Mu + EG");
	tmp->GetXaxis()->SetTitle("Mu P_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("EG E_{T}^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);	
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,30,40);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());
	
// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2_SingleMu_CJet_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
   tmp->SetTitle("Single Mu + Central Jet");
	tmp->GetXaxis()->SetTitle("Mu P_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("CJet E_{T}^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);		
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,70,51);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());	

// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2_SingleMu_ETM_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
	tmp->SetTitle("Single Mu + ETM");
	tmp->GetXaxis()->SetTitle("Mu P_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("ETM^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);		
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,70,161);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());

// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2_SingleEG_CJet_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
	tmp->SetTitle("Single EG + Central Jet");
	tmp->GetXaxis()->SetTitle("EG E_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("CJet E_{T}^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);		
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,-1,-1);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());	

// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2_SingleEG_ETM_byThreshold";	
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
   tmp->SetTitle("Single EG + ETM");
	tmp->GetXaxis()->SetTitle("EG E_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("ETM^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);			
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,-1,-1);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());

// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2_DoubleJet_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
	tmp->SetTitle("Double Jet");
	tmp->GetXaxis()->SetTitle("Leading Jet E_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("2nd Jet E_{T}^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);		
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,-1,-1);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());

// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2_DoubleTau_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
	tmp->SetTitle("Double Tau");
	tmp->GetXaxis()->SetTitle("Leading Tau E_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("2nd Tau E_{T}^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);		
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,26,26);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());

// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2A_QuadJetCentral_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
	tmp->SetTitle("Quad Central Jet 1,3");
	tmp->GetXaxis()->SetTitle("Leading Jet E_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("2nd, 3rd, and 4th Jet E_{T}^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);			
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,51,21);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());

// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2B_QuadJetCentral_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
	tmp->SetTitle("Quad Central Jet 2,2");
	tmp->GetXaxis()->SetTitle("Leading & 2nd Jet E_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("3rd & 4th Jet E_{T}^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);			
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,51,21);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());


// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2_SingleCJet_ETM_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
	tmp->SetTitle("Single CenJet + ETM");
	tmp->GetXaxis()->SetTitle(" Jet E_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("ETM^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);			
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,51,-1);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());

// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2_DoubleCJet_ETM_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
	tmp->SetTitle("Double CenJet + ETM");
	tmp->GetXaxis()->SetTitle(" Jets E_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("ETM^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);			
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,51,-1);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());

// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2_HTT_ETM_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
	tmp->SetTitle("HTT + ETM");
	tmp->GetXaxis()->SetTitle(" HTT^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("ETM^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);			
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,-1,-1);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());

// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2_SingleTau_ETM_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
	tmp->SetTitle("Single Tau + ETM");
	tmp->GetXaxis()->SetTitle(" Tau E_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("ETM^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);			
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,51,-1);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());

// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2_SingleTau_CJet_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
	tmp->SetTitle("Single Tau + CenJet");
	tmp->GetXaxis()->SetTitle("Tau E_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("Jet E_{T}^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);			
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,51,-1);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());

// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2_Mu_Tau_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
	tmp->SetTitle("Mu + Tau");
	tmp->GetXaxis()->SetTitle("Mu P_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("Tau E_{T}^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);			
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,101,21);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());

// Draw one of the 2D plots; set the ranges if  needed	
   TString hname = "h2_EG_Tau_byThreshold";
	TH2F* tmp = (TH2F*)(f->Get(hname))->Clone();
	tmp->SetTitle("EG + Tau");
	tmp->GetXaxis()->SetTitle("EG E_{T}^{Th}"); tmp->GetXaxis()->SetTitleColor(kBlack); tmp->GetXaxis()->SetTitleFont(62);
	tmp->GetYaxis()->SetTitle("Tau E_{T}^{Th}"); tmp->GetYaxis()->SetTitleColor(kBlack); tmp->GetYaxis()->SetTitleFont(62);			
	plotCOLZ(tmp,maxVal,minVal,nCont,breakTh,-1,36);
	hname += code;
	hname += ".png";
	if(save) Plots->SaveAs(hname.Data());
		
}

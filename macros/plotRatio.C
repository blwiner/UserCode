plotRatio (TH1F* dataHist, TH1F* mcHist, int nbins=-1, TString lbl1 = "Data", TString lbl2 = "MC", bool logSc = false, float min_tot=-1., float max_tot = -1., float ratioMax = -1., float ratioMin =-1.) {


//   TH1F* dataHist = (TH1F*)f->Get(histName)->Clone();
// 	TH1F* mcHist   = (TH1F*)f2->Get(histName)->Clone();

//mcHist->Scale(1./857000.);
//dataHist->Scale(1./357000.);  

//double ratioMax = 2.3;
//double ratioMin = 0.0;

//Hack to get it plotted with ratio plot
TCanvas* myC = new TCanvas("myC", "myC", 1000,700);
gStyle->SetPadBorderMode(0);
gStyle->SetFrameBorderMode(0);
Float_t small = 1.e-5;
myC->Divide(1,2,small,small);
const float padding=1e-5; const float ydivide=0.3;
myC->GetPad(1)->SetPad( padding, ydivide + padding , 1-padding, 1-padding);
myC->GetPad(2)->SetPad( padding, padding, 1-padding, ydivide-padding);
myC->GetPad(1)->SetRightMargin(.05);
myC->GetPad(2)->SetRightMargin(.05);
myC->GetPad(1)->SetBottomMargin(.4);
myC->GetPad(2)->SetBottomMargin(.4);
myC->GetPad(2)->SetGridy(true);
myC->GetPad(1)->SetGridy(true);
myC->GetPad(1)->SetGridx(true);
myC->GetPad(1)->Modified();
myC->GetPad(2)->Modified();
myC->cd(1);
gPad->SetBottomMargin(small);
gPad->Modified();



TH1F* myRatio = (TH1F*)dataHist->Clone("Ratio");
myRatio->Reset("ICE");
myRatio->SetTitle("");

 dataHist->SetStats(0);
 mcHist->SetStats(0); 
 if(nbins>0) {
   dataHist->GetXaxis()->SetRange(1,nbins);
   myRatio->GetXaxis()->SetRange(1,nbins);
   mcHist->GetXaxis()->SetRange(1,nbins);
 }
 mcHist->SetLineColor(kRed);
 mcHist->SetLineWidth(2.);
 
TLegend *mleg = new TLegend(0.3,0.80,0.6,0.9); 
mleg->AddEntry(dataHist,lbl1,"p");
mleg->AddEntry(mcHist,lbl2,"l");
mleg->SetFillColor(0);
mleg->SetTextSize(0.05);


myRatio->SetStats(0);
myRatio->SetLineColor(kBlue);
myRatio->SetMarkerColor(kBlue);

myRatio->Divide(dataHist,mcHist);
//for(int i=1; i<myRatio->GetNbinsX(); i++) if(mcHist->GetBinContent(i)>0.)myRatio->SetBinContent(i,(dataHist->GetBinContent(i)/mcHist->GetBinContent(i)));

// x-axis parameters
myRatio->GetXaxis()->SetLabelSize(0.12);
myRatio->GetXaxis()->SetLabelOffset(0.04);
myRatio->GetXaxis()->SetTitleOffset(1.);
myRatio->GetXaxis()->SetTitle(dataHist->GetXaxis()->GetTitle()); //make y label bigger
myRatio->GetXaxis()->SetTitleColor(kBlack);
myRatio->GetXaxis()->SetTitleFont(62);
myRatio->GetXaxis()->SetTitleSize(0.12);
//y axis parameters
myRatio->GetYaxis()->SetNdivisions(50000+404);
myRatio->GetYaxis()->SetLabelSize(0.1); //make y label bigger
myRatio->GetXaxis()->SetLabelSize(0.1); //make y label bigger
myRatio->GetYaxis()->SetTitle("Ratio");
myRatio->GetYaxis()->SetTitleSize(0.1);
myRatio->GetYaxis()->SetTitleOffset(.43);

myRatio->SetMarkerStyle(20);
myRatio->SetMarkerSize(1.0);

if(ratioMin<0.) ratioMin = 0.;
if(ratioMax<0.) {
  ratioMax = 1.15*myRatio->GetMaximum();
  if(ratioMax<2.0) ratioMax = 2.0;
}
myRatio->SetMinimum(ratioMin);
myRatio->SetMaximum(ratioMax);

myC->cd(2);
gPad->SetTopMargin(small);
gPad->SetTickx();
gPad->Modified();

myRatio->Draw("p");

if(min_tot==-1.) min_tot = 0.;
if(max_tot==-1.) (dataHist->GetMaximum() > mcHist->GetMaximum()) ? max_tot=dataHist->GetMaximum()*1.15 : max_tot=mcHist->GetMaximum()*1.15;
dataHist->GetYaxis()->SetRangeUser(min_tot,max_tot);

myC->cd(1);
if(logSc) myC_1->SetLogy(true);

//mcHist_1sig->SetFillStyle(3354);
//mcHist_1sig->SetFillColor(kBlack);

dataHist->GetYaxis()->SetTitle("Num/event");
dataHist->GetYaxis()->SetTitle("Rate (kHz)");
dataHist->GetYaxis()->SetTitleSize(0.04);
dataHist->GetYaxis()->SetTitleOffset(1.0);
dataHist->SetMarkerSize(1.0);
dataHist->SetMarkerStyle(20);


dataHist->Draw("pe");
mcHist->Draw("histsame");
//mcHist_1sig->Draw("e2same");
//dataHist->Draw("pe1same");
mleg->Draw();

}
void plotAll(int plotCode = 0xFFF) {

TString labelg = "_Data_JS_0.529_PreSc60";
TString label1 = "5GeV JetSeed";
TString label2 = "No JetSeed";
bool logSc = false;
minValue = 0.0;
TString nameH = "";

  
  
if(plotCode&0x001 > 0) {
  nameH  = "h_Egamma";
  plotRatio((TH1F*)f->Get(nameH),(TH1F*)f2->Get(nameH),12,label1,label2,logSc,minValue);
  nameH += labelg;
  nameH += ".png";
  myC->SaveAs(nameH);
}

  printf("plotCode %i  %i \n",plotCode,(plotCode&0x002));
  
if( (plotCode&0x002) > 0) {
	nameH  = "h_MultiEgamma";
	plotRatio((TH1F*)f->Get(nameH),(TH1F*)f2->Get(nameH),5,label1,label2,logSc,minValue);
	nameH += labelg;
	nameH += ".png";
	myC->SaveAs(nameH);
}

if( (plotCode&0x004) > 0) {
	nameH  = "h_Muons";
	plotRatio((TH1F*)f->Get(nameH),(TH1F*)f2->Get(nameH),13,label1,label2,logSc,minValue);
	nameH += labelg;
	nameH += ".png";
	myC->SaveAs(nameH);
}

if( (plotCode&0x008) > 0) {
	nameH  = "h_MultiMuons";
	plotRatio((TH1F*)f->Get(nameH),(TH1F*)f2->Get(nameH),13,label1,label2,logSc,minValue);
	nameH += labelg;
	nameH += ".png";
	myC->SaveAs(nameH);
}

if( (plotCode&0x010) > 0) {
	nameH  = "h_Jets";
	plotRatio((TH1F*)f->Get(nameH),(TH1F*)f2->Get(nameH),9,label1,label2,logSc,minValue);
	nameH += labelg;
	nameH += ".png";
	myC->SaveAs(nameH);
}

if( (plotCode&0x020) > 0) {
	nameH  = "h_MultiJets";
	plotRatio((TH1F*)f->Get(nameH),(TH1F*)f2->Get(nameH),15,label1,label2,logSc,minValue);
	nameH += labelg;
	nameH += ".png";
	myC->SaveAs(nameH);
}

if( (plotCode&0x040) > 0) {
	nameH  = "h_Sums";
	plotRatio((TH1F*)f->Get(nameH),(TH1F*)f2->Get(nameH),14,label1,label2,logSc,minValue);
	nameH += labelg;
	nameH += ".png";
	myC->SaveAs(nameH);
}

if( (plotCode&0x080) > 0) {
	nameH  = "h_Cross";
	plotRatio((TH1F*)f->Get(nameH),(TH1F*)f2->Get(nameH),22,label1,label2,logSc,minValue);
	nameH += labelg;
	nameH += ".png";
	myC->SaveAs(nameH);
}

if( (plotCode&0x100) > 0) {
	nameH  = "h_MultiCross";
	plotRatio((TH1F*)f->Get(nameH),(TH1F*)f2->Get(nameH),10,label1,label2,logSc,minValue);
	nameH += labelg;
	nameH += ".png";
	myC->SaveAs(nameH);
}

}


void plotTrigThr(TString nameH, bool save=false, double minValue = 1., int nBins = -1, bool logSc=true) {

/*
TString labelg = "_FastvsFullSim";
TString label1 = " 14 TeV FastSim PU 66";
TString label2 = " 14 TeV FullSim PU 50";
*/
TString labelg = "_45v66PU";
TString label1 = "8 TeV Data PU 66";
TString label2 = "8 TeV Data PU 45";

TString labelx = "Threshold (GeV)";
TString histN = "h_";
histN += nameH;
histN += "_byThreshold";


  ((TH1F*)f->Get(histN))->GetXaxis()->SetTitle(labelx);
  ((TH1F*)f->Get(histN))->SetTitle(nameH);
  plotRatio((TH1F*)f->Get(histN),(TH1F*)f2->Get(histN),nBins,label1,label2,logSc,minValue);
  nameH += labelg;
  nameH += ".png";
  if(save) myC->SaveAs(nameH);
  
  

}
void plotTrigPrim(int code) {

  switch(code) {
    case 1:
	   plotTrigThr("h_CJet_Njet",0.0,-1,false);
		break;
    case 2:
	   plotTrigThr("h_CJet_Et",1e-7,-1,true);
		break;
    case 3:
	   plotTrigThr("h_CJet_Eta",0.0,-1,false);
		break; 
    case 4:
	   plotTrigThr("h_CJet_Phi",0.,-1,false);
		break; 

    case 21:
	   plotTrigThr("h_FJet_Njet",0.0,-1,false);
		break;		
    case 22:
	   plotTrigThr("h_FJet_Et",1e-7,-1,true);
		break;		
    case 23:
	   plotTrigThr("h_FJet_Eta",0.,-1,false);
		break; 
    case 24:
	   plotTrigThr("h_FJet_Phi",0.,-1,false);
		break; 

    case 31:
	   plotTrigThr("h_TJet_Njet",0.0,-1,false);
		break;		
    case 32:
	   plotTrigThr("h_TJet_Et",1e-7,-1,true);
		break;		
    case 33:
	   plotTrigThr("h_TJet_Eta",0.,-1,false);
		break; 
    case 34:
	   plotTrigThr("h_TJet_Phi",0.,-1,false);
		break; 

    case 41:
	   plotTrigThr("h_isoEG_Nele",0.0,-1,false);
		break;		
    case 42:
	   plotTrigThr("h_isoEG_Et",1e-7,-1,true);
		break;
    case 43:
	   plotTrigThr("h_isoEG_Eta",0.,-1,false);
		break; 
    case 44:
	   plotTrigThr("h_isoEG_Phi",0.,-1,false);
		break; 

    case 51:
	   plotTrigThr("h_nIsoEG_Nele",0.0,-1,false);
		break;		
    case 52:
	   plotTrigThr("h_nIsoEG_Et",1e-7,-1,true);
		break;
    case 53:
	   plotTrigThr("h_nIsoEG_Eta",0.,-1,false);
		break; 
    case 54:
	   plotTrigThr("h_nIsoEG_Phi",0.,-1,false);
		break;
		 
    case 61:
	   plotTrigThr("h_Mu_Nmu",0.0,-1,false);
		break;		
    case 62:
	   plotTrigThr("h_Mu_Et",1e-5,-1,true);
		break;
    case 63:
	   plotTrigThr("h_Mu_Eta",0.,-1,false);
		break; 
    case 64:
	   plotTrigThr("h_Mu_Phi",0.,-1,false);
		break; 


    case 71:
	   plotTrigThr("h_Sum_ETT",1e-7,1000,true);
		break;
    case 72:
	   plotTrigThr("h_Sum_ETM",1e-7,-1,true);
		break;					
    case 73:
	   plotTrigThr("h_Sum_PhiETM",0,-1,false);
		break;	
    case 74:
	   plotTrigThr("h_Sum_HTT",1e-7,1000,true);
		break;
    case 75:
	   plotTrigThr("h_Sum_HTM",1e-7,-1,true);
		break;					
    case 76:
	   plotTrigThr("h_Sum_PhiHTM",0.,-1,false);
		break;	
	} 
	 	

}

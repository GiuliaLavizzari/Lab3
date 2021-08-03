#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include "TMath.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TColor.h"
#include "TLegend.h"


using namespace std;

int main (int argc, char** argv)
{
	
	TApplication* myApp = new TApplication("myApp", NULL, NULL);
	int nBins = 35500;
	
	double gB =  pow(10, 9/20.);//gain in ADC
	double gC =  pow(10, 9/20.);
	double gL =  pow(10, 5/20.);
		
	// ########## LYSO ######################
	TH1D histoLS ("LYSOwB", "LYSOwB", nBins, -995.5/gL, 353995.0/gL);
	double binLS;
	double contentLS;
		
	int counterLS = 0;
	TString nomeLS = Form("segnaleNa5_lyso_histo.txt");
	ifstream datafileLS (nomeLS.Data());
	while (datafileLS.good()) 
	{
		datafileLS >> binLS >> contentLS;
		histoLS.SetBinContent(counterLS, contentLS);
		counterLS ++;
	}
	
	TH1D histoLF ("Fondo LYSO", "Fondo LYSO", nBins, -995.5/gL, 353995.0/gL);
	double binLF;
	double contentLF;
	int counterLF = 0;
	TString nomeLF = Form ("fondoNa5_lyso_histo.txt");
	ifstream datafileLF(nomeLF.Data());
	while (datafileLF.good()) 
	{
		datafileLF >> binLF >> contentLF;
		histoLF.SetBinContent(counterLF, contentLF);
		counterLF ++;
	}
	
	TH1D histoL ("LYSO con sottrazione del fondo", "LYSO con sottrazione del fondo", nBins, -995.5/gL, 353995.0/gL);
	histoL.Add(&histoLS, &histoLF, 1, -1);
	histoL.Rebin(10);
	histoL.SetTitle("Spettri del Sodio:");
	//histoL.SetFillColor(kMagenta);
	histoL.SetLineColor(kMagenta);
	//histoL.SetFillStyle(3004); 
	histoL.GetXaxis()->SetTitle("ADC channels");
	histoL.GetYaxis()->SetTitle("counts");
	histoL.GetYaxis()->SetRangeUser(0, 1500);
	
	
	
	//########### BGO ################
	TH1D histoBS ("BGO", "BGO", nBins, -995.5/gB, 353995.0/gB);
	double binBS;
	double contentBS;
	int counterBS = 0;
	TString nomeBS = Form("segnaleNa_bgo.txt");
	ifstream datafileBS (nomeBS.Data());
	while (datafileBS.good()) 
	{
		datafileBS >> binBS >> contentBS;
		histoBS.SetBinContent(counterBS, contentBS);
		counterBS ++;
	}
	
	//histoBS.SetFillColor(kBlue);
	histoBS.SetLineColor(kBlue); 
	//histoBS.SetFillStyle(3004);
	histoBS.GetXaxis()->SetTitle("ADC channels");
	histoBS.GetYaxis()->SetTitle("counts");
	histoBS.Rebin(10); 
	
	
	
	
	
	//########### CSI ################
	TH1D histoCS ("CsI", "CsI", nBins, -995.5/gC, 353995.0/gC);
	double binCS;
	double contentCS;
	int counterCS = 0;
	TString nomeCS = Form("segnaleNa_csi.txt");
	ifstream datafileCS (nomeCS.Data());
	while (datafileCS.good()) 
	{
		datafileCS >> binCS >> contentCS;
		histoCS.SetBinContent(counterCS, contentCS);
		counterCS ++;
	}
	//histoCS.SetFillColor(kBlack);
	histoCS.SetLineColor(kBlack); 
	//histoCS.SetFillStyle(3004);
	histoCS.GetXaxis()->SetTitle("ADC");
	histoCS.GetYaxis()->SetTitle("conteggi");
	histoCS.Rebin(10); 
	
	
	
	TCanvas* myC0 = new TCanvas("myC0","myC0");
	
	myC0->cd();
	histoL.Draw("hist, SAME");
	histoBS.Draw("SAME");
	histoCS.Draw("SAME");
	
	TLegend* legend1 = new TLegend (0.1, 0.1, 0.3, 0.3); 
	legend1->AddEntry(&histoL,"LYSO","f");
	legend1->AddEntry(&histoBS,"BGO","f");
	legend1->AddEntry(&histoCS,"CsI", "f");
	legend1->Draw();
	
	myC0->Modified();
	myC0->Update();
	
	
	
	
	TH1D histo2LS ("LYSO2wB", "LYSO2wB", nBins, -995.5/gL, 353995.0/gL);
	double bin2LS;
	double content2LS;
	int counter2LS = 0;
	TString nome2LS = Form("segnaleCo_lyso.txt");
	ifstream datafile2LS (nome2LS.Data());
	while (datafile2LS.good()) 
	{
		datafile2LS >> bin2LS >> content2LS;
		histo2LS.SetBinContent(counter2LS, content2LS);
		counter2LS ++;
	}
	TH1D histo2L ("LYSO2 con sottrazione del fondo", "LYSO2 con sottrazione del fondo", nBins, -995.5/gL, 353995.0/gL);
	histo2L.Add(&histo2LS, &histoLF, 1, -1);
	histo2L.Rebin(10);
	histo2L.SetTitle("Spettri del Cobalto:");
	//histo2L.SetFillColor(kMagenta);
	histo2L.SetLineColor(kMagenta);
	//histo2L.SetFillStyle(3004); 
	histo2L.GetXaxis()->SetTitle("ADC");
	histo2L.GetYaxis()->SetTitle("conteggi");
	
	
	TCanvas* myC1 = new TCanvas("myC1","myC1");
	myC1->cd();
	histo2L.Draw("hist, SAME");
	//histo2BS.Draw("SAME");
	//histo2CS.Draw("SAME");
	
	TLegend* legend2 = new TLegend (0.1, 0.1, 0.3, 0.3); 
	legend2->AddEntry(&histoL,"LYSO","f");
	//legend2->AddEntry(&histoBS,"BGO","f");
	//legend2->AddEntry(&histoCS,"CsI", "f");
	legend2->Draw();
	
	myC1->Modified();
	myC1->Update();
	
	
	/*
	
	TCanvas* myC2 = new TCanvas("myC2","myC2");
	myC2->cd();
	myC2->Modified();
	myC2->Update();
	
	*/
	
	

	
	
	
	myApp->Run();
	return 0;
}

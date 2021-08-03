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
		
	// ########## LYSO ######################
	TH1D histoLS ("LYSO", "LYSO", nBins, -995.5, 353995.0);
	double binLS;
	double contentLS;
	vector<double> vecLYSO;
	vector<double> vecfLYSO;
	vector<double> vecBGO;
	vector<double> vecfBGO;
	vector<double> vecCSI;
	vector<double> vecfCSI;
	
	int counterLS = 0;
	TString nomeLS = Form("segnaleNa5_lyso_histo.txt");
	ifstream datafileLS (nomeLS.Data());
	while (datafileLS.good()) 
	{
		datafileLS >> binLS >> contentLS;
		histoLS.SetBinContent(counterLS, contentLS);
		vecLYSO.push_back(contentLS);
		counterLS ++;
	}
	
	TH1D histoLF ("Fondo LYSO", "Fondo LYSO", nBins, -995.5, 353995.0);
	double binLF;
	double contentLF;
	int counterLF = 0;
	TString nomeLF = Form ("fondoNa5_lyso_histo.txt");
	ifstream datafileLF(nomeLF.Data());
	while (datafileLF.good()) 
	{
		datafileLF >> binLF >> contentLF;
		histoLF.SetBinContent(counterLF, contentLF);
		vecfLYSO.push_back(contentLF);
		counterLF ++;
	}
	
	double sumL = 0;
	double sumLf = 0;
	for (int i = 3101; i < 4801; i++) 
		{
		sumL = sumL + vecLYSO.at(i);
		sumLf = sumLf + vecfLYSO.at(i);
		}
	cout << "sumL: " << sumL << endl;
	cout << "sumLf: " << sumLf << endl;
	
	
	//########### BGO ################
	TH1D histoBS ("BGO", "BGO", nBins, -995.5, 353995.0);
	double binBS;
	double contentBS;
	int counterBS = 0;
	TString nomeBS = Form("segnaleNa_bgo.txt");
	ifstream datafileBS (nomeBS.Data());
	while (datafileBS.good()) 
	{
		datafileBS >> binBS >> contentBS;
		histoBS.SetBinContent(counterBS, contentBS);
		vecBGO.push_back(contentBS);
		counterBS ++;
	}
	
	TH1D histoBF ("Fondo BGO", "Fondo BGO", nBins, -995.5, 353995.0);
	double binBF;
	double contentBF;
	int counterBF = 0;
	//TString nomeBF = Form ("fondo_bgo_histo.txt");
	TString nomeBF = Form ("fondoNa_bgo_histo.txt");
	ifstream datafileBF(nomeBF.Data());
	while (datafileBF.good()) 
	{
		datafileBF >> binBF >> contentBF;
		histoBF.SetBinContent(counterBF, contentBF);
		vecfBGO.push_back(contentBF);
		counterBF ++;
	}
	
	double sumB = 0;
	double sumBf = 0;
	for (int i = 1001; i < 2101; i++) 
		{
		sumB = sumB + vecBGO.at(i);
		sumBf = sumBf + vecfBGO.at(i);
		}
	cout << "sumB: " << sumB << endl;
	cout << "sumBf: " << sumBf << endl;
	
	//########### CSI ################
	TH1D histoCS ("CsI", "CsI", nBins, -995.5, 353995.0);
	double binCS;
	double contentCS;
	int counterCS = 0;
	TString nomeCS = Form("segnaleNa_csi.txt");
	ifstream datafileCS (nomeCS.Data());
	while (datafileCS.good()) 
	{
		datafileCS >> binCS >> contentCS;
		histoCS.SetBinContent(counterCS, contentCS);
		vecCSI.push_back(contentCS);
		counterCS ++;
	}
	
	TH1D histoCF ("Fondo CsI", "Fondo CsI", nBins, -995.5, 353995.0);
	double binCF;
	double contentCF;
	int counterCF = 0;
	//TString nomeCF = Form ("fondo_csi_histo.txt"); 
	TString nomeCF = Form ("fondoNa_csi_histo.txt"); 
	ifstream datafileCF(nomeCF.Data());
	while (datafileCF.good()) 
	{
		datafileCF >> binCF >> contentCF;
		histoCF.SetBinContent(counterCF, contentCF);
		vecfCSI.push_back(contentCF);
		counterCF ++;
	}
	
	double sumC = 0;
	double sumCf = 0;
	for (int i = 71020; i < 10601; i++) 
		{
		sumC = sumC + vecCSI.at(i);
		sumCf = sumCf + vecfCSI.at(i);
		}
	cout << "sumC: " << sumC << endl;
	cout << "sumCf: " << sumCf << endl;
	
	
	TCanvas* myC0 = new TCanvas("myC0","myC0");
	
	myC0->cd();
	
	//histoLS.SetFillColor(kMagenta);
	histoLS.SetLineColor(kMagenta); //kOrange+5
	//histoLS.SetFillStyle(3004);
	histoLS.GetXaxis()->SetTitle("ADC");
	histoLS.GetYaxis()->SetTitle("conteggi");
	histoLS.SetTitle("LYSO: fondo e segnale sodio");
	histoLS.Rebin(20); //rebin = merge degi bin a dieci a dieci
	histoLS.Draw("hist, SAME");
	
	histoLF.SetFillColor(kMagenta+3);
	histoLF.SetLineColor(kMagenta+3); //kCyan+3
	histoLF.SetFillStyle(3004); //3004
	histoLF.GetXaxis()->SetTitle("ADC");
	histoLF.GetYaxis()->SetTitle("contggi");
	histoLF.Rebin(20); //rebin = merge degi bin a dieci a dieci
	histoLF.Draw("hist, SAME");
	
	TLegend* legend1 = new TLegend (0.1, 0.1, 0.3, 0.3); 
	legend1->AddEntry(&histoLS,"segnale","f");
	legend1->AddEntry(&histoLF,"fondo","f");
	legend1->Draw();
	
	myC0->Modified();
	myC0->Update();
	
		
		
		
		
	TCanvas* myC1 = new TCanvas("myC1","myC1");
	myC1->cd();
	
	//histoBS.SetFillColor(kBlue);
	histoBS.SetLineColor(kBlue); //kOrange+5
	//histoBS.SetFillStyle(3004);
	histoBS.GetXaxis()->SetTitle("ADC");
	histoBS.GetYaxis()->SetTitle("conteggi");
	histoBS.Rebin(5); //rebin = merge degi bin a dieci a dieci
	histoBS.SetTitle("BGO: fondo e segnale sodio");
	histoBS.Draw("SAME");
	
	histoBF.SetFillColor(kBlue+2);
	histoBF.SetLineColor(kBlue+2); //kCyan+3
	histoBF.SetFillStyle(3004); //3004
	histoBF.GetXaxis()->SetTitle("ADC");
	histoBF.GetYaxis()->SetTitle("conteggi");
	histoBF.Rebin(10); //rebin = merge degi bin a dieci a dieci
	histoBF.Draw("hist, SAME");
	
	TLegend* legend2 = new TLegend (0.1, 0.1, 0.3, 0.3); 
	legend2->AddEntry(&histoBS,"segnale","f");
	legend2->AddEntry(&histoBF,"fondo","f");
	legend2->Draw();
	
	myC1->Modified();
	myC1->Update();
	
	
	
	
	TCanvas* myC2 = new TCanvas("myC2","myC2");
	myC2->cd();
	
	//histoCS.SetFillColor(kBlack);
	histoCS.SetLineColor(kBlack); //kOrange+5
	//histoCS.SetFillStyle(3004);
	histoCS.GetXaxis()->SetTitle("ADC");
	histoCS.GetYaxis()->SetTitle("conteggi");
	histoCS.Rebin(20); 
	histoCS.SetTitle("CsI: fondo e segnale sodio");
	histoCS.Draw("SAME");
	
	histoCF.SetFillColor(kGray+2);
	histoCF.SetLineColor(kGray+2); //kCyan+3
	histoCF.SetFillStyle(3004); //3004
	histoCF.GetXaxis()->SetTitle("ADC");
	histoCF.GetYaxis()->SetTitle("conteggi");
	histoCF.Rebin(40); //rebin = merge degi bin a dieci a dieci
	histoCF.Draw("SAME");
	
	TLegend* legend3 = new TLegend (0.1, 0.1, 0.3, 0.3); 
	legend3->AddEntry(&histoCS,"segnale","f");
	legend3->AddEntry(&histoCF,"fondo","f");
	legend3->Draw();
	
	myC2->Modified();
	myC2->Update();
	
	
	
	
	TCanvas* myC3 = new TCanvas("myC3","myC3");
	myC3->cd();
	
	histoLF.SetName("Fondo LYSO");
	histoLF.SetFillColor(kMagenta+3);
	histoLF.SetLineColor(kMagenta+3); //kOrange+5
	histoLF.SetFillStyle(3004);
	histoLF.Draw("hist");
	
	myC3->Modified();
	myC3->Update();
	
	
	TCanvas* myC4 = new TCanvas("myC4","myC4");	
	myC4->cd();
	
	histoBF.SetName("Fondo BGO");
	histoBF.SetFillColor(kBlue+2);
	histoBF.SetLineColor(kBlue+2); //kOrange+5
	histoBF.SetFillStyle(3004);
	histoBF.Draw("hist");
	
	myC4->Modified();
	myC4->Update();
	
	
	TCanvas* myC5 = new TCanvas("myC5","myC5");
	myC5->cd();
	
	histoCF.SetName("Fondo CsI");
	histoCF.SetFillColor(kGray+2);
	histoCF.SetLineColor(kGray+2); //kOrange+5
	histoCF.SetFillStyle(3004);
	histoCF.Draw("hist");
	
	myC5->Modified();
	myC5->Update();
	
	
	
	
	myApp->Run();
	return 0;
}

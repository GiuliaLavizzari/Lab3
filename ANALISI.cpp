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


double gaussiana (double* x, double* par)
{
	return par[0]*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]));
}

double polinomio (double* x, double* par)
{
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

double somma (double* x, double* par)
{
	return gaussiana(x, par) + polinomio(x, &par[3]);
}

double weightedmean (std::vector<double> val, std::vector<double> err)
{
	double num = 0;
	double den = 0;
	for (int i = 0; i < val.size(); i++)
	{
		num = num + val.at(i)/pow(err.at(i),2);
		den = den + 1/pow(err.at(i), 2);
	}
	return num/den;
}

double weightederror (std::vector<double> err)
{
	double den = 0;
	for (int i = 0; i < err.size(); i++) den = den + 1/pow(err.at(i), 2);
	return sqrt(1/den);
}



using namespace std;

int main (int argc, char** argv)
{
	
	TApplication* myApp = new TApplication("myApp", NULL, NULL);
	gStyle->SetOptFit(1112);
	
	int nBins = 35500;
	int nBins2 = 25500;
	
	double gB =  pow(10, 9/20.);//gain in ADC
	double gC =  pow(10, 9/20.);
	double gL =  pow(10, 5/20.);
		
	// ########## LYSO ######################
	TH1D histoLS ("LYSO_Na", "LYSO_Na", nBins, -995.5, 353995.0);
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
		counterLF ++;
	}
	
	TH1D histoL ("LYSO_BS", "LYSO_BS", nBins, -995.5, 353995.0);
	histoL.Add(&histoLS, &histoLF, 1, -1);
	//histoL.Rebin(10);
	histoL.SetTitle("Sodio con LYSO:");
	histoL.GetXaxis()->SetTitle("ADC channels");
	histoL.GetYaxis()->SetTitle("counts");

	
	TH1D CohistoLS ("LYSOCo", "LYSOCo", nBins, -995.5, 353995.0);
	double binLSCo;
	double contentLSCo;
	int counterLSCo = 0;
	TString nomeLSCo = Form ("segnaleCo_lyso.txt");
	ifstream datafileLSCo(nomeLSCo.Data());
	while (datafileLSCo.good()) 
	{
		datafileLSCo >> binLSCo >> contentLSCo;
		CohistoLS.SetBinContent(counterLSCo, contentLSCo);
		counterLSCo ++;
	}
	
	TH1D CohistoL ("LYSO_Co", "LYSO_Co", nBins, -995.5, 353995.0);
	CohistoL.Add(&CohistoLS, &histoLF, 1, -1);
	//CohistoL.Rebin(10);
	CohistoL.SetTitle("Cobalto con LYSO:");
	CohistoL.GetXaxis()->SetTitle("ADC channels");
	CohistoL.GetYaxis()->SetTitle("counts");
	
	
	
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
		counterBS ++;
	}
	
	histoBS.GetXaxis()->SetTitle("ADC channels");
	histoBS.GetYaxis()->SetTitle("counts");
	histoBS.SetTitle("Sodio con BGO:");
	//histoBS.Rebin(10); 
	
	TH1D histoBSCo ("BGOCos", "BGOCos", nBins2, -995.5, 381492.500);
	double binBSCo;
	double contentBSCo;
	int counterBSCo = 0;
	TString nomeBSCo = Form("segnaleCo9_bgo_histo.txt");
	ifstream datafileBSCo (nomeBSCo.Data());
	while (datafileBSCo.good()) 
	{
		datafileBSCo >> binBSCo >> contentBSCo;
		histoBSCo.SetBinContent(counterBSCo, contentBSCo);
		counterBSCo ++;
	}
	
	histoBSCo.GetXaxis()->SetTitle("ADC channels");
	histoBSCo.GetYaxis()->SetTitle("counts");
	histoBSCo.SetTitle("Cobalto con BGO:");
	
	TH1D histoBFCo ("BGOCof", "BGOCof", nBins2, -995.5, 381492.500);
	double binBFCo;
	double contentBFCo;
	int counterBFCo = 0;
	TString nomeBFCo = Form("fondoCo9_bgo_histo.txt");
	ifstream datafileBFCo (nomeBFCo.Data());
	while (datafileBFCo.good()) 
	{
		datafileBFCo >> binBFCo >> contentBFCo;
		histoBFCo.SetBinContent(counterBFCo, contentBFCo);
		counterBFCo ++;
	}
	
	histoBFCo.GetXaxis()->SetTitle("ADC channels");
	histoBFCo.GetYaxis()->SetTitle("counts");
	histoBFCo.SetTitle("Cobalto con BGO:");
	
	TH1D CohistoB ("BGO_Co", "BGO_Co", nBins2, -995.5, 381492.500);
	CohistoB.Add(&histoBSCo, &histoBFCo, 1, -1);
	CohistoB.Rebin(20);
	CohistoB.SetTitle("Cobalto con BGO:");
	CohistoB.GetXaxis()->SetTitle("ADC channels");
	CohistoB.GetYaxis()->SetTitle("counts");
	//CohistoB.Rebin(1);
	
	
	
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
		counterCS ++;
	}
	histoCS.GetXaxis()->SetTitle("ADC");
	histoCS.GetYaxis()->SetTitle("conteggi");
	histoCS.SetTitle("Sodio con CSI:");
	//histoCS.Rebin(10); 
	
	TH1D histoCSCo ("CSICos", "CSICos", nBins2, -995.5, 381492.500);
	double binCSCo;
	double contentCSCo;
	int counterCSCo = 0;
	TString nomeCSCo = Form("segnaleCo9_csi_histo.txt");
	ifstream datafileCSCo (nomeCSCo.Data());
	while (datafileCSCo.good()) 
	{
		datafileCSCo >> binCSCo >> contentCSCo;
		histoCSCo.SetBinContent(counterCSCo, contentCSCo);
		counterCSCo ++;
	}
	
	histoCSCo.GetXaxis()->SetTitle("ADC channels");
	histoCSCo.GetYaxis()->SetTitle("counts");
	histoCSCo.SetTitle("Cobalto con CSI:");
	
	TH1D histoCFCo ("CSICof", "CSICof", nBins2, -995.5, 381492.500);
	double binCFCo;
	double contentCFCo;
	int counterCFCo = 0;
	TString nomeCFCo = Form("fondoCo9_csi_histo.txt");
	ifstream datafileCFCo (nomeCFCo.Data());
	while (datafileCFCo.good()) 
	{
		datafileCFCo >> binCFCo >> contentCFCo;
		histoCFCo.SetBinContent(counterCFCo, contentCFCo);
		counterCFCo ++;
	}
	
	histoCFCo.GetXaxis()->SetTitle("ADC channels");
	histoCFCo.GetYaxis()->SetTitle("counts");
	histoCFCo.SetTitle("Cobalto con CSI:");
	
	TH1D CohistoC ("CSI_Co", "CSI_Co", nBins2, -995.5, 381492.500);
	CohistoC.Add(&histoCSCo, &histoCFCo, 1, -1);
	CohistoC.Rebin(20);
	CohistoC.SetTitle("Cobalto con CSI:"); 
	CohistoC.GetXaxis()->SetTitle("ADC channels");
	CohistoC.GetYaxis()->SetTitle("counts");
	
	//CohistoC.Rebin(1);
	
	
	//########################################################################
	//############################ 22Na #####################################
	//########################################################################
	
	//############################## LYSO ####################################
	
	TF1* gausL1 = new TF1 ("gausL1", gaussiana, 20000, 51000, 3);
	gausL1->SetLineColor(kRed-7);
	gausL1->SetLineStyle(5);
	gausL1->SetLineWidth(2.);
	
	TF1* polL1 = new TF1 ("polL1", polinomio, 20000, 50000, 3);
	polL1->SetLineColor(kRed);
	
	TF1* sumL1 = new TF1 ("sumL1", somma, 20000, 50000, 6);
	sumL1->SetParameter(0, 5.35116e+02);
	sumL1->SetParName(0, "A");
	sumL1->SetParameter(1, 3.99569e+04);
	sumL1->SetParName(1, "Mean");
	sumL1->SetParameter(2, 5000);
	sumL1->SetParName(2, "Sigma");
	sumL1->SetParameter(3, 0.2*1.83002e+02);
	sumL1->SetParName(3, "p0");
	sumL1->SetParameter(4, 0.2*(-3.54445e-03));
	sumL1->SetParName(4, "p1");
	sumL1->SetParameter(5, 0.2*1.75785e-08);
	sumL1->SetParName(5, "p2");
	sumL1->SetLineColor(kPink-3);
	sumL1->SetLineWidth(2.);
	//sumL1->SetLineStyle(5);
	
	
	
	TF1* gausL2 = new TF1 ("gausL2", gaussiana, 60700, 70000, 3);
	gausL2->SetLineColor(kAzure+10);
	gausL2->SetLineStyle(5);
	gausL2->SetLineWidth(2.);
	
	TF1* polL2 = new TF1 ("polL2", polinomio, 61000, 70000, 3);
	polL2->SetLineColor(kPink-2);
	
	TF1* sumL2 = new TF1 ("sumL2", somma, 61000, 70000, 6);
	sumL2->SetParameter(0, 200);
	sumL2->SetParName(0, "A");
	sumL2->SetParameter(1, 64000);
	sumL2->SetParName(1, "Mean");
	sumL2->SetParameter(2, 5000);
	sumL2->SetParName(2, "Sigma");
	sumL2->SetParameter(3, 0.2*1.83002e+02);
	sumL2->SetParName(3, "p0");
	sumL2->SetParameter(4, 0.2*(3.54445e-03));
	sumL2->SetParName(4, "p1");
	sumL2->SetParameter(5, 0.2*1.75785e-08);
	sumL2->SetParName(5, "p2");
	sumL2->SetLineColor(kAzure+8);
	sumL2->SetLineWidth(2.);
	//sumL2->SetLineStyle(5);
	
	TLegend* legend2 = new TLegend (0.1, 0.1, 0.3, 0.3); 
	legend2->AddEntry(sumL1,"gaus+pol 511 keV","l");
	legend2->AddEntry(gausL1,"gaussiana 511 keV","l");
	legend2->AddEntry(sumL2,"gaus+pol 1274 keV", "l");
	legend2->AddEntry(gausL2, "gaussiana 1274 keV", "l");
	
	TCanvas* myC0 = new TCanvas("myC0","myC0");
	myC0->cd();
	histoL.Rebin(10);
	histoL.Fit("sumL1", "R+");
	gausL1-> SetParameter(0, sumL1->GetParameter(0));
	gausL1-> SetParameter(1, sumL1->GetParameter(1));
	gausL1-> SetParameter(2, sumL1->GetParameter(2));
	gausL1->Draw("same");
	cout << "chi: " << sumL1->GetChisquare()/sumL1->GetNDF() << endl;
	histoL.Fit("sumL2", "R+");
	gausL2-> SetParameter(0, sumL2->GetParameter(0));
	gausL2-> SetParameter(1, sumL2->GetParameter(1));
	gausL2-> SetParameter(2, sumL2->GetParameter(2));
	gausL2->Draw("same");
	legend2->Draw();
	cout << "\nchi: " << sumL2->GetChisquare()/sumL2->GetNDF() << endl;
	myC0->Modified();
	myC0->Update();
	
	
	//########################### BGO ################################
	
	TF1* gausB1 = new TF1 ("gausB1", gaussiana, 7000, 22000, 3);
	gausB1->SetLineColor(kOrange);
	
	TF1* polB1 = new TF1 ("polB1", polinomio, 7000, 22000, 3);
	polB1->SetLineColor(kRed);
	
	TF1* sumB1 = new TF1 ("sumB1", somma, 7000, 22000, 6);
	sumB1->SetParameter(0, 3.75984e+02);
	sumB1->SetParName(0, "A");
	sumB1->SetParameter(1, 15000);
	sumB1->SetParName(1, "Mean");
	sumB1->SetParameter(2, 2000);
	sumB1->SetParName(2, "Sigma");
	sumB1->SetParameter(3, 2.29553e+02);
	sumB1->SetParName(3, "p0");
	sumB1->SetParameter(4, -1.60452e-02);
	sumB1->SetParName(4, "p1");
	sumB1->SetParameter(5, 2.85057e-07);
	sumB1->SetParName(5, "p2");
	sumB1->SetLineColor(kGreen);

	TF1* gausB2 = new TF1 ("gausB2", gaussiana, 30000, 45000, 3);
	
	TF1* polB2 = new TF1 ("polB2", polinomio, 30000, 45000, 3);
	
	TF1* sumB2 = new TF1 ("sumB2", somma, 30000, 45000, 6);
	sumB2->SetParameter(0, 100);
	sumB2->SetParName(0, "A");
	sumB2->SetParameter(1, 37000);
	sumB2->SetParName(1, "Mean");
	sumB2->SetParameter(2, 5000);
	sumB2->SetParName(2, "Sigma");
	sumB2->SetParameter(3, 2.29553e+02);
	sumB2->SetParName(3, "p0");
	sumB2->SetParameter(4, 1.60452e-02);
	sumB2->SetParName(4, "p1");
	sumB2->SetParameter(5, 2.85057e-07);
	sumB2->SetParName(5, "p2");
	sumB2->SetLineColor(kMagenta);

	TCanvas* myC1 = new TCanvas("myC1","myC1");
	myC1->cd();
	histoBS.Rebin(10);
	histoBS.Fit("sumB1", "R+");
	cout << "chi: " << sumB1->GetChisquare()/sumB1->GetNDF() << endl;
	histoBS.Fit("sumB2", "R+");
	cout << "chi: " << sumB2->GetChisquare()/sumB2->GetNDF() << endl;
	myC1->Modified();
	myC1->Update();
	
	//########################### CsI #################################
	
	
	TF1* gausC1 = new TF1 ("gausC1", gaussiana, 40000, 110000, 3);
	gausC1->SetLineColor(kOrange);
	
	TF1* polC1 = new TF1 ("polC1", polinomio, 40000, 110000, 3);
	polC1->SetLineColor(kRed);
	
	TF1* sumC1 = new TF1 ("sumC1", somma, 40000, 110000, 6);
	sumC1->SetParameter(0, 20);
	sumC1->SetParName(0, "A");
	sumC1->SetParameter(1, 85000);
	sumC1->SetParName(1, "Mean");
	sumC1->SetParameter(2, 10000);
	sumC1->SetParName(2, "Sigma");
	sumC1->SetParameter(3, 0.2*1.83002e+02);
	sumC1->SetParName(3, "p0");
	sumC1->SetParameter(4, 0.2*(-3.54445e-03));
	sumC1->SetParName(4, "p1");
	sumC1->SetParameter(5, 0.2*1.75785e-08);
	sumC1->SetParName(5, "p2");
	sumC1->SetLineColor(kGreen);
	
	TF1* gausC2 = new TF1 ("gausC2", gaussiana, 160000, 220000, 3);
	
	TF1* polC2 = new TF1 ("polC2", polinomio, 160000, 220000, 3);
	
	TF1* sumC2 = new TF1 ("sumC2", somma, 160000, 220000, 6);
	sumC2->SetParameter(0, 14);
	sumC2->SetParName(0, "A");
	sumC2->SetParameter(1, 184000);
	sumC2->SetParName(1, "Mean");
	sumC2->SetParameter(2, 5000);
	sumC2->SetParName(2, "Sigma");
	sumC2->SetParameter(3, 0.2*1.83002e+02);
	sumC2->SetParName(3, "p0");
	sumC2->SetParameter(4, 0.2*(3.54445e-03));
	sumC2->SetParName(4, "p1");
	sumC2->SetParameter(5, 0.2*1.75785e-08);
	sumC2->SetParName(5, "p2");
	sumC2->SetLineColor(kMagenta);
	
	TCanvas* myC2 = new TCanvas("myC2","myC2");
	myC2->cd();
	histoCS.Rebin(10); 
	histoCS.Fit("sumC1", "R+");
	cout << "chi: " << sumC1->GetChisquare()/sumC1->GetNDF() << endl;
	histoCS.Fit("sumC2", "R+");
	cout << "chi: " << sumC2->GetChisquare()/sumC2->GetNDF() << endl;
	myC2->Modified();
	myC2->Update();
	
	
	//########################################################################
	//############################ 57Co #####################################
	//########################################################################
	
	//############################## LYSO ####################################
	
	TF1* gausL3 = new TF1 ("gausL3", gaussiana, 5500, 17000, 3);
	gausL3->SetLineColor(kOrange-3);
	gausL3->SetLineStyle(5);
	
	TF1* polL3 = new TF1 ("polL3", polinomio, 5500, 17000, 3);
	polL3->SetLineColor(kRed);
	
	TF1* sumL3 = new TF1 ("sumL3", somma, 5500, 17000, 6);
	sumL3->SetParameter(0, 300);
	sumL3->SetParName(0, "A");
	sumL3->SetParameter(1, 9900);
	sumL3->SetParName(1, "Mean");
	sumL3->SetParameter(2, 2000);
	sumL3->SetParName(2, "Sigma");
	sumL3->SetParameter(3, 0.2*1.83002e+02);
	sumL3->SetParName(3, "p0");
	sumL3->SetParameter(4, 0.2*(3.54445e-03));
	sumL3->SetParName(4, "p1");
	sumL3->SetParameter(5, 0.2*1.75785e-08);
	sumL3->SetParName(5, "p2");
	sumL3->SetLineColor(kOrange+7);
	
	TLegend* legend3 = new TLegend (0.1, 0.1, 0.3, 0.3); 
	legend3->AddEntry(sumL3,"gaus+pol 122 keV","l");
	legend3->AddEntry(gausL3,"gaussiana 122 keV","l");
	
	TCanvas* myC3 = new TCanvas("myC3","myC3");
	myC3->cd();
	CohistoL.Rebin(25);
	CohistoL.Fit("sumL3", "R+");
	gausL3-> SetParameter(0, sumL3->GetParameter(0));
	gausL3-> SetParameter(1, sumL3->GetParameter(1));
	gausL3-> SetParameter(2, sumL3->GetParameter(2));
	gausL3->Draw("same");
	legend3->Draw();
	cout << "chi: " << sumL3->GetChisquare()/sumL3->GetNDF() << endl;
	myC3->Modified();
	myC3->Update();
	
	
	
	
	//############################## BGO ####################################
	
	TF1* gausB3 = new TF1 ("gausB3", gaussiana, 2000, 12000, 3);
	gausB3->SetLineColor(kOrange);
	
	TF1* polB3 = new TF1 ("polB3", polinomio, 2000, 12000, 3);
	polB3->SetLineColor(kRed);
	
	TF1* sumB3 = new TF1 ("sumB3", somma, 2000, 12000, 6);
	sumB3->SetParameter(0, 600);
	sumB3->SetParName(0, "A");
	sumB3->SetParameter(1, 6000);
	sumB3->SetParName(1, "Mean");
	sumB3->SetParameter(2, 2000);
	sumB3->SetParName(2, "Sigma");
	sumB3->SetParameter(3, 0.2*1.83002e+02);
	sumB3->SetParName(3, "p0");
	sumB3->SetParameter(4, 0.2*(3.54445e-03));
	sumB3->SetParName(4, "p1");
	sumB3->SetParameter(5, 0.2*1.75785e-08);
	sumB3->SetParName(5, "p2");
	sumB3->SetLineColor(kGreen);
	
	TCanvas* myC4 = new TCanvas("myC4","myC4");
	myC4->cd();
	CohistoB.Rebin(1);
	CohistoB.Fit("sumB3", "R+");
	cout << "chi: " << sumB3->GetChisquare()/sumB3->GetNDF() << endl;
	myC4->Modified();
	myC4->Update();
	
	
	//############################## CsI ####################################
	TF1* gausC3 = new TF1 ("gausC3", gaussiana, 20000, 60000, 3);
	gausC3->SetLineColor(kOrange);
	
	TF1* polC3 = new TF1 ("polC3", polinomio, 20000, 60000, 3);
	polC3->SetLineColor(kRed);
	
	TF1* sumC3 = new TF1 ("sumC3", somma, 20000, 60000, 6);
	sumC3->SetParameter(0, 140);
	sumC3->SetParName(0, "A");
	sumC3->SetParameter(1, 35000);
	sumC3->SetParName(1, "Mean");
	sumC3->SetParameter(2, 10000);
	sumC3->SetParName(2, "Sigma");
	sumC3->SetParameter(3, 0.2*1.83002e+02);
	sumC3->SetParName(3, "p0");
	sumC3->SetParameter(4, 0.2*(3.54445e-03));
	sumC3->SetParName(4, "p1");
	sumC3->SetParameter(5, 0.2*1.75785e-08);
	sumC3->SetParName(5, "p2");
	sumC3->SetLineColor(kGreen);
	
	
	TCanvas* myC5 = new TCanvas("myC5","myC5");
	myC5->cd();
	CohistoC.Rebin(1);
	CohistoC.Fit("sumC3", "R+");
	cout << "chi: " << sumC3->GetChisquare()/sumC3->GetNDF() << endl;
	myC5->Modified();
	myC5->Update();
	
	
	//#####################################################################
	vector<double> veclyso; // [0]=511keV [1]=1274keV [2]=122keV
	vector<double> vecdevstL;
	vector<double> vecerrlyso; // calcolato come FWHM media/gain*2*sqrt(2*log(2))
	vector<double> errFWHMlyso; //serve per la risoluzione
	veclyso.push_back(sumL1->GetParameter(1)/gL);
	veclyso.push_back(sumL2->GetParameter(1)/gL);
	veclyso.push_back(sumL3->GetParameter(1)/gL);
	vecdevstL.push_back(sumL1->GetParameter(2)/gL);
	vecdevstL.push_back(sumL2->GetParameter(2)/gL);
	vecdevstL.push_back(sumL3->GetParameter(2)/gL);
	errFWHMlyso.push_back((sumL1->GetParError(2)/gL)*2*sqrt(2*log(2)));
	errFWHMlyso.push_back((sumL2->GetParError(2)/gL)*2*sqrt(2*log(2)));
	errFWHMlyso.push_back((sumL3->GetParError(2)/gL)*2*sqrt(2*log(2)));
	for (int i = 0; i < veclyso.size(); i++) vecerrlyso.push_back(abs(vecdevstL.at(i)*2*sqrt(2*log(2))));
	for (int i = 0; i < veclyso.size(); i++) cout << veclyso.at(i) << " +- " << vecerrlyso.at(i) << endl;
	
	vector<double> vecbgo;
	vector<double> vecdevstB;
	vector<double> vecerrbgo;
	vector<double> errFWHMbgo;
	vecbgo.push_back(sumB1->GetParameter(1)/gB);
	vecbgo.push_back(sumB2->GetParameter(1)/gB);
	vecbgo.push_back(sumB3->GetParameter(1)/gB);
	vecdevstB.push_back(sumB1->GetParameter(2)/gB);
	vecdevstB.push_back(sumB2->GetParameter(2)/gB);
	vecdevstB.push_back(sumB3->GetParameter(2)/gB);
	errFWHMbgo.push_back((sumB1->GetParError(2)/gB)*2*sqrt(2*log(2)));
	errFWHMbgo.push_back((sumB2->GetParError(2)/gB)*2*sqrt(2*log(2)));
	errFWHMbgo.push_back((sumB3->GetParError(2)/gB)*2*sqrt(2*log(2)));
	for (int i = 0; i < vecbgo.size(); i++) vecerrbgo.push_back(abs(vecdevstB.at(i)*2*sqrt(2*log(2))));
	for (int i = 0; i < vecbgo.size(); i++) cout << vecbgo.at(i) << " +- " << vecerrbgo.at(i) << endl;
	vector<double> veccsi;
	vector<double> vecdevstC;
	vector<double> vecerrcsi;
	vector<double> errFWHMcsi;
	veccsi.push_back(sumC1->GetParameter(1)/gC);
	veccsi.push_back(sumC2->GetParameter(1)/gC);
	veccsi.push_back(sumC3->GetParameter(1)/gC);
	vecdevstC.push_back(sumC1->GetParameter(2)/gC);
	vecdevstC.push_back(sumC2->GetParameter(2)/gC);
	vecdevstC.push_back(sumC3->GetParameter(2)/gC);
	errFWHMcsi.push_back((sumC1->GetParError(2)/gC)*2*sqrt(2*log(2)));
	errFWHMcsi.push_back((sumC2->GetParError(2)/gC)*2*sqrt(2*log(2)));
	errFWHMcsi.push_back((sumC3->GetParError(2)/gC)*2*sqrt(2*log(2)));
	for (int i = 0; i < veccsi.size(); i++) vecerrcsi.push_back(abs(vecdevstC.at(i)*2*sqrt(2*log(2))));
	for (int i = 0; i < veccsi.size(); i++) cout << veccsi.at(i) << " +- " << vecerrcsi.at(i) << endl;
	
	
	
	
	
	//#########################################################################
	//############################# LIGHT YIELD ###############################
	//#########################################################################
	vector<double> ly1;
	vector<double> errly1;
	vector<double> ly2;
	vector<double> errly2;
	vector<double> ly3; 
	vector<double> errly3;
	
	double A, B;
	double errA, errB;
	
	cout << "\n\n\nlight yield bgo/csi (expected = 0.16):\n" << endl;
	for (int i = 0; i < veclyso.size(); i++) 
	{
		A = vecbgo.at(i);
		B = veccsi.at(i);
		errA = vecerrbgo.at(i);
		errB = vecerrcsi.at(i);
		ly1.push_back(A/B);
		errly1.push_back(sqrt((1/pow(B,2))*pow(errA,2) + (pow(A,2)/pow(B,4))*pow(errB,2)));
		cout << ly1.at(i) << " +- " << errly1.at(i) << endl;
	}
	cout << "\n\n\nlight yield lyso/csi (expected = 0.52):\n" << endl;
	for (int i = 0; i < veclyso.size(); i++) 
	{
		A = veclyso.at(i);
		B = veccsi.at(i);
		errA = vecerrlyso.at(i);
		errB = vecerrcsi.at(i);
		ly2.push_back(A/B);
		errly2.push_back(sqrt((1/pow(B,2))*pow(errA,2) + (pow(A,2)/pow(B,4))*pow(errB,2)));
		cout << ly2.at(i) << " +- " << errly2.at(i) << endl;
	}
	cout << "\n\n\nlight yield lyso/bgo (expected = 3.29):\n" << endl;
	for (int i = 0; i < veclyso.size(); i++) 
	{
		A = veclyso.at(i);
		B = vecbgo.at(i);
		errA = vecerrlyso.at(i);
		errB = vecerrbgo.at(i);
		ly3.push_back(A/B);
		errly3.push_back(sqrt((1/pow(B,2))*pow(errA,2) + (pow(A,2)/pow(B,4))*pow(errB,2)));
		cout << ly3.at(i) << " +- " << errly3.at(i) << endl;
	}

	cout << "\n\nLIGHT YIELDS: \n" << endl;
	double ly1mean = weightedmean(ly1, errly1);
	double errly1mean = weightederror(errly1);
	cout << "weighted mean " << ly1mean << " +- " << errly1mean << endl;
	cout << "comp: " << abs(ly1mean - 0.16)/errly1mean << endl << endl;
	double ly2mean = weightedmean(ly2, errly2);
	double errly2mean = weightederror(errly2);
	cout << "weighted mean " << ly2mean << " +- " << errly2mean << endl;
	cout << "comp: " << abs(ly2mean - 0.52)/errly2mean << endl << endl;
	double ly3mean = weightedmean(ly3, errly3);
	double errly3mean = weightederror(errly3);
	cout << "weighted mean " << ly3mean << " +- " << errly3mean << endl;
	cout << "comp: " << abs(ly3mean - 3.29)/errly3mean << endl << endl;
	
	
	
	//#########################################################################
	//############################# RISOLUZIONE ###############################
	//#########################################################################
	
	vector<double> risL, errrisL;
	vector<double> risB, errrisB;
	vector<double> risC, errrisC;
	
	cout << "\n\n\nRISOLUZIONE: \n\n" << endl;
	cout << "lyso: \n" << endl;
	for (int i = 0; i < veclyso.size(); i++) 
	{
		A = vecerrlyso.at(i);
		B = veclyso.at(i);
		errA = errFWHMlyso.at(i);
		errB = vecdevstL.at(i);
		//errB = vecerrlyso.at(i);
		risL.push_back(A/B);
		errrisL.push_back(sqrt((1/pow(B,2))*pow(errA,2) + (pow(A,2)/pow(B,4))*pow(errB,2)));
		cout << risL.at(i) << " +- " << errrisL.at(i) << endl;
	}
	
	cout << "\nbgo: \n" << endl;
	for (int i = 0; i < vecbgo.size(); i++) 
	{
		A = vecerrbgo.at(i);
		B = vecbgo.at(i);
		errA = errFWHMbgo.at(i);
		errB = vecdevstB.at(i);
		//errB = vecerrbgo.at(i);
		risB.push_back(A/B);
		errrisB.push_back(sqrt((1/pow(B,2))*pow(errA,2) + (pow(A,2)/pow(B,4))*pow(errB,2)));
		cout << risB.at(i) << " +- " << errrisB.at(i) << endl;
	}
	
	cout << "\ncsi: \n" << endl;
	for (int i = 0; i < veccsi.size(); i++) 
	{
		A = vecerrcsi.at(i);
		B = veccsi.at(i);
		errA = errFWHMcsi.at(i);
		errB = vecdevstC.at(i);
		//errB = vecerrcsi.at(i);
		risC.push_back(A/B);
		errrisC.push_back(sqrt((1/pow(B,2))*pow(errA,2) + (pow(A,2)/pow(B,4))*pow(errB,2)));
		cout << risC.at(i) << " +- " << errrisC.at(i) << endl;
	}
	
	
	//############## curve di calibrazione (aka risoluzione vs energia) ###############
	
	vector<double> energ = {122., 511., 1274.};
	vector<int> ind = {2, 0, 1}; // ho avuto la brillante idea di metterli nel vettore in un altro ordine e non ho voglia di cambiare tutto sopra.
	
	TCanvas* myC6 = new TCanvas("myC6","myC6");
	
	
	TGraphErrors graphL;
	for (int i = 0; i < risL.size(); i++) 
	{
		graphL.SetPoint(graphL.GetN(), energ.at(i), risL.at(ind.at(i)));
		graphL.SetPointError(graphL.GetN()-1, 0, errrisL.at(ind.at(i)));
	}
	graphL.SetTitle("Risoluzione energetica");
	graphL.GetXaxis()->SetTitle("Energia (keV)");
	graphL.GetXaxis()->SetRangeUser(80., 1330.);
	graphL.GetYaxis()->SetTitle("Risoluzione");
	graphL.GetYaxis()->SetRangeUser(0., 0.4);
	graphL.SetMarkerStyle(8);
	graphL.SetMarkerSize(0.8);
	graphL.SetMarkerColor(kMagenta);
	graphL.SetLineColor(kMagenta);
	
	TGraphErrors graphB;
	for (int i = 0; i < risB.size(); i++) 
	{
		graphB.SetPoint(graphB.GetN(), energ.at(i), risB.at(ind.at(i)));
		graphB.SetPointError(graphB.GetN()-1, 0, errrisB.at(ind.at(i)));
	}
	graphB.SetMarkerStyle(8);
	graphB.SetMarkerSize(0.8);
	graphB.SetMarkerColor(kBlue);
	graphB.SetLineColor(kBlue);
	
	
	TGraphErrors graphC;
	for (int i = 0; i < risC.size(); i++) 
	{
		graphC.SetPoint(graphC.GetN(), energ.at(i), risC.at(ind.at(i)));
		graphC.SetPointError(graphC.GetN()-1, 0, errrisC.at(ind.at(i)));
	}
	graphC.SetMarkerStyle(8);
	graphC.SetMarkerSize(0.8);
	graphC.SetMarkerColor(kBlack);
	graphC.SetLineColor(kBlack);
	
	TLegend* legend1 = new TLegend (0.1, 0.1, 0.3, 0.3); 
	legend1->AddEntry(&graphL,"curva LYSO","pl");
	legend1->AddEntry(&graphB,"curva BGO","pl");
	legend1->AddEntry(&graphC,"curva CsI", "pl");
	
	myC6->cd();
	graphL.Draw("APL");
	graphB.Draw("same");
	graphC.Draw("same");
	legend1->Draw();
	myC6-> Modified();
	myC6-> Update();
	
		
	
	myApp->Run();
	return 0;
}

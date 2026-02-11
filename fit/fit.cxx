#include <map>
#include <unistd.h>
#include <string>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "TStyle.h"
#include "TDirectory.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "RooRealVar.h"
#include "RooNovosibirsk.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooArgusBG.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooMCStudy.h"
#include "RooGaussian.h"
#include "RooKeysPdf.h"
#include "RooGaussian.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
using namespace RooFit;
using namespace std;
void draw_fit(std::string target, std::string file_name, Double_t down, Double_t up, std::string cut);
void draw_inmc(std::string target, std::string file_name, double down, double up, std::string cut);
void fit()
{
	draw_inmc("m_Kstar",      "m_Kstar", 0, 3, "m_Kshort!=0");
	draw_inmc("m_Kshort",     "m_Kshort", 0.4, 0.6, "m_Kshort!=0");
	draw_inmc("chisq_kmfit",  "chisq_kmfit", 0, 300, "m_Kshort!=0");
	draw_inmc("ratio_Kshort", "ratio_Kshort", 0, 10, "m_Kshort!=0");
	draw_inmc("m_Kstar",      "m_Kstar_cutted", 0, 3, "ratio_Kshort>2");
	draw_inmc("m_Kshort",     "m_Kshort_cutted", 0.4, 0.6, "ratio_Kshort>2");
	draw_inmc("chisq_kmfit",  "chisq_kmfit_cutted", 0, 300, "ratio_Kshort>2");
	draw_inmc("ratio_Kshort", "ratio_Kshort_cutted", 0, 10, "ratio_Kshort>2");
	Double_t lo = 0.3;
	Double_t up = 0.7;
	draw_fit("m_Kshort", "m_Kshort_fitted", lo, up, "m_Kshort!=0");
}


void draw_fit(std::string target, std::string file_name, Double_t lo, Double_t up, std::string cut)
{
	int n_bin(100);
	const unsigned int n_hist=6;
	TChain chain_data("chi2omg");
	TChain chain_inmc("chi2omg");
	chain_inmc.Add("/scratchfs/bes/cloudsail/JpsiToPiKKshort/master/run/topo/topoana.root");
	chain_data.Add("/scratchfs/bes/cloudsail/JpsiToPiKKshort/master/run/data.root");
	
	RooRealVar mass("mass", "mass", lo, up);
	RooRealVar mean = RooRealVar("mean", "mean",0, -0.5, 0.5);
	RooRealVar sigma = RooRealVar("sigma", "sigma", 0.1, 1e-4, 0.1);
	RooGaussian *gaus = new RooGaussian("gaus", "gaus", mass, mean, sigma);
	
	//int i_dcy_tr[n_hist] = {6,0,27,1,29,4};
	int i_dcy_tr[n_hist] = {4,29,1,27,0,6};
	int color[n_hist] = {3,4,5,6,7,2};
	int line_style[n_hist] = {10,10,10,10,10,9};

	char hs_title[200];
	sprintf(hs_title, "#bf{Distribution of %s accumulated}", target.c_str());
	THStack *hs = new THStack("hs",hs_title);
	TH1D *h[n_hist];
	RooDataHist *rh[n_hist];
	RooHistPdf  *rp[n_hist];
	RooRealVar  *rv[n_hist];
	RooExtendPdf *rep[n_hist];
	char hname[20], rhname[20], rpname[20], rvname[20], repname[20];
	char varexp[20];
	char selection[100];
	for (int i=0;i<n_hist;i++) {
		sprintf(hname, "h%i", i);
		sprintf(varexp, "%s>>h%i", target.c_str(), i);
		sprintf(selection, "iDcyTr==%i&&%s", i_dcy_tr[i], cut.c_str());
		h[i] = new TH1D(hname, hname, n_bin, lo, up);
		chain_inmc.Draw(varexp, selection);
		h[i]->SetFillColor(color[i]);
		hs->Add(h[i]);
		sprintf(rhname, "rh%i", i);
		rh[i] = new RooDataHist(rhname, rhname, mass, h[i]);
		sprintf(rpname, "rp%i", i);
		rp[i] = new RooHistPdf(rpname, rpname, mass, *rh[i], 3);
		sprintf(rvname, "rv%i", i);
		rv[i] = new RooRealVar(rvname, rvname, 0.01*h[i]->Integral(), 0, 3*h[i]->Integral());
		sprintf(repname, "rep%i", i);
		rep[i] = new RooExtendPdf(repname, repname, *rp[i], *rv[i]);
	}
	RooAbsPdf *sig = new RooFFTConvPdf("sig", "sig", mass, *rp[5], *gaus);
	rep[5] = new RooExtendPdf("rep5", "rep5", *sig, *rv[5]);
	RooAddPdf *model = new RooAddPdf("model", "model", RooArgList(*rep[0], *rep[1], *rep[2], *rep[3], *rep[4], *rep[5]));
	RooFitResult *fr;
	TH1D *hdata = new TH1D("hdata", "hdata", n_bin, lo, up);
	chain_data.Draw((target+">>hdata").c_str(), cut.c_str());
	RooDataHist *rhdata = new RooDataHist("rhdata", "rhdata", mass, hdata);
	fr= model->fitTo(*rhdata, RooFit::Extended(kTRUE),RooFit::Save(kTRUE),RooFit::Minos(kFALSE));
	RooPlot *frame = mass.frame();
	rhdata->plotOn(frame, MarkerStyle(20), LineWidth(3.0), MarkerSize(1.0));
	for (int i=0;i<n_hist;i++) {
		std::string components_name = "rep"+to_string(i);
		model->plotOn(frame, Components(components_name.c_str()), LineColor(color[i]), LineStyle(line_style[i]), LineWidth(3));
	}
	model->plotOn(frame, Components(RooArgList(*rep[0], *rep[1], *rep[2], *rep[3], *rep[4], *rep[5])),LineColor(color[5]), LineStyle(1));
	TCanvas* cc = new TCanvas("cc","cc",4800,1200);
	cc->Divide(4);

	RooHist *hresid = frame->residHist();
	RooHist *hpull = frame->pullHist();
	RooPlot *frame2 = mass.frame(Title("Residual Distribution"));
	frame2->addPlotable(hresid, "P");
	RooPlot *frame3 = mass.frame(Title("Pull Distribution"));
	frame3->addPlotable(hpull, "P");
	
	cc->cd(1);
	gPad->SetLeftMargin(0.15);
	hs->Draw();
	cc->cd(2);
	gPad->SetLeftMargin(0.15);
	frame->Draw();
	cc->cd(3);
	gPad->SetLeftMargin(0.15);
	frame2->Draw();
	cc->cd(4);
	gPad->SetLeftMargin(0.15);
	frame3->Draw();

	cc->cd(1);
	TLegend *legend = new TLegend(0.6175,0.3725,0.8825,0.8825);
	TLegendEntry *h_legend_entry[n_hist];
	char h_legend_entry_title[n_hist][200] = {
		"J/#psi #rightarrow #pi^{#plus} K^{0}_{S} K^{#minus}",
		"J/#psi #rightarrow #pi^{0} #pi^{#plus} #pi^{#plus} #pi^{#minus} #pi^{#minus}",
		"J/#psi #rightarrow #pi^{#plus} #pi^{#minus} #omega",
		"J/#psi #rightarrow #pi^{#plus} #pi^{#plus} #pi^{#minus} #pi^{#minus}",
		"J/#psi #rightarrow #rho^{0} a^{0}_{2}",
		"J/#psi #rightarrow #pi^{#minus} #bar{K^{*}} K^{* #plus}"
	};
	for(int i=n_hist-1;i>=0;i--) {
		h_legend_entry[i] = legend->AddEntry(h[i],h_legend_entry_title[i],"f");
		h_legend_entry[i]->SetTextColor(color[i]);
	}
	legend->Draw();
	cc->Print((file_name+".png").c_str());
	std::cout << "Chisq: " << frame->chiSquare() << std::endl;
}

void draw_inmc(std::string target, std::string file_name, double down, double up, std::string cut)
{
	TCanvas* cc = new TCanvas("cc","cc",1);
	int Nbin(100), Npars(5);

	TChain chain_data("chi2omg");
	TChain chain_inmc("chi2omg");

	chain_inmc.Add("/scratchfs/bes/cloudsail/JpsiToPiKKshort/master/run/topo/topoana.root");
	
	chain_data.Add("/scratchfs/bes/cloudsail/JpsiToPiKKshort/master/run/data/data_master_000.root");
	chain_data.Add("/scratchfs/bes/cloudsail/JpsiToPiKKshort/master/run/data/data_master_001.root");
	chain_data.Add("/scratchfs/bes/cloudsail/JpsiToPiKKshort/master/run/data/data_master_002.root");
	chain_data.Add("/scratchfs/bes/cloudsail/JpsiToPiKKshort/master/run/data/data_master_003.root");
	chain_data.Add("/scratchfs/bes/cloudsail/JpsiToPiKKshort/master/run/data/data_master_004.root");
	chain_data.Add("/scratchfs/bes/cloudsail/JpsiToPiKKshort/master/run/data/data_master_005.root");
	chain_data.Add("/scratchfs/bes/cloudsail/JpsiToPiKKshort/master/run/data/data_master_006.root");
	chain_data.Add("/scratchfs/bes/cloudsail/JpsiToPiKKshort/master/run/data/data_master_007.root");
	chain_data.Add("/scratchfs/bes/cloudsail/JpsiToPiKKshort/master/run/data/data_master_008.root");
	chain_data.Add("/scratchfs/bes/cloudsail/JpsiToPiKKshort/master/run/data/data_master_009.root");
	

	char hs_title[200];
	sprintf(hs_title, "#bf{Distribution of %s accumulated}", target.c_str());
	THStack *hs = new THStack("hs",hs_title);
	const unsigned int n_hist=6;
	TH1D *h[n_hist];
	char hname[20];
	char varexp[20];
	char selection[100];
	int i_dcy_tr[n_hist] = {6,0,27,1,29,4};
	//int i_dcy_tr[n_hist] = {4,29,1,27,0,6};
	int h_fill_color[n_hist] = {2,3,4,6,7,8};
	TLegend *legend = new TLegend(0.6175,0.3725,0.8825,0.8825);
	TLegendEntry *h_legend_entry[n_hist];
	char h_legend_entry_title[n_hist][200] = {
		"J/#psi #rightarrow #pi^{#plus} K^{0}_{S} K^{#minus}",
		"J/#psi #rightarrow #pi^{0} #pi^{#plus} #pi^{#plus} #pi^{#minus} #pi^{#minus}",
		"J/#psi #rightarrow #pi^{#plus} #pi^{#minus} #omega",
		"J/#psi #rightarrow #pi^{#plus} #pi^{#plus} #pi^{#minus} #pi^{#minus}",
		"J/#psi #rightarrow #rho^{0} a^{0}_{2}",
		"J/#psi #rightarrow #pi^{#minus} #bar{K^{*}} K^{* #plus}"
	};
	for(int i=n_hist-1;i>=0;i--) {
		sprintf(hname, "h%i", i);
		sprintf(varexp, "%s>>h%i", target.c_str(), i);
		sprintf(selection, "iDcyTr==%i&&%s", i_dcy_tr[i], cut.c_str());
		h[i] = new TH1D(hname, hname, 120, down, up);
		chain_inmc.Draw(varexp, selection);
		h[i]->SetFillColor(h_fill_color[i]);
		hs->Add(h[i]);

		h_legend_entry[i] = legend->AddEntry(h[i],h_legend_entry_title[i],"f");
		h_legend_entry[i]->SetTextColor(h_fill_color[i]);
	}
	hs->Draw();
	legend->Draw();

	cc->Print((file_name+".png").c_str());
}


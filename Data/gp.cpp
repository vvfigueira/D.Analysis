#include "TCanvas.h"
#include "TGraphErrors.h"
#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TF1.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TRatioPlot.h"
#include "TMath.h"
#include "TMathText.h"
#include "TPad.h"

int npar = 4;
int q = 1;

Double_t xmin = 6000, xmax = 8300;

Double_t Malus(Double_t *x, Double_t *par) {
   return TMath::Power(par[0]*TMath::Cos(TMath::DegToRad()*(x[0]-par[1]))+par[2],2)+par[3];
}

Double_t Per(Double_t *x, Double_t *par){
    return par[0] + par[1]*TMath::Cos(x[0]*par[2]+par[3]);
}

void gp(){
    
    TLatex *latex = new TLatex();

    TCanvas *c1 = new TCanvas();
    auto pad1 = new TPad("graph","",0.,0.3,1,0.95);

    auto pad2 = new TPad("rezid","",0.,0.,1,0.32);

    pad1->SetFillStyle(4000);
    pad2->SetFillStyle(4000);
    pad1->SetMargin(0.08,0.08,0.01,0);
    pad2->SetMargin(0.08,0.08,0.21,0.01);
    pad1->SetGridx();
    pad1->SetGridy();
    pad2->SetGridx();
    pad2->SetGridy();
    pad2->Draw();
    pad1->Draw();
    

    c1->SetTickx();
    c1->SetTicky();
    c1->SetGridx();
    c1->SetGridy();

    TGraphErrors *gr = new TGraphErrors();

    //gr->GetXaxis()->SetTitle("#bf{[Graus]}");
    gr->GetYaxis()->SetTitle("Residuos");//#bf{[]}");
    gr->GetXaxis()->SetLabelFont(132);
    gr->GetXaxis()->SetTitleFont(132);
    gr->GetYaxis()->SetLabelFont(132);
    gr->GetYaxis()->SetTitleFont(132);
    gr->GetYaxis()->SetLabelSize(0.05);
    gr->GetYaxis()->SetTitleSize(0.05);
    gr->GetXaxis()->SetLabelSize(0.05);
    gr->GetXaxis()->SetTitleSize(0.05);
    gr->GetYaxis()->SetTitleOffset(0.5);
    //gr->SetMarkerStyle(kFullCircle);
    //gr->SetMarkerSize(0.7);

    std::fstream file;
    file.open("DLP2-1.6.tsv", std::ios::in);

    double x1[2],ex = 20, ey = 0.003;

    int n = 0;

    while(1){
        file >> x1[0] >> x1[1]; //>> x1[2] >> x1[3];
    
        n = gr->GetN();

        gr->SetPoint(n, x1[0], x1[1]);
        gr->SetPointError(n, ex, ey);

        if(file.eof()) break;
    
    }

    file.close();

    gr->GetXaxis()->SetLimits(xmin, xmax);
    //gr->GetYaxis()->SetRangeUser(6, 102);
    pad1->cd();
    gr->Draw("AP");

    TF1 *mal = new TF1("mal", Per, xmin, xmax, npar);
    mal->SetLineColor(kRed);
    mal->SetParameter(0, 0);
    mal->SetParameter(1, 0.01);
    mal->SetParameter(2, TMath::TwoPi()/365);
    mal->SetParameter(3, -152*TMath::TwoPi()/365);

    gr->Fit("mal","","",xmin,xmax);

    Double_t par[npar];

    mal->GetParameters(par);
    
    //std::cout << gr->Chisquare(mal) << "\n";

    auto *red = new TGraphErrors();

    int m = 0;

    std::fstream file2;
    file2.open("DLP2-1.6.tsv", std::ios::in);

    Double_t x2[2];

    while(1){
        
        file2 >> x2[0] >> x2[1];// >> x2[2] >> x2[3];
    
        m = red->GetN();

        red->SetPoint(m, x2[0], -(Malus(x2,par) - x2[1]));
        red->SetPointError(m, ex, ey);

        if(file2.eof()) break;

    }

    file2.close();

    red->GetXaxis()->SetLimits(xmin, xmax);
    /*auto c2 = new TCanvas();

    c2->SetTickx();
    c2->SetTicky();
    c2->SetGridx();
    c2->SetGridy();*/

    red->GetXaxis()->SetTitle("Dias");//#bf{[D]}");
    //red->GetYaxis()->SetTitle("#bf{[V]}");
    red->GetXaxis()->SetLabelFont(132);
    red->GetXaxis()->SetTitleFont(132);
    red->GetYaxis()->SetLabelFont(132);
    red->GetYaxis()->SetTitleFont(132);
    red->GetYaxis()->SetTitleOffset(1.3);
    red->GetYaxis()->SetLabelSize(0.1);
    red->GetYaxis()->SetTitleSize(0.1);
    red->GetXaxis()->SetLabelSize(0.1);
    red->GetXaxis()->SetTitleSize(0.1);
    //red->GetYaxis()->SetRangeUser(-0.05, 0.058);

    auto fnc = new TF1("fnc",Malus,xmin,xmax,npar);

    fnc->SetParameters(par);
    pad2->cd();
    red->Draw("AP");

    //fnc->Draw("SAME");

    /*TLegend *leg = new TLegend(0.68, 0.705, 0.82, 0.805);
    leg->SetBorderSize(0);
    leg->AddEntry(gr, "Dados", "p");
    leg->AddEntry(mal, "Ajuste", "l");
    leg->Draw();*/

}
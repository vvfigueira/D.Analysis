#include "TCanvas.h"
#include "TGraphErrors.h"
#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TRatioPlot.h"
#include "TMath.h"
#include "TMathText.h"
#include "TPad.h"
#include "TRatioPlot.h"

int npar = 4;
int q = 1;
char* name = "DamaP2.2-6Error.tsv";
double confi = 0.05;

Double_t xmin = 3000, xmax = 8500;

Double_t Cnst(Double_t* x,Double_t* par){
    return confi;
}

Double_t Cnst2(Double_t* x,Double_t* par){
    return -confi;
}

Double_t Malus(Double_t *x, Double_t *par) {
   return TMath::Power(par[0]*TMath::Cos(TMath::DegToRad()*(x[0]-par[1]))+par[2],2)+par[3];
}

Double_t Per(Double_t *x, Double_t *par){
    return par[0] + par[1]*TMath::Cos(x[0]*par[2]+par[3]);
}

void gp(){

    std::fstream file;
    double sumsq = 0, prod = 0, autocorr = 0, x1[4], x2[4];
    Double_t par[npar];
    int n = 0, m = 0;
    
    TLatex *Latex = new TLatex();

    TCanvas *Princ = new TCanvas();

    Princ->SetTickx();
    Princ->SetTicky();
    Princ->SetGridx();
    Princ->SetGridy();

    auto GPad = new TPad("Graph","",0.,0.3,1,0.95);

    GPad->SetFillStyle(4000);
    GPad->SetMargin(0.08,0.08,0.01,0.05);
    GPad->SetGridx();
    GPad->SetGridy();
    GPad->Draw();

    auto RPad = new TPad("Rezid","",0.,0.,1,0.32);
    
    RPad->SetFillStyle(4000);
    RPad->SetMargin(0.08,0.08,0.21,0.01);
    RPad->SetGridx();
    RPad->SetGridy();
    RPad->Draw();

    TGraphErrors *Graph = new TGraphErrors();

    Graph->GetYaxis()->SetTitle("Residuos");
    Graph->GetXaxis()->SetLabelFont(132);
    Graph->GetXaxis()->SetTitleFont(132);
    Graph->GetYaxis()->SetLabelFont(132);
    Graph->GetYaxis()->SetTitleFont(132);
    Graph->GetYaxis()->SetLabelSize(0.05);
    Graph->GetYaxis()->SetTitleSize(0.05);
    Graph->GetXaxis()->SetLabelSize(0.05);
    Graph->GetXaxis()->SetTitleSize(0.05);
    Graph->GetYaxis()->SetTitleOffset(0.65);
    Graph->SetMarkerStyle(kOpenCircle);
    Graph->SetMarkerSize(0.6);
    Graph->GetXaxis()->SetLimits(xmin, xmax);
    Graph->GetYaxis()->SetRangeUser(-0.07, 0.07);

    file.open(name, std::ios::in);

    while(1){
        file >> x1[0] >> x1[1]>> x1[2] >> x1[3];
    
        n = Graph->GetN();

        Graph->SetPoint(n, x1[0], x1[2]);
        Graph->SetPointError(n, x1[1], x1[3]);

        if(file.eof()) break;
    
    }

    file.close();

    GPad->cd();
    Graph->Draw("AP");

    for(int i=0;i<Graph->GetN();i++){
        sumsq = sumsq + pow(Graph->GetPointY(i)-Graph->GetMean(2),2);
    }

    TF1 *Ajt = new TF1("Ajt", Per, xmin, xmax, npar);

    Ajt->SetLineColor(kRed);
    Ajt->SetParameter(0, 0);
    Ajt->SetParameter(1, 0.01);
    Ajt->SetParameter(2, TMath::TwoPi()/365);
    Ajt->SetParameter(3, -152*TMath::TwoPi()/365);
    Ajt->SetNpx(1.e5);

    Graph->Fit("Ajt","","",xmin,xmax);

    Ajt->GetParameters(par);

    auto *Rezid = new TGraphErrors();

    Rezid->GetXaxis()->SetTitle("Dias");
    Rezid->GetXaxis()->SetLabelFont(132);
    Rezid->GetXaxis()->SetTitleFont(132);
    Rezid->GetYaxis()->SetLabelFont(132);
    Rezid->GetYaxis()->SetTitleFont(132);
    Rezid->GetYaxis()->SetTitleOffset(1.3);
    Rezid->GetYaxis()->SetLabelSize(0.1);
    Rezid->GetYaxis()->SetTitleSize(0.1);
    Rezid->GetXaxis()->SetLabelSize(0.1);
    Rezid->GetXaxis()->SetTitleSize(0.1);
    Rezid->SetMarkerStyle(kOpenCircle);
    Rezid->SetMarkerSize(0.6);
    Rezid->GetYaxis()->SetNdivisions(8);
    Rezid->GetXaxis()->SetLimits(xmin, xmax);

    file.open(name, std::ios::in);

    while(1){
        
        file >> x2[0] >> x2[1] >> x2[2] >> x2[3];
    
        m = Rezid->GetN();

        Rezid->SetPoint(m, x2[0], -(Per(x2,par) - x2[2])/x2[3]);
        Rezid->SetPointError(m, x2[1], x2[3]);

        if(file.eof()) ;

    }

    file.close();

    RPad->cd();
    Rezid->Draw("AL");

    TCanvas *RP = new TCanvas();

    auto h1 = new TH1D("h1","h1",200,xmin,xmax);
    double x,y;
    for(int i=0; i < Graph->GetN(); ++i) {
        Graph->GetPoint(i, x, y);
        h1->Fill(x,y);
    }
    h1->GetXaxis()->SetLabelFont(132);
    h1->GetXaxis()->SetTitleFont(132);
    h1->GetYaxis()->SetLabelFont(132);
    h1->GetYaxis()->SetTitleFont(132);
    //h1->GetYaxis()->SetTitleOffset(1.3);
    //h1->GetYaxis()->SetLabelSize(0.1);
    //h1->GetYaxis()->SetTitleSize(0.1);
    //h1->GetXaxis()->SetLabelSize(0.1);
    //h1->GetXaxis()->SetTitleSize(0.1);
    h1->SetMarkerStyle(kOpenCircle);
    h1->SetMarkerSize(0.6);
    h1->SetMarkerColor(kBlack);

    h1->Fit("Ajt","","",xmin,xmax);

    auto rp1 = new TRatioPlot(h1,"divsym");
    RP->cd();
    rp1->Draw();
    rp1->GetLowerRefYaxis()->SetTitle("ratio");
    rp1->GetUpperRefYaxis()->SetTitle("entries");

    TCanvas *CorDat = new TCanvas();

    CorDat->SetTickx();
    CorDat->SetTicky();
    CorDat->SetGridx();
    CorDat->SetGridy();
    
    TGraphErrors *GCorDat = new TGraphErrors();

    int N =Graph->GetN();

    GCorDat->SetTitle("Auto-Correlacao dos Dados");
    GCorDat->GetYaxis()->SetTitle("Auto-Correlacao");
    GCorDat->GetXaxis()->SetTitle("Defasagem");
    GCorDat->GetXaxis()->SetLabelFont(132);
    GCorDat->GetXaxis()->SetTitleFont(132);
    GCorDat->GetYaxis()->SetLabelFont(132);
    GCorDat->GetYaxis()->SetTitleFont(132);
    GCorDat->GetYaxis()->SetLabelSize(0.035);
    GCorDat->GetYaxis()->SetTitleSize(0.035);
    GCorDat->GetXaxis()->SetLabelSize(0.035);
    GCorDat->GetXaxis()->SetTitleSize(0.035);
    GCorDat->GetYaxis()->SetTitleOffset(1);
    GCorDat->SetMarkerStyle(kFullCircle);
    GCorDat->SetMarkerSize(0.6);
    

    for(int k=-Graph->GetN();k<Graph->GetN();k++){
        
        prod = 0;
        for(int j=0;j<Graph->GetN();j++){
            
            if(j-k >= 0 && j-k < n){
                prod = prod + (Graph->GetPointY(j-k)-Graph->GetMean(2))*(Graph->GetPointY(j)-Graph->GetMean(2));
            }else{prod = prod + (-Graph->GetMean(2))*(Graph->GetPointY(j)-Graph->GetMean(2));}
        }

        autocorr = prod/sumsq;

        GCorDat->SetPoint(k+n, k, autocorr);
        GCorDat->SetPointError(n+k, 0, 0);
    
    }

    GCorDat->GetXaxis()->SetLimits(-Graph->GetN(),Graph->GetN());
    CorDat->cd();
    GCorDat->Draw("AL");

    TCanvas *CorRez = new TCanvas();

    CorRez->SetTickx();
    CorRez->SetTicky();
    CorRez->SetGridx();
    CorRez->SetGridy();

    TGraphErrors *GCorRez = new TGraphErrors();

    GCorRez->SetTitle("Auto-Correlacao dos Residuos");
    GCorRez->GetYaxis()->SetTitle("Auto-Correlacao");
    GCorRez->GetXaxis()->SetTitle("Defasagem");
    GCorRez->GetXaxis()->SetLabelFont(132);
    GCorRez->GetXaxis()->SetTitleFont(132);
    GCorRez->GetYaxis()->SetLabelFont(132);
    GCorRez->GetYaxis()->SetTitleFont(132);
    GCorRez->GetYaxis()->SetLabelSize(0.035);
    GCorRez->GetYaxis()->SetTitleSize(0.035);
    GCorRez->GetXaxis()->SetLabelSize(0.035);
    GCorRez->GetXaxis()->SetTitleSize(0.035);
    GCorRez->GetYaxis()->SetTitleOffset(1);
    GCorRez->SetMarkerStyle(kFullCircle);
    GCorRez->SetMarkerSize(0.6);

    sumsq = 0, prod=0,autocorr=0;
    for(int i=0;i<Rezid->GetN();i++){
        sumsq = sumsq + pow(Graph->GetPointY(i)-Rezid->GetMean(2),2);
    }

    for(int k=-Rezid->GetN();k<Rezid->GetN();k++){
        
        prod = 0;
        for(int j=0;j<Rezid->GetN();j++){
            
            if(j-k >= 0 && j-k < Rezid->GetN()){
                prod = prod + (Rezid->GetPointY(j-k)-Rezid->GetMean(2))*(Rezid->GetPointY(j)-Rezid->GetMean(2));
            }else{prod = prod + (-Rezid->GetMean(2))*(Rezid->GetPointY(j)-Rezid->GetMean(2));}
        }

        autocorr = prod/sumsq;

        GCorRez->SetPoint(k+n, k, autocorr);
        GCorRez->SetPointError(n+k, 0, 0);
    
    }
    
    GCorRez->GetXaxis()->SetLimits(-Rezid->GetN(), Rezid->GetN());
    CorRez->cd();
    GCorRez->Draw("AL");

    TF1 *Confi1 = new TF1("Confi1", Cnst, -10*Rezid->GetN(), 10*Rezid->GetN(), 0);

    Confi1->SetFillColor(kRed);
    Confi1->Draw("same");

    TF1 *Confi2 = new TF1("Confi2", Cnst2, -10*Rezid->GetN(), 10*Rezid->GetN(), 0);
    
    Confi2->SetFillColor(kRed);
    Confi2->Draw("same");

    /*TLegend *leg = new TLegend(0.68, 0.705, 0.82, 0.805);
    leg->SetBorderSize(0);
    leg->AddEntry(Graph, "Dados", "p");
    leg->AddEntry(Ajt, "Ajuste", "l");
    leg->Draw();*/

}
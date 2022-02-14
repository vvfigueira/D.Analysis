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

char* name2 = "D1.tsv";
char* name1 = "D2.tsv";
double confi = 0.05;

Double_t Cnst(Double_t* x,Double_t* par){
    return confi;
}

Double_t Cnst2(Double_t* x,Double_t* par){
    return -confi;
}

void ac(){
    
    Double_t x1[4];
    int n,m;
    std::fstream file;
    double sumsq1 = 0, sumsq2 = 0, prod = 0, autocorr = 0;
    double Nmin;

    TLatex *latex = new TLatex();

    TCanvas *C1 = new TCanvas();

    C1->SetTickx();
    C1->SetTicky();
    C1->SetGridx();
    C1->SetGridy();

    TGraphErrors *GTemp1 = new TGraphErrors();

    file.open(name1, std::ios::in);

    while(1){
        file >> x1[0] >> x1[1]>> x1[2] >> x1[3];
    
        n = GTemp1->GetN();

        GTemp1->SetPoint(n, x1[0], x1[2]);
        GTemp1->SetPointError(n, x1[1], x1[3]);

        if(file.eof()) break;
    
    }

    file.close();

    TGraphErrors *GTemp2 = new TGraphErrors();

    file.open(name2, std::ios::in);

    while(1){
        file >> x1[0] >> x1[1]>> x1[2] >> x1[3];
    
        m = GTemp2->GetN();

        GTemp2->SetPoint(m, x1[0], x1[2]);
        GTemp2->SetPointError(m, x1[1], x1[3]);

        if(file.eof()) break;
    
    }

    file.close();

    for(int i=0;i<GTemp1->GetN();i++){
        sumsq1 = sumsq1 + pow(GTemp1->GetPointY(i)-GTemp1->GetMean(2),2);
    }

    sumsq1 = sqrt(sumsq1);

    for(int i=0;i<GTemp2->GetN();i++){
        sumsq2 = sumsq2 + pow(GTemp2->GetPointY(i)-GTemp2->GetMean(2),2);
    }

    sumsq2 = sqrt(sumsq2);

    TGraphErrors *GAcr = new TGraphErrors();

    for(int k= -GTemp1->GetN();k < GTemp1->GetN();k++){
        
        prod = 0;
        for(int j = 0;j < GTemp2->GetN();j++){
            
            if(j-k >= 0 && j-k < GTemp1->GetN()){
                prod = prod + (GTemp1->GetPointY(j-k)- GTemp1->GetMean(2))*(GTemp2->GetPointY(j) - GTemp2->GetMean(2));
            }else{prod = prod + ( -GTemp1->GetMean(2))*(GTemp2->GetPointY(j) -GTemp2->GetMean(2));}
        }

        autocorr = prod/(sumsq1*sumsq2);

        GAcr->SetPoint(k+GTemp1->GetN(), k, autocorr);
        GAcr->SetPointError(GTemp1->GetN()+k, 0, 0);
    
    }

    GAcr->GetYaxis()->SetTitle("Auto-Correlacao");
    GAcr->GetXaxis()->SetTitle("Defasagem");
    GAcr->GetXaxis()->SetLabelFont(132);
    GAcr->GetXaxis()->SetTitleFont(132);
    GAcr->GetYaxis()->SetLabelFont(132);
    GAcr->GetYaxis()->SetTitleFont(132);
    GAcr->GetYaxis()->SetLabelSize(0.035);
    GAcr->GetYaxis()->SetTitleSize(0.035);
    GAcr->GetXaxis()->SetLabelSize(0.035);
    GAcr->GetXaxis()->SetTitleSize(0.035);
    GAcr->GetYaxis()->SetTitleOffset(1);
    GAcr->SetMarkerStyle(kFullCircle);
    GAcr->SetMarkerSize(0.6);
    GAcr->GetXaxis()->SetLimits(-GTemp1->GetN(), GTemp1->GetN());

    C1->cd();
    GAcr->Draw("AP");

    TF1 *R1 = new TF1("R1", Cnst, -GTemp1->GetN(), GTemp1->GetN(), 0);

    R1->SetFillColor(kRed);
    R1->Draw("same");

    TF1 *R2 = new TF1("R2", Cnst2, -GTemp1->GetN(), GTemp1->GetN(), 0);
    
    R2->SetFillColor(kRed);
    R2->Draw("same");
}
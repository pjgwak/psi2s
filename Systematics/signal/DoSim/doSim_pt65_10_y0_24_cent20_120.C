#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <Fit/Fitter.h>
#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Math/WrappedMultiTF1.h>
#include <HFitInterface.h>
#include "../../../headers/commonUtility.h"
#include "../../../headers/cutsAndBin.h"
#include "../../../headers/HiEvtPlaneList.h"
#include "../../../headers/Style.h"
#include "../../../headers/tdrstyle.C"
#include "../../../headers/CMS_lumi_v2mass.C"
using namespace std;

const int nParmM = 11;
const int nParmV = 16;
Int_t iparmass[nParmM] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
Int_t iparvn[nParmV] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13,14,15};

struct GlobalChi2_width
{
    GlobalChi2_width(ROOT::Math::IMultiGenFunction & f1, ROOT::Math::IMultiGenFunction & f2):
        fChi2_1(&f1), fChi2_2(&f2) {}

    Double_t operator() (const double *par) const
    {
        Double_t p1[nParmM];
        for(Int_t i = 0; i < nParmM; i++) p1[i] = par[iparmass[i]];
        Double_t p2[nParmV];
        for(Int_t i = 0; i < nParmV; i++) p2[i] = par[iparvn[i]];
        return (*fChi2_1)(p1) + (*fChi2_2)(p2);
    }
    const ROOT::Math::IMultiGenFunction * fChi2_1;
    const ROOT::Math::IMultiGenFunction * fChi2_2;
};

//totalYield{{{
Double_t TotalYield(Double_t* x, Double_t* par)
{
    Double_t N1 = par[0]; //Number of Psi2S yield
    Double_t Nbkg = par[1]; //Nuber of Bkg
    Double_t mean = par[2]; //Crystall Ball mean
    Double_t sigma = par[3]; //Crystall Ball sigma
    Double_t alpha = par[4]; //crystall ball alpha
    Double_t n = par[5]; //Crystall ball n
    Double_t ratio = par[6]; //For fraction of Double Crystall ball
    Double_t frac = par[7]; //Crystall ball f
    Double_t cheb0 = par[8];
    Double_t cheb1 = par[9];
    Double_t cheb2 = par[10];
    Double_t sigma1_2 = sigma*ratio; //

    //t2 > t1
    Double_t JPsi_t1 = (x[0]-mean)/sigma;
    Double_t JPsi_t2 = (x[0]-mean)/sigma1_2;

    if (alpha < 0)
    {
        cout << "ERROR ::: alpha variable negative!!!! " << endl;
        return -1;
    }

    //Crystall ball fucntion
    Double_t absAlpha = TMath::Abs(alpha);
    Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
    Double_t b = n/absAlpha - absAlpha;

    Double_t JPsi_1 = -1;
    Double_t JPsi_2 = -1;

    if(JPsi_t1 > -alpha){
        JPsi_1 = exp(-JPsi_t1*JPsi_t1/2.);
    }
    if(JPsi_t1 <= -alpha){
        JPsi_1 = a*TMath::Power((b-JPsi_t1),-n);
    }

    JPsi_2=1/(sigma1_2*TMath::Sqrt(2*TMath::Pi()))*exp(-0.5*JPsi_t2*JPsi_t2);
    Double_t fC = n/absAlpha*1/(n-1)*exp(-absAlpha*absAlpha/2);
    Double_t fD = TMath::Sqrt(TMath::Pi()/2)*(1+TMath::Erf(absAlpha/TMath::Sqrt(2)));
    Double_t fN_1 = 1./(sigma*(fC+fD));

    double shx = (10*x[0]-37)/3;
    Double_t BkgM = Nbkg*(1+cheb0*shx + cheb1*(2*shx*shx-1) + cheb2*(4*shx*shx*shx-3*shx));

    return N1*(fN_1*frac*JPsi_1 + (1-frac)*JPsi_2) + BkgM;
}

//totalYieldSig{{{
Double_t TotalYieldSig(Double_t* x, Double_t* par)
{
    Double_t N1 = par[0];
    Double_t mean = par[1];
    Double_t sigma = par[2];
    Double_t alpha = par[3];
    Double_t n = par[4];
    Double_t ratio = par[5];
    Double_t frac = par[6];
    Double_t normMin = par[7];
    Double_t sigma1_2 = sigma*ratio;

    //t2 > t1
    Double_t JPsi_t1 = (x[0]-mean)/sigma;
    Double_t JPsi_t2 = (x[0]-mean)/sigma1_2;

    if (alpha < 0)
    {
        cout << "ERROR ::: alpha variable negative!!!! " << endl;
        return -1;
    }

    Double_t absAlpha = fabs((Double_t)alpha);
    Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
    Double_t b = n/absAlpha - absAlpha;
    Double_t JPsi_1 = -1;
    Double_t JPsi_2 = -1;

    if(JPsi_t1 > -alpha){
        JPsi_1 = exp(-JPsi_t1*JPsi_t1/2.);
    }
    else if(JPsi_t1 <= -alpha){
        JPsi_1 = a*TMath::Power((b-JPsi_t1),-n);
    }

    JPsi_2=1/(sigma1_2*TMath::Sqrt(2*TMath::Pi()))*exp(-0.5*JPsi_t2*JPsi_t2);
    Double_t fC = n/absAlpha*1/(n-1)*exp(-absAlpha*absAlpha/2);
    Double_t fD = TMath::Sqrt(TMath::Pi()/2)*(1+TMath::Erf(absAlpha/TMath::Sqrt(2)));
    Double_t fN_1 = 1./(sigma*(fC+fD));

    normMin = 0;  // y axis cosmetic
    return normMin+N1*(fN_1*frac*JPsi_1 + (1-frac)*JPsi_2);
}

//totalvn pol2 bkg Psi2S{{{
Double_t Totalvnpol1JPsi(Double_t* x, Double_t* par)
{
    Double_t N1 = par[0];
    Double_t Nbkg = par[1];
    Double_t mean = par[2];
    Double_t sigma = par[3];
    Double_t alpha = par[4];
    Double_t n = par[5];
    Double_t ratio = par[6];
    Double_t frac = par[7];
    Double_t cheb0 = par[8];
    Double_t cheb1 = par[9];
    Double_t cheb2 = par[10];
    Double_t c = par[11];
    Double_t c1 = par[12];
    Double_t c2 = par[13];
    Double_t c3 = par[14];

    Double_t sigma1_2 = sigma*ratio;

    //t2 > t1
    Double_t JPsi_t1 = (x[0]-mean)/sigma;
    Double_t JPsi_t2 = (x[0]-mean)/sigma1_2;
    if (alpha < 0)
    {
        cout << "ERROR ::: alpha variable negative!!!! " << endl;
        return -1;
    }

    Double_t absAlpha = fabs((Double_t)alpha);
    Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
    Double_t b = n/absAlpha - absAlpha;

    Double_t JPsi_1 = -1;
    Double_t JPsi_2 = -1;

    if(JPsi_t1 > -alpha){
        JPsi_1 = exp(-JPsi_t1*JPsi_t1/2.);
    }
    else if(JPsi_t1 <= -alpha){
        JPsi_1 = a*TMath::Power((b-JPsi_t1),-n);
    }

    JPsi_2=1/(sigma1_2*TMath::Sqrt(2*TMath::Pi()))*exp(-0.5*JPsi_t2*JPsi_t2);
    Double_t fC = n/absAlpha*1/(n-1)*exp(-absAlpha*absAlpha/2);
    Double_t fD = TMath::Sqrt(TMath::Pi()/2)*(1+TMath::Erf(absAlpha/TMath::Sqrt(2)));
    Double_t fN_1 = 1./(sigma*(fC+fD));
    Double_t SigM = N1*(fN_1*frac*JPsi_1 + (1-frac)*JPsi_2);
    
    double shx = (10*x[0]-37)/3;
    Double_t BkgM = Nbkg*(1+cheb0*shx + cheb1*(2*shx*shx-1) + cheb2*(4*shx*shx*shx-3*shx));

    return c*(SigM/(SigM+BkgM)) + (1 - SigM/(SigM+BkgM))*(c2 + c1*x[0]);

}
//}}}

//totalvn pol2 bkg Psi2S{{{
Double_t Totalvnpol2JPsi(Double_t* x, Double_t* par)
{
    Double_t N1 = par[0];
    Double_t Nbkg = par[1];
    Double_t mean = par[2];
    Double_t sigma = par[3];
    Double_t alpha = par[4];
    Double_t n = par[5];
    Double_t ratio = par[6];
    Double_t frac = par[7];
    Double_t cheb0 = par[8];
    Double_t cheb1 = par[9];
    Double_t cheb2 = par[10];
    Double_t c = par[11];
    Double_t c1 = par[12];
    Double_t c2 = par[13];
    Double_t c3 = par[14];

    Double_t sigma1_2 = sigma*ratio;

    //t2 > t1
    Double_t JPsi_t1 = (x[0]-mean)/sigma;
    Double_t JPsi_t2 = (x[0]-mean)/sigma1_2;
    if (alpha < 0)
    {
        cout << "ERROR ::: alpha variable negative!!!! " << endl;
        return -1;
    }

    Double_t absAlpha = fabs((Double_t)alpha);
    Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
    Double_t b = n/absAlpha - absAlpha;

    Double_t JPsi_1 = -1;
    Double_t JPsi_2 = -1;

    if(JPsi_t1 > -alpha){
        JPsi_1 = exp(-JPsi_t1*JPsi_t1/2.);
    }
    else if(JPsi_t1 <= -alpha){
        JPsi_1 = a*TMath::Power((b-JPsi_t1),-n);
    }

    JPsi_2=1/(sigma1_2*TMath::Sqrt(2*TMath::Pi()))*exp(-0.5*JPsi_t2*JPsi_t2);
    Double_t fC = n/absAlpha*1/(n-1)*exp(-absAlpha*absAlpha/2);
    Double_t fD = TMath::Sqrt(TMath::Pi()/2)*(1+TMath::Erf(absAlpha/TMath::Sqrt(2)));
    Double_t fN_1 = 1./(sigma*(fC+fD));
    Double_t SigM = N1*(fN_1*frac*JPsi_1 + (1-frac)*JPsi_2);

    double shx = (10*x[0]-37)/3;
    Double_t BkgM = Nbkg*(1+cheb0*shx + cheb1*(2*shx*shx-1) + cheb2*(4*shx*shx*shx-3*shx));

    return c*(SigM/(SigM+BkgM)) + (1 - SigM/(SigM+BkgM))*(c3 + c2*x[0] + c1*x[0]*x[0]);

}
//}}}

//totalYieldBkg{{{
Double_t TotalYieldBkg(Double_t* x, Double_t* par)
{
    Double_t Nbkg = par[0]; //Nuber of Bkg
    Double_t cheb0 = par[1];
    Double_t cheb1 = par[2];
    Double_t cheb2 = par[3];

    double shx = (10*x[0]-37)/3;
    Double_t BkgM = Nbkg*(1+cheb0*shx + cheb1*(2*shx*shx-1) + cheb2*(4*shx*shx*shx-3*shx));

    return BkgM;
}


//totalvn pol3 bkg Psi2S{{{
Double_t Totalvnpol3JPsi(Double_t* x, Double_t* par)
{
    Double_t N1 = par[0];
    Double_t Nbkg = par[1];
    Double_t mean = par[2];
    Double_t sigma = par[3];
    Double_t alpha = par[4];
    Double_t n = par[5];
    Double_t ratio = par[6];
    Double_t frac = par[7];
    Double_t cheb0 = par[8];
    Double_t cheb1 = par[9];
    Double_t cheb2 = par[10];
    Double_t c = par[11];
    Double_t c1 = par[12];
    Double_t c2 = par[13];
    Double_t c3 = par[14];
    Double_t c4 = par[15];

    Double_t sigma1_2 = sigma*ratio;

    //t2 > t1
    Double_t JPsi_t1 = (x[0]-mean)/sigma;
    Double_t JPsi_t2 = (x[0]-mean)/sigma1_2;
    if (alpha < 0)
    {
        cout << "ERROR ::: alpha variable negative!!!! " << endl;
        return -1;
    }

    Double_t absAlpha = fabs((Double_t)alpha);
    Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
    Double_t b = n/absAlpha - absAlpha;

    Double_t JPsi_1 = -1;
    Double_t JPsi_2 = -1;

    if(JPsi_t1 > -alpha){
        JPsi_1 = exp(-JPsi_t1*JPsi_t1/2.);
    }
    else if(JPsi_t1 <= -alpha){
        JPsi_1 = a*TMath::Power((b-JPsi_t1),-n);
    }

    JPsi_2=1/(sigma1_2*TMath::Sqrt(2*TMath::Pi()))*exp(-0.5*JPsi_t2*JPsi_t2);
    Double_t fC = n/absAlpha*1/(n-1)*exp(-absAlpha*absAlpha/2);
    Double_t fD = TMath::Sqrt(TMath::Pi()/2)*(1+TMath::Erf(absAlpha/TMath::Sqrt(2)));
    Double_t fN_1 = 1./(sigma*(fC+fD));
    Double_t SigM = N1*(fN_1*frac*JPsi_1 + (1-frac)*JPsi_2);

    double shx = (10*x[0]-37)/3;
    Double_t BkgM = Nbkg*(1+cheb0*shx + cheb1*(2*shx*shx-1) + cheb2*(4*shx*shx*shx-3*shx));

    return c*(SigM/(SigM+BkgM)) + (1 - SigM/(SigM+BkgM))*(c4 + c3*x[0] + c2*x[0]*x[0] + c1*x[0]*x[0]*x[0]);

}
//}}}

Double_t vnPol1BkgAlpha(Double_t* x, Double_t* par)
{
    Double_t N1 = par[0];
    Double_t Nbkg = par[1];
    Double_t mean = par[2];
    Double_t sigma = par[3];
    Double_t alpha = par[4];
    Double_t n = par[5];
    Double_t ratio = par[6];
    Double_t frac = par[7];
    Double_t cheb0 = par[8];
    Double_t cheb1 = par[9];
    Double_t cheb2 = par[10];
    Double_t c1 = par[11];
    Double_t c2 = par[12];
    Double_t c3 = par[13];

    Double_t sigma1_2 = sigma*ratio;

    //t2 > t1
    Double_t JPsi_t1 = (x[0]-mean)/sigma;
    Double_t JPsi_t2 = (x[0]-mean)/sigma1_2;
    if (alpha < 0)
    {
        cout << "ERROR ::: alpha variable negative!!!! " << endl;
        return -1;
    }

    Double_t absAlpha = fabs((Double_t)alpha);
    Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
    Double_t b = n/absAlpha - absAlpha;

    Double_t JPsi_1 = -1;
    Double_t JPsi_2 = -1;

    if(JPsi_t1 > -alpha){
        JPsi_1 = exp(-JPsi_t1*JPsi_t1/2.);
    }
    else if(JPsi_t1 <= -alpha){
        JPsi_1 = a*TMath::Power((b-JPsi_t1),-n);
    }

    JPsi_2=1/(sigma1_2*TMath::Sqrt(2*TMath::Pi()))*exp(-0.5*JPsi_t2*JPsi_t2);
    Double_t fC = n/absAlpha*1/(n-1)*exp(-absAlpha*absAlpha/2);
    Double_t fD = TMath::Sqrt(TMath::Pi()/2)*(1+TMath::Erf(absAlpha/TMath::Sqrt(2)));
    Double_t fN_1 = 1./(sigma*(fC+fD));
    Double_t SigM = N1*(fN_1*frac*JPsi_1 + (1-frac)*JPsi_2);

    double shx = (10*x[0]-37)/3;
    Double_t BkgM = Nbkg*(1+cheb0*shx + cheb1*(2*shx*shx-1) + cheb2*(4*shx*shx*shx-3*shx));

    return (1 - SigM/(SigM+BkgM))*(c2 + c1*x[0]);

}

//totalvn pol2 bkg Psi2S * (1-alpha){{{
Double_t vnPol2BkgAlpha(Double_t* x, Double_t* par)
{
    Double_t N1 = par[0];
    Double_t Nbkg = par[1];
    Double_t mean = par[2];
    Double_t sigma = par[3];
    Double_t alpha = par[4];
    Double_t n = par[5];
    Double_t ratio = par[6];
    Double_t frac = par[7];
    Double_t cheb0 = par[8];
    Double_t cheb1 = par[9];
    Double_t cheb2 = par[10];
    Double_t c1 = par[11];
    Double_t c2 = par[12];
    Double_t c3 = par[13];

    Double_t sigma1_2 = sigma*ratio;

    //t2 > t1
    Double_t JPsi_t1 = (x[0]-mean)/sigma;
    Double_t JPsi_t2 = (x[0]-mean)/sigma1_2;
    if (alpha < 0)
    {
        cout << "ERROR ::: alpha variable negative!!!! " << endl;
        return -1;
    }

    Double_t absAlpha = fabs((Double_t)alpha);
    Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
    Double_t b = n/absAlpha - absAlpha;

    Double_t JPsi_1 = -1;
    Double_t JPsi_2 = -1;

    if(JPsi_t1 > -alpha){
        JPsi_1 = exp(-JPsi_t1*JPsi_t1/2.);
    }
    else if(JPsi_t1 <= -alpha){
        JPsi_1 = a*TMath::Power((b-JPsi_t1),-n);
    }

    JPsi_2=1/(sigma1_2*TMath::Sqrt(2*TMath::Pi()))*exp(-0.5*JPsi_t2*JPsi_t2);
    Double_t fC = n/absAlpha*1/(n-1)*exp(-absAlpha*absAlpha/2);
    Double_t fD = TMath::Sqrt(TMath::Pi()/2)*(1+TMath::Erf(absAlpha/TMath::Sqrt(2)));
    Double_t fN_1 = 1./(sigma*(fC+fD));
    Double_t SigM = N1*(fN_1*frac*JPsi_1 + (1-frac)*JPsi_2);

    double shx = (10*x[0]-37)/3;
    Double_t BkgM = Nbkg*(1+cheb0*shx + cheb1*(2*shx*shx-1) + cheb2*(4*shx*shx*shx-3*shx));

    return (1 - SigM/(SigM+BkgM))*(c3 + c2*x[0] + c1*x[0]*x[0]);

}
//}}}

//totalvn pol3 bkg Psi2S * (1-alpha){{{
Double_t vnPol3BkgAlpha(Double_t* x, Double_t* par)
{
    Double_t N1 = par[0];
    Double_t Nbkg = par[1];
    Double_t mean = par[2];
    Double_t sigma = par[3];
    Double_t alpha = par[4];
    Double_t n = par[5];
    Double_t ratio = par[6];
    Double_t frac = par[7];
    Double_t cheb0 = par[8];
    Double_t cheb1 = par[9];
    Double_t cheb2 = par[10];
    Double_t c1 = par[11];
    Double_t c2 = par[12];
    Double_t c3 = par[13];
    Double_t c4 = par[14];

    Double_t sigma1_2 = sigma*ratio;

    //t2 > t1
    Double_t JPsi_t1 = (x[0]-mean)/sigma;
    Double_t JPsi_t2 = (x[0]-mean)/sigma1_2;
    if (alpha < 0)
    {
        cout << "ERROR ::: alpha variable negative!!!! " << endl;
        return -1;
    }

    Double_t absAlpha = fabs((Double_t)alpha);
    Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
    Double_t b = n/absAlpha - absAlpha;

    Double_t JPsi_1 = -1;
    Double_t JPsi_2 = -1;

    if(JPsi_t1 > -alpha){
        JPsi_1 = exp(-JPsi_t1*JPsi_t1/2.);
    }
    else if(JPsi_t1 <= -alpha){
        JPsi_1 = a*TMath::Power((b-JPsi_t1),-n);
    }

    JPsi_2=1/(sigma1_2*TMath::Sqrt(2*TMath::Pi()))*exp(-0.5*JPsi_t2*JPsi_t2);
    Double_t fC = n/absAlpha*1/(n-1)*exp(-absAlpha*absAlpha/2);
    Double_t fD = TMath::Sqrt(TMath::Pi()/2)*(1+TMath::Erf(absAlpha/TMath::Sqrt(2)));
    Double_t fN_1 = 1./(sigma*(fC+fD));
    Double_t SigM = N1*(fN_1*frac*JPsi_1 + (1-frac)*JPsi_2);
    
    double shx = (10*x[0]-37)/3;
    Double_t BkgM = Nbkg*(1+cheb0*shx + cheb1*(2*shx*shx-1) + cheb2*(4*shx*shx*shx-3*shx));

    return (1 - SigM/(SigM+BkgM))*(c4 + c3*x[0] + c2*x[0]*x[0] + c1*x[0]*x[0]*x[0]);

}
//}}}

//Alpha Funct{{{
Double_t alphaFunct(Double_t* x, Double_t* par)
{
    Double_t N1 = par[0];
    Double_t Nbkg = par[1];
    Double_t mean = par[2];
    Double_t sigma = par[3];
    Double_t alpha = par[4];
    Double_t n = par[5];
    Double_t ratio = par[6];
    Double_t frac = par[7];
    Double_t cheb0 = par[8];
    Double_t cheb1 = par[9];
    Double_t cheb2 = par[10];
    Double_t sigma1_2 = sigma*ratio;

    //t2 > t1
    Double_t JPsi_t1 = (x[0]-mean)/sigma;
    Double_t JPsi_t2 = (x[0]-mean)/sigma1_2;
    if (alpha < 0)
    {

        cout << "ERROR ::: alpha variable negative!!!! " << endl;
        return -1;
    }

    Double_t absAlpha = fabs((Double_t)alpha);
    Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
    Double_t b = n/absAlpha - absAlpha;

    Double_t JPsi_1 = -1;
    Double_t JPsi_2 = -1;

    if(JPsi_t1 > -alpha){
        JPsi_1 = exp(-JPsi_t1*JPsi_t1/2.);
    }
    else if(JPsi_t1 <= -alpha){
        JPsi_1 = a*TMath::Power((b-JPsi_t1),-n);
    }

    JPsi_2=1/(sigma1_2*TMath::Sqrt(2*TMath::Pi()))*exp(-0.5*JPsi_t2*JPsi_t2);
    Double_t fC = n/absAlpha*1/(n-1)*exp(-absAlpha*absAlpha/2);
    Double_t fD = TMath::Sqrt(TMath::Pi()/2)*(1+TMath::Erf(absAlpha/TMath::Sqrt(2)));
    Double_t fN_1 = 1./(sigma*(fC+fD));
    Double_t SigM = N1*(fN_1*frac*JPsi_1 + (1-frac)*JPsi_2);
    
    double shx = (10*x[0]-37)/3;
    Double_t BkgM = Nbkg*(1+cheb0*shx + cheb1*(2*shx*shx-1) + cheb2*(4*shx*shx*shx-3*shx));

    return (SigM/(SigM+BkgM));

}
//}}}

//pol1 bkg{{{
Double_t pol1bkg(Double_t* x, Double_t* par)
{
    Double_t c1 = par[0];
    Double_t c2 = par[1];

    return c1*x[0]+c2;
}
//}}}


//pol2 bkg{{{
Double_t pol2bkg(Double_t* x, Double_t* par)
{
    Double_t c1 = par[0];
    Double_t c2 = par[1];
    Double_t c3 = par[2];

    return c1*x[0]*x[0]+c2*x[0]+c3;
}
//}}}

//pol2 bkg{{{
Double_t pol3bkg(Double_t* x, Double_t* par)
{
    Double_t c1 = par[0];
    Double_t c2 = par[1];
    Double_t c3 = par[2];
    Double_t c4 = par[3];

    return c1*x[0]*x[0]*x[0]+c2*x[0]*x[0]+c3*x[0] + c4;
}
//}}}

void doSim_pt65_10_y0_24_cent20_120(int cLow = 20, int cHigh = 120,
		float ptLow =  6.5, float ptHigh = 10.0,
		float yLow = 0.0, float yHigh = 2.4,
		float SiMuPtCut = 0, float massLow = 3.3, float massHigh =4.1, bool dimusign=true, 
		int ibkg_vn_sel = fpol1, bool fixSigPar=true)
{
    TString DATE = "210928";
    gSystem->mkdir(Form("roots/%s",DATE.Data()),kTRUE);
    gSystem->mkdir(Form("figs/%s",DATE.Data()),kTRUE);
    
	setTDRStyle();
	gStyle->SetOptFit(0000);
	writeExtraText= true;
	int iPeriod = 2;
	int iPos = 33;

	TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, SiMuPtCut, cLow, cHigh) ;
	TString dimusignString;
	if(dimusign) dimusignString = "OS";
	else if(!dimusign) dimusignString = "SS";

	kineLabel = kineLabel + Form("_m%.1f-%.1f",massLow,massHigh) + "_" + dimusignString;

	TF1* fyieldtot;
	TF1* fvntot;
	Double_t v2;
	Double_t v2e;
	Double_t v2_bkg;

	int nParmV_;
	int nParBkg;
	if(ibkg_vn_sel == fpol2) {nParmV_ = 15; nParBkg = 3;} 
	else if(ibkg_vn_sel == fpol3) {nParmV_ = 16; nParBkg = 4;}
	else if(ibkg_vn_sel == fpol1) {nParmV_ = 14, nParBkg = 2;}
	else{
		cout << "ERROR ::: No Selection for v2 background function!!!!" << endl;
		return;
	}

	//Get yield distribution{{{
	//TFile* rf = new TFile(Form("roots/v2mass_hist/Psi2S_NonPrompt_%s.root",kineLabel.Data()),"read");
    TFile* rf = new TFile(Form("../make_v2_mass_hist/roots/Psi2S_Prompt_%s_Eff1_Acc1_PtW1_TnP1.root",kineLabel.Data()),"read");
	TH1D* h_v2_SplusB = (TH1D*) rf->Get("h_v2_SplusB");
	TGraphAsymmErrors* g_mass = (TGraphAsymmErrors*) rf->Get("g_mass");

    TFile *wf = new TFile(Form("roots/%s/SimFitResult_Prompt_%s.root",DATE.Data(),kineLabel.Data()),"recreate");

	//define function for simultaneous fitting{{{
	TF1* fmass_total = new TF1("fmass_total", TotalYield, massLow, massHigh, nParmM);
	TF1* fvn_simul;
	if(ibkg_vn_sel == fpol2) fvn_simul = new TF1("fvn_simul", Totalvnpol2JPsi, massLow, massHigh, nParmV_);
	else if(ibkg_vn_sel == fpol3) fvn_simul = new TF1("fvn_simul", Totalvnpol3JPsi, massLow, massHigh, nParmV_);
	else if(ibkg_vn_sel == fpol1) fvn_simul = new TF1("fvn_simul", Totalvnpol1JPsi, massLow, massHigh, nParmV_);
	//}}}

	//combine functions{{{
	fmass_total->SetLineColor(2);
	fmass_total->SetLineWidth(1);

	fvn_simul->SetLineColor(2);
	fvn_simul->SetLineWidth(1);

	ROOT::Math::WrappedMultiTF1 wmass(*fmass_total, 1);
	ROOT::Math::WrappedMultiTF1 wvn(*fvn_simul, 1);

	ROOT::Fit::DataOptions opt;
	ROOT::Fit::DataRange massrange;

	massrange.SetRange(massLow,massHigh);
	ROOT::Fit::BinData datamass(opt, massrange);
	ROOT::Fit::FillData(datamass, g_mass);

	ROOT::Fit::DataRange vnrange;
	vnrange.SetRange(massLow,massHigh);
	ROOT::Fit::BinData datavn(opt, vnrange);
	ROOT::Fit::FillData(datavn, h_v2_SplusB);

	ROOT::Fit::Chi2Function mass_chi2(datamass, wmass);
	ROOT::Fit::Chi2Function vn_chi2(datavn, wvn);

	GlobalChi2_width globalChi2(mass_chi2, vn_chi2);

	ROOT::Fit::Fitter fitter;
	//}}}

	TString kineLabel_ = getKineLabel (ptLow, ptHigh, yLow, yHigh, SiMuPtCut, cLow, cHigh) ;
    TFile* f_mass = new TFile(Form("../MassFit/data_fit/roots/MassFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1.root",kineLabel_.Data()),"read");
	RooWorkspace *ws = new RooWorkspace("workspace");
	RooDataSet *datasetMass = (RooDataSet*)f_mass->Get("datasetMass");
	ws->import(*datasetMass);
	f_mass->cd();

	Double_t N1_ = ws->var("N_Jpsi")->getVal();
	Double_t Nbkg_ = ws->var("N_Bkg")->getVal();
	Double_t mean_ = pdgMass.Psi2S;
	Double_t sigma_ = ws->var("sigma_1_A")->getVal();
	Double_t alpha_ = ws->var("alpha_1_A")->getVal();
	Double_t n_ = ws->var("n_1_A")->getVal();
	Double_t ratio_ = ws->var("x_A")->getVal();
	Double_t frac_ = ws->var("f")->getVal();
    Double_t cheb0_ = ws->var("sl1")->getVal();
    Double_t cheb1_ = ws->var("sl2")->getVal();
    Double_t cheb2_ = ws->var("sl3")->getVal();
	Double_t c_  = 0.190;
	Double_t c1_ = 0.6210;
	Double_t c2_ = 0.3706;



	Double_t c3_ = 4.0410;
	Double_t c4_ = -0.0010;

	// Double_t cheb0_ = 0.0121;
	// Double_t cheb1_ = 0.0135;
	// Double_t cheb2_ = 0.4226;
	// Double_t c_  = 0.00921;
	// Double_t c1_ = 9.1210;
	// Double_t c2_ = 3.4006;

	// Double_t cheb0_ = 2.0121;
	// Double_t cheb1_ = 2.0135;
	// Double_t cheb2_ = 1.0226;
	// Double_t c_  = 0.35021;
	// Double_t c1_ = 2.1210;
	// Double_t c2_ = 3.4006;
	//

	std::cout << "----- OK? ------" << std::endl;
	Double_t par0[nParmV];
	par0[0] = N1_;
	par0[1] = Nbkg_;
	par0[2] = mean_;
	par0[3] = sigma_;
	par0[4] = alpha_;
	par0[5] = n_;
	par0[6] = ratio_;
	par0[7] = frac_;
	par0[8] = cheb0_;
	par0[9] = cheb1_;
	par0[10] = cheb2_;
	par0[11] = c_;
	par0[12] = c1_;
	par0[13] = c2_;
	par0[14] = c3_;
	par0[15] = c4_;


    Double_t parLimitLow[nParmV]  = {    0,       0, mean_-0.005,     0,   0.,   0.,    0,     0,  -10, -10, -10,  0, -20, -20, -20,-20};
    Double_t parLimitHigh[nParmV] = {N1_*1.2, Nbkg_*1.2, mean_+0.005,      0.4,    5.,   5.,  5.,  1.,   10,  10,  10, 0.3,  20,  20,  20, 20};


	fitter.Config().SetParamsSettings(nParmV_, par0);
	for(int ipar = 0; ipar<nParmV_; ipar++){
		fitter.Config().ParSettings(ipar).SetLimits(parLimitLow[ipar],parLimitHigh[ipar]);
	}

	int parmid_sigi = 3;
	int parmid_sigf = 7;
	if(fixSigPar){
		for(int iparmsig=parmid_sigi; iparmsig<=parmid_sigf;iparmsig++) fitter.Config().ParSettings(iparmsig).Fix();
	}
	fitter.Config().MinimizerOptions().SetPrintLevel(0);
	fitter.Config().SetMinimizer("Minuit2","Migrad");
	//}}}

	fitter.FitFCN(nParmV_, globalChi2, 0, datamass.Size()+datavn.Size(), true);
	ROOT::Fit::FitResult result = fitter.Result();
	result.Print(std::cout);


	cout << "****************************************" << endl;
	cout << "	" << "	Upper Limit" << "	" << "Lower Limit" << endl;
	// Print Parameter and Parameter Limits
	for(int i=0; i<nParmV_; i++){
		cout << "par[" << i << "]	:	" << parLimitHigh[i] << "	" <<"	" << parLimitLow[i] << endl;
	}
	cout << " " << endl;

	//Yield fitting result{{{
	fmass_total->SetFitResult(result, iparmass);
	fmass_total->SetRange(massrange().first, massrange().second);
	g_mass->GetListOfFunctions()->Add(fmass_total);
	h_v2_SplusB->GetListOfFunctions()->Add(fvn_simul);


	int prmid_bkgyield = 1;
	int prmid_bkgcheb0   = 8;
	int prmid_bkgcheb1   = 9;
	int prmid_bkgcheb2   = 10;
	int prmid_vnbkg1   = 12;
	int prmid_vnbkg2   = 13;
	int prmid_vnbkg3   = 14;
	int prmid_vnbkg4   = 15;
	int prmid_v2val    = 11;
	int nprm_sigf      = 8;
	int nprm_bkgf      = 4;
	int nprm_alpha     = 11;
	double massYMin=0;//0
	// double massYMax=7000;//1000

	TF1* fyield_bkg = new TF1("fyield_bkg", TotalYieldBkg, massLow, massHigh,nprm_bkgf);
	fyield_bkg->FixParameter(0, fmass_total->GetParameter(prmid_bkgyield));
	fyield_bkg->FixParameter(1, fmass_total->GetParameter(prmid_bkgcheb0));
	fyield_bkg->FixParameter(2, fmass_total->GetParameter(prmid_bkgcheb1));
	fyield_bkg->FixParameter(3, fmass_total->GetParameter(prmid_bkgcheb2));
	g_mass->GetListOfFunctions()->Add(fyield_bkg);
	//}}}

	//vn fitting result{{{
	fvn_simul->SetFitResult(result, iparvn);
	fvn_simul->SetRange(vnrange().first, vnrange().second);

	TF1* fyield_sig = new TF1("fyield_sig",TotalYieldSig, massLow, massHigh, nprm_sigf);
	for(int iparm=0;iparm<nprm_sigf-1; iparm++){
		if(iparm==0) fyield_sig->FixParameter(iparm, fvn_simul->GetParameter(iparm));
		else if(iparm!=0) fyield_sig->FixParameter(iparm, fvn_simul->GetParameter(iparm+1));
	}
	fyield_sig->FixParameter(nprm_sigf-1,massYMin);
	g_mass->GetListOfFunctions()->Add(fyield_sig);
	//}}}

	//Alpha function{{{
	TF1* fAlpha = new TF1("fAlpha",alphaFunct,massLow,massHigh,nprm_alpha);
	for(int iparm=0;iparm<nprm_alpha; iparm++){
		fAlpha->FixParameter(iparm,fvn_simul->GetParameter(iparm));
	}
	//h_v2_SplusB->GetListOfFunctions()->Add(fAlpha);
	//}}}


	//vn bkg {{{
	TF1* fvn_bkg;
	TF1* fvn_bkg_alpha;
	if(ibkg_vn_sel == fpol2){
		fvn_bkg = new TF1("fvn_bkg",pol2bkg, massLow, massHigh, nParBkg);
		fvn_bkg->FixParameter(0, fvn_simul->GetParameter(prmid_vnbkg1));
		fvn_bkg->FixParameter(1, fvn_simul->GetParameter(prmid_vnbkg2));
		fvn_bkg->FixParameter(2, fvn_simul->GetParameter(prmid_vnbkg3));
		fvn_bkg_alpha = new TF1("fvn_bkg_alpha",vnPol2BkgAlpha,massLow,massHigh,nParmV_ -1);
		for(int iparm=0;iparm<nParmV_ -1; iparm++){
			if(iparm<=prmid_v2val) fvn_bkg_alpha->FixParameter(iparm,fvn_simul->GetParameter(iparm));
			else if(iparm>=prmid_v2val) fvn_bkg_alpha->FixParameter(iparm,fvn_simul->GetParameter(iparm+1));
		}
	}
	else if(ibkg_vn_sel == fpol3){
		fvn_bkg = new TF1("fvn_bkg",pol3bkg, massLow, massHigh, nParBkg);
		fvn_bkg->FixParameter(0, fvn_simul->GetParameter(prmid_vnbkg1));
		fvn_bkg->FixParameter(1, fvn_simul->GetParameter(prmid_vnbkg2));
		fvn_bkg->FixParameter(2, fvn_simul->GetParameter(prmid_vnbkg3));
		fvn_bkg->FixParameter(3, fvn_simul->GetParameter(prmid_vnbkg4));
		fvn_bkg_alpha = new TF1("fvn_bkg_alpha",vnPol3BkgAlpha,massLow,massHigh,nParmV_ -1);
		for(int iparm=0;iparm<nParmV_ -1; iparm++){
			if(iparm<=prmid_v2val) fvn_bkg_alpha->FixParameter(iparm,fvn_simul->GetParameter(iparm));
			else if(iparm>=prmid_v2val) fvn_bkg_alpha->FixParameter(iparm,fvn_simul->GetParameter(iparm+1));
		}
	}
	else if(ibkg_vn_sel == fpol1){
		fvn_bkg = new TF1("fvn_bkg",pol1bkg, massLow, massHigh, nParBkg);
		fvn_bkg->FixParameter(0, fvn_simul->GetParameter(prmid_vnbkg1));
		fvn_bkg->FixParameter(1, fvn_simul->GetParameter(prmid_vnbkg2));
		fvn_bkg->FixParameter(2, fvn_simul->GetParameter(prmid_vnbkg3));
		fvn_bkg_alpha = new TF1("fvn_bkg_alpha",vnPol1BkgAlpha,massLow,massHigh,nParmV_ -1);
		for(int iparm=0;iparm<nParmV_ -1; iparm++){
			if(iparm<=prmid_v2val) fvn_bkg_alpha->FixParameter(iparm,fvn_simul->GetParameter(iparm));
			else if(iparm>=prmid_v2val) fvn_bkg_alpha->FixParameter(iparm,fvn_simul->GetParameter(iparm+1));
		}
	}


	h_v2_SplusB->GetListOfFunctions()->Add(fvn_bkg);
	//}}}

	unsigned int nfpxl = 2000;
	fvn_simul->SetNpx(nfpxl);
	fmass_total->SetNpx(nfpxl);
	fAlpha->SetNpx(nfpxl);
	fvn_bkg_alpha->SetNpx(nfpxl);

	v2 = fvn_simul->GetParameter(prmid_v2val);
	v2e = fvn_simul->GetParError(prmid_v2val);

	// Drawing 
	fyield_bkg->SetLineColor(kBlue);
	fyield_bkg->SetLineStyle(kDashed);
	fyield_bkg->SetLineWidth(3);
	fyield_sig->SetFillStyle(3004);
	fyield_sig->SetFillColor(kRed);
	fmass_total->SetLineColor(kBlack);
	fmass_total->SetLineWidth(3);
	fvn_simul->SetLineColor(kRed+2);
	fvn_simul->SetLineWidth(3);
	fvn_bkg->SetLineColor(kRed+2);
	fvn_bkg->SetLineWidth(2);
	fvn_bkg->SetLineStyle(kDashed);

	fvn_bkg_alpha->SetLineStyle(4);
	fvn_bkg_alpha->SetLineColor(kOrange+2);
	fvn_bkg_alpha->SetLineWidth(2);
	fAlpha->SetLineStyle(6);
	fAlpha->SetLineColor(kGreen+2);
	fAlpha->SetLineWidth(2);

	h_v2_SplusB->GetYaxis()->SetRangeUser(0.05,0.26);
	h_v2_SplusB->GetYaxis()->SetTitle("v_{2}^{S+B}");
	h_v2_SplusB->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
	h_v2_SplusB->GetYaxis()->SetLabelSize(0.055);
	h_v2_SplusB->GetXaxis()->SetLabelSize(0.055);
	h_v2_SplusB->GetXaxis()->SetTitleSize(0.07);
	h_v2_SplusB->GetYaxis()->SetTitleSize(0.07);
	h_v2_SplusB->GetYaxis()->SetTitleOffset(1.2);
	h_v2_SplusB->GetXaxis()->SetTitleOffset(1.0);
	h_v2_SplusB->GetXaxis()->SetLabelOffset(0.011);
	h_v2_SplusB->GetXaxis()->SetRangeUser(massLow,massHigh);
	h_v2_SplusB->GetXaxis()->SetLimits(massLow,massHigh);
	SetHistStyle(h_v2_SplusB,0,0);
	SetGraphStyle2(g_mass,0,0);

	g_mass->SetMarkerSize(1);
	//g_mass->SetMinimum(massYMin);
	//g_mass->SetMaximum(massYMax);
	g_mass->GetXaxis()->SetLimits(massLow,massHigh);
	g_mass->GetXaxis()->SetRangeUser(massLow,massHigh);
	g_mass->GetYaxis()->SetTitleOffset(1.7);
	g_mass->GetYaxis()->SetTitle("Events/(0.02GeV/c^{2})");
	g_mass->GetYaxis()->SetLabelSize(0.055);
	g_mass->GetXaxis()->SetLabelSize(0.055);
	g_mass->GetXaxis()->SetTitleSize(0.07);
	g_mass->GetYaxis()->SetTitleSize(0.07);
	g_mass->GetYaxis()->SetTitleOffset(1.2);
	g_mass->GetXaxis()->SetNdivisions(510);

	double sizeTick = 12;

	float pos_x = 0.03;
	float pos_x_mass = 0.23;
	float pos_y = 0.76;
	float pos_y_diff = 0.071;
	int text_color = 1;
	float text_size = 16;
	TString perc = "%";

	TLegend *leg1 = new TLegend(0.73,0.48,0.98,0.70);
	SetLegendStyle(leg1);
	leg1->SetTextSize(0.05);
	leg1->AddEntry(g_mass,"Data","lp");
	leg1->AddEntry(fmass_total,"Fit","l");
	leg1->AddEntry(fyield_sig,"#psi(2S) Signal","lf");
	leg1->AddEntry(fyield_bkg,"Background","lf");

	TCanvas* c_mass_v2 = new TCanvas("c_mass_v2","",590,750);
	TPad* pad1 = new TPad("pad1","pad1",0,0.5,1.0,1.0);
	TPad* pad2 = new TPad("pad2","pad2",0,0.0,1.0,0.5);
	c_mass_v2->cd();
	pad1->SetTicks(1,1);
	pad1->SetBottomMargin(0);
	pad1->SetLeftMargin(0.19);
	pad1->SetTopMargin(0.08);
	pad1->Draw();
	pad1->cd();
	// pad1->SetLogy(1);
	double pad1W = pad1->GetWw()*pad1->GetAbsWNDC();
	double pad1H = pad1->GetWh()*pad1->GetAbsHNDC();
	double tickScaleX = (pad1->GetUxmax() - pad1->GetUxmin())/(pad1->GetX2()-pad1->GetX1())*pad1H;
	g_mass->GetXaxis()->SetTickLength(sizeTick/tickScaleX);   
	g_mass->GetXaxis()->SetTitleSize(0);
	g_mass->GetXaxis()->SetLabelSize(0);
	g_mass->GetXaxis()->SetNdivisions(510);
	g_mass->Draw("AP");
	leg1->Draw("same");
	//fyield_sig->Draw("same");
	if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_x_mass,pos_y,text_color,text_size);
	else if(ptLow!=0) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
	if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
	else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
	drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
	drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
	drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);

	pad1->Update();

	c_mass_v2->cd();
	pad2->SetTicks(1,1);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.17);
	pad2->SetLeftMargin(0.19);
	pad2->Draw();
	pad2->cd();
	double pad2W = pad2->GetWw()*pad2->GetAbsWNDC();
	double pad2H = pad2->GetWh()*pad2->GetAbsHNDC();
	h_v2_SplusB->Draw("P");
	pad1->Update();
	pad2->Update();
	h_v2_SplusB->GetXaxis()->SetTickLength(sizeTick/tickScaleX);
	h_v2_SplusB->Draw("P");
	jumSun(massLow,0,massHigh,0,1,1);
	drawText(Form("v_{2}^{S} = %.3f #pm %.3f",v2,v2e),pos_x_mass,pos_y,text_color,text_size);

	TLegend *leg2 = new TLegend(0.73,0.6,0.95,0.8);
	SetLegendStyle(leg2);
	leg2->SetTextSize(0.05);
	leg2->AddEntry(fvn_simul,"v_{2}^{S+B}","l");
	leg2->AddEntry(fvn_bkg,"v_{2}^{B}","l");
	//leg2->AddEntry(fAlpha,"#alpha","l");
	leg2->Draw("same");

	CMS_lumi_v2mass(pad1,iPeriod,iPos);  
	pad1->Update();
	pad2->Update();
	c_mass_v2->cd();
	pad1->Draw();
	pad2->Draw();
	c_mass_v2->Update();
	c_mass_v2->SaveAs(Form("figs/%s/v2Mass_Prompt_%s.pdf",DATE.Data(),kineLabel.Data()));
	wf->cd();
	//store individual function{{{
	fyieldtot = (TF1*) fmass_total->Clone();
	fyieldtot->SetName("massfit");
	fyieldtot->Write();

	fvntot = (TF1*) fvn_simul->Clone();
	fvntot->SetName("vnfit");
	fvntot->Write();

	fyield_bkg->Write();
	fyield_sig->Write();
	fvn_bkg->Write();
	fAlpha->Write();
	fvn_bkg_alpha->Write();
	//}}}

	//get Chi2{{{
	Double_t ptar = ptHigh+ptLow/2;
	TGraphErrors* v2plot = new TGraphErrors();
	v2plot->SetPoint(0,ptar,v2);
	v2plot->SetPointError(0,0,v2e);
	v2plot->SetTitle("");
	v2plot->SetName("v2vspt");
	v2plot->Write();
}


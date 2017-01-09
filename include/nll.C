#include "TMath.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "TROOT.h"

#include <iostream>
using namespace std;


RooWorkspace*  f_w = 0;

RooAbsPdf* zptpdf = 0;
RooAbsPdf* metpdf = 0;
RooAbsPdf* ldppdf = 0;
RooAbsPdf* mlbpdf = 0;
		 
RooArgSet* obszpt = 0;
RooArgSet* obsmet = 0;
RooArgSet* obsldp = 0;
RooArgSet* obsmlb = 0;

void SetW(RooWorkspace* w){ f_w = w; }

void Setzptpdf(RooAbsPdf* pdf){ zptpdf = pdf; }
void Setmetpdf(RooAbsPdf* pdf){	metpdf = pdf; }
void Setldppdf(RooAbsPdf* pdf){	ldppdf = pdf; }
void Setmlbpdf(RooAbsPdf* pdf){ mlbpdf = pdf; }

void Setobszpt(RooArgSet*obs){ obszpt = obs; }
void Setobsmet(RooArgSet*obs){ obsmet = obs; }
void Setobsldp(RooArgSet*obs){ obsldp = obs; }
void Setobsmlb(RooArgSet*obs){ obsmlb = obs; }



float nll(float met, float zpt, float mlb, float ldp)
{


  f_w->var("lepsZPt_Edge")->setVal(zpt);
  float zptPdfVal = zptpdf->getVal(obszpt);
  
  f_w->var("met_Edge")->setVal(met);
  float metPdfVal = metpdf->getVal(obsmet);
  
  f_w->var("lepsDPhi_Edge")->setVal(ldp);
  float ldpPdfVal = ldppdf->getVal(obsldp);

  f_w->var("sum_mlb_Edge")->setVal(mlb);
  float mlbPdfVal = mlbpdf->getVal(obsmlb);

  return -1.*TMath::Log(zptPdfVal*metPdfVal*ldpPdfVal*mlbPdfVal);
  
}

#include "TMath.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "TROOT.h"
#include "TVector3.h"

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


float deltaPhi(float phi1, float phi2)
{
  float result = phi2-phi1;
  while (result > TMath::Pi()) result -= 2*TMath::Pi();
  while (result < -TMath::Pi())result += 2*TMath::Pi();
  return result;
}

float deltaPhiZZ_Jet(float ptl1, float ptl2, float phil1, float phil2, float met_pt, float phiMET, float phiJet, float ptJet)
{
  TVector3 l1, l2, met, jet;
  l1 .SetPtEtaPhi(ptl1, 0., phil1);
  l2 .SetPtEtaPhi(ptl2, 0., phil2);
  met.SetPtEtaPhi(met_pt,  0., phiMET);
  jet.SetPtEtaPhi(ptJet, 0., phiJet); // ptJet not really needed but ok 

  return (l1+l2+met).DeltaPhi(jet);
}

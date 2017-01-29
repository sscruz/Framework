#include<iostream>
#include"TH1F.h"

using namespace std;


TH1F* puHist = 0;
TH1F* puHistUp = 0;
TH1F* puHistDn = 0;

TH1D* OldPUHist = 0;
void SetPUhist(TH1F* h)  { puHist = h; }
void SetPUhistUp(TH1F* h){ puHistUp = h; }
void SetPUhistDn(TH1F* h){ puHistDn = h; }
void SetOldPUHist(TH1D* h){ OldPUHist = h; }

// float PUWeight( float ntrue){
//   return puHist->GetBinContent( puHist->FindBin(ntrue)) 
// }

// float PUWeightUp( float ntrue){
//   return puHistUp->GetBinContent( puHistUp->FindBin(ntrue)) 
// }
// float PUWeightDn( float ntrue){
//   return puHistDn->GetBinContent( puHistDn->FindBin(ntrue)) 
// }


float PUWeight( float oldweight){
  float ntrue = -1;
  for (int bin = 0; bin < OldPUHist->GetNbinsX()+1; ++bin){
    if (TMath::Abs( OldPUHist->GetBinContent(bin+1) - oldweight ) < 1e-4){
      ntrue = OldPUHist->GetBinCenter(bin+1);
    }
  }
  if (ntrue < 0){
    cout << "ntrue int not found" << endl;
    for (int bin = 0; bin < OldPUHist->GetNbinsX()+1; ++bin){
      cout << oldweight - OldPUHist->GetBinContent(bin+1) << endl;
    }
  }
  return puHist->GetBinContent( puHist->FindBin(ntrue));
}

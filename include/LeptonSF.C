#include<iostream>
#include"TH1.h"
#include"TH2.h"
#include"TGraphAsymmErrors.h"

using namespace std;

vector<TH2*> hElec;
vector<TH2*> hMuon;


void AddElec (TH2* h){ hElec.push_back(h); }
void AddMuon (TH2* h){ hMuon.push_back(h); }



Double_t LepSF(Double_t pt, Double_t eta, Int_t pdgId, TString sys="")
{
  // sys can be 
  // "": nominal
  // "ElUp": obvious
  // "ElDn": obvious
  // "MuUp": obvious
  // "MuDn": obvious
  // cout << "[LepSF]" << pt << " " << eta << " " << pdgId << " " << sys << endl;

  vector<TH2*> hists = abs(pdgId)==11 ? hElec : hMuon;
  float out = 1;
  int var = 0;

  if (sys.Contains("Up")) var = 1;
  else if (sys.Contains("Dn")) var = 1;
  else var = 0;

  for (auto hist : hists){
    int ptbin  = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt)));
    int etabin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(fabs(eta))));
    out *=  hist->GetBinContent(ptbin,etabin)+var*hist->GetBinError(ptbin,etabin);
    if (abs(pdgId) == 13)
      out = out+var*out*0.03;
  }
  return out;
}


// Double_t LepSFElUp(Double_t pt, Double_t eta, Int_t pdgId)
// {
//   return LepSF(pt, eta, pdgId, "ElUp");
// }
// Double_t LepSFElDn(Double_t pt, Double_t eta, Int_t pdgId)
// {
//   return LepSF(pt, eta, pdgId, "ElDn");
// }
// Double_t LepSFMuUp(Double_t pt, Double_t eta, Int_t pdgId)
// {
//   return LepSF(pt, eta, pdgId, "MuUp");
// }
// Double_t LepSFMuDn(Double_t pt, Double_t eta, Int_t pdgId)
// {
//   return LepSF(pt, eta, pdgId, "MuDn");
// }


int njnb(int j, int b){
  if (j == 0 && b == 0) return 0;
  if (j == 1 && b == 0) return 1;
  if (j == 1 && b == 1) return 2;  
  if (j == 2 && b == 0) return 3;
  if (j == 2 && b == 1) return 4;
  if (j == 2 && b == 2) return 5;
  else return 6;
}

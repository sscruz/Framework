#include<iostream>
#include"TH1.h"
#include"TH2.h"


using namespace std;

TH2D* hElecID  = 0;
TH2D* hElecIP  = 0;
TH2D* hElecISO = 0;
TH2D* hElecTrk = 0;
TH2D* hMuonID  = 0;
TH2D* hMuonIP1  = 0;
TH2D* hMuonIP2  = 0;
TH2D* hMuonISO = 0;
TH2D* hMuonTrk = 0;


void SetFSElecID (TH2D* h){ hElecID  = h; }
void SetFSElecIP (TH2D* h){ hElecIP  = h; }
void SetFSElecISO(TH2D* h){ hElecISO = h; }
void SetFSElecTrk(TH2D* h){ hElecTrk = h; }
void SetFSMuonID (TH2D* h){ hMuonID  = h; }
void SetFSMuonIP1 (TH2D* h){ hMuonIP1  = h; }
void SetFSMuonIP2 (TH2D* h){ hMuonIP2  = h; }
void SetFSMuonISO(TH2D* h){ hMuonISO = h; }
void SetFSMuonTrk(TH2D* h){ hMuonTrk = h; }

Double_t _getFSIDoverRECO(Double_t pt, Double_t eta, Int_t pdgId, TString sys)
{

  eta = TMath::Abs(eta);
  Double_t sf   = 0.;
  Double_t sf_e = 0.;
  if (TMath::Abs(pdgId) == 13){
    if (pt > 200) pt = 199.9;
    sf   = hMuonID->GetBinContent( hMuonID->FindBin(pt,eta));
    sf_e = 0;
  }
  else{
    if (pt > 200) pt = 199.9;
    sf   = hElecID->GetBinContent( hElecID->FindBin(pt,eta));
    sf_e = 0.;
  }
  // error hardcoded to 2 later on
  // sf_e = (sys.Contains("Mu")) ? hMuonID->GetBinError  ( hMuonID->FindBin(pt,eta)) : 0.;

  if (sys.Contains("Up"))     return sf+sf_e;
  else if(sys.Contains("Dn")) return sf-sf_e;
  else return sf;
}


Double_t _getFSIPoverID(Double_t pt, Double_t eta, Int_t pdgId, TString sys)
{
  eta = TMath::Abs(eta);
  Double_t sf   = 0.;
  Double_t sf_e = 0.;
  if (TMath::Abs(pdgId) == 13){
    if (pt > 200) pt = 199.9;
    sf    = hMuonIP1->GetBinContent( hMuonIP1->FindBin(pt,eta));
    sf   *= hMuonIP2->GetBinContent( hMuonIP2->FindBin(pt,eta));
    sf_e = 0;
  }
  else{
    if (pt > 200) pt = 199.9;
    sf   = hElecIP->GetBinContent( hElecIP->FindBin(pt,eta));
    sf_e = 0;
  }
  // error hardcoded to 2 later on

  if (sys.Contains("Up"))     return sf+sf_e;
  else if(sys.Contains("Dn")) return sf-sf_e;
  else return sf;
}


Double_t _getFSISOoverIP(Double_t pt, Double_t eta, Int_t pdgId, TString sys)
{
  eta = TMath::Abs(eta);
  Double_t sf   = 0.;
  Double_t sf_e = 0.;
  if (TMath::Abs(pdgId) == 13){
    if (pt > 200) pt = 199.9;
    sf   = hMuonISO->GetBinContent( hMuonISO->FindBin(pt,eta));
    sf_e = 0;
  }
  else{
    if (pt > 200) pt = 199.9;
    sf   = hElecISO->GetBinContent( hElecISO->FindBin(pt,eta));
    sf_e = 0.;
  }
  // error hardcoded to 3 later on

  if (sys.Contains("Up"))     return sf+sf_e;
  else if(sys.Contains("Dn")) return sf-sf_e;
  else return sf;
}


Double_t LepSFFastSim(Double_t pt, Double_t eta, Int_t pdgId, TString sys="")
{
  // sys can be 
  // "": nominal
  // "ElUp": obvious
  // "ElDn": obvious
  // "MuUp": obvious
  // "MuDn": obvious
  // cout << "[LepSFFastSim]" << pt << " " << eta << " " << pdgId << " " << sys << endl;

  Double_t sf = _getFSIDoverRECO(pt,eta,pdgId,sys)*_getFSIPoverID(pt,eta,pdgId,sys)*_getFSISOoverIP(pt,eta,pdgId,sys);
  if (TMath::Abs(pdgId) == 13){
    if (sys.Contains("MuUp")){return TMath::Max(1e-100,sf*1.02);}
    else if (sys.Contains("MuDn")){return TMath::Max(1e-100,sf*0.98);}
    else return TMath::Max(1e-100,sf);
  }
  else{
    if (sys.Contains("ElUp")){return TMath::Max(1e-100,sf*1.02);}
    else if (sys.Contains("ElDn")){return TMath::Max(1e-100,sf*0.98);}
    else return TMath::Max(1e-100,sf);
  }

}


Double_t LepSFFastSimElUp(Double_t pt, Double_t eta, Int_t pdgId)
{
  return LepSFFastSim(pt, eta, pdgId, "ElUp");
}
Double_t LepSFFastSimElDn(Double_t pt, Double_t eta, Int_t pdgId)
{
  return LepSFFastSim(pt, eta, pdgId, "ElDn");
}
Double_t LepSFFastSimMuUp(Double_t pt, Double_t eta, Int_t pdgId)
{
  return LepSFFastSim(pt, eta, pdgId, "MuUp");
}
Double_t LepSFFastSimMuDn(Double_t pt, Double_t eta, Int_t pdgId)
{
  return LepSFFastSim(pt, eta, pdgId, "MuDn");
}



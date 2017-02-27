#include<iostream>
#include"TH1.h"
#include"TH2.h"


using namespace std;

TH2D* hSFElecID  = 0;
TH2D* hSFElecIP  = 0;
TH2D* hSFElecISO = 0;
TH2D* hSFElecTrk = 0;
TH2D* hSFMuonID  = 0;
TH2D* hSFMuonIP1  = 0;
TH2D* hSFMuonIP2  = 0;
TH2D* hSFMuonISO = 0;
TH2D* hSFMuonTrk = 0;


void SetFSElecID (TH2D* h){ hSFElecID  = h; }
void SetFSElecIP (TH2D* h){ hSFElecIP  = h; }
void SetFSElecISO(TH2D* h){ hSFElecISO = h; }
void SetFSElecTrk(TH2D* h){ hSFElecTrk = h; }
void SetFSMuonID (TH2D* h){ hSFMuonID  = h; }
void SetFSMuonIP1 (TH2D* h){ hSFMuonIP1  = h; }
void SetFSMuonIP2 (TH2D* h){ hSFMuonIP2  = h; }
void SetFSMuonISO(TH2D* h){ hSFMuonISO = h; }


Double_t _getFSIDoverRECO(Double_t pt, Double_t eta, Int_t pdgId, TString sys)
{

  eta = TMath::Abs(eta);
  Double_t sf   = 0.;
  Double_t sf_e = 0.;
  if (TMath::Abs(pdgId) == 13){
    if (pt > 200) pt = 199.9;
    sf   = hSFMuonID->GetBinContent( hSFMuonID->FindBin(pt,eta));
    sf_e = 0;
  }
  else{
    if (pt > 200) pt = 199.9;
    sf   = hSFElecID->GetBinContent( hSFElecID->FindBin(pt,eta));
    sf_e = 0.;
  }
  // error hardcoded to 2 later on
  // sf_e = (sys.Contains("Mu")) ? hSFMuonID->GetBinError  ( hSFMuonID->FindBin(pt,eta)) : 0.;

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
    sf    = hSFMuonIP1->GetBinContent( hSFMuonIP1->FindBin(pt,eta));
    sf   *= hSFMuonIP2->GetBinContent( hSFMuonIP2->FindBin(pt,eta));
    sf_e = 0;
  }
  else{
    if (pt > 200) pt = 199.9;
    sf   = hSFElecIP->GetBinContent( hSFElecIP->FindBin(pt,eta));
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
    sf   = hSFMuonISO->GetBinContent( hSFMuonISO->FindBin(pt,eta));
    sf_e = 0;
  }
  else{
    if (pt > 200) pt = 199.9;
    sf   = hSFElecISO->GetBinContent( hSFElecISO->FindBin(pt,eta));
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



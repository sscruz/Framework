//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                              //
// This program skims our trees making a basic selection that already contains at least 2 jet and good leptons  // 
//                                                                                                              //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Compile as: g++ Skim.C -o Skim `root-config --glibs --libs --cflags --ldflags` -lFoam -lTreePlayer -g//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include <math.h>
#include "TLorentzVector.h"
#include "TDirectory.h"
#include <string>
#include <cstdlib>
//---------------------------------------------------//
//              Some global variables                //
//---------------------------------------------------//
Int_t isDATA;

Int_t isData;

UInt_t run, lumi;

ULong64_t evt;

Int_t nLepTight_Edge, nPairLep_Edge, Lep1_pdgId_Edge, Lep2_pdgId_Edge, nJetSel_Edge, nPFHad10_Egde, nPFLep5_Edge, Flag_hbheFilterIso_Edge, Flag_hbheFilterNew25ns_Edge, Flag_HBHENoiseFilter_Edge,Flag_HBHENoiseIsoFilter_Edge, Flag_EcalDeadCellTriggerPrimitiveFilter_Edge, Flag_goodVertices_Edge, Flag_eeBadScFilter_Edge, Flag_globalTightHalo2016Filter_Edge, Flag_CSCTightHalo2016Filter_Edge;

Float_t Flag_eeBadScFilter, HLT_mu17mu8_dz, HLT_mu27tkmu8, HLT_mu17tkmu8_dz, HLT_el17el12_dz, HLT_ele33ele33, HLT_mu8el17, HLT_mu17el12, HLT_mu30ele30, HLT_DoubleMu, Flag_badMuonFilter;

Float_t HLT_pfht200, HLT_pfht250, HLT_pfht300, HLT_pfht400, HLT_pfht475, HLT_pfht600, HLT_pfht800, HLT_at51, HLT_at52, HLT_at53, HLT_at55;

Int_t nBJetLoose35_Edge, nBJetMedium35_Edge; 

Float_t mZ1_Edge, mZ2_Edge, mt2bb_Edge, mt2bb_jecUp_Edge ,  mt2bb_jecDown_Edge , mbb_Edge, mbb_jecUp_Edge,  mt2bb_jecDown_Edge ,met_pt, met_phi, lepsZPt_Edge, lepsJZB_Edge, lepsJZB_raw_Edge, lepsJZB_recoil_Edge, Lep1_phi_Edge, Lep2_phi_Edge, Lep1_eta_Edge, Lep2_eta_Edge, Lep1_pt_Edge, Lep2_pt_Edge, lepsDR_Edge, lepsMll_Edge;
    
Int_t nTrueInt, nVert;

Float_t puWeight, genWeight, PileupW_Edge;

TFile *f;

TFile *SkimFile, *SkimFileFriend;

TTree *tree;

TTree *SkimTree, *SkimFriendTree;

//---------------------------------------------------//
//              Some methods                         //
//---------------------------------------------------//
void setSourceBranches();
void setOutputBranches();


int main(int argc, char *argv[]) {

    if(argc != 6) {
      std::cout << "Usage: ./skim SourceTreeFile SourceTreeFriendFile DestinyTreeFile DestinyTreeFriendFile" << std::endl;
      return 1;
    }

    std::string fileName(argv[1]);
    std::string friendTree(argv[2]);
    std::string Skim(argv[3]);
    std::string SkimFriend(argv[4]);
    isDATA = std::atoi(argv[5]);

    //Opening the input file and getting the tree and the 3 counter histograms
    f = new TFile(fileName.c_str());
    tree = (TTree *)f->Get("tree");
    TH1F *Count = (TH1F *)f->Get("Count");
    TH1F *CountLHE, *SumGenWeights;
    if(!isDATA) { 
    	CountLHE = (TH1F *)f->Get("CountLHE");
    	SumGenWeights = (TH1F *)f->Get("SumGenWeights");
    }
    tree->AddFriend("sf/t", friendTree.c_str());
    
    setSourceBranches();

    //This will be the new Tree output file    
    SkimFile = new TFile(Skim.c_str(), "RECREATE");
    SkimTree = new TTree("tree", "tree");
   
    //This will be the new friend Tree output file    
    SkimFileFriend = new TFile(SkimFriend.c_str(), "RECREATE");
    TDirectory *dir = SkimFileFriend->mkdir("sf");
    SkimFriendTree = new TTree("t", "t");

    setOutputBranches();

    Int_t nentries = (Int_t)tree->GetEntries();
    for(Int_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        if(! (fabs(Lep1_eta_Edge)<2.4 && fabs(Lep2_eta_Edge)<2.4 && Lep1_pt_Edge > 20. && Lep2_pt_Edge > 20. )) continue;
        if(! (nJetSel_Edge > 1)) continue;
        if(! (nPairLep_Edge > 0)) continue;
        if(! ((lepsDR_Edge > 0.1 && (fabs(fabs(Lep1_eta_Edge) - 1.5) > 0.1 && fabs(fabs(Lep2_eta_Edge) - 1.5) > 0.1) && lepsMll_Edge > 20 ))) continue;
        SkimTree->Fill();
        SkimFriendTree->Fill();
    }

    SkimFile->cd();
    Count->Write();
    if(!isDATA) { 
        CountLHE->Write();
    	SumGenWeights->Write();
    }
    SkimTree->Write();
    SkimFile->Close();
    SkimFileFriend->cd();
    dir->cd(); 
    SkimFriendTree->Write();
    SkimFileFriend->Close();
    f->Close();

    return 0;

}




void setSourceBranches() {

    tree->SetBranchAddress("run",&run);
    tree->SetBranchAddress("lumi",&lumi);
    tree->SetBranchAddress("evt",&evt);
    tree->SetBranchAddress("isData",&isData);
    tree->SetBranchAddress("nPairLep_Edge",&nPairLep_Edge);
    tree->SetBranchAddress("nLepTight_Edge",&nLepTight_Edge);
    tree->SetBranchAddress("nPFHad10_Edge",&nPFHad10_Edge);
    tree->SetBranchAddress("nPFLep5_Edge",&nPFLep5_Edge);
    tree->SetBranchAddress("Lep1_phi_Edge",&Lep1_phi_Edge);
    tree->SetBranchAddress("Lep2_phi_Edge",&Lep2_phi_Edge);
    tree->SetBranchAddress("Lep1_pdgId_Edge",&Lep1_pdgId_Edge);
    tree->SetBranchAddress("Lep2_pdgId_Edge",&Lep2_pdgId_Edge);
    tree->SetBranchAddress("Flag_HBHENoiseFilter_Edge", &Flag_HBHENoiseFilter_Edge);                              
    tree->SetBranchAddress("Flag_HBHENoiseIsoFilter_Edge", &Flag_HBHENoiseIsoFilter_Edge);                          
    tree->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter_Edge", &Flag_EcalDeadCellTriggerPrimitiveFilter_Edge);        
    tree->SetBranchAddress("Flag_goodVertices_Edge",&Flag_goodVertices_Edge);                              
    tree->SetBranchAddress("Flag_eeBadScFilter_Edge",&Flag_eeBadScFilter_Edge);                            
    tree->SetBranchAddress("Flag_globalTightHalo2016Filter_Edge", &Flag_globalTightHalo2016Filter_Edge);               
    tree->SetBranchAddress("Flag_CSCTightHalo2016Filter_Edge", &Flag_CSCTightHalo2016Filter_Edge);                 
    tree->SetBranchAddress("Flag_badMuonFilter_Edge", &Flag_badMuonFilter_Edge);                                                    
    tree->SetBranchAddress("HLT_mu17mu8_dz",&HLT_mu17mu8_dz);
    tree->SetBranchAddress("HLT_mu27tkmu8",&HLT_mu27tkmu8);
    tree->SetBranchAddress("HLT_mu17tkmu8_dz",&HLT_mu17tkmu8_dz);
    tree->SetBranchAddress("HLT_el17el12_dz",&HLT_el17el12_dz);
    tree->SetBranchAddress("HLT_ele33ele33",&HLT_ele33ele33);
    tree->SetBranchAddress("HLT_mu8el17",&HLT_mu8el17);
    tree->SetBranchAddress("HLT_mu17el12",&HLT_mu17el12);
    tree->SetBranchAddress("HLT_mu30ele30",&HLT_mu30ele30);
    tree->SetBranchAddress("HLT_DoubleMu",&HLT_DoubleMu);
    tree->SetBranchAddress("HLT_pfht200",&HLT_pfht200);
    tree->SetBranchAddress("HLT_pfht250",&HLT_pfht250);
    tree->SetBranchAddress("HLT_pfht300",&HLT_pfht300);
    tree->SetBranchAddress("HLT_pfht400",&HLT_pfht400);
    tree->SetBranchAddress("HLT_pfht475",&HLT_pfht475);
    tree->SetBranchAddress("HLT_pfht600",&HLT_pfht600);
    tree->SetBranchAddress("HLT_pfht800",&HLT_pfht800);
    tree->SetBranchAddress("HLT_at51",&HLT_at51);
    tree->SetBranchAddress("HLT_at52",&HLT_at52);
    tree->SetBranchAddress("HLT_at53",&HLT_at53);
    tree->SetBranchAddress("HLT_at55",&HLT_at55);
    tree->SetBranchAddress("met_pt",&met_pt);
    tree->SetBranchAddress("met_phi",&met_phi);
    tree->SetBranchAddress("mZ1_Edge",&mZ1_Edge);
    tree->SetBranchAddress("mZ2_Edge",&mZ2_Edge);
    tree->SetBranchAddress("mt2bb_Edge",&mt2bb_Edge);
    tree->SetBranchAddress("mt2bb_jecUp_Edge",&mt2bb_jecUp_Edge);
    tree->SetBranchAddress("mt2bb_jecDown_Edge",&mt2bb_jecDown_Edge);
    tree->SetBranchAddress("mbb_Edge",&mbb_Edge);
    tree->SetBranchAddress("mbb_jecUp_Edge",&mbb_jecUp_Edge);
    tree->SetBranchAddress("mbb_jecDown_Edge",&mbb_jecDown_Edge);     
    tree->SetBranchAddress("nBJetLoose35_Edge",&nBJetLoose35_Edge);
    tree->SetBranchAddress("nBJetMedium35_Edge",&nBJetMedium35_Edge);
    tree->SetBranchAddress("lepsZPt_Edge",&lepsZPt_Edge);
    tree->SetBranchAddress("lepsJZB_Edge",&lepsJZB_Edge);
    tree->SetBranchAddress("lepsJZB_raw_Edge",&lepsJZB_raw_Edge);
    tree->SetBranchAddress("lepsJZB_recoil_Edge",&lepsJZB_recoil_Edge);
    tree->SetBranchAddress("Lep1_eta_Edge",&Lep1_eta_Edge);
    tree->SetBranchAddress("Lep2_eta_Edge",&Lep2_eta_Edge);
    tree->SetBranchAddress("Lep1_pt_Edge",&Lep1_pt_Edge);
    tree->SetBranchAddress("Lep2_pt_Edge",&Lep2_pt_Edge);
    tree->SetBranchAddress("lepsDR_Edge",&lepsDR_Edge);
    tree->SetBranchAddress("nJetSel_Edge",&nJetSel_Edge);
    tree->SetBranchAddress("lepsMll_Edge",&lepsMll_Edge);
    tree->SetBranchAddress("hbheFilterIso", &hbheFilterIso);
    tree->SetBranchAddress("hbheFilterNew25ns", &hbheFilterNew25ns);
    tree->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter);
    tree->SetBranchAddress("PileupW_Edge", &PileupW_Edge);
    tree->SetBranchAddress("nVert", &nVert);
    if(!isDATA) {
    	tree->SetBranchAddress("puWeight", &puWeight);
	tree->SetBranchAddress("nTrueInt", &nTrueInt);
        tree->SetBranchAddress("genWeight", &genWeight);
    }
    
}


void setOutputBranches() {

    SkimTree->Branch("run",&run);
    SkimTree->Branch("lumi",&lumi);
    SkimTree->Branch("evt",&evt);
    SkimTree->Branch("isData",&isData);
    SkimFriendTree->Branch("nPairLep_Edge",&nPairLep_Edge);
    SkimFriendTree->Branch("nLepTight_Edge",&nLepTight_Edge);
    SkimFriendTree->Branch("nPFHad10_Edge",&nPFHad10_Edge);
    SkimFriendTree->Branch("nPFLep5_Edge",&nPFLep5_Edge);
    SkimFriendTree->Branch("Lep1_phi_Edge",&Lep1_phi_Edge);
    SkimFriendTree->Branch("Lep2_phi_Edge",&Lep2_phi_Edge);
    SkimFriendTree->Branch("Lep1_pdgId_Edge",&Lep1_pdgId_Edge);
    SkimFriendTree->Branch("Lep2_pdgId_Edge",&Lep2_pdgId_Edge);
    SkimFriendTree->Branch("Flag_HBHENoiseFilter_Edge", &Flag_HBHENoiseFilter_Edge);                              
    SkimFriendTree->Branch("Flag_HBHENoiseIsoFilter_Edge", &Flag_HBHENoiseIsoFilter_Edge);                          
    SkimFriendTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter_Edge", &Flag_EcalDeadCellTriggerPrimitiveFilter_Edge);        
    SkimFriendTree->Branch("Flag_goodVertices_Edge",&Flag_goodVertices_Edge);                              
    SkimFriendTree->Branch("Flag_eeBadScFilter_Edge",&Flag_eeBadScFilter_Edge);                            
    SkimFriendTree->Branch("Flag_globalTightHalo2016Filter_Edge", &Flag_globalTightHalo2016Filter_Edge);               
    SkimFriendTree->Branch("Flag_CSCTightHalo2016Filter_Edge", &Flag_CSCTightHalo2016Filter_Edge);                 
    SkimFriendTree->Branch("Flag_badMuonFilter_Edge", &Flag_badMuonFilter_Edge);                                                   
    SkimTree->Branch("HLT_mu17mu8_dz",&HLT_mu17mu8_dz);
    SkimTree->Branch("HLT_mu27tkmu8",&HLT_mu27tkmu8);
    SkimTree->Branch("HLT_mu17tkmu8_dz",&HLT_mu17tkmu8_dz);
    SkimTree->Branch("HLT_el17el12_dz",&HLT_el17el12_dz);
    SkimTree->Branch("HLT_ele33ele33",&HLT_ele33ele33);
    SkimTree->Branch("HLT_mu8el17",&HLT_mu8el17);
    SkimTree->Branch("HLT_mu17el12",&HLT_mu17el12);
    SkimTree->Branch("HLT_mu30ele30",&HLT_mu30ele30);
    SkimTree->Branch("HLT_DoubleMu",&HLT_DoubleMu);
    SkimTree->Branch("HLT_pfht200",&HLT_pfht200);
    SkimTree->Branch("HLT_pfht250",&HLT_pfht250);
    SkimTree->Branch("HLT_pfht300",&HLT_pfht300);
    SkimTree->Branch("HLT_pfht400",&HLT_pfht400);
    SkimTree->Branch("HLT_pfht475",&HLT_pfht475);
    SkimTree->Branch("HLT_pfht600",&HLT_pfht600);
    SkimTree->Branch("HLT_pfht800",&HLT_pfht800);
    SkimTree->Branch("HLT_at51",&HLT_at51);
    SkimTree->Branch("HLT_at52",&HLT_at52);
    SkimTree->Branch("HLT_at53",&HLT_at53);
    SkimTree->Branch("HLT_at55",&HLT_at55);
    SkimTree->Branch("met_pt",&met_pt);
    SkimTree->Branch("met_phi",&met_phi);
    SkimTree->Branch("met_phi",&met_phi);
    SkimTree->Branch("mZ1_Edge",&mZ1_Edge);
    SkimTree->Branch("mZ2_Edge",&mZ2_Edge);
    SkimTree->Branch("mt2bb_Edge",&mt2bb_Edge);
    SkimTree->Branch("mt2bb_jecUp_Edge",&mt2bb_jecUp_Edge);
    SkimTree->Branch("mt2bb_jecDown_Edge",&mt2bb_jecDown_Edge);
    SkimTree->Branch("mbb_Edge",&mbb_Edge);
    SkimTree->Branch("mbb_jecUp_Edge",&mbb_jecUp_Edge);
    SkimTree->Branch("mbb_jecDown_Edge",&mbb_jecDown_Edge);     
    SkimFriendTree->Branch("nBJetLoose35_Edge",&nBJetLoose35_Edge);
    SkimFriendTree->Branch("nBJetMedium35_Edge",&nBJetMedium35_Edge);
    SkimFriendTree->Branch("lepsZPt_Edge",&lepsZPt_Edge);
    SkimFriendTree->Branch("lepsJZB_Edge",&lepsJZB_Edge);
    SkimFriendTree->Branch("lepsJZB_raw_Edge",&lepsJZB_raw_Edge);
    SkimFriendTree->Branch("lepsJZB_recoil_Edge",&lepsJZB_recoil_Edge);
    SkimFriendTree->Branch("Lep1_eta_Edge",&Lep1_eta_Edge);
    SkimFriendTree->Branch("Lep2_eta_Edge",&Lep2_eta_Edge);
    SkimFriendTree->Branch("Lep1_pt_Edge",&Lep1_pt_Edge);
    SkimFriendTree->Branch("Lep2_pt_Edge",&Lep2_pt_Edge);
    SkimFriendTree->Branch("lepsDR_Edge",&lepsDR_Edge);
    SkimFriendTree->Branch("nJetSel_Edge",&nJetSel_Edge);
    SkimFriendTree->Branch("lepsMll_Edge",&lepsMll_Edge);
    SkimFriendTree->Branch("PileupW_Edge", &PileupW_Edge);
    SkimTree->Branch("hbheFilterIso", &hbheFilterIso);
    SkimTree->Branch("hbheFilterNew25ns", &hbheFilterNew25ns);
    SkimTree->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter);
    SkimTree->Branch("nVert", &nVert);
    if(!isDATA) {
    	SkimTree->Branch("puWeight", &puWeight);
    	SkimTree->Branch("nTrueInt", &nTrueInt);
    	SkimTree->Branch("genWeight", &genWeight);
    }
}














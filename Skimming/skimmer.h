//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec  7 19:00:27 2016 by ROOT version 6.06/08
// from TTree t/t
// found on file: /pool/ciencias/HeppyTrees/EdgeZ/Moriond17/finalFriends/evVarFriend_DoubleEG_Run2016F_23Sep2016_v1_runs_271036_284044_part1.root
//////////////////////////////////////////////////////////

#ifndef skimmer_h
#define skimmer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>

// Header file for the classes stored in the TTree if any.

class skimmer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TTree          *outputtree;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TH1D           *counts;
   TH1D           *genWeights;
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        evt_Edge;
   Int_t           run_Edge;
   Int_t           lumi_Edge;
   Int_t           nVert_Edge;
   Int_t           nPFHad10_Edge;
   Int_t           nPFLep5_Edge;
   Int_t           Flag_HBHENoiseFilter_Edge;
   Int_t           Flag_HBHENoiseIsoFilter_Edge;
   Int_t           Flag_EcalDeadCellTriggerPrimitiveFilter_Edge;
   Int_t           Flag_goodVertices_Edge;
   Int_t           Flag_eeBadScFilter_Edge;
   Int_t           Flag_globalTightHalo2016Filter_Edge;
   Int_t           Flag_CSCTightHalo2016Filter_Edge;                      
   Float_t         Flag_badMuonFilter_Edge;                      
   Float_t           Flag_badChargedHadronFilter_Edge;                      
   Int_t           Flag_badCloneMuonMoriond2017_Edge;                      
   Int_t           Flag_badMuonMoriond2017_Edge;                      
   Float_t         badCloneMuonMoriond2017_maxPt_Edge;                      
   Float_t         badNotCloneMuonMoriond2017_maxPt_Edge;                      
   Int_t           nLepTight_Edge;
   Int_t           nLepLoose_Edge;
   Int_t           nJetSel_Edge;
   Int_t           nJetSel_jecUp_Edge;
   Int_t           nJetSel_jecDn_Edge;
   Float_t         bestMjj_Edge;
   Float_t         minMjj_Edge;
   Float_t         maxMjj_Edge;
   Float_t         hardMjj_Edge;
   Float_t         dphiMjj_Edge;
   Float_t         drMjj_Edge;
   Float_t         hardJJDphi_Edge;
   Float_t         hardJJDR_Edge;
   Float_t         j1MetDPhi_Edge;
   Float_t         j2MetDPhi_Edge;
   Int_t           nPairLep_Edge;
   Int_t           iLT_Edge[3];   //[nLepTight_Edge]
   Int_t           iJ_Edge[13];   //[nJetSel_Edge]
   Int_t           nLepGood20_Edge;
   Int_t           nLepGood20T_Edge;
   Int_t           nJet35_Edge;
   Int_t           nJet35_jecUp_Edge;
   Int_t           nJet35_jecDn_Edge;
   Float_t         htJet35j_Edge;
   Float_t         htJet35j_jecUp_Edge;
   Float_t         htJet35j_jecDn_Edge;
   Int_t           nBJetMedium25_Edge;
   Int_t           nBJetMedium25_jecUp_Edge;
   Int_t           nBJetMedium25_jecDn_Edge;
   Int_t           nBJetLoose35_Edge;
   Int_t           nBJetLoose35_jecUp_Edge;
   Int_t           nBJetLoose35_jecDn_Edge;
   Int_t           nBJetMedium35_Edge;
   Int_t           nBJetMedium35_jecUp_Edge;
   Int_t           nBJetMedium35_jecDn_Edge;
   Int_t           iL1T_Edge;
   Int_t           iL2T_Edge;
   Float_t         lepsMll_Edge;
   Float_t         lepsJZB_Edge;
   Float_t         lepsJZB_raw_Edge;
   Float_t         lepsJZB_recoil_Edge;
   Float_t         lepsDR_Edge;
   Float_t         lepsMETRec_Edge;
   Float_t         lepsZPt_Edge;
   Float_t         metl1DPhi_Edge;
   Float_t         metl2DPhi_Edge;
   Float_t         met_Edge;
   Float_t         met_phi_Edge;
   Float_t         met_jecUp_Edge;
   Float_t         met_jecDn_Edge;
   Float_t         met_raw_Edge;
   Float_t         mZ1_Edge;
   Float_t         mZ2_Edge;
   Float_t         mt2bb_Edge;
   Float_t         mt2bb_jecUp_Edge;
   Float_t         mt2bb_jecDn_Edge;
   Float_t         mbb_Edge;
   Float_t         mbb_jecUp_Edge;
   Float_t         mbb_jecDn_Edge;     
   Float_t         genMet_Edge;
   Float_t         genMet_phi_Edge;
   Float_t         lepsDPhi_Edge;
   Float_t         Lep1_pt_Edge;
   Float_t         Lep1_eta_Edge;
   Float_t         Lep1_phi_Edge;
   Float_t         Lep1_miniRelIso_Edge;
   Float_t         Lep1_relIso03_Edge;
   Float_t         Lep1_relIso04_Edge;
   Float_t         Lep1_dxy_Edge;
   Float_t         Lep1_dz_Edge;
   Float_t         Lep1_sip3d_Edge;
   Int_t           Lep1_pdgId_Edge;
   Float_t         Lep1_tightCharge_Edge;
   Float_t         Lep1_mvaIdSpring15_Edge;
   Float_t         Lep1_mcMatchId_Edge;
   Float_t         Lep1_mcMatchTau_Edge;
   Float_t         Lep1_minTauDR_Edge;
   Float_t         Lep2_pt_Edge;
   Float_t         Lep2_eta_Edge;
   Float_t         Lep2_phi_Edge;
   Float_t         Lep2_miniRelIso_Edge;
   Float_t         Lep2_relIso03_Edge;
   Float_t         Lep2_relIso04_Edge;
   Float_t         Lep2_dxy_Edge;
   Float_t         Lep2_dz_Edge;
   Float_t         Lep2_sip3d_Edge;
   Int_t           Lep2_pdgId_Edge;
   Float_t         Lep2_tightCharge_Edge;
   Float_t         Lep2_mvaIdSpring15_Edge;
   Float_t         Lep2_mcMatchId_Edge;
   Float_t         Lep2_mcMatchTau_Edge;
   Float_t         Lep2_minTauDR_Edge;
   Float_t         PileupW_Edge;
   Float_t         min_mlb1_Edge;
   Float_t         min_mlb2_Edge;
   Float_t         min_mlb1Up_Edge;
   Float_t         min_mlb2Up_Edge;
   Float_t         min_mlb1Dn_Edge;
   Float_t         min_mlb2Dn_Edge;
   Float_t         sum_mlb_Edge;
   Float_t         sum_mlbUp_Edge;
   Float_t         sum_mlbDn_Edge;
   Float_t         st_Edge;
   Int_t           srID_Edge;
   Float_t         mt2_Edge;
   Float_t         mt2_jecUp_Edge;
   Float_t         mt2_jecDn_Edge;
   Float_t         weight_trigger_Edge;
   Float_t         weight_btagsf_Edge;
   Float_t         weight_btagsf_heavy_UP_Edge;
   Float_t         weight_btagsf_heavy_DN_Edge;
   Float_t         weight_btagsf_light_UP_Edge;
   Float_t         weight_btagsf_light_DN_Edge;
   Float_t         d3D_Edge;
   Float_t         parPt_Edge;
   Float_t         ortPt_Edge;
   Float_t         dTheta_Edge;
   Float_t         genWeight_Edge;
   Float_t         HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v_Edge;
   Float_t         HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_Edge;
   Float_t         HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v_Edge;
   Float_t         HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v_Edge;
   Float_t         HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge;
   Float_t         HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge;
   Float_t         HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v_Edge;
   Float_t         HLT_BIT_HLT_Mu27_TkMu8_v_Edge;
   Float_t         HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge;
   Float_t         HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge;
   Float_t         HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v_Edge;
   Float_t         HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v_Edge;
   Float_t         HLT_BIT_HLT_Mu30_TkMu11_v_Edge;
   Float_t         HLT_BIT_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v_Edge;
   Float_t         HLT_BIT_HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v_Edge;
   Float_t         HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v_Edge;
   Float_t         HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v_Edge;
   Float_t         HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge;
   Float_t         HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge;
   Float_t         HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge;    
   Float_t         HLT_BIT_HLT_PFHT200_v_Edge;
   Float_t         HLT_BIT_HLT_PFHT250_v_Edge;
   Float_t         HLT_BIT_HLT_PFHT300_v_Edge;
   Float_t         HLT_BIT_HLT_PFHT350_v_Edge;
   Float_t         HLT_BIT_HLT_PFHT400_v_Edge;
   Float_t         HLT_BIT_HLT_PFHT475_v_Edge;
   Float_t         HLT_BIT_HLT_PFHT600_v_Edge;
   Float_t         HLT_BIT_HLT_PFHT650_v_Edge;
   Float_t         HLT_BIT_HLT_PFHT800_v_Edge;
   Float_t         HLT_BIT_HLT_PFHT300_PFMET110_v_Edge;
   Float_t         GenSusyMScan1_Edge;
   Float_t         GenSusyMScan2_Edge;
   Float_t         GenSusyMScan3_Edge;
   Float_t         GenSusyMScan4_Edge;
   Float_t         GenSusyMGluino_Edge;
   Float_t         GenSusyMGravitino_Edge;
   Float_t         GenSusyMStop_Edge;
   Float_t         GenSusyMSbottom_Edge;
   Float_t         GenSusyMStop2_Edge;
   Float_t         GenSusyMSbottom2_Edge;
   Float_t         GenSusyMSquark_Edge;
   Float_t         GenSusyMNeutralino_Edge;
   Float_t         GenSusyMNeutralino2_Edge;
   Float_t         GenSusyMNeutralino3_Edge;
   Float_t         GenSusyMNeutralino4_Edge;
   Float_t         GenSusyMChargino_Edge;
   Float_t         GenSusyMChargino2_Edge;
   Float_t         JetSel_Edge_pt[13];   //[nJetSel_Edge]
   Float_t         JetSel_Edge_eta[13];   //[nJetSel_Edge]
   Float_t         JetSel_Edge_phi[13];   //[nJetSel_Edge]
   Float_t         JetSel_Edge_mass[13];   //[nJetSel_Edge]
   Float_t         JetSel_Edge_btagCSV[13];   //[nJetSel_Edge]
   Float_t         JetSel_Edge_rawPt[13];   //[nJetSel_Edge]
   Float_t         JetSel_Edge_mcPt[13];   //[nJetSel_Edge]
   Int_t           JetSel_Edge_mcFlavour[13];   //[nJetSel_Edge]
   Int_t           JetSel_Edge_mcMatchId[13];   //[nJetSel_Edge]

   // List of branches
   TBranch        *b_evt_Edge;   //!
   TBranch        *b_run_Edge;   //!
   TBranch        *b_lumi_Edge;   //!
   TBranch        *b_nVert_Edge;   //!
   TBranch        *b_nLepTight_Edge;   //!
   TBranch        *b_nPFHad10_Edge;   //!
   TBranch        *b_nPFLep5_Edge;   //!
   TBranch        *b_Flag_HBHENoiseFilter_Edge;
   TBranch        *b_Flag_HBHENoiseIsoFilter_Edge;
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter_Edge;
   TBranch        *b_Flag_goodVertices_Edge;
   TBranch        *b_Flag_eeBadScFilter_Edge;
   TBranch        *b_Flag_globalTightHalo2016Filter_Edge;
   TBranch        *b_Flag_CSCTightHalo2016Filter_Edge;             
   TBranch        *b_Flag_badMuonFilter_Edge;                                    
   TBranch        *b_Flag_badChargedHadronFilter_Edge;                                    
   TBranch        *b_Flag_badCloneMuonMoriond2017_Edge;                                    
   TBranch        *b_Flag_badMuonMoriond2017_Edge;                                    
   TBranch        *b_badCloneMuonMoriond2017_maxPt_Edge;                                    
   TBranch        *b_badNotCloneMuonMoriond2017_maxPt_Edge;                                    
   TBranch        *b_nLepLoose_Edge;   //!
   TBranch        *b_nJetSel_Edge;   //!
   TBranch        *b_nJetSel_jecUp_Edge;   //!
   TBranch        *b_nJetSel_jecDn_Edge;   //!
   TBranch        *b_bestMjj_Edge;   //!
   TBranch        *b_minMjj_Edge;   //!
   TBranch        *b_maxMjj_Edge;   //!
   TBranch        *b_hardMjj_Edge;   //!
   TBranch        *b_dphiMjj_Edge;   //!
   TBranch        *b_drMjj_Edge;   //!
   TBranch        *b_hardJJDphi_Edge;   //!
   TBranch        *b_hardJJDR_Edge;   //!
   TBranch        *b_j1MetDPhi_Edge;   //!
   TBranch        *b_j2MetDPhi_Edge;   //!
   TBranch        *b_nPairLep_Edge;   //!
   TBranch        *b_iLT_Edge;   //!
   TBranch        *b_iJ_Edge;   //!
   TBranch        *b_nLepGood20_Edge;   //!
   TBranch        *b_nLepGood20T_Edge;   //!
   TBranch        *b_nJet35_Edge;   //!
   TBranch        *b_nJet35_jecUp_Edge;   //!
   TBranch        *b_nJet35_jecDn_Edge;   //!
   TBranch        *b_htJet35j_Edge;   //!
   TBranch        *b_htJet35j_jecUp_Edge;   //!
   TBranch        *b_htJet35j_jecDn_Edge;   //!
   TBranch        *b_nBJetMedium25_Edge;   //!
   TBranch        *b_nBJetMedium25_jecUp_Edge;   //!
   TBranch        *b_nBJetMedium25_jecDn_Edge;   //!
   TBranch        *b_nBJetLoose35_Edge;   //!
   TBranch        *b_nBJetLoose35_jecUp_Edge;   //!
   TBranch        *b_nBJetLoose35_jecDn_Edge;   //!
   TBranch        *b_nBJetMedium35_Edge;   //!
   TBranch        *b_nBJetMedium35_jecUp_Edge;   //!
   TBranch        *b_nBJetMedium35_jecDn_Edge;   //!
   TBranch        *b_iL1T_Edge;   //!
   TBranch        *b_iL2T_Edge;   //!
   TBranch        *b_lepsMll_Edge;   //!
   TBranch        *b_lepsJZB_Edge;   //!
   TBranch        *b_lepsJZB_raw_Edge;   //!
   TBranch        *b_lepsJZB_recoil_Edge;   //!
   TBranch        *b_lepsDR_Edge;   //!
   TBranch        *b_lepsMETRec_Edge;   //!
   TBranch        *b_lepsZPt_Edge;   //!
   TBranch        *b_metl1DPhi_Edge;   //!
   TBranch        *b_metl2DPhi_Edge;   //!
   TBranch        *b_met_Edge;   //!
   TBranch        *b_met_phi_Edge;   //!
   TBranch        *b_met_jecUp_Edge;   //!
   TBranch        *b_met_jecDn_Edge;   //!
   TBranch        *b_met_raw_Edge;   //!
   TBranch        *b_mZ1_Edge;
   TBranch        *b_mZ2_Edge;
   TBranch        *b_mt2bb_Edge;
   TBranch        *b_mt2bb_jecUp_Edge;
   TBranch        *b_mt2bb_jecDn_Edge;
   TBranch        *b_mbb_Edge;
   TBranch        *b_mbb_jecUp_Edge;
   TBranch        *b_mbb_jecDn_Edge;     
   TBranch        *b_genMet_Edge;   //!
   TBranch        *b_genMet_phi_Edge;   //!
   TBranch        *b_lepsDPhi_Edge;   //!
   TBranch        *b_Lep1_pt_Edge;   //!
   TBranch        *b_Lep1_eta_Edge;   //!
   TBranch        *b_Lep1_phi_Edge;   //!
   TBranch        *b_Lep1_miniRelIso_Edge;   //!
   TBranch        *b_Lep1_relIso03_Edge;   //!
   TBranch        *b_Lep1_relIso04_Edge;   //!
   TBranch        *b_Lep1_dxy_Edge;   //!
   TBranch        *b_Lep1_dz_Edge;   //!
   TBranch        *b_Lep1_sip3d_Edge;   //!
   TBranch        *b_Lep1_pdgId_Edge;   //!
   TBranch        *b_Lep1_tightCharge_Edge;   //!
   TBranch        *b_Lep1_mvaIdSpring15_Edge;   //!
   TBranch        *b_Lep1_mcMatchId_Edge;   //!
   TBranch        *b_Lep1_mcMatchTau_Edge;   //!
   TBranch        *b_Lep1_minTauDR_Edge;   //!
   TBranch        *b_Lep2_pt_Edge;   //!
   TBranch        *b_Lep2_eta_Edge;   //!
   TBranch        *b_Lep2_phi_Edge;   //!
   TBranch        *b_Lep2_miniRelIso_Edge;   //!
   TBranch        *b_Lep2_relIso03_Edge;   //!
   TBranch        *b_Lep2_relIso04_Edge;   //!
   TBranch        *b_Lep2_dxy_Edge;   //!
   TBranch        *b_Lep2_dz_Edge;   //!
   TBranch        *b_Lep2_sip3d_Edge;   //!
   TBranch        *b_Lep2_pdgId_Edge;   //!
   TBranch        *b_Lep2_tightCharge_Edge;   //!
   TBranch        *b_Lep2_mvaIdSpring15_Edge;   //!
   TBranch        *b_Lep2_mcMatchId_Edge;   //!
   TBranch        *b_Lep2_mcMatchTau_Edge;   //!
   TBranch        *b_Lep2_minTauDR_Edge;   //!
   TBranch        *b_PileupW_Edge;   //!
   TBranch        *b_min_mlb1_Edge;   //!
   TBranch        *b_min_mlb2_Edge;   //!
   TBranch        *b_min_mlb1Up_Edge;   //!
   TBranch        *b_min_mlb2Up_Edge;   //!
   TBranch        *b_min_mlb1Dn_Edge;   //!
   TBranch        *b_min_mlb2Dn_Edge;   //!
   TBranch        *b_sum_mlb_Edge;   //!
   TBranch        *b_sum_mlbUp_Edge;   //!
   TBranch        *b_sum_mlbDn_Edge;   //!
   TBranch        *b_st_Edge;   //!
   TBranch        *b_srID_Edge;   //!
   TBranch        *b_mt2_Edge;   //!
   TBranch        *b_mt2_jecUp_Edge;   //!
   TBranch        *b_mt2_jecDn_Edge;   //!
   TBranch        *b_weight_trigger_Edge;   //!
   TBranch        *b_weight_btagsf_Edge;   //!
   TBranch        *b_weight_btagsf_heavy_UP_Edge;   //!
   TBranch        *b_weight_btagsf_heavy_DN_Edge;   //!
   TBranch        *b_weight_btagsf_light_UP_Edge;   //!
   TBranch        *b_weight_btagsf_light_DN_Edge;   //!
   TBranch        *b_d3D_Edge;   //!
   TBranch        *b_parPt_Edge;   //!
   TBranch        *b_ortPt_Edge;   //!
   TBranch        *b_dTheta_Edge;   //!
   TBranch        *b_genWeight_Edge;   //!
   TBranch        *b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Mu27_TkMu8_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Mu30_TkMu11_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v_Edge;
   TBranch        *b_HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v_Edge;
   TBranch        *b_HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge;
   TBranch        *b_HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge;    
   TBranch        *b_HLT_BIT_HLT_PFHT200_v_Edge;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT250_v_Edge;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT300_v_Edge;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT350_v_Edge;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT400_v_Edge;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT475_v_Edge;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT600_v_Edge;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT650_v_Edge;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT800_v_Edge;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT300_PFMET110_v_Edge;   //!
   TBranch        *b_GenSusyMScan1_Edge;   //!
   TBranch        *b_GenSusyMScan2_Edge;   //!
   TBranch        *b_GenSusyMScan3_Edge;   //!
   TBranch        *b_GenSusyMScan4_Edge;   //!
   TBranch        *b_GenSusyMGluino_Edge;   //!
   TBranch        *b_GenSusyMGravitino_Edge;   //!
   TBranch        *b_GenSusyMStop_Edge;   //!
   TBranch        *b_GenSusyMSbottom_Edge;   //!
   TBranch        *b_GenSusyMStop2_Edge;   //!
   TBranch        *b_GenSusyMSbottom2_Edge;   //!
   TBranch        *b_GenSusyMSquark_Edge;   //!
   TBranch        *b_GenSusyMNeutralino_Edge;   //!
   TBranch        *b_GenSusyMNeutralino2_Edge;   //!
   TBranch        *b_GenSusyMNeutralino3_Edge;   //!
   TBranch        *b_GenSusyMNeutralino4_Edge;   //!
   TBranch        *b_GenSusyMChargino_Edge;   //!
   TBranch        *b_GenSusyMChargino2_Edge;   //!
   TBranch        *b_JetSel_Edge_pt;   //!
   TBranch        *b_JetSel_Edge_eta;   //!
   TBranch        *b_JetSel_Edge_phi;   //!
   TBranch        *b_JetSel_Edge_mass;   //!
   TBranch        *b_JetSel_Edge_btagCSV;   //!
   TBranch        *b_JetSel_Edge_rawPt;   //!
   TBranch        *b_JetSel_Edge_mcPt;   //!
   TBranch        *b_JetSel_Edge_mcFlavour;   //!
   TBranch        *b_JetSel_Edge_mcMatchId;   //!

   skimmer(TString, TString);
   virtual ~skimmer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     SetOutVariables();
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   TString sample;
   TString path;
};

#endif

#ifdef skimmer_cxx
skimmer::skimmer(TString sampleName, TString pathString) : fChain(0) 
{
  sample = sampleName;
  path   = pathString;
  TFile* f = TFile::Open(path + "evVarFriend_" + sample + ".root","READ");
  TTree* tree = (TTree*) f->Get("sf/t");
  counts = (TH1D *) f->Get("Count");
  genWeights = (TH1D *) f->Get("SumGenWeights");

   /* if (tree == 0) { */
   /*    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/pool/ciencias/HeppyTrees/EdgeZ/Moriond17/finalFriends/evVarFriend_DoubleEG_Run2016F_23Sep2016_v1_runs_271036_284044_part1.root"); */
   /*    if (!f || !f->IsOpen()) { */
   /*       f = new TFile("/pool/ciencias/HeppyTrees/EdgeZ/Moriond17/finalFriends/evVarFriend_DoubleEG_Run2016F_23Sep2016_v1_runs_271036_284044_part1.root"); */
   /*    } */
   /*    TDirectory * dir = (TDirectory*)f->Get("/pool/ciencias/HeppyTrees/EdgeZ/Moriond17/finalFriends/evVarFriend_DoubleEG_Run2016F_23Sep2016_v1_runs_271036_284044_part1.root:/sf"); */
   /*    dir->GetObject("t",tree); */

   /* } */
   Init(tree);
   Loop();
}

skimmer::~skimmer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t skimmer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t skimmer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void skimmer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evt_Edge", &evt_Edge, &b_evt_Edge);
   fChain->SetBranchAddress("run_Edge", &run_Edge, &b_run_Edge);
   fChain->SetBranchAddress("lumi_Edge", &lumi_Edge, &b_lumi_Edge);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter_Edge", &Flag_HBHENoiseFilter_Edge, &b_Flag_HBHENoiseFilter_Edge); 
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter_Edge", &Flag_HBHENoiseIsoFilter_Edge, &b_Flag_HBHENoiseIsoFilter_Edge);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter_Edge", &Flag_EcalDeadCellTriggerPrimitiveFilter_Edge, &b_Flag_EcalDeadCellTriggerPrimitiveFilter_Edge);
   fChain->SetBranchAddress("Flag_goodVertices_Edge", &Flag_goodVertices_Edge, &b_Flag_goodVertices_Edge);
   fChain->SetBranchAddress("Flag_eeBadScFilter_Edge", &Flag_eeBadScFilter_Edge, &b_Flag_eeBadScFilter_Edge);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter_Edge", &Flag_globalTightHalo2016Filter_Edge, &b_Flag_globalTightHalo2016Filter_Edge);
   fChain->SetBranchAddress("Flag_CSCTightHalo2016Filter_Edge", &Flag_CSCTightHalo2016Filter_Edge, &b_Flag_CSCTightHalo2016Filter_Edge);
   fChain->SetBranchAddress("Flag_badMuonFilter_Edge", &Flag_badMuonFilter_Edge, &b_Flag_badMuonFilter_Edge);                                                                  
   fChain->SetBranchAddress("Flag_badChargedHadronFilter_Edge", &Flag_badChargedHadronFilter_Edge, &b_Flag_badChargedHadronFilter_Edge);                                                       
   fChain->SetBranchAddress("Flag_badCloneMuonMoriond2017_Edge", &Flag_badCloneMuonMoriond2017_Edge, &b_Flag_badCloneMuonMoriond2017_Edge);                                                       
   fChain->SetBranchAddress("Flag_badMuonMoriond2017_Edge", &Flag_badMuonMoriond2017_Edge, &b_Flag_badMuonMoriond2017_Edge);                                                       
   fChain->SetBranchAddress("badCloneMuonMoriond2017_maxPt_Edge", &badCloneMuonMoriond2017_maxPt_Edge, &b_badCloneMuonMoriond2017_maxPt_Edge);                                                       
   fChain->SetBranchAddress("badNotCloneMuonMoriond2017_maxPt_Edge", &badNotCloneMuonMoriond2017_maxPt_Edge, &b_badNotCloneMuonMoriond2017_maxPt_Edge);                                            
   fChain->SetBranchAddress("nLepTight_Edge", &nLepTight_Edge, &b_nLepTight_Edge);
   fChain->SetBranchAddress("nLepLoose_Edge", &nLepLoose_Edge, &b_nLepLoose_Edge);
   fChain->SetBranchAddress("nPFHad10_Edge", &nPFHad10_Edge, &b_nPFHad10_Edge);
   fChain->SetBranchAddress("nPFLep5_Edge", &nPFLep5_Edge, &b_nPFLep5_Edge);
   fChain->SetBranchAddress("nJetSel_Edge", &nJetSel_Edge, &b_nJetSel_Edge);
   fChain->SetBranchAddress("nJetSel_jecUp_Edge", &nJetSel_jecUp_Edge, &b_nJetSel_jecUp_Edge);
   fChain->SetBranchAddress("nJetSel_jecDn_Edge", &nJetSel_jecDn_Edge, &b_nJetSel_jecDn_Edge);
   fChain->SetBranchAddress("bestMjj_Edge", &bestMjj_Edge, &b_bestMjj_Edge);
   fChain->SetBranchAddress("minMjj_Edge", &minMjj_Edge, &b_minMjj_Edge);
   fChain->SetBranchAddress("maxMjj_Edge", &maxMjj_Edge, &b_maxMjj_Edge);
   fChain->SetBranchAddress("hardMjj_Edge", &hardMjj_Edge, &b_hardMjj_Edge);
   fChain->SetBranchAddress("dphiMjj_Edge", &dphiMjj_Edge, &b_dphiMjj_Edge);
   fChain->SetBranchAddress("drMjj_Edge", &drMjj_Edge, &b_drMjj_Edge);
   fChain->SetBranchAddress("hardJJDphi_Edge", &hardJJDphi_Edge, &b_hardJJDphi_Edge);
   fChain->SetBranchAddress("hardJJDR_Edge", &hardJJDR_Edge, &b_hardJJDR_Edge);
   fChain->SetBranchAddress("j1MetDPhi_Edge", &j1MetDPhi_Edge, &b_j1MetDPhi_Edge);
   fChain->SetBranchAddress("j2MetDPhi_Edge", &j2MetDPhi_Edge, &b_j2MetDPhi_Edge);
   fChain->SetBranchAddress("nPairLep_Edge", &nPairLep_Edge, &b_nPairLep_Edge);
   fChain->SetBranchAddress("iLT_Edge", iLT_Edge, &b_iLT_Edge);
   fChain->SetBranchAddress("iJ_Edge", iJ_Edge, &b_iJ_Edge);
   fChain->SetBranchAddress("nLepGood20_Edge", &nLepGood20_Edge, &b_nLepGood20_Edge);
   fChain->SetBranchAddress("nLepGood20T_Edge", &nLepGood20T_Edge, &b_nLepGood20T_Edge);
   fChain->SetBranchAddress("nJet35_Edge", &nJet35_Edge, &b_nJet35_Edge);
   fChain->SetBranchAddress("nJet35_jecUp_Edge", &nJet35_jecUp_Edge, &b_nJet35_jecUp_Edge);
   fChain->SetBranchAddress("nJet35_jecDn_Edge", &nJet35_jecDn_Edge, &b_nJet35_jecDn_Edge);
   fChain->SetBranchAddress("htJet35j_Edge", &htJet35j_Edge, &b_htJet35j_Edge);
   fChain->SetBranchAddress("htJet35j_jecUp_Edge", &htJet35j_jecUp_Edge, &b_htJet35j_jecUp_Edge);
   fChain->SetBranchAddress("htJet35j_jecDn_Edge", &htJet35j_jecDn_Edge, &b_htJet35j_jecDn_Edge);
   fChain->SetBranchAddress("nBJetMedium25_Edge", &nBJetMedium25_Edge, &b_nBJetMedium25_Edge);
   fChain->SetBranchAddress("nBJetMedium25_jecUp_Edge", &nBJetMedium25_jecUp_Edge, &b_nBJetMedium25_jecUp_Edge);
   fChain->SetBranchAddress("nBJetMedium25_jecDn_Edge", &nBJetMedium25_jecDn_Edge, &b_nBJetMedium25_jecDn_Edge);
   fChain->SetBranchAddress("nBJetLoose35_Edge", &nBJetLoose35_Edge, &b_nBJetLoose35_Edge);
   fChain->SetBranchAddress("nBJetLoose35_jecUp_Edge", &nBJetLoose35_jecUp_Edge, &b_nBJetLoose35_jecUp_Edge);
   fChain->SetBranchAddress("nBJetLoose35_jecDn_Edge", &nBJetLoose35_jecDn_Edge, &b_nBJetLoose35_jecDn_Edge);
   fChain->SetBranchAddress("nBJetMedium35_Edge", &nBJetMedium35_Edge, &b_nBJetMedium35_Edge);
   fChain->SetBranchAddress("nBJetMedium35_jecUp_Edge", &nBJetMedium35_jecUp_Edge, &b_nBJetMedium35_jecUp_Edge);
   fChain->SetBranchAddress("nBJetMedium35_jecDn_Edge", &nBJetMedium35_jecDn_Edge, &b_nBJetMedium35_jecDn_Edge);
   fChain->SetBranchAddress("iL1T_Edge", &iL1T_Edge, &b_iL1T_Edge);
   fChain->SetBranchAddress("iL2T_Edge", &iL2T_Edge, &b_iL2T_Edge);
   fChain->SetBranchAddress("lepsMll_Edge", &lepsMll_Edge, &b_lepsMll_Edge);
   fChain->SetBranchAddress("lepsJZB_Edge", &lepsJZB_Edge, &b_lepsJZB_Edge);
   fChain->SetBranchAddress("lepsJZB_raw_Edge", &lepsJZB_raw_Edge, &b_lepsJZB_raw_Edge);
   fChain->SetBranchAddress("lepsJZB_recoil_Edge", &lepsJZB_recoil_Edge, &b_lepsJZB_recoil_Edge);
   fChain->SetBranchAddress("lepsDR_Edge", &lepsDR_Edge, &b_lepsDR_Edge);
   fChain->SetBranchAddress("lepsMETRec_Edge", &lepsMETRec_Edge, &b_lepsMETRec_Edge);
   fChain->SetBranchAddress("lepsZPt_Edge", &lepsZPt_Edge, &b_lepsZPt_Edge);
   fChain->SetBranchAddress("metl1DPhi_Edge", &metl1DPhi_Edge, &b_metl1DPhi_Edge);
   fChain->SetBranchAddress("metl2DPhi_Edge", &metl2DPhi_Edge, &b_metl2DPhi_Edge);
   fChain->SetBranchAddress("met_Edge", &met_Edge, &b_met_Edge);
   fChain->SetBranchAddress("met_phi_Edge", &met_phi_Edge, &b_met_phi_Edge);
   fChain->SetBranchAddress("met_jecUp_Edge", &met_jecUp_Edge, &b_met_jecUp_Edge);
   fChain->SetBranchAddress("met_jecDn_Edge", &met_jecDn_Edge, &b_met_jecDn_Edge);
   fChain->SetBranchAddress("met_raw_Edge", &met_raw_Edge, &b_met_raw_Edge);
   fChain->SetBranchAddress("genMet_Edge", &genMet_Edge, &b_genMet_Edge);
   fChain->SetBranchAddress("genMet_phi_Edge", &genMet_phi_Edge, &b_genMet_phi_Edge);
   fChain->SetBranchAddress("lepsDPhi_Edge", &lepsDPhi_Edge, &b_lepsDPhi_Edge);
   fChain->SetBranchAddress("Lep1_pt_Edge", &Lep1_pt_Edge, &b_Lep1_pt_Edge);
   fChain->SetBranchAddress("Lep1_eta_Edge", &Lep1_eta_Edge, &b_Lep1_eta_Edge);
   fChain->SetBranchAddress("Lep1_phi_Edge", &Lep1_phi_Edge, &b_Lep1_phi_Edge);
   fChain->SetBranchAddress("Lep1_miniRelIso_Edge", &Lep1_miniRelIso_Edge, &b_Lep1_miniRelIso_Edge);
   fChain->SetBranchAddress("Lep1_relIso03_Edge", &Lep1_relIso03_Edge, &b_Lep1_relIso03_Edge);
   fChain->SetBranchAddress("Lep1_relIso04_Edge", &Lep1_relIso04_Edge, &b_Lep1_relIso04_Edge);
   fChain->SetBranchAddress("Lep1_dxy_Edge", &Lep1_dxy_Edge, &b_Lep1_dxy_Edge);
   fChain->SetBranchAddress("Lep1_dz_Edge", &Lep1_dz_Edge, &b_Lep1_dz_Edge);
   fChain->SetBranchAddress("Lep1_sip3d_Edge", &Lep1_sip3d_Edge, &b_Lep1_sip3d_Edge);
   fChain->SetBranchAddress("Lep1_pdgId_Edge", &Lep1_pdgId_Edge, &b_Lep1_pdgId_Edge);
   fChain->SetBranchAddress("Lep1_tightCharge_Edge", &Lep1_tightCharge_Edge, &b_Lep1_tightCharge_Edge);
   fChain->SetBranchAddress("Lep1_mvaIdSpring15_Edge", &Lep1_mvaIdSpring15_Edge, &b_Lep1_mvaIdSpring15_Edge);
   fChain->SetBranchAddress("Lep1_mcMatchId_Edge", &Lep1_mcMatchId_Edge, &b_Lep1_mcMatchId_Edge);
   fChain->SetBranchAddress("Lep1_mcMatchTau_Edge", &Lep1_mcMatchTau_Edge, &b_Lep1_mcMatchTau_Edge);
   fChain->SetBranchAddress("Lep1_minTauDR_Edge", &Lep1_minTauDR_Edge, &b_Lep1_minTauDR_Edge);
   fChain->SetBranchAddress("Lep2_pt_Edge", &Lep2_pt_Edge, &b_Lep2_pt_Edge);
   fChain->SetBranchAddress("Lep2_eta_Edge", &Lep2_eta_Edge, &b_Lep2_eta_Edge);
   fChain->SetBranchAddress("Lep2_phi_Edge", &Lep2_phi_Edge, &b_Lep2_phi_Edge);
   fChain->SetBranchAddress("Lep2_miniRelIso_Edge", &Lep2_miniRelIso_Edge, &b_Lep2_miniRelIso_Edge);
   fChain->SetBranchAddress("Lep2_relIso03_Edge", &Lep2_relIso03_Edge, &b_Lep2_relIso03_Edge);
   fChain->SetBranchAddress("Lep2_relIso04_Edge", &Lep2_relIso04_Edge, &b_Lep2_relIso04_Edge);
   fChain->SetBranchAddress("Lep2_dxy_Edge", &Lep2_dxy_Edge, &b_Lep2_dxy_Edge);
   fChain->SetBranchAddress("Lep2_dz_Edge", &Lep2_dz_Edge, &b_Lep2_dz_Edge);
   fChain->SetBranchAddress("Lep2_sip3d_Edge", &Lep2_sip3d_Edge, &b_Lep2_sip3d_Edge);
   fChain->SetBranchAddress("Lep2_pdgId_Edge", &Lep2_pdgId_Edge, &b_Lep2_pdgId_Edge);
   fChain->SetBranchAddress("Lep2_tightCharge_Edge", &Lep2_tightCharge_Edge, &b_Lep2_tightCharge_Edge);
   fChain->SetBranchAddress("Lep2_mvaIdSpring15_Edge", &Lep2_mvaIdSpring15_Edge, &b_Lep2_mvaIdSpring15_Edge);
   fChain->SetBranchAddress("Lep2_mcMatchId_Edge", &Lep2_mcMatchId_Edge, &b_Lep2_mcMatchId_Edge);
   fChain->SetBranchAddress("Lep2_mcMatchTau_Edge", &Lep2_mcMatchTau_Edge, &b_Lep2_mcMatchTau_Edge);
   fChain->SetBranchAddress("Lep2_minTauDR_Edge", &Lep2_minTauDR_Edge, &b_Lep2_minTauDR_Edge);
   fChain->SetBranchAddress("PileupW_Edge", &PileupW_Edge, &b_PileupW_Edge);
   fChain->SetBranchAddress("min_mlb1_Edge", &min_mlb1_Edge, &b_min_mlb1_Edge);
   fChain->SetBranchAddress("min_mlb2_Edge", &min_mlb2_Edge, &b_min_mlb2_Edge);
   fChain->SetBranchAddress("min_mlb1Up_Edge", &min_mlb1Up_Edge, &b_min_mlb1Up_Edge);
   fChain->SetBranchAddress("min_mlb2Up_Edge", &min_mlb2Up_Edge, &b_min_mlb2Up_Edge);
   fChain->SetBranchAddress("min_mlb1Dn_Edge", &min_mlb1Dn_Edge, &b_min_mlb1Dn_Edge);
   fChain->SetBranchAddress("min_mlb2Dn_Edge", &min_mlb2Dn_Edge, &b_min_mlb2Dn_Edge);
   fChain->SetBranchAddress("sum_mlb_Edge", &sum_mlb_Edge, &b_sum_mlb_Edge);
   fChain->SetBranchAddress("sum_mlbUp_Edge", &sum_mlbUp_Edge, &b_sum_mlbUp_Edge);
   fChain->SetBranchAddress("sum_mlbDn_Edge", &sum_mlbDn_Edge, &b_sum_mlbDn_Edge);
   fChain->SetBranchAddress("st_Edge", &st_Edge, &b_st_Edge);
   fChain->SetBranchAddress("srID_Edge", &srID_Edge, &b_srID_Edge);
   fChain->SetBranchAddress("mt2_Edge", &mt2_Edge, &b_mt2_Edge);
   fChain->SetBranchAddress("mt2_jecUp_Edge", &mt2_jecUp_Edge, &b_mt2_jecUp_Edge);
   fChain->SetBranchAddress("mZ1_Edge", &mZ1_Edge, &b_mZ1_Edge);
   fChain->SetBranchAddress("mZ2_Edge", &mZ2_Edge, &b_mZ2_Edge);
   fChain->SetBranchAddress("mt2bb_Edge", &mt2bb_Edge, &b_mt2bb_Edge);
   fChain->SetBranchAddress("mt2bb_jecUp_Edge", &mt2bb_jecUp_Edge, &b_mt2bb_jecUp_Edge);
   fChain->SetBranchAddress("mt2bb_jecDn_Edge", &mt2bb_jecDn_Edge, &b_mt2bb_jecDn_Edge);
   fChain->SetBranchAddress("mbb_Edge", &mbb_Edge, &b_mbb_Edge);
   fChain->SetBranchAddress("mbb_jecUp_Edge", &mbb_jecUp_Edge, &b_mbb_jecUp_Edge);
   fChain->SetBranchAddress("mbb_jecDn_Edge", &mbb_jecDn_Edge, &b_mbb_jecDn_Edge);          
   fChain->SetBranchAddress("weight_trigger_Edge", &weight_trigger_Edge, &b_weight_trigger_Edge);
   fChain->SetBranchAddress("weight_btagsf_Edge", &weight_btagsf_Edge, &b_weight_btagsf_Edge);
   fChain->SetBranchAddress("weight_btagsf_heavy_UP_Edge", &weight_btagsf_heavy_UP_Edge, &b_weight_btagsf_heavy_UP_Edge);
   fChain->SetBranchAddress("weight_btagsf_heavy_DN_Edge", &weight_btagsf_heavy_DN_Edge, &b_weight_btagsf_heavy_DN_Edge);
   fChain->SetBranchAddress("weight_btagsf_light_UP_Edge", &weight_btagsf_light_UP_Edge, &b_weight_btagsf_light_UP_Edge);
   fChain->SetBranchAddress("weight_btagsf_light_DN_Edge", &weight_btagsf_light_DN_Edge, &b_weight_btagsf_light_DN_Edge);
   fChain->SetBranchAddress("d3D_Edge", &d3D_Edge, &b_d3D_Edge);
   fChain->SetBranchAddress("parPt_Edge", &parPt_Edge, &b_parPt_Edge);
   fChain->SetBranchAddress("ortPt_Edge", &ortPt_Edge, &b_ortPt_Edge);
   fChain->SetBranchAddress("dTheta_Edge", &dTheta_Edge, &b_dTheta_Edge);
   fChain->SetBranchAddress("genWeight_Edge", &genWeight_Edge, &b_genWeight_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v_Edge", &HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v_Edge, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_Edge", &HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_Edge, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v_Edge", &HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v_Edge, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v_Edge", &HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v_Edge, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge", &HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge", &HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge, &b_HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v_Edge", &HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v_Edge, &b_HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu27_TkMu8_v_Edge", &HLT_BIT_HLT_Mu27_TkMu8_v_Edge, &b_HLT_BIT_HLT_Mu27_TkMu8_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge", &HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge, &b_HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge", &HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge, &b_HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v_Edge", &HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v_Edge, &b_HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v_Edge", &HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v_Edge, &b_HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu30_TkMu11_v_Edge", &HLT_BIT_HLT_Mu30_TkMu11_v_Edge, &b_HLT_BIT_HLT_Mu30_TkMu11_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v_Edge", &HLT_BIT_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v_Edge, &b_HLT_BIT_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v_Edge", &HLT_BIT_HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v_Edge, &b_HLT_BIT_HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v_Edge", &HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v_Edge, &b_HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v_Edge", &HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v_Edge, &b_HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge", &HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge, &b_HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge", &HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge, &b_HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge", &HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge, &b_HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT200_v_Edge", &HLT_BIT_HLT_PFHT200_v_Edge, &b_HLT_BIT_HLT_PFHT200_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT250_v_Edge", &HLT_BIT_HLT_PFHT250_v_Edge, &b_HLT_BIT_HLT_PFHT250_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT300_v_Edge", &HLT_BIT_HLT_PFHT300_v_Edge, &b_HLT_BIT_HLT_PFHT300_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT350_v_Edge", &HLT_BIT_HLT_PFHT350_v_Edge, &b_HLT_BIT_HLT_PFHT350_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT400_v_Edge", &HLT_BIT_HLT_PFHT400_v_Edge, &b_HLT_BIT_HLT_PFHT400_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT475_v_Edge", &HLT_BIT_HLT_PFHT475_v_Edge, &b_HLT_BIT_HLT_PFHT475_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT600_v_Edge", &HLT_BIT_HLT_PFHT600_v_Edge, &b_HLT_BIT_HLT_PFHT600_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT650_v_Edge", &HLT_BIT_HLT_PFHT650_v_Edge, &b_HLT_BIT_HLT_PFHT650_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT800_v_Edge", &HLT_BIT_HLT_PFHT800_v_Edge, &b_HLT_BIT_HLT_PFHT800_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT300_PFMET110_v_Edge", &HLT_BIT_HLT_PFHT300_PFMET110_v_Edge, &b_HLT_BIT_HLT_PFHT300_PFMET110_v_Edge);
   fChain->SetBranchAddress("GenSusyMScan1_Edge", &GenSusyMScan1_Edge, &b_GenSusyMScan1_Edge);
   fChain->SetBranchAddress("GenSusyMScan2_Edge", &GenSusyMScan2_Edge, &b_GenSusyMScan2_Edge);
   fChain->SetBranchAddress("GenSusyMScan3_Edge", &GenSusyMScan3_Edge, &b_GenSusyMScan3_Edge);
   fChain->SetBranchAddress("GenSusyMScan4_Edge", &GenSusyMScan4_Edge, &b_GenSusyMScan4_Edge);
   fChain->SetBranchAddress("GenSusyMGluino_Edge", &GenSusyMGluino_Edge, &b_GenSusyMGluino_Edge);
   fChain->SetBranchAddress("GenSusyMGravitino_Edge", &GenSusyMGravitino_Edge, &b_GenSusyMGravitino_Edge);
   fChain->SetBranchAddress("GenSusyMStop_Edge", &GenSusyMStop_Edge, &b_GenSusyMStop_Edge);
   fChain->SetBranchAddress("GenSusyMSbottom_Edge", &GenSusyMSbottom_Edge, &b_GenSusyMSbottom_Edge);
   fChain->SetBranchAddress("GenSusyMStop2_Edge", &GenSusyMStop2_Edge, &b_GenSusyMStop2_Edge);
   fChain->SetBranchAddress("GenSusyMSbottom2_Edge", &GenSusyMSbottom2_Edge, &b_GenSusyMSbottom2_Edge);
   fChain->SetBranchAddress("GenSusyMSquark_Edge", &GenSusyMSquark_Edge, &b_GenSusyMSquark_Edge);
   fChain->SetBranchAddress("GenSusyMNeutralino_Edge", &GenSusyMNeutralino_Edge, &b_GenSusyMNeutralino_Edge);
   fChain->SetBranchAddress("GenSusyMNeutralino2_Edge", &GenSusyMNeutralino2_Edge, &b_GenSusyMNeutralino2_Edge);
   fChain->SetBranchAddress("GenSusyMNeutralino3_Edge", &GenSusyMNeutralino3_Edge, &b_GenSusyMNeutralino3_Edge);
   fChain->SetBranchAddress("GenSusyMNeutralino4_Edge", &GenSusyMNeutralino4_Edge, &b_GenSusyMNeutralino4_Edge);
   fChain->SetBranchAddress("GenSusyMChargino_Edge", &GenSusyMChargino_Edge, &b_GenSusyMChargino_Edge);
   fChain->SetBranchAddress("GenSusyMChargino2_Edge", &GenSusyMChargino2_Edge, &b_GenSusyMChargino2_Edge);
   fChain->SetBranchAddress("JetSel_Edge_pt", JetSel_Edge_pt, &b_JetSel_Edge_pt);
   fChain->SetBranchAddress("JetSel_Edge_eta", JetSel_Edge_eta, &b_JetSel_Edge_eta);
   fChain->SetBranchAddress("JetSel_Edge_phi", JetSel_Edge_phi, &b_JetSel_Edge_phi);
   fChain->SetBranchAddress("JetSel_Edge_mass", JetSel_Edge_mass, &b_JetSel_Edge_mass);
   fChain->SetBranchAddress("JetSel_Edge_btagCSV", JetSel_Edge_btagCSV, &b_JetSel_Edge_btagCSV);
   fChain->SetBranchAddress("JetSel_Edge_rawPt", JetSel_Edge_rawPt, &b_JetSel_Edge_rawPt);
   fChain->SetBranchAddress("JetSel_Edge_mcPt", JetSel_Edge_mcPt, &b_JetSel_Edge_mcPt);
   fChain->SetBranchAddress("JetSel_Edge_mcFlavour", JetSel_Edge_mcFlavour, &b_JetSel_Edge_mcFlavour);
   fChain->SetBranchAddress("JetSel_Edge_mcMatchId", JetSel_Edge_mcMatchId, &b_JetSel_Edge_mcMatchId);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge", &HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge,&b_HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge", &HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge, &b_HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge", &HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge, &b_HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge);
   Notify();
}

void skimmer::SetOutVariables()
{
   outputtree->SetBranchAddress("evt_Edge", &evt_Edge, &b_evt_Edge);
   outputtree->SetBranchAddress("run_Edge", &run_Edge, &b_run_Edge);
   outputtree->SetBranchAddress("lumi_Edge", &lumi_Edge, &b_lumi_Edge);
   outputtree->SetBranchAddress("nVert_Edge", &nVert_Edge, &b_nVert_Edge);
   outputtree->SetBranchAddress("Flag_HBHENoiseFilter_Edge", &Flag_HBHENoiseFilter_Edge, &b_Flag_HBHENoiseFilter_Edge); 
   outputtree->SetBranchAddress("Flag_HBHENoiseIsoFilter_Edge", &Flag_HBHENoiseIsoFilter_Edge, &b_Flag_HBHENoiseIsoFilter_Edge);
   outputtree->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter_Edge", &Flag_EcalDeadCellTriggerPrimitiveFilter_Edge, &b_Flag_EcalDeadCellTriggerPrimitiveFilter_Edge);
   outputtree->SetBranchAddress("Flag_goodVertices_Edge", &Flag_goodVertices_Edge, &b_Flag_goodVertices_Edge);
   outputtree->SetBranchAddress("Flag_eeBadScFilter_Edge", &Flag_eeBadScFilter_Edge, &b_Flag_eeBadScFilter_Edge);
   outputtree->SetBranchAddress("Flag_globalTightHalo2016Filter_Edge", &Flag_globalTightHalo2016Filter_Edge, &b_Flag_globalTightHalo2016Filter_Edge);
   outputtree->SetBranchAddress("Flag_CSCTightHalo2016Filter_Edge", &Flag_CSCTightHalo2016Filter_Edge, &b_Flag_CSCTightHalo2016Filter_Edge);
   outputtree->SetBranchAddress("Flag_badMuonFilter_Edge", &Flag_badMuonFilter_Edge, &b_Flag_badMuonFilter_Edge);                                                                  
   outputtree->SetBranchAddress("Flag_badChargedHadronFilter_Edge", &Flag_badChargedHadronFilter_Edge, &b_Flag_badChargedHadronFilter_Edge);                                              
   outputtree->SetBranchAddress("Flag_badCloneMuonMoriond2017_Edge", &Flag_badCloneMuonMoriond2017_Edge, &b_Flag_badCloneMuonMoriond2017_Edge);                  
   outputtree->SetBranchAddress("Flag_badMuonMoriond2017_Edge", &Flag_badMuonMoriond2017_Edge, &b_Flag_badMuonMoriond2017_Edge);                                 
   outputtree->SetBranchAddress("badCloneMuonMoriond2017_maxPt_Edge", &badCloneMuonMoriond2017_maxPt_Edge, &b_badCloneMuonMoriond2017_maxPt_Edge);               
   outputtree->SetBranchAddress("badNotCloneMuonMoriond2017_maxPt_Edge", &badNotCloneMuonMoriond2017_maxPt_Edge, &b_badNotCloneMuonMoriond2017_maxPt_Edge);      
   outputtree->SetBranchAddress("nLepTight_Edge", &nLepTight_Edge, &b_nLepTight_Edge);
   outputtree->SetBranchAddress("nLepLoose_Edge", &nLepLoose_Edge, &b_nLepLoose_Edge);
   outputtree->SetBranchAddress("nPFHad10_Edge", &nPFHad10_Edge, &b_nPFHad10_Edge);
   outputtree->SetBranchAddress("nPFLep5_Edge", &nPFLep5_Edge, &b_nPFLep5_Edge);
   outputtree->SetBranchAddress("nJetSel_Edge", &nJetSel_Edge, &b_nJetSel_Edge);
   outputtree->SetBranchAddress("nJetSel_jecUp_Edge", &nJetSel_jecUp_Edge, &b_nJetSel_jecUp_Edge);
   outputtree->SetBranchAddress("nJetSel_jecDn_Edge", &nJetSel_jecDn_Edge, &b_nJetSel_jecDn_Edge);
   outputtree->SetBranchAddress("bestMjj_Edge", &bestMjj_Edge, &b_bestMjj_Edge);
   outputtree->SetBranchAddress("minMjj_Edge", &minMjj_Edge, &b_minMjj_Edge);
   outputtree->SetBranchAddress("maxMjj_Edge", &maxMjj_Edge, &b_maxMjj_Edge);
   outputtree->SetBranchAddress("hardMjj_Edge", &hardMjj_Edge, &b_hardMjj_Edge);
   outputtree->SetBranchAddress("dphiMjj_Edge", &dphiMjj_Edge, &b_dphiMjj_Edge);
   outputtree->SetBranchAddress("drMjj_Edge", &drMjj_Edge, &b_drMjj_Edge);
   outputtree->SetBranchAddress("hardJJDphi_Edge", &hardJJDphi_Edge, &b_hardJJDphi_Edge);
   outputtree->SetBranchAddress("hardJJDR_Edge", &hardJJDR_Edge, &b_hardJJDR_Edge);
   outputtree->SetBranchAddress("j1MetDPhi_Edge", &j1MetDPhi_Edge, &b_j1MetDPhi_Edge);
   outputtree->SetBranchAddress("j2MetDPhi_Edge", &j2MetDPhi_Edge, &b_j2MetDPhi_Edge);
   outputtree->SetBranchAddress("nPairLep_Edge", &nPairLep_Edge, &b_nPairLep_Edge);
   outputtree->SetBranchAddress("iLT_Edge", iLT_Edge, &b_iLT_Edge);
   outputtree->SetBranchAddress("iJ_Edge", iJ_Edge, &b_iJ_Edge);
   outputtree->SetBranchAddress("nLepGood20_Edge", &nLepGood20_Edge, &b_nLepGood20_Edge);
   outputtree->SetBranchAddress("nLepGood20T_Edge", &nLepGood20T_Edge, &b_nLepGood20T_Edge);
   outputtree->SetBranchAddress("nJet35_Edge", &nJet35_Edge, &b_nJet35_Edge);
   outputtree->SetBranchAddress("nJet35_jecUp_Edge", &nJet35_jecUp_Edge, &b_nJet35_jecUp_Edge);
   outputtree->SetBranchAddress("nJet35_jecDn_Edge", &nJet35_jecDn_Edge, &b_nJet35_jecDn_Edge);
   outputtree->SetBranchAddress("htJet35j_Edge", &htJet35j_Edge, &b_htJet35j_Edge);
   outputtree->SetBranchAddress("htJet35j_jecUp_Edge", &htJet35j_jecUp_Edge, &b_htJet35j_jecUp_Edge);
   outputtree->SetBranchAddress("htJet35j_jecDn_Edge", &htJet35j_jecDn_Edge, &b_htJet35j_jecDn_Edge);
   outputtree->SetBranchAddress("nBJetMedium25_Edge", &nBJetMedium25_Edge, &b_nBJetMedium25_Edge);
   outputtree->SetBranchAddress("nBJetMedium25_jecUp_Edge", &nBJetMedium25_jecUp_Edge, &b_nBJetMedium25_jecUp_Edge);
   outputtree->SetBranchAddress("nBJetMedium25_jecDn_Edge", &nBJetMedium25_jecDn_Edge, &b_nBJetMedium25_jecDn_Edge);
   outputtree->SetBranchAddress("nBJetLoose35_Edge", &nBJetLoose35_Edge, &b_nBJetLoose35_Edge);
   outputtree->SetBranchAddress("nBJetLoose35_jecUp_Edge", &nBJetLoose35_jecUp_Edge, &b_nBJetLoose35_jecUp_Edge);
   outputtree->SetBranchAddress("nBJetLoose35_jecDn_Edge", &nBJetLoose35_jecDn_Edge, &b_nBJetLoose35_jecDn_Edge);
   outputtree->SetBranchAddress("nBJetMedium35_Edge", &nBJetMedium35_Edge, &b_nBJetMedium35_Edge);
   outputtree->SetBranchAddress("nBJetMedium35_jecUp_Edge", &nBJetMedium35_jecUp_Edge, &b_nBJetMedium35_jecUp_Edge);
   outputtree->SetBranchAddress("nBJetMedium35_jecDn_Edge", &nBJetMedium35_jecDn_Edge, &b_nBJetMedium35_jecDn_Edge);
   outputtree->SetBranchAddress("iL1T_Edge", &iL1T_Edge, &b_iL1T_Edge);
   outputtree->SetBranchAddress("iL2T_Edge", &iL2T_Edge, &b_iL2T_Edge);
   outputtree->SetBranchAddress("lepsMll_Edge", &lepsMll_Edge, &b_lepsMll_Edge);
   outputtree->SetBranchAddress("lepsJZB_Edge", &lepsJZB_Edge, &b_lepsJZB_Edge);
   outputtree->SetBranchAddress("lepsJZB_raw_Edge", &lepsJZB_raw_Edge, &b_lepsJZB_raw_Edge);
   outputtree->SetBranchAddress("lepsJZB_recoil_Edge", &lepsJZB_recoil_Edge, &b_lepsJZB_recoil_Edge);
   outputtree->SetBranchAddress("lepsDR_Edge", &lepsDR_Edge, &b_lepsDR_Edge);
   outputtree->SetBranchAddress("lepsMETRec_Edge", &lepsMETRec_Edge, &b_lepsMETRec_Edge);
   outputtree->SetBranchAddress("lepsZPt_Edge", &lepsZPt_Edge, &b_lepsZPt_Edge);
   outputtree->SetBranchAddress("metl1DPhi_Edge", &metl1DPhi_Edge, &b_metl1DPhi_Edge);
   outputtree->SetBranchAddress("metl2DPhi_Edge", &metl2DPhi_Edge, &b_metl2DPhi_Edge);
   outputtree->SetBranchAddress("met_Edge", &met_Edge, &b_met_Edge);
   outputtree->SetBranchAddress("met_phi_Edge", &met_phi_Edge, &b_met_phi_Edge);
   outputtree->SetBranchAddress("met_jecUp_Edge", &met_jecUp_Edge, &b_met_jecUp_Edge);
   outputtree->SetBranchAddress("met_jecDn_Edge", &met_jecDn_Edge, &b_met_jecDn_Edge);
   outputtree->SetBranchAddress("met_raw_Edge", &met_raw_Edge, &b_met_raw_Edge);
   outputtree->SetBranchAddress("genMet_Edge", &genMet_Edge, &b_genMet_Edge);
   outputtree->SetBranchAddress("genMet_phi_Edge", &genMet_phi_Edge, &b_genMet_phi_Edge);
   outputtree->SetBranchAddress("lepsDPhi_Edge", &lepsDPhi_Edge, &b_lepsDPhi_Edge);
   outputtree->SetBranchAddress("Lep1_pt_Edge", &Lep1_pt_Edge, &b_Lep1_pt_Edge);
   outputtree->SetBranchAddress("Lep1_eta_Edge", &Lep1_eta_Edge, &b_Lep1_eta_Edge);
   outputtree->SetBranchAddress("Lep1_phi_Edge", &Lep1_phi_Edge, &b_Lep1_phi_Edge);
   outputtree->SetBranchAddress("Lep1_miniRelIso_Edge", &Lep1_miniRelIso_Edge, &b_Lep1_miniRelIso_Edge);
   outputtree->SetBranchAddress("Lep1_relIso03_Edge", &Lep1_relIso03_Edge, &b_Lep1_relIso03_Edge);
   outputtree->SetBranchAddress("Lep1_relIso04_Edge", &Lep1_relIso04_Edge, &b_Lep1_relIso04_Edge);
   outputtree->SetBranchAddress("Lep1_dxy_Edge", &Lep1_dxy_Edge, &b_Lep1_dxy_Edge);
   outputtree->SetBranchAddress("Lep1_dz_Edge", &Lep1_dz_Edge, &b_Lep1_dz_Edge);
   outputtree->SetBranchAddress("Lep1_sip3d_Edge", &Lep1_sip3d_Edge, &b_Lep1_sip3d_Edge);
   outputtree->SetBranchAddress("Lep1_pdgId_Edge", &Lep1_pdgId_Edge, &b_Lep1_pdgId_Edge);
   outputtree->SetBranchAddress("Lep1_tightCharge_Edge", &Lep1_tightCharge_Edge, &b_Lep1_tightCharge_Edge);
   outputtree->SetBranchAddress("Lep1_mvaIdSpring15_Edge", &Lep1_mvaIdSpring15_Edge, &b_Lep1_mvaIdSpring15_Edge);
   outputtree->SetBranchAddress("Lep1_mcMatchId_Edge", &Lep1_mcMatchId_Edge, &b_Lep1_mcMatchId_Edge);
   outputtree->SetBranchAddress("Lep1_mcMatchTau_Edge", &Lep1_mcMatchTau_Edge, &b_Lep1_mcMatchTau_Edge);
   outputtree->SetBranchAddress("Lep1_minTauDR_Edge", &Lep1_minTauDR_Edge, &b_Lep1_minTauDR_Edge);
   outputtree->SetBranchAddress("Lep2_pt_Edge", &Lep2_pt_Edge, &b_Lep2_pt_Edge);
   outputtree->SetBranchAddress("Lep2_eta_Edge", &Lep2_eta_Edge, &b_Lep2_eta_Edge);
   outputtree->SetBranchAddress("Lep2_phi_Edge", &Lep2_phi_Edge, &b_Lep2_phi_Edge);
   outputtree->SetBranchAddress("Lep2_miniRelIso_Edge", &Lep2_miniRelIso_Edge, &b_Lep2_miniRelIso_Edge);
   outputtree->SetBranchAddress("Lep2_relIso03_Edge", &Lep2_relIso03_Edge, &b_Lep2_relIso03_Edge);
   outputtree->SetBranchAddress("Lep2_relIso04_Edge", &Lep2_relIso04_Edge, &b_Lep2_relIso04_Edge);
   outputtree->SetBranchAddress("Lep2_dxy_Edge", &Lep2_dxy_Edge, &b_Lep2_dxy_Edge);
   outputtree->SetBranchAddress("Lep2_dz_Edge", &Lep2_dz_Edge, &b_Lep2_dz_Edge);
   outputtree->SetBranchAddress("Lep2_sip3d_Edge", &Lep2_sip3d_Edge, &b_Lep2_sip3d_Edge);
   outputtree->SetBranchAddress("Lep2_pdgId_Edge", &Lep2_pdgId_Edge, &b_Lep2_pdgId_Edge);
   outputtree->SetBranchAddress("Lep2_tightCharge_Edge", &Lep2_tightCharge_Edge, &b_Lep2_tightCharge_Edge);
   outputtree->SetBranchAddress("Lep2_mvaIdSpring15_Edge", &Lep2_mvaIdSpring15_Edge, &b_Lep2_mvaIdSpring15_Edge);
   outputtree->SetBranchAddress("Lep2_mcMatchId_Edge", &Lep2_mcMatchId_Edge, &b_Lep2_mcMatchId_Edge);
   outputtree->SetBranchAddress("Lep2_mcMatchTau_Edge", &Lep2_mcMatchTau_Edge, &b_Lep2_mcMatchTau_Edge);
   outputtree->SetBranchAddress("Lep2_minTauDR_Edge", &Lep2_minTauDR_Edge, &b_Lep2_minTauDR_Edge);
   outputtree->SetBranchAddress("PileupW_Edge", &PileupW_Edge, &b_PileupW_Edge);
   outputtree->SetBranchAddress("min_mlb1_Edge", &min_mlb1_Edge, &b_min_mlb1_Edge);
   outputtree->SetBranchAddress("min_mlb2_Edge", &min_mlb2_Edge, &b_min_mlb2_Edge);
   outputtree->SetBranchAddress("min_mlb1Up_Edge", &min_mlb1Up_Edge, &b_min_mlb1Up_Edge);
   outputtree->SetBranchAddress("min_mlb2Up_Edge", &min_mlb2Up_Edge, &b_min_mlb2Up_Edge);
   outputtree->SetBranchAddress("min_mlb1Dn_Edge", &min_mlb1Dn_Edge, &b_min_mlb1Dn_Edge);
   outputtree->SetBranchAddress("min_mlb2Dn_Edge", &min_mlb2Dn_Edge, &b_min_mlb2Dn_Edge);
   outputtree->SetBranchAddress("sum_mlb_Edge", &sum_mlb_Edge, &b_sum_mlb_Edge);
   outputtree->SetBranchAddress("sum_mlbUp_Edge", &sum_mlbUp_Edge, &b_sum_mlbUp_Edge);
   outputtree->SetBranchAddress("sum_mlbDn_Edge", &sum_mlbDn_Edge, &b_sum_mlbDn_Edge);
   outputtree->SetBranchAddress("st_Edge", &st_Edge, &b_st_Edge);
   outputtree->SetBranchAddress("srID_Edge", &srID_Edge, &b_srID_Edge);
   outputtree->SetBranchAddress("mt2_Edge", &mt2_Edge, &b_mt2_Edge);
   outputtree->SetBranchAddress("mt2_jecUp_Edge", &mt2_jecUp_Edge, &b_mt2_jecUp_Edge);
   outputtree->SetBranchAddress("mt2_jecDn_Edge", &mt2_jecDn_Edge, &b_mt2_jecDn_Edge);
   outputtree->SetBranchAddress("mZ1_Edge", &mZ1_Edge, &b_mZ1_Edge);
   outputtree->SetBranchAddress("mZ2_Edge", &mZ2_Edge, &b_mZ2_Edge);
   outputtree->SetBranchAddress("mt2bb_Edge", &mt2bb_Edge, &b_mt2bb_Edge);
   outputtree->SetBranchAddress("mt2bb_jecUp_Edge", &mt2bb_jecUp_Edge, &b_mt2bb_jecUp_Edge);
   outputtree->SetBranchAddress("mt2bb_jecDn_Edge", &mt2bb_jecDn_Edge, &b_mt2bb_jecDn_Edge);
   outputtree->SetBranchAddress("mbb_Edge", &mbb_Edge, &b_mbb_Edge);
   outputtree->SetBranchAddress("mbb_jecUp_Edge", &mbb_jecUp_Edge, &b_mbb_jecUp_Edge);
   outputtree->SetBranchAddress("mbb_jecDn_Edge", &mbb_jecDn_Edge, &b_mbb_jecDn_Edge);          
   outputtree->SetBranchAddress("weight_trigger_Edge", &weight_trigger_Edge, &b_weight_trigger_Edge);
   outputtree->SetBranchAddress("weight_btagsf_Edge", &weight_btagsf_Edge, &b_weight_btagsf_Edge);
   outputtree->SetBranchAddress("weight_btagsf_heavy_UP_Edge", &weight_btagsf_heavy_UP_Edge, &b_weight_btagsf_heavy_UP_Edge);
   outputtree->SetBranchAddress("weight_btagsf_heavy_DN_Edge", &weight_btagsf_heavy_DN_Edge, &b_weight_btagsf_heavy_DN_Edge);
   outputtree->SetBranchAddress("weight_btagsf_light_UP_Edge", &weight_btagsf_light_UP_Edge, &b_weight_btagsf_light_UP_Edge);
   outputtree->SetBranchAddress("weight_btagsf_light_DN_Edge", &weight_btagsf_light_DN_Edge, &b_weight_btagsf_light_DN_Edge);
   outputtree->SetBranchAddress("d3D_Edge", &d3D_Edge, &b_d3D_Edge);
   outputtree->SetBranchAddress("parPt_Edge", &parPt_Edge, &b_parPt_Edge);
   outputtree->SetBranchAddress("ortPt_Edge", &ortPt_Edge, &b_ortPt_Edge);
   outputtree->SetBranchAddress("dTheta_Edge", &dTheta_Edge, &b_dTheta_Edge);
   outputtree->SetBranchAddress("genWeight_Edge", &genWeight_Edge, &b_genWeight_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v_Edge", &HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v_Edge, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_Edge", &HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_Edge, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v_Edge", &HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v_Edge, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v_Edge", &HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v_Edge, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge", &HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge", &HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge, &b_HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v_Edge", &HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v_Edge, &b_HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu27_TkMu8_v_Edge", &HLT_BIT_HLT_Mu27_TkMu8_v_Edge, &b_HLT_BIT_HLT_Mu27_TkMu8_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge", &HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge, &b_HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge", &HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge, &b_HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v_Edge", &HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v_Edge, &b_HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v_Edge", &HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v_Edge, &b_HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu30_TkMu11_v_Edge", &HLT_BIT_HLT_Mu30_TkMu11_v_Edge, &b_HLT_BIT_HLT_Mu30_TkMu11_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v_Edge", &HLT_BIT_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v_Edge, &b_HLT_BIT_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v_Edge", &HLT_BIT_HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v_Edge, &b_HLT_BIT_HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v_Edge", &HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v_Edge, &b_HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v_Edge", &HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v_Edge, &b_HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge", &HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge, &b_HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge", &HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge, &b_HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge", &HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge, &b_HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_PFHT200_v_Edge", &HLT_BIT_HLT_PFHT200_v_Edge, &b_HLT_BIT_HLT_PFHT200_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_PFHT250_v_Edge", &HLT_BIT_HLT_PFHT250_v_Edge, &b_HLT_BIT_HLT_PFHT250_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_PFHT300_v_Edge", &HLT_BIT_HLT_PFHT300_v_Edge, &b_HLT_BIT_HLT_PFHT300_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_PFHT350_v_Edge", &HLT_BIT_HLT_PFHT350_v_Edge, &b_HLT_BIT_HLT_PFHT350_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_PFHT400_v_Edge", &HLT_BIT_HLT_PFHT400_v_Edge, &b_HLT_BIT_HLT_PFHT400_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_PFHT475_v_Edge", &HLT_BIT_HLT_PFHT475_v_Edge, &b_HLT_BIT_HLT_PFHT475_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_PFHT600_v_Edge", &HLT_BIT_HLT_PFHT600_v_Edge, &b_HLT_BIT_HLT_PFHT600_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_PFHT650_v_Edge", &HLT_BIT_HLT_PFHT650_v_Edge, &b_HLT_BIT_HLT_PFHT650_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_PFHT800_v_Edge", &HLT_BIT_HLT_PFHT800_v_Edge, &b_HLT_BIT_HLT_PFHT800_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_PFHT300_PFMET110_v_Edge", &HLT_BIT_HLT_PFHT300_PFMET110_v_Edge, &b_HLT_BIT_HLT_PFHT300_PFMET110_v_Edge);
   outputtree->SetBranchAddress("GenSusyMScan1_Edge", &GenSusyMScan1_Edge, &b_GenSusyMScan1_Edge);
   outputtree->SetBranchAddress("GenSusyMScan2_Edge", &GenSusyMScan2_Edge, &b_GenSusyMScan2_Edge);
   outputtree->SetBranchAddress("GenSusyMScan3_Edge", &GenSusyMScan3_Edge, &b_GenSusyMScan3_Edge);
   outputtree->SetBranchAddress("GenSusyMScan4_Edge", &GenSusyMScan4_Edge, &b_GenSusyMScan4_Edge);
   outputtree->SetBranchAddress("GenSusyMGluino_Edge", &GenSusyMGluino_Edge, &b_GenSusyMGluino_Edge);
   outputtree->SetBranchAddress("GenSusyMGravitino_Edge", &GenSusyMGravitino_Edge, &b_GenSusyMGravitino_Edge);
   outputtree->SetBranchAddress("GenSusyMStop_Edge", &GenSusyMStop_Edge, &b_GenSusyMStop_Edge);
   outputtree->SetBranchAddress("GenSusyMSbottom_Edge", &GenSusyMSbottom_Edge, &b_GenSusyMSbottom_Edge);
   outputtree->SetBranchAddress("GenSusyMStop2_Edge", &GenSusyMStop2_Edge, &b_GenSusyMStop2_Edge);
   outputtree->SetBranchAddress("GenSusyMSbottom2_Edge", &GenSusyMSbottom2_Edge, &b_GenSusyMSbottom2_Edge);
   outputtree->SetBranchAddress("GenSusyMSquark_Edge", &GenSusyMSquark_Edge, &b_GenSusyMSquark_Edge);
   outputtree->SetBranchAddress("GenSusyMNeutralino_Edge", &GenSusyMNeutralino_Edge, &b_GenSusyMNeutralino_Edge);
   outputtree->SetBranchAddress("GenSusyMNeutralino2_Edge", &GenSusyMNeutralino2_Edge, &b_GenSusyMNeutralino2_Edge);
   outputtree->SetBranchAddress("GenSusyMNeutralino3_Edge", &GenSusyMNeutralino3_Edge, &b_GenSusyMNeutralino3_Edge);
   outputtree->SetBranchAddress("GenSusyMNeutralino4_Edge", &GenSusyMNeutralino4_Edge, &b_GenSusyMNeutralino4_Edge);
   outputtree->SetBranchAddress("GenSusyMChargino_Edge", &GenSusyMChargino_Edge, &b_GenSusyMChargino_Edge);
   outputtree->SetBranchAddress("GenSusyMChargino2_Edge", &GenSusyMChargino2_Edge, &b_GenSusyMChargino2_Edge);
   outputtree->SetBranchAddress("JetSel_Edge_pt", JetSel_Edge_pt, &b_JetSel_Edge_pt);
   outputtree->SetBranchAddress("JetSel_Edge_eta", JetSel_Edge_eta, &b_JetSel_Edge_eta);
   outputtree->SetBranchAddress("JetSel_Edge_phi", JetSel_Edge_phi, &b_JetSel_Edge_phi);
   outputtree->SetBranchAddress("JetSel_Edge_mass", JetSel_Edge_mass, &b_JetSel_Edge_mass);
   outputtree->SetBranchAddress("JetSel_Edge_btagCSV", JetSel_Edge_btagCSV, &b_JetSel_Edge_btagCSV);
   outputtree->SetBranchAddress("JetSel_Edge_rawPt", JetSel_Edge_rawPt, &b_JetSel_Edge_rawPt);
   outputtree->SetBranchAddress("JetSel_Edge_mcPt", JetSel_Edge_mcPt, &b_JetSel_Edge_mcPt);
   outputtree->SetBranchAddress("JetSel_Edge_mcFlavour", JetSel_Edge_mcFlavour, &b_JetSel_Edge_mcFlavour);
   outputtree->SetBranchAddress("JetSel_Edge_mcMatchId", JetSel_Edge_mcMatchId, &b_JetSel_Edge_mcMatchId);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge", &HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge,&b_HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge", &HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge, &b_HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge);
   outputtree->SetBranchAddress("HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge", &HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge, &b_HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge);


}

Bool_t skimmer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void skimmer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t skimmer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef skimmer_cxx

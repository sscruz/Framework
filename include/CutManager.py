

class CutManager:
   'This class serves as an on-demand cut server'

   def __init__(self):

      self.protection = '((run_Edge<=276811)  ||  (278820<=run_Edge && run_Edge<=279931))'

      ########################################################################
      ######Basic Lepton Cuts ################################################
      ########################################################################
      self.twoLeptons = "nPairLep_Edge > 0 &&  Flag_HBHENoiseFilter_Edge ==1 && Flag_HBHENoiseIsoFilter_Edge ==1 && Flag_badChargedHadronFilter_Edge == 1 && Flag_EcalDeadCellTriggerPrimitiveFilter_Edge == 1 && Flag_goodVertices_Edge == 1 && Flag_globalTightHalo2016Filter_Edge ==1 "
      self.tightCharge = 'Lep1_tightCharge_Edge > 0 && Lep2_tightCharge_Edge > 0'
      self.leptonPt = "Lep1_pt_Edge > 25 && Lep2_pt_Edge > 20."
      self.diLeptonPt = "lepsZPt_Edge > 25"
      self.leptonDR = "lepsDR_Edge > 0.1"       
      self.leptonsMll = "lepsMll_Edge > 20"
      self.ee = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -121)"
      self.mm = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -169)"
      self.OF = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -143)"
      self.SF = "(" + self.ee + " || " +  self.mm + ")"
      self.AF = "(" + self.SF + " || " +  self.OF + ")"
      self.trigMMNEW = "( HLT_BIT_HLT_Mu30_TkMu11_v_Edge ==1 ||HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v_Edge == 1  || HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v_Edge ==1 || HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_Edge == 1 || HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v_Edge == 1 ||HLT_BIT_HLT_Mu27_TkMu8_v_Edge==1 )"
      self.trigEENEW = "(HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge == 1 || HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge==1 || HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v_Edge==1 || HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v_Edge==1)"
      self.trigEMNEW = "(HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge ==1 || HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge == 1 || HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v_Edge ==1 || HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v_Edge == 1 || HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v_Edge ==1 || HLT_BIT_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v_Edge ==1 ||  HLT_BIT_HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v_Edge ==1  || HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge == 1 || HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge == 1 || HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge == 1)"
      self.trigger = "((" + self.trigMMNEW + " && " + self.mm + ") || (" + self.trigEENEW + " && " + self.ee + ") || (" + self.trigEMNEW + " && " + self.OF + "))"
      self.triggerForCR = "((" + self.trigMMNEW + " && abs(Lep1_pdgId_Edge * Lep2_pdgId_Edge) == 169) || (" + self.trigEENEW + " && abs(Lep1_pdgId_Edge * Lep2_pdgId_Edge) == 121 ) || (" + self.trigEMNEW + " && abs(Lep1_pdgId_Edge * Lep2_pdgId_Edge) == 143))"
      self.goodLeptonInc = "("+self.trigger + "&&" +self.twoLeptons +  "&&" + self.leptonDR + ")"
      self.goodLepton = "("+self.trigger + "&&" +self.twoLeptons + "&&" + self.diLeptonPt +"&&" + self.leptonPt + "&&" + self.leptonDR + "&&"  + self.leptonsMll +  ")"
      self.goodLeptonNoTrigger = "(" +self.twoLeptons + "&&" + self.diLeptonPt +"&&" + self.leptonPt + "&&" + self.leptonDR + "&&"  + self.leptonsMll +  ")"
      self.ThirdLeptonVeto = '(nLepLoose_Edge == 2 && nPFHad10_Edge == 0 && nPFLep5_Edge <= 2)'
      self.tightIso = 'max(Lep1_miniRelIso_Edge, Lep2_miniRelIso_Edge) < 0.05'
      self.threeTightLeptons = 'nLepTight_Edge == 3'
      self.fourTightLeptons = 'nLepTight_Edge == 4'
      self.twoLooseLeptons = 'nLepLoose_Edge == 2'
      self.fourLooseLeptons = 'nLepLoose_Edge > 3'
      ########################################################################
      ######Basic Jets and bJets Cuts###############################################
      ########################################################################
      self.nj2 = "(nJetSel_Edge >= 2)"
      self.njless2 = "(nJetSel_Edge < 2)"
      self.nbjless2 = "(nBJetMedium25_Edge < 2)"
      self.nj1 = "(nJetSel_Edge >= 1)"
      self.nj0 = "(nJetSel_Edge >= 0)"
      self.nbj2 = "(nBJetMedium25_Edge >= 2)"
      self.nbj1 = "(nBJetMedium25_Edge >= 1)"
      self.nbj0 = "(nBJetMedium25_Edge >= 0)"
      self.bveto = "(nBJetMedium25_Edge == 0)"
      self.binv = "(!(nBJetMedium25_Edge == 0))"
      self.njExact2 = '(nJetSel_Edge == 2)'
      self.nbExact2 = '(nBJetMedium25_Edge == 2)'
      self.mjj110  = 'dphiMjj_Edge < 110.'
      self.mbb150  = 'mbb_Edge < 150.'
      self.mZ2g20  = 'mZ2_Edge > 20.'
      self.mZ1inZ  = 'mZ1_Edge >= 86 && mZ1_Edge < 96'
      self.ZmassInc  = 'mZ1_Edge >= 80 && mZ1_Edge < 100'
      self.lepsFromZ  = 'abs(Lep1_mcMatchId_Edge) == 23 && abs(Lep2_mcMatchId_Edge) == 23'
      
      ########################################################################
      ######Basic MET cuts##########################################################
      ########################################################################
      self.METl50  = "(met_Edge < 50)"
      self.METl100 = "(met_Edge < 100)"
      self.METl150 = "(met_Edge < 150)"
      self.METg30 = "(met_Edge > 30)"
      self.METg60 = "(met_Edge > 60)"
      self.METg50 = "(met_Edge > 50)"
      self.METg80 = "(met_Edge >  80)"
      self.METg100 = "(met_Edge >= 100)"
      self.METg150 = "(met_Edge >= 150)"
      self.METg200 = "(met_Edge >= 200)"
      self.MET50_100 = "(met_Edge >= 50 && met_Edge < 100)"
      self.MET100_150 = "(met_Edge >= 100 && met_Edge < 150)"
      self.MET150_250 = "(met_Edge >= 150 && met_Edge < 250)"
      self.MET250 = "(met_Edge >= 250)"
      self.dPhiJETMET = " (abs(j1MetDPhi_Edge)>= 0.4)&& (abs(j2MetDPhi_Edge)>= 0.4)"
      
      ########################################################################
      ######Basic NLL cut#####################################################
      ########################################################################
      self.NLL = '(nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) > 21)'

      ########################################################################
      ######Basic MT2 cuts####################################################
      ########################################################################
      self.mT2_80 = "(mt2_Edge >= 80)"
      self.mT2_100 = "(mt2_Edge >= 100)"
      self.mT2b200 = '(mt2bb_Edge >= 200.)'

      self.FSCentralJetCleaning = 'FS_central_jets_Edge > 0'
      
      ########################################################################
      ######Basic Mass cuts###################################################
      ########################################################################
      self.wideZ = "lepsMll_Edge > 76 && lepsMll_Edge < 106"
      self.mZ140 =  "mZ1_Edge < 140"
      self.mZ4_120 = "mZ2_Edge < 120 && mZ1_Edge > 4 && mZ2_Edge > 4 && lepsMll_Edge > 4 "
      self.Zmass = "lepsMll_Edge >= 86 && lepsMll_Edge < 96"
      self.mll20_60= "lepsMll_Edge < 60. && lepsMll_Edge >= 20 "
      self.mll60_86= "lepsMll_Edge < 86. && lepsMll_Edge >= 60 "
      self.mll81_101= "lepsMll_Edge < 101. && lepsMll_Edge >= 81 "
      self.mll86_96= "lepsMll_Edge < 96. && lepsMll_Edge >= 86 "
      self.mll70_110= "lepsMll_Edge < 110. && lepsMll_Edge >= 70 "
      self.mll96_150= "lepsMll_Edge < 150. && lepsMll_Edge >= 96 "
      self.mll150_200= "lepsMll_Edge < 200. && lepsMll_Edge >= 150 "
      self.mll200_300= "lepsMll_Edge < 300. && lepsMll_Edge >= 200 "
      self.mll300_400= "lepsMll_Edge < 400. && lepsMll_Edge >= 300 "
      self.mll400= "lepsMll_Edge >= 400 "
      self.ZmassExtended = "lepsMll_Edge >= 60 && lepsMll_Edge < 120"
      self.ZmassExtendedRsfof = "lepsMll_Edge >= 70 && lepsMll_Edge < 110"
      self.Zveto = "!(lepsMll_Edge >= 86 && lepsMll_Edge < 96)"
      self.baseline = self.AddList([self.nj2, self.mT2_80, self.dPhiJETMET,self.goodLepton])
      self.baselineNoTrigger = self.AddList([self.nj2, self.mT2_80, self.dPhiJETMET,self.goodLeptonNoTrigger])
      self.baselineNoMT2 = self.AddList([self.nj2,  self.dPhiJETMET,self.goodLepton])
      self.baselineNoMT2NoTrigger = self.AddList([self.nj2,  self.dPhiJETMET,self.goodLeptonNoTrigger])

      ########################################################################
      ######Edge regions #####################################################
      ########################################################################
      self.JETMETBaseline = self.AddList([self.METg150, self.nj2, self.mT2_80])
      self.JETMETBaselineNoMT2 = self.AddList([self.METg150, self.nj2])
      self.JETMETBaselineNoMT2SF_86_96 = self.AddList([self.METg150, self.nj2, self.SF, self.mll86_96])
      self.SignalRegionNoDPhi = self.AddList([self.goodLepton, self.JETMETBaseline]) 
      self.SignalRegion = self.AddList([self.goodLepton, self.JETMETBaseline, self.dPhiJETMET]) 
      self.SignalRegionNoMT2 = self.AddList([self.goodLepton, self.JETMETBaselineNoMT2, self.dPhiJETMET]) 
      self.SignalRegionNLL = self.AddList([self.goodLepton, self.JETMETBaseline, self.dPhiJETMET, self.NLL]) 

      self.RSFOFDirectControlRegion = self.AddList([self.METg100, self.METl150, self.njExact2, self.donot(self.mll70_110)])
      self.RSFOFDirectControlRegionOF = self.AddList([self.METg100, self.METl150, self.njExact2, self.donot(self.mll70_110), self.OF])
      self.RSFOFDirectControlRegionNoJet = self.AddList([self.goodLepton, self.METg100, self.METl150, self.donot(self.mll70_110)])
      self.RSFOFDirectControlRegionNoMll = self.AddList([self.goodLepton,  self.METg100, self.METl150, self.njExact2])
      self.RSFOFDirectControlRegionNoMET = self.AddList([self.goodLepton,  self.njExact2, self.donot(self.mll70_110)])
      self.RSFOFDirectSignalRegion = self.AddList([self.goodLepton, self.nj2, self.METg150, self.donot(self.mll70_110)])
      self.RSFOFDirectSignalRegionNoMET = self.AddList([self.goodLepton, self.nj2, self.donot(self.mll70_110)])
      self.RSFOFDirectSignalRegionNoJet = self.AddList([self.goodLepton, self.METg150, self.donot(self.mll70_110)])
 
      self.DYControlRegion = self.AddList([self.goodLepton, self.METl50, self.nj2, self.ZmassExtended])
      self.DYControlRegionNoJet = self.AddList([self.goodLepton, self.METl50, self.ZmassExtended])
      self.DYControlRegionNoMllNoNJet = self.AddList([self.goodLepton, self.METl50])
      self.DYControlRegionNoMll = self.AddList([self.goodLepton, self.METl50, self.nj2])
      self.DYControlRegionNoMET = self.AddList([self.goodLepton, self.nj2, self.ZmassExtended])
      self.DYControlRegionNoMllNoMET = self.AddList([self.goodLepton, self.nj2])

      self.Edge20Mll60   = '(lepsMll_Edge < 60)'                                              
      self.Edge60Mll86   = '(lepsMll_Edge > 60 ) && (lepsMll_Edge < 86 )'
      self.Edge96Mll150  = '(lepsMll_Edge > 96 ) && (lepsMll_Edge < 150)'
      self.Edge150Mll200 = '(lepsMll_Edge > 150) && (lepsMll_Edge < 200)'
      self.Edge200Mll300 = '(lepsMll_Edge > 200) && (lepsMll_Edge < 300)'
      self.Edge300Mll400 = '(lepsMll_Edge > 300) && (lepsMll_Edge < 400)'
      self.Edge400MllInf = '(lepsMll_Edge > 400)'
                                                                                            
      self.ttBarLike    = ' nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) < 21'
      self.NonttBarLike = ' nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) > 21'
      self.loNLL = self.ttBarLike
      self.hiNLL = self.NonttBarLike                                                        

      ########################################################################
      ######EWK signal regions ###############################################
      self.region3l  = self.AddList([self.goodLepton, self.threeTightLeptons, self.dPhiJETMET, self.METg60, self.bveto, self.nj2])
      self.region4l  = self.AddList([self.goodLepton, self.fourTightLeptons,  self.dPhiJETMET, self.mZ2g20, self.bveto, self.nj2])
      self.region4lInc  = self.AddList([self.goodLeptonInc, self.fourLooseLeptons, self.ZmassInc, self.mZ140, self.mZ4_120 ])
      self.regionttZ = self.AddList([self.goodLepton, self.threeTightLeptons, self.dPhiJETMET, self.nbj2, self.METg30, self.nj2])
      self.regionEdge3l  = self.AddList([self.goodLepton, self.threeTightLeptons, self.njless2, self.dPhiJETMET, self.METg60, self.bveto])  
      self.regionEdge4l  = self.AddList([self.goodLepton, self.fourTightLeptons,  self.njless2, self.dPhiJETMET, self.mZ2g20, self.bveto])
      self.regionEdgettZ = self.AddList([self.goodLepton, self.threeTightLeptons, self.nbjless2, self.njless2, self.dPhiJETMET, self.METg30])

      ########################################################################
      ######EWK signal regions ###############################################
      ########################################################################
      self.Baseline = self.AddList([self.nj2, self.METg100,self.dPhiJETMET,self.goodLepton])
      self.BaselineNoTrigger = self.AddList([self.nj2, self.METg100,self.dPhiJETMET,self.goodLeptonNoTrigger])
      self.EdgeBaseline = self.AddList( [self.baseline,self.METg150, self.mT2_80])
      self.ewinoWZ      = self.AddList([self.baselineNoTrigger, self.bveto, self.mjj110, self.ThirdLeptonVeto, self.Zmass])
      self.ewinoWZNoTrigger = self.AddList([self.baselineNoTrigger, self.bveto, self.mjj110, self.ThirdLeptonVeto, self.Zmass])
      self.ewinoWZExtMll  = self.AddList([self.baseline, self.bveto, self.mjj110, self.ThirdLeptonVeto])
      self.ewinoZH        = self.AddList([self.baselineNoMT2NoTrigger,self.nbExact2, self.mT2b200, self.mbb150, self.ThirdLeptonVeto, self.Zmass])
      self.ewinoZHExtMll  = self.AddList([self.baselineNoMT2,self.nbExact2, self.mT2b200, self.mbb150, self.ThirdLeptonVeto])
      self.strongOnZBVeto     = self.AddList([self.baseline,self.ThirdLeptonVeto,self.mT2_80,  self.bveto])
      self.strongOnZWithB     = self.AddList([self.baseline,self.ThirdLeptonVeto,self.mT2_100, self.nbj1])
      self.strongOnZBase      = self.AddList([self.baseline,self.Zmass, self.ThirdLeptonVeto, self.OR(self.strongOnZBVeto,self.strongOnZWithB)])
      ########################################################################
      ######Basic HT cuts#####################################################
      ########################################################################
      self.HT = "(htJet35j_Edge > 200) "
      self.triggerHT = "( HLT_BIT_HLT_PFHT200_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT250_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT300_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT350_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT400_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT475_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT600_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT650_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT800_v_Edge > 0.5  )" 
      self.numerator =   self.AddList([self.goodLeptonNoTrigger,  self.donot(self.JETMETBaselineNoMT2SF_86_96), self.donot(self.RSFOFDirectControlRegionOF), self.HT, self.triggerHT, self.trigger])
      self.denominator = self.AddList([self.goodLeptonNoTrigger,  self.donot(self.JETMETBaselineNoMT2SF_86_96), self.donot(self.RSFOFDirectControlRegionOF), self.HT, self.triggerHT])  

   def donot(self, cut):
     return '(!' + self.brackets(cut) + ')'

   def brackets(self, cut):
      if cut != '':
          return '('+cut+')'

   def AddList(self, cutlist):
      returncut = ''
      for cut in cutlist:
          if cut != '':
              returncut += cut
              if not cutlist.index(cut) == len(cutlist)-1:
                  returncut += ' && '
      return self.brackets(returncut)
  
   def Add(self, cut1, cut2):
      if cut1 == '':
          return cut2 
      if cut2 == '':
          return cut1
      return self.brackets(cut1 + " && " + cut2 )
  
   def OR(self, cut1, cut2):

      return self.brackets(cut1 + " || " + cut2 )

   def MaxRun(self, run):
      
      return self.brackets("run <= %d"%(run))

   def MinRun(self, run):
      
      return self.brackets("run >  %d"%(run))

      

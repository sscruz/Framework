

class CutManager:
   'This class serves as an on-demand cut server'

   def __init__(self):

      self.twoLeptons = "nPairLep_Edge > 0 && passesFilters_Edge > 0 "
      self.tightCharge = 'Lep1_tightCharge_Edge > 0 && Lep2_tightCharge_Edge > 0'
      self.trigMMc = "(HLT_mu17mu8_dz_Edge  > 0 || HLT_mu17mu8_Edge > 0 || HLT_mu17tkmu8_Edge > 0 || HLT_mu17tkmu8_dz_Edge > 0 || HLT_mu27tkmu8_Edge > 0 || HLT_mu30tkmu11_noniso_Edge > 0 )"
      self.trigMMNEW = "( HLT_BIT_HLT_Mu30_TkMu11_v ==1 ||HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v == 1 || HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v == 1 || HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v == 1 ||HLT_BIT_HLT_Mu27_TkMu8_v==1 )"
      self.trigEEc = "(HLT_el17el12_dz_Edge > 0 || HLT_el23el12_dz_Edge > 0 || HLT_doubleele33_noniso_Edge > 0)"
      self.trigEENEW = "(HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v == 1 || HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v==1 || HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v==1 || HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v==1)"
      self.trigEMNEW = "(HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v ==1 || HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v == 1 || HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v ==1 || HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v == 1 || HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v ==1 || HLT_BIT_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v ==1 )"
      self.trigEMc = "(HLT_mu8el17_Edge > 0 || HLT_mu8el23_Edge > 0 || HLT_mu17el12_Edge > 0 || HLT_mu23el12_Edge > 0 || HLT_mu30ele30_noniso_Edge > 0)"
      self.leptonPt = "Lep1_pt_Edge > 25 && Lep2_pt_Edge > 20."
      self.leptonDR = "lepsDR_Edge > 0.1"       
      self.ECALCrack = "abs(abs(Lep1_eta_Edge) - 1.5) > 0.1 && abs(abs(Lep2_eta_Edge) - 1.5) > 0.1"
      self.leptonsMll = "lepsMll_Edge > 20"

      self.goodLepton = "("+self.twoLeptons + "&&" + self.leptonPt + "&&" + self.leptonDR + "&&" + self.ECALCrack + "&&" + self.leptonsMll +  ")"

      self.ee = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -121)"
      self.mm = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -169)"
      self.OF = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -143)"
      self.SF = "(" + self.ee + " || " +  self.mm + ")"
      self.AF = "(" + self.SF + " || " +  self.OF + ")"
      self.nj2 = "(nJetSel_Edge >= 2)"
      self.nj1 = "(nJetSel_Edge >= 1)"
      self.nj0 = "(nJetSel_Edge >= 0)"
      self.nbj2 = "(nBJetMedium25_Edge >= 2)"
      self.nbj1 = "(nBJetMedium25_Edge >= 1)"
      self.nbj0 = "(nBJetMedium25_Edge >= 0)"
      self.MET50  = "(met_Edge < 50)"
      self.MET100 = "(met_Edge < 100)"
      self.METg100 = "(met_Edge > 100)"
      self.MET150 = "(met_Edge > 150)"
      self.MET200 = "(met_Edge > 200)"
      self.JetMETBaseline = "(met_Edge > 150 && nJetSel_Edge >= 2 && mt2_Edge > 80.)"
      self.dPhiJetMET = " (abs(j1MetDPhi_Edge)> 0.4)&& (abs(j2MetDPhi_Edge)> 0.4)"
      self.Zmass = "lepsMll_Edge > 86 && lepsMll_Edge < 96"
      self.mll20_86= "lepsMll_Edge < 86. && lepsMll_Edge > 20 "
      self.mll86_96= "lepsMll_Edge < 96. && lepsMll_Edge > 86 "
      self.mll96_150= "lepsMll_Edge < 150. && lepsMll_Edge > 96 "
      self.mll150_200= "lepsMll_Edge < 200. && lepsMll_Edge > 150 "
      self.mll200_300= "lepsMll_Edge < 300. && lepsMll_Edge > 200 "
      self.mll300_400= "lepsMll_Edge < 400. && lepsMll_Edge > 300 "
      self.mll400= "lepsMll_Edge > 400 "
      self.ZmassExtended = "lepsMll_Edge > 61 && lepsMll_Edge < 121"
      self.Zveto = "!(lepsMll_Edge > 86 && lepsMll_Edge < 96)"
      self.bveto = "(nBJetMedium25_Edge == 0)"
      self.binv = "!(nBJetMedium25_Edge == 0)"
      self.trigger = "((" + self.trigMMc + " && " + self.mm + ") || (" + self.trigEEc + " && " + self.ee + ") || (" + self.trigEMc + " && " + self.OF + "))"
      self.SignalRegionBaseLine          = self.AddList([self.goodLepton, self.JetMETBaseline, self.dPhiJetMET, self.trigger ]) 
      self.SignalRegionBaseLineNoTrigger = self.AddList([self.goodLepton, self.JetMETBaseline, self.dPhiJetMET]) 
      self.region3l = '('+self.dPhiJetMET +' && nLepTight_Edge == 3 && met_Edge > 60 && nBJetMedium25_Edge == 0)'
      self.region4l = '('+self.dPhiJetMET +' && nLepTight_Edge == 4)'
      self.regionttZ ='('+self.dPhiJetMET +' && nLepTight_Edge == 3 &&  nBJetMedium25_Edge == 2)'
      self.ewino =  '('+self.goodLepton+'&&'+self.dPhiJetMET+'&&'+self.SF+'&&'+ self.nj2 +'&&'+ self.Zmass +     ' && nLepLoose_Edge == 2 )'
      self.ewinoClosureReg =  '('+self.goodLepton+'&&'+self.dPhiJetMET+'&&'+ self.nj2 +'&&'+ self.Zmass +        ' && nLepLoose_Edge == 2  && mt2_Edge > 80. && nBJetMedium25_Edge == 0  )'
      self.ewinoClosureExt =  '('+self.goodLepton+'&&'+self.dPhiJetMET+'&&'+ self.nj2 +'&&'+ self.ZmassExtended +' && nLepLoose_Edge == 2  && mt2_Edge > 80.   )'
      self.ewinoCR ='('+self.goodLepton+'&&'+self.dPhiJetMET+'&&'+ self.nj2 +' && nLepLoose_Edge == 2 )'
      self.ewinoSR ='('+self.goodLepton+'&&'+self.dPhiJetMET+'&&'+ self.nj2 +'&&'+ self.Zmass +' && nBJetMedium25_Edge == 0 && mt2_Edge > 80. && nLepLoose_Edge == 2 && met_Edge > 150 )'

      self.nbExact2 = '(nBJetMedium25_Edge == 2)'
      self.mT2_80 = "(mt2_Edge > 80)"
      self.mT2_100 = "(mt2_Edge > 100)"
      self.bVeto  = "(nBJetMedium25_Edge == 0)"
      self.ThirdLeptonVeto = '(nLepLoose_Edge == 2)'
      self.JetMETPhi04 = "abs(j1MetDPhi_Edge) >  0.4 && abs(j2MetDPhi_Edge) > 0.4"
      # to do cuts #####################################
      print 'still to do this'
      self.mjj110  = 'dphiMjj_Edge < 110.'
      self.mbb150  = 'mbb_Edge < 150.'
      self.mT2b200 = 'mt2bb_Edge > 200.'
######################
      self.narrowZMass = '(lepsMll_Edge > 86 && lepsMll_Edge < 96)'
      self.trigger = "((" + self.trigMMc + " && " + self.mm + ") || (" + self.trigEEc + " && " + self.ee + ") || (" + self.trigEMc + " && " + self.OF + "))"
      self.Baseline = self.AddList([self.nj2,self.METg100,self.JetMETPhi04,self.goodLepton, self.trigger])
      self.BaselineNoTrigger = self.AddList([self.nj2,self.METg100,self.JetMETPhi04,self.goodLepton])
      self.EdgeBaseline = self.AddList( [self.MET150, self.mT2_80])
      self.ewinoCharNeu = self.AddList( [self.MET150, self.bVeto,  self.mjj110, self.narrowZMass, self.ThirdLeptonVeto]); 
      print 'has this a third lepton veto?'
      self.ewinoNeuNeu  = self.AddList( [self.nbExact2, self.narrowZMass, self.mT2b200, self.mbb150, self.ThirdLeptonVeto])
      self.strongOnZBVeto    = self.AddList( [ self.mT2_80  , self.bVeto ])
      self.strongOnZWithB    = self.AddList( [ self.mT2_100 , self.nbj1   ])
      self.strongOnZBase     = self.AddList( [self.narrowZMass, self.ThirdLeptonVeto, self.OR(self.strongOnZBVeto,self.strongOnZWithB)])

      self.lowmass = "lepsMll_Edge > 20 && lepsMll_Edge < 81"
      self.loMass= "lepsMll_Edge <  81."
      self.hiMass= "lepsMll_Edge > 101."
      self.highmass = "lepsMll_Edge > 101"
      self.SignalRegionBaseLine          = self.AddList([self.goodLepton, self.JetMETBaseline, self.trigger ]) 
      self.SignalRegionBaseLineNoTrigger = self.AddList([self.goodLepton, self.JetMETBaseline]) 
#      self.SignalRegionBaseLine = self.AddList([self.goodLepton, self.trigger, self.JetMETBaseline]) 
#      self.SignalRegionBaseLineNoTrigger = self.AddList([self.goodLepton, self.JetMETBaseline]) 

      self.region3l = '(nLepTight_Edge == 3 && met_Edge > 60 && nBJetMedium25_Edge == 0)'
      self.region4l = '(nLepTight_Edge == 4)'

      
      ##### Needed by RSFOF direct calculation ########
      self.RSFOFDirectControlRegion = "((met_Edge > 100 && met_Edge < 150 && nJetSel_Edge == 2) && ((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110))"
      self.RSFOFDirectControlRegionNoMll = "((met_Edge > 100 && met_Edge < 150 && nJetSel_Edge == 2) && (lepsMll_Edge>20)&& (abs(j1MetDPhi_Edge)> 0.4)&& (abs(j2MetDPhi_Edge)> 0.4))"
      self.RSFOFDirectControlRegionNoMET = "((nJetSel_Edge == 2) && (abs(j1MetDPhi_Edge)> 0.4)&& (abs(j2MetDPhi_Edge)> 0.4)&&((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110))"
      self.RSFOFDirectSignalRegion = "((met_Edge > 150 && nJetSel_Edge >= 2) &&(abs(j1MetDPhi_Edge)> 0.4)&& (abs(j2MetDPhi_Edge)> 0.4)&&((lepsMll_Edge > 20&& lepsMll_Edge < 70)||lepsMll_Edge >110))"
      self.RSFOFDirectSignalRegionNoMET = "((nJetSel_Edge >= 2) && (abs(j1MetDPhi_Edge)> 0.4)&& (abs(j2MetDPhi_Edge)> 0.4)&&((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110))"
      self.RSFOFDirectSignalRegionNoJet = "((nJetSel_Edge >= 2) && (abs(j1MetDPhi_Edge)> 0.4)&& (abs(j2MetDPhi_Edge)> 0.4)&&((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110) && met_Edge > 150.)"
      self.RSFOFDirectSignalRegionNLL   = "((nJetSel_Edge >= 2) && (abs(j1MetDPhi_Edge)> 0.4)&& (abs(j2MetDPhi_Edge)> 0.4)&& met_Edge > 150. && ((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110))"
      ##### Needed by rmue calculation################ 
      self.DYControlRegion =            "(nBJetMedium25_Edge == 0 &&"+self.MET50 +"&& nJetSel_Edge >= 2 && "+ self.Zmass + ")"
      self.DYControlRegionNoMll =       "(nBJetMedium25_Edge == 0 &&"+self.MET50 +"&& nJetSel_Edge >= 2 )"
      self.DYControlRegionNoNjetNoMll = "(nBJetMedium25_Edge == 0 &&"+self.MET50 +")"
      self.DYControlRegionNoJet =       "(nBJetMedium25_Edge == 0 &&"+self.MET50 +"&&"+ self.Zmass + ")"
      self.DYControlRegionNoMETNoMll =  "(nBJetMedium25_Edge == 0 && nJetSel_Edge >= 2 )"
      self.DYControlRegionNoMET =       "(nBJetMedium25_Edge == 0 && nJetSel_Edge >= 2 && "+ self.Zmass + ")"
      ##### Needed by RT calculation################ 
      self.HT = "(htJet35j_Edge > 200)"
      self.triggerHT = "(HLT_htall_Edge > 0 )" 
      self.triggerHTNEW = "( HLT_BIT_HLT_PFHT200_v ==1 ||HLT_BIT_HLT_PFHT250_v ==1 ||HLT_BIT_HLT_PFHT300_v ==1 ||HLT_BIT_HLT_PFHT350_v ==1 ||HLT_BIT_HLT_PFHT400_v ==1 ||HLT_BIT_HLT_PFHT475_v ==1 ||HLT_BIT_HLT_PFHT600_v ==1 ||HLT_BIT_HLT_PFHT650_v ==1 ||HLT_BIT_HLT_PFHT2800_v ==1 || HLT_BIT_HLT_PFHT300_PFMET110_v ==1 )" 
      self.numerator =   self.AddList([self.goodLepton, self.dPhiJetMET, self.donot(self.JetMETBaseline), self.donot(self.RSFOFDirectControlRegion), self.HT, self.triggerHT, self.trigger])
      self.denominator = self.AddList([self.goodLepton, self.dPhiJetMET, self.donot(self.JetMETBaseline), self.donot(self.RSFOFDirectControlRegion), self.HT, self.triggerHT])  

      ## for checks of the excess
      self.tightIso = 'max(Lep1_miniRelIso_Edge, Lep2_miniRelIso_Edge) < 0.05'

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

      

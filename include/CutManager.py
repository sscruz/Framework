

class CutManager:
   'This class serves as an on-demand cut server'

   def __init__(self):

      self.protection = '((run_Edge<=276811)  ||  (278820<=run_Edge && run_Edge<=279931))'

      ########################################################################
      ######Basic Lepton Cuts ################################################
      ########################################################################
      self.twoLeptons = "nPairLep_Edge > 0 && passesFilters_Edge > 0 "
      self.tightCharge = 'Lep1_tightCharge_Edge > 0 && Lep2_tightCharge_Edge > 0'
      self.leptonPt = "Lep1_pt_Edge > 25 && Lep2_pt_Edge > 20."
      self.dilepPt = "lepsZPt_Edge > 25"
      self.leptonDR = "lepsDR_Edge > 0.1"       
      self.ECALCrack = "abs(abs(Lep1_eta_Edge) - 1.5) > 0.1 && abs(abs(Lep2_eta_Edge) - 1.5) > 0.1"
      self.leptonsMll = "lepsMll_Edge > 20"
      self.goodLepton = "("+self.twoLeptons + "&&" + self.dilepPt + "&&" + self.leptonPt + "&&" + self.leptonDR + "&&" + self.ECALCrack + "&&" + self.leptonsMll +  ")"
      self.ee = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -121)"
      self.mm = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -169)"
      self.OF = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -143)"
      self.SF = "(" + self.ee + " || " +  self.mm + ")"
      self.AF = "(" + self.SF + " || " +  self.OF + ")"
      self.trigMMNEW = "( HLT_BIT_HLT_Mu30_TkMu11_v_Edge ==1 ||HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v_Edge == 1 || HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_Edge == 1 || HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v_Edge == 1 ||HLT_BIT_HLT_Mu27_TkMu8_v_Edge==1 )"
      self.trigEENEW = "(HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge == 1 || HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge==1 || HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v_Edge==1 || HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v_Edge==1)"
      self.trigEMNEW = "(HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge ==1 || HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge == 1 || HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v_Edge ==1 || HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v_Edge == 1 || HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v_Edge ==1 || HLT_BIT_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v_Edge ==1 || HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge == 1 || HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge == 1 || HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge == 1)"

      self.trigger = "((" + self.trigMMNEW + " && " + self.mm + ") || (" + self.trigEENEW + " && " + self.ee + ") || (" + self.trigEMNEW + " && " + self.OF + "))"
      self.ThirdLeptonVeto = '(nLepLoose_Edge == 2)'
      self.tightIso = 'max(Lep1_miniRelIso_Edge, Lep2_miniRelIso_Edge) < 0.05'
      self.threeTightLeptons = 'nLepTight_Edge == 3'
      self.fourTightLeptons = 'nLepTight_Edge == 4'
      self.twoLooseLeptons = 'nLepLoose_Edge == 2'
      ########################################################################
      ######Basic Jets and bJets Cuts###############################################
      ########################################################################
      self.nj2 = "(nJetSel_Edge >= 2)"
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
      
      ########################################################################
      ######Basic MET cuts##########################################################
      ########################################################################
      self.METl50  = "(met_Edge < 50)"
      self.METl100 = "(met_Edge < 100)"
      self.METl150 = "(met_Edge < 150)"
      self.METg30 = "(met_Edge > 30)"
      self.METg60 = "(met_Edge > 60)"
      self.METg80 = "(met_Edge >  80)"
      self.METg100 = "(met_Edge >= 100)"
      self.METg150 = "(met_Edge >= 150)"
      self.METg200 = "(met_Edge >= 200)"
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

      ########################################################################
      ######Basic Mass cuts###################################################
      ########################################################################
      self.wideZ = "lepsMll_Edge > 76 && lepsMll_Edge < 106"
      self.Zmass = "lepsMll_Edge >= 86 && lepsMll_Edge < 96"
      self.mll20_60= "lepsMll_Edge < 60. && lepsMll_Edge >= 20 "
      self.mll60_86= "lepsMll_Edge < 86. && lepsMll_Edge >= 60 "
      self.mll86_96= "lepsMll_Edge < 96. && lepsMll_Edge >= 86 "
      self.mll96_150= "lepsMll_Edge < 150. && lepsMll_Edge >= 96 "
      self.mll150_200= "lepsMll_Edge < 200. && lepsMll_Edge >= 150 "
      self.mll200_300= "lepsMll_Edge < 300. && lepsMll_Edge >= 200 "
      self.mll300_400= "lepsMll_Edge < 400. && lepsMll_Edge >= 300 "
      self.mll400= "lepsMll_Edge >= 400 "
      self.ZmassExtended = "lepsMll_Edge >= 61 && lepsMll_Edge < 121"
      self.Zveto = "!(lepsMll_Edge >= 86 && lepsMll_Edge < 96)"
      self.baseline = self.AddList([self.dilepPt, self.nj2, self.mT2_80, self.dPhiJETMET,self.goodLepton])
      self.baselineNoMT2 = self.AddList([self.dilepPt, self.nj2,  self.dPhiJETMET,self.goodLepton])
 

      ########################################################################
      ######Edge regions #####################################################
      ########################################################################
      self.JETMETBaseline = self.AddList([self.METg150, self.nj2, self.mT2_80])
      self.JETMETBaselineNoMT2 = self.AddList([self.METg150, self.nj2])
      self.SignalRegion = self.AddList([self.goodLepton, self.JETMETBaseline, self.dPhiJETMET]) 
      self.SignalRegionNoMT2 = self.AddList([self.goodLepton, self.JETMETBaselineNoMT2, self.dPhiJETMET]) 
      self.SignalRegionNLL = self.AddList([self.goodLepton, self.JETMETBaseline, self.dPhiJETMET, self.NLL]) 

      self.RSFOFDirectControlRegion = self.AddList([self.goodLepton, self.dPhiJETMET, self.METg100, self.METl150, self.njExact2, self.donot(self.Zmass)])
      self.RSFOFDirectControlRegionNoMll = self.AddList([self.goodLepton, self.dPhiJETMET, self.METg100, self.METl150, self.njExact2])
      self.RSFOFDirectControlRegionNoMET = self.AddList([self.goodLepton, self.dPhiJETMET, self.njExact2, self.donot(self.Zmass)])
      self.RSFOFDirectSignalRegion = "((met_Edge > 150 && nJetSel_Edge >= 2) &&(abs(j1MetDPhi_Edge)> 0.4)&& (abs(j2MetDPhi_Edge)> 0.4)&&((lepsMll_Edge > 20&& lepsMll_Edge < 70)||lepsMll_Edge >110))"
      self.RSFOFDirectSignalRegionNoMET = "((nJetSel_Edge >= 2) && (abs(j1MetDPhi_Edge)> 0.4)&& (abs(j2MetDPhi_Edge)> 0.4)&&((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110))"
      self.RSFOFDirectSignalRegionNoJet = "((nJetSel_Edge >= 2) && (abs(j1MetDPhi_Edge)> 0.4)&& (abs(j2MetDPhi_Edge)> 0.4)&&((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110) && met_Edge > 150.)"
      self.RSFOFDirectSignalRegionNLL   = "((nJetSel_Edge >= 2) && (abs(j1MetDPhi_Edge)> 0.4)&& (abs(j2MetDPhi_Edge)> 0.4)&& met_Edge > 150. && ((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110))"
 
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
      ########################################################################
      self.region3l = self.AddList([self.threeTightLeptons, self.dPhiJETMET, self.METg60, self.bveto])
      self.region4l = self.AddList([self.fourTightLeptons, self.dPhiJETMET])
      self.regionttZ = self.AddList([self.threeTightLeptons, self.dPhiJETMET, self.nbj2, self.METg30 ])

      self.ewinoExtMll  = self.AddList([self.baseline,self.goodLepton,self.dPhiJETMET, self.nj2, self.twoLooseLeptons, self.mT2_80, self.bveto, self.mjj110, self.ThirdLeptonVeto])
      self.ewinoSR      = self.AddList([self.baseline,self.goodLepton,self.dPhiJETMET, self.nj2, self.twoLooseLeptons, self.mT2_80, self.bveto, self.mjj110, self.ThirdLeptonVeto, self.Zmass])

      ########################################################################
      ######EWK signal regions ###############################################
      ########################################################################
      ##### Warning: these regions I didn't touch ############################
      self.Baseline = self.AddList([self.nj2, self.METg100,self.dPhiJETMET,self.goodLepton])
      self.BaselineNoTrigger = self.AddList([self.nj2,self.METg100,self.dPhiJETMET,self.goodLepton])
      self.EdgeBaseline = self.AddList( [self.baseline,self.METg150, self.mT2_80])
      self.ewinoNeuNeu        = self.AddList([self.baseline,self.nbExact2, self.mT2b200, self.mbb150, self.ThirdLeptonVeto, self.Zmass])
      self.ewinoNeuNeuExtMll  = self.AddList([self.baseline,self.nbExact2, self.mT2b200, self.mbb150, self.ThirdLeptonVeto])
      self.strongOnZBVeto     = self.AddList([self.baseline,self.ThirdLeptonVeto,self.mT2_80, self.bveto + "&& ( ((htJet35j_Edge >= 500)&& ( nJetSel_Edge < 6))|| (nJetSel_Edge >=6))  " ])
      self.strongOnZWithB     = self.AddList([self.baseline,self.ThirdLeptonVeto,self.mT2_100, self.nbj1  + "&& ( ((htJet35j_Edge >= 200)&& ( nJetSel_Edge < 6))|| (nJetSel_Edge >=6 ))  " ])
      self.strongOnZBase      = self.AddList([self.baseline,self.Zmass, self.ThirdLeptonVeto, self.OR(self.strongOnZBVeto,self.strongOnZWithB)])

      ########################################################################
      ######Basic HT cuts#####################################################
      ########################################################################
      self.HT = "(htJet35j_Edge > 200) "
      self.triggerHT = "( HLT_BIT_HLT_PFHT200_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT250_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT300_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT350_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT400_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT475_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT600_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT650_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT800_v_Edge > 0.5 || HLT_BIT_HLT_PFHT300_PFMET110_v_Edge > 0.5 )" 
      self.numerator =   self.AddList([self.goodLepton, self.donot(self.JETMETBaseline), self.donot(self.RSFOFDirectControlRegion), self.HT, self.triggerHT, self.trigger])
      self.denominator = self.AddList([self.goodLepton, self.donot(self.JETMETBaseline), self.donot(self.RSFOFDirectControlRegion), self.HT, self.triggerHT])  



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

      

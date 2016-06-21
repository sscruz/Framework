

class CutManager:
   'This class serves as an on-demand cut server'

   def __init__(self):

      self.twoLeptons = "nPairLep_Edge > 0 && hbheFilterIso_Edge > 0 && hbheFilterNew25ns_Edge > 0 && Flag_eeBadScFilter_Edge > 0 "
      self.trigMMc = "(HLT_mu17mu8_dz  > 0 || HLT_mu30tkmu11_noniso > 0)"
      self.trigEEc = "(HLT_el17el12_dz > 0 || HLT_el23el12_dz > 0 || HLT_mu30el30_noniso > 0)"
      self.trigEMc = "(HLT_mu8el17 > 0 || HLT_mu8el23 > 0)"
      self.leptonPt = "Lep1_pt_Edge > 25. && Lep2_pt_Edge > 20."
      self.leptonDR = "lepsDR_Edge > 0.1"       
      self.ECALCrack = "abs(abs(Lep1_eta_Edge) - 1.5) > 0.1 && abs(abs(Lep2_eta_Edge) - 1.5) > 0.1"
      self.leptonsMll = "lepsMll_Edge > 20"
      self.goodLepton = self.twoLeptons + "&&" + self.leptonPt + "&&" + self.leptonDR + "&&" + self.ECALCrack + "&&" + self.leptonsMll
      self.ee = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -121)"
      self.mm = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -169)"
      self.OF = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -143)"
      self.SF = "(" + self.ee + " || " +  self.mm + ")"
      self.AF = "(" + self.SF + " || " +  self.OF + ")"
      self.nj2 = "(nJetSel_Edge >= 2)"
      self.nj1 = "(nJetSel_Edge >= 0)"
      self.nj0 = "(nJetSel_Edge >= 0)"
      self.nbj2 = "(nbJetSel_Edge >= 2)"
      self.nbj1 = "(nbJetSel_Edge >= 0)"
      self.nbj0 = "(nbJetSel_Edge >= 0)"
      self.MET100 = "(met_Edge > 100)"
      self.MET150 = "(met_Edge > 150)"
      self.MET200 = "(met_Edge > 200)"
      self.JetMETBaseline = "(met_Edge > 150 && nJetSel_Edge >= 2)"
      self.lowmass = "lepsMll_Edge > 20 && lepsMll_Edge < 81"
      self.Zmass = "lepsMll_Edge > 81 && lepsMll_Edge < 101"
      self.highmass = "lepsMll_Edge > 101"
      self.trigger = "((" + self.trigMMc + " && " + self.mm + ") || (" + self.trigEEc + " && " + self.ee + ") || (" + self.trigEMc + " && " + self.OF + "))"
      self.SignalRegionBaseLine = self.AddList([self.goodLepton, self.trigger, self.JetMETBaseline]) 
      
      ##### Needed by RSFOF direct calculation ########
      self.RSFOFDirectControlRegion = "((met_Edge > 100 && met_Edge < 150 && nJetSel_Edge == 2) && ((lepsMll_Edge > 20 && lepsMll_Edge < 81) || lepsMll_Edge > 101))"
      self.RSFOFDirectControlRegionRestrictive = "((met_Edge > 100 && met_Edge < 150 && nJetSel_Edge == 2) && ((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110))"
      self.RSFOFDirectSignalRegionRestrictive = "((met_Edge > 150 && nJetSel_Edge >= 2) && ((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110))"
      self.RSFOFDirectSignalRegionRestrictiveNoMET = "((nJetSel_Edge >= 2) && ((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110))"
      self.RSFOFDirectSignalRegionRestrictiveNoJet = "((nJetSel_Edge >= 2) && ((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110))"
      self.RSFOFDirectControlRegionNoMll = "((met_Edge > 100 && met_Edge < 150 && nJetSel_Edge == 2) && (lepsMll_Edge))"
      self.RSFOFDirectControlRegionNoMET = "((nJetSel_Edge == 2) && ((lepsMll_Edge > 20 && lepsMll_Edge < 81) || lepsMll_Edge > 101))"
      ##### Needed by rmue calculation################ 
      self.DYControlRegion = "(met_Edge < 50 && nJetSel_Edge >= 2 && lepsMll_Edge > 81 && lepsMll_Edge < 101)"
      self.DYControlRegionNoMll = "(met_Edge < 50 && nJetSel_Edge >= 2)"
      self.DYControlRegionNoMET = "(nJetSel_Edge >= 2 && lepsMll_Edge > 81 && lepsMll_Edge < 101)"
      self.DYControlRegionNoJet = "(met_Edge < 50 && lepsMll_Edge > 81 && lepsMll_Edge < 101)"
      ##### Needed by RT calculation################ 
      self.HT = "(htJet35j_Edge > 200)"
      self.triggerHT = "(HLT_htall > 0 || HLT_htmet > 0 || HLT_atall > 0)"
      self.numerator = self.AddList([self.goodLepton, self.donot(self.JetMETBaseline), self.donot(self.RSFOFDirectControlRegion), self.HT, self.triggerHT])


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
 
      

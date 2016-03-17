

class CutManager:
   'This class serves as an on-demand cut server'

   def __init__(self):

      self.twoLeptons = "t.nPairLep_Edge > 0 && hbheFilterIso > 0 && hbheFilterNew25ns > 0 && Flag_eeBadScFilter > 0 "
      self.trigMMc = "(HLT_DoubleMu > 0 || HLT_mu27tkmu8 > 0)"
      self.trigEEc = "(HLT_el17el12_dz > 0 || HLT_ele33ele33 > 0)"
      self.trigEMc = "(HLT_mu8el17 > 0 || HLT_mu17el12 > 0 || HLT_mu30ele30 > 0)"
      self.leptonPt = "t.Lep1_pt_Edge > 20. && t.Lep2_pt_Edge > 20."
      self.leptonDR = "t.lepsDR_Edge > 0.3"       
      self.ECALCrack = "abs(abs(Lep1_eta_Edge) - 1.5) > 0.1 && abs(abs(Lep2_eta_Edge) - 1.5) > 0.1"
      self.leptonsMll = "t.lepsMll_Edge > 20"
      self.goodLepton = self.twoLeptons + "&&" + self.leptonPt + "&&" + self.leptonDR + "&&" + self.ECALCrack + "&&" + self.leptonsMll
      self.ee = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -121)"
      self.mm = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -169)"
      self.OF = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -143)"
      self.SF = "(" + self.ee + " || " +  self.mm + ")"
      self.AF = "(" + self.SF + " || " +  self.OF + ")"
      self.nj2 = "(t.nJetSel_Edge >= 2)"
      self.nj0 = "(t.nJetSel_Edge >= 0)"
      self.METJetsSignalRegion = "((met_pt > 150 && t.nJetSel_Edge > 1) || (met_pt > 100 && t.nJetSel_Edge > 2))"
      self.METJetsSignalRegion2J = "(((met_pt > 150 && t.nJetSel_Edge > 1) || (met_pt > 100 && t.nJetSel_Edge > 2)) && t.nJetSel_Edge == 2)"
      self.METJetsSignalRegion3J = "(((met_pt > 150 && t.nJetSel_Edge > 1) || (met_pt > 100 && t.nJetSel_Edge > 2)) && t.nJetSel_Edge >  2)"
      self.METJetsSignalRegionMET100 = "(met_pt > 100 && t.nJetSel_Edge > 1)"
      self.METJetsSignalRegionMET150 = "(met_pt > 150 && t.nJetSel_Edge > 1)"
      self.METJetsControlRegion = "(met_pt > 100 && met_pt < 150 && t.nJetSel_Edge == 2)"
      self.RSFOFControlAlternative = "(met_pt > 50 && met_pt < 150 && t.nJetSel_Edge == 2 && t.nBJetLoose35_Edge >=1 )"
      self.DYControlRegion = "(met_pt < 50 && t.nJetSel_Edge >= 2)"
      self.blinded = "!"+self.METJetsSignalRegion
      self.DYmet = "(met_pt < 50)"
      self.DYmass = "t.lepsMll_Edge > 60 && t.lepsMll_Edge < 120"
      self.ZmassVeto = "(t.lepsMll_Edge < 81 || t.lepsMll_Edge > 101)"
      self.DYmassVeto = "(t.lepsMll_Edge < 70 || t.lepsMll_Edge > 110)"
      self.lowmass = "t.lepsMll_Edge > 20 && t.lepsMll_Edge < 70"
      self.Zmass = "t.lepsMll_Edge > 81 && t.lepsMll_Edge < 101"
      self.highmass = "t.lepsMll_Edge > 120"
      self.central = "(abs(t.Lep1_eta_Edge)<1.4 && abs(t.Lep2_eta_Edge)<1.4)"
      self.forward = "(abs(t.Lep1_eta_Edge)>1.4 || abs(t.Lep2_eta_Edge)>1.4)"
      self.trigger = "((" + self.trigMMc + " && " + self.mm + ") || (" + self.trigEEc + " && " + self.ee + ") || (" + self.trigEMc + " && " + self.OF + "))"
      self.HT = "(t.htJet35j_Edge > 200)"
      self.triggerHT = "(HLT_pfht200 > 0 || HLT_pfht250 > 0 || HLT_pfht300 > 0 || HLT_pfht350 > 0 || HLT_pfht400>0 || HLT_pfht475>0 || HLT_pfht600>0 || HLT_pfht800>0)"
      self.triggerCalculation = self.AddList([self.triggerHT, self.HT, self.goodLepton, self.donot(self.METJetsSignalRegion), self.donot(self.AddList([self.METJetsControlRegion, self.OR(self.lowmass, self.highmass)]))])


   def donot(self, cut):
     return '(!' + self.brackets(cut) + ')'

   def brackets(self, cut):
      return '('+cut+')'

   def AddList(self, cutlist):
      returncut = ''
      for cut in cutlist:
          returncut += cut
          if not cutlist.index(cut) == len(cutlist)-1:
            returncut += ' && '
      return self.brackets(returncut)
  
   def Add(self, cut1, cut2):

      return self.brackets(cut1 + " && " + cut2 )
  
   def OR(self, cut1, cut2):

      return self.brackets(cut1 + " || " + cut2 )


   def MaxRun(self, run):
      
      return self.brackets("run <= %d"%(run))

   def MinRun(self, run):
      
      return self.brackets("run >  %d"%(run))

   def Central(self):
      
      return self.brackets(self.central)

   def trigMM(self):

      return self.brackets(self.trigMMc)
 
   def trigEE(self):

      return self.brackets(self.trigEEc)
 
   def trigEM(self):

      return self.brackets(self.trigEMc)
 
   def Forward(self):
      
      return self.brackets(self.forward)

   def DYMass(self):
 
      return self.brackets(self.DYmass)
 
   def GoodLeptonNoTriggerSF(self):

      return self.brackets(self.goodLepton + " && " + self.SF)

   def GoodLeptonNoTriggerOF(self):

      return self.brackets(self.goodLepton + " && " + self.OF) 

   def GoodLeptonNoTriggeree(self):

      return self.brackets(self.goodLepton + " && " + self.ee)

   def GoodLeptonNoTriggermm(self):

      return self.brackets(self.goodLepton + " && " + self.mm) 

   def GoodLepton(self):

      return self.brackets(self.goodLepton + " && " + self.trigger)

   def GoodLeptonAF(self):

      return self.brackets(self.goodLepton + " && " + self.AF + " && " + self.trigger)

   def GoodLeptonSFNoTrigger(self):

      return self.brackets(self.goodLepton + " && " + self.SF )

   def GoodLeptonSF(self):

      return self.brackets(self.goodLepton + " && " + self.SF + " && " + self.trigger)

   def GoodLeptonOF(self):

      return self.brackets(self.goodLepton + " && " + self.OF + " && " + self.trigger)

   def GoodLeptonee(self):

      return self.brackets(self.goodLepton + " && " + self.ee + " && " + self.trigger)

   def GoodLeptonmm(self):

      return self.brackets(self.goodLepton + " && " + self.mm + " && " + self.trigger)

   def SignalNoMassLeptonSF(self):

      return self.brackets(self.GoodLeptonSF() + " && " + self.METJetsSignalRegion)
   
   def SignalNoMassLeptonOF(self):

      return self.brackets(self.GoodLeptonOF() + " && " + self.METJetsSignalRegion)
   
   def SignalNoMassLeptonee(self):

      return self.brackets(self.GoodLeptonee() + " && " + self.METJetsSignalRegion)

   def SignalNoMassLeptonmm(self):

      return self.brackets(self.GoodLeptonmm() + " && " + self.METJetsSignalRegion)

   def ControlNoMassLeptonSF(self):

      return self.brackets(self.GoodLeptonSF() + " && " + self.METJetsControlRegion  + " && " + self.trigger)
   
   def ControlNoMassLeptonOF(self):

      return self.brackets(self.GoodLeptonOF() + " && " + self.METJetsControlRegion  + " && " + self.trigger)
   
   def ControlNoMassLeptonee(self):

      return self.brackets(self.GoodLeptonee() + " && " + self.METJetsControlRegion)

   def ControlNoMassLeptonmm(self):

      return self.brackets(self.GoodLeptonmm() + " && " + self.METJetsControlRegion)

   def DYControlNoMassLeptonSF(self):

      return self.brackets(self.GoodLeptonSF() + " && " + self.DYControlRegion  + " && " + self.trigger)
   
   def DYControlNoMassLeptonOF(self):

      return self.brackets(self.GoodLeptonOF() + " && " + self.DYControlRegion  + " && " + self.trigger)
   
   def DYControlNoMassLeptonee(self):

      return self.brackets(self.GoodLeptonee() + " && " + self.DYControlRegion)

   def DYControlNoMassLeptonmm(self):

      return self.brackets(self.GoodLeptonmm() + " && " + self.DYControlRegion)

   def SignalLowMassSF(self):

      return self.brackets(self.SignalNoMassLeptonSF() + " && " + self.lowmass)
   
   def SignalLowMassOF(self):

      return self.brackets(self.SignalNoMassLeptonOF() + " && " + self.lowmass)
   
   def SignalLowMassee(self):

      return self.brackets(self.SignalNoMassLeptonee() + " && " + self.lowmass)
   
   def SignalLowMassmm(self):

      return self.brackets(self.SignalNoMassLeptonmm() + " && " + self.lowmass)

   def SignalZMassSF(self):

      return self.brackets(self.SignalNoMassLeptonSF() + " && " + self.Zmass)
   
   def SignalZMassOF(self):

      return self.brackets(self.SignalNoMassLeptonOF() + " && " + self.Zmass)
   
   def SignalZMassee(self):

      return self.brackets(self.SignalNoMassLeptonee() + " && " + self.Zmass)
   
   def SignalZMassmm(self):

      return self.brackets(self.SignalNoMassLeptonmm() + " && " + self.Zmass)

   def SignalHighMassSF(self):

      return self.brackets(self.SignalNoMassLeptonSF() + " && " + self.highmass)
   
   def SignalHighMassOF(self):

      return self.brackets(self.SignalNoMassLeptonOF() + " && " + self.highmass)
   
   def SignalHighMassee(self):

      return self.brackets(self.SignalNoMassLeptonee() + " && " + self.highmass)
   
   def SignalHighMassmm(self):

      return self.brackets(self.SignalNoMassLeptonmm() + " && " + self.highmass)
 
   def ControlLowMassSF(self):

      return self.brackets(self.ControlNoMassLeptonSF() + " && " + self.lowmass)
   
   def ControlLowMassOF(self):

      return self.brackets(self.ControlNoMassLeptonOF() + " && " + self.lowmass)
   
   def ControlLowMassee(self):

      return self.brackets(self.ControlNoMassLeptonee() + " && " + self.lowmass)
   
   def ControlLowMassmm(self):

      return self.brackets(self.ControlNoMassLeptonmm() + " && " + self.lowmass)

   def ControlZMassSF(self):

      return self.brackets(self.ControlZMassLeptonSF() + " && " + self.Zmass)
   
   def ControlZMassOF(self):

      return self.brackets(self.ControlZMassLeptonOF() + " && " + self.Zmass)
   
   def ControlZMassee(self):

      return self.brackets(self.ControlZMassLeptonee() + " && " + self.Zmass)
   
   def ControlZMassmm(self):

      return self.brackets(self.ControlNoMassLeptonmm() + " && " + self.Zmass)

   def ControlHighMassSF(self):

      return self.brackets(self.ControlNoMassLeptonSF() + " && " + self.highmass)
   
   def ControlHighMassOF(self):

      return self.brackets(self.ControlNoMassLeptonOF() + " && " + self.highmass)
   
   def ControlHighMassee(self):

      return self.brackets(self.ControlNoMassLeptonee() + " && " + self.highmass)
   
   def Control2JetsSF(self):

      return self.brackets(self.nj2 + " && " + self.GoodLeptonSF()  + " && " + self.trigger)
   
   def Control2JetsOF(self):

      return self.brackets(self.nj2 + " && " + self.GoodLeptonOF()  + " && " + self.trigger)
   
   def Control2JetsAF(self):

      return self.brackets(self.nj2 + " && " + self.GoodLeptonAF()  + " && " + self.trigger)
   
   def RSFOFControlRegion(self):

      return self.brackets(self.RSFOFControlAlternative + " && " + self.trigger)
   
   def Control2Jetsee(self):

      return self.brackets(self.nj2 + " && " + self.GoodLeptonee())
   
   def Control2Jetsmm(self):

	return self.brackets(self.nj2 + " && " + self.GoodLeptonmm())
      
   def Control2JetsMETSF(self):

      return self.brackets(self.nj2 + " && " + self.DYmet + " && " + self.GoodLeptonSF())
   
   def Control2JetsMETOF(self):

      return self.brackets(self.nj2 + " && " + self.DYmet + " && " + self.GoodLeptonOF())
   
   def Control2JetsMETee(self):

      return self.brackets(self.nj2 + " && " + self.DYmet + " && " + self.GoodLeptonee())
   
   def Control2JetsMETmm(self):

	return self.brackets(self.nj2 + " && " + self.DYmet + " && " + self.GoodLeptonmm())
      

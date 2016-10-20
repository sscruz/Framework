import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TH1F, TGraph
import math, sys, optparse, copy, re, array, itertools, os


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables
import cPickle as pickle


class  darkmatterSys:
    '''lightweight version of the dark matter class, in order to implement all the systematics studies'''

    def __init__(self,name, sys):
        global cuts, lumi
        cuts = CutManager.CutManager()
        lumi = 12.9
        self.sys = sys
        self.name = name+'sys'
        self.setSRsCuts()
        self.setTrees()

    def setTrees(self):
        mcDatasetsHT  = ['TT_pow_ext34','DYJetsToLL_M50','DYJetsToLL_M10to50',
                         'WWTo2L2Nu', 'WZTo2L2Q', 'TTZ_LO','TTW_LO','TBar_tWch','ZZ','T_tWch']
        mcDatasetsNoTT  = ['DYJetsToLL_M50','WWTo2L2Nu', 'WZTo2L2Q', 'TTZ_LO','TTW_LO','TBar_tWch','ZZ','T_tWch']
        mcDatasetsNoTTZ = ['DYJetsToLL_M50_HT100to200_ext', 'DYJetsToLL_M50_HT200to400_ext',
                           'DYJetsToLL_M50_HT400to600_ext','DYJetsToLL_M50_HT600toInf_ext',
                           'WWTo2L2Nu', 'WZTo2L2Q', 'TTW_LO','TBar_tWch','ZZ','T_tWch']
        ttDatasets       = ['TT_pow_ext34']
        daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2_runs_273150_275376',
                      'DoubleEG_Run2016B-PromptReco-v2_runs_273150_275376',
                      'MuonEG_Run2016B-PromptReco-v2_runs_273150_275376',
                      'DoubleMuon_Run2016C-PromptReco-v2_runs_275420_276283',
                      'DoubleEG_Run2016C-PromptReco-v2_runs_275420_276283',
                      'MuonEG_Run2016C-PromptReco-v2_runs_275420_276283',
                      'DoubleMuon_Run2016D-PromptReco-v2_runs_276315_276811',
                      'DoubleEG_Run2016D-PromptReco-v2_runs_276315_276811',
                      'MuonEG_Run2016D-PromptReco-v2_runs_276315_276811']
        ttZDatasets = ['TTZ_LO']
        siDatasets   = ['TTDM_scalar_Mchi-1_Mphi-10']
        siDatasets300= ['TTDM_scalar_Mchi-1_Mphi-300']
        self.treeMCHT   = Sample.Tree(helper.selectSamples('samples.dat', mcDatasetsHT,   'MC'), 'MC', 0, isScan = 0, extraWeight=self.extraWeightsForSys[self.sys])
        self.treeTT     = Sample.Tree(helper.selectSamples('samples.dat', ttDatasets,   'MC'), 'MC', 0, isScan = 0, extraWeight=self.extraWeightsForSys[self.sys])
        self.treeDA   = Sample.Tree(helper.selectSamples('samples.dat', daDatasets,   'DA'), 'DA', 1, isScan = 0)
        self.treeSI   = Sample.Tree(helper.selectSamples('samples.dat', siDatasets,   'MC'), 'MC', 1, isScan = 0, extraWeight=self.extraWeightsForSys[self.sys])
        self.treeSI300= Sample.Tree(helper.selectSamples('samples.dat', siDatasets300,'MC'), 'MC', 1, isScan = 0, extraWeight=self.extraWeightsForSys[self.sys])
# TreeSI for data driven studies. choose treeDA when studies are done (choose treeMC for closure tests)
        self.treeDataDriven = self.treeDA
        self.treeMCNoTT     = Sample.Tree(helper.selectSamples('samples.dat', mcDatasetsNoTT,   'MC'),
                                          'MC', 0, isScan = 0, extraWeight=self.extraWeightsForSys[self.sys])
        self.treeTTZ        = Sample.Tree(helper.selectSamples('samples.dat', ttZDatasets,   'MC'),
                                          'MC', 0, isScan = 0, extraWeight=self.extraWeightsForSys[self.sys])
        self.treeMCNoTTZNoTT= Sample.Tree(helper.selectSamples('samples.dat', mcDatasetsNoTTZ,   'MC'),
                                          'MC', 0, isScan = 0, extraWeight=self.extraWeightsForSys[self.sys])
    def setSRsCuts(self):
        theCuts = []
        self.baselineCuts = [cuts.goodLepton, 'nBJetMedium35_Edge > 0','mt2_Edge > 120',
                             'TMath::Abs(phi_ptboost) < 1.','nLepLoose_Edge < 3', cuts.ZVeto]
        self.ttControCuts = [cuts.goodLepton, 'nBJetMedium35_Edge > 0', 'TMath::Abs(phi_ptboost) < 1.',
                             'nLepLoose_Edge < 3','mt2_Edge > 60','mt2_Edge < 100','nJetSel_Edge>1',cuts.ZVeto]
        self.ttzControCuts = [cuts.goodLepton, 'nBJetMedium35_Edge > 0','mt2_Edge > 90', 'nLepTight_Edge > 2',
                              'nJetSel_Edge>1'] 
        self.baseCuts      = [ cuts.goodLepton, cuts.ZVeto]

        self.replaceCutsForSys = {'':      [],
                                  'jecUp': [['nBJetMedium35_Edge','nBJetMedium35_jecUp_Edge'],
                                            ['mt2_Edge','mt2_jecUp_Edge'],
                                            ['phi_ptboost','phi_ptboost_jecUp'],
                                            ['nJetSel_Edge','nJetSel_jecUp_Edge']],
                                  'jecDn': [['nBJetMedium35_Edge','nBJetMedium35_jecDn_Edge'],
                                            ['mt2_Edge','mt2_jecDn_Edge'],
                                            ['phi_ptboost','phi_ptboost_jecDn'],
                                            ['nJetSel_Edge','nJetSel_jecDn_Edge']],
                                  'bHeUp': [],
                                  'bHeDn': [],
                                  'bLiUp': [],
                                  'bLiDn': [],
                                  'ElUp' : [],
                                  'ElDn' : [],
                                  'MuUp' : [],
                                  'MuDn' : []}
        

        self.extraWeightsForSys = {'': '1',
                                   'jecUp': '1',
                                   'jecDn': '1',
                                   'bHeUp': 'weight_btagsf_heavy_UP_Edge/weight_btagsf_Edge',
                                   'bHeDn': 'weight_btagsf_heavy_DN_Edge/weight_btagsf_Edge',
                                   'bLiUp': 'weight_btagsf_light_UP_Edge/weight_btagsf_Edge',
                                   'bLiDn': 'weight_btagsf_light_DN_Edge/weight_btagsf_Edge',
                                   'ElUp' : 'sfs_ElUp / sfs',
                                   'ElDn' : 'sfs_ElDn / sfs',
                                   'MuUp' : 'sfs_MuUp / sfs',
                                   'MuDn' : 'sfs_MuDn / sfs'}


        self.baselineCuts = self.getTheCuts(self.baselineCuts)
        self.ttControCuts = self.getTheCuts(self.ttControCuts)
        self.ttzControCuts= self.getTheCuts(self.ttzControCuts)

        self.baselineCuts = self.makeReplacement(self.baselineCuts)
        self.ttControCuts = self.makeReplacement(self.ttControCuts)
        self.ttzControCuts= self.makeReplacement(self.ttzControCuts)

        # print 'jes up' 
        # self.baselineCutsjecUp = [cuts.goodLepton, 'nBJetMedium35_jecUp_Edge > 0','mt2_jecUp_Edge > 120',
        #                      'TMath::Abs(phi_ptboost_jecUp) < 1.','nLepLoose_Edge < 3']
        # self.ttControCutsjecUp = [cuts.goodLepton, 'nBJetMedium35_jecUp_Edge > 0', 'TMath::Abs(phi_ptboost_jecUp) < 1.',
        #                      'nLepLoose_Edge < 3','mt2_jecUp_Edge > 60','mt2_Edge < 100','nJetSel_jecUp_Edge>1']
        # self.ttzControCutsjecUp = [cuts.goodLepton, 'nBJetMedium35_jecUp_Edge > 0','mt2_jecUp_Edge > 90', 'nLepTight_Edge > 2',
        #                       'nJetSel_jecUp_Edge>1'] 


        # print 'jec dn'
        # self.baselineCutsjecDn = [cuts.goodLepton, 'nBJetMedium35_jecDn_Edge > 0','mt2_jecDn_Edge > 120',
        #                      'TMath::Abs(phi_ptboost_jecDn) < 1.','nLepLoose_Edge < 3']
        # self.ttControCutsjecDn = [cuts.goodLepton, 'nBJetMedium35_jecDn_Edge > 0', 'TMath::Abs(phi_ptboost_jecDn) < 1.',
        #                      'nLepLoose_Edge < 3','mt2_jecDn_Edge > 60','mt2_Edge < 100','nJetSel_jecDn_Edge>1']
        # self.ttzControCutsjecDn = [cuts.goodLepton, 'nBJetMedium35_jecDn_Edge > 0','mt2_jecDn_Edge > 90', 'nLepTight_Edge > 2',
        #                       'nJetSel_jecDn_Edge>1'] 



#         for sys in ['','jecUp','jecDn']:
# #            setattr(self,'baseCuts'+sys     , self.getTheCuts(getattr(self,'baseCuts'+sys)))
#             setattr(self,'baselineCuts'+sys , self.getTheCuts(getattr(self,'baselineCuts'+sys)))
#             setattr(self,'ttControCuts'+sys , self.getTheCuts(getattr(self,'ttControCuts'+sys)))
#             setattr(self,'ttzControCuts'+sys, self.getTheCuts(getattr(self,'ttzControCuts'+sys),False))


    def makeReplacement(self, cut):
        for rpl in self.replaceCutsForSys[self.sys]:
            cut = cut.replace(rpl[0],rpl[1])
        return cut

    def getTheCuts(self,theCuts):
        finalcuts = ''
        finalcuts = cuts.AddList(theCuts)
        return finalcuts

    def GetDataMCFactorTT(self):
        # estimation = observed in CR x Transfer factor
        _den = self.treeTT.getYields(lumi,'0.5',0,1,"("+self.ttControCuts+")")
        den = _den[0]; den_e = _den[1]
        _num = self.treeDA.getYields(lumi,'0.5',0,1,"("+self.ttControCuts+") && (%s)"%cuts.trigger)
        num = _num[0]; num_e = _num[1]
        _bkg = self.treeMCNoTT.getYields(lumi,'0.5',0,1,"("+self.ttControCuts+")")
        bkg = _bkg[0]; bkg_e = _bkg[1]
        DataMCFactor   = (num-bkg)/den
        DataMCFactor_e = num_e/den + bkg_e/den + num*den_e/den**2
        self.TTScaleFactor   = DataMCFactor
        self.TTScaleFactor_e = DataMCFactor_e
        print 'ttbar sf is', self.TTScaleFactor, '+/-', self.TTScaleFactor_e, self.sys

    def GetDataMCFactorTTZ(self):
        # tt bar is subtracted using the data/mc correction
        # estimation = observed in CR x Transfer factor
        if not hasattr(self,'TTScaleFactor'):
            print 'GetDataMCFactorTTZ shouldnt be called like that. GetDataMCFactorTT should be called first'
            print 'Calling GetDataMCFactorTT...'
            self.GetDataMCFactorTT()



        _den = self.treeTTZ.getYields(lumi,'0.5',0,1,self.ttzControCuts)
        den = _den[0]; den_e = _den[1]
        _num = self.treeDataDriven.getYields(lumi,'0.5',0,1,self.ttzControCuts)
        num = _num[0]; num_e = _num[1]
        _bkg = self.treeMCNoTTZNoTT.getYields(lumi,'0.5',0,1,self.ttzControCuts)
        bkg = _bkg[0]; bkg_e = _bkg[1]
        _tt  = self.treeTT.getYields(lumi,'0.5',0,1,self.ttzControCuts)
        tt = _tt[0]*self.TTScaleFactor; tt_e = _tt[1]*self.TTScaleFactor
        bkg = bkg + tt; bkg_e = bkg_e + tt_e


        DataMCFactor   = (num-bkg)/den
        DataMCFactor_e = num_e/den + bkg_e/den + num*den_e/den**2
        self.TTZScaleFactor  =DataMCFactor
        self.TTZScaleFactor_e=DataMCFactor_e
        print 'ttz sf is', self.TTZScaleFactor, '+/-', self.TTZScaleFactor_e
            
    def GetCorrectionFactors(self):
        self.GetDataMCFactorTT()
        self.GetDataMCFactorTTZ()
        

    def GetDataDriven(self):
        searchReg = self.treeMCNoTTZNoTT.getTH1F(lumi, "stack_some", 'met_Edge',10 ,100,500,
                                                  "(" + self.baselineCuts+")", "", "MET [GeV]")
        tt       = self.treeTT .getTH1F(lumi, "tt_sr",  'met_Edge',10 ,100,500,
                                        "(" + self.baselineCuts+")", "", "MET [GeV]")
        ttz      = self.treeTTZ.getTH1F(lumi, "ttz_sr", 'met_Edge',10 ,100,500,
                                        "(" + self.baselineCuts+")", "", "MET [GeV]")
        
        for bin in range(1,tt.GetNbinsX()):
            con_tt   = tt.GetBinContent(bin)
            con_tt_e = tt.GetBinError(bin)
            tt.SetBinContent(bin,con_tt*self.TTScaleFactor)
            tt.SetBinError(bin, con_tt_e*self.TTScaleFactor+con_tt*self.TTScaleFactor_e)
        for bin in range(1,ttz.GetNbinsX()):
            con_ttz  = ttz.GetBinContent(bin)
            con_ttz_e= ttz.GetBinError(bin)
            ttz.SetBinContent(bin,con_ttz*self.TTZScaleFactor)
            ttz.SetBinError(bin,con_ttz_e*self.TTZScaleFactor+con_ttz*self.TTZScaleFactor_e)

        sig     = self.treeSI.getTH1F(lumi, "sig",  'met_Edge',10 ,100,500,
                                      "(" + self.baselineCuts+")", "", "MET [GeV]")
        sig300  = self.treeSI300.getTH1F(lumi, "sig300",  'met_Edge',10 ,100,500,
                                      "(" + self.baselineCuts+")", "", "MET [GeV]")
        sig300.Scale(3.5*3.5)


        self.tt    = tt
        self.ttz   = ttz
        self.other = searchReg
        self.sig   = sig
        self.sig300= sig300




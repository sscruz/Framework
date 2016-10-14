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


class  darkmatter:
    def __init__(self,name):
        self.setTrees()
        self.setSRsCuts()
    def setTrees(self):
        mcDatasetsHT  = ['TT_pow_ext34','DYJetsToLL_M50_HT100to200_ext', 'DYJetsToLL_M50_HT200to400_ext',
                         'DYJetsToLL_M50_HT400to600_ext','DYJetsToLL_M50_HT600toInf_ext',
                         'WWTo2L2Nu', 'WZTo2L2Q', 'TTZ_LO','TTW_LO','TBar_tWch','ZZ','T_tWch']
        mcDatasetsNoTT  = ['DYJetsToLL_M50_HT100to200_ext', 'DYJetsToLL_M50_HT200to400_ext',
                           'DYJetsToLL_M50_HT400to600_ext','DYJetsToLL_M50_HT600toInf_ext',
                           'WWTo2L2Nu', 'WZTo2L2Q', 'TTZ_LO','TTW_LO','TBar_tWch','ZZ','T_tWch']
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
        
        self.treeMCHT   = Sample.Tree(helper.selectSamples('samples.dat', mcDatasetsHT,   'MC'), 'MC', 0, isScan = 0)
        self.treeTT     = Sample.Tree(helper.selectSamples('samples.dat', ttDatasets,   'MC'), 'MC', 0, isScan = 0)
        self.treeDA   = Sample.Tree(helper.selectSamples('samples.dat', daDatasets,   'DA'), 'DA', 1, isScan = 0)

# tree for data driven studies. choose treeDA when studies are done (choose treeMC for closure tests)
        self.treeDataDriven = self.treeDA
        self.treeMCNoTT     = Sample.Tree(helper.selectSamples('samples.dat', mcDatasetsNoTT,   'MC'),
                                          'MC', 0, isScan = 0)
        self.treeTTZ        = Sample.Tree(helper.selectSamples('samples.dat', ttZDatasets,   'MC'),
                                          'MC', 0, isScan = 0)
        self.treeMCNoTTZNoTT= Sample.Tree(helper.selectSamples('samples.dat', mcDatasetsNoTTZ,   'MC'),
                                          'MC', 0, isScan = 0)

    def setSRsCuts(self):
        theCuts = []
        self.baselineCuts = [cuts.goodLepton, 'nBJetMedium35_Edge > 0','mt2_Edge > 120',
                             'TMath::Abs(phi_ptboost) < 1.','nLepLoose_Edge < 3']
        self.ttControCuts = [cuts.goodLepton, 'nBJetMedium35_Edge > 0', 'TMath::Abs(phi_ptboost) < 1.',
                             'nLepLoose_Edge < 3','mt2_Edge > 60','mt2_Edge < 100','nJetSel_Edge>1']
        self.ttzControCuts = [cuts.goodLepton, 'nBJetMedium35_Edge > 0','mt2_Edge > 90', 'nLepTight_Edge > 2',
                              'nJetSel_Edge>1'] 

        self.baselineCuts  = self.getTheCuts(self.baselineCuts)
        self.ttControCuts  = self.getTheCuts(self.ttControCuts)
        self.ttzControCuts = self.getTheCuts(self.ttzControCuts,False)

    def getTheCuts(self,theCuts,doZVeto=True):
        finalcuts = ''
        for fl in ['ee','OF','mm']:
            if doZVeto:
                if len(finalcuts) == 0: 
                    finalcuts = "("+ cuts.AddList(theCuts+[getattr(cuts,fl)]+([] if fl == 'OF' else ['TMath::Abs(lepsMll_Edge - 91) > 20'])) + ")"
                else: 
                    finalcuts += "|| (" + cuts.AddList(theCuts+[getattr(cuts,fl)]+([] if fl == 'OF' else ['TMath::Abs(lepsMll_Edge - 91) > 20'])) + ")"
            else:
                if len(finalcuts) == 0: 
                    finalcuts = "("+ cuts.AddList(theCuts+[getattr(cuts,fl)]) + ")"
                else: 
                    finalcuts += "|| (" + cuts.AddList(theCuts+[getattr(cuts,fl)]) + ")"

        return finalcuts

    def doTTControlRegionStack(self):
        ttControl = self.treeMCHT.getStack(lumi, "ttCont_all", 'mt2_Edge',20 ,0,200,
                                           "(" + self.ttControCuts+")", "", "MET [GeV]")
        data      = self.treeDA.getTH1F(lumi,    "ttCont_data", 'mt2_Edge',20 ,0,200,
                                        "(" + self.ttControCuts+") && (%s)"%cuts.trigger, "", "MET [GeV]")
        data.SetMarkerStyle(r.kFullCircle)
        plot = Canvas.Canvas('DMStuff/TT_Control','png,pdf,cxx', 0.6, 0.60, 0.8, 0.8)
        plot.addStack(ttControl,'HIST',True,1)
        plot.addHisto(data,'P,SAME',"Data (%4.0f)"%data.Integral(),"P",r.kBlack,True,0)
        plot.save(True, True, False, lumi,0.1,5e4)

    def doTTZControlRegionStack(self):
        ttzControl = self.treeMCHT.getStack(lumi, "ttzCont_all", 'nJetSel_Edge',6 ,-0.5,5.5,
                                            "(" + self.ttzControCuts+")", "", "MET [GeV]")

        plot = Canvas.Canvas('DMStuff/TTZ_Control_nBJet','png,pdf,cxx', 0.6, 0.60, 0.8, 0.8)
        plot.addStack(ttzControl,'HIST',True,1)
        plot.save(True, True, False, lumi,0.1,5e4)

    def DoMCStackSR(self):
        searchReg = self.treeMCHT.getStack(lumi, "stack_all", 'met_Edge',10 ,100,500,
                                           "(" + self.baselineCuts+")", "", "MET [GeV]")
        plot = Canvas.Canvas('DMStuff/ATLAS_met_all','png,pdf,cxx', 0.6, 0.60, 0.8, 0.8)
        plot.addStack(stack,'HIST',True,1)
        plot.save(True, True, False, lumi,0.1,5e4)

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
        print 'ttbar sf is', self.TTScaleFactor, '+/-', self.TTScaleFactor_e

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

        searchReg = self.treeMCNoTTZNoTT.getStack(lumi, "stack_some", 'met_Edge',10 ,100,500,
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


        tt.SetFillColor(r.kCyan)
        ttz.SetFillColor(r.kRed)

        searchReg.Add(ttz);         searchReg.Add(tt)
        data     = self.treeDA.getTH1F(lumi, "data",  'met_Edge',10 ,100,500,
                                       "(" + self.baselineCuts+") && (%s)"%cuts.trigger, "", "MET [GeV]")
        plot = Canvas.Canvas('DMStuff/SR_data_driven','png,pdf,cxx', 0.6, 0.60, 0.8, 0.8)
        plot.addHisto(data,'P',"Data (%4.0f)"%data.Integral(),"P",r.kBlack,True,0)
        plot.addStack(searchReg,'HIST,SAME',True,1)
#        plot.addHisto(data,'P,SAME',"Data (%4.0f)"%data.Integral(),"P",r.kBlack,True,0)
        plot.save(True, True, False, lumi,0.1,5e4)
    
    def dump(self,filename):
        f = open(filename,'wb')
        for argname in ['TTZScaleFactor','TTZScaleFactor_e','TTScaleFactor','TTScaleFactor_e']:
            f.write('%s      %f\n'%(argname,getattr(self,argname)))
        f.close()

    def load(self,filename):
        f = open(filename,'rb')
        for line in f.read().splitlines():
            setattr(self,line.split()[0],float(line.split()[1]))

if __name__ == "__main__":
    
    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option("-r", action='store_true', dest="reload")
    (opts,arg) = parser.parse_args()
    
    gROOT.SetBatch(True)
    gROOT.ProcessLine('.L include/tdrstyle.C')
    r.setTDRStyle() 
    global lumi
    global cuts
    cuts = CutManager.CutManager()
    lumi= 12.9

    if os.path.exists('darkmatter.txt') and not opts.reload:
        print 'loading pickle file'
        analysis = darkmatter('theAnalysis')
        analysis.load('darkmatter.txt')

    else:
        print 'recreating darkmatter class. this will take time' 
        analysis = darkmatter('theAnalysis')
        analysis.GetCorrectionFactors()
        analysis.dump('darkmatter.txt')

    print 'done.'

    analysis.GetDataDriven()

import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TH1F, TGraph
import math, sys, optparse, copy, re, array, itertools, os


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables

gROOT.SetBatch(True)
gROOT.ProcessLine('.L include/tdrstyle.C')
r.setTDRStyle() 

mcDatasetsHT  = ['TT_pow_ext34','DYJetsToLL_M50_HT100to200_ext', 'DYJetsToLL_M50_HT200to400_ext',
                 'DYJetsToLL_M50_HT400to600_ext','DYJetsToLL_M50_HT600toInf_ext',
                 'WWTo2L2Nu', 'WZTo2L2Q', 'TTZ_LO','TTW_LO','TBar_tWch','ZZ','T_tWch']

daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2_runs_271036_276097',
              'DoubleEG_Run2016B-PromptReco-v2_runs_271036_276097',
              'MuonEG_Run2016B-PromptReco-v2_runs_271036_276097',
              'DoubleMuon_Run2016C-PromptReco-v2_runs_271036_276811',
              'DoubleEG_Run2016C-PromptReco-v2_runs_271036_276811',
              'MuonEG_Run2016C-PromptReco-v2_runs_271036_276811',
              'DoubleMuon_Run2016D-PromptReco-v2_runs_271036_276811',
              'DoubleEG_Run2016D-PromptReco-v2_runs_271036_276811',
              'MuonEG_Run2016D-PromptReco-v2_runs_271036_276811']

treeMCHT   = Sample.Tree(helper.selectSamples('samples.dat', mcDatasetsHT,   'MC'), 'MC', 0, isScan = 0)
#treeDA   = Sample.Tree(helper.selectSamples('samples.dat', daDatasets,   'DA'), 'DA', 1, isScan = 0)

cuts = CutManager.CutManager()
theCuts = [cuts.goodLepton]
#baselineCuts = [ 'nBJetMedium35_Edge > 0','nJetSel_Edge > 1', 'met_Edge > 50', 'TMath::Abs(phi_metzpt) > 1.2']
baselineCuts = [ 'nBJetMedium35_Edge > 0','mt2_Edge > 120', 'TMath::Abs(phi_ptboost) < 1.','nLepLoose_Edge < 3']
lumi = 12.9

finalcuts = ''

for fl in ['ee','OF','mm']:
    if len(finalcuts) == 0: 
        finalcuts = "("+ cuts.AddList(theCuts+baselineCuts+[getattr(cuts,fl)]+([] if fl == 'OF' else ['TMath::Abs(lepsMll_Edge - 91) > 20'])) + ")"
    else: 
        finalcuts += "|| (" + cuts.AddList(theCuts+baselineCuts+[getattr(cuts,fl)]+([] if fl == 'OF' else ['TMath::Abs(lepsMll_Edge - 91) > 20'])) + ")"


# [50,75,100,125,150,175,200,250,300]
print finalcuts 
stack = treeMCHT.getStack(lumi, "stack_%s"%('all'), 'met_Edge',10 ,100,500, "(" + finalcuts+")", "", "MET [GeV]")
#datos = treeDA.getTH1F(lumi, "th1_%s"%('all'), 'met_Edge', 10, 100,500, cuts.AddList([finalcuts,cuts.trigger]), "", "MET [GeV]")


plot = Canvas.Canvas('DMStuff/ATLAS_met_all','png,pdf,cxx', 0.6, 0.60, 0.8, 0.8)
plot.addStack(stack,'HIST',True,1)
#datos.SetMarkerStyle(r.kFullCircle)
#plot.addHisto(datos,'E,SAME', 'Data','P',r.kBlack,True,0)
plot.save(True, True, False, lumi,0.1,5e4)


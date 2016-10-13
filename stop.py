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


mcDatasetsHT  = ['TT_pow_ext3','TT_pow_ext4','DYJetsToLL_M50_HT100to200_ext', 'DYJetsToLL_M50_HT200to400_ext',
                 'DYJetsToLL_M50_HT400to600_ext','DYJetsToLL_M50_HT600toInf_ext',
                 'WWTo2L2Nu', 'WZTo2L2Q', 'TTZ_LO','TTW_LO','TBar_tWch']
                               
print 'poner zz','T_tWch'
siDatasets = ['SMS-T2tt_mStop-150to250']

#treeSI     = Sample.Tree(helper.selectSamples('samples_forstop.dat', siDatasets,     'SI'), 'SI', 0, isScan = 1)
treeMCHT   = Sample.Tree(helper.selectSamples('samples.dat', mcDatasetsHT,   'MC'), 'MC', 0, isScan = 0)

cuts = CutManager.CutManager()
theCuts = [cuts.goodLepton]
baselineCuts = ['nJetSel_Edge > 1', 'nBJetMedium35_Edge > 0', 'TMath::Abs(lepsMll_Edge - 91) > 15','met_Edge > 80', 'met_Edge/TMath::Sqrt(htJet35j_Edge) > 5', 'TMath::Cos(JetSel_Edge_phi[0]-met_phi_Edge) < 0.8', 'TMath::Cos(JetSel_Edge_phi[1]-met_phi_Edge)<TMath::Cos(0.25)']
lumi = 12.9
histos = []
plot = Canvas.Canvas('stopStuff/sig','png,pdf,cxx', 0.17, 0.70, 0.37, 0.8)
for fl in ['OF','ee','mm']:
    stack = treeMCHT.getStack(lumi, "stack_%s"%(fl), 'mt2_Edge', [0.,20.,40.,60.,80,100.,140.,240.], 0., 0., cuts.AddList(theCuts+baselineCuts+[getattr(cuts,fl)]), "", "m_{T2} [GeV]")
#    signal = treeSI.getTH1F(lumi, "stack_%s"%(fl), 'mt2_Edge', [0.,20.,40.,60.,80,100.,140.,240.], 0., 0., (cuts.AddList(theCuts+baselineCuts+[getattr(cuts,fl)]+['GenSusyMStop_Edge == 650 && GenSusyMNeutralino_Edge ==1'])).replace(cuts.twoLeptons, 't.nPairLep_Edge > 0'), "", "m_{T2} [GeV]")
#    signal.Scale( 1 / signal.Integral() )
#    plot = Canvas.Canvas('stopStuff/mt2_%s'%fl,'png,pdf,cxx', 0.6, 0.60, 0.8, 0.8)
    plot.addStack(stack,'HIST',True,0)

    color = r.kBlue if fl == 'ee' else r.kRed if fl == 'mm' else r.kBlack
##    plot.addHisto(signal,'HIST,SAME','T2tt (650,1) %s'%fl,'L',color,True,True)

plot.save(True, False, True, lumi)

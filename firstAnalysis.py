#####################################################################
######                                                              #
###### 88888888888         88                        888888888888   #  
###### 88                  88                                 ,88   #
###### 88                  88                               ,88"    #  
###### 88aaaaa     ,adPPYb,88  ,adPPYb,d8  ,adPPYba,      ,88"      #
###### 88"""""    a8"    `Y88 a8"    `Y88 a8P_____88    ,88"        #
###### 88         8b       88 8b       88 8PP"""""""  ,88"          #
###### 88         "8a,   ,d88 "8a,   ,d88 "8b,   ,aa 88"            #
###### 88888888888 `"8bbdP"Y8  `"YbbdP"Y8  `"Ybbd8"' 888888888888   #
######                       aa,    ,88                             #
######                         "Y8bbdP"                             #
######                                                              #
#####################################################################

import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors
import math, sys, optparse, array
import Rounder as rounder

import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample




if __name__ == "__main__":


    parser = optparse.OptionParser(usage="usage: %prog [options] FilenameWithSamples", version="%prog 1.0")
    parser.add_option("-m", "--mode", action="store", dest="mode", default="rmue", help="Operation mode")
    (options, args) = parser.parse_args()

    inputFileName = args[0]

    print 'Going to load DATA and MC trees...'
    mcDatasets = ['TTJets_LO', 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    dyDatasets = ['DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    daDatasets = ['DoubleMuon_Run2015C', 'DoubleEG_Run2015C', 'MuonEG_Run2015C',
                  'DoubleMuon_Run2015D', 'DoubleEG_Run2015D', 'MuonEG_Run2015D']
    treeMC = Sample.Tree(helper.selectSamples(inputFileName, mcDatasets, 'MC'), 'MC'  , 0)
    treeDY = Sample.Tree(helper.selectSamples(inputFileName, dyDatasets, 'DY'), 'DY'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(inputFileName, daDatasets, 'DA'), 'DATA', 1)
    #tree = treeMC
    print 'Trees successfully loaded...'

    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    ####Cuts needed by rmue
    cuts = CutManager.CutManager()

    lumi = 0.225


    regions = []
    setLog = []
    Control2JetsSF = Region.region('Control2JetsSF',
                       [cuts.Control2JetsSF()],
                       ['mll', 'met'],
                       [range(10,310,10), range(10,310,10)],
                       True)
    regions.append(Control2JetsSF)
    setLog.append(True)
#    Control2JetsOF = Region.region('Control2JetsOF',
#                        [cuts.Control2JetsOF()],
#                        ['mll', 'met'],
#                        [range(10,310,10), range(10,310,10),],
#                        True)
#    regions.append(Control2JetsOF) 
#    setLog.append(False)                      
#    DYControlNoMassLeptonSF = Region.region('DYControlNoMassLeptonSF',
#                        [cuts.DYControlNoMassLeptonSF()],
#                        ['mll', 'met'],
#                        [range(10,310,10), range(10,310,10),],
#                        True)
#    regions.append(DYControlNoMassLeptonSF)
#    setLog.append(True)
#    DYControlNoMassLeptonOF = Region.region('DYControlNoMassLeptonOF',
#                        [cuts.DYControlNoMassLeptonOF()],
#                        ['mll', 'met'],
#                        [range(10,310,10), range(10,310,10),],
#                        True)
#    regions.append(DYControlNoMassLeptonOF)                       
#    setLog.append(False)
#    ControlNoMassLeptonSF = Region.region('ControlNoMassLeptonSF',
#                        [cuts.ControlNoMassLeptonSF()],
#                        ['mll', 'met'],
#                        [range(10,310,10), range(10,310,10),],
#                        True)
#    regions.append(ControlNoMassLeptonSF)
#    setLog.append(True)
#    ControlNoMassLeptonOF = Region.region('ControlNoMassLeptonOF',
#                        [cuts.ControlNoMassLeptonOF()],
#                        ['mll', 'met'],
#                        [range(10,310,10), range(10,310,10),],
#                        True)
#    regions.append(ControlNoMassLeptonOF)                       
#    setLog.append(False)
#    Control2JetsMETee = Region.region('Control2JetsMETee',
#                        [cuts.Control2JetsMETee()],
#                        ['mll', 'met'],
#                        [range(10,310,10), range(10,310,10),],
#                        True)
#    regions.append(Control2JetsMETee)
#    setLog.append(True)				   
#    Control2JetsMETmm = Region.region('Control2JetsMETmm',
#                        [cuts.Control2JetsMETmm()],
#                        ['mll', 'met'],
#                        [range(10,310,10), range(10,310,10),],
#                        True)
#    regions.append(Control2JetsMETmm)     
#    setLog.append(True)    
    for reg in regions:
        print 'i am at region', reg.name
        for eta in ['central', 'forward']:
					
		my_cuts = cuts.AddList([cuts.goodLepton]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts)
		for tree in ([treeMC, treeDA] if reg.doData else [treeMC]):
		
#			dataMC = 'DATA' if tree == treeDA else 'MC'  
		
			if tree == treeDA: 
				dataMC = 'DATA'
			elif tree == treeMC: 
				dataMC = 'MC' 
			else: 
				dataMC = 'DY' 

			if 'mll' in reg.rvars:
				reg.mll.setHisto(tree.getTH1F(lumi, "mll_"+eta+reg.name+str(dataMC), "t.lepsMll_Edge", 	reg.bins[reg.rvars.index('mll')], 1, 1, my_cuts, "", "m_{ll} (GeV)"), dataMC, eta)
				reg.mll_dy.setHisto(treeDY.getTH1F(lumi, "mll_dy_"+eta+reg.name+str(dataMC), "t.lepsMll_Edge", 	reg.bins[reg.rvars.index('mll')], 1, 1, my_cuts, "", "m_{ll} (GeV)"), 'DY', eta)
			if 'met' in reg.rvars:
				reg.met.setHisto(tree.getTH1F(lumi, "met_"+eta+reg.name+str(dataMC), "met_pt", 	reg.bins[reg.rvars.index('met')], 1, 1, my_cuts, "", "ME_{T} (GeV)"), dataMC, eta)

    for reg in regions:
		print 'plotting region', reg.name

		for eta in ['central', 'forward']:
			plot_mll = Canvas.Canvas("test/mll_"+eta+"_"+reg.name, "png,pdf", 0.6, 0.6, 0.8, 0.8)
			print 'mll making some canvas' 

			plot_mll.addHisto(reg.mll.getHisto('MC', eta), "HIST", "MC"  , "PL", r.kRed+1 , 1, 0)		
			plot_mll.addHisto(reg.mll_dy.getHisto('DY', eta), "SAME", "DY"  , "PL", r.kRed+1 , 1, 0)		
			print 'mll adding MC histo'
			plot_mll.addHisto(reg.mll.getHisto('DATA', eta), "E,SAME", "Data", "PL", r.kBlack , 1, 1)
			print 'mll adding data histo'
			plot_mll.saveRatio(1, 1, 0, lumi, reg.mll.getHisto('DATA', eta), reg.mll.getHisto('MC', eta))
			print 'mll saving ratio'
			del plot_mll  
			plot_met = Canvas.Canvas("test/met_"+eta+"_"+reg.name, "png,pdf", 0.6, 0.6, 0.8, 0.8)
			print 'met making some canvas' 
			plot_met.addHisto(reg.met.getHisto('MC', eta), "HIST", "MC"  , "PL", r.kRed+1 , 1, 0)		
			print 'met adding MC histo'
			plot_met.addHisto(reg.met.getHisto('DATA', eta), "E,SAME", "Data", "PL", r.kBlack , 1, 1)
			print 'met adding data histo'
			plot_met.saveRatio(1, 1, 0, lumi, reg.met.getHisto('DATA', eta), reg.met.getHisto('MC', eta))
			del plot_met    

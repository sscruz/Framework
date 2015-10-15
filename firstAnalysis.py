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
import math, sys, optparse, array, time
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
 #   daDatasets = ['DoubleMuon_Run2015C', 'DoubleEG_Run2015C', 'MuonEG_Run2015C',
 #                 'DoubleMuon_Run2015D', 'DoubleEG_Run2015D', 'MuonEG_Run2015D']
    treeMC = Sample.Tree(helper.selectSamples(inputFileName, mcDatasets, 'MC'), 'MC'  , 0)
    treeDY = Sample.Tree(helper.selectSamples(inputFileName, dyDatasets, 'DY'), 'DY'  , 0)
 #   treeDA = Sample.Tree(helper.selectSamples(inputFileName, daDatasets, 'DA'), 'DATA', 1)
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
    Control2JetsMETSF = Region.region('Control2JetsMETSF',
                       [cuts.Control2JetsMETSF()],
                       ['min_mlb', 'max_mlb'],
                       [range(10,310,10), range(10,310,10)],
                       True)
    regions.append(Control2JetsMETSF)
    setLog.append(True)
    Control2JetsMETOF = Region.region('Control2JetsMETOF',
                        [cuts.Control2JetsMETOF()],
                        ['min_mlb', 'max_mlb'],
                        [range(10,310,10), range(10,310,10)],
                        True)
    regions.append(Control2JetsMETOF) 
    setLog.append(False)                      
#    RSFOFControlRegionAlt = Region.region('RSFOFControlRegionAlternative',
#                        [cuts.RSFOFControlRegion()],
#                        ['mll', 'met', 'nb', 'nj'],
#                        [range(10,310,10), range(10,310,10), range(0,8,1), range(0,8,1)],
#                        True)
#    regions.append(RSFOFControlRegionAlt) 
#    setLog.append(False)                      
#    RSFOFControlRegion = Region.region('RSFOFControlRegion',
#                        [cuts.RSFOFControlRegion()],
#                        ['mll', 'met', 'nb', 'nj'],
#                        [range(10,310,10), range(10,310,10), range(0,8,1), range(0,8,1)],
#                        True)
#    regions.append(RSFOFControlRegion) 
#    setLog.append(False)                      
    DYControlNoMassLeptonSF = Region.region('DYControlNoMassLeptonSF',
                        [cuts.DYControlNoMassLeptonSF()],
                        ['min_mlb', 'max_mlb'],
                        [range(10,310,10), range(10,310,10),],
                        True)
    regions.append(DYControlNoMassLeptonSF)
    setLog.append(True)
    DYControlNoMassLeptonOF = Region.region('DYControlNoMassLeptonOF',
                        [cuts.DYControlNoMassLeptonOF()],
                        ['min_mlb', 'max_mlb'],
                        [range(10,310,10), range(10,310,10),],
                        True)
    regions.append(DYControlNoMassLeptonOF)                       
    setLog.append(False)
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
#    Control2Jetsee = Region.region('Control2Jetsee',
#                        [cuts.Control2Jetsee()],
#                        ['met', 'nvtx'],
#                        [range(10,310,10), range(0,40,1)],
#                        True)
#    regions.append(Control2Jetsee)
#    setLog.append(True)                   
#    Control2Jetsmm = Region.region('Control2Jetsmm',
#                        [cuts.Control2Jetsmm()],
#                        ['met', 'nvtx'],
#                        [range(10,310,10), range(0,40,1)],
#                        True)
#    regions.append(Control2Jetsmm)     
#    setLog.append(True)    
#    InclusiveCROF = Region.region('InclusiveCROF',
#                       [cuts.InclusiveCROF()],
#                       ['mll', 'met', 'nb', 'nj'],
#                       [range(10,310,10), range(10,310,10), range(0,8,1), range(0,8,1)],
#                       True)
#    regions.append(InclusiveCROF)
#    setLog.append(True)                   
#    InclusiveCRSF = Region.region('InclusiveCRSF',
#                       [cuts.InclusiveCRSF()],
#                       ['mll', 'met', 'nb', 'nj'],
#                       [range(10,310,10), range(10,310,10), range(0,8,1), range(0,8,1)],
#                       True)
#    regions.append(InclusiveCRSF)     
#    setLog.append(True)    
    for reg in regions:
        print 'i am at region', reg.name
        for eta in ['central', 'forward']:
                    
            my_cuts = cuts.AddList([cuts.goodLepton]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts)
       #     for tree in ([treeMC, treeDA] if reg.doData else [treeMC]):
            for tree in [treeMC]:
 #               if tree == treeDA: 
 #                   dataMC = 'DATA'
                if tree == treeMC: 
                    dataMC = 'MC' 
                else: 
                    dataMC = 'DY' 

                if 'mll' in reg.rvars:
                    reg.mll.setHisto(tree.getTH1F(lumi, "mll_"+eta+reg.name+str(dataMC), "t.lepsMll_Edge",     reg.bins[reg.rvars.index('mll')], 1, 1, my_cuts, "", "m_{ll} (GeV)"), dataMC, eta)
                    if dataMC == "MC":
                         reg.mll_dy.setHisto(treeDY.getTH1F(lumi, "mll_dy_"+eta+reg.name+str(dataMC), "t.lepsMll_Edge",     reg.bins[reg.rvars.index('mll')], 1, 1, my_cuts, "", "m_{ll} (GeV)"), 'MC', eta)
                if 'met' in reg.rvars:
                    reg.met.setHisto(tree.getTH1F(lumi, "met_"+eta+reg.name+str(dataMC), "met_pt",     reg.bins[reg.rvars.index('met')], 1, 1, my_cuts, "", "ME_{T} (GeV)"), dataMC, eta)
                    if dataMC == "MC":
                        reg.met_dy.setHisto(treeDY.getTH1F(lumi, "met_dy_"+eta+reg.name+str(dataMC), "met_pt",     reg.bins[reg.rvars.index('met')], 1, 1, my_cuts, "", "ME_{T} (GeV)"), 'MC', eta)
                if 'nb' in reg.rvars:
                    reg.nb.setHisto(tree.getTH1F(lumi, "nb_"+eta+reg.name+str(dataMC), "t.nBJetMedium35_Edge",     reg.bins[reg.rvars.index('nb')], 1, 1, my_cuts, "", "nBJets, med35"), dataMC, eta)
                    if dataMC == "MC":
                        reg.nb_dy.setHisto(treeDY.getTH1F(lumi, "nb_dy_"+eta+reg.name+str(dataMC), "t.nBJetMedium35_Edge",     reg.bins[reg.rvars.index('nb')], 1, 1, my_cuts, "", "nBJets, med35"), 'MC', eta)
                if 'nj' in reg.rvars:
                    reg.nj.setHisto(tree.getTH1F(lumi, "nj_"+eta+reg.name+str(dataMC), "t.nJetSel_Edge",     reg.bins[reg.rvars.index('nj')], 1, 1, my_cuts, "", "nJets"), dataMC, eta)
                    if dataMC == "MC":
                        reg.nj_dy.setHisto(treeDY.getTH1F(lumi, "nj_dy_"+eta+reg.name+str(dataMC), "t.nJetSel_Edge",     reg.bins[reg.rvars.index('nj')], 1, 1, my_cuts, "", "nJets"), 'MC', eta)
                if 'min_mlb' in reg.rvars:
                    reg.min_mlb.setHisto(tree.getTH1F(lumi, "min_mlb_"+eta+reg.name+str(dataMC), "t.min_mlb1_Edge ",     reg.bins[reg.rvars.index('min_mlb')], 1, 1, my_cuts, "", "min Mlb"), dataMC, eta)
                    if dataMC == "MC":
                        reg.min_mlb_dy.setHisto(treeDY.getTH1F(lumi, "min_mlb_dy_"+eta+reg.name+str(dataMC), "t.min_mlb1_Edge",     reg.bins[reg.rvars.index('min_mlb')], 1, 1, my_cuts, "", "min Mlb"), 'MC', eta)
                if 'max_mlb' in reg.rvars:
                    reg.max_mlb.setHisto(tree.getTH1F(lumi, "max_mlb_"+eta+reg.name+str(dataMC), "t.min_mlb2_Edge ",     reg.bins[reg.rvars.index('max_mlb')], 1, 1, my_cuts, "", "next to min Mlb"), dataMC, eta)
                    if dataMC == "MC":
                        reg.max_mlb_dy.setHisto(treeDY.getTH1F(lumi, "max_mlb_dy_"+eta+reg.name+str(dataMC), "t.min_mlb2_Edge",     reg.bins[reg.rvars.index('max_mlb')], 1, 1, my_cuts, "", "next to min Mlb"), 'MC', eta)

                    





    for reg in regions:
        print 'plotting region', reg.name

        for eta in ['central', 'forward']:
#            reg.mll_dy.getHisto('MC', eta).SetFillColor(r.kBlue+1)
#            reg.mll_dy.getHisto('MC', eta).SetFillColorAlpha(r.kBlue+1, 0.5)
#            plot_mll = Canvas.Canvas("mll/mll_"+eta+"_"+reg.name, "png,pdf", 0.6, 0.6, 0.8, 0.8)
#            plot_mll.addHisto(reg.mll   .getHisto('MC', eta), "HIST", "DY+TTJets"  , "PL", r.kRed+1 , 1, 0)        
#            plot_mll.addHisto(reg.mll_dy.getHisto('MC', eta), "HIST SAME", "DY"  , "PL", r.kBlue+1 , 1, 0)        
#            plot_mll.addHisto(reg.mll   .getHisto('DATA', eta), "E,SAME", "Data", "PL", r.kBlack , 1, 1)
#            plot_mll.saveRatio(1, 1, 1, lumi, reg.mll.getHisto('DATA', eta), reg.mll.getHisto('MC', eta))
#            #del plot_mll  
#            time.sleep(0.1)

#            reg.met_dy.getHisto('MC', eta).SetFillColor(r.kBlue+1)
#            reg.met_dy.getHisto('MC', eta).SetFillColorAlpha(r.kBlue+1, 0.5)
#            plot_met = Canvas.Canvas("met/met_"+eta+"_"+reg.name, "png,pdf", 0.6, 0.6, 0.8, 0.8)
#            plot_met.addHisto(reg.met.getHisto('MC', eta), "HIST", "DY+TTJets"  , "PL", r.kRed+1 , 1, 0)        
#            plot_met.addHisto(reg.met_dy.getHisto('MC', eta), "HIST SAME", "DY"  , "PL", r.kBlue+1 , 1, 0)        
#            plot_met.addHisto(reg.met.getHisto('DATA', eta), "E,SAME", "Data", "PL", r.kBlack , 1, 1)
#            plot_met.saveRatio(1, 1, 1, lumi, reg.met.getHisto('DATA', eta), reg.met.getHisto('MC', eta))
#            #del plot_met    
#            time.sleep(0.1)

#            reg.nb_dy.getHisto('MC', eta).SetFillColor(r.kBlue+1)
#            reg.nb_dy.getHisto('MC', eta).SetFillColorAlpha(r.kBlue+1, 0.5)
#            plot_nb = Canvas.Canvas("nb/nb_"+eta+"_"+reg.name, "png,pdf", 0.6, 0.6, 0.8, 0.8)
#            plot_nb.addHisto(reg.nb.getHisto('MC', eta), "HIST", "MC"  , "PL", r.kRed+1 , 1, 0)        
#            plot_nb.addHisto(reg.nb_dy.getHisto('MC', eta), "HIST SAME", "DY"  , "PL", r.kBlue+1 , 1, 0)        
#            plot_nb.addHisto(reg.nb.getHisto('DATA', eta), "E,SAME", "Data", "PL", r.kBlack , 1, 1)
#            plot_nb.saveRatio(1, 1, 1, lumi, reg.nb.getHisto('DATA', eta), reg.nb.getHisto('MC', eta))
#            #del plot_nb    
#            time.sleep(0.1)

#            reg.nj_dy.getHisto('MC', eta).SetFillColor(r.kBlue+1)
#            reg.nj_dy.getHisto('MC', eta).SetFillColorAlpha(r.kBlue+1, 0.5)
#            plot_nj = Canvas.Canvas("nj/nj_"+eta+"_"+reg.name, "png,pdf", 0.6, 0.6, 0.8, 0.8)
#            plot_nj.addHisto(reg.nj.getHisto('MC', eta), "HIST", "MC"  , "PL", r.kRed+1 , 1, 0)        
#            plot_nj.addHisto(reg.nj_dy.getHisto('MC', eta), "HIST SAME", "DY"  , "PL", r.kBlue+1 , 1, 0)        
#            plot_nj.addHisto(reg.nj.getHisto('DATA', eta), "E,SAME", "Data", "PL", r.kBlack , 1, 1)
#            plot_nj.saveRatio(1, 1, 1, lumi, reg.nj.getHisto('DATA', eta), reg.nj.getHisto('MC', eta))
#            #del plot_nj    
#            time.sleep(0.1)

            reg.min_mlb_dy.getHisto('MC', eta).SetFillColor(r.kBlue+1)
            reg.min_mlb_dy.getHisto('MC', eta).SetFillColorAlpha(r.kBlue+1, 0.5)
            plot_min_mlb = Canvas.Canvas("min_mlb/min_mlb_"+eta+"_"+reg.name, "png,pdf", 0.6, 0.6, 0.8, 0.8)
            plot_min_mlb.addHisto(reg.min_mlb.getHisto('MC', eta), "HIST", "MC"  , "PL", r.kRed+1 , 1, 0)        
            plot_min_mlb.addHisto(reg.min_mlb_dy.getHisto('MC', eta), "HIST SAME", "DY"  , "PL", r.kBlue+1 , 1, 0)        
#            plot_min_mlb.addHisto(reg.min_mlb.getHisto('DATA', eta), "E,SAME", "Data", "PL", r.kBlack , 1, 1)
            plot_min_mlb.saveRatio(1, 1, 0, lumi, reg.min_mlb.getHisto('MC', eta), reg.min_mlb.getHisto('MC', eta))
            #del plot_min_mlb1    
            time.sleep(0.1)

 
            reg.max_mlb_dy.getHisto('MC', eta).SetFillColor(r.kBlue+1)  
            reg.max_mlb_dy.getHisto('MC', eta).SetFillColorAlpha(r.kBlue+1, 0.5)
            plot_max_mlb = Canvas.Canvas("min_mlb/max_mlb_"+eta+"_"+reg.name, "png,pdf", 0.6, 0.6, 0.8, 0.8)
            plot_max_mlb.addHisto(reg.max_mlb.getHisto('MC', eta), "HIST", "MC"  , "PL", r.kRed+1 , 1, 0)        
            plot_max_mlb.addHisto(reg.max_mlb_dy.getHisto('MC', eta), "HIST SAME", "DY"  , "PL", r.kBlue+1 , 1, 0)     
#            plot_max_mlb.addHisto(reg.max_mlb.getHisto('DATA', eta), "E,SAME", "Data", "PL", r.kBlack , 1, 1)
            plot_max_mlb.saveRatio(1, 1, 0, lumi, reg.max_mlb.getHisto('MC', eta), reg.max_mlb.getHisto('MC', eta))
            time.sleep(0.1)  


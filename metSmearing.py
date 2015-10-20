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
    mcDatasets = ['TTLep_pow', 'DYJetsToLL_M50', 'DYJetsToLL_M10to50', 'WWTo2L2Nu', 'WZTo2L2Q', 'ZZTo2L2Q']
    dyDatasets = ['DYJetsToLL_M50', 'DYJetsToLL_M10to50']
    daDatasets = ['DoubleMuon_Run2015C', 'DoubleEG_Run2015C', 'MuonEG_Run2015C',
                  'DoubleMuon_Run2015D', 'DoubleEG_Run2015D', 'MuonEG_Run2015D']
    treeMC = Sample.Tree(helper.selectSamples(inputFileName, mcDatasets, 'MC'), 'MC'  , 0)
    treeDY = Sample.Tree(helper.selectSamples(inputFileName, dyDatasets, 'DY'), 'DY'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(inputFileName, daDatasets, 'DA'), 'DATA', 1)
    print 'Trees successfully loaded...'

    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    ####Cuts needed by rmue
    cuts = CutManager.CutManager()

    lumi = 0.225


                       
   # InclusiveCRSF = Region.region('InclusiveCRSF',
   #                     [cuts.InclusiveCRSF()],
   #                     [ 'met', 'jzb', 'met_x', 'met_y'],
   #                     [range(10,310,10), range(-260,260,10), range(-260,260,10),range(-260,260,10)],
   #                     True)
   # regions.append(InclusiveCRSF) 
   # setLog.append(True)
                        
   # for tree in ([treeMC, treeDA] if reg.doData else [treeMC]):
    #    if tree == treeDA: 
    #        dataMC = 'DATA'
    #    else:  
    #        dataMC = 'MC' 

    met_InclusiveSF_data = treeDA.getTH1F(lumi, "met_InclusiveSF_data", "met_pt", 40, 20, 300, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} (GeV)")
    met_InclusiveSF_mc = treeMC.getTH1F(lumi, "met_InclusiveSF_mc", "met_pt", 40, 20, 300, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} (GeV)")   
    met_InclusiveSF_dy = treeDY.getTH1F(lumi, "met_InclusiveSF_dy", "met_pt", 40, 20, 300, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} (GeV)")   
              

    met_x_InclusiveSF_data = treeDA.getTH1F(lumi, "met_x_InclusiveSF_data", "met_pt*cos(met_phi)", 50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} x (GeV)")
    met_x_InclusiveSF_mc = treeMC.getTH1F(lumi, "met_x_InclusiveSF_mc", "met_pt*cos(met_phi)", 50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} x (GeV)")   
    met_x_InclusiveSF_dy = treeDY.getTH1F(lumi, "met_x_InclusiveSF_dy", "met_pt*cos(met_phi)",  50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} x (GeV)")   


    met_y_InclusiveSF_data = treeDA.getTH1F(lumi, "met_y_InclusiveSF_data", "met_pt*sin(met_phi)", 50, -260, 260 , cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} y (GeV)")
    met_y_InclusiveSF_mc = treeMC.getTH1F(lumi, "met_y_InclusiveSF_mc", "met_pt*sin(met_phi)",  50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} y (GeV)")   
    met_y_InclusiveSF_dy = treeDY.getTH1F(lumi, "met_y_InclusiveSF_dy", "met_pt*sin(met_phi)",  50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} y (GeV)")   
       
       
    jzb_InclusiveSF_data = treeDA.getTH1F(lumi, "jzb_InclusiveSF_data", "lepsJZB_Edge",  50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "JZB (GeV)")
    jzb_InclusiveSF_mc = treeMC.getTH1F(lumi, "jzb_InclusiveSF_mc", "lepsJZB_Edge",  50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "JZB (GeV)")   
    jzb_InclusiveSF_dy = treeDY.getTH1F(lumi, "jzb_InclusiveSF_dy", "lepsJZB_Edge",  50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "JZB (GeV)")   



    met_InclusiveSF_dy.SetFillColor(r.kBlue+1)                     
    met_InclusiveSF_dy.SetFillColorAlpha(r.kBlue+1, 0.5)            
    met_InclusiveSF = Canvas.Canvas("metSmearing/met_InclusiveSF", "png", 0.6, 0.6, 0.8, 0.8)                    
    met_InclusiveSF.addHisto(met_InclusiveSF_mc, "HIST", "DY+TTJets", "P", r.kRed+1, 1, 0)
    met_InclusiveSF.addHisto(met_InclusiveSF_dy, "HIST SAME", "DY"  , "PL", r.kBlue+1 , 1, 0) 
    met_InclusiveSF.addHisto(met_InclusiveSF_data, "E,SAME", "Data"  , "PL", r.kBlack , 1, 1) 
    met_InclusiveSF.saveRatio(1, 1, 1, lumi, met_InclusiveSF_data, met_InclusiveSF_mc)
    time.sleep(0.1)                                                                       
        
    met_x_InclusiveSF_dy.SetFillColor(r.kBlue+1)                             
    met_x_InclusiveSF_dy.SetFillColorAlpha(r.kBlue+1, 0.5)         
    met_x_InclusiveSF = Canvas.Canvas("metSmearing/met_x_InclusiveSF", "png", 0.6, 0.6, 0.8, 0.8)                 
    met_x_InclusiveSF.addHisto(met_x_InclusiveSF_mc, "HIST", "DY+TTJets", "P", r.kRed+1, 1, 0)
    met_x_InclusiveSF.addHisto(met_x_InclusiveSF_dy, "HIST SAME", "DY"  , "PL", r.kBlue+1 , 1, 0) 
    met_x_InclusiveSF.addHisto(met_x_InclusiveSF_data, "E,SAME", "Data"  , "PL", r.kBlack , 1, 1) 
    met_x_InclusiveSF.saveRatio(1, 1, 1, lumi, met_x_InclusiveSF_data, met_x_InclusiveSF_mc)
    time.sleep(0.1)                                                                  

    met_y_InclusiveSF_dy.SetFillColor(r.kBlue+1)                     
    met_y_InclusiveSF_dy.SetFillColorAlpha(r.kBlue+1, 0.5)         
    met_y_InclusiveSF = Canvas.Canvas("metSmearing/met_y_InclusiveSF", "png", 0.6, 0.6, 0.8, 0.8)                 
    met_y_InclusiveSF.addHisto(met_y_InclusiveSF_mc, "HIST", "DY+TTJets", "P", r.kRed+1, 1, 0)
    met_y_InclusiveSF.addHisto(met_y_InclusiveSF_dy, "HIST SAME", "DY"  , "PL", r.kBlue+1 , 1, 0) 
    met_y_InclusiveSF.addHisto(met_y_InclusiveSF_data, "E,SAME", "Data"  , "PL", r.kBlack , 1, 1) 
    met_y_InclusiveSF.saveRatio(1, 1, 1, lumi, met_y_InclusiveSF_data, met_y_InclusiveSF_mc)
    time.sleep(0.1)                                                                  

    jzb_InclusiveSF_dy.SetFillColor(r.kBlue+1)                     
    jzb_InclusiveSF_dy.SetFillColorAlpha(r.kBlue+1, 0.5)         
    jzb_InclusiveSF = Canvas.Canvas("metSmearing/jzb_InclusiveSF", "png", 0.6, 0.6, 0.8, 0.8)                 
    jzb_InclusiveSF.addHisto(jzb_InclusiveSF_mc, "HIST", "DY+TTJets", "P", r.kRed+1, 1, 0)
    jzb_InclusiveSF.addHisto(jzb_InclusiveSF_dy, "HIST SAME", "DY"  , "PL", r.kBlue+1 , 1, 0) 
    jzb_InclusiveSF.addHisto(jzb_InclusiveSF_data, "E,SAME", "Data"  , "PL", r.kBlack , 1, 1) 
    jzb_InclusiveSF.saveRatio(1, 1, 1, lumi,  jzb_InclusiveSF_data, jzb_InclusiveSF_mc)
    time.sleep(0.1)                                                                  






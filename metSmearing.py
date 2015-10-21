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
    mcDatasets = ['TTLep_pow']
  #  mcDatasets = ['TTLep_pow', 'DYJetsToLL_M50', 'DYJetsToLL_M10to50', 'WWTo2L2Nu', 'WZTo2L2Q', 'ZZTo2L2Q']
  #  dyDatasets = ['DYJetsToLL_M50', 'DYJetsToLL_M10to50']
  #  daDatasets = ['DoubleEG_Run2015D_05Oct_v1_runs_246908_258751', 'DoubleEG_Run2015D_v4_runs_246908_258751', 'DoubleMuon_Run2015D_05Oct_v1_runs_246908_258751', 'DoubleMuon_Run2015D_v4_runs_246908_258751', 
  #              'MuonEG_Run2015D_05Oct_v2_runs_246908_258751', 'MuonEG_Run2015D_v4_runs_246908_258751']
    treeMC = Sample.Tree(helper.selectSamples(inputFileName, mcDatasets, 'MC'), 'MC'  , 0)
  #  treeDY = Sample.Tree(helper.selectSamples(inputFileName, dyDatasets, 'DY'), 'DY'  , 0)
  #  treeDA = Sample.Tree(helper.selectSamples(inputFileName, daDatasets, 'DA'), 'DATA', 1)
    print 'Trees successfully loaded...'

    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()

    lumi = 1.264

 #   met_InclusiveSF_data = treeDA.getTH1F(lumi, "met_InclusiveSF_data", "met_pt", 40, 20, 300, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} (GeV)")
 #   met_InclusiveSF_mc = treeMC.getTH1F(lumi, "met_InclusiveSF_mc", "met_pt", 40, 20, 300, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} (GeV)")   
 #   met_InclusiveSF_dy = treeDY.getTH1F(lumi, "met_InclusiveSF_dy", "met_pt", 40, 20, 300, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} (GeV)")   

    met_InclusiveOF = treeMC.getTH1F(lumi, "met_InclusiveOF", "met_pt", 40, 20, 300, cuts.AddList([cuts.InclusiveCROF()]), "", "ME_{T} (GeV)")
    met_InclusiveEE = treeMC.getTH1F(lumi, "met_InclusiveEE", "met_pt", 40, 20, 300, cuts.AddList([cuts.InclusiveCRee()]), "", "ME_{T} (GeV)")   
    met_InclusiveMM = treeMC.getTH1F(lumi, "met_InclusiveMM", "met_pt", 40, 20, 300, cuts.AddList([cuts.InclusiveCRmm()]), "", "ME_{T} (GeV)")   
    met_InclusiveSF = treeMC.getTH1F(lumi, "met_InclusiveSF", "met_pt", 40, 20, 300, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} (GeV)") 
#    met_raw_InclusiveSF_data = treeDA.getTH1F(lumi, "met_raw_InclusiveSF_data", "met_rawPt", 40, 20, 300, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} (Raw) (GeV)")
#    met_raw_InclusiveSF_mc = treeMC.getTH1F(lumi, "met_raw_InclusiveSF_mc", "met_rawPt", 40, 20, 300, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} (Raw) (GeV)")   
#    met_raw_InclusiveSF_dy = treeDY.getTH1F(lumi, "met_raw_InclusiveSF_dy", "met_rawPt", 40, 20, 300, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} (Raw) (GeV)")   

 #   met_x_InclusiveSF_data = treeDA.getTH1F(lumi, "met_x_InclusiveSF_data", "met_pt*cos(met_phi)", 50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} x (GeV)")
 #   met_x_InclusiveSF_mc = treeMC.getTH1F(lumi, "met_x_InclusiveSF_mc", "met_pt*cos(met_phi)", 50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} x (GeV)")   
 #   met_x_InclusiveSF_dy = treeDY.getTH1F(lumi, "met_x_InclusiveSF_dy", "met_pt*cos(met_phi)",  50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} x (GeV)")   
       
  #  met_y_InclusiveSF_data = treeDA.getTH1F(lumi, "met_y_InclusiveSF_data", "met_pt*sin(met_phi)", 50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} y (GeV)")
  #  met_y_InclusiveSF_mc = treeMC.getTH1F(lumi, "met_y_InclusiveSF_mc", "met_pt*sin(met_phi)", 50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} y (GeV)")   
  #  met_y_InclusiveSF_dy = treeDY.getTH1F(lumi, "met_y_InclusiveSF_dy", "met_pt*sin(met_phi)",  50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} y (GeV)")   


  #  met_phi_InclusiveSF_data = treeDA.getTH1F(lumi, "met_y_InclusiveSF_data", "met_phi", 40, -3.2, 3.2 , cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} phi")
  #  met_phi_InclusiveSF_mc = treeMC.getTH1F(lumi, "met_y_InclusiveSF_mc", "met_phi",  40, -3.2, 3.2, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} phi")   
  #  met_phi_InclusiveSF_dy = treeDY.getTH1F(lumi, "met_y_InclusiveSF_dy", "met_phi",  40, -3.2, 3.2, cuts.AddList([cuts.InclusiveCRSF()]), "", "ME_{T} phi")   



  #  jzb_InclusiveSF_data = treeDA.getTH1F(lumi, "jzb_InclusiveSF_data", "lepsJZB_Edge",  50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "JZB (GeV)")
  #  jzb_InclusiveSF_mc = treeMC.getTH1F(lumi, "jzb_InclusiveSF_mc", "lepsJZB_Edge",  50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "JZB (GeV)")   
  #  jzb_InclusiveSF_dy = treeDY.getTH1F(lumi, "jzb_InclusiveSF_dy", "lepsJZB_Edge",  50, -260, 260, cuts.AddList([cuts.InclusiveCRSF()]), "", "JZB (GeV)")   



   # met_InclusiveSF_dy.SetFillColor(r.kBlue+1)                     
   # met_InclusiveSF_dy.SetFillColorAlpha(r.kBlue+1, 0.5)            
   # met_InclusiveSF = Canvas.Canvas("metSmearing/met_InclusiveSF", "png", 0.6, 0.6, 0.8, 0.8)                    
   # met_InclusiveSF.addHisto(met_InclusiveSF_mc, "HIST", "MC", "P", r.kRed+1, 1, 0)
   # met_InclusiveSF.addHisto(met_InclusiveSF_dy, "HIST SAME", "DY"  , "PL", r.kBlue+1 , 1, 0) 
   # met_InclusiveSF.addHisto(met_InclusiveSF_data, "E,SAME", "Data"  , "PL", r.kBlack , 1, 1) 
   # met_InclusiveSF.saveRatio(1, 1, 1, lumi, met_InclusiveSF_data, met_InclusiveSF_mc)
  #  del met_InclusiveSF
   # time.sleep(0.1)                                                                       

 #   met_InclusiveOF.SetFillColor(r.kGreen+1)                     
 #   met_InclusiveOF.SetFillColorAlpha(r.kGreen+1, 0.5)            
 #   met_InclusiveEE.SetFillColor(r.kRed+1)               
 #   met_InclusiveEE.SetFillColorAlpha(r.kRed+1, 0.5)     
 #   met_InclusiveMM.SetFillColor(r.kBlue+1)               
 #   met_InclusiveMM.SetFillColorAlpha(r.kBlue+1, 0.5)     


    met_RSFOF = Canvas.Canvas("met_studies/met_RSFOF", "png", 0.6, 0.6, 0.8, 0.8)
    met_RSFOF.addHisto(met_InclusiveOF, "HIST, SAME", "OF", "P", r.kRed+1, 1, 0)
    met_RSFOF.addHisto(met_InclusiveSF, "HIST, SAME", "SF","P", r.kBlue+1, 1, 0)
    met_RSFOF.saveRatio(1,1,1, lumi, met_InclusiveSF, met_InclusiveOF)
  #  del met_InclusiveSF
    time.sleep(0.1)                                                                       

    met_REEOF = Canvas.Canvas("met_studies/met_REEOF", "png", 0.6, 0.6, 0.8, 0.8)
    met_REEOF.addHisto(met_InclusiveOF, "HIST, SAME", "OF", "P", r.kRed+1, 1, 0)
    met_REEOF.addHisto(met_InclusiveEE, "HIST, SAME", "ee","P", r.kBlue+1, 1, 0)
    met_REEOF.saveRatio(1,1,1, lumi, met_InclusiveEE, met_InclusiveOF)
  #  del met_InclusiveSF
    time.sleep(0.1)                                                                    

    met_RMMOF = Canvas.Canvas("met_studies/met_RMMOF", "png", 0.6, 0.6, 0.8, 0.8)
    met_RMMOF.addHisto(met_InclusiveOF, "HIST, SAME", "OF", "P", r.kRed+1, 1, 0)
    met_RMMOF.addHisto(met_InclusiveMM, "HIST, SAME", "mumu","P", r.kBlue+1, 1, 0)
    met_RMMOF.saveRatio(1,1,1, lumi, met_InclusiveMM, met_InclusiveOF)
  #  del met_InclusiveSF
    time.sleep(0.1)                                                                    

   # plot_met_InclusiveEE = Canvas.Canvas("metSmearing/met_InclusiveEE", "png", 0.6, 0.6, 0.8, 0.8)
   # plot_met_InclusiveEE.addHisto(met_InclusiveEE, "HIST", "ee", "P", r.kRed+1 , 1, 0)
   # plot_met_InclusiveEE.save(1,1,0, lumi)
  #  del met_InclusiveSF
   # time.sleep(0.1)

   # plot_met_InclusiveMM = Canvas.Canvas("metSmearing/met_InclusiveMM", "png", 0.6, 0.6, 0.8, 0.8)
   # plot_met_InclusiveMM.addHisto(met_InclusiveMM, "HIST", "mumu", "P", r.kRed+1 , 1, 0)
   # plot_met_InclusiveMM.save(1,1,0, lumi)
  #  del met_InclusiveSF
   # time.sleep(0.1)                                                                       
#   met_raw_InclusiveSF_dy.SetFillColor(r.kBlue+1)                     
 #   met_raw_InclusiveSF_dy.SetFillColorAlpha(r.kBlue+1, 0.5)            
 #   met_raw_InclusiveSF = Canvas.Canvas("metSmearing/met_raw_InclusiveSF", "png", 0.6, 0.6, 0.8, 0.8)
 #   met_raw_InclusiveSF.addHisto(met_raw_InclusiveSF_mc, "HIST", "MC", "P", r.kRed+1, 1, 0)
 #   met_raw_InclusiveSF.addHisto(met_raw_InclusiveSF_dy, "HIST SAME", "DY"  , "PL", r.kBlue+1 , 1, 0) 
 #   met_raw_InclusiveSF.addHisto(met_raw_InclusiveSF_data, "E,SAME", "Data"  , "PL", r.kBlack , 1, 1)
 #   met_raw_InclusiveSF.saveRatio(1, 1, 1, lumi, met_raw_InclusiveSF_data, met_raw_InclusiveSF_mc)
  #  del met_InclusiveSF
 #   time.sleep(0.1)                                                                       


   # met_x_InclusiveSF_dy.SetFillColor(r.kBlue+1)                             
   # met_x_InclusiveSF_dy.SetFillColorAlpha(r.kBlue+1, 0.5)         
   # met_x_InclusiveSF = Canvas.Canvas("metSmearing/met_x_InclusiveSF", "png", 0.6, 0.6, 0.8, 0.8)                 
   # met_x_InclusiveSF.addHisto(met_x_InclusiveSF_mc, "HIST", "MC", "P", r.kRed+1, 1, 0)
   # met_x_InclusiveSF.addHisto(met_x_InclusiveSF_dy, "HIST SAME", "DY"  , "PL", r.kBlue+1 , 1, 0) 
   # met_x_InclusiveSF.addHisto(met_x_InclusiveSF_data, "E,SAME", "Data"  , "PL", r.kBlack , 1, 1) 
   # met_x_InclusiveSF.saveRatio(1, 1, 1, lumi, met_x_InclusiveSF_data, met_x_InclusiveSF_mc)
  #  del met_x_InclusiveSF
   # time.sleep(0.1)                                                                  

   # met_y_InclusiveSF_dy.SetFillColor(r.kBlue+1)                             
   # met_y_InclusiveSF_dy.SetFillColorAlpha(r.kBlue+1, 0.5)         
   # met_y_InclusiveSF = Canvas.Canvas("metSmearing/met_y_InclusiveSF", "png", 0.6, 0.6, 0.8, 0.8)                 
   # met_y_InclusiveSF.addHisto(met_y_InclusiveSF_mc, "HIST", "MC", "P", r.kRed+1, 1, 0)
   # met_y_InclusiveSF.addHisto(met_y_InclusiveSF_dy, "HIST SAME", "DY"  , "PL", r.kBlue+1 , 1, 0) 
   # met_y_InclusiveSF.addHisto(met_y_InclusiveSF_data, "E,SAME", "Data"  , "PL", r.kBlack , 1, 1) 
   # met_y_InclusiveSF.saveRatio(1, 1, 1, lumi, met_y_InclusiveSF_data, met_y_InclusiveSF_mc)
  #  del met_y_InclusiveSF
   # time.sleep(0.1)                                                                  


   # met_phi_InclusiveSF_dy.SetFillColor(r.kBlue+1)                                                
   # met_phi_InclusiveSF_dy.SetFillColorAlpha(r.kBlue+1, 0.5)                                      
   # met_phi_InclusiveSF = Canvas.Canvas("metSmearing/met_phi_InclusiveSF", "png", 0.6, 0.6, 0.8, 0.8)  
   # met_phi_InclusiveSF.addHisto(met_phi_InclusiveSF_mc, "HIST", "MC", "P", r.kRed+1, 1, 0)   
   # met_phi_InclusiveSF.addHisto(met_phi_InclusiveSF_dy, "HIST SAME", "DY"  , "PL", r.kBlue+1 , 1, 0) 
   # met_phi_InclusiveSF.addHisto(met_phi_InclusiveSF_data, "E,SAME", "Data"  , "PL", r.kBlack , 1, 1) 
   # met_phi_InclusiveSF.saveRatio(1, 1, 1, lumi, met_phi_InclusiveSF_data, met_phi_InclusiveSF_mc)     
  #  del met_phi_InclusiveSF
   # time.sleep(0.1)                                                                               



   # jzb_InclusiveSF_dy.SetFillColor(r.kBlue+1)                     
   # jzb_InclusiveSF_dy.SetFillColorAlpha(r.kBlue+1, 0.5)         
   # jzb_InclusiveSF = Canvas.Canvas("metSmearing/jzb_InclusiveSF", "png", 0.6, 0.6, 0.8, 0.8)                 
   # jzb_InclusiveSF.addHisto(jzb_InclusiveSF_mc, "HIST", "MC", "P", r.kRed+1, 1, 0)
   # jzb_InclusiveSF.addHisto(jzb_InclusiveSF_dy, "HIST SAME", "DY"  , "PL", r.kBlue+1 , 1, 0) 
   # jzb_InclusiveSF.addHisto(jzb_InclusiveSF_data, "E,SAME", "Data"  , "PL", r.kBlack , 1, 1) 
   # jzb_InclusiveSF.saveRatio(1, 1, 1, lumi,  jzb_InclusiveSF_data, jzb_InclusiveSF_mc)
  #  del jzb_InclusiveSF
   # time.sleep(0.1)                                                                  






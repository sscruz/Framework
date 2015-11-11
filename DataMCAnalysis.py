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
import math as math

from optparse import OptionParser
from ROOT import gROOT, TCanvas, TFile
from Sample import Sample, Block, Tree
from CutManager import CutManager
from Canvas import Canvas




if __name__ == "__main__":

    parser = OptionParser(usage="usage: %prog [options] FilenameWithSamples", version="%prog 1.0")
    parser.add_option("-m", "--mode", action="store", dest="mode", default="rmue", help="Operation mode")
    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.error("wrong number of arguments")

    inputFileNameMC = args[0]
    inputFileNameData = args[1]
  
    treeMC = Tree(inputFileNameMC, "MC", 0)
    treeData = Tree(inputFileNameData, "Data", 0)
   
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 

    cuts = CutManager()

    lumi = 0.042

    #inclusive region
    mll_SF_central_data = treeData.getTH1F(lumi, "mll_SF_central_MC", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.GoodLeptonSF(), cuts.Central()]), "", "m_{ll} [GeV]")
    mll_SF_central_mc = treeMC.getStack(lumi, "mll_SF_central_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.GoodLeptonSF(), cuts.Central()]), "", "m_{ll} [GeV]")
    mll_SF_central_mc_h = treeMC.getTH1F(lumi, "mll_SF_central_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.GoodLeptonSF(), cuts.Central()]), "", "m_{ll} [GeV]")
    
    mll_OF_central_data = treeData.getTH1F(lumi, "mll_OF_central_MC", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.GoodLeptonOF(), cuts.Central()]), "", "m_{ll} [GeV]")
    mll_OF_central_mc = treeMC.getStack(lumi, "mll_OF_central_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.GoodLeptonOF(), cuts.Central()]), "", "m_{ll} [GeV]")
    mll_OF_central_mc_h = treeMC.getTH1F(lumi, "mll_OF_central_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.GoodLeptonOF(), cuts.Central()]), "", "m_{ll} [GeV]")
 
    plot_SF_central = Canvas("plot_SF_central", "png", 0.6, 0.6, 0.8, 0.8)
    plot_SF_central.addStack(mll_SF_central_mc, "HIST", 1, 1)
    plot_SF_central.addHisto(mll_SF_central_data, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plot_SF_central.saveRatio(1, 1, 1, lumi, mll_SF_central_data, mll_SF_central_mc_h)
    
    plot_OF_central = Canvas("plot_OF_central", "png", 0.6, 0.6, 0.8, 0.8)
    plot_OF_central.addStack(mll_OF_central_mc, "HIST", 1, 1)
    plot_OF_central.addHisto(mll_OF_central_data, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plot_OF_central.saveRatio(1, 1, 1, lumi, mll_OF_central_data, mll_OF_central_mc_h)
   

    mll_SF_forward_data = treeData.getTH1F(lumi, "mll_SF_forward_MC", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.GoodLeptonSF(), cuts.Forward()]), "", "m_{ll} [GeV]")
    mll_SF_forward_mc = treeMC.getStack(lumi, "mll_SF_forward_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.GoodLeptonSF(), cuts.Forward()]), "", "m_{ll} [GeV]")
    mll_SF_forward_mc_h = treeMC.getTH1F(lumi, "mll_SF_forward_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.GoodLeptonSF(), cuts.Forward()]), "", "m_{ll} [GeV]")
    
    mll_OF_forward_data = treeData.getTH1F(lumi, "mll_OF_forward_MC", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.GoodLeptonOF(), cuts.Forward()]), "", "m_{ll} [GeV]")
    mll_OF_forward_mc = treeMC.getStack(lumi, "mll_OF_forward_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.GoodLeptonOF(), cuts.Forward()]), "", "m_{ll} [GeV]")
    mll_OF_forward_mc_h = treeMC.getTH1F(lumi, "mll_OF_forward_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.GoodLeptonOF(), cuts.Forward()]), "", "m_{ll} [GeV]")
 
    plot_SF_forward = Canvas("plot_SF_forward", "png", 0.6, 0.6, 0.8, 0.8)
    plot_SF_forward.addStack(mll_SF_forward_mc, "HIST", 1, 1)
    plot_SF_forward.addHisto(mll_SF_forward_data, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plot_SF_forward.saveRatio(1, 1, 1, lumi, mll_SF_forward_data, mll_SF_forward_mc_h)
    
    plot_OF_forward = Canvas("plot_OF_forward", "png", 0.6, 0.6, 0.8, 0.8)
    plot_OF_forward.addStack(mll_OF_forward_mc, "HIST", 1, 1)
    plot_OF_forward.addHisto(mll_OF_forward_data, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plot_OF_forward.saveRatio(1, 1, 1, lumi, mll_OF_forward_data, mll_OF_forward_mc_h)

    #DY region

    mllDY_SF_central_data = treeData.getTH1F(lumi, "mllDY_SF_central_MC", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.DYControlNoMassLeptonSF(), cuts.Central()]), "", "m_{ll} [GeV]")
    mllDY_SF_central_mc = treeMC.getStack(lumi, "mllDY_SF_central_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.DYControlNoMassLeptonSF(), cuts.Central()]), "", "m_{ll} [GeV]")
    mllDY_SF_central_mc_h = treeMC.getTH1F(lumi, "mllDY_SF_central_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.DYControlNoMassLeptonSF(), cuts.Central()]), "", "m_{ll} [GeV]")
    
    mllDY_OF_central_data = treeData.getTH1F(lumi, "mllDY_OF_central_MC", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.DYControlNoMassLeptonOF(), cuts.Central()]), "", "m_{ll} [GeV]")
    mllDY_OF_central_mc = treeMC.getStack(lumi, "mllDY_OF_central_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.DYControlNoMassLeptonOF(), cuts.Central()]), "", "m_{ll} [GeV]")
    mllDY_OF_central_mc_h = treeMC.getTH1F(lumi, "mllDY_OF_central_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.DYControlNoMassLeptonOF(), cuts.Central()]), "", "m_{ll} [GeV]")
 
    plotDY_SF_central = Canvas("plotDY_SF_central", "png", 0.6, 0.6, 0.8, 0.8)
    plotDY_SF_central.addStack(mllDY_SF_central_mc, "HIST", 1, 1)
    plotDY_SF_central.addHisto(mllDY_SF_central_data, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plotDY_SF_central.saveRatio(1, 1, 1, lumi, mllDY_SF_central_data, mllDY_SF_central_mc_h)
    
    plotDY_OF_central = Canvas("plotDY_OF_central", "png", 0.6, 0.6, 0.8, 0.8)
    plotDY_OF_central.addStack(mllDY_OF_central_mc, "HIST", 1, 1)
    plotDY_OF_central.addHisto(mllDY_OF_central_data, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plotDY_OF_central.saveRatio(1, 1, 1, lumi, mllDY_OF_central_data, mllDY_OF_central_mc_h)
   

    mllDY_SF_forward_data = treeData.getTH1F(lumi, "mllDY_SF_forward_MC", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.DYControlNoMassLeptonSF(), cuts.Forward()]), "", "m_{ll} [GeV]")
    mllDY_SF_forward_mc = treeMC.getStack(lumi, "mllDY_SF_forward_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.DYControlNoMassLeptonSF(), cuts.Forward()]), "", "m_{ll} [GeV]")
    mllDY_SF_forward_mc_h = treeMC.getTH1F(lumi, "mllDY_SF_forward_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.DYControlNoMassLeptonSF(), cuts.Forward()]), "", "m_{ll} [GeV]")
    
    mllDY_OF_forward_data = treeData.getTH1F(lumi, "mllDY_OF_forward_MC", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.DYControlNoMassLeptonOF(), cuts.Forward()]), "", "m_{ll} [GeV]")
    mllDY_OF_forward_mc = treeMC.getStack(lumi, "mllDY_OF_forward_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.DYControlNoMassLeptonOF(), cuts.Forward()]), "", "m_{ll} [GeV]")
    mllDY_OF_forward_mc_h = treeMC.getTH1F(lumi, "mllDY_OF_forward_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.DYControlNoMassLeptonOF(), cuts.Forward()]), "", "m_{ll} [GeV]")
 
    plotDY_SF_forward = Canvas("plotDY_SF_forward", "png", 0.6, 0.6, 0.8, 0.8)
    plotDY_SF_forward.addStack(mllDY_SF_forward_mc, "HIST", 1, 1)
    plotDY_SF_forward.addHisto(mllDY_SF_forward_data, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plotDY_SF_forward.saveRatio(1, 1, 1, lumi, mllDY_SF_forward_data, mllDY_SF_forward_mc_h)
    
    plotDY_OF_forward = Canvas("plotDY_OF_forward", "png", 0.6, 0.6, 0.8, 0.8)
    plotDY_OF_forward.addStack(mllDY_OF_forward_mc, "HIST", 1, 1)
    plotDY_OF_forward.addHisto(mllDY_OF_forward_data, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plotDY_OF_forward.saveRatio(1, 1, 1, lumi, mllDY_OF_forward_data, mllDY_OF_forward_mc_h)
  

 
    #ttbar Control region
 
    mllCR_SF_central_data = treeData.getTH1F(lumi, "mllCR_SF_central_MC", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.ControlNoMassLeptonSF(), cuts.Central()]), "", "m_{ll} [GeV]")
    mllCR_SF_central_mc = treeMC.getStack(lumi, "mllCR_SF_central_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.ControlNoMassLeptonSF(), cuts.Central()]), "", "m_{ll} [GeV]")
    mllCR_SF_central_mc_h = treeMC.getTH1F(lumi, "mllCR_SF_central_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.ControlNoMassLeptonSF(), cuts.Central()]), "", "m_{ll} [GeV]")
    
    mllCR_OF_central_data = treeData.getTH1F(lumi, "mllCR_OF_central_MC", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.ControlNoMassLeptonOF(), cuts.Central()]), "", "m_{ll} [GeV]")
    mllCR_OF_central_mc = treeMC.getStack(lumi, "mllCR_OF_central_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.ControlNoMassLeptonOF(), cuts.Central()]), "", "m_{ll} [GeV]")
    mllCR_OF_central_mc_h = treeMC.getTH1F(lumi, "mllCR_OF_central_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.ControlNoMassLeptonOF(), cuts.Central()]), "", "m_{ll} [GeV]")
 
    plotCR_SF_central = Canvas("plotCR_SF_central", "png", 0.6, 0.6, 0.8, 0.8)
    plotCR_SF_central.addStack(mllCR_SF_central_mc, "HIST", 1, 1)
    plotCR_SF_central.addHisto(mllCR_SF_central_data, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plotCR_SF_central.saveRatio(1, 1, 1, lumi, mllCR_SF_central_data, mllCR_SF_central_mc_h)
    
    plotCR_OF_central = Canvas("plotCR_OF_central", "png", 0.6, 0.6, 0.8, 0.8)
    plotCR_OF_central.addStack(mllCR_OF_central_mc, "HIST", 1, 1)
    plotCR_OF_central.addHisto(mllCR_OF_central_data, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plotCR_OF_central.saveRatio(1, 1, 1, lumi, mllCR_OF_central_data, mllCR_OF_central_mc_h)
   

    mllCR_SF_forward_data = treeData.getTH1F(lumi, "mllCR_SF_forward_MC", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.ControlNoMassLeptonSF(), cuts.Forward()]), "", "m_{ll} [GeV]")
    mllCR_SF_forward_mc = treeMC.getStack(lumi, "mllCR_SF_forward_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.ControlNoMassLeptonSF(), cuts.Forward()]), "", "m_{ll} [GeV]")
    mllCR_SF_forward_mc_h = treeMC.getTH1F(lumi, "mllCR_SF_forward_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.ControlNoMassLeptonSF(), cuts.Forward()]), "", "m_{ll} [GeV]")
    
    mllCR_OF_forward_data = treeData.getTH1F(lumi, "mllCR_OF_forward_MC", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.ControlNoMassLeptonOF(), cuts.Forward()]), "", "m_{ll} [GeV]")
    mllCR_OF_forward_mc = treeMC.getStack(lumi, "mllCR_OF_forward_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.ControlNoMassLeptonOF(), cuts.Forward()]), "", "m_{ll} [GeV]")
    mllCR_OF_forward_mc_h = treeMC.getTH1F(lumi, "mllCR_OF_forward_data", "t.lepsMll_Edge", 40, 20, 300, cuts.AddList([cuts.ControlNoMassLeptonOF(), cuts.Forward()]), "", "m_{ll} [GeV]")
 
    plotCR_SF_forward = Canvas("plotCR_SF_forward", "png", 0.6, 0.6, 0.8, 0.8)
    plotCR_SF_forward.addStack(mllCR_SF_forward_mc, "HIST", 1, 1)
    plotCR_SF_forward.addHisto(mllCR_SF_forward_data, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plotCR_SF_forward.saveRatio(1, 1, 1, lumi, mllCR_SF_forward_data, mllCR_SF_forward_mc_h)
    
    plotCR_OF_forward = Canvas("plotCR_OF_forward", "png", 0.6, 0.6, 0.8, 0.8)
    plotCR_OF_forward.addStack(mllCR_OF_forward_mc, "HIST", 1, 1)
    plotCR_OF_forward.addHisto(mllCR_OF_forward_data, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plotCR_OF_forward.saveRatio(1, 1, 1, lumi, mllCR_OF_forward_data, mllCR_OF_forward_mc_h)
      

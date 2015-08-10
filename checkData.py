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

    if len(args) != 1:
        parser.error("wrong number of arguments")

    inputFileName = args[0]
    tree = Tree(inputFileName, "Data", 0)
   
    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager()

    mll_SF_central_data = tree.getTH1F(4, "mll_SF_central_data", "t.lepsMll_Edge", 40, 20, 300, cuts.Add(cuts.DYControlNoMassLeptonSF(), cuts.Central()), "", "m_{ll} [GeV]")
    
    plot_SFvsOF_central = Canvas("plot_SF_DataMC_central", "png", 0.6, 0.6, 0.8, 0.8)
    plot_SFvsOF_central.addHisto(mll_SF_central_data, "E1", "Data", "P", r.kBlack, 1, 0)
    plot_SFvsOF_central.save(1, 0, 0, 4.0)
    

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
from ROOT import TGraphErrors
from array import array


if __name__ == "__main__":

    parser = OptionParser(usage="usage: %prog [options] FilenameWithSamples", version="%prog 1.0")
    parser.add_option("-m", "--mode", action="store", dest="mode", default="rmue", help="Operation mode")
    (options, args) = parser.parse_args()

    if len(args) != 2:
      parser.error("wrong number of arguments")

    inputFileName = args[0]
    tree23 = Tree(inputFileName, "MC", 0)
    inputFileName2 = args[1]
    tree43 = Tree(inputFileName2, "MC", 0)
   
    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager()

    bins = [50, 60, 70, 81, 101, 120, 150, 180, 220, 260, 300]
    #####Attention: setting the first bin in 50 because we only have DY mll>50 at this point
    ht23 = tree23.getTH1F(4, "ht23", "t.htJet35j_Edge", 20, 0, 500, cuts.Add(cuts.DYControlNoMassLeptonee(), cuts.Central()), "", "H_T [GeV]")
    ht43 = tree43.getTH1F(4, "ht43", "t.htJet35j_Edge", 20, 0, 500, cuts.Add(cuts.DYControlNoMassLeptonee(), cuts.Central()), "", "H_T [GeV]")
    ht23_40 = tree23.getTH1F(4, "ht23_40", "htJet40", 20, 0, 500, cuts.Add(cuts.DYControlNoMassLeptonee(), cuts.Central()), "", "H_T [GeV]")
    ht43_40 = tree43.getTH1F(4, "ht43_40", "htJet40", 20, 0, 500, cuts.Add(cuts.DYControlNoMassLeptonee(), cuts.Central()), "", "H_T [GeV]")
 
    plot_ht_central = Canvas("plot_ht_central", "png", 0.6, 0.6, 0.8, 0.8)
    plot_ht_central.addHisto(ht23, "HIST", "23", "L", r.kRed, 1, 0)
    plot_ht_central.addHisto(ht43, "HIST, SAME", "43", "L", r.kBlue, 1, 0)
    plot_ht_central.saveRatio(1, 0, 0, 4.0, ht23, ht43)
    

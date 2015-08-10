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
    tree = Tree(inputFileName, "MC", 0)
   
    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager()

    mll_SF = tree.getTH1F(4, "mll_SF", "t.lepsMll_Edge", 40, 20, 300, cuts.DYControlNoMassLeptonSF(), "", "m_{ll} (GeV)")
    mll_OF = tree.getTH1F(4, "mll_OF", "t.lepsMll_Edge", 40, 20, 300, cuts.DYControlNoMassLeptonOF(), "", "m_{ll} (GeV)")

    jzb_SF = tree.getTH1F(4, "jzb_SF", "t.lepsJZB_Edge", 40, -100, 100, cuts.DYControlNoMassLeptonSF(), "", "JZB (GeV)")

    plot = Canvas("mll", "png", 0.6, 0.6, 0.8, 0.8)
    plot.addHisto(mll_OF, "HISTO", "OF", "L", r.kBlack)
    plot.addHisto(mll_SF, "E1,SAME", "SF", "P", r.kBlue)
    plot.saveRatio(1, 0, 1, 4.0, 0.8, 1.2)
    
    plot2 = Canvas("jzb", "png,pdf", 0.6, 0.6, 0.8, 0.8)
    plot2.addHisto(jzb_SF, "PE", "SF", "F", r.kBlack)
    plot2.save(0, 0, 0, 4.0)
    












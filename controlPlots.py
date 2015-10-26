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
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership
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
    ##mcDatasets = ['TTLep_pow', 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    ##dyDatasets = ['DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    ##daDatasets = ['DoubleMuon_Run2015C', 'DoubleEG_Run2015C', 'MuonEG_Run2015C',
    ##              'DoubleMuon_Run2015D', 'DoubleEG_Run2015D', 'MuonEG_Run2015D',
    ##              'DoubleMuon_Run2015D_v4', 'DoubleEG_Run2015D_v4', 'MuonEG_Run2015D_v4']
    #mcDatasets = ['TTLep_pow', 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    ttDatasets = ['TTLep_pow']
    dyDatasets = ['DYJetsToLL_M10to50', 'DYJetsToLL_M50']
   # stDatasets = ['TToLeptons_tch_amcatnlo']
    raDatasets = ['WWTo2L2Nu', 'ZZTo2L2Q','WZTo2L2Q']
   # mcDatasets = ['TTLep_pow', 'DYJetsToLL_M10to50', 'DYJetsToLL_M50','WWTo2L2Nu', 'ZZTo2L2Q','WZTo2L2Q']
    daDatasets = ['DoubleMuon_Run2015D_05Oct_v1_runs_246908_258751', 'DoubleEG_Run2015D_05Oct_v1_runs_246908_258751', 'MuonEG_Run2015D_05Oct_v2_runs_246908_258751',
                  'DoubleMuon_Run2015D_v4_runs_246908_258751', 'DoubleEG_Run2015D_v4_runs_246908_258751', 'MuonEG_Run2015D_v4_runs_246908_258751']

    ##treeMC = Sample.Tree(helper.selectSamples(inputFileName, mcDatasets, 'MC'), 'MC'  , 0)
  #  treeMC = Sample.Tree(helper.selectSamples(inputFileName, mcDatasets, 'MC'), 'MC'  , 0)
    treeTT = Sample.Tree(helper.selectSamples(inputFileName, ttDatasets, 'TT'), 'TT'  , 0)
    treeDY = Sample.Tree(helper.selectSamples(inputFileName, dyDatasets, 'DY'), 'DY'  , 0)
    treeRA = Sample.Tree(helper.selectSamples(inputFileName, raDatasets, 'RA'), 'RA'  , 0)
   # treeST = Sample.Tree(helper.selectSamples(inputFileName, stDatasets, 'ST'), 'ST',   0)
    treeDA = Sample.Tree(helper.selectSamples(inputFileName, daDatasets, 'DA'), 'DATA', 1)

   # mcTrees = [treeDY, treeRA, treeTT]
    mcTrees = [treeRA, treeTT, treeRA]
   # mcTrees = [treeTT]
    print 'Trees successfully loaded...'

    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    ####Cuts needed by rmue
    cuts = CutManager.CutManager()

    #lumi = 0.849
    lumi = 1.3
    lumi_str = 'lumi'+str(lumi).replace('.', 'p')


    regions = []
    setLog = []
    Control2JetsSF = Region.region('Control2JetsSF',
                       [cuts.Control2JetsSF()],
                       #['mll', 'met', 'nb', 'nj', 'nvtx'],
                       #[range(20,310,10), range(0,310,10), range(0,5,1), range(0,8,1), range(0,35)],
                       ['mll'], [range(0,310,10)],
                       True)
    regions.append(Control2JetsSF)
    setLog.append(True)

#    Control2JetsOF = Region.region('Control2JetsOF',
#                       [cuts.Control2JetsOF()],
#                       #['mll', 'met', 'nb', 'nj', 'nvtx'],
#                       #[range(10,310,10), range(10,310,10), range(0,5,1), range(0,8,1), range(0,35)],
#                       ['mll'], [range(0,310,10)],
#                       True)
#    regions.append(Control2JetsOF) 
#    setLog.append(False)                      

    for reg in regions:
        print 'i am at region', reg.name
        for eta in ['central', 'forward']:         
            my_cuts = cuts.AddList([cuts.goodLepton, cuts.trigger]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts)
            for tree in ( (mcTrees+[treeDA]) if reg.doData else mcTrees):
           
                if tree == treeDA: 
                    dataMC = 'DATA'
                else: 
                    dataMC = 'MC' 

                for var in reg.rvars:

                    if   var == 'mll':
                        varTitle    = 'm_{ll} (GeV)'
                        varVariable = 't.lepsMll_Edge'
                    elif var == 'met':
                        varTitle    = 'ME_{T} (GeV)'
                        varVariable = 'met_pt'
                    elif var == 'nj':
                        varTitle    = 'n-{jets}'
                        varVariable = 't.nJetSel_Edge'
                    elif var == 'nb':
                        varTitle    = 'n_{b-jets}'
                        varVariable = 't.nBJetMedium35_Edge'
                    elif var == 'nvtx':
                        varTitle    = 'n_{vertices}'
                        varVariable = 'nVert'
                    elif var == 'min_mlb':
                        varTitle    = 'min(m_{lb})'
                        varVariable = 't.min_mlb1_Edge'
                    elif var == 'max_mlb':
                        varTitle    = 'max(m_{lb})'
                        varVariable = 't.min_mlb2_Edge'


                    print 'loading variable %s in %s for %s'%(var, eta, tree.name)

                    attr = var+('' if tree.name in ['DATA', 'MC'] else '_'+tree.name.lower())
                    tmp_full= tree.getTH1F(lumi, var+"_"+eta+reg.name+tree.name, varVariable, reg.bins[reg.rvars.index(var)], 1, 1, my_cuts, "", varTitle)
                    getattr(reg, attr).setHisto(tmp_full, dataMC, eta) 


    for reg in regions:
        for var in reg.rvars:
            for eta in ['central', 'forward']:
                print 'in region ', eta
                stack = r.THStack()
            ## add the MCs up in the MC histo
                mc_hist = 0
                vhists = []
                for tree in mcTrees:
                    treename = tree.name.lower()
                    tmp_hist = getattr(reg, var+'_'+treename).getHisto('MC', eta).Clone(var+eta+reg.name)
                    getattr(reg, var).setHisto(tmp_hist, 'MC', eta)
                    if treename == 'dy':
                        tmp_hist.SetFillColorAlpha(r.kBlue, 0.5)
                    elif treename == 'tt':
                        tmp_hist.SetFillColorAlpha(r.kRed, 0.5)
                    elif treename == 'ra':
                        tmp_hist.SetFillColorAlpha(r.kYellow, 0.5)
                    vhists.append(tmp_hist)
                    if not mc_hist:
                        mc_hist = tmp_hist
                    else:
                        mc_hist.Add(tmp_hist, 1.) 
                    #print 'mc_hist is', mc_hist.GetEntries()
                    #print "tmp_hist is", tmp_hist.GetName(), tmp_hist.GetEntries()
                    #stack.Add(tmp_hist, "")
                    #del tmp_hist
                    del tmp_hist
                getattr(reg, var).setHisto(mc_hist, 'MC', eta)

                stack.Add(vhists[0])
                stack.Add(vhists[1])
                stack.Add(vhists[2])
                # done adding the MC hist
                print 'plotting %s in region %s in %s' %(var, reg.name, eta) 
                
                plot_var = Canvas.Canvas("test/%s/%s_%s_%s"%(lumi_str, var, eta, reg.name), "png,pdf", 0.6, 0.6, 0.8, 0.8)
                SetOwnership(stack, False )
                #plot_var.addHisto(getattr(reg, var).getHisto('DATA', eta)  , "E"     , "Data", "PL", r.kBlack , 1, 1)
                plot_var.addStack(stack, "HIST", 1, 1)  
                plot_var.addHisto(getattr(reg, var).getHisto('DATA', eta), "E, SAME"   , "Data"  , "PL", r.kBlack , 1, 1)   
                print 'data is (getEntries)', getattr(reg, var).getHisto('DATA', eta).GetEntries()
                print 'data is (integral)', getattr(reg, var).getHisto('DATA', eta).Integral()
                print 'total mc is (getEntries)', getattr(reg, var).getHisto('MC', eta).GetEntries()
                print 'total mc is (integral)', getattr(reg, var).getHisto('MC', eta).Integral()
                plot_var.saveRatio(1, 1, 1, lumi, getattr(reg, var).getHisto('DATA', eta), getattr(reg, var).getHisto('MC', eta))
                del plot_var  
                del stack
                time.sleep(0.1)









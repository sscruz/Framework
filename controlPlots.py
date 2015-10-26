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
import math, sys, optparse, array, time, copy
import Rounder as rounder

import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample




if __name__ == "__main__":


    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option("-m", "--mode", action="store", dest="mode", default="rmue", help="Operation mode")
    parser.add_option("-t", "--test", action="store_true", dest="test", default=False, help="just do a testrun. takes one variable in one eta for one region")
    (opts, args) = parser.parse_args()

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
    raDatasets = ['WWTo2L2Nu', 'ZZTo2L2Q', 'TToLeptons_tch_amcatnlo']#'WZTo2L2Q']
    daDatasets = ['DoubleMuon_Run2015D_05Oct_v1_runs_246908_258751', 'DoubleEG_Run2015D_05Oct_v1_runs_246908_258751', 'MuonEG_Run2015D_05Oct_v2_runs_246908_258751',
                  'DoubleMuon_Run2015D_v4_runs_246908_258751', 'DoubleEG_Run2015D_v4_runs_246908_258751', 'MuonEG_Run2015D_v4_runs_246908_258751']

    ##treeMC = Sample.Tree(helper.selectSamples(inputFileName, mcDatasets, 'MC'), 'MC'  , 0)
    treeTT = Sample.Tree(helper.selectSamples(inputFileName, ttDatasets, 'TT'), 'TT'  , 0)
    treeDY = Sample.Tree(helper.selectSamples(inputFileName, dyDatasets, 'DY'), 'DY'  , 0)
    treeRA = Sample.Tree(helper.selectSamples(inputFileName, raDatasets, 'RA'), 'RA'  , 0)

    treeDA = Sample.Tree(helper.selectSamples(inputFileName, daDatasets, 'DA'), 'DATA', 1)

    mcTrees = [treeRA, treeDY, treeTT]
    #tree = treeMC
    print 'Trees successfully loaded...'

    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    ####Cuts needed by rmue
    cuts = CutManager.CutManager()

    color = helper.color

    #lumi = 0.849
    #lumi = 1.28
    lumi = 1.3
    lumi_str = 'lumi'+str(lumi).replace('.', 'p')+'_blind'


    regions = []
    setLog = []

    doVariables = ['mll', 'met', 'nb', 'nj', 'nvtx', 'mlb', 'l1pt', 'l2pt', 'ht']
    binnings    = [range(20,310,10), range(0,165, 5), range(0,5,1), range(0,8,1), range(0,35), range(0,310,10), range(20, 205, 5), range(20, 155, 5), range(40, 520, 20)]

    if opts.test:
        doVariables = [doVariables[1]]
        binnings    = [binnings[1]]
    doVariables =doVariables[-4:]
    binnings =binnings[-4:]

    Control2JetsAF_blind = Region.region('Control2JetsAF_blind',
                       [cuts.Control2JetsAF(), cuts.blinded], doVariables, binnings, True)
    regions.append(Control2JetsAF_blind) 

    ##Control0JetsSF = Region.region('Control0JetsSF',
    ##                   [cuts.nj0, cuts.SF], doVariables, binnings, True)
    ##regions.append(Control0JetsSF)

    ##Control0JetsOF = Region.region('Control0JetsOF',
    ##                   [cuts.nj0, cuts.OF], doVariables, binnings, True)
    ##regions.append(Control0JetsOF)
    ##
    ##Control2JetsSF = Region.region('Control2JetsSF',
    ##                   [cuts.Control2JetsSF()], doVariables, binnings, True)
    ##regions.append(Control2JetsSF)
    ##
    ##Control2JetsOF = Region.region('Control2JetsOF',
    ##                   [cuts.Control2JetsOF()], doVariables, binnings, True)
    ##regions.append(Control2JetsOF) 


    etas = ['central' , 'forward']
    if opts.test:
        regions = regions[:1]
        etas = etas[:1]
        

    for reg in regions:
        print color.bold+color.red+'=========================================='
        print 'i am at region', reg.name
        print '=========================================='+color.end
        #for eta in ['central', 'forward']:
        for eta in etas:
                    

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
                elif var == 'mlb':
                    varTitle    = 'min(m_{lb})'
                    varVariable = 't.min_mlb1_Edge'
                elif var == 'l1pt':
                    varTitle    = 'p_{T, l1}'
                    varVariable = 't.Lep1_pt_Edge'
                elif var == 'l2pt':
                    varTitle    = 'p_{T, l2}'
                    varVariable = 't.Lep2_pt_Edge'
                elif var == 'ht':
                    varTitle    = 'H_{T}'
                    varVariable = 't.htJet35j_Edge'

                print 'loading variable %s in %s'%(var, eta)
                mc_histo = 0
                mc_stack = r.THStack()

                for tree in ( (mcTrees+[treeDA]) if reg.doData else mcTrees):
                    ind = 0
                    treename = tree.name.lower()
                    print '... doing tree %s' %treename

                    my_cuts = cuts.AddList([cuts.GoodLepton()]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts)

                    block = tree.blocks[0]
            
                    dataMC = ('DATA' if tree == treeDA else 'MC'); isData = dataMC == 'DATA';

                    attr = var+('' if tree.name in ['DATA', 'MC'] else '_'+treename)
                    tmp_full= tree.getTH1F(lumi, var+"_"+eta+reg.name+treename, varVariable, reg.bins[reg.rvars.index(var)], 1, 1, my_cuts, "", varTitle)
                    if tree == treeTT:
                        tmp_full.SetLineColor(r.kRed+1)
                        tmp_full.SetFillColor(block.color)
                    else:
                        tmp_full.SetFillColorAlpha(block.color, 0.5)
                    tmp_full.SetTitle(block.name)

                    getattr(reg, attr).setHisto(tmp_full, dataMC, eta)

                    ## don't do any adding for data
                    if isData: 
                        continue

                    tmp_histo = copy.deepcopy(tmp_full.Clone(var+eta+reg.name))
                    mc_stack.Add(tmp_histo)
                    if not ind: mc_stack.Draw()
                    ind+=1
                    mc_stack.GetXaxis().SetTitle(tmp_histo.GetXaxis().GetTitle())
                    print 'tmp_histo has integral %.2f'%tmp_histo.Integral()
                    if not mc_histo:
                        mc_histo = copy.deepcopy(tmp_histo)
                    else:
                        mc_histo.Add(tmp_histo, 1.)
                setattr(reg, '%s_mc_histo_%s'%(var, eta), mc_histo)
                setattr(reg, '%s_mc_stack_%s'%(var, eta), mc_stack)

                print 'plotting %s in region %s in %s' %(var, reg.name, eta)
                plot_var = Canvas.Canvas("test/%s/%s_%s_%s"%(lumi_str, var, eta, reg.name), "png,pdf", 0.6, 0.7, 0.8, 0.9)
                plot_var.addStack(mc_stack  , "hist" , 1, 1)
                plot_var.addHisto(getattr(reg, var).getHisto('DATA', eta), "E,SAME"   , "Data"  , "PL", r.kBlack , 1, 0)
                plot_var.saveRatio(1, 1, 1, lumi, getattr(reg, var).getHisto('DATA', eta), mc_histo )
                time.sleep(0.1)


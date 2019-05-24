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
import math, sys, optparse, array, copy
import gc, inspect
from prettytable import PrettyTable

import include.LeptonSF
import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
from itertools import product
from include.tableFormats import latexTable


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def makePlot(lumi, lumi_str, treeDA, treeMC, var, name, nbin, xmin, xmax, theCut, cuts, labelx, logx, protection):

    print 'getting stack'
    MCS  = treeMC.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, theCut, "", labelx,'1',0)
    print' getting mc'
    DATA = treeDA.getTH1F(lumi, "hDATA_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx,'1',0)
    MC = MCS.GetStack().Last()


    maxValmc = MC.GetMaximum() #GetBinContent(MC.GetMaximumBin())
    maxValdata = DATA.GetMaximum() #GetBinContent(DATA.GetMaximumBin())
    maxVal = max(maxValmc, maxValdata)
    if not logx:
        MC.GetYaxis().SetRangeUser(0.1, 1.3*maxVal)
        MCS.SetMaximum(1.3*maxVal)
        DATA.GetYaxis().SetRangeUser(0.1, 1.3*maxVal)
    else:
        MC.GetYaxis().SetRangeUser(0.1, 2.0*maxVal)
        MCS.SetMaximum(2.0*maxVal)
        DATA.GetYaxis().SetRangeUser(0.1, 2.0*maxVal)
    print "DATA integral ", DATA.Integral()
    print "MC   integral ", MC.Integral()
    SetOwnership(MC, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(MCS, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(DATA, 0 )   # 0 = release (not keep), 1 = kee
    plot = Canvas.Canvas('DataMC/%s/plot_%s'%(lumi_str,name), 'png,pdf,root', 0.65, 0.6, 0.85, 0.9)
    plot.addStack(MCS, "HIST", 1, 1)
    DATA.SetMarkerStyle(r.kFullCircle)
    plot.addHisto(DATA, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plot.saveRatio(1, 1, logx, lumi, DATA, MC, r_ymin=0.6, r_ymax=1.4)


def makeTable(lumi, treeDA, treeMC, theCut, saveName=''):
    print 'mctable'
    mctable   = treeMC.getYieldTable( lumi, theCut, "1", '1')
    print 'data table'
    datatable = treeDA.getYieldTable( lumi, theCut, "1", '1')
    print 'finish'
    t = PrettyTable(['Process', 'Yield'])
    total = 0; total_e = 0;
    for key, val in mctable.iteritems():
        total   += val[0]
        total_e += val[1]**2
        t.add_row([key, '%4.2f +/- %4.2f'%(val[0],val[1])])
    t.add_row(['Total mc', '%4.2f +/- %4.2f'%(total,math.sqrt(total_e))])
    for key, val in datatable.iteritems():
        t.add_row(['Data', '%d'%int(val[0])])
    print t 
    if saveName != "":
        out = open('plots/DataMC/yields_%s.txt'%saveName,'w')
        xt = latexTable(t)
        out.write(xt.get_string())
        out.close()

    return

if __name__ == "__main__":

    print bcolors.HEADER
    print '#######################################################################'
    print '                  Starting r_SFOF analysis...                          '
    print '#######################################################################' + bcolors.ENDC

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    parser.add_option('-y', '--year', action='store', type=int, dest='year', default=0, help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()

    print bcolors.HEADER + '[Data - MC comparisons] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    #DYDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    #ttDatasets = ['TTJets_DiLepton']
    #mcDatasets = ['TTJets_SingleLeptonFromT',  'T_tch_powheg', 'TBar_tch_powheg', 'WWTo2L2Nu',  'WZTo3LNu', 'ZZTo4L', 'ZZTo2L2Nu', 'WWW', 'WZZ', 'ZZZ', 'TWZ', 'tZq_ll', 'TTZToLLNuNu' , 'TTLLJets_m1to10', 'TTWToLNu', 'TTTT', 'TTHnobb_pow', 'GGHZZ4L',  'WJetsToLNu_LO']
    #mcDatasets += ttDatasets
    #mcDatasets += DYDatasets

    dyDatasets = ['DYJetsToLL_M50_part1+DYJetsToLL_M50_part2','DYJetsToLL_M10to50']
    ttDatasets = ['TTTo2L2Nu_part1+TTTo2L2Nu_part2','TTJets_SingleLeptonFromT+TTJets_SingleLeptonFromT_ext1','TTJets_SingleLeptonFromTbar+TTJets_SingleLeptonFromTbar_ext1'] 
    stDatasets = ['TW','TbarW'] # t and s (lol) channel missing 
    ttzDatasets = ['TTZToLLNuNu_ext+TTZToLLNuNu_ext3','TTZToLLNuNu_m1to10' ,'TTWToLNu','TTWW']#,'TTWZ','TTGJets_newpmx']
    zz2lDatasets = ['ZZTo2L2Nu']
    zz4lDatasets = ['ZZTo4L',
                    'GluGluToContinToZZTo2e2mu',
                    'GluGluToContinToZZTo2e2nu',
                    'GluGluToContinToZZTo2mu2nu',
                    'GluGluToContinToZZTo4e',
                    'GluGluToContinToZZTo4mu'
    ]
    wwDatasets = ['WWTo2L2Nu']
    wzDatasets = ['WZTo3LNu','WZTo2L2Q']
    raDatasets = [ 'TTHnobb']
    mcDatasets = zz4lDatasets + zz2lDatasets + ttzDatasets + raDatasets + wwDatasets +wzDatasets + stDatasets+  ttDatasets + dyDatasets
    

    
    if opts.year == 2016:
        lumi = 35.9
        lumi_str = '35.9fb'
        pds  = 'DoubleEG,DoubleMuon,SingleElectron,SingleMuon,MuonEG'.split(',')
        eras = 'C,D,E,F,G,H'.split(',')
        mcDatasets = helper.replaceSamp(mcDatasets,'DYJetsToLL_M10to50','DYJetsToLL_M10to50_LO')
    elif opts.year == 2017:
        lumi = 41.4
        lumi_str = '41.1fb'
        pds  = 'DoubleEG,DoubleMuon,SingleElectron,SingleMuon,MuonEG'.split(',')
        eras = 'B,C,D,E,F'.split(',')
        mcDatasets = helper.replaceSamp(mcDatasets,'TTZToLLNuNu_ext+TTZToLLNuNu_ext3','TTZToLLNuNu')
        mcDatasets = helper.replaceSamp(mcDatasets,'TTJets_SingleLeptonFromT+TTJets_SingleLeptonFromT_ext1','TTToSemiLeptonic+TTToSemiLeptonic_PSweights')
        mcDatasets = helper.removeSamp(mcDatasets,'TTJets_SingleLeptonFromTbar+TTJets_SingleLeptonFromTbar_ext1')
        mcDatasets = helper.replaceSamp(mcDatasets,'DYJetsToLL_M50_part1+DYJetsToLL_M50_part2','DYJetsToLL_M50_ext_part1+DYJetsToLL_M50_ext_part2+DYJetsToLL_M50_ext_part3+DYJetsToLL_M50_part1+DYJetsToLL_M50_part2+DYJetsToLL_M50_part3')
        mcDatasets = helper.removeSamp(mcDatasets,'TTHnobb')

    elif opts.year == 2018:
        lumi = 58.8
        lumi_str = '58.8fb'
        pds  = 'DoubleMuon,EGamma,SingleMuon,MuonEG'.split(',')
        eras = 'A,B,C,D'.split(',')
        mcDatasets = helper.replaceSamp(mcDatasets,'TTZToLLNuNu_ext+TTZToLLNuNu_ext3','TTZToLLNuNu')
        mcDatasets = helper.replaceSamp(mcDatasets,'TTJets_SingleLeptonFromT+TTJets_SingleLeptonFromT_ext1','TTToSemiLeptonic')
        mcDatasets = helper.removeSamp(mcDatasets,'TTJets_SingleLeptonFromTbar+TTJets_SingleLeptonFromTbar_ext1')
        mcDatasets = helper.replaceSamp(mcDatasets,'DYJetsToLL_M50_part1+DYJetsToLL_M50_part2','DYJetsToLL_M50_LO_part1+DYJetsToLL_M50_LO_part2')
        mcDatasets = helper.removeSamp(mcDatasets,'TTZToLLNuNu_m1to10')
        mcDatasets = helper.removeSamp(mcDatasets,'TTHnobb')

    daDatasets = [ '%s_Run%d%s_%d'%(x,opts.year,y,opts.year) for x,y in product(pds,eras) ]
    mcDatasets = [ '+'.join([(co.replace('_part','_%d_part'%opts.year) if '_part' in co else co+'_%d'%opts.year) for co in da.split('+') ]) for da in mcDatasets ]  # easy to debug
    

    treeMC  = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0, isScan = 0)
    treeDA  = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1, isScan = 0)
    
    
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()
    labelmll = 'm_{ll} [GeV]'
    labelmet = 'E_{T}^{miss} [GeV]'
    labelmt2 = 'MT2 [GeV]'
    labelpt2 = 'p_{T} (2nd) [GeV]'
    labelpt1 = 'p_{T} (1nd) [GeV]'
    labelnjet = "N. Jets"
    labelnbjet = "N. B Jets"
    labelnbjetnjet = "(N_{jets},N_{b-jets})"
    labelfatjet = "N. AK8 jets" 




    makePlot(lumi, lumi_str, treeDA, treeMC, "MET_pt_Edge", "MET_DYJets_OF", 20, 0, 500, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.diLeptonPt]), cuts, labelmet,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "MET_pt_Edge", "MET_DYJets_ee", 20, 0, 500, cuts.AddList([cuts.ee, cuts.goodLepton, cuts.diLeptonPt]), cuts, labelmet,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "MET_pt_Edge", "MET_DYJets_mm", 20, 0, 500, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.diLeptonPt]), cuts, labelmet,  1, 0)

    makePlot(lumi, lumi_str, treeDA, treeMC, "Lep1_pt_Edge", "lep1_pt_DYJets_OF", 20, 0, 500, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.diLeptonPt]), cuts, labelpt1,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "Lep1_pt_Edge", "lep1_pt_DYJets_ee", 20, 0, 500, cuts.AddList([cuts.ee, cuts.goodLepton, cuts.diLeptonPt]), cuts, labelpt1,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "Lep1_pt_Edge", "lep1_pt_DYJets_mm", 20, 0, 500, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.diLeptonPt]), cuts, labelpt1,  1, 0)

    makePlot(lumi, lumi_str, treeDA, treeMC, "Lep2_pt_Edge", "lep2_pt_DYJets_OF", 20, 0, 500, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.diLeptonPt]), cuts, labelpt2,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "Lep2_pt_Edge", "lep2_pt_DYJets_ee", 20, 0, 500, cuts.AddList([cuts.ee, cuts.goodLepton, cuts.diLeptonPt]), cuts, labelpt2,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "Lep2_pt_Edge", "lep2_pt_DYJets_mm", 20, 0, 500, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.diLeptonPt]), cuts, labelpt2,  1, 0)

    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_DYJets_OF", 20, 0, 500, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.diLeptonPt]), cuts, labelnjet,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_DYJets_ee", 20, 0, 500, cuts.AddList([cuts.ee, cuts.goodLepton, cuts.diLeptonPt]), cuts, labelnjet,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_DYJets_mm", 20, 0, 500, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.diLeptonPt]), cuts, labelnjet,  1, 0)




#####################################################################
######                                                              #
###### 8=========D         OO                        8==========D   #  
###### OO                  ||                                ,88   #
###### ||                  ||                              ,88'    #  
###### ||O---\     ,adPPYb,|| ,adPPYb,d8  ,adPPYba,      ,88'      #
###### ||O---/    a8'    `Y||a8'    `Y88 a8P_____88    ,88'        #
###### ||         8b       ||8b       88 8PP'''''''  ,88'          #
###### \/         '8a,   ,d||'8a,   ,d88 '8b,   ,aa 88'            #
###### 8=========D `'8bbdP'\/  `'YbbdP'Y8  `'Ybbd8'' 8==========D   #
######                       aa,    ,88                             #
######                         'Y8bbdP'                             #
######                                                              #
#####################################################################

import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats
import math, sys, optparse, copy, re, array


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables

class onZResult:
    def __init__(self, infile):
        self.infile = infile
        self.readOnZResults()
        self.makeFinalResult()

    def readOnZResults(self):
        f = open(self.infile, 'r')
        for line in f.readlines():
            if '#' in line: continue
            if not len(line.rstrip('\n')): continue
            setattr(self, line.split()[0]       , float(line.split()[1]))
            setattr(self, line.split()[0]+'_err', float(line.split()[2]))
        f.close()

    def makeFinalResult(self):
        self.cen_incb   = self.cen_incbtag_2jet + self.cen_incbtag_3jet
        self.fwd_incb   = self.fwd_incbtag_2jet + self.fwd_incbtag_3jet
        self.cen_incb_e = math.sqrt(self.cen_incbtag_2jet_err**2 + self.cen_incbtag_3jet_err**2)
        self.fwd_incb_e = math.sqrt(self.fwd_incbtag_2jet_err**2 + self.fwd_incbtag_3jet_err**2)

        self.cen_0b   = self.cen_0btag_2jet + self.cen_0btag_3jet
        self.fwd_0b   = self.fwd_0btag_2jet + self.fwd_0btag_3jet
        self.cen_0b_e = math.sqrt(self.cen_0btag_2jet_err**2 + self.cen_0btag_3jet_err**2)
        self.fwd_0b_e = math.sqrt(self.fwd_0btag_2jet_err**2 + self.fwd_0btag_3jet_err**2)

        self.cen_1b   = self.cen_1btag_2jet + self.cen_1btag_3jet
        self.fwd_1b   = self.fwd_1btag_2jet + self.fwd_1btag_3jet
        self.cen_1b_e = math.sqrt(self.cen_1btag_2jet_err**2 + self.cen_1btag_3jet_err**2)
        self.fwd_1b_e = math.sqrt(self.fwd_1btag_2jet_err**2 + self.fwd_1btag_3jet_err**2)

        self.cen_2b   = self.cen_2btag_2jet + self.cen_2btag_3jet
        self.fwd_2b   = self.fwd_2btag_2jet + self.fwd_2btag_3jet
        self.cen_2b_e = math.sqrt(self.cen_2btag_2jet_err**2 + self.cen_2btag_3jet_err**2)
        self.fwd_2b_e = math.sqrt(self.fwd_2btag_2jet_err**2 + self.fwd_2btag_3jet_err**2)


def makeResultsTable(binnedSR, dyShapes, dataMC, eta, nbs):
    line0 = ' %30s &           & '  %('')
    line1 = ' %30s & OF pred.  & ' %('\multirow{3}{*}{%s}' %(binnedSR.name))
    line2 = ' %30s & DY pred.  & ' %('')
    line3 = ' %30s & total     & ' %('')
    line4 = ' %30s & obs.      & ' %('')
    tmp_histo_obs  = binnedSR.mll     .getHisto(dataMC, eta)
    tmp_histo_pred = binnedSR.mll_pred.getHisto(dataMC, eta)
    tmp_histo_dy   = dyShapes['%db_%s_%s_binned'%(nbs, dataMC[:2].lower(), eta)]

    my_range = range(1,tmp_histo_obs.GetNbinsX()+1)
    for i in my_range:
        tmp_dy   = tmp_histo_dy.GetBinContent(i)  ; tmp_dy_e   = tmp_histo_dy.GetBinError(i)
        tmp_of   = tmp_histo_pred.GetBinContent(i); tmp_of_e   = tmp_histo_pred.GetBinError(i)
        tmp_full = tmp_dy + tmp_of                ; tmp_full_e = math.sqrt(tmp_dy_e**2 + tmp_of_e**2)
        tmp_obs  = tmp_histo_obs.GetBinContent(i) ; tmp_obs_e  = tmp_histo_obs.GetBinError(i)
        mll_low , mll_high = tmp_histo_pred.GetXaxis().GetBinLowEdge(i), tmp_histo_pred.GetXaxis().GetBinUpEdge(i)

        line0 += '%.0f $<$ \\mll $<$ %.0f %s' %(mll_low, mll_high   , ' & ' if i != max(my_range) else '\\\\')
        line1 += '  %.2f $\\pm$ %.2f      %s' %(tmp_of  , tmp_of_e  , ' & ' if i != max(my_range) else '\\\\')
        line2 += '  %.2f $\\pm$ %.2f      %s' %(tmp_dy  , tmp_dy_e  , ' & ' if i != max(my_range) else '\\\\')
        line3 += '  %.2f $\\pm$ %.2f      %s' %(tmp_full, tmp_full_e, ' & ' if i != max(my_range) else '\\\\')
        line4 += '  %.2f $\\pm$ %.2f      %s' %(tmp_obs , tmp_obs_e , ' & ' if i != max(my_range) else '\\\\')
    line0 += '\\hline'; line2 += '\\hline'; line3 += '\\hline \\hline'

    return line0, line1, line2, line3, line4

    

def makeRatioTable(myregion, dataMC, eta, nbs): ## for this to make sense the region should be properly binned!!
    header= 'THIS IS THE TABLE FOR %s in %s for %s b-tags'%(dataMC, eta, str(nbs))
    line0 = ' %30s &           & '  %('')
    line1 = ' %30s & %s pred.  & ' %('\multirow{3}{*}{%s}' %(myregion.name), dataMC)
    line2 = ' %30s & %s obs.   & ' %('', dataMC)
    line3 = ' %30s &  ratio    & ' %('')
    tmp_histo_obs  = myregion.mll     .getHisto(dataMC, eta)
    tmp_histo_pred = myregion.mll_pred.getHisto(dataMC, eta)
    my_range = range(1,tmp_histo_obs.GetNbinsX()+1)
    for i in my_range:
        tmp_ratio   = tmp_histo_obs.GetBinContent(i) / tmp_histo_pred.GetBinContent(i)
        tmp_ratio_e = math.sqrt( (tmp_histo_obs.GetBinError(i)/tmp_histo_obs.GetBinContent(i))**2  +  (tmp_histo_pred.GetBinError(i)/tmp_histo_obs.GetBinContent(i))**2)  * tmp_ratio
        line0 += '%.0f $<$ \\mll $<$ %.0f ' %(tmp_histo_pred.GetXaxis().GetBinLowEdge(i), tmp_histo_pred.GetXaxis().GetBinUpEdge(i))
        line1 += '  %.2f $\\pm$ %.2f      ' %(tmp_histo_pred.GetBinContent(i), tmp_histo_pred.GetBinError(i))
        line2 += '  %.2f $\\pm$ %.2f      ' %(tmp_histo_obs.GetBinContent(i), tmp_histo_obs.GetBinError(i))
        line3 += '  %.2f $\\pm$ %.2f      ' %(tmp_ratio, tmp_ratio_e)
        if i != max(my_range):
            line0+=' & '
            line1+=' & '
            line2+=' & '
            line3+=' & '
        else:
            line0+=' \\\\ '
            line1+=' \\\\ '
            line2+=' \\\\ '
            line3+=' \\\\ '
    return header, line0, line1, line2, line3

def scaleZ(histo, eta, nbs):
    bin81  = histo.FindBin(81)
    bin101 = histo.FindBin(100.5)
    nb_pred = getattr(onZ, '%s_%sb'%( ('cen' if eta == 'central' else 'fwd'), str(nbs)))
    ret_hist = copy.deepcopy(histo)
    integral = ret_hist.Integral(bin81, bin101)
    if not integral:
        print 'WARNING: the DY shape histogram has no content around the Z-mass !!!'
        print 'WARNING: this happens for nb %s' %(str(nbs))
        print 'WARNING: putting the whole prediction into one bin at the Z-mass'
        ret_hist.SetBinContent(ret_hist.FindBin(91.), nb_pred)
    else:
        scale  = nb_pred / histo.Integral(bin81, bin101)
        ret_hist.Scale(lumi/2.1*scale)
    return ret_hist

def makePrediction(of_histo, ing, eta):
    central = (eta == 'central')
    sf_pred = copy.deepcopy(of_histo)
    scale = ing.rsfof_final_cen_val  if central else ing.rsfof_final_fwd_val
    s_err = ing.rsfof_final_cen_err  if central else ing.rsfof_final_fwd_err
    print '========================================='
    print 'using rsfof for %s in %s: %.4f +- %.4f' %(ing.dType, eta, scale, s_err)
    print '========================================='
    for _bin in range(1,of_histo.GetNbinsX()+1):
        tmp_cont  = of_histo.GetBinContent(_bin)
        tmp_err   = of_histo.GetBinError  (_bin)
        tmp_mass  = of_histo.GetXaxis().GetBinCenter(_bin)
        tmp_pred  = scale*tmp_cont
        tmp_prede = math.sqrt(tmp_err*tmp_err + s_err*s_err*tmp_cont*tmp_cont)
        sf_pred.SetBinContent(_bin, tmp_pred )
        sf_pred.SetBinError  (_bin, tmp_prede)
    return sf_pred

def getDYShapes(binsPlot, binsFine, binsSR):
    dy_reg = Region.region('dy_region_shortBin' , [cuts.DYControlRegion], ['mll'], [range(20,310,1)], True)
    dy_dict = {}
    ## open the file with the dy shapes in READ or RECREATE depending on if you want to reload it
    dy_shapes_file = r.TFile('dy_shapes.root', 'READ')
    dy_shapes_file.cd()

    for eta in ['central', 'forward']:
        fileHasShapes = bool((dy_shapes_file.Get('dy_shape_da_%s'%(eta)) and dy_shapes_file.Get('dy_shape_mc_%s'%(eta))))
        print 'the DY-shape file already has the shapes for %s : %s' %(eta, str(fileHasShapes))
        if opts.loadShapes or not fileHasShapes:
            ## ========================================
            ## get the dy shapes from the rinout region
            ## ========================================
            print 're-openining the file with UPDATE'
            dy_shapes_file.ReOpen('UPDATE')
            etacut = cuts.Central() if eta == 'central' else cuts.Forward()
            thecutsSF_nomass = cuts.AddList([cuts.GoodLeptonSF()] + [etacut] + dy_reg.cuts)
            thecutsOF_nomass = cuts.AddList([cuts.GoodLeptonOF()] + [etacut] + dy_reg.cuts)
            print 'getting SF data'
            setattr(dy_reg, 'mll_SF_da_'+eta, treeDA.getTH1F(lumi, "mll_SF_in_" +eta+dy_reg.name+'_data' , "t.lepsMll_Edge", dy_reg.bins[0], 1, 1, thecutsSF_nomass, "", "m_{ll} (GeV)"))
            print 'getting OF data'
            setattr(dy_reg, 'mll_OF_da_'+eta, treeDA.getTH1F(lumi, "mll_OF_in_" +eta+dy_reg.name+'_data' , "t.lepsMll_Edge", dy_reg.bins[0], 1, 1, thecutsOF_nomass, "", "m_{ll} (GeV)"))
            print 'getting SF mc'
            setattr(dy_reg, 'mll_SF_mc_'+eta, treeMC.getTH1F(lumi, "mll_SF_in_" +eta+dy_reg.name+'_mc'   , "t.lepsMll_Edge", dy_reg.bins[0], 1, 1, thecutsSF_nomass, "", "m_{ll} (GeV)"))
            print 'getting OF mc'
            setattr(dy_reg, 'mll_OF_mc_'+eta, treeMC.getTH1F(lumi, "mll_OF_in_" +eta+dy_reg.name+'_mc'   , "t.lepsMll_Edge", dy_reg.bins[0], 1, 1, thecutsOF_nomass, "", "m_{ll} (GeV)"))

            # take the SF histogram of the finely binned distribution
            setattr(dy_reg, 'dy_shape_da_'+eta, getattr(dy_reg, 'mll_SF_da_'+eta).Clone('dy_shape_da_%s'%(eta)))
            setattr(dy_reg, 'dy_shape_mc_'+eta, getattr(dy_reg, 'mll_SF_mc_'+eta).Clone('dy_shape_mc_%s'%(eta)))
            # subtract OF from the SF yield
            getattr(dy_reg, 'dy_shape_da_'+eta).Add(getattr(dy_reg, 'mll_OF_da_'+eta), -1.)
            getattr(dy_reg, 'dy_shape_mc_'+eta).Add(getattr(dy_reg, 'mll_OF_mc_'+eta), -1.)
            maxbin = getattr(dy_reg, 'dy_shape_mc_'+eta).GetXaxis().GetXmax() #axis().GetBinUpEdge(getattr(dy_reg, 'dy_shape_mc_'+eta).GetXaxis().GetLast())
            print 'maxbin is', maxbin

            print  'we now have the shapes with bin-width = 1. now let\'s rebin them'
            # NOTE THAT THE SHAPES ARE **NOT** NORMALIZED TO ANYTHING USEFUL!!
            dy_dict['da_%s'%(eta)] = getattr(dy_reg, 'dy_shape_da_%s'%eta)
            dy_dict['mc_%s'%(eta)] = getattr(dy_reg, 'dy_shape_mc_%s'%eta)
            ## print 'DA xmin  %.2f, xmax  %.2f, nbins  %.2f ' %(dy_dict['da_%s'%(eta)].GetXaxis().GetXmin(), dy_dict['da_%s'%(eta)].GetXaxis().GetXmax(), dy_dict['da_%s'%(eta)].GetNbinsX())
            ## print 'MC xmin  %.2f, xmax  %.2f, nbins  %.2f ' %(dy_dict['mc_%s'%(eta)].GetXaxis().GetXmin(), dy_dict['mc_%s'%(eta)].GetXaxis().GetXmax(), dy_dict['mc_%s'%(eta)].GetNbinsX())

            # binning relevant for the plots
            dy_dict['da_%s_binsPlot'%(eta)] = dy_dict['da_%s'%(eta)].Rebin(len(binsPlot), dy_dict['da_%s'%(eta)].GetName()+'_binsPlot', array.array('d', binsPlot+[maxbin]))
            dy_dict['mc_%s_binsPlot'%(eta)] = dy_dict['mc_%s'%(eta)].Rebin(len(binsPlot), dy_dict['mc_%s'%(eta)].GetName()+'_binsPlot', array.array('d', binsPlot+[maxbin]))
            # binning relevant for the plots, but finer to show the shape
            dy_dict['da_%s_binsFine'%(eta)] = dy_dict['da_%s'%(eta)].Rebin(len(binsFine), dy_dict['da_%s'%(eta)].GetName()+'_binsFine', array.array('d', binsFine+[maxbin]))
            dy_dict['mc_%s_binsFine'%(eta)] = dy_dict['mc_%s'%(eta)].Rebin(len(binsFine), dy_dict['mc_%s'%(eta)].GetName()+'_binsFine', array.array('d', binsFine+[maxbin]))
            # binning for the binned signal regions
            dy_dict['da_%s_binsSR'  %(eta)] = dy_dict['da_%s'%(eta)].Rebin(len(binsSR  ), dy_dict['da_%s'%(eta)].GetName()+'_binsSR'  , array.array('d', binsSR  +[maxbin]))
            dy_dict['mc_%s_binsSR'  %(eta)] = dy_dict['mc_%s'%(eta)].Rebin(len(binsSR  ), dy_dict['mc_%s'%(eta)].GetName()+'_binsSR'  , array.array('d', binsSR  +[maxbin]))

            returndict = copy.deepcopy(dy_dict)
            print 'writing shapes to file'
            for i in ['', '_binsPlot', '_binsFine', '_binsSR']:
                copy.deepcopy(dy_dict['mc_%s%s'%(eta, i)]).Write()
                copy.deepcopy(dy_dict['da_%s%s'%(eta, i)]).Write()

        else:
            print 'taking pre-computed shapes'
            for i in ['', '_binsFine', '_binsPlot', '_binsSR']:
                dy_dict['mc_%s%s'%(eta, i)] = copy.deepcopy(dy_shapes_file.Get('dy_shape_da_%s%s'%(eta, i)))
                dy_dict['da_%s%s'%(eta, i)] = copy.deepcopy(dy_shapes_file.Get('dy_shape_mc_%s%s'%(eta, i)))
            returndict = copy.deepcopy(dy_dict)

    dy_shapes_file.Close()
    return returndict

def makePlot(pHisto, pHistoName, oHisto, oHistoName, plotname, eta, nbs, nbstr):
    #nb = re.sub('[^0-9]', '', nbstring)
    pGraph = r.TGraphErrors(pHisto)
    pGraph.SetFillColorAlpha(r.kBlue+1, 0.2)
    pGraph.SetFillStyle(3001)

    maxCont = pHisto.GetBinContent(pHisto.GetMaximumBin())
    pHisto.GetYaxis().SetRangeUser(0., 1.7*maxCont)

    plot = Canvas.Canvas('results/%s/plot_results_mll_%s_%s%s_nb%s'%(lumi_str, plotname, eta, ('' if not opts.onlyTTbar else '_onlyTT'), str(nbs)), 'png,pdf', 0.6, 0.65, 0.85, 0.82)
    plot.addHisto(pHisto, 'hist'     , pHistoName, 'L' , r.kBlue+1, 1,  0)
    plot.addGraph(pGraph, '2,same'   , pHistoName, 'L' , r.kBlue+1, 1, -1)
    plot.addHisto(oHisto, 'PE,SAME'  , oHistoName, 'PL', r.kBlack , 1,  1)
    if not opts.onlyTTbar: 
        dy_histo = scaleZ(dy_shapes['da_%s_binsPlot'%(eta)], eta, nbs)
        dy_graph = r.TGraphErrors(dy_histo)
        print 'adding dy histo to the canvas'
        dy_histo.SetLineStyle(2)
        dy_histo.SetLineColor(r.kGreen+2)
        plot.addHisto(dy_histo, 'HIST,SAME', 'Drell-Yan ', 'L', r.kGreen+2, 1,  2)
        print 'adding dy graph to the canvas'
        dy_graph.SetFillColorAlpha(r.kGreen+2, 0.2)
        dy_graph.SetFillStyle(3001)
        plot.addGraph(dy_graph, '2,same   ', 'Drell-Yan ', 'L', r.kGreen+2, 1, -1)
    plot.addLatex(0.7, 0.53, eta  , 62)
    plot.addLatex(0.7, 0.45, nbstr, 62)
    plot.saveRatio(1, 0, 0 , lumi, oHisto, pHisto, 0.5, 1.5)


if __name__ == '__main__':

    print 'Starting to produce some good ol\' results...'
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples FilenameWithIngredients', version='%prog 1.0')
    parser.add_option('-t', action='store_true', dest='onlyTTbar', default=False, help='just do OF closure test')
    parser.add_option('-c', action='store_true', dest='onlyClosure', default=False, help='just do the closure test. don\'t bother with data')
    parser.add_option('-b', '--nbs', action='store', type=int, dest='nbs', default=0, help='do this for different numbers of b\'s')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-l', '--loadShapes', action='store_true', dest='loadShapes', default=False, help='reload dy shapes. default is off since this takes a while')
    parser.add_option('-M', '--maxRun', action='store', type=int, dest='maxRun', default=999999, help='max run to use for analysis (run is included)')
    parser.add_option('-m', '--minRun', action='store', type=int, dest='minRun', default=-1    , help='min run to use for analysis (run not included)')

    ## make the options globa.. also the lumi
    global opts, lumi, lumi_str, dy_shapes, nbstring
    global ingMC, ingDA, onZ, treeDA, treeMC
    (opts, args) = parser.parse_args()

    print 'running with these options \n'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)
    print ' \n\n'

    ingMC = helper.ingredients(opts.ingredientFile, 'MC'  )
    ingDA = helper.ingredients(opts.ingredientFile, 'DATA')

    rsfofTable = Tables.makeRSFOFTable(ingDA, ingMC)

    ##print asdf
    print 'Going to load DATA and MC trees...'
    mcDatasets = ['TTLep_pow'] + ([] if opts.onlyTTbar else [ 'DYJetsToLL_M10to50', 'DYJetsToLL_M50'])
    lumi = 2.1
    if lumi == 1.3:
        daDatasets = ['DoubleMuon_Run2015C_25ns_05Oct_v1_runs_246908_258714' , 'DoubleEG_Run2015C_25ns_05Oct_v1_runs_246908_258714' , 'MuonEG_Run2015C_25ns_05Oct_v1_runs_246908_258714' ,
                      'DoubleMuon_Run2015D_05Oct_v1_runs_246908_258751'      , 'DoubleEG_Run2015D_05Oct_v1_runs_246908_258751'      , 'MuonEG_Run2015D_05Oct_v2_runs_246908_258751'      ,
                      'DoubleMuon_Run2015D_v4_runs_246908_258751'            , 'DoubleEG_Run2015D_v4_runs_246908_258751'            , 'MuonEG_Run2015D_v4_runs_246908_258751'            ]
    elif lumi == 2.1:
        daDatasets = ['DoubleMuon_Run2015C_25ns-05Oct_v1_runs_246908_260627' , 'DoubleEG_Run2015C_25ns-05Oct_v1_runs_246908_260627' , 'MuonEG_Run2015C_25ns-05Oct_v1_runs_246908_260627' ,
                      'DoubleMuon_Run2015D-05Oct_v1_runs_246908_260627'      , 'DoubleEG_Run2015D-05Oct_v1_runs_246908_260627'      , 'MuonEG_Run2015D-05Oct_v2_runs_246908_260627'      ,
                      'DoubleMuon_Run2015D_v4_runs_246908_260627'            , 'DoubleEG_Run2015D_v4_runs_246908_260627'            , 'MuonEG_Run2015D_v4_runs_246908_260627'            ]

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print 'Trees successfully loaded...'


    #lumi_str = 'lumi'+str(lumi).replace('.', 'p')+('_inclB' if inclusiveB else '_exclB')
    lumi_str = 'lumi'+str(lumi).replace('.', 'p')+'_forApproval'
    print 'Running with an integrated luminosity of %.2f fb-1' %(lumi)

    isBlinded =False

    ## load the on-Z results from the MET templates
    onZ = onZResult('vinceResults.txt')

   
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()


    finalBinning = range(20,310,10) # + range(81,111,10) + range(110,310,10))
    nbcut_0 = 't.nBJetMedium35_Edge == 0'
    nbcut_1 = 't.nBJetMedium35_Edge >= 1'
    nbcut_2 = 't.nBJetMedium35_Edge >= 2'


    #dy_shapes = getDYShapes(finalBinning, range(20,302,2), [20., 70., 81., 101., 120., 300.])

    #print asdf  
    regions = []
    signalRegionincb = Region.region('signalRegionincb', 
                                 [cuts.METJetsSignalRegion],
                                 ['mll', 'met'],
                                 #[range(20,310,10), range(0,210,10)],
                                 [finalBinning, range(0,210,10)],
                                 True if not opts.onlyClosure else False)
    binnedSRincb     = Region.region('SR8TeVincb', 
                                 [cuts.METJetsSignalRegion],
                                 ['mll'],
                                 [ [20., 70., 81., 101., 120., 13000.] ],
                                 True if not opts.onlyClosure else False)
    signalRegion0b = Region.region('signalRegion0b', 
                                 [cuts.METJetsSignalRegion, nbcut_0],
                                 ['mll', 'met'],
                                 #[range(20,310,10), range(0,210,10)],
                                 [finalBinning, range(0,210,10)],
                                 True if not opts.onlyClosure else False)
    binnedSR0b     = Region.region('SR8TeV0b', 
                                 [cuts.METJetsSignalRegion, nbcut_0],
                                 ['mll'],
                                 [ [20., 70., 81., 101., 120., 13000.] ],
                                 True if not opts.onlyClosure else False)
    signalRegion1b = Region.region('signalRegion1b', 
                                 [cuts.METJetsSignalRegion, nbcut_1],
                                 ['mll', 'met'],
                                 #[range(20,310,10), range(0,210,10)],
                                 [finalBinning, range(0,210,10)],
                                 True if not opts.onlyClosure else False)
    binnedSR1b     = Region.region('SR8TeV1b', 
                                 [cuts.METJetsSignalRegion, nbcut_1],
                                 ['mll'],
                                 [ [20., 70., 81., 101., 120., 13000.] ],
                                 True if not opts.onlyClosure else False)
    signalRegion2b = Region.region('signalRegion2b', 
                                 [cuts.METJetsSignalRegion, nbcut_1],
                                 ['mll', 'met'],
                                 #[range(20,310,10), range(0,210,10)],
                                 [finalBinning, range(0,210,10)],
                                 True if not opts.onlyClosure else False)
    binnedSR2b     = Region.region('SR8TeV2b', 
                                 [cuts.METJetsSignalRegion, nbcut_2],
                                 ['mll'],
                                 [ [20., 70., 81., 101., 120., 13000.] ],
                                 True if not opts.onlyClosure else False)

    #regions.append(signalRegionincb)
    #regions.append(signalRegion0b)
    #regions.append(signalRegion1b)
    #regions.append(signalRegion2b)
    regions.append(binnedSRincb)
    regions.append(binnedSR0b)
    regions.append(binnedSR1b)
    regions.append(binnedSR2b)

    ## ==================================================
    ## look for DY shapes in file, if not there make them
    ## ==================================================
    dy_shapes = getDYShapes(finalBinning, range(20,302,2), [20., 70., 81., 101., 120., 300.])
    ## ==================================================
    ## done with DY shapes. moving on with life          
    ## ==================================================

    #print asdf

    global binnedSR
    binnedSR = copy.deepcopy(binnedSR0b)

    for reg in regions:
        print 'i am at region %s' %(reg.name)
        for eta in ['central', 'forward']:
            print '... in %s' %(eta)

            cuts_sf = cuts.AddList([cuts.MaxRun(opts.maxRun), cuts.MinRun(opts.minRun), cuts.GoodLeptonSF()]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts+[cuts.trigger])
            cuts_of = cuts.AddList([cuts.MaxRun(opts.maxRun), cuts.MinRun(opts.minRun), cuts.GoodLeptonOF()]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts+[cuts.trigger])

            for tree in ([treeMC, treeDA] if reg.doData else [treeMC]):
                dataMC = 'DATA' if tree == treeDA else 'MC'

                if 'mll' in reg.rvars:
                    reg.mll_sf = tree.getTH1F(lumi, 'mll_sf_'+eta+reg.name+dataMC, 't.lepsMll_Edge', reg.bins[reg.rvars.index('mll')], 1, 1, cuts_sf, '', "m_{ll} (GeV)")
                    reg.mll_of = tree.getTH1F(lumi, 'mll_of_'+eta+reg.name+dataMC, 't.lepsMll_Edge', reg.bins[reg.rvars.index('mll')], 1, 1, cuts_of, '', "m_{ll} (GeV)")

                    isData = (dataMC == 'DATA')
                    reg.mll     .setHisto(reg.mll_sf, dataMC, eta)
                    reg.mll_pred.setHisto(makePrediction(reg.mll_of, ingDA if isData else ingMC, eta), dataMC, eta)

                    del reg.mll_sf, reg.mll_of




    ## for signalRegion in [signalRegionincb, signalRegion0b, signalRegion1b]:
    ##     nb = 2 if '2b' in signalRegion.name else 1 if '1b' in signalRegion.name else 'inc' if 'incb' in signalRegion.name else 0
    ##     #nbstr = 'n_{b} '+(' #geq ' if not nb in [0,1] else ' = ')+(str(nb) if nb != 'inc' else str(0))
    ##     if   nb == 'inc': nbstr = 'n_{b} #geq 0'
    ##     elif nb == 0    : nbstr = 'n_{b} = 0'
    ##     elif nb == 1    : nbstr = 'n_{b} #geq 1'
    ##     elif nb == 2    : nbstr = 'n_{b} #geq 2'
    ##     for eta in ['central']:

    ##         ## get the mc prediction and observation
    ##         mcObsHisto  = signalRegion.mll     .getHisto('MC', eta)
    ##         mcPredHisto = signalRegion.mll_pred.getHisto('MC', eta).Clone('mcPredHisto_'+eta)
    ##         #add the DY contribution to the MC prediction
    ##         if not opts.onlyTTbar: 
    ##             print 'adding dy shape to the mcPredHisto'
    ##             mcPredHisto.Add(scaleZ(dy_shapes['mc_%s'%(eta)+'_binsPlot'], eta, nb), 1.)
    ##         ## get the data observed and predicted
    ##         daObsHisto  = signalRegion.mll.getHisto('DATA', eta)
    ##         daPredHisto = copy.deepcopy( signalRegion.mll_pred.getHisto('DATA', eta) )
    ##         if not opts.onlyTTbar:
    ##             print 'adding dy shape to the daPredHisto'
    ##             daPredHisto.Add(scaleZ(dy_shapes['da_%s'%(eta)+'_binsPlot'], eta, nb), 1.)

    ##         ## finally, make the plots.

    ##         makePlot(mcPredHisto, 'pred. (MC)'  , mcObsHisto , 'obs. (MC)'   , 'MCClosure'    , eta, nb, nbstr) ## closure plot
    ##         makePlot(mcPredHisto, 'pred. (MC)'  , daPredHisto, 'pred. (DATA)', 'dataMCPred'   , eta, nb, nbstr) ## prediction comparison
    ##         makePlot(daPredHisto, 'pred. (DATA)', mcObsHisto , 'obs. (MC)'   , 'predDataObsMC', eta, nb, nbstr) ## predData, obsMC
    ##         makePlot(mcObsHisto , 'obs. (MC)'   , daObsHisto , 'obs. (DATA)' , 'dataMCobs'    , eta, nb, nbstr) ## observed data/mc
    ##         makePlot(daPredHisto, 'pred. (DATA)', daObsHisto , 'obs. (DATA)' , 'data'         , eta, nb, nbstr) ## money plot


    # a = Tables.makeConciseTable(binnedSRincb, binnedSR0b, binnedSR1b, ingDA, ingMC, onZ)
    # b = Tables.makeConciseTable(binnedSRincb, binnedSR0b, binnedSR1b, binnedSR2b, ingDA, ingMC, onZ)
    # for i in a: print i
    # for i in b: print i

    ## a_cen    = makeResultsTable(binnedSR, dy_shapes, 'MC', 'central', opts.nbs)
    ## a_fwd    = makeResultsTable(binnedSR, dy_shapes, 'MC', 'forward', opts.nbs)
    ## a_cen_da = makeResultsTable(binnedSR, dy_shapes, 'DATA', 'central', opts.nbs)
    ## a_fwd_da = makeResultsTable(binnedSR, dy_shapes, 'DATA', 'forward', opts.nbs)
    ## print 'CENTRAL TABLE'
    ## for i in a_cen_da: print i
    ## print 'FORWARD TABLE'
    ## for i in a_fwd_da: print i


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
import math, sys, optparse, copy


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample

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

def scaleZ(histo, onZ_pred, eta, nbs):
    bin81  = histo.FindBin(81)
    bin101 = histo.FindBin(100)
    nb_pred = getattr(onZ_pred, '%s_%db'%( ('cen' if eta == 'central' else 'fwd'), nbs))
    ret_hist = copy.deepcopy(histo)
    integral = ret_hist.Integral(bin81, bin101)
    if not integral:
        print 'WARNING: the DY shape histogram has no content around the Z-mass !!!'
        print 'WARNING: this happens for nb %s' %(str(nbs))
        print 'WARNING: putting the whole prediction into one bin at the Z-mass'
        ret_hist.SetBinContent(ret_hist.FindBin(91.), nb_pred)
    else:
        scale  = nb_pred / histo.Integral(bin81, bin101)
        ret_hist.Scale(scale)
    return ret_hist

def makePrediction(of_histo, ing, eta):
    central = (eta == 'central')
    sf_pred = copy.deepcopy(of_histo)
    scale = ing.rsfof_ttcr_lm.cen_val   if central else ing.rsfof_ttcr_lm.fwd_val
    s_err = ing.rsfof_ttcr_lm.cen_stat  if central else ing.rsfof_ttcr_lm.fwd_stat
    print '========================================='
    print 'using rsfof: %.4f' %(scale)
    print '========================================='
    for _bin in range(1,of_histo.GetNbinsX()+1):
        tmp_cont  = of_histo.GetBinContent(_bin)
        tmp_err   = of_histo.GetBinError  (_bin)
        tmp_mass  = of_histo.GetXaxis().GetBinCenter(_bin)
        ##if   20 <= tmp_mass <=  70.:
        ##    scale = ing.rsfof_ttcr_lm.cen_val   if central else ing.rsfof_ttcr_lm.fwd_val
        ##    s_err = ing.rsfof_ttcr_lm.cen_stat  if central else ing.rsfof_ttcr_lm.fwd_stat
        ##elif 70 <  tmp_mass <= 120.:
        ##    scale = ing.rsfof_ttcr_onZ.cen_val  if central else ing.rsfof_ttcr_onZ.fwd_val
        ##    s_err = ing.rsfof_ttcr_onZ.cen_stat if central else ing.rsfof_ttcr_onZ.fwd_stat
        ##elif 120 < tmp_mass:
        ##    scale = ing.rsfof_ttcr_hm.cen_val   if central else ing.rsfof_ttcr_hm.fwd_val
        ##    s_err = ing.rsfof_ttcr_hm.cen_stat  if central else ing.rsfof_ttcr_hm.fwd_stat
        tmp_pred  = scale*tmp_cont
        tmp_prede = math.sqrt(tmp_err*tmp_err + s_err*s_err*tmp_cont*tmp_cont)
        sf_pred.SetBinContent(_bin, tmp_pred )
        sf_pred.SetBinError  (_bin, tmp_prede)
    return sf_pred



if __name__ == '__main__':

    print 'Starting to produce some good ol\' results...'
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples FilenameWithIngredients', version='%prog 1.0')
    parser.add_option('-t', action='store_true', dest='onlyTTbar', default=False, help='just do OF closure test')
    parser.add_option('-c', action='store_true', dest='onlyClosure', default=False, help='just do the closure test. don\'t bother with data')
    parser.add_option('-b', '--nbs', action='store', type=int, dest='nbs', default=0, help='do this for different numbers of b\'s')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-l', '--loadShapes', action='store_true', dest='loadShapes', default=False, help='reload dy shapes. default is off since this takes a while')
    parser.add_option('-m', '--maxRun', action='store', type=int, dest='maxRun', default=999999, help='max run to use for analysis (run is included)')

    (opts, args) = parser.parse_args()

    print 'running with these options \n'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)
    print ' \n\n'

    print 'Going to load DATA and MC trees...'
    mcDatasets = ['TTLep_pow'] + ([] if opts.onlyTTbar else [ 'DYJetsToLL_M10to50', 'DYJetsToLL_M50'])
    daDatasets = ['DoubleMuon_Run2015C_25ns_05Oct_v1_runs_246908_258714' , 'DoubleEG_Run2015C_25ns_05Oct_v1_runs_246908_258714' , 'MuonEG_Run2015C_25ns_05Oct_v1_runs_246908_258714' ,
                  'DoubleMuon_Run2015D_05Oct_v1_runs_246908_258751'      , 'DoubleEG_Run2015D_05Oct_v1_runs_246908_258751'      , 'MuonEG_Run2015D_05Oct_v2_runs_246908_258751'      ,
                  'DoubleMuon_Run2015D_v4_runs_246908_258751'            , 'DoubleEG_Run2015D_v4_runs_246908_258751'            , 'MuonEG_Run2015D_v4_runs_246908_258751'            ]

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print 'Trees successfully loaded...'


    lumi = 1.3 
    #lumi_str = 'lumi'+str(lumi).replace('.', 'p')+('_inclB' if inclusiveB else '_exclB')
    #lumi_str = 'lumi'+str(lumi).replace('.', 'p')+'_e0bge1bge2b'
    lumi_str = 'lumi'+str(lumi).replace('.', 'p')+'_ge0Inclusive_normal'
    print 'Running with an integrated luminosity of %.2f fb-1' %(lumi)

    isBlinded =False

    ## load the on-Z results from the MET templates
    onZ = onZResult('vinceResults.txt')

   
    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()


    ingMC = helper.ingredients(opts.ingredientFile, 'MC'  )
    ingDA = helper.ingredients(opts.ingredientFile, 'DATA')

    finalBinning = (range(20,80,10) + range(81,111,10) + range(110,310,10))
    regions = []
    nbcut = 't.nBJetMedium35_Edge %s 0'%( '>=' )
    if   opts.nbs == 1:
        nbcut = 't.nBJetMedium35_Edge %s 1'%( '>=' )
    elif opts.nbs == 2:
        nbcut = 't.nBJetMedium35_Edge %s 2'%( '>=' )

    signalRegion = Region.region('signalRegion', 
                                 [cuts.METJetsSignalRegion, nbcut],
                                 ['mll', 'met'],
                                 #[range(20,310,10), range(0,210,10)],
                                 [finalBinning, range(0,210,10)],
                                 True if not opts.onlyClosure else False)
    binnedSR     = Region.region('SR8TeV', 
                                 [cuts.METJetsSignalRegion, nbcut],
                                 ['mll'],
                                 [ [20., 70., 81., 101., 120., 500.] ],
                                 True if not opts.onlyClosure else False)

    regions.append(signalRegion)
    regions.append(binnedSR)


    for reg in regions:
        print 'i am at region %s' %(reg.name)
        for eta in ['central', 'forward']:
            print '... in %s' %(eta)

            cuts_sf = cuts.AddList([cuts.MaxRun(opts.maxRun), cuts.GoodLeptonSF()]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts+[cuts.trigger])
            cuts_of = cuts.AddList([cuts.MaxRun(opts.maxRun), cuts.GoodLeptonOF()]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts+[cuts.trigger])

            for tree in ([treeMC, treeDA] if reg.doData else [treeMC]):
                dataMC = 'DATA' if tree == treeDA else 'MC'

                if 'mll' in reg.rvars:
                    reg.mll_sf = tree.getTH1F(lumi, 'mll_sf_'+eta+reg.name+dataMC, 't.lepsMll_Edge', reg.bins[reg.rvars.index('mll')], 1, 1, cuts_sf, '', "m_{ll} (GeV)")
                    reg.mll_of = tree.getTH1F(lumi, 'mll_of_'+eta+reg.name+dataMC, 't.lepsMll_Edge', reg.bins[reg.rvars.index('mll')], 1, 1, cuts_of, '', "m_{ll} (GeV)")

                    isData = (dataMC == 'DATA')
                    reg.mll     .setHisto(reg.mll_sf, dataMC, eta)
                    reg.mll_pred.setHisto(makePrediction(reg.mll_of, ingDA if isData else ingMC, eta), dataMC, eta)

                    del reg.mll_sf, reg.mll_of




    ## ==================================================
    ## look for DY shapes in file, if not there make them
    ## ==================================================
    if not opts.onlyTTbar:
        dy_nomass          = Region.region(''         , [cuts.DYControlRegion], ['mll']       , [range(20,302,2)], True)
        dy_nomass_rightBin = Region.region('_rightBin', [cuts.DYControlRegion], ['mll']       , [finalBinning]   , True)
        dy_nomass_binned   = Region.region('_binned'  , [cuts.DYControlRegion], binnedSR.rvars, binnedSR.bins    , True)
        dy_nomasses = [dy_nomass, dy_nomass_rightBin, dy_nomass_binned]
        dy_shapes = {}
        ## open the file with the dy shapes in READ or RECREATE depending on if you want to reload it
        dy_shapes_file = r.TFile('dy_shapes.root', 'READ')
        dy_shapes_file.cd()

        for dy_reg in dy_nomasses:
            for eta in ['central', 'forward']:
                fileHasShapes = (dy_shapes_file.Get('dy_shape_da_%s%s_%d'%(eta, dy_reg.name, opts.nbs)) and dy_shapes_file.Get('dy_shape_mc_%s%s_%d'%(eta, dy_reg.name, opts.nbs)))
                print 'the file has the shapes:', fileHasShapes
                if opts.loadShapes or not fileHasShapes:
                    ## ========================================
                    ## get the dy shapes from the rinout region
                    ## ========================================
                    print 're-openining the file with UPDATE'
                    dy_shapes_file.ReOpen('UPDATE')
                    etacut = cuts.Central() if eta == 'central' else cuts.Forward()
                    thecutsSF_nomass = cuts.AddList([cuts.GoodLeptonSF(), cuts.trigger] + [etacut] + dy_reg.cuts)
                    thecutsOF_nomass = cuts.AddList([cuts.GoodLeptonOF(), cuts.trigger] + [etacut] + dy_reg.cuts)
                    setattr(dy_reg, 'mll_SF_da_'+eta, treeDA.getTH1F(lumi, "mll_SF_in_" +eta+dy_reg.name+'_data' , "t.lepsMll_Edge", dy_reg.bins[0], 1, 1, thecutsSF_nomass, "", "m_{ll} (GeV)"))
                    setattr(dy_reg, 'mll_OF_da_'+eta, treeDA.getTH1F(lumi, "mll_OF_in_" +eta+dy_reg.name+'_data' , "t.lepsMll_Edge", dy_reg.bins[0], 1, 1, thecutsOF_nomass, "", "m_{ll} (GeV)"))
                    setattr(dy_reg, 'mll_SF_mc_'+eta, treeMC.getTH1F(lumi, "mll_SF_in_" +eta+dy_reg.name+'_mc'   , "t.lepsMll_Edge", dy_reg.bins[0], 1, 1, thecutsSF_nomass, "", "m_{ll} (GeV)"))
                    setattr(dy_reg, 'mll_OF_mc_'+eta, treeMC.getTH1F(lumi, "mll_OF_in_" +eta+dy_reg.name+'_mc'   , "t.lepsMll_Edge", dy_reg.bins[0], 1, 1, thecutsOF_nomass, "", "m_{ll} (GeV)"))

                    setattr(dy_reg, 'dy_shape_da_'+eta, getattr(dy_reg, 'mll_SF_da_'+eta).Clone('dy_shape_da_%s%s_%d'%(eta, dy_reg.name, opts.nbs)))
                    setattr(dy_reg, 'dy_shape_mc_'+eta, getattr(dy_reg, 'mll_SF_mc_'+eta).Clone('dy_shape_mc_%s%s_%d'%(eta, dy_reg.name, opts.nbs)))
                    getattr(dy_reg, 'dy_shape_da_'+eta).Add(getattr(dy_reg, 'mll_OF_da_'+eta), -1.)
                    getattr(dy_reg, 'dy_shape_mc_'+eta).Add(getattr(dy_reg, 'mll_OF_mc_'+eta), -1.)
            
                    dy_shapes['%db_mc_%s%s'%(opts.nbs, eta, dy_reg.name)] = scaleZ(getattr(dy_reg, 'dy_shape_mc_%s'%eta), onZ, eta, opts.nbs)
                    dy_shapes['%db_da_%s%s'%(opts.nbs, eta, dy_reg.name)] = scaleZ(getattr(dy_reg, 'dy_shape_da_%s'%eta), onZ, eta, opts.nbs)

                    print 'writing shape to file'
                    dy_shapes['%db_mc_%s%s'%(opts.nbs, eta, dy_reg.name)].Write()
                    dy_shapes['%db_da_%s%s'%(opts.nbs, eta, dy_reg.name)].Write()

                else:
                    print 'taking pre-computed shapes'
                    dy_shapes['%db_mc_%s%s'%(opts.nbs, eta, dy_reg.name)] = copy.deepcopy(dy_shapes_file.Get('dy_shape_mc_%s%s_%d'%(eta, dy_reg.name, opts.nbs)))
                    dy_shapes['%db_da_%s%s'%(opts.nbs, eta, dy_reg.name)] = copy.deepcopy(dy_shapes_file.Get('dy_shape_da_%s%s_%d'%(eta, dy_reg.name, opts.nbs)))

        dy_shapes_file.Close()

        for key,value in dy_shapes.items():
            value.SetLineStyle(2)
            value.SetLineColor(r.kRed+2)
            if 'rightBin' in key:
                value.SetLineColor(r.kGreen+2)
    ## ==================================================
    ## done with DY shapes. moving on with life          
    ## ==================================================

    tables = []
    for eta in ['central', 'forward']:
        yscale = 75. if eta == 'central' else 40.
        ## ===============================
        ## MAKE THE Mll PLOT for MC only
        ## ===============================

        mcObsHisto  = signalRegion.mll.getHisto('MC', eta)
        mcPredHisto = signalRegion.mll_pred.getHisto('MC'  , eta).Clone('mcPredHisto_'+eta)
        if not opts.onlyTTbar: 
            mcPredHisto.Add(dy_shapes['%db_mc_%s'%(opts.nbs,eta)+'_rightBin'], 1.)

        mcPredGraph = r.TGraphErrors(mcPredHisto)
        mcPredGraph.SetFillColorAlpha(r.kBlue+1, 0.2)
        mcPredGraph.SetFillStyle(3001)

        #mcPredHisto.GetYaxis().SetRangeUser(0.,  yscale*lumi)
        mcPredHisto.GetYaxis().SetRangeUser(0.,  1.6*mcPredHisto.GetMaximum() )

        plot_closure = Canvas.Canvas('results/%s/plot_results_mll_MCClosure_%s%s_nb%d'%(lumi_str, eta, ('' if not opts.onlyTTbar else '_onlyTT'), opts.nbs), 'png,pdf', 0.6, 0.5, 0.85, 0.7)
        plot_closure.addHisto(mcPredHisto, 'hist'     , 'pred. (MC)', 'L', r.kBlue+1 , 1, 0)
        plot_closure.addGraph(mcPredGraph, '2,same'   , 'pred. (MC)', 'L', r.kBlue+1  , 1, -1)
        plot_closure.addHisto(mcObsHisto , 'PE,SAME'  , 'obs.  (MC)', 'PL', r.kBlack, 1, 1)
        if not opts.onlyTTbar: plot_closure.addHisto(dy_shapes['%db_mc_%s'%(opts.nbs,eta)], 'HIST,SAME', 'DY shape  (MC)', 'L', r.kRed+2, 1, 1)
        plot_closure.addLatex(0.7, 0.8, eta, 62)
        plot_closure.addLatex(0.7, 0.45, 'n_{b} '+('#geq ' if opts.nbs != 0 else '= ')+str(opts.nbs), 62)
        plot_closure.saveRatio(1, 0, 0, lumi, mcObsHisto, mcPredHisto, 0.5, 1.5)

        if not opts.onlyClosure:
            daPredHisto = signalRegion.mll_pred.getHisto('DATA', eta)
            daPredHisto.GetYaxis().SetRangeUser(0.,  1.8*daPredHisto.GetMaximum() )
            if not opts.onlyTTbar: daPredHisto.Add(dy_shapes['%db_da_%s'%(opts.nbs,eta)+'_rightBin'], 1.)
            ## ============================================================
            ## MAKE THE Mll PLOT for DATA MC PREDICTION COMPARISON
            ## ============================================================
            plot_dataMC = Canvas.Canvas('results/%s/plot_results_mll_dataMCPred_%s%s_nb%d'%(lumi_str, eta, ('' if not opts.onlyTTbar else '_onlyTT'), opts.nbs), 'png,pdf', 0.6, 0.5, 0.85, 0.7)
            plot_dataMC.addHisto(mcPredHisto, 'hist'   , 'pred. (MC)'  , 'L', r.kBlue+1 , 1, 0)
            plot_dataMC.addGraph(mcPredGraph, '2,same' , 'pred. (MC)'  , 'L', r.kBlue+1  , 1, -1)
            plot_dataMC.addHisto(daPredHisto, 'PE,SAME', 'pred. (DATA)', 'PL', r.kBlack, 1, 1)
            if not opts.onlyTTbar: plot_dataMC.addHisto(dy_shapes['%db_da_%s'%(opts.nbs,eta)], 'HIST,SAME', 'DY shape  (DATA)', 'L', r.kRed+2, 1, 1)
            plot_dataMC.addLatex(0.7, 0.8, eta, 62)
            plot_dataMC.addLatex(0.7, 0.45, 'n_{b} '+('#geq ' if opts.nbs != 0 else '= ')+str(opts.nbs), 62)
            plot_dataMC.saveRatio(1, 0, 0, lumi, daPredHisto, mcPredHisto, 0.5, 1.5)

            tables.append(makeRatioTable(binnedSR, 'MC', eta, opts.nbs))

            ## ============================================================
            ## MAKE THE Mll PLOT for DATA PREDICTION vs MC OBSERVATION
            ## ============================================================
            plot_predDataMC = Canvas.Canvas('results/%s/plot_results_mll_predDataObsMC_%s%s_nb%d'%(lumi_str, eta, ('' if not opts.onlyTTbar else '_onlyTT'), opts.nbs), 'png,pdf', 0.6, 0.5, 0.85, 0.7)
            plot_predDataMC.addHisto(daPredHisto, 'hist'   , 'SR - pred. (DATA)', 'PL', r.kBlack, 1, 1)
            plot_predDataMC.addGraph(mcPredHisto, 'PE,same', 'SR - obs. (MC)'   , 'L' , r.kBlue+1  , 1, -1)
            if not opts.onlyTTbar: plot_predDataMC.addHisto(dy_shapes['%db_da_%s'%(opts.nbs,eta)], 'HIST,SAME', 'DY shape  (DATA)', 'L', r.kRed+2, 1, 1)
            plot_predDataMC.addLatex(0.7, 0.8, eta, 62)
            plot_predDataMC.addLatex(0.7, 0.45, 'n_{b} '+('#geq ' if opts.nbs != 0 else '= ')+str(opts.nbs), 62)
            plot_predDataMC.saveRatio(1, 0, 0, lumi, mcPredHisto, daPredHisto, 0.5, 1.5)


            ## these two only when not blinded!!!
            if not isBlinded:
                ## =================
                ## MAKE THE Mll PLOT for DATA only
                ## =================
                daObsHisto = signalRegion.mll.getHisto('DATA', eta)
                daPredHisto.GetXaxis().SetTitle('m_{ll}')
                daPredGraph = r.TGraphErrors(daPredHisto)
                daPredGraph.SetFillColorAlpha(r.kBlue+1, 0.2)
                daPredGraph.SetFillStyle(3001)
                plot_result = Canvas.Canvas('results/%s/plot_results_mll_data_%s%s_nb%d'%(lumi_str, eta, ('' if not opts.onlyTTbar else '_onlyTT'), opts.nbs), 'png,pdf', 0.6, 0.5, 0.80, 0.7)
                plot_result.addHisto(daPredHisto, 'hist'   , 'pred. (DATA)', 'L', r.kBlue+1 , 1, 0)
                plot_result.addGraph(daPredGraph, '2,same' , 'pred. (DATA)', 'L', r.kBlue+1  , 1, -1)
                plot_result.addHisto(daObsHisto , 'PE,SAME', 'obs.  (DATA)', 'PL', r.kBlack, 1, 1)
                if not opts.onlyTTbar: plot_result.addHisto(dy_shapes['%db_da_%s'%(opts.nbs,eta)], 'HIST,SAME', 'DY shape  (DATA)', 'L', r.kRed+2, 1, 1)
                plot_result.addLatex(0.7, 0.8, eta, 62)
                plot_result.addLatex(0.7, 0.45, 'n_{b} '+('#geq ' if opts.nbs != 0 else '= ')+str(opts.nbs), 62)
                plot_result.saveRatio(1, 0, 0, lumi, daObsHisto, daPredHisto, 0.5, 1.5)
            else:
                print 'you sneaky bastard. stop this!!!'
    a_cen = makeResultsTable(binnedSR, dy_shapes, 'MC', 'central', opts.nbs)
    a_fwd = makeResultsTable(binnedSR, dy_shapes, 'MC', 'forward', opts.nbs)
    a_cen_da = makeResultsTable(binnedSR, dy_shapes, 'DATA', 'central', opts.nbs)
    a_fwd_da = makeResultsTable(binnedSR, dy_shapes, 'DATA', 'forward', opts.nbs)
    print 'CENTRAL TABLE'
    for i in a_cen_da: print i
    print 'FORWARD TABLE'
    for i in a_fwd_da: print i


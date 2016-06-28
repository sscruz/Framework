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
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats, TStyle
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
        self.DYcen_incb   = self.DYcen_incbtag
        self.DYfwd_incb   = self.DYfwd_incbtag
        self.DYcen_incb_e = self.DYcen_incbtag_err
        self.DYfwd_incb_e = self.DYfwd_incbtag_err

        self.DYcen_0b   = self.DYcen_0btag
        self.DYfwd_0b   = self.DYfwd_0btag
        self.DYcen_0b_e = self.DYcen_0btag_err
        self.DYfwd_0b_e = self.DYfwd_0btag_err

        self.DYcen_1b   = self.DYcen_1btag
        self.DYfwd_1b   = self.DYfwd_1btag
        self.DYcen_1b_e = self.DYcen_1btag_err
        self.DYfwd_1b_e = self.DYfwd_1btag_err


        self.OTHERcen_incb   = self.OTHERcen_incbtag
        self.OTHERfwd_incb   = self.OTHERfwd_incbtag
        self.OTHERcen_incb_e = self.OTHERcen_incbtag_err
        self.OTHERfwd_incb_e = self.OTHERfwd_incbtag_err

        self.OTHERcen_0b   = self.OTHERcen_0btag
        self.OTHERfwd_0b   = self.OTHERfwd_0btag
        self.OTHERcen_0b_e = self.OTHERcen_0btag_err
        self.OTHERfwd_0b_e = self.OTHERfwd_0btag_err

        self.OTHERcen_1b   = self.OTHERcen_1btag
        self.OTHERfwd_1b   = self.OTHERfwd_1btag
        self.OTHERcen_1b_e = self.OTHERcen_1btag_err
        self.OTHERfwd_1b_e = self.OTHERfwd_1btag_err

        self.cen_incb   = self.cen_incbtag_2jet + self.cen_incbtag_3jet
        self.fwd_incb   = self.fwd_incbtag_2jet + self.fwd_incbtag_3jet
        self.cen_incb_e = math.sqrt(self.cen_incbtag_2jet_err**2 + self.cen_incbtag_3jet_err**2)
        self.fwd_incb_e = math.sqrt(self.fwd_incbtag_2jet_err**2 + self.fwd_incbtag_3jet_err**2)

        self.cen_0b   = self.cen_0btag_2jet + self.cen_0btag_3jet
        self.fwd_0b   = self.fwd_0btag_2jet + self.fwd_0btag_3jet
        self.cen_0b_e = math.sqrt(self.cen_0btag_2jet_err**2 + self.cen_0btag_3jet_err**2)
        self.fwd_0b_e = math.sqrt(self.fwd_0btag_2jet_err**2 + self.fwd_0btag_3jet_err**2)

        self.cen_1b   = self.cen_1btag_2jet + self.cen_1btag_3jet + self.cen_2btag_2jet + self.cen_2btag_3jet
        self.fwd_1b   = self.fwd_1btag_2jet + self.fwd_1btag_3jet + self.fwd_2btag_2jet + self.fwd_2btag_3jet
        self.cen_1b_e = math.sqrt(self.cen_1btag_2jet_err**2 + self.cen_1btag_3jet_err**2 + self.cen_2btag_2jet_err**2 + self.cen_2btag_3jet_err**2)
        self.fwd_1b_e = math.sqrt(self.fwd_1btag_2jet_err**2 + self.fwd_1btag_3jet_err**2 + self.fwd_2btag_2jet_err**2 + self.fwd_2btag_3jet_err**2)

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

def makePrediction(of_histo, ing):
    #central = (eta == 'central')
    sf_pred = copy.deepcopy(of_histo.Clone(of_histo.GetName().replace('of','pred')))
    scale = ing.rsfof_final_cen_val  ##if central else ing.rsfof_final_fwd_val
    s_err = ing.rsfof_final_cen_err  ##if central else ing.rsfof_final_fwd_err
    print '========================================='
    print 'using rsfof for %s : %.4f +- %.4f' %(ing.dType, scale, s_err)
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

def makePlot(pHisto, pHistoName, oHisto, oHistoName, plotname):#, eta, nbs, nbstr):

    colors = helper.createMyColors()

    #nb = re.sub('[^0-9]', '', nbstring)
    pGraph = r.TGraphErrors(pHisto)
    pGraph.SetFillColorAlpha(r.kBlue+3, 0.2)
    pGraph.SetFillStyle(3001)
    pGraph.SetFillColor(helper.myColors["MyBlue"])
    pGraph.SetMarkerStyle(20)
    pHisto.SetMarkerStyle(20)
    pHisto.SetLineColor(r.kBlack)
    oHisto.SetLineColor(r.kBlue+3)
    oHisto.SetLineWidth(2)

    maxCont = pHisto.GetBinContent(pHisto.GetMaximumBin())
    pHisto.GetYaxis().SetRangeUser(0., 1.7*maxCont)

    plot = Canvas.Canvas('2016/results/%s/plot_results_%s%s'%(lumi_str, plotname, ('' if not opts.onlyTTbar else '_onlyTT')), 'png,pdf', 0.6, 0.65, 0.85, 0.82)
    plot.addHisto(pHisto, 'hist'     , pHistoName, 'L' , r.kBlack, 1,  0)
    plot.addGraph(pGraph, '2,same'   , pHistoName, 'L' , r.kBlack, 1, -1)
    plot.addHisto(oHisto, 'PE,SAME'  , oHistoName, 'PL', r.kBlue+3 , 1,  1)


    if not opts.onlyTTbar: 
        dy_histo = scaleZ(dy_shapes['da_%s_binsPlot'%(eta)], eta, nbs)
        dy_histo.SetLineColor(r.kGreen+3)
        dy_histo.SetFillColor(r.kGreen+3)
        dy_histo.SetFillStyle(3002)
        dy_graph = r.TGraphErrors(dy_histo)
        print 'adding dy histo to the canvas'
        dy_histo.SetLineStyle(2)
        dy_histo.SetLineColor(r.kGreen+2)
        plot.addHisto(dy_histo, 'HIST,SAME', 'Drell-Yan ', 'L', r.kGreen+3, 1,  2)
        print 'adding dy graph to the canvas'
        dy_graph.SetFillColorAlpha(r.kGreen+3, 0.2)
        dy_graph.SetFillStyle(3001)
        plot.addGraph(dy_graph, '2,same   ', 'Drell-Yan ', 'L', r.kGreen+3, 1, -1)
    #plot.addLatex(0.7, 0.53, eta  , 62)
    #plot.addLatex(0.7, 0.45, nbstr, 62)
    plot.saveRatio(1, 0, 0 , lumi, oHisto, pHisto, 0.5, 1.5)

def makePlots(srlist):
    for sR in srlist:
        print 'making plots for region', sR.name
        for v in sR.rvars:
            print '.. and variable',v
            #if not getattr(sR, v).getHisto('MC'): continue
            ## get the mc prediction and observation
            mcObsHisto  = getattr(sR, v        ).getHisto('MC')
            mcPredHisto = getattr(sR, v+'_pred').getHisto('MC').Clone('mcPredHisto')
            #add the DY contribution to the MC prediction
            if not opts.onlyTTbar: 
                print 'adding dy shape to the mcPredHisto'
                mcPredHisto.Add(scaleZ(dy_shapes['mc_%s'%(eta)+'_binsPlot'], eta, nb), 1.)
            ## get the data observed and predicted
            daObsHisto  = getattr(sR, v).getHisto('DATA')
            daPredHisto = copy.deepcopy( getattr(sR, v+'_pred').getHisto('DATA') )
            if not opts.onlyTTbar:
                print 'adding dy shape to the daPredHisto'
                daPredHisto.Add(scaleZ(dy_shapes['da_%s'%(eta)+'_binsPlot'], eta, nb), 1.)

            ## finally, make the plots.
            makePlot(mcPredHisto, 'pred. (MC)'  , mcObsHisto , 'obs. (MC)'   , '{sn}_{vn}_MCClosure'    .format(sn=sR.name, vn=v) )## , eta, nb, nbstr) ## closure plot
            if opts.onlyClosure: continue
            makePlot(mcPredHisto, 'pred. (MC)'  , daPredHisto, 'pred. (DATA)', '{sn}_{vn}_dataMCPred'   .format(sn=sR.name, vn=v) )## , eta, nb, nbstr) ## prediction comparison
            makePlot(daPredHisto, 'pred. (DATA)', mcObsHisto , 'obs. (MC)'   , '{sn}_{vn}_predDataObsMC'.format(sn=sR.name, vn=v) )## , eta, nb, nbstr) ## predData, obsMC
            makePlot(mcObsHisto , 'obs. (MC)'   , daObsHisto , 'obs. (DATA)' , '{sn}_{vn}_dataMCobs'    .format(sn=sR.name, vn=v) )## , eta, nb, nbstr) ## observed data/mc
            makePlot(daPredHisto, 'pred. (DATA)', daObsHisto , 'obs. (DATA)' , '{sn}_{vn}_data'         .format(sn=sR.name, vn=v) )## , eta, nb, nbstr) ## money plot

def scaleByRSFOF(histo, rsfof, rsfof_err):
    h_rsfof = copy.deepcopy(histo)
    h_rsfof.SetName('h_rsfof')
    for i in range(1, h_rsfof.GetNbinsX()+1):
        h_rsfof.SetBinContent(i,rsfof)
        h_rsfof.SetBinError  (i,rsfof_err)
    histo.Multiply(h_rsfof)
    return histo

def makeResultTable(resultPlotLoNLL, resultplotHiNLL):
    h_obsLoNLL = resultPlotLoNLL.histos[0]
    h_preLoNLL = resultPlotLoNLL.histos[1]
    h_obsHiNLL = resultPlotHiNLL.histos[0]
    h_preHiNLL = resultPlotHiNLL.histos[1]
    bin01 = 1
    bin90 = h_obsLoNLL.FindBin(90.)
    binZZ = h_obsLoNLL.GetNbinsX()+1

    dyTotal = 11.8; dyTotal_e = 3.1
    rinoutLoM = 0.109; rinoutLoM_e = math.sqrt(0.001**2 + 0.027**2)
    rinoutHiM = 0.063; rinoutHiM_e = math.sqrt(0.001**2 + 0.016**2)
    dyEffLoNLL = 0.65
    dyEffHiNLL = 1.-dyEffLoNLL

    prFSloNloM_e, prFSloNhiM_e = r.Double(),  r.Double()
    prFShiNloM_e, prFShiNhiM_e = r.Double(),  r.Double()

    prFSloNloM = h_preLoNLL.IntegralAndError(bin01, bin90, prFSloNloM_e)
    prDYloNloM = dyTotal*rinoutLoM*dyEffLoNLL; prDYloNloM_e = math.sqrt(prDYloNloM)
    prTOloNloM = prFSloNloM+prDYloNloM; prTOloNloM_e = math.sqrt(prFSloNloM_e**2 + prDYloNloM_e**2)
    obTOloNloM = h_obsLoNLL.Integral        (bin01, bin90              )
    prFSloNhiM = h_preLoNLL.IntegralAndError(bin90, binZZ, prFSloNhiM_e)
    prDYloNhiM = dyTotal*rinoutHiM*dyEffLoNLL; prDYloNhiM_e = math.sqrt(prDYloNhiM)
    prTOloNhiM = prFSloNhiM+prDYloNhiM; prTOloNhiM_e = math.sqrt(prFSloNhiM_e**2 + prDYloNhiM_e**2)
    obTOloNhiM = h_obsLoNLL.Integral        (bin90, binZZ              )

    prFShiNloM = h_preHiNLL.IntegralAndError(bin01, bin90, prFShiNloM_e)
    prDYhiNloM = dyTotal*rinoutLoM*dyEffHiNLL; prDYhiNloM_e = math.sqrt(prDYhiNloM)
    prTOhiNloM = prFShiNloM+prDYhiNloM; prTOhiNloM_e = math.sqrt(prFShiNloM_e**2 + prDYhiNloM_e**2)
    obTOhiNloM = h_obsHiNLL.Integral        (bin01, bin90              )
    prFShiNhiM = h_preHiNLL.IntegralAndError(bin90, binZZ, prFShiNhiM_e)
    prDYhiNhiM = dyTotal*rinoutHiM*dyEffHiNLL; prDYhiNhiM_e = math.sqrt(prDYhiNhiM)
    prTOhiNhiM = prFShiNhiM+prDYhiNhiM; prTOhiNhiM_e = math.sqrt(prFShiNhiM_e**2 + prDYhiNhiM_e**2)
    obTOhiNhiM = h_obsHiNLL.Integral        (bin90, binZZ              )
    
    resultTable = '''\\documentclass[12pt,a4paper]{{article}}
\\usepackage{{multirow}}
\\begin{{document}}
\\begin{{table}}[hbtp] 
\\begin{{center}} 
\\bgroup 
\\def\\arraystretch{{1.2}} 
\\caption{{Predicted and observed yields for 0.8 fb$^{{-1}}$ of 2016 data.}} 
\\label{{tab:resultTableData}} 
\\begin{{tabular}}{{r l c c }} 
              &                 & ttbar-like  & non-ttbar-like\\\\ \\hline
\\multirow{{4}}{{*}}{{mll $<$ 81 GeV}}       & pred. FS        & {prFSloNloM:.2f}  $\\pm$  {prFSloNloM_e:.2f}    & {prFShiNloM:.2f}     $\\pm$  {prFShiNloM_e:.2f}  \\\\
                                             & pred. DY        & {prDYloNloM:.2f}  $\\pm$  {prDYloNloM_e:.2f}    & {prDYhiNloM:.2f}     $\\pm$  {prDYhiNloM_e:.2f}  \\\\
                                             & pred. total     & {prTOloNloM:.2f}  $\\pm$  {prTOloNloM_e:.2f}    & {prTOhiNloM:.2f}     $\\pm$  {prTOhiNloM_e:.2f}  \\\\
                                             & \\textbf{{obs}} & \\textbf{{{obTOloNloM}}}                        & \\textbf{{{obTOhiNloM}}}                         \\\\ \\hline
\\multirow{{4}}{{*}}{{mll $>$ 101 GeV}}      & pred. FS     & {prFSloNhiM:.2f}  $\\pm$  {prFSloNhiM_e:.2f}       & {prFShiNhiM:.2f}     $\\pm$  {prFShiNhiM_e:.2f}  \\\\
                                             & pred. DY     & {prDYloNhiM:.2f}  $\\pm$  {prDYloNhiM_e:.2f}       & {prDYhiNhiM:.2f}     $\\pm$  {prDYhiNhiM_e:.2f}  \\\\
                                             & pred. total  & {prTOloNhiM:.2f}  $\\pm$  {prTOloNhiM_e:.2f}       & {prTOhiNhiM:.2f}     $\\pm$  {prTOhiNhiM_e:.2f}  \\\\
                                             & \\textbf{{obs}} & \\textbf{{{obTOloNhiM}}}                        & \\textbf{{{obTOhiNhiM}}}                         \\\\
\\end{{tabular}} 
\\egroup 
\\end{{center}} 
\\end{{table}} 
\\end{{document}}'''.format(
    prFSloNloM_e=prFSloNloM_e, prDYloNloM_e=prDYloNloM_e, prTOloNloM_e=prTOloNloM_e, prFSloNhiM_e=prFSloNhiM_e, prDYloNhiM_e=prDYloNhiM_e, prTOloNhiM_e=prTOloNhiM_e,
    prFShiNloM_e=prFShiNloM_e, prDYhiNloM_e=prDYhiNloM_e, prTOhiNloM_e=prTOhiNloM_e, prFShiNhiM_e=prFShiNhiM_e, prDYhiNhiM_e=prDYhiNhiM_e, prTOhiNhiM_e=prTOhiNhiM_e,

    prFSloNloM = prFSloNloM,      prFShiNloM = prFShiNloM,
    prDYloNloM = prDYloNloM,      prDYhiNloM = prDYhiNloM,
    prTOloNloM = prTOloNloM,      prTOhiNloM = prTOhiNloM,
    obTOloNloM = int(obTOloNloM), obTOhiNloM = int(obTOhiNloM),
    prFSloNhiM = prFSloNhiM,      prFShiNhiM = prFShiNhiM,
    prDYloNhiM = prDYloNhiM,      prDYhiNhiM = prDYhiNhiM,
    prTOloNhiM = prTOloNhiM,      prTOhiNhiM = prTOhiNhiM,
    obTOloNhiM = int(obTOloNhiM), obTOhiNhiM = int(obTOhiNhiM))

    compTableFile = open('plots/results/tables/resultTable_800invpb.tex','w')
    compTableFile.write(resultTable)
    compTableFile.close()

def makeClosureTests(var, specialcut = '', scutstring = '', doCumulative = False):
    rsfof_mc = 1.060; rsfof_mc_err = 0.046

    if var == 'mll':
        treevar = 'lepsMll_Edge'
        nbins, xmin, xmax = 23, 20, 250
        xlabel = 'm_{ll} (GeV)'
    elif var == 'nll':
        treevar = 'nll_Edge'
        nbins, xmin, xmax = 52, 10, 36
        xlabel = 'NLL'
        

    ## ## mll ditributions
    mc_OF = treeMC.getTH1F(lumi, var+"mc_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.OF, cuts.Zveto]), '', xlabel)
    mc_SF = treeMC.getTH1F(lumi, var+"mc_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.SF, cuts.Zveto]), '', xlabel)
    da_OF = treeDA.getTH1F(lumi, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.OF, cuts.Zveto]), '', xlabel)
    da_SF = treeDA.getTH1F(lumi, var+"da_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.SF, cuts.Zveto]), '', xlabel)
    dy_SF = treeDY.getTH1F(lumi, var+"dy_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.SF, cuts.Zveto]), '', xlabel)
    mc_OF_err = copy.deepcopy(mc_OF)
    mc_OF_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    mc_OF_err.SetFillStyle(3004); mc_OF_err.SetMarkerSize(0.)
    dy_SF.SetFillColorAlpha(r.kGreen+2,0.5)

    ## da_mll_nllInc_OF_err = copy.deepcopy(da_mll_nllInc_OF)
    ## da_mll_nllInc_OF_err.SetFillColorAlpha(r.kBlack, 0.8)
    ## da_mll_nllInc_OF_err.SetFillStyle(3004); da_mll_nllInc_OF_err.SetMarkerSize(0.)

    mc_SF.GetYaxis().SetRangeUser(0., 1.3*mc_SF.GetMaximum())
    print helper.bcolors.HEADER + '[MC only closure test not scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    plot_closure_noRSFOF = Canvas.Canvas('closure/%s/plot_closure_%s_mcPredmcObs%s_noRSFOF'%(lumi_str, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.6, 0.6, 0.75, 0.8)
    plot_closure_noRSFOF.addHisto(mc_SF    , 'PE'       , 'MC - SF', 'PL', r.kRed+1  , 1,  0)
    plot_closure_noRSFOF.addHisto(mc_OF_err, 'e2,same'  , ''       , 'PL', r.kBlue+1 , 1, -1)
    plot_closure_noRSFOF.addHisto(mc_OF    , 'hist,SAME', 'MC - OF', 'L' , r.kBlue+1 , 1,  1)
    plot_closure_noRSFOF.addHisto(dy_SF    , 'hist,SAME', 'DY - SF', 'FL', r.kGreen+2, 1,  2)
    plot_closure_noRSFOF.addLatex (0.2, 0.8, 'no R_{SFOF} scaling')
    plot_closure_noRSFOF.saveRatio(1, 0, 0, lumi, mc_SF, mc_OF, 0.2, 1.8)

    mc_OF_rsfofScaled = copy.deepcopy(mc_OF)
    mc_OF_rsfofScaled = scaleByRSFOF(mc_OF_rsfofScaled, rsfof_mc, rsfof_mc_err)

    mc_OF_rsfofScaled_err = copy.deepcopy(mc_OF_rsfofScaled)
    mc_OF_rsfofScaled_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    mc_OF_rsfofScaled_err.SetFillStyle(3004); mc_OF_rsfofScaled_err.SetMarkerSize(0.)

    print helper.bcolors.HEADER + '[MC only closure test scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    plot_closure = Canvas.Canvas('closure/%s/plot_closure_%s_mcPredmcObs%s'%(lumi_str, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.6, 0.6, 0.75, 0.8)
    plot_closure.addHisto(mc_SF                , 'PE'       , 'MC - SF', 'PL', r.kRed+1  , 1,  0)
    plot_closure.addHisto(mc_OF_rsfofScaled_err, 'e2,same'  , ''       , 'PL', r.kBlue+1 , 1, -1)
    plot_closure.addHisto(mc_OF_rsfofScaled    , 'hist,SAME', 'MC - OF', 'L' , r.kBlue+1 , 1,  1)
    plot_closure.addHisto(dy_SF                , 'hist,SAME', 'DY - SF', 'FL', r.kGreen+2, 1,  2)
    plot_closure.addLatex (0.2, 0.8, 'R_{SFOF} scaled')
    plot_closure.saveRatio(1, 0, 0, lumi, mc_SF, mc_OF_rsfofScaled, 0.2, 1.8)

    ## if True:
    ##     mc_OF_rsfofScaled    .Scale(0.8/10.)
    ##     mc_OF_rsfofScaled_err.Scale(0.8/10.)
    ##     dy_SF                .Scale(0.8/10.)
    ##     plot_closure = Canvas.Canvas('closure/%s/plot_closure_%s_mcPreddaObs%s'%(lumi_str, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.6, 0.6, 0.75, 0.8)
    ##     plot_closure.addHisto(da_SF                , 'PE'       , 'data-SF', 'PL', r.kRed+1  , 1,  0)
    ##     plot_closure.addHisto(mc_OF_rsfofScaled_err, 'e2,same'  , ''       , 'PL', r.kBlue+1 , 1, -1)
    ##     plot_closure.addHisto(mc_OF_rsfofScaled    , 'hist,SAME', 'MC - OF', 'L' , r.kBlue+1 , 1,  1)
    ##     plot_closure.addHisto(dy_SF                , 'hist,SAME', 'DY - SF', 'FL', r.kGreen+2, 1,  2)
    ##     plot_closure.addLatex (0.2, 0.8, 'R_{SFOF} scaled')
    ##     plot_closure.saveRatio(1, 0, 0, 0.8 , mc_SF, mc_OF_rsfofScaled, 0.2, 1.8)

    ## make cumulative distributions
    if doCumulative:
        mc_SF_cum                 = copy.deepcopy(mc_SF                ).GetCumulative()
        mc_OF_rsfofScaled_err_cum = copy.deepcopy(mc_OF_rsfofScaled_err).GetCumulative()
        mc_OF_rsfofScaled_cum     = copy.deepcopy(mc_OF_rsfofScaled    ).GetCumulative()
        dy_SF_cum                 = copy.deepcopy(dy_SF                ).GetCumulative()
        mc_SF_cum            .Scale(1./mc_SF            .Integral()); mc_SF_cum            .SetLineWidth(2)
        mc_OF_rsfofScaled_cum.Scale(1./mc_OF_rsfofScaled.Integral()); mc_OF_rsfofScaled_cum.SetLineWidth(2)
        dy_SF_cum            .Scale(1./dy_SF            .Integral()); dy_SF_cum            .SetLineWidth(2); dy_SF_cum.SetFillColor(r.kWhite)
        mc_SF_cum.GetYaxis() .SetRangeUser(0.9, 1.01)
        print helper.bcolors.HEADER + '[MC only cumulative distribution, scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
        plot_cumulative = Canvas.Canvas('closure/%s/plot_cumulative_%s_mcPredmcObs%s'%(lumi_str, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.8, 0.2, 0.95, 0.4)
        plot_cumulative.addHisto(mc_SF_cum                 , 'hist,same', 'MC - SF', 'L' , r.kRed+1  , 1,  0)
        plot_cumulative.addHisto(mc_OF_rsfofScaled_cum     , 'hist,SAME', 'MC - OF', 'L', r.kBlue+1 , 1,  1)
        #plot_cumulative.addHisto(dy_SF_cum                 , 'hist,SAME', 'DY - SF', 'L', r.kGreen+2, 1,  2)
        plot_cumulative.addLatex (0.2, 0.8, 'R_{SFOF} scaled')
        plot_cumulative.saveRatio(1, 0, 0, lumi, mc_SF_cum, mc_OF_rsfofScaled_cum, 0.2, 1.8)
        return mc_SF_cum

def makeResultData(var, maxrun = 274240, lint = 0.864, specialcut = '', scutstring = '', returnplot = False):
    rsfof_da = 1.080 ; rsfof_da_err = 0.050
    if var == 'mll':
        treevar = 'lepsMll_Edge'
        nbins, xmin, xmax = 23, 20, 250
        xlabel = 'm_{ll} (GeV)'
    elif var == 'nll':
        treevar = 'nll_Edge'
        nbins, xmin, xmax = 26, 10, 36
        xlabel = 'NLL'
        

    newLumiString = str(lint).replace('.','p')+'invfb'
    if not specialcut:
        specialcut = 'run_Edge < {run}'.format(run=maxrun)
    else:
        specialcut += ' && run_Edge < {run}'.format(run=maxrun)
    ## ## mll ditributions
    da_SF = treeDA.getTH1F(lint, var+"da_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.SF, cuts.Zveto]), '', xlabel)
    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.OF, cuts.Zveto]), '', xlabel)
    da_OF_err = copy.deepcopy(da_OF)
    da_OF_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    da_OF_err.SetFillStyle(3004); da_OF_err.SetMarkerSize(0.)

    da_OF_rsfofScaled = copy.deepcopy(da_OF)
    da_OF_rsfofScaled = scaleByRSFOF(da_OF_rsfofScaled, rsfof_da, rsfof_da_err)

    da_OF_rsfofScaled_err = copy.deepcopy(da_OF_rsfofScaled)
    da_OF_rsfofScaled_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    da_OF_rsfofScaled_err.SetFillStyle(3004); da_OF_rsfofScaled_err.SetMarkerSize(0.)

    print helper.bcolors.HEADER + '[result scaled by RSFOF for DATA] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    plot_result = Canvas.Canvas('results/%s/plot_result_%s_daPreddaObs%s'%(newLumiString, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.6, 0.6, 0.75, 0.8)
    plot_result.addHisto(da_SF                , 'PE'       , 'data - SF', 'PL', r.kRed+1  , 1,  0)
    plot_result.addHisto(da_OF_rsfofScaled_err, 'e2,same'  , ''         , 'PL', r.kBlue+1 , 1, -1)
    plot_result.addHisto(da_OF_rsfofScaled    , 'hist,SAME', 'data - OF', 'L' , r.kBlue+1 , 1,  1)
    plot_result.addLatex (0.2, 0.80, 'R_{SFOF} scaled')
    plot_result.addLatex (0.2, 0.85, 'max. run {run}'.format(run=maxrun))
    plot_result.saveRatio(1, 0, 0, lint, da_SF, da_OF_rsfofScaled, 0.2, 1.8)
    if returnplot:
        return plot_result

def makePlotsCombinedSR(srlist):
    for i,sR in enumerate(srlist):
        if not i: newSR = copy.deepcopy(sR)
        else:
            newSR = newSR.add(sR)
    print 'this is newSR.mll', newSR.mll
    print 'this is newSR.mll_pred', newSR.mll_pred
    print 'this is newSR.mll_pred.mc', newSR.mll_pred.mc
    makePlots([newSR])

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
    global ingMC, ingDA, onZ, treeDA, treeMC, treeDY, treeTT
    (opts, args) = parser.parse_args()

    print 'running with these options \n'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)
    print ' \n\n'

    ## redp this for 2016 ingMC = helper.ingredients(opts.ingredientFile, 'MC'  )
    ## redp this for 2016 ingDA = helper.ingredients(opts.ingredientFile, 'DATA')

    ## rsfofTable = Tables.makeRSFOFTable(ingDA, ingMC)

    ##print asdf
    print 'Going to load DATA and MC trees...'
    dyDatasets = ['DYJetsToLL_M10to50', 'DYJetsToLL_M50_HT100to200_ext', 'DYJetsToLL_M50_HT100to200_ext', 'DYJetsToLL_M50_HT200to400_ext', 'DYJetsToLL_M50_HT400to600_ext', 'DYJetsToLL_M50_HT600toInf_ext']
    ttDatasets = ['TTJets_DiLepton_total']
    mcDatasets = ttDatasets + ([] if opts.onlyTTbar else dyDatasets)
    daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2_runs_271036_275125', 'DoubleEG_Run2016B-PromptReco-v2_runs_271036_275125', 'MuonEG_Run2016B-PromptReco-v2_runs_271036_275125']

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets, 'DY'), 'DY'  , 0)
    treeTT = Sample.Tree(helper.selectSamples(opts.sampleFile, ttDatasets, 'TT'), 'TT'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print 'Trees successfully loaded...'

    lumi     = 10.
    lumi_str = str(lumi).replace('.', 'p')+'invfb'
    print 'Running with an integrated luminosity of %.2f fb-1' %(lumi)

    isBlinded = True

    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()

    #resultPlot = makeResultData('mll', maxrun = 274240, lint = 0.864, specialcut = '' , scutstring = '', returnplot = True)
    ##makeResultData('nll', maxrun = 274240, lint = 0.864)
    ## resultPlotLoNLL = makeResultData('mll',          274240,        0.864, 'nll_Edge < 21.', 'nllBelow21', returnplot = True)
    ## resultPlotHiNLL = makeResultData('mll',          274240,        0.864, 'nll_Edge > 21.', 'nllAbove21', returnplot = True)
    ## makeResultTable(resultPlotLoNLL, resultPlotHiNLL)
    ## print addsf
    #####print asdfsadf
    ##makeClosureTests('mll','nll_Edge > 21. && run_Edge <= 274240', 'nllAbove21_0p8fb-1')
    #makeClosureTests('nll','lepsMll_Edge > 101. && run_Edge <= 274240', 'highMass_0p8fb-1')

    ##a = makeClosureTests('nll', '', '', True)
    ##makeClosureTests('mll')
    ##makeClosureTests('mll','nll_Edge > 21.', 'nllAbove21')
    ##makeClosureTests('mll','nll_Edge < 21.', 'nllBelow21')
    makeClosureTests('nll', '', '', True)
    ##makeClosureTests('nll','lepsMll_Edge < 81.' , 'lowMass')
    ##makeClosureTests('nll','lepsMll_Edge > 101.', 'highMass')
    print asdfasdfsdf

    finalBinning = range(20,310,10) # + range(81,111,10) + range(110,310,10))

    #dy_shapes = getDYShapes(finalBinning, range(20,302,2), [20., 70., 81., 101., 120., 300.])

    #print asdf  
    regions = []
    loNLL = Region.region('loNLL', 
                           [cuts.AddList([cuts.METJetsSignalRegionMET150, cuts.MaxNLL(20.), '(lepsMll_Edge < 81. || lepsMll_Edge > 101.)']) ],
                           ['mll', 'met', 'nll'],
                           [finalBinning, range(150,310,10), range(10,31,1)],
                           True if not opts.onlyClosure else False)
    hiNLL = Region.region('hiNLL', 
                           [cuts.AddList([cuts.METJetsSignalRegionMET150, cuts.MinNLL(20.), '(lepsMll_Edge < 81. || lepsMll_Edge > 101.)']) ],
                           ['mll', 'met', 'nll'],
                           [finalBinning, range(150,310,10), range(10,31,1)],
                           True if not opts.onlyClosure else False)
    regions.append(loNLL)
    regions.append(hiNLL)

    ## ==================================================
    ## look for DY shapes in file, if not there make them
    ## ==================================================
    dy_shapes = getDYShapes(finalBinning, range(20,302,2), [20., 70., 81., 101., 120., 300.])
    print dy_shapes
    ## ==================================================
    ## done with DY shapes. moving on with life          
    ## ==================================================

    #print asdf

    ## global binnedSR
    ## binnedSR = copy.deepcopy(binnedSR0b)

    for reg in regions:
        print 'i am at region %s' %(reg.name)

        cuts_sf = cuts.AddList([cuts.GoodLeptonSF()]+reg.cuts+[cuts.trigger])
        cuts_of = cuts.AddList([cuts.GoodLeptonOF()]+reg.cuts+[cuts.trigger])

        for tree in ([treeMC, treeDA] if reg.doData else [treeMC]):
            dataMC = 'DATA' if tree == treeDA else 'MC'
            isData = (dataMC == 'DATA')
            for var in reg.rvars:
                v = helper.vParams(var)
                #if not v.name in ['mll']: continue

                setattr(reg, v.name+'_sf'  , tree.getTH1F(lumi, v.name+'_sf_'+reg.name+dataMC, v.vtree, reg.bins[reg.rvars.index(v.name)], 1, 1, cuts_sf, '', v.title) )
                setattr(reg, v.name+'_of'  , tree.getTH1F(lumi, v.name+'_of_'+reg.name+dataMC, v.vtree, reg.bins[reg.rvars.index(v.name)], 1, 1, cuts_of, '', v.title) )
                getattr(reg, v.name        ).setHisto(getattr(reg, v.name+'_sf'), dataMC)
                getattr(reg, v.name+'_pred').setHisto(makePrediction(getattr(reg, v.name+'_of'), ingDA if isData else ingMC), dataMC)

                #reg.mll     .setHisto(reg.mll_sf, dataMC)
                #reg.mll_pred.setHisto(makePrediction(reg.mll_of, ingDA if isData else ingMC), dataMC)

                #del reg.mll_sf, reg.mll_of


    ## make some plots!!
    ## =====================

    #makePlots([loNLL, hiNLL])
    makePlotsCombinedSR([loNLL, hiNLL])


    ## make some tables!!
    ## =====================

    ## a = Tables.makeConciseTable(binnedSRincb, binnedSR0b, binnedSR1b, ingDA, ingMC, onZ)
    #b = Tables.makeConciseTableWith2b(binnedSRincb, binnedSR0b, binnedSR1b, binnedSR2b, ingDA, ingMC, onZ)
    # for i in a: print i
    #for i in a: print i

    ## make datacards
    # Tables.makeDataCards([binnedSRincb, binnedSR0b, binnedSR1b, binnedSR2b], '', ingDA, onZ)

    #a_cen    = makeResultsTable(binnedSR, dy_shapes, 'MC', 'central', opts.nbs)
    #a_fwd    = makeResultsTable(binnedSR, dy_shapes, 'MC', 'forward', opts.nbs)
    #a_cen_da = makeResultsTable(binnedSR, dy_shapes, 'DATA', 'central', opts.nbs)
    #a_fwd_da = makeResultsTable(binnedSR, dy_shapes, 'DATA', 'forward', opts.nbs)
    #print 'CENTRAL TABLE'
    #for i in a_cen_da: print i
    #print 'FORWARD TABLE'
    #for i in a_fwd_da: print i


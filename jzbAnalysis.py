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
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats, TH1F
import math, sys, optparse, copy

import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample


def make_jzbDependencies(histo2d, xvar, opt = 'dep'):

    print '==========================================================='
    print '===== at variable %s ====================================' %(xvar)
    print '==========================================================='

    means = []
    meanes = []
    bin_div = 1
    if xvar == 'met':
        bin_div = 1
    elif xvar == 'zpt':
        bin_div = 1
    elif xvar == 'nvx':
        bin_div = 1
    i = 1

    c = TCanvas()
    c.cd()
    while i+bin_div-1 <= histo2d.GetXaxis().GetNbins():
        tmp_histo = histo2d.ProjectionY(histo2d.GetName()+'_py'+str(i), i, i+bin_div-1, 'e')
        if opt == 'dep': tmp_histo.Sumw2()
        if opt == 'dep':
            fit_tf1 = r.TF1('mygaus','gaus(0)')
            fit_tf1.SetParameter(1, tmp_histo.GetMean())
            fit_tf1.SetParameter(2, 10.)
            tmp_histo.Fit('mygaus')#,'','', -15., 15)
            tmp_histo.Draw()
            c.SaveAs('plots/jzb/controlHistos/'+tmp_histo.GetName()+'.png')
            ff = tmp_histo.GetFunction('mygaus')

            if ff and ff.GetNDF() and ff.GetChisquare()/ff.GetNDF() < 2.5:
                means .append(ff.GetParameter(1))
                meanes.append(ff.GetParError(1))
            else:
                print '==========================================='
                print 'ATTENTION FIT DID NOT CONVERGE OR SOMETHING'
                print '==========================================='
                means .append(tmp_histo.GetMean())
                meanes.append(tmp_histo.GetRMS())

        elif opt == 'response':
            fit_tf1 = r.TF1('mygaus','gaus(0)')
            fit_tf1.SetParLimits(0, 0., 2.*tmp_histo.GetMaximum())
            fit_tf1.SetParameter(1, tmp_histo.GetMean())
            fit_tf1.SetParameter(2,  2.)
            tmp_histo.Fit('mygaus','','',  -5., 5.)
            tmp_histo.Draw()
            c.SaveAs('plots/jzb/controlHistos/'+tmp_histo.GetName()+'.png')
            ff = tmp_histo.GetFunction('mygaus')

            if ff and ff.GetNDF() and ff.GetChisquare()/ff.GetNDF() < 2.5:
                means .append(ff.GetParameter(1))
                meanes.append(ff.GetParError(1))

            else:
                print '==========================================='
                print 'ATTENTION FIT DID NOT CONVERGE OR SOMETHING'
                print '==========================================='
                means .append(tmp_histo.GetMean())
                meanes.append(tmp_histo.GetRMS())

            ## means.append(tmp_histo.GetMean())
            ## meanes.append(tmp_histo.GetRMS())

        i+= bin_div

        del tmp_histo

    h = TH1F(histo2d.GetName(),histo2d.GetName(), len(means), histo2d.GetXaxis().GetXmin(), histo2d.GetXaxis().GetXmax())
    h.GetXaxis().SetTitle(histo2d.GetXaxis().GetTitle())
    
    for mean in means:
        idx = means.index(mean)
        h.SetBinContent(idx+1, mean)
        h.SetBinError(idx+1, meanes[idx])
    if opt == 'response':
        fitfunc = r.TF1('myfunc', ' [0] + [1]/x + [2]/(x*x)')
        fitfunc.SetParLimits(0, -0.1, 0.1)
        #fitfunc.SetParLimits(1, )
        fitfunc.SetParLimits(2,  0., 100.)
        h.Fit('myfunc', '', '', 0., 400.)

    return h

def getVertexCorrection(tree, cuts, eta):
    print '----------------------------------------------------------------'
    print '------ getting vertex correction factors -----------------------'
    print '----------------------------------------------------------------'
    nvxSlope = 0.

    if   eta == 'central':
        etaCut = cuts.Central()
    elif eta == 'forward':
        etaCut = cuts.Forward()

    for var in ['nvx']:#, 'met', 'mll', 'zpt']:

        if   var == 'nvx':
            varTree = 't.lepsJZB_Edge:nVert'
            varName = 'n_{vertices}'
            varBins = [ 5, 10.,  20.]

        elif var == 'mll':
            varTree = 't.lepsJZB_Edge:t.lepsMll_Edge'
            varName = 'm_{ll}'
            varBins = [20, 80., 100.]

        elif var == 'zpt':
            varTree = 't.lepsJZB_Edge:t.lepsZPt_Edge'
            varName = 'm_{ll}'
            varBins = [15,  0., 150.]

        elif var == 'nj':
            varTree = 't.lepsJZB_Edge:t.nJetSel_Edge'
            varName = 'm_{ll}'
            varBins = [10,  0.,  10.]

        jzb_dep  = tree.getTH2F(lumi, 'jzb_vs_'+var+'_'+eta, varTree, varBins[0], varBins[1],  varBins[2], 20, -100., 100., cuts.Add(cuts.DYControlNoMassLeptonSF(), etaCut), '', varName, 'JZB '+eta)
        jzb_dep_proj = make_jzbDependencies(jzb_dep, var)
        jzb_dep_proj.GetYaxis().SetRangeUser(-5., 10.)
        if   var == 'nvx':
            jzb_dep_proj.Fit('pol1', '')
            nvxSlope = jzb_dep_proj.GetFunction('pol1').GetParameter(1)
        plot_jzb_dep_proj = Canvas.Canvas('jzb/plot_jzb_vs_'+var+'_'+eta, 'png,pdf', 0.6, 0.6, 0.8, 0.8)
        plot_jzb_dep_proj.addHisto(jzb_dep_proj, 'E1', 'OF', 'L', r.kBlack, 1, 0)
        plot_jzb_dep_proj.addLine(jzb_dep_proj.GetXaxis().GetXmin(), 0., jzb_dep_proj.GetXaxis().GetXmax(),0., 3)
        plot_jzb_dep_proj.save(0, 0, 0, lumi)

        del plot_jzb_dep_proj, jzb_dep, jzb_dep_proj

    #return the values for central and forward
    return nvxSlope


def getResponseCorrectionFunction(tree, cuts):
    # response correction
    # response_histo  = tree.getTH2F(lumi, 'response_vs_zpt'+'_'+eta, 't.lepsMETRec_Edge/t.lepsZPt_Edge:t.lepsZPt_Edge',  50, 0.,  500., 100, 0., 10., cuts.Add(cuts.DYControlNoMassLeptonSF(), etaCut), '', 'p_{T} (Z) '+eta, 'MET recoil '+eta)
    response_histo  = tree.getTH2F(lumi, 'response_vs_zpt', 't.lepsJZB_Edge/t.lepsZPt_Edge:t.lepsZPt_Edge',  50, 0.,  500., 100, -1.,  2., cuts.AddList([cuts.DYControlNoMassLeptonSF()]), '', 'p_{T} (Z)', 'JZB/(p_{T} (Z))')
    response_histo.GetZaxis().SetRangeUser(0., response_histo.GetMaximum()*1.2)
    response_proj = make_jzbDependencies(response_histo, 'zpt', 'response')
    plot_response = Canvas.Canvas('jzb/response'+'_inclusive', 'png,pdf', 0.6, 0.6, 0.8, 0.8)
    plot_response.addHisto(response_histo, 'COLZ', 'OF', 'L', r.kBlack, 1, 0)
    plot_response.addLine (response_proj.GetXaxis().GetXmin(), 0., response_proj.GetXaxis().GetXmax(),  0., r.kGreen)
    plot_response.addHisto(response_proj , 'PE,SAMES', 'OF', 'L', r.kGray, 1, 0)
    plot_response.save(0, 0, 0, lumi)
    
    responseCorrFunc = copy.deepcopy(response_proj.GetFunction('myfunc'))
    return responseCorrFunc

if __name__ == '__main__':

    print 'Starting JZB analysis...'
    parser = optparse.OptionParser(usage='usage: %prog [options] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    parser.add_option('-t', '--trigger', action='store', type='int', dest='triggerFlag', default='1', help='Trigger cut. Set to 0 if you want to run without trigger')
    (options, args) = parser.parse_args()


    if len(args) != 1:
      parser.error('wrong number of arguments')

    inputFileName = args[0]

    print 'Going to load DATA and MC trees...'
    mcDatasets = ['DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    daDatasets = ['DoubleMuon_Run2015C', 'DoubleEG_Run2015C', 'MuonEG_Run2015C',
                  'DoubleMuon_Run2015D', 'DoubleEG_Run2015D', 'MuonEG_Run2015D']
    treeMC = Sample.Tree(helper.selectSamples(inputFileName, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(inputFileName, daDatasets, 'DA'), 'DATA', 1)
    print 'Trees successfully loaded...'


    global lumi 
    lumi = 0.225

    print 'Running with an integrated luminosity of', lumi,'fb-1'
   
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()

    responseCorrFunc = getResponseCorrectionFunction(treeMC, cuts)
        
    responseCorr = responseCorrFunc.GetParameter(0)

    const = str(responseCorrFunc.GetParameter(0))
    lin   = str(responseCorrFunc.GetParameter(1))
    quad  = str(responseCorrFunc.GetParameter(2))

    responseCorr = '('+const+'+'+lin+'/t.lepsZPt_Edge + '+quad+'/(t.lepsZPt_Edge*t.lepsZPt_Edge) )'
    print 'this is the response correction', responseCorr
    resCor = '(t.lepsZPt_Edge*('+str(responseCorr)+'))'
    rawJZB = 't.lepsJZB_Edge'
    resJZB = rawJZB+'-'+resCor

    response_histo_corr  = treeMC.getTH2F(lumi, 'response_vs_zpt', '('+resJZB+')/t.lepsZPt_Edge:t.lepsZPt_Edge',  50, 0.,  500., 100, -1.,  2., cuts.AddList([cuts.DYControlNoMassLeptonSF()]), '', 'p_{T} (Z)', 'JZB/(p_{T} (Z)) corrected ')
    response_histo_corr.GetZaxis().SetRangeUser(0., response_histo_corr.GetMaximum()*1.2)
    response_proj = make_jzbDependencies(response_histo_corr, 'zpt', 'response')
    plot_response_corr = Canvas.Canvas('jzb/response'+'_inclusive_corrected', 'png,pdf', 0.6, 0.6, 0.8, 0.8)
    plot_response_corr.addHisto(response_histo_corr, 'colz', 'OF', 'L', r.kBlack, 1, 0)
    plot_response_corr.addHisto(response_proj , 'e1 sames', 'OF', 'L', r.kGray, 1, 0)
    plot_response_corr.addLine (response_proj.GetXaxis().GetXmin(), 0., response_proj.GetXaxis().GetXmax(),  0., r.kGreen)
    plot_response_corr.save(0, 0, 0, lumi)

    for eta in ['central', 'forward']:

        if   eta == 'central':
            etaCut = cuts.Central()
        elif eta == 'forward':
            etaCut = cuts.Forward()

        if eta == 'central':
            nvxSlope = 0.
        if eta == 'forward':
            nvxSlope = 0.
        
        ## redo the vertex slope fitting only if initial values are set to 0.
        if nvxSlope == 0.:
            nvxSlope = getVertexCorrection(treeMC, cuts, eta)


        rawJZB = 't.lepsJZB_Edge'
        resCor = '(t.lepsZPt_Edge*('+str(responseCorr)+'))'
        nvxCor = '(nVert*'+str(nvxSlope)+')'
        corJZB = rawJZB+'-'+resCor+'-'+nvxCor
        resJZB = rawJZB+'-'+resCor
        nvxJZB = rawJZB+'-'+nvxCor

        jzb_devel_before = treeMC.getTH1F(lumi, 'jzb_devel_before'+'_'+eta, rawJZB, 40, -60., 60., cuts.Add(cuts.DYControlNoMassLeptonSF(), etaCut), '', 'JZB '+eta)
        jzb_devel_after  = treeMC.getTH1F(lumi, 'jzb_devel_after' +'_'+eta, corJZB, 40, -60., 60., cuts.Add(cuts.DYControlNoMassLeptonSF(), etaCut), '', 'JZB '+eta)
        jzb_devel_resCor = treeMC.getTH1F(lumi, 'jzb_devel_resCor'+'_'+eta, resJZB, 40, -60., 60., cuts.Add(cuts.DYControlNoMassLeptonSF(), etaCut), '', 'JZB '+eta)
        jzb_devel_nvxCor = treeMC.getTH1F(lumi, 'jzb_devel_nvxCor'+'_'+eta, nvxJZB, 40, -60., 60., cuts.Add(cuts.DYControlNoMassLeptonSF(), etaCut), '', 'JZB '+eta)

        ## calculate the residual correction
        print '====================================================================='
        print '============= this is the fit of the residual ======================='
        print '====================================================================='
        jzb_devel_after.Fit('gaus','', '', -15, 15)
        mu = jzb_devel_after.GetFunction('gaus').GetParameter(1)
        print 'residual correction', mu
        print '====================================================================='
        print '====================================================================='
        
        jzb_devel_corrected  = treeMC.getTH1F(lumi, 'jzb_devel_corrected'+'_'+eta , corJZB+'-'+str(mu), 40, -60., 60., cuts.Add(cuts.DYControlNoMassLeptonSF(), etaCut), '', 'JZB '+eta)
        jzb_devel_corrected.Fit('gaus','')
        mufinal = jzb_devel_corrected.GetFunction('gaus').GetParameter(1)
        jzb_devel_corrected.GetFunction('gaus').SetLineColor(r.kMagenta)

        plot_jzb_devel = Canvas.Canvas('jzb/plot_jzb_devel'+'_'+eta, 'png,pdf', 0.4, 0.2, 0.7, 0.4)
        plot_jzb_devel.addHisto(jzb_devel_before    , 'e1'       , 'uncor.'        , 'PEL' , r.kBlack , 1 , 0)
        plot_jzb_devel.addHisto(jzb_devel_resCor    , 'e1 same'  , 'response cor.' , 'PEL' , r.kGray  , 1 , 1)
        plot_jzb_devel.addHisto(jzb_devel_nvxCor    , 'e1 same'  , 'nvx cor.'      , 'PEL' , r.kGreen , 1 , 2)
        plot_jzb_devel.addHisto(jzb_devel_after     , 'e1 same'  , 'resnvx cor.'   , 'PEL' , r.kBlue  , 1 , 3)
        plot_jzb_devel.addHisto(jzb_devel_corrected , 'e1 sames' , 'full cor'      , 'PEL' , r.kRed   , 1 , 4)
        plot_jzb_devel.save(1, 0, 0, lumi)

        jzb_closure_pos = treeMC.getTH1F(lumi, 'jzb_closure_pos'+'_'+eta, '1.*('+corJZB+')' , 20, 0., 60., cuts.AddList([cuts.DYControlNoMassLeptonSF(), etaCut, corJZB+'>0']), '', 'JZB '+eta)
        jzb_closure_neg = treeMC.getTH1F(lumi, 'jzb_closure_neg'+'_'+eta, '-1.*('+corJZB+')', 20, 0., 60., cuts.AddList([cuts.DYControlNoMassLeptonSF(), etaCut, corJZB+'<0']), '', 'JZB '+eta)
        plot_jzb_closure = Canvas.Canvas('jzb/plot_jzb_closure'+'_'+eta, 'png,pdf', 0.7, 0.8, 0.8, 0.9)
        plot_jzb_closure.addHisto(jzb_closure_pos, 'e1 same'  , 'positive JZB', 'PEL' , r.kBlack, 1 , 0)
        plot_jzb_closure.addHisto(jzb_closure_neg, 'e1 same'  , 'negative JZB', 'PEL' , r.kBlue , 1 , 1)
        plot_jzb_closure.saveRatio(1, 0, 0, lumi, jzb_closure_pos, jzb_closure_neg, 0.5, 1.5)

        corMET = '(met_pt + (1.00*'+resCor+'))'
        met_closure_pos = treeMC.getTH1F(lumi, 'met_closure_pos'+'_'+eta, corMET, 20, 0., 60., cuts.AddList([cuts.DYControlNoMassLeptonSF(), etaCut, corJZB+'>0']), '', 'MET '+eta)
        met_closure_neg = treeMC.getTH1F(lumi, 'met_closure_neg'+'_'+eta, corMET, 20, 0., 60., cuts.AddList([cuts.DYControlNoMassLeptonSF(), etaCut, corJZB+'<0']), '', 'MET '+eta)
        plot_met_closure = Canvas.Canvas('jzb/plot_met_closure'+'_'+eta, 'png,pdf', 0.7, 0.8, 0.8, 0.9)
        plot_met_closure.addHisto(met_closure_pos, 'e1 same'  , 'positive JZB', 'PEL' , r.kBlack, 1 , 0)
        plot_met_closure.addHisto(met_closure_neg, 'e1 same'  , 'negative JZB', 'PEL' , r.kBlue , 1 , 1)
        plot_met_closure.saveRatio(1, 0, 0, lumi, met_closure_pos, met_closure_neg, 0.5, 1.5)


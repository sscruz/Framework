import ROOT as r
import math as math
import copy

from optparse import OptionParser
from ROOT import gROOT, TCanvas, TFile, TF1, TH1F, TLine
from Sample import Sample, Block, Tree
from CutManager import CutManager
from Canvas import Canvas


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
            c.SaveAs('plots/controlHistos/'+tmp_histo.GetName()+'.png')
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
            c.SaveAs('plots/controlHistos/'+tmp_histo.GetName()+'.png')
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

        jzb_dep  = tree.getTH2F(4, 'jzb_vs_'+var+'_'+eta, varTree, varBins[0], varBins[1],  varBins[2], 20, -100., 100., cuts.Add(cuts.DYControlNoMassLeptonSF(), etaCut), '', varName, 'JZB '+eta)
        jzb_dep_proj = make_jzbDependencies(jzb_dep, var)
        jzb_dep_proj.GetYaxis().SetRangeUser(-5., 10.)
        if   var == 'nvx':
            jzb_dep_proj.Fit('pol1', '')
            nvxSlope = jzb_dep_proj.GetFunction('pol1').GetParameter(1)
        plot_jzb_dep_proj = Canvas('plot_jzb_vs_'+var+'_'+eta, 'png,pdf', 0.6, 0.6, 0.8, 0.8)
        plot_jzb_dep_proj.addHisto(jzb_dep_proj, 'E1', 'OF', 'L', r.kBlack, 1, 0)
        plot_jzb_dep_proj.addLine(jzb_dep_proj.GetXaxis().GetXmin(), 0., jzb_dep_proj.GetXaxis().GetXmax(),0., 3)
        plot_jzb_dep_proj.save(0, 0, 0, 4.0)

        del plot_jzb_dep_proj, jzb_dep, jzb_dep_proj

    #return the values for central and forward
    return nvxSlope


def getResponseCorrectionFunction(tree, cuts):
    # response correction
    # response_histo  = tree.getTH2F(4, 'response_vs_zpt'+'_'+eta, 't.lepsMETRec_Edge/t.lepsZPt_Edge:t.lepsZPt_Edge',  50, 0.,  500., 100, 0., 10., cuts.Add(cuts.DYControlNoMassLeptonSF(), etaCut), '', 'p_{T} (Z) '+eta, 'MET recoil '+eta)
    response_histo  = tree.getTH2F(4, 'response_vs_zpt', 't.lepsJZB_Edge/t.lepsZPt_Edge:t.lepsZPt_Edge',  50, 0.,  500., 100, -1.,  2., cuts.AddList([cuts.DYControlNoMassLeptonSF()]), '', 'p_{T} (Z)', 'JZB/(p_{T} (Z))')
    response_proj = make_jzbDependencies(response_histo, 'zpt', 'response')
    plot_response = Canvas('response'+'_inclusive', 'png,pdf', 0.6, 0.6, 0.8, 0.8)
    plot_response.addHisto(response_histo, 'colz', 'OF', 'L', r.kBlack, 1, 0)
    plot_response.addHisto(response_proj , 'e1 sames', 'OF', 'L', r.kGray, 1, 0)
    plot_response.save(0, 0, 0, 1.0)
    
    responseCorrFunc = copy.deepcopy(response_proj.GetFunction('myfunc'))
    return responseCorrFunc

if __name__ == '__main__':

    parser = OptionParser(usage='usage: %prog [options] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    (options, args) = parser.parse_args()

    if len(args) != 1:
      parser.error('wrong number of arguments')

    inputFileName = args[0]
    tree = Tree(inputFileName, 'MC', 0)
   
    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    r.gStyle.SetStatX(0.7)
    r.gStyle.SetStatY(0.85)
    r.gStyle.SetStatH(0.1)
    r.gStyle.SetStatW(0.2)

    r.gStyle.SetPalette(53)
    r.gStyle.SetNumberContours( 999 )

    cuts = CutManager()

    responseCorrFunc = getResponseCorrectionFunction(tree, cuts)
        
    responseCorr = responseCorrFunc.GetParameter(0)

    const = str(responseCorrFunc.GetParameter(0))
    lin   = str(responseCorrFunc.GetParameter(1))
    quad  = str(responseCorrFunc.GetParameter(2))

    responseCorr = '('+const+'+'+lin+'/t.lepsZPt_Edge + '+quad+'/(t.lepsZPt_Edge*t.lepsZPt_Edge) )'
    print 'this is the response correction', responseCorr
    resCor = '(t.lepsZPt_Edge*('+str(responseCorr)+'))'
    rawJZB = 't.lepsJZB_Edge'
    resJZB = rawJZB+'-'+resCor

    response_histo_corr  = tree.getTH2F(4, 'response_vs_zpt', '('+resJZB+')/t.lepsZPt_Edge:t.lepsZPt_Edge',  50, 0.,  500., 100, -1.,  2., cuts.AddList([cuts.DYControlNoMassLeptonSF()]), '', 'p_{T} (Z)', 'JZB/(p_{T} (Z)) corrected ')
    response_proj = make_jzbDependencies(response_histo_corr, 'zpt', 'response')
    plot_response_corr = Canvas('response'+'_inclusive_corrected', 'png,pdf', 0.6, 0.6, 0.8, 0.8)
    plot_response_corr.addHisto(response_histo_corr, 'colz', 'OF', 'L', r.kBlack, 1, 0)
    plot_response_corr.addHisto(response_proj , 'e1 sames', 'OF', 'L', r.kGray, 1, 0)
    plot_response_corr.save(0, 0, 0, 1.0)

    for eta in ['central', 'forward']:

        if   eta == 'central':
            etaCut = cuts.Central()
        elif eta == 'forward':
            etaCut = cuts.Forward()

        if eta == 'central':
            nvxSlope = 0.135
        if eta == 'forward':
            nvxSlope = 0.189
        
        ## redo the vertex slope fitting only if initial values are set to 0.
        if nvxSlope == 0.:
            nvxSlope = getVertexCorrection(tree, cuts, eta)


        rawJZB = 't.lepsJZB_Edge'
        resCor = '(t.lepsZPt_Edge*('+str(responseCorr)+'))'
        nvxCor = '(nVert*'+str(nvxSlope)+')'
        corJZB = rawJZB+'-'+resCor+'-'+nvxCor
        resJZB = rawJZB+'-'+resCor
        nvxJZB = rawJZB+'-'+nvxCor

        jzb_devel_before = tree.getTH1F(4, 'jzb_devel_before'+'_'+eta, rawJZB, 40, -60., 60., cuts.Add(cuts.DYControlNoMassLeptonSF(), etaCut), '', 'JZB '+eta)
        jzb_devel_after  = tree.getTH1F(4, 'jzb_devel_after' +'_'+eta, corJZB, 40, -60., 60., cuts.Add(cuts.DYControlNoMassLeptonSF(), etaCut), '', 'JZB '+eta)
        jzb_devel_resCor = tree.getTH1F(4, 'jzb_devel_resCor'+'_'+eta, resJZB, 40, -60., 60., cuts.Add(cuts.DYControlNoMassLeptonSF(), etaCut), '', 'JZB '+eta)
        jzb_devel_nvxCor = tree.getTH1F(4, 'jzb_devel_nvxCor'+'_'+eta, nvxJZB, 40, -60., 60., cuts.Add(cuts.DYControlNoMassLeptonSF(), etaCut), '', 'JZB '+eta)

        ## calculate the residual correction
        print '====================================================================='
        print '============= this is the fit of the residual ======================='
        print '====================================================================='
        jzb_devel_after.Fit('gaus','', '', -15, 15)
        mu = jzb_devel_after.GetFunction('gaus').GetParameter(1)
        print 'residual correction', mu
        print '====================================================================='
        print '====================================================================='
        
        jzb_devel_corrected  = tree.getTH1F(4, 'jzb_devel_corrected'+'_'+eta , corJZB+'-'+str(mu), 40, -60., 60., cuts.Add(cuts.DYControlNoMassLeptonSF(), etaCut), '', 'JZB '+eta)
        jzb_devel_corrected.Fit('gaus','')
        mufinal = jzb_devel_corrected.GetFunction('gaus').GetParameter(1)
        jzb_devel_corrected.GetFunction('gaus').SetLineColor(r.kMagenta)

        plot_jzb_devel = Canvas('plot_jzb_devel'+'_'+eta, 'png,pdf', 0.4, 0.2, 0.7, 0.4)
        plot_jzb_devel.addHisto(jzb_devel_before    , 'e1'       , 'uncor.'        , 'PEL' , r.kBlack , 1 , 0)
        plot_jzb_devel.addHisto(jzb_devel_resCor    , 'e1 same'  , 'response cor.' , 'PEL' , r.kGray  , 1 , 1)
        plot_jzb_devel.addHisto(jzb_devel_nvxCor    , 'e1 same'  , 'nvx cor.'      , 'PEL' , r.kGreen , 1 , 2)
        plot_jzb_devel.addHisto(jzb_devel_after     , 'e1 same'  , 'resnvx cor.'   , 'PEL' , r.kBlue  , 1 , 3)
        plot_jzb_devel.addHisto(jzb_devel_corrected , 'e1 sames' , 'full cor'      , 'PEL' , r.kRed   , 1 , 4)
        plot_jzb_devel.save(1, 0, 0, 4.0)

        jzb_closure_pos = tree.getTH1F(4, 'jzb_closure_pos'+'_'+eta, '1.*('+corJZB+')' , 20, 0., 60., cuts.AddList([cuts.DYControlNoMassLeptonSF(), etaCut, corJZB+'>0']), '', 'JZB '+eta)
        jzb_closure_neg = tree.getTH1F(4, 'jzb_closure_neg'+'_'+eta, '-1.*('+corJZB+')', 20, 0., 60., cuts.AddList([cuts.DYControlNoMassLeptonSF(), etaCut, corJZB+'<0']), '', 'JZB '+eta)
        plot_jzb_closure = Canvas('plot_jzb_closure'+'_'+eta, 'png,pdf', 0.7, 0.8, 0.8, 0.9)
        plot_jzb_closure.addHisto(jzb_closure_pos, 'e1 same'  , 'positive JZB', 'PEL' , r.kBlack, 1 , 0)
        plot_jzb_closure.addHisto(jzb_closure_neg, 'e1 same'  , 'negative JZB', 'PEL' , r.kBlue , 1 , 1)
        plot_jzb_closure.saveRatio(1, 0, 0, 4.0, jzb_closure_pos, jzb_closure_neg, 0.5, 1.5)

        corMET = '(met_pt + (1.00*'+resCor+'))'
        met_closure_pos = tree.getTH1F(4, 'met_closure_pos'+'_'+eta, corMET, 20, 0., 60., cuts.AddList([cuts.DYControlNoMassLeptonSF(), etaCut, corJZB+'>0']), '', 'MET '+eta)
        met_closure_neg = tree.getTH1F(4, 'met_closure_neg'+'_'+eta, corMET, 20, 0., 60., cuts.AddList([cuts.DYControlNoMassLeptonSF(), etaCut, corJZB+'<0']), '', 'MET '+eta)
        plot_met_closure = Canvas('plot_met_closure'+'_'+eta, 'png,pdf', 0.7, 0.8, 0.8, 0.9)
        plot_met_closure.addHisto(met_closure_pos, 'e1 same'  , 'positive JZB', 'PEL' , r.kBlack, 1 , 0)
        plot_met_closure.addHisto(met_closure_neg, 'e1 same'  , 'negative JZB', 'PEL' , r.kBlue , 1 , 1)
        plot_met_closure.saveRatio(1, 0, 0, 4.0, met_closure_pos, met_closure_neg, 0.5, 1.5)


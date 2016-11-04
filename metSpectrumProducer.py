import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TH1F, TGraph
import math, sys, optparse, copy, re, array, itertools, random


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables

import metSpectrumProducer

colorlist = [r.kBlack   , r.kRed    , r.kGreen , r.kBlue  , r.kGray    , r.kOrange  ,
             r.kMagenta , r.kCyan   , r.kRed-3 , r.kYellow, r.kGreen-6 , r.kOrange-5, 
             r.kSpring-7, r.kAzure-6, r.kGray+2, 8        , r.kOrange+3, r.kBlue-6  ,
             r.kRed+3   , r.kAzure+3, r.kGreen+3, r.kYellow+3, r.kBlue+2, r.kMagenta+2 ] ## 24 colors

def getColor():
    global colorlist
    retcolor = colorlist[0]
    colorlist = colorlist[1:]
    return retcolor


def makeRSFOFPlots(ttt, regionCut):
    cuts = CutManager.CutManager() ## True means ewino. adapts the pts and the dR. and requires == 2 tight leptons
    regionCut = cuts.AddList([cuts.GoodLeptonSF(), 'nBJetMedium25_Edge >= 0', cuts.nj2, 'met_Edge > 50.'])
    bcof = regionCut.replace('121', '143')
    bcof = bcof.replace('169', '143')
    r.gROOT.SetBatch()
    lumi = 7.65
    for var in ['met']:
        v = helper.vParams(var)
        sfhisto = ttt.getTH1F(lumi, 'sfhisto', v.vtree, v.nbins, v.xmin, v.xmax, regionCut, '', v.title, ofBin = True)
        ofhisto = ttt.getTH1F(lumi, 'ofhisto', v.vtree, v.nbins, v.xmin, v.xmax, bcof     , '', v.title, ofBin = True)
    
        div = copy.deepcopy(sfhisto)
        div.Divide(ofhisto)
        div.GetYaxis().SetRangeUser(0.5, 1.5)
        div.SetMarkerStyle(20)
        div.SetMarkerSize (0.8)
        div.SetMarkerColor(r.kRed+1)

        div.GetYaxis().SetTitle('r_{SF/OF}')
        
        c_tmp = r.TCanvas()
        div.Fit('pol0')
        div.Draw('pe')
        c_tmp.SaveAs('plots/ewkino/rsfof/'+var+'_rsfof.pdf')
        c_tmp.SaveAs('plots/ewkino/rsfof/'+var+'_rsfof.png')

class metObject:
    def __init__(self, cutlist, bkgTrees, daTree, onlyTT):
        self.cl = cutlist
        self.allHistos = []
        self.bkg  = bkgTrees
        self.data = daTree
        self.treett =  filter(lambda x: x.name == 'TT', self.bkg)
        self.onlyTT = onlyTT
        self.allHistos = []
        #self.getHistos('met')

    def getHistos(self, var):
        met_v = helper.vParams(var)
        for tmp_tree in (self.bkg if not self.onlyTT else self.treett):
            for ci, (cn, cut) in enumerate(self.cl.items()):
                tmp_h = tmp_tree.getTH1F(lumi, cn+'_%s_%s'%(tmp_tree.name, var), met_v.vtree, met_v.nbins, met_v.xmin, met_v.xmax, cut, '', met_v.title, ofBin = True)
                tmp_color = getColor()
                tmp_h.SetMarkerColor(tmp_color)
                tmp_h.SetLineColor  (tmp_color)
                tmp_h.SetLineWidth  (2)
                setattr(self, cn+'_%s_%s'%(tmp_tree.name, var), copy.deepcopy(tmp_h))
                if tmp_h.Integral(): tmp_h.Scale(1./tmp_h.Integral())
                tmp_h.SetName(tmp_h.GetName()+'_norm')
                setattr(self, cn+'_%s_%s_%s'%(tmp_tree.name, var, 'norm'), copy.deepcopy(tmp_h))
                self.allHistos.append(getattr(self, cn+'_%s_%s'   %(tmp_tree.name, var)        ))
                self.allHistos.append(getattr(self, cn+'_%s_%s_%s'%(tmp_tree.name, var, 'norm')))
        for ci, (cn, cut) in enumerate(self.cl.items()):
            tmp_h_da = self.data.getTH1F(lumi, cn+'_DA', met_v.vtree, met_v.nbins, met_v.xmin, met_v.xmax, cut, '', met_v.title, ofBin = True)
            tmp_style = 20 if 'sf' in cn and '0b' in cn else 21 if 'of' in cn and '0b' in cn else 24 if 'of' in cn and '0b' in cn else 25
            tmp_color = r.kBlack if 'sf' in cn else r.kGray
            tmp_h_da.SetLineColor  (tmp_color)
            tmp_h_da.SetMarkerColor(tmp_color)
            tmp_h_da.SetMarkerStyle(tmp_style)
            setattr(self, cn+'_DA_%s'%(var), copy.deepcopy(tmp_h_da))
            if tmp_h_da.Integral(): tmp_h_da.Scale(1./tmp_h_da.Integral())
            tmp_h_da.SetName(tmp_h_da.GetName()+'_norm')
            setattr(self, cn+'_DA_%s_%s'%(var, 'norm'), copy.deepcopy(tmp_h_da))
            self.allHistos.append(getattr(self, cn+'_DA_%s'   %(var)   ))
            self.allHistos.append(getattr(self, cn+'_DA_%s_%s'%(var, 'norm')))

    def makeMETPlotSimple(self):
        ratio_h = []
        plot = Canvas.Canvas('ewkino/met_comparison_TT', 'png,pdf', 0.6, 0.50, 0.85, 0.82)
        _i = 0
        histolist = filter(lambda x: 'norm' in x.GetName() and not any(i in x.GetName() for i in ['mt285', 'mt20', 'DA']), self.allHistos)
        histolist_da = filter(lambda x: 'norm' in x.GetName() and not any(i in x.GetName() for i in ['mt285', 'mt20', 'TT', 'sf']), self.allHistos)
        for _i,_h in enumerate(histolist):
            ratio_h.append(_h)
            plot.addHisto(_h, 'hist same', _h.GetName().replace('_',' '), 'L' , _h.GetLineColor(), 1,  _i); _i+=1;
        for _i,_h in enumerate(histolist_da):
            ratio_h.append(_h)
            plot.addHisto(_h, 'pe same', _h.GetName().replace('_',' '), 'PL' , _h.GetLineColor(), 1,  _i); _i+=1;
        plot.saveRatio(1, 0, 0 , lumi, self.sf_0b_nomt2_TT_norm, ratio_h, 0.0, 2.0)
        plot.saveRatio(1, 0, 1 , lumi, self.sf_0b_nomt2_TT_norm, ratio_h, 0.0, 2.0)
        
def controlPlots(mo):
    tt_of_1b = copy.deepcopy(mo.of_1b_mt280_TT_met)
    tt_sf_1b = copy.deepcopy(mo.sf_1b_mt280_TT_met)
    tt_sf_0b = copy.deepcopy(mo.sf_0b_mt280_TT_met)
    tt_of_0b = copy.deepcopy(mo.of_0b_mt280_TT_met)
    helper.makeRatioPlot([tt_sf_0b,tt_of_0b,tt_sf_1b,tt_of_1b], 'plots/ewkino/controlPlots/metComparisons', {'norm': True, 'legco': [0.65,0.60,0.88,0.90]})

def ratioError(num, num_e, den, den_e, opt=''):
    tmp_n = r.TH1F('n', 'n', 1,0,1)
    tmp_d = r.TH1F('d', 'd', 1,0,1)
    tmp_r = r.TH1F('r', 'r', 1,0,1)
    tmp_n.SetBinContent(1,num); tmp_n.SetBinError(1,num_e)
    tmp_d.SetBinContent(1,den); tmp_d.SetBinError(1,den_e)
    #tmp_n.Divide(tmp_d)
    tmp_r.Divide(tmp_n, tmp_d, 1., 1., opt)
    return tmp_r.GetBinContent(1), tmp_r.GetBinError(1)

def getFmll(mo,withData=False):
    tt_of_1b = copy.deepcopy(mo.of_1b_mt280_61mll121_TT_mll)
    tt_sf_1b = copy.deepcopy(mo.sf_1b_mt280_61mll121_TT_mll)
    tt_sf_0b = copy.deepcopy(mo.sf_0b_mt280_61mll121_TT_mll)
    tt_of_0b = copy.deepcopy(mo.of_0b_mt280_61mll121_TT_mll)
    da_of_1b = copy.deepcopy(mo.of_1b_mt280_61mll121_DA_mll)
    da_sf_1b = copy.deepcopy(mo.sf_1b_mt280_61mll121_DA_mll)
    da_sf_0b = copy.deepcopy(mo.sf_0b_mt280_61mll121_DA_mll)
    da_of_0b = copy.deepcopy(mo.of_0b_mt280_61mll121_DA_mll)


    bin61  = tt_of_1b.FindBin(61.)
    bin81  = tt_of_1b.FindBin(81.)
    bin101 = tt_of_1b.FindBin(101.)
    bin121 = tt_of_1b.FindBin(121.)

    tt_frac_of_1b_eN = r.Double(-1.); tt_frac_of_1b_eD = r.Double(-1.)
    tt_frac_sf_1b_eN = r.Double(-1.); tt_frac_sf_1b_eD = r.Double(-1.)
    tt_frac_of_0b_eN = r.Double(-1.); tt_frac_of_0b_eD = r.Double(-1.)
    tt_frac_sf_0b_eN = r.Double(-1.); tt_frac_sf_0b_eD = r.Double(-1.)
    tt_frac_of_1b_N = tt_of_1b.IntegralAndError(bin81, bin101-1, tt_frac_of_1b_eN); tt_frac_of_1b_D = tt_of_1b.IntegralAndError(bin61, bin121-1, tt_frac_of_1b_eD)
    tt_frac_sf_1b_N = tt_sf_1b.IntegralAndError(bin81, bin101-1, tt_frac_sf_1b_eN); tt_frac_sf_1b_D = tt_sf_1b.IntegralAndError(bin61, bin121-1, tt_frac_sf_1b_eD)
    tt_frac_of_0b_N = tt_of_0b.IntegralAndError(bin81, bin101-1, tt_frac_of_0b_eN); tt_frac_of_0b_D = tt_of_0b.IntegralAndError(bin61, bin121-1, tt_frac_of_0b_eD)
    tt_frac_sf_0b_N = tt_sf_0b.IntegralAndError(bin81, bin101-1, tt_frac_sf_0b_eN); tt_frac_sf_0b_D = tt_sf_0b.IntegralAndError(bin61, bin121-1, tt_frac_sf_0b_eD)
    
    da_frac_of_1b_eN = r.Double(-1.); da_frac_of_1b_eD = r.Double(-1.)
    da_frac_sf_1b_eN = r.Double(-1.); da_frac_sf_1b_eD = r.Double(-1.)
    da_frac_of_0b_eN = r.Double(-1.); da_frac_of_0b_eD = r.Double(-1.)
    da_frac_sf_0b_eN = r.Double(-1.); da_frac_sf_0b_eD = r.Double(-1.)
    da_frac_of_1b_N = da_of_1b.IntegralAndError(bin81, bin101-1, da_frac_of_1b_eN); da_frac_of_1b_D = da_of_1b.IntegralAndError(bin61, bin121-1, da_frac_of_1b_eD)
    da_frac_sf_1b_N = da_sf_1b.IntegralAndError(bin81, bin101-1, da_frac_sf_1b_eN); da_frac_sf_1b_D = da_sf_1b.IntegralAndError(bin61, bin121-1, da_frac_sf_1b_eD)
    da_frac_of_0b_N = da_of_0b.IntegralAndError(bin81, bin101-1, da_frac_of_0b_eN); da_frac_of_0b_D = da_of_0b.IntegralAndError(bin61, bin121-1, da_frac_of_0b_eD)
    da_frac_sf_0b_N = da_sf_0b.IntegralAndError(bin81, bin101-1, da_frac_sf_0b_eN); da_frac_sf_0b_D = da_sf_0b.IntegralAndError(bin61, bin121-1, da_frac_sf_0b_eD)
    
    print '======= ttbar ====================='
    print 'numerator: frac_of_1b: %.3f +- %.3f'%( tt_frac_of_1b_N , tt_frac_of_1b_eN)
    print 'numerator: frac_sf_1b: %.3f +- %.3f'%( tt_frac_sf_1b_N , tt_frac_sf_1b_eN)
    print 'numerator: frac_of_0b: %.3f +- %.3f'%( tt_frac_of_0b_N , tt_frac_of_0b_eN)
    print 'numerator: frac_sf_0b: %.3f +- %.3f'%( tt_frac_sf_0b_N , tt_frac_sf_0b_eN)
    print 'denom    : frac_of_1b: %.3f +- %.3f'%( tt_frac_of_1b_D , tt_frac_of_1b_eD)
    print 'denom    : frac_sf_1b: %.3f +- %.3f'%( tt_frac_sf_1b_D , tt_frac_sf_1b_eD)
    print 'denom    : frac_of_0b: %.3f +- %.3f'%( tt_frac_of_0b_D , tt_frac_of_0b_eD)
    print 'denom    : frac_sf_0b: %.3f +- %.3f'%( tt_frac_sf_0b_D , tt_frac_sf_0b_eD)

    tt_ratio_of_1b = ratioError(tt_frac_of_1b_N , tt_frac_of_1b_eN, tt_frac_of_1b_D , tt_frac_of_1b_eD, 'B')   
    tt_ratio_sf_1b = ratioError(tt_frac_sf_1b_N , tt_frac_sf_1b_eN, tt_frac_sf_1b_D , tt_frac_sf_1b_eD, 'B')
    tt_ratio_of_0b = ratioError(tt_frac_of_0b_N , tt_frac_of_0b_eN, tt_frac_of_0b_D , tt_frac_of_0b_eD, 'B')
    tt_ratio_sf_0b = ratioError(tt_frac_sf_0b_N , tt_frac_sf_0b_eN, tt_frac_sf_0b_D , tt_frac_sf_0b_eD, 'B')
    
    print 'ratio: frac_of_1b: %.3f +- %.3f'%( tt_ratio_of_1b[0], tt_ratio_of_1b[1])
    print 'ratio: frac_sf_1b: %.3f +- %.3f'%( tt_ratio_sf_1b[0], tt_ratio_sf_1b[1])
    print 'ratio: frac_of_0b: %.3f +- %.3f'%( tt_ratio_of_0b[0], tt_ratio_of_0b[1])
    print 'ratio: frac_sf_0b: %.3f +- %.3f'%( tt_ratio_sf_0b[0], tt_ratio_sf_0b[1])

    print '======= data ======================'
    print 'numerator: frac_of_1b: %.3f +- %.3f'%( da_frac_of_1b_N , da_frac_of_1b_eN)
    print 'numerator: frac_sf_1b: %.3f +- %.3f'%( da_frac_sf_1b_N , da_frac_sf_1b_eN)
    print 'numerator: frac_of_0b: %.3f +- %.3f'%( da_frac_of_0b_N , da_frac_of_0b_eN)
    print 'denom    : frac_of_1b: %.3f +- %.3f'%( da_frac_of_1b_D , da_frac_of_1b_eD)
    print 'denom    : frac_sf_1b: %.3f +- %.3f'%( da_frac_sf_1b_D , da_frac_sf_1b_eD)
    print 'denom    : frac_of_0b: %.3f +- %.3f'%( da_frac_of_0b_D , da_frac_of_0b_eD)

    da_ratio_of_1b = ratioError(da_frac_of_1b_N , da_frac_of_1b_eN, da_frac_of_1b_D , da_frac_of_1b_eD, 'B')   
    da_ratio_sf_1b = ratioError(da_frac_sf_1b_N , da_frac_sf_1b_eN, da_frac_sf_1b_D , da_frac_sf_1b_eD, 'B')
    da_ratio_of_0b = ratioError(da_frac_of_0b_N , da_frac_of_0b_eN, da_frac_of_0b_D , da_frac_of_0b_eD, 'B')
    
    print 'ratio: frac_of_1b: %.3f +- %.3f'%( da_ratio_of_1b[0], da_ratio_of_1b[1])
    print 'ratio: frac_sf_1b: %.3f +- %.3f'%( da_ratio_sf_1b[0], da_ratio_sf_1b[1])
    print 'ratio: frac_of_0b: %.3f +- %.3f'%( da_ratio_of_0b[0], da_ratio_of_0b[1])

    c1 = r.TCanvas()
    l1 = r.TLegend(0.4,0.8,0.6,0.89)
    c1.SetLeftMargin(0.15)
    c1.SetRightMargin(0.05)
    tmp_h = r.TH1F('fracs','f_{mll}',3,0,3)
    tmp_h.SetBinContent(1, tt_ratio_of_1b[0]); tmp_h.SetBinError(1, tt_ratio_of_1b[1]); tmp_h.GetXaxis().SetBinLabel(1, 'OF #geq 1b');
    tmp_h.SetBinContent(2, tt_ratio_of_0b[0]); tmp_h.SetBinError(2, tt_ratio_of_0b[1]); tmp_h.GetXaxis().SetBinLabel(2, 'OF = 0b');
    tmp_h.SetBinContent(3, tt_ratio_sf_0b[0]); tmp_h.SetBinError(3, tt_ratio_sf_0b[1]); tmp_h.GetXaxis().SetBinLabel(3, 'SF = 0b');
    #tmp_h.SetBinContent(2, tt_ratio_sf_1b[0]); tmp_h.SetBinError(2, tt_ratio_sf_1b[1]); tmp_h.GetXaxis().SetBinLabel(2, 'SF #geq 1b');
    tmp_h.GetYaxis().SetRangeUser(0.10, 0.60)
    tmp_h.GetYaxis().SetTitle('f_{mll}')
    tmp_h.GetYaxis().SetTitleOffset(1.2*tmp_h.GetYaxis().GetTitleOffset())
    tmp_h.SetMarkerColor(r.kBlue+2)
    tmp_h.SetLineColor  (r.kBlue-7)
    tmp_h.SetLineWidth  (2)
    tmp_h.SetMarkerSize (1.2*tmp_h.GetMarkerSize())
    r.gStyle.SetPaintTextFormat('.3f')
    tmp_h.Draw('e1')
    tmp_h.Draw('text00 same')
    l1.AddEntry(tmp_h, 'ttbar MC', 'pl')
    if withData:
        tmp_h_da = tmp_h.Clone('data fracs')
        tmp_h_da.SetBinContent(1, da_ratio_of_1b[0]); tmp_h_da.SetBinError(1, da_ratio_of_1b[1]); tmp_h_da.GetXaxis().SetBinLabel(1, 'OF #geq 1b');
        tmp_h_da.SetBinContent(2, da_ratio_of_0b[0]); tmp_h_da.SetBinError(2, da_ratio_of_0b[1]); tmp_h_da.GetXaxis().SetBinLabel(2, 'OF = 0b');
        tmp_h_da.SetBinContent(3, 0.               ); tmp_h_da.SetBinError(3, 0.               ); tmp_h_da.GetXaxis().SetBinLabel(3, 'SF = 0b');
        #tmp_h_da.SetBinContent(2, da_ratio_sf_1b[0]); tmp_h_da.SetBinError(2, da_ratio_sf_1b[1]); tmp_h_da.GetXaxis().SetBinLabel(2, 'SF #geq 1b');
        tmp_h_da.SetMarkerStyle(22)
        tmp_h_da.SetMarkerColor(r.kBlack)
        tmp_h_da.SetLineColor  (r.kGray+1)
        tmp_h_da.SetLineWidth  (2)
        tmp_h_da.Draw('e1 ')
        tmp_h_da.Draw('text00 same ')
        tmp_h.Draw('e1 same')
        tmp_h.Draw('text00 same')
        l1.AddEntry(tmp_h_da, 'data', 'pl')
    l1.Draw('same')
    lat = r.TLatex(); lat.SetNDC(); lat.SetTextFont(tmp_h.GetYaxis().GetTitleFont()); lat.SetTextSize(tmp_h.GetYaxis().GetTitleSize())
    lumi = 12.9
    lumistr = 'lumi12p9invfb'
    lat.DrawLatex(0.76, 0.93, '%.2f fb^{-1}'%(lumi))
    #lat.Draw('same')
    c1.SaveAs('plots/ttbarMET/%s/fmll_highStats%s.png' %(lumistr, '_withData' if withData else ''))
    c1.SaveAs('plots/ttbarMET/%s/fmll_highStats%s.pdf' %(lumistr, '_withData' if withData else ''))
    c1.SaveAs('plots/ttbarMET/%s/fmll_highStats%s.root'%(lumistr, '_withData' if withData else ''))

    ## calculate also r0b1b
    print '================================================== '
    print '============================================ r0b1b '
    print '================================================== '
    tt_r0b1b_sf_incM = ratioError(tt_frac_sf_0b_D, tt_frac_sf_0b_eD, tt_frac_sf_1b_D, tt_frac_sf_1b_eD) 
    tt_r0b1b_of_incM = ratioError(tt_frac_of_0b_D, tt_frac_of_0b_eD, tt_frac_of_1b_D, tt_frac_of_1b_eD) 
    tt_r0b1b_sf_onZ  = ratioError(tt_frac_sf_0b_N, tt_frac_sf_0b_eN, tt_frac_sf_1b_N, tt_frac_sf_1b_eN) 
    tt_r0b1b_of_onZ  = ratioError(tt_frac_of_0b_N, tt_frac_of_0b_eN, tt_frac_of_1b_N, tt_frac_of_1b_eN) 
    print '\n\n ttbar ==================================='
    print 'ratio 0b1b: sf_incM: %.3f +- %.3f (%.3f +- %.3f / %.3f +- %.3f)'%( tt_r0b1b_sf_incM[0], tt_r0b1b_sf_incM[1], tt_frac_sf_0b_D, tt_frac_sf_0b_eD, tt_frac_sf_1b_D, tt_frac_sf_1b_eD)
    print 'ratio 0b1b: of_incM: %.3f +- %.3f (%.3f +- %.3f / %.3f +- %.3f)'%( tt_r0b1b_of_incM[0], tt_r0b1b_of_incM[1], tt_frac_of_0b_D, tt_frac_of_0b_eD, tt_frac_of_1b_D, tt_frac_of_1b_eD)
    print 'ratio 0b1b: sf_onZ : %.3f +- %.3f (%.3f +- %.3f / %.3f +- %.3f)'%( tt_r0b1b_sf_onZ [0], tt_r0b1b_sf_onZ [1], tt_frac_sf_0b_N, tt_frac_sf_0b_eN, tt_frac_sf_1b_N, tt_frac_sf_1b_eN)
    print 'ratio 0b1b: of_onZ : %.3f +- %.3f (%.3f +- %.3f / %.3f +- %.3f)'%( tt_r0b1b_of_onZ [0], tt_r0b1b_of_onZ [1], tt_frac_of_0b_N, tt_frac_of_0b_eN, tt_frac_of_1b_N, tt_frac_of_1b_eN)

    da_r0b1b_sf_incM = ratioError(da_frac_sf_0b_D, da_frac_sf_0b_eD, da_frac_sf_1b_D, da_frac_sf_1b_eD) 
    da_r0b1b_of_incM = ratioError(da_frac_of_0b_D, da_frac_of_0b_eD, da_frac_of_1b_D, da_frac_of_1b_eD) 
    da_r0b1b_sf_onZ  = ratioError(da_frac_sf_0b_N, da_frac_sf_0b_eN, da_frac_sf_1b_N, da_frac_sf_1b_eN) 
    da_r0b1b_of_onZ  = ratioError(da_frac_of_0b_N, da_frac_of_0b_eN, da_frac_of_1b_N, da_frac_of_1b_eN) 
    print '\n\n data ===================================='
    print 'ratio 0b1b: sf_incM: %.3f +- %.3f (%.3f +- %.3f / %.3f +- %.3f)'%( da_r0b1b_sf_incM[0], da_r0b1b_sf_incM[1], da_frac_sf_0b_D, da_frac_sf_0b_eD, da_frac_sf_1b_D, da_frac_sf_1b_eD)
    print 'ratio 0b1b: of_incM: %.3f +- %.3f (%.3f +- %.3f / %.3f +- %.3f)'%( da_r0b1b_of_incM[0], da_r0b1b_of_incM[1], da_frac_of_0b_D, da_frac_of_0b_eD, da_frac_of_1b_D, da_frac_of_1b_eD)
    print 'ratio 0b1b: sf_onZ : %.3f +- %.3f (%.3f +- %.3f / %.3f +- %.3f)'%( da_r0b1b_sf_onZ [0], da_r0b1b_sf_onZ [1], da_frac_sf_0b_N, da_frac_sf_0b_eN, da_frac_sf_1b_N, da_frac_sf_1b_eN)
    print 'ratio 0b1b: of_onZ : %.3f +- %.3f (%.3f +- %.3f / %.3f +- %.3f)'%( da_r0b1b_of_onZ [0], da_r0b1b_of_onZ [1], da_frac_of_0b_N, da_frac_of_0b_eN, da_frac_of_1b_N, da_frac_of_1b_eN)

    l2 = r.TLegend(0.2,0.8,0.4,0.89)
    tmp_r = r.TH1F('r0b1b','r_{0b1b}',3,0,3)
    tmp_r.SetBinContent(1, tt_r0b1b_of_incM[0]); tmp_r.SetBinError(1, tt_r0b1b_of_incM[1]); tmp_r.GetXaxis().SetBinLabel(1, 'OF m_{ll} in 61-121');
    tmp_r.SetBinContent(2, tt_r0b1b_of_onZ [0]); tmp_r.SetBinError(2, tt_r0b1b_of_onZ [1]); tmp_r.GetXaxis().SetBinLabel(2, 'OF on-Z');
    tmp_r.SetBinContent(3, tt_r0b1b_sf_incM[0]); tmp_r.SetBinError(3, tt_r0b1b_sf_incM[1]); tmp_r.GetXaxis().SetBinLabel(3, 'SF m_{ll} in 61-121');
    #tmp_r.SetBinContent(2, tt_r0b1b_sf_onZ [0]); tmp_r.SetBinError(2, tt_r0b1b_sf_onZ [1]); tmp_r.GetXaxis().SetBinLabel(2, 'SF on-Z');
    tmp_r.GetYaxis().SetRangeUser(0.00, 0.50)
    tmp_r.GetYaxis().SetTitle('r_{0b1b}')
    tmp_r.GetYaxis().SetTitleOffset(1.2*tmp_r.GetYaxis().GetTitleOffset())
    tmp_r.SetMarkerColor(r.kRed+2)
    tmp_r.SetLineColor  (r.kRed-9)
    tmp_r.SetLineWidth  (2)
    tmp_r.SetMarkerSize (1.2*tmp_r.GetMarkerSize())
    tmp_r.Draw('e1')
    tmp_r.Draw('text00 same')
    l2.AddEntry(tmp_r, 'ttbar MC', 'pl')
    if withData:
        tmp_r_da = tmp_r.Clone('data r0b1b')
        tmp_r_da.SetBinContent(1, da_r0b1b_of_incM[0]); tmp_r_da.SetBinError(1, da_r0b1b_of_incM[1]); tmp_r_da.GetXaxis().SetBinLabel(1, 'OF m_{ll} in 61-121');
        tmp_r_da.SetBinContent(2, da_r0b1b_of_onZ [0]); tmp_r_da.SetBinError(2, da_r0b1b_of_onZ [1]); tmp_r_da.GetXaxis().SetBinLabel(2, 'OF on-Z');
        tmp_r_da.SetBinContent(3, 0.                 ); tmp_r_da.SetBinError(3, 0.                 ); tmp_r_da.GetXaxis().SetBinLabel(3, 'SF m_{ll} in 61-121');
        #tmp_r_da.SetBinContent(3, 0.                 ); tmp_r_da.SetBinError(3, 0.                 ); tmp_r_da.GetXaxis().SetBinLabel(3, 'SF on-Z');
        tmp_r_da.SetMarkerStyle(22)
        tmp_r_da.SetMarkerColor(r.kBlack)
        tmp_r_da.SetLineColor  (r.kGray+1)
        tmp_r_da.SetLineWidth  (2)
        tmp_r_da.Draw('e1 ')
        tmp_r_da.Draw('text00 same')
        tmp_r.Draw('e1 same')
        tmp_r.Draw('text00 same')
        l2.AddEntry(tmp_r_da, 'data', 'pl')
    l2.Draw('same')
    lat.DrawLatex(0.76, 0.93, '%.2f fb^{-1}'%(lumi))
    c1.SaveAs('plots/ttbarMET/%s/r0b1b_highStats%s.png' %(lumistr, '_withData' if withData else ''))
    c1.SaveAs('plots/ttbarMET/%s/r0b1b_highStats%s.pdf' %(lumistr, '_withData' if withData else ''))
    c1.SaveAs('plots/ttbarMET/%s/r0b1b_highStats%s.root'%(lumistr, '_withData' if withData else ''))
    
    
    return tmp_h,tmp_r, tmp_h_da, tmp_r_da
    

def closureTest(mo, withData=False):
    mt2 = 80;
    fmll  = 0.330; fmll_e  = 0.020
    r0b1b = 0.279; r0b1b_e = 0.021
    rsfof = 1.085; rsfof_e = 0.050
    met_sf_0b_onZ     = copy.deepcopy(getattr(mo, 'sf_0b_mt2{mt2}_TT_met'         .format(mt2=mt2))) ## final histogram to closure from
    met_of_1b_incM    = copy.deepcopy(getattr(mo, 'of_1b_mt2{mt2}_61mll121_TT_met'.format(mt2=mt2))) ## initial histo for estimation
    met_of_1b_incM_da = copy.deepcopy(getattr(mo, 'of_1b_mt2{mt2}_61mll121_DA_met'.format(mt2=mt2))) ## initial DATA histogram
    #met_of_1b_incM_da.Scale(10./2.6) ## simple scaling

    h_fmll  = met_sf_0b_onZ.Clone('h_fmll' ); h_fmll  .Reset(); h_fmll  .Sumw2()
    h_r0b1b = met_sf_0b_onZ.Clone('h_r0b1b'); h_r0b1b .Reset(); h_r0b1b .Sumw2()
    h_rsfof = met_sf_0b_onZ.Clone('h_rsfof'); h_rsfof .Reset(); h_rsfof .Sumw2()
    
    for i in range(1,h_fmll.GetNbinsX()+1):
        h_fmll .SetBinContent(i, fmll ); h_fmll .SetBinError(i, fmll_e )
        h_r0b1b.SetBinContent(i, r0b1b); h_r0b1b.SetBinError(i, r0b1b_e)
        h_rsfof.SetBinContent(i, rsfof); h_rsfof.SetBinError(i, rsfof_e)

    ## scale the estimator by fmll
    met_of_1b_onZ_est  = copy.deepcopy(met_of_1b_incM)
    met_of_1b_onZ_est.Multiply(h_fmll)
    met_of_1b_onZ_est_da  = copy.deepcopy(met_of_1b_incM_da)
    met_of_1b_onZ_est_da.Multiply(h_fmll)
        
    ## scale the estimator by r0b1b
    met_of_0b_onZ_est  = copy.deepcopy(met_of_1b_onZ_est)
    met_of_0b_onZ_est.Multiply(h_r0b1b)
    met_of_0b_onZ_est_da  = copy.deepcopy(met_of_1b_onZ_est_da)
    met_of_0b_onZ_est_da.Multiply(h_r0b1b)
    
    ## scale the estimator by rsfof
    met_sf_0b_onZ_est  = copy.deepcopy(met_of_0b_onZ_est)
    met_sf_0b_onZ_est.Multiply(h_rsfof)
    met_sf_0b_onZ_est_da  = copy.deepcopy(met_of_0b_onZ_est_da)
    met_sf_0b_onZ_est_da.Multiply(h_rsfof)

    ## cosmetics
    met_sf_0b_onZ       .SetName('observed  ttbar')
    met_sf_0b_onZ_est   .SetName('estimated ttbar')
    met_sf_0b_onZ_est_da.SetName('estimated data')
    met_sf_0b_onZ       .SetMarkerColor(r.kBlue-1)
    met_sf_0b_onZ_est   .SetMarkerColor(r.kRed)
    met_sf_0b_onZ_est_da.SetMarkerColor(r.kBlack)
    met_sf_0b_onZ_est_da.SetMarkerStyle(22)

    helper.ensureDirectory('plots/ttbarMET/%s/'%(lumistr))
    retval = helper.makeRatioPlot([met_sf_0b_onZ, met_sf_0b_onZ_est,met_sf_0b_onZ_est_da], 'plots/ttbarMET/%s/closure_highStats%s'%(lumistr,'_withData' if withData else ''), {'fitratio': [1,'pol1']})
    print 'this is retval', retval
    print 'this is retval.GetParamter(1)', retval.GetParameter(1)
    print 'this is retval.GetParError(1)', retval.GetParError(1)
    doCorrection = abs(retval.GetParameter(1)) > abs(retval.GetParError(1))
    print 'i am going to do the correction', doCorrection
    met_sf_0b_onZ_est_tt_corr = copy.deepcopy(met_sf_0b_onZ_est.Clone('tt corrected'))
    if doCorrection: met_sf_0b_onZ_est_tt_corr.Divide(retval)
    met_sf_0b_onZ_est_tt_corr.SetMarkerColor(r.kAzure)
    met_sf_0b_onZ_est_tt_corr.SetLineColor  (r.kGray+2)
    met_sf_0b_onZ_est_tt_corr.SetMarkerStyle(24)
    met_sf_0b_onZ_est_da_corr = copy.deepcopy(met_sf_0b_onZ_est_da.Clone('data est. corrected'))
    if doCorrection: met_sf_0b_onZ_est_da_corr.Divide(retval)
    met_sf_0b_onZ_est_da_corr.SetMarkerColor(r.kGreen+2)
    helper.makeRatioPlot([met_sf_0b_onZ, met_sf_0b_onZ_est,met_sf_0b_onZ_est_da, met_sf_0b_onZ_est_da_corr, met_sf_0b_onZ_est_tt_corr], 'plots/ttbarMET/%s/closure_highStats%s'%(lumistr,'_withDataCorrected' if withData else ''))

    bin100 = met_sf_0b_onZ.FindBin(100.)
    bin150 = met_sf_0b_onZ.FindBin(150.)
    nbins  = met_sf_0b_onZ.GetNbinsX()

    obs100p_e = r.Double(); exp100p_e = r.Double();
    obs150p_e = r.Double(); exp150p_e = r.Double();
    obs100p   = met_sf_0b_onZ    .IntegralAndError(bin100, nbins+1, obs100p_e)
    obs150p   = met_sf_0b_onZ    .IntegralAndError(bin150, nbins+1, obs150p_e)
    exp100p   = met_sf_0b_onZ_est.IntegralAndError(bin100, nbins+1, exp100p_e)
    exp150p   = met_sf_0b_onZ_est.IntegralAndError(bin150, nbins+1, exp150p_e)
    ratio100 = ratioError(obs100p, obs100p_e, exp100p, exp100p_e)
    ratio150 = ratioError(obs150p, obs150p_e, exp150p, exp150p_e)

    print 'inclusive closure for MET 100+'
    print 'obs: {obs:.3f} +- {obse:.3f}\t\t{exp:.3f} +- {expe:.3f}'.format(obs=obs100p, obse=obs100p_e, exp=exp100p, expe=exp100p_e)
    print 'ratio: {rat:.3f} +- {rate:.3f}'.format(rat=ratio100[0],rate=ratio100[1])
    print '================================='
    print 'inclusive closure for MET 150+'
    print 'obs: {obs:.3f} +- {obse:.3f}\t\t{exp:.3f} +- {expe:.3f}'.format(obs=obs150p, obse=obs150p_e, exp=exp150p, expe=exp150p_e)
    print 'ratio: {rat:.3f} +- {rate:.3f}'.format(rat=ratio150[0],rate=ratio150[1])

    obs100p_tt_corr_e = r.Double(); exp100p_tt_corr_e = r.Double();
    obs150p_tt_corr_e = r.Double(); exp150p_tt_corr_e = r.Double();
    obs100p_da_corr_e = r.Double(); exp100p_da_corr_e = r.Double();
    obs150p_da_corr_e = r.Double(); exp150p_da_corr_e = r.Double();
    exp100p_tt_corr   = met_sf_0b_onZ_est_tt_corr.IntegralAndError(bin100, nbins+1, exp100p_tt_corr_e)
    exp150p_tt_corr   = met_sf_0b_onZ_est_tt_corr.IntegralAndError(bin150, nbins+1, exp150p_tt_corr_e)
    exp100p_da_corr   = met_sf_0b_onZ_est_da_corr.IntegralAndError(bin100, nbins+1, exp100p_da_corr_e)
    exp150p_da_corr   = met_sf_0b_onZ_est_da_corr.IntegralAndError(bin150, nbins+1, exp150p_da_corr_e)
    ratio100_tt_corr = ratioError(obs100p, obs100p_e, exp100p_tt_corr, exp100p_tt_corr_e)
    ratio150_tt_corr = ratioError(obs150p, obs150p_e, exp150p_tt_corr, exp150p_tt_corr_e)
    ratio100_da_corr = ratioError(obs100p, obs100p_e, exp100p_da_corr, exp100p_da_corr_e)
    ratio150_da_corr = ratioError(obs150p, obs150p_e, exp150p_da_corr, exp150p_da_corr_e)

    print '\n\ncorrected closure for MET 100+ ttbar'
    print 'obs: {obs:.3f} +- {obse:.3f}\t\t{exp:.3f} +- {expe:.3f}'.format(obs=obs100p, obse=obs100p_e, exp=exp100p_tt_corr, expe=exp100p_tt_corr_e)
    print 'ratio: {rat:.3f} +- {rate:.3f}'.format(rat=ratio100_tt_corr[0],rate=ratio100_tt_corr[1])
    print 'corrected closure for MET 100+ data'
    print 'obs: {obs:.3f} +- {obse:.3f}\t\t{exp:.3f} +- {expe:.3f}'.format(obs=obs100p, obse=obs100p_e, exp=exp100p_da_corr, expe=exp100p_da_corr_e)
    print 'ratio: {rat:.3f} +- {rate:.3f}'.format(rat=ratio100_da_corr[0],rate=ratio100_da_corr[1])
    print '================================='
    print 'corrected closure for MET 150+ ttbar'
    print 'obs: {obs:.3f} +- {obse:.3f}\t\t{exp:.3f} +- {expe:.3f}'.format(obs=obs150p, obse=obs150p_e, exp=exp150p_tt_corr, expe=exp150p_tt_corr_e)
    print 'ratio: {rat:.3f} +- {rate:.3f}'.format(rat=ratio150_tt_corr[0],rate=ratio150_tt_corr[1])
    print 'corrected closure for MET 150+ data'
    print 'obs: {obs:.3f} +- {obse:.3f}\t\t{exp:.3f} +- {expe:.3f}'.format(obs=obs150p, obse=obs150p_e, exp=exp150p_da_corr, expe=exp150p_da_corr_e)
    print 'ratio: {rat:.3f} +- {rate:.3f}'.format(rat=ratio150_da_corr[0],rate=ratio150_da_corr[1])

    tabClosureNonCorrString = '''\\documentclass[12pt,a4paper]{{article}}
\\begin{{document}}
\\begin{{table}}[hbtp] 
\\begin{{center}} 
\\bgroup 
\\def\\arraystretch{{1.2}} 
\\caption{{Closure for E$_{{T}}^{{miss.}}$ for the electroweak baseline signal region before slope correction.}} 
\\label{{tab:ewkClosureTTbar}} 
\\begin{{tabular}}{{l| c c }} 
              & E$_{{T}}^{{miss.}}$ $>$ 100 GeV & E$_{{T}}^{{miss.}}$ $>$ 150 GeV \\\\ \\hline 
obs.  (tt)   &  {obs100:.2f}     $\\pm$  {obs100e:.2f}    & {obs150:.2f}     $\\pm$  {obs150e:.2f}       \\\\ \\hline 
pred. (tt)   &  {pred100tt:.2f}  $\\pm$  {pred100tte:.2f} & {pred150tt:.2f}  $\\pm$  {pred150tte:.2f}    \\\\ \\hline 
obs./pred. (tt)   &  {rat100tt:.2f}   $\\pm$  {rat100tte:.2f}  & {rat150tt:.2f}   $\\pm$  {rat150tte:.2f}     \\\\ \\hline 

\\end{{tabular}} 
\\egroup 
\\end{{center}} 
\\end{{table}} 
\\end{{document}}'''.format(obs100=obs100p, obs100e=obs100p_e, obs150=obs150p, obs150e=obs150p_e,
                          pred100tt=exp100p, pred100tte=exp100p_e, pred150tt=exp150p, pred150tte=exp150p_e,
                          rat100tt=ratio100[0], rat100tte=ratio100[1], rat150tt=ratio150[0], rat150tte=ratio150[1])
    helper.ensureDirectory('plots/ttbarMET/%s/tables/'%lumistr)
    nonCorrTableFile = open('plots/ttbarMET/%s/tables/table_closureNonCorrected.tex'%(lumistr),'w')
    nonCorrTableFile.write(tabClosureNonCorrString)
    nonCorrTableFile.close()
    
        
    tabClosureCorrString = '''\\documentclass[12pt,a4paper]{{article}}
\\begin{{document}}
\\begin{{table}}[hbtp] 
\\begin{{center}} 
\\bgroup 
\\def\\arraystretch{{1.2}} 
\\caption{{Closure for E$_{{T}}^{{miss.}}$ for the electroweak baseline signal region with slope correction applied.}} 
\\label{{tab:ewkClosureTTbarCorr}} 
\\begin{{tabular}}{{l| c c }} 
              & E$_{{T}}^{{miss.}}$ $>$ 100 GeV & E$_{{T}}^{{miss.}}$ $>$ 150 GeV \\\\ \\hline 
obs.  (tt)   &  {obs100:.2f}     $\\pm$  {obs100e:.2f}    & {obs150:.2f}     $\\pm$  {obs150e:.2f}       \\\\ \\hline 
pred. (tt)   &  {pred100tt:.2f}  $\\pm$  {pred100tte:.2f} & {pred150tt:.2f}  $\\pm$  {pred150tte:.2f}    \\\\ \\hline 
pred. (data) &  {pred100da:.2f}  $\\pm$  {pred100dae:.2f} & {pred150da:.2f}  $\\pm$  {pred150dae:.2f}    \\\\ \\hline \\hline 
obs./pred. (tt)   &  {rat100tt:.2f}   $\\pm$  {rat100tte:.2f}  & {rat150tt:.2f}   $\\pm$  {rat150tte:.2f}     \\\\ \\hline 
obs./pred. (data) &  {rat100da:.2f}   $\\pm$  {rat100dae:.2f}  & {rat150da:.2f}   $\\pm$  {rat150dae:.2f}     \\\\ \\hline 

\\end{{tabular}} 
\\egroup 
\\end{{center}} 
\\end{{table}}
\\end{{document}}'''.format(obs100=obs100p, obs100e=obs100p_e, obs150=obs150p, obs150e=obs150p_e,
                          pred100tt=exp100p_tt_corr, pred100tte=exp100p_tt_corr_e, pred150tt=exp150p_tt_corr, pred150tte=exp150p_tt_corr_e,
                          pred100da=exp100p_da_corr, pred100dae=exp100p_da_corr_e, pred150da=exp150p_da_corr, pred150dae=exp150p_da_corr_e,
                          rat100tt=ratio100_tt_corr[0], rat100tte=ratio100_tt_corr[1], rat150tt=ratio150_tt_corr[0], rat150tte=ratio150_tt_corr[1],
                          rat100da=ratio100_da_corr[0], rat100dae=ratio100_da_corr[1], rat150da=ratio150_da_corr[0], rat150dae=ratio150_da_corr[1] )
    corrTableFile = open('plots/ttbarMET/%s/tables/table_closureCorrected.tex'%(lumistr),'w')
    corrTableFile.write(tabClosureCorrString)
    corrTableFile.close()

def makeEWKinoSummaryPlot(fileMarc, fileVince):
    return 0
    
        
 
if __name__ == "__main__":
    global opts

    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option('-s', '--samples' , action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-b', '--batch'   , action='store_true', dest='batchMode' , default=False            , help='set batch mode. default false')
    parser.add_option('-v', '--variable', action='store', type=str, dest='variable' , default=''            , help='do certain variable')
    parser.add_option('-l', '--sleptons', action='store_true', dest='isSlepton' , default=False   , help='is this a slepton scan')
    (opts, args) = parser.parse_args()


    if opts.batchMode:
        print 'i am setting this stupid batch mode now'
        gROOT.SetBatch(True)
    gROOT.ProcessLine('.L include/tdrstyle.C')
    r.setTDRStyle() 

    print 'running with these options'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)

    global lumi, lumistr, cutstring, cuts
    #lumi = 10.
    #lumi = 0.8; maxrun = 274240
    #lumi = 4.0; maxrun = 999999; lumistr = 'lumi4invfbTTPOW'
    #lumi = 7.65; maxrun = 999999; lumistr = 'lumi7.65invfbTTPOW'
    lumi = 12.9; maxrun = 999999; lumistr = 'lumi12.9invfbTTPOW'
    


    cuts = CutManager.CutManager() ## True means ewino. adapts the pts and the dR. and requires == 2 tight leptons
    cutlist = {}
    ## cutlist['sf_0b_nomt2'] = cuts.AddList([cuts.goodLepton, cuts.SF, 'nBJetMedium25_Edge == 0', cuts.nj2, 'j1MetDPhi_Edge > 1.', cuts.Zmass, 'run_Edge <= %d'%maxrun])
    ## cutlist['of_0b_nomt2'] = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge == 0', cuts.nj2, 'j1MetDPhi_Edge > 1.', cuts.Zmass, 'run_Edge <= %d'%maxrun])
    ## cutlist['sf_1b_nomt2'] = cuts.AddList([cuts.goodLepton, cuts.SF, 'nBJetMedium25_Edge >= 1', cuts.nj2, 'j1MetDPhi_Edge > 1.', cuts.Zmass, 'run_Edge <= %d'%maxrun])
    ## cutlist['of_1b_nomt2'] = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge >= 1', cuts.nj2, 'j1MetDPhi_Edge > 1.', cuts.Zmass, 'run_Edge <= %d'%maxrun])
                                                                                                                                                                 
    cutlist['sf_0b_mt280'] = cuts.AddList([cuts.goodLepton, cuts.SF, 'nBJetMedium25_Edge == 0', cuts.nj2, 'j1MetDPhi_Edge > 1.', cuts.Zmass, 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    cutlist['of_0b_mt280'] = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge == 0', cuts.nj2, 'j1MetDPhi_Edge > 1.', cuts.Zmass, 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    cutlist['sf_1b_mt280'] = cuts.AddList([cuts.goodLepton, cuts.SF, 'nBJetMedium25_Edge >= 1', cuts.nj2, 'j1MetDPhi_Edge > 1.', cuts.Zmass, 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    cutlist['of_1b_mt280'] = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge >= 1', cuts.nj2, 'j1MetDPhi_Edge > 1.', cuts.Zmass, 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])

    ## z-mass extended   
    ## cutlist['sf_0b_nomt2_61mll121'] = cuts.AddList([cuts.goodLepton, cuts.SF, 'nBJetMedium25_Edge == 0', 'j1MetDPhi_Edge > 1.', cuts.nj2, 'run_Edge <= %d'%maxrun, 'lepsMll_Edge > 61 && lepsMll_Edge < 121'])
    ## cutlist['of_0b_nomt2_61mll121'] = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge == 0', 'j1MetDPhi_Edge > 1.', cuts.nj2, 'run_Edge <= %d'%maxrun, 'lepsMll_Edge > 61 && lepsMll_Edge < 121'])
    ## cutlist['sf_1b_nomt2_61mll121'] = cuts.AddList([cuts.goodLepton, cuts.SF, 'nBJetMedium25_Edge >= 1', 'j1MetDPhi_Edge > 1.', cuts.nj2, 'run_Edge <= %d'%maxrun, 'lepsMll_Edge > 61 && lepsMll_Edge < 121'])
    ## cutlist['of_1b_nomt2_61mll121'] = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge >= 1', 'j1MetDPhi_Edge > 1.', cuts.nj2, 'run_Edge <= %d'%maxrun, 'lepsMll_Edge > 61 && lepsMll_Edge < 121'])
                                                                                                                                                       
    cutlist['sf_0b_mt280_61mll121'] = cuts.AddList([cuts.goodLepton, cuts.SF, 'nBJetMedium25_Edge == 0', 'j1MetDPhi_Edge > 1.', cuts.nj2, 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121'])
    cutlist['of_0b_mt280_61mll121'] = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge == 0', 'j1MetDPhi_Edge > 1.', cuts.nj2, 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121'])
    cutlist['sf_1b_mt280_61mll121'] = cuts.AddList([cuts.goodLepton, cuts.SF, 'nBJetMedium25_Edge >= 1', 'j1MetDPhi_Edge > 1.', cuts.nj2, 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121'])
    cutlist['of_1b_mt280_61mll121'] = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge >= 1', 'j1MetDPhi_Edge > 1.', cuts.nj2, 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121'])
       
    ## dataset selection
    # BACKGROUNDS
    # ===============================
    #ttDatasets = ['TTJets_DiLepton_total']
    ttDatasets = ['TT_pow_ext34']
    wzDatasets = ['WZTo3LNu']
    dyDatasets = ['DYJetsToLL_M50_HT100to200_ext', 'DYJetsToLL_M50_HT200to400_ext', 'DYJetsToLL_M50_HT400to600_ext', 'DYJetsToLL_M50_HT600toInf_ext']
    zzDatasets = ['ZZ']
    wwDatasets = ['WWTo2L2Nu']
    treeTT  = Sample.Tree(helper.selectSamples(opts.sampleFile, ttDatasets, 'TT' ), 'TT' , 0, isScan = 0)
    treeWZ  = Sample.Tree(helper.selectSamples(opts.sampleFile, wzDatasets, 'WZ' ), 'WZ' , 0, isScan = 0)
    treeZZ  = Sample.Tree(helper.selectSamples(opts.sampleFile, zzDatasets, 'ZZ' ), 'ZZ' , 0, isScan = 0)
    treeWW  = Sample.Tree(helper.selectSamples(opts.sampleFile, wwDatasets, 'WW' ), 'WW' , 0, isScan = 0)
    treeDY  = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets, 'DY' ), 'DY' , 0, isScan = 0)
    #bkgTrees = [treeWZ, treeZZ, treeWW, treeDY, treeTT]
    bkgTrees = [treeTT]

    # DATA
    # ===============================
    #daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2_runs_271036_275125' , 'DoubleEG_Run2016B-PromptReco-v2_runs_271036_275125' , 'MuonEG_Run2016B-PromptReco-v2_runs_271036_275125']
    #daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2_runs_271036_276097', 'DoubleEG_Run2016B-PromptReco-v2_runs_271036_276097', 'MuonEG_Run2016B-PromptReco-v2_runs_271036_276097',
    #              'DoubleMuon_Run2016C-PromptReco-v2_runs_275420_276098', 'DoubleEG_Run2016C-PromptReco-v2_runs_275420_276098', 'MuonEG_Run2016C-PromptReco-v2_runs_275420_276098']
    daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2_runs_271036_276097', 'DoubleEG_Run2016B-PromptReco-v2_runs_271036_276097', 'MuonEG_Run2016B-PromptReco-v2_runs_271036_276097',
                  'DoubleMuon_Run2016C-PromptReco-v2_runs_271036_276811', 'DoubleEG_Run2016C-PromptReco-v2_runs_271036_276811', 'MuonEG_Run2016C-PromptReco-v2_runs_271036_276811',
                  'DoubleMuon_Run2016D-PromptReco-v2_runs_271036_276811', 'DoubleEG_Run2016D-PromptReco-v2_runs_271036_276811', 'MuonEG_Run2016D-PromptReco-v2_runs_271036_276811']
    treeDA  = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA' ), 'DATA' , 1)

    ## marc 80x # SIGNAL
    ## marc 80x # ===============================
    ## marc 80x sig350_20  = ['TChiNeu_WZ_mCha350_mLSP20' ]
    ## marc 80x sig300_20  = ['TNeuCha300_20']
    ## marc 80x sig400_20  = ['TNeuCha400_20']
    ## marc 80x sig350_100 = ['TChiNeu_WZ_mCha350_mLSP100']
    ## marc 80x sig200_100 = ['TChiNeu_WZ_mCha200_mLSP100']
    ## marc 80x tree35020  = Sample.Tree(helper.selectSamples(opts.sampleFile, sig350_20 , 'SIG35020' ), 'SIG35020' , 0, isScan = 1)
    ## marc 80x tree30020  = Sample.Tree(helper.selectSamples(opts.sampleFile, sig300_20 , 'SIG30020' ), 'SIG30020' , 0, isScan = 1)
    ## marc 80x tree40020  = Sample.Tree(helper.selectSamples(opts.sampleFile, sig400_20 , 'SIG40020' ), 'SIG40020' , 0, isScan = 1)
    ## marc 80x tree350100 = Sample.Tree(helper.selectSamples(opts.sampleFile, sig350_100, 'SIG350100'), 'SIG350100', 0, isScan = 1)
    ## marc 80x tree200100 = Sample.Tree(helper.selectSamples(opts.sampleFile, sig200_100, 'SIG200100'), 'SIG200100', 0, isScan = 1)

    #controlPlots(bkgTrees, treeDA, cutlist)


    mo = metObject(cutlist, bkgTrees, treeDA, onlyTT = True)


## what we need is a MET spectrum for ttbar in high MT2 and SF, 0b, on the Z mass
## we have all other combinations available. namely, the met spectrum at MT2 == 0

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
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership, TLatex
import math, sys, optparse, array, copy
import gc, inspect

import include.LeptonSF
import include.nll
import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample



class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def dumpObjects():
    gc.collect()
    oo = gc.get_objects()
    for o in oo:
        if getattr(o, "__class__", None):
            name = o.__class__.__name__
            if name not in exclude:
                filename = inspect.getabsfile(o.__class__)            
                print "Object of class:", name, "...",
                print "defined in file:", filename                


def makePlot(lumi, lumi_str, treeMC, var, name, nbin, xmin, xmax, theCut, cuts,labelx, logx):
    normalized = False
    theCutMC = cuts.AddList([theCut, cuts.goodLepton])
    #theCutSIG = cuts.AddList([theCut, cuts.goodLeptonSignal])
    MC   = treeMC.getTH1F(lumi, "hMC_%s"%(name), var, nbin, xmin, xmax, theCutMC, '', labelx)
    MCS  = treeMC.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, theCutMC, "", labelx)
   # SIG1 = sigs[0].getTH1F(lumi, "hSIG1_%s"%(name), var, nbin, xmin, xmax, theCutSIG, '', labelx)
   # SIG2 = sigs[1].getTH1F(lumi, "hSIG2_%s"%(name), var, nbin, xmin, xmax, theCutSIG, '', labelx)
   # SIG3 = sigs[2].getTH1F(lumi, "hSIG3_%s"%(name), var, nbin, xmin, xmax, theCutSIG, '', labelx)
   # SIG4 = sigs[3].getTH1F(lumi, "hSIG4_%s"%(name), var, nbin, xmin, xmax, theCutSIG, '', labelx)
   # SIG5 = sigs[4].getTH1F(lumi, "hSIG5_%s"%(name), var, nbin, xmin, xmax, theCutSIG, '', labelx)

    for _i,_h in enumerate(MCS.GetHists()):
        print "sample ", _h.GetName(), " has integral ", _h.Integral()
        
    if normalized:
        hists = []
        for _i,_h in enumerate(MCS.GetHists()):
            yie_e = r.Double()
            yie = _h.IntegralAndError(1, _h.GetNbinsX()+1, yie_e)
            _h.SetName(_h.GetName().split('_')[-2]+' %.1f +- %.1f'%(yie, yie_e))
            _h.Scale(1./_h.Integral())
            _h.SetLineColor(_h.GetFillColor())
            _h.SetFillColor(0)
            _h.SetLineWidth(2)
            hists.append(_h)
        yie_e = r.Double()
        yie = MC.IntegralAndError(1, MC.GetNbinsX()+1, yie_e)
        MC.SetName('Full MC')
        MC  .Scale(1./MC.Integral())
        MC  .SetLineColor(r.kRed)
        MC  .SetLineWidth(2)         
        hists.append(MC)
    #SIG1.SetLineStyle(2);SIG2.SetLineStyle(2);SIG3.SetLineStyle(2);SIG4.SetLineStyle(2);SIG5.SetLineStyle(2);
    #SIG1.SetLineWidth(2);SIG2.SetLineWidth(2);SIG3.SetLineWidth(2);SIG4.SetLineWidth(2);SIG5.SetLineWidth(2);
    #print "SIG1 ", SIG1.Integral()
    maxVal = MC.GetMaximum() #GetBinContent(MC.GetMaximumBin())
    if not logx:
        MCS.SetMaximum(1.5*maxVal)
    else:
        MCS.SetMaximum(10000*maxVal)

    print name
    SetOwnership(MC, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(MCS, 0 )   # 0 = release (not keep), 1 = keep
    plot = Canvas.Canvas('DataMC/%s/plot_%s'%(lumi_str,name), 'png,pdf,root', 0.6, 0.5, 0.85, 0.9)
    if not normalized:
        plot.addStack(MCS, "HIST", 1, 1)
    else:
        for h in hists:
            plot.addHisto(h, 'hist,same' if not h.GetName() == 'data' else 'p,same', h.GetName(), 'PL', h.GetLineColor(), 1, 0)
    plot.save(1, 1, logx, lumi)                                                                                                                                           


def makeTable(lumi, lumi_str, treeMC, sigs,var, name, nbin, xmin, xmax, theCut, cuts,labelx, logx, names):
    normalized = False
    err1 = r.Double()    
    err2 = r.Double()    
    err3 = r.Double()    
    theCut1 = cuts.AddList([theCut, cuts.goodLepton, 'met_Edge > 100 && met_Edge < 175'])
    theCut2 = cuts.AddList([theCut, cuts.goodLepton, 'met_Edge > 175 && met_Edge < 250'])
    theCut3 = cuts.AddList([theCut, cuts.goodLepton, 'met_Edge > 250'])
    theCutSIG1 = cuts.AddList([theCut, cuts.goodLeptonSignal, 'met_Edge > 100 && met_Edge < 175'])
    theCutSIG2 = cuts.AddList([theCut, cuts.goodLeptonSignal, 'met_Edge > 175 && met_Edge < 250'])
    theCutSIG3 = cuts.AddList([theCut, cuts.goodLeptonSignal, 'met_Edge > 250'])
    metbins = [100, 175, 250]
    MC1   = treeMC.getTH1F(lumi, "hMC1_%s"%(name), var, nbin, xmin, xmax, theCut1, '', labelx)
    MC2   = treeMC.getTH1F(lumi, "hMC2_%s"%(name), var, nbin, xmin, xmax, theCut2, '', labelx)
    MC3   = treeMC.getTH1F(lumi, "hMC3_%s"%(name), var, nbin, xmin, xmax, theCut3, '', labelx)
    MCS1  = treeMC.getStack(lumi, "hMCS1_%s"%(name), var, nbin, xmin, xmax, theCut1, "", labelx)
    MCS2  = treeMC.getStack(lumi, "hMCS2_%s"%(name), var, nbin, xmin, xmax, theCut2, "", labelx)
    MCS3  = treeMC.getStack(lumi, "hMCS3_%s"%(name), var, nbin, xmin, xmax, theCut3, "", labelx)

    SIGlowMet1 = sigs[0].getTH1F(lumi, "hSIG11_%s"%(name), var, nbin, xmin, xmax, theCutSIG1, '', labelx)
    SIGlowMet2 = sigs[1].getTH1F(lumi, "hSIG12_%s"%(name), var, nbin, xmin, xmax, theCutSIG1, '', labelx)
    SIGlowMet3 = sigs[2].getTH1F(lumi, "hSIG13_%s"%(name), var, nbin, xmin, xmax, theCutSIG1, '', labelx)
    SIGlowMet4 = sigs[3].getTH1F(lumi, "hSIG14_%s"%(name), var, nbin, xmin, xmax, theCutSIG1, '', labelx)
    SIGmidMet1 = sigs[0].getTH1F(lumi, "hSIG21_%s"%(name), var, nbin, xmin, xmax, theCutSIG2, '', labelx)
    SIGmidMet2 = sigs[1].getTH1F(lumi, "hSIG22_%s"%(name), var, nbin, xmin, xmax, theCutSIG2, '', labelx)
    SIGmidMet3 = sigs[2].getTH1F(lumi, "hSIG23_%s"%(name), var, nbin, xmin, xmax, theCutSIG2, '', labelx)
    SIGmidMet4 = sigs[3].getTH1F(lumi, "hSIG24_%s"%(name), var, nbin, xmin, xmax, theCutSIG2, '', labelx)
    SIGhiMet1 = sigs[0].getTH1F(lumi, "hSIG31_%s"%(name), var, nbin, xmin, xmax, theCutSIG3, '', labelx)
    SIGhiMet2 = sigs[1].getTH1F(lumi, "hSIG32_%s"%(name), var, nbin, xmin, xmax, theCutSIG3, '', labelx)
    SIGhiMet3 = sigs[2].getTH1F(lumi, "hSIG33_%s"%(name), var, nbin, xmin, xmax, theCutSIG3, '', labelx)
    SIGhiMet4 = sigs[3].getTH1F(lumi, "hSIG34_%s"%(name), var, nbin, xmin, xmax, theCutSIG3, '', labelx)

    line0 ='       & 100 $<$ MET $<$ 175 & 175 $<$ MET $<$ 250 & MET $>$ 250 \\\\'
    line1 =' ttbar ' 
    line2 =' WW2l  ' 
    line3 =' DY    ' 
    line4 =' ZZ2l  ' 
    line5 =' WZ3l  ' 
    line6 =' rares ' 
    #line7 =' VVV   ' 
    line8 =' TTV   ' 
    line9 =' Total MC & %.1f $\pm$ %.1f& %.1f $\pm$ %.1f& %.1f $\pm$ %.1f' %(MC1.IntegralAndError(1, MC1.GetNbinsX()+1, err1), err1, MC2.IntegralAndError(1, MC2.GetNbinsX()+1, err2), err2,MC3.IntegralAndError(1, MC3.GetNbinsX()+1, err3), err3) 
    line10 ='m_{l}: %s, m_{LSP} : %s & %.1f $\pm$ %.1f& %.1f $\pm$ %.1f& %.1f $\pm$ %.1f \\\\' %(names[0][0], names[0][1], SIGlowMet1.IntegralAndError(1, SIGlowMet1.GetNbinsX()+1, err1), err1, SIGmidMet1.IntegralAndError(1, SIGmidMet1.GetNbinsX()+1, err2), err2,SIGhiMet1.IntegralAndError(1, SIGhiMet1.GetNbinsX()+1, err3), err3) 
    line11 ='m_{l}: %s, m_{LSP} : %s & %.1f $\pm$ %.1f& %.1f $\pm$ %.1f& %.1f $\pm$ %.1f \\\\' %(names[1][0], names[1][1], SIGlowMet2.IntegralAndError(1, SIGlowMet2.GetNbinsX()+1, err1), err1, SIGmidMet2.IntegralAndError(1, SIGmidMet2.GetNbinsX()+1, err2), err2,SIGhiMet2.IntegralAndError(1, SIGhiMet2.GetNbinsX()+1, err3), err3) 
    line12 ='m_{l}: %s, m_{LSP} : %s & %.1f $\pm$ %.1f& %.1f $\pm$ %.1f& %.1f $\pm$ %.1f \\\\' %(names[2][0], names[2][1], SIGlowMet3.IntegralAndError(1, SIGlowMet3.GetNbinsX()+1, err1), err1, SIGmidMet3.IntegralAndError(1, SIGmidMet3.GetNbinsX()+1, err2), err2,SIGhiMet3.IntegralAndError(1, SIGhiMet3.GetNbinsX()+1, err3), err3) 
    line13 ='m_{l}: %s, m_{LSP} : %s & %.1f $\pm$ %.1f& %.1f $\pm$ %.1f& %.1f $\pm$ %.1f \\\\' %(names[3][0], names[3][1], SIGlowMet4.IntegralAndError(1, SIGlowMet4.GetNbinsX()+1, err1), err1, SIGmidMet4.IntegralAndError(1, SIGmidMet4.GetNbinsX()+1, err2), err2,SIGhiMet4.IntegralAndError(1, SIGhiMet4.GetNbinsX()+1, err3), err3) 
    #line10 ='m_{l}: %s, m_{LSP} : %s & %.1f $\pm$ %.1f& %.1f $\pm$ %.1f& %.1f $\pm$ %.1f \\\\' %(names[0][0], names[0][1], SIG11.IntegralAndError(1, SIG11.GetNbinsX()+1, err1), err1, SIG12.IntegralAndError(1, SIG12.GetNbinsX()+1, err2), err2,SIG13.IntegralAndError(1, SIG13.GetNbinsX()+1, err3), err3) 
    #line11 ='m_{l}: %s, m_{LSP} : %s & %.1f $\pm$ %.1f& %.1f $\pm$ %.1f& %.1f $\pm$ %.1f \\\\' %(names[1][0], names[1][1], SIG21.IntegralAndError(1, SIG21.GetNbinsX()+1, err1), err1, SIG22.IntegralAndError(1, SIG22.GetNbinsX()+1, err2), err2,SIG23.IntegralAndError(1, SIG23.GetNbinsX()+1, err3), err3) 
    #line12 ='m_{l}: %s, m_{LSP} : %s & %.1f $\pm$ %.1f& %.1f $\pm$ %.1f& %.1f $\pm$ %.1f \\\\' %(names[2][0], names[2][1], SIG31.IntegralAndError(1, SIG31.GetNbinsX()+1, err1), err1, SIG32.IntegralAndError(1, SIG32.GetNbinsX()+1, err2), err2,SIG33.IntegralAndError(1, SIG33.GetNbinsX()+1, err3), err3) 

    err = r.Double()    
    for _i1, _h1 in enumerate(MCS1.GetHists()):
        print "sample ", _h1.GetName(), " has integral ", _h1.Integral()
        if  _h1.GetName() == 'auxStack_block_hMCS1_mt2_SR_TTJets_blockHisto': line1 += ' & %.1f $\pm$ %.1f'  %(_h1.IntegralAndError(1, _h1.GetNbinsX()+1, err), err)
        if  _h1.GetName() == 'auxStack_block_hMCS1_mt2_SR_WW2l_blockHisto'  : line2 += ' & %.1f $\pm$ %.1f'  %(_h1.IntegralAndError(1, _h1.GetNbinsX()+1, err), err)
        if  _h1.GetName() == 'auxStack_block_hMCS1_mt2_SR_DYJets_blockHisto': line3 += ' & %.1f $\pm$ %.1f'  %(_h1.IntegralAndError(1, _h1.GetNbinsX()+1, err), err)
        if  _h1.GetName() == 'auxStack_block_hMCS1_mt2_SR_ZZ2l_blockHisto'  : line4 += ' & %.1f $\pm$ %.1f'  %(_h1.IntegralAndError(1, _h1.GetNbinsX()+1, err), err)
        if  _h1.GetName() == 'auxStack_block_hMCS1_mt2_SR_WZ3l_blockHisto'  : line5 += ' & %.1f $\pm$ %.1f'  %(_h1.IntegralAndError(1, _h1.GetNbinsX()+1, err), err)
        if  _h1.GetName() == 'auxStack_block_hMCS1_mt2_SR_rares_blockHisto' : line6 += ' & %.1f $\pm$ %.1f'  %(_h1.IntegralAndError(1, _h1.GetNbinsX()+1, err), err)
        #if  _h1.GetName() == 'auxStack_block_hMCS1_mt2_SR_VVV_blockHisto'   : line7 += ' & %.1f $\pm$ %.1f'  %(_h1.IntegralAndError(1, _h1.GetNbinsX()+1, err), err)
        if  _h1.GetName() == 'auxStack_block_hMCS1_mt2_SR_TTV_blockHisto'   : line8 += ' & %.1f $\pm$ %.1f'  %(_h1.IntegralAndError(1, _h1.GetNbinsX()+1, err), err)    
    
    for _i2, _h2 in enumerate(MCS2.GetHists()):
        if  _h2.GetName() == 'auxStack_block_hMCS2_mt2_SR_TTJets_blockHisto': line1 += ' & %.1f $\pm$ %.1f'  %(_h2.IntegralAndError(1, _h2.GetNbinsX()+1, err), err)
        if  _h2.GetName() == 'auxStack_block_hMCS2_mt2_SR_WW2l_blockHisto'  : line2 += ' & %.1f $\pm$ %.1f'  %(_h2.IntegralAndError(1, _h2.GetNbinsX()+1, err), err)
        if  _h2.GetName() == 'auxStack_block_hMCS2_mt2_SR_DYJets_blockHisto': line3 += ' & %.1f $\pm$ %.1f'  %(_h2.IntegralAndError(1, _h2.GetNbinsX()+1, err), err)
        if  _h2.GetName() == 'auxStack_block_hMCS2_mt2_SR_ZZ2l_blockHisto'  : line4 += ' & %.1f $\pm$ %.1f'  %(_h2.IntegralAndError(1, _h2.GetNbinsX()+1, err), err)
        if  _h2.GetName() == 'auxStack_block_hMCS2_mt2_SR_WZ3l_blockHisto'  : line5 += ' & %.1f $\pm$ %.1f'  %(_h2.IntegralAndError(1, _h2.GetNbinsX()+1, err), err)
        if  _h2.GetName() == 'auxStack_block_hMCS2_mt2_SR_rares_blockHisto' : line6 += ' & %.1f $\pm$ %.1f'  %(_h2.IntegralAndError(1, _h2.GetNbinsX()+1, err), err)
        #if  _h2.GetName() == 'auxStack_block_hMCS2_mt2_SR_VVV_blockHisto'   : line7 += ' & %.1f $\pm$ %.1f'  %(_h2.IntegralAndError(1, _h2.GetNbinsX()+1, err), err)
        if  _h2.GetName() == 'auxStack_block_hMCS2_mt2_SR_TTV_blockHisto'   : line8 += ' & %.1f $\pm$ %.1f'  %(_h2.IntegralAndError(1, _h2.GetNbinsX()+1, err), err)    
        
    for _i3, _h3 in enumerate(MCS3.GetHists()):
        if  _h3.GetName() == 'auxStack_block_hMCS3_mt2_SR_TTJets_blockHisto': line1 += ' & %.1f $\pm$ %.1f \\\\'  %(_h3.IntegralAndError(1, _h3.GetNbinsX()+1, err), err)
        if  _h3.GetName() == 'auxStack_block_hMCS3_mt2_SR_WW2l_blockHisto'  : line2 += ' & %.1f $\pm$ %.1f \\\\'  %(_h3.IntegralAndError(1, _h3.GetNbinsX()+1, err), err)
        if  _h3.GetName() == 'auxStack_block_hMCS3_mt2_SR_DYJets_blockHisto': line3 += ' & %.1f $\pm$ %.1f \\\\'  %(_h3.IntegralAndError(1, _h3.GetNbinsX()+1, err), err)
        if  _h3.GetName() == 'auxStack_block_hMCS3_mt2_SR_ZZ2l_blockHisto'  : line4 += ' & %.1f $\pm$ %.1f \\\\'  %(_h3.IntegralAndError(1, _h3.GetNbinsX()+1, err), err)
        if  _h3.GetName() == 'auxStack_block_hMCS3_mt2_SR_WZ3l_blockHisto'  : line5 += ' & %.1f $\pm$ %.1f \\\\'  %(_h3.IntegralAndError(1, _h3.GetNbinsX()+1, err), err)
        if  _h3.GetName() == 'auxStack_block_hMCS3_mt2_SR_rares_blockHisto' : line6 += ' & %.1f $\pm$ %.1f \\\\'  %(_h3.IntegralAndError(1, _h3.GetNbinsX()+1, err), err)
        #if  _h3.GetName() == 'auxStack_block_hMCS3_mt2_SR_VVV_blockHisto'   : line7 += ' & %.1f $\pm$ %.1f \\\\'  %(_h3.IntegralAndError(1, _h3.GetNbinsX()+1, err), err)
        if  _h3.GetName() == 'auxStack_block_hMCS3_mt2_SR_TTV_blockHisto'   : line8 += ' & %.1f $\pm$ %.1f \\\\'  %(_h3.IntegralAndError(1, _h3.GetNbinsX()+1, err), err) 
        
    line0 += '\\hline'; line8 += '\\hline'; line9 += '\\\\\hline' ;  line13 += '\\hline\hline'                                                       
    

    print line0
    print line1
    print line2
    print line3
    print line4
    print line5
    print line6
    #print line7
    print line8
    print line9
    print line10
    print line11
    print line12
    print line13


if __name__ == "__main__":

    print bcolors.HEADER
    print '#######################################################################'
    print '                  Starting slepton analysis...                          '
    print '#######################################################################' + bcolors.ENDC
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()
     
    print bcolors.HEADER + '[Data - MC comparisons] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    dyDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    zzDatasets = ['ZZTo2L2Nu',  'GluGluToContinToZZTo2mu2nu', 'GluGluToContinToZZTo2e2nu']
    wzDatasets = ['WZTo3LNu']
    othersDatasets = ['WWZ', 'WZZ', 'ZZZ', 'TWZ', 'tZq_ll', 'TTZToLLNuNu_ext2', 'TTZToQQ', 'TTLLJets_m1to10', 'TTHnobb_pow', 'VHToNonbb']
    fsDatasets = ['TTTT', 'TTTo2L2Nu', 'TBar_tch_powheg', 'T_tch_powheg', 'WWTo2L2Nu', 'WWW', 'WWDouble', 'WpWpJJ', 'TTWToLNu_ext2',  'TTWToQQ', 'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT'] 


    mcDatasets = fsDatasets+dyDatasets + othersDatasets + zzDatasets + wzDatasets
    #mcDatasets = ['GGHWWTo2L2Nu', 'GGWWTo2L2Nu', 'WZTo3LNu','WWZ',  'WZZ', 'ZZZ', 'DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50', 'WWTo2L2Nu', 'WWDouble', 'WpWpJJ', 'TTTo2L2Nu', 'ZZTo2L2Nu', 'tZq_ll', 'TTZToLLNuNu_ext2', 'TWZ',  'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT','TTZToQQ', 'TTLLJets_m1to10', 'GluGluToContinToZZTo2e2nu', 'GluGluToContinToZZTo2mu2nu']
    #mcDatasets = ['GGHWWTo2L2Nu', 'GGWWTo2L2Nu', 'GGHZZ4L', 'ZZTo4L', 'WZTo3LNu','WWG',  'WWW', 'WWZ',  'WZZ', 'ZZZ', 'ZGTo2LG', 'DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50', 'WWTo2L2Nu', 'WWDouble', 'WpWpJJ', 'WGToLNuG', 'TTTo2L2Nu', 'ZZTo2L2Nu', 'TTHnobb_pow', 'VHToNonbb', 'TWZ', 'WZTo2L2Q',  'TBar_tch_powheg', 'T_tch_powheg', 'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT', 'TTTT',  'TTZToQQ', 'TTWToQQ', 'TTLLJets_m1to10', 'TTWToLNu_ext2', 'GluGluToContinToZZTo2e2tau', 'GluGluToContinToZZTo2mu2tau', 'GluGluToContinToZZTo2mu2nu', 'GluGluToContinToZZTo4e', 'GluGluToContinToZZTo4mu', 'GluGluToContinToZZTo4tau']

    signals = ['SMS_450_50', 'SMS_350_150', 'SMS_250_180', 'SMS_100_1']
    signames = [[450, 50], [350, 150], [250, 180], [100, 1]]
    #xsec = [ 2*270.79, 2*270.79,  2*270.79]
    xsec = [2*(0.77+ 0.3),2*(2.33+0.89), 2*(9.210+3.470), 2*270.79]


    treeSIG1 = Sample.Tree(helper.selectSamples("samplesSlepton.dat", [signals[0]], 'SI'), 'SI', 0, isScan = xsec[0])
    treeSIG2 = Sample.Tree(helper.selectSamples("samplesSlepton.dat", [signals[1]], 'SI'), 'SI', 0, isScan = xsec[1])
    treeSIG3 = Sample.Tree(helper.selectSamples("samplesSlepton.dat", [signals[2]], 'SI'), 'SI', 0, isScan = xsec[2])
    treeSIG4 = Sample.Tree(helper.selectSamples("samplesSlepton.dat", [signals[3]], 'SI'), 'SI', 0, isScan = xsec[3])
    sigs = [treeSIG1, treeSIG2, treeSIG3, treeSIG4]                                                               

    treeMC = Sample.Tree(helper.selectSamples("samplesSlepton.dat", mcDatasets, 'MC'), 'MC'  , 0, isScan = 0)
    print bcolors.HEADER + '[Data - MC comparisons] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    lumi = 35.9 ; maxrun = 999999; lumi_str = '35.9invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()
    labelmll = 'm_{ll} [GeV]'
    labelmet = 'E_{T}^{miss} [GeV]'
    labelnjet = "N. Jets"
    cut  = "nJetSel_Edge >= 0"
    #cut  = "lepsMll_Edge < 101 && lepsMll_Edge > 81"
    name = "SR"
    metbins = [100, 175, 250]
    makeTable(lumi, lumi_str, treeMC, sigs, "mt2_Edge","mt2_"+name, 40, 0, 400, cuts.AddList([cuts.SF, cuts.bvetoLoose25, cuts.ThirdLeptonVeto, cuts.ZvetoExt, ' htJet25j_Edge  == 0 && mt2_Edge > 90']),cuts, "mt2 [GeV] (SF, njet = 0)",  0, signames)

    makePlot(lumi, lumi_str, treeMC, "mt2_Edge","mt2_incMET", 40, 0, 400, cuts.AddList([cuts.SF, cuts.bvetoLoose25, cuts.ThirdLeptonVeto, cuts.ZvetoExt, ' htJet25j_Edge  == 0 && met_Edge > 100 ']),cuts, "MT_{2} (E_{T}^{miss} > 100) [GeV]",  0)
    makePlot(lumi, lumi_str, treeMC, "mt2_Edge","mt2_lowMET", 40, 0, 400, cuts.AddList([cuts.SF, cuts.bvetoLoose25, cuts.ThirdLeptonVeto, cuts.ZvetoExt, ' htJet25j_Edge  == 0 && met_Edge > 100 && met_Edge < 175']),cuts, "MT_{2} (100 < E_{T}^{miss} < 175) [GeV]",  0)
    makePlot(lumi, lumi_str, treeMC, "mt2_Edge","mt2_medMET", 40, 0, 400, cuts.AddList([cuts.SF, cuts.bvetoLoose25, cuts.ThirdLeptonVeto, cuts.ZvetoExt, ' htJet25j_Edge  == 0 && met_Edge > 175 && met_Edge < 250']),cuts, "MT_{2} (175 < E_{T}^{miss} < 250) [GeV]",  0)
    makePlot(lumi, lumi_str, treeMC, "mt2_Edge","mt2_hiMET", 40, 0, 400, cuts.AddList([cuts.SF, cuts.bvetoLoose25, cuts.ThirdLeptonVeto, cuts.ZvetoExt, ' htJet25j_Edge  == 0 && met_Edge > 250']),cuts, "MT_{2} (E_{T}^{miss} > 250) [GeV]",  0)
    makePlot(lumi, lumi_str, treeMC, "met_Edge","met_SR", 40, 0, 400, cuts.AddList([cuts.SF, cuts.bvetoLoose25, cuts.ThirdLeptonVeto, cuts.ZvetoExt, ' htJet25j_Edge  == 0 && mt2_Edge > 90']),cuts, "E_{T}^{miss} [GeV] ",  0)                                                  



















### makePlot(lumi, lumi_str, treeMC, treeDA, "lepsMll_Edge","mll_AF_debug"+name, 40, 0, 400,cuts.AddList([ cuts.bvetoLoose25, cuts.ThirdLeptonVeto, ' htJet25j_Edge  == 0 &&  mt2_Edge > 90 && met_Edge > 100 && met_Edge < 175']),cuts, "m_{ll} [GeV] (AF, njet = 0)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "nBJetLoose25_Edge","nBJetLoose25_debug"+name, 4, 0, 4, cuts.AddList([cuts.OF,cuts.ThirdLeptonVeto, cuts.ZvetoExt, ' htJet25j_Edge  == 0 &&  mt2_Edge > 90 && met_Edge > 100 && met_Edge < 175']),cuts, "nbjets (OF, njet = 0)",  0)
   ### makePlot(lumi, lumi_str, treeMC, treeDA, "nJet25_Edge","nJet25_debug"+name, 4, 0, 4, cuts.AddList([cuts.OF, cuts.bvetoLoose25, cuts.ThirdLeptonVeto, cuts.ZvetoExt, '   mt2_Edge > 90 && met_Edge > 100 && met_Edge < 175']),cuts, "nbjets (OF, njet = 0)",  0)
   ### makePlot(lumi, lumi_str, treeMC, treeDA, "lepsMll_Edge","mll_debug"+name, 40, 0, 400,cuts.AddList([cuts.OF, cuts.bvetoLoose25, cuts.ThirdLeptonVeto, ' htJet25j_Edge  == 0 &&  mt2_Edge > 90 && met_Edge > 100 && met_Edge < 175']),cuts, "m_{ll} [GeV] (OF, njet = 0)",  0)
   ### makePlot(lumi, lumi_str, treeMC, treeDA, "mt2_Edge","mt2_debug"+name, 40, 0, 400, cuts.AddList([cuts.OF, cuts.bvetoLoose25, cuts.ThirdLeptonVeto, cuts.ZvetoExt, ' htJet25j_Edge  == 0 && met_Edge > 100 && met_Edge < 175']),cuts, "mt2 [GeV] (OF, njet = 0)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "JetSel_Edge_pt","JetPt_SF"+name, 30, 0, 300, cuts.AddList([cuts.SF,  cut]),cuts, "pt jets > 35 GeV  [GeV] (SF)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "JetSel_Edge_eta","JetEta_SF"+name, 30, -5, 5, cuts.AddList([cuts.SF,  cut]),cuts, "eta jets > 35 GeV  [GeV] (SF)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "nJet35_Edge","nJet35_SF"+name, 8, 0, 8, cuts.AddList([cuts.SF,  cut]),cuts, "number of jets > 35 GeV  [GeV] (SF)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "nJet35_Edge","nJet35_OF"+name, 8, 0, 8, cuts.AddList([cuts.OF,  cut]),cuts, "number of jets > 35 GeV  [GeV] (OF)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "lepsMll_Edge","mll_DYCR"+name, 40, 0, 400, cuts.AddList([cuts.SF,  cut, "met_Edge < 50"]),cuts, "m_{ll} [GeV] (DY CR)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "lepsMll_Edge","mll_ttbarCR"+name, 40, 0, 400, cuts.AddList([cuts.OF,  cut, "met_Edge > 100 && met_Edge < 150 && nJetSel_Edge == 2"]),cuts, "m_{ll} [GeV] (ttbar CR)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "lepsZPt_Edge","zPt_SF_njet0"+name, 40, 0, 400, cuts.AddList([cuts.SF, cut, cuts.slep0jet ]),cuts, "Z p_{T} GeV (SF, njet = 0)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "lepsZPt_Edge","zPt_SF_njet0"+name, 40, 0, 400, cuts.AddList([cuts.SF, cut ]),cuts, "Z p_{T} GeV (OF, njet = 0)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "lepsZPt_Edge","zPt_OF_njet0"+name, 40, 0, 400, cuts.AddList([cuts.OF, cut, cuts.slep0jet ]),cuts, "Z p_{T} GeV (OF, njet = 0)",  0)
   #### makePlot(lumi, lumi_str, treeMC, treeDA, "nJet35_Edge","nJet35_SF"+name, 4, 0, 4, cuts.AddList([cuts.SF, cut ]),cuts, "number of jets > 35 GeV (SF)",  0)
   #### makePlot(lumi, lumi_str, treeMC, treeDA, "nJet35_Edge","nJet35_OF"+name, 4, 0, 4, cuts.AddList([cuts.OF, cut]),cuts, "number of jets > 35 GeV (OF)",  0)
   #### makePlot(lumi, lumi_str, treeMC, treeDA, "nJet35_Edge","nJet35_SF"+name, 4, 0, 4, cuts.AddList([cuts.SF, cut ]),cuts, "number of jets > 35 GeV (SF)",  0)
   #### makePlot(lumi, lumi_str, treeMC, treeDA, "nJet35_Edge","nJet35_ee"+name, 4, 0, 4, cuts.AddList([cuts.ee, cut ]),cuts, "number of jets > 35 GeV (SF)",  0)
   #### makePlot(lumi, lumi_str, treeMC, treeDA, "nJet35_Edge","nJet35_mm"+name, 4, 0, 4, cuts.AddList([cuts.mm, cut ]),cuts, "number of jets > 35 GeV (SF)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "lepsMll_Edge","mll_OF_njet1"+name, 40, 0, 400, cuts.AddList([cuts.OF, cuts.nj25eq1, cut]),cuts, "m_{ll} [GeV] (OF, njet = 1)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "lepsMll_Edge","mll_SF_njet1"+name, 40, 0, 400, cuts.AddList([cuts.SF, cuts.nj25eq1, cut]),cuts, "m_{ll} [GeV] (SF, njet = 1)",  0)
   ### makePlot(lumi, lumi_str, treeMC, treeDA, "lepsMll_Edge","mll_OF_njet0"+name, 40, 0, 400, cuts.AddList([cuts.OF, cuts.slep0jet, cut]),cuts, "m_{ll} [GeV] (OF, njet = 0)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "lepsMll_Edge","mll_SF_njet0"+name, 40, 0, 400, cuts.AddList([cuts.SF, cuts.slep0jet, cut]),cuts, "m_{ll} [GeV] (SF, njet = 0)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "met_Edge","met_SR_OF_njet0"+name, 30, 100, 400, cuts.AddList([cuts.OF, cuts.slep0jet]),cuts, "E_{T}^{miss} [GeV] (OF, njet = 0)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "met_Edge","met_SR_OF_njet1"+name, 30, 100, 400, cuts.AddList([cuts.OF, cuts.nj25eq1]),cuts, "E_{T}^{miss} [GeV] (OF, njet = 1)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "mt2_Edge","mt2_OF_njet1"+name, 40, 0, 400, cuts.AddList([cuts.OF, cuts.nj25eq1]),cuts, "mt2 [GeV] (OF, njet = 1)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "mt2_Edge","mt2_SF_njet1"+name, 40, 0, 400, cuts.AddList([cuts.SF, cuts.nj25eq1]),cuts, "mt2 [GeV] (SF, njet = 1)",  0)
   ### makePlot(lumi, lumi_str, treeMC, treeDA, "mt2_Edge","mt2_OF_njet0"+name, 40, 0, 400, cuts.AddList([cuts.OF, cuts.slep0jet]),cuts, "mt2 [GeV] (OF, njet = 0)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "mt2_Edge","mt2_SF_njet0"+name, 40, 0, 400, cuts.AddList([cuts.SF, cuts.slep0jet]),cuts, "mt2 [GeV] (SF, njet = 0)",  0)
   ### makePlot(lumi, lumi_str, treeMC, treeDA, "met_Edge","met_OF_njet0"+name, 40, 0, 200, cuts.AddList([cuts.OF, cuts.slep0jet]),cuts, "E_{T}^{miss} [GeV] (OF, njet = 0)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "met_Edge","met_SF_njet0"+name, 40, 0, 200, cuts.AddList([cuts.SF, cuts.slep0jet]),cuts, "E_{T}^{miss} [GeV] (SF, njet = 0)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "Lep1_pt_Edge","Lep1_pt_SF_njet0"+name, 20, 0, 400, cuts.AddList([cuts.SF, cuts.slep0jet]),cuts, "Leading lepton p_{T} [GeV] (SF, njet = 0)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "Lep2_pt_Edge","Lep2_pt_SF_njet0"+name, 20, 0, 400, cuts.AddList([cuts.SF, cuts.slep0jet]),cuts, "Subleading lepton p_{T} [GeV] (SF, njet = 0)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "Lep1_eta_Edge","Lep1_eta_SF_njet0"+name, 20, -4, 4, cuts.AddList([cuts.SF, cuts.slep0jet]),cuts, "Leading lepton #eta (SF, njet = 0)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "Lep2_eta_Edge","Lep2_eta_SF_njet0"+name, 20, -4, 4, cuts.AddList([cuts.SF, cuts.slep0jet]),cuts, "Subleading lepton #eta (SF, njet = 0)",  0)          
   ### makePlot(lumi, lumi_str, treeMC, treeDA, "Lep1_pt_Edge","Lep1_pt_OF_njet0"+name, 20, 0, 400, cuts.AddList([cuts.OF, cuts.slep0jet]),cuts, "Leading lepton p_{T} [GeV] (OF, njet = 0)",  0)
   ### makePlot(lumi, lumi_str, treeMC, treeDA, "Lep2_pt_Edge","Lep2_pt_OF_njet0"+name, 20, 0, 400, cuts.AddList([cuts.OF, cuts.slep0jet]),cuts, "Subleading lepton p_{T} [GeV] (OF, njet = 0)",  0)
   ### makePlot(lumi, lumi_str, treeMC, treeDA, "Lep1_eta_Edge","Lep1_eta_OF_njet0"+name, 20, -4, 4, cuts.AddList([cuts.OF, cuts.slep0jet]),cuts, "Leading lepton #eta (OF, njet = 0)",  0)
   ### makePlot(lumi, lumi_str, treeMC, treeDA, "Lep2_eta_Edge","Lep2_eta_OF_njet0"+name, 20, -4, 4, cuts.AddList([cuts.OF, cuts.slep0jet]),cuts, "Subleading lepton #eta (OF, njet = 0)",  0)          
   ###
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "met_Edge","met_OF_njet1", 40, 0, 200, cuts.AddList([cuts.OF, cuts.nj25eq1]),cuts, "E_{T}^{miss} [GeV] (OF, njet = 1)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "met_Edge","met_SF_njet1", 40, 0, 200, cuts.AddList([cuts.SF, cuts.nj25eq1]),cuts, "E_{T}^{miss} [GeV] (SF, njet = 1)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "Lep1_pt_Edge","Lep1_pt_SF_njet1", 20, 0, 400, cuts.AddList([cuts.SF, cuts.nj25eq1]),cuts, "Leading lepton p_{T} [GeV] (SF, njet = 1)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "Lep2_pt_Edge","Lep2_pt_SF_njet1", 20, 0, 400, cuts.AddList([cuts.SF, cuts.nj25eq1]),cuts, "Subleading lepton p_{T} [GeV] (SF, njet = 1)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "Lep1_eta_Edge","Lep1_eta_SF_njet1", 20, -4, 4, cuts.AddList([cuts.SF, cuts.nj25eq1]),cuts, "Leading lepton #eta (SF, njet = 1)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "Lep2_eta_Edge","Lep2_eta_SF_njet1", 20, -4, 4, cuts.AddList([cuts.SF, cuts.nj25eq1]),cuts, "Subleading lepton #eta (SF, njet = 1)",  0)          
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "Lep1_pt_Edge","Lep1_pt_OF_njet1", 20, 0, 400, cuts.AddList([cuts.OF, cuts.nj25eq1]),cuts, "Leading lepton p_{T} [GeV] (OF, njet = 1)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "Lep2_pt_Edge","Lep2_pt_OF_njet1", 20, 0, 400, cuts.AddList([cuts.OF, cuts.nj25eq1]),cuts, "Subleading lepton p_{T} [GeV] (OF, njet = 1)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "Lep1_eta_Edge","Lep1_eta_OF_njet1", 20, -4, 4, cuts.AddList([cuts.OF, cuts.nj25eq1]),cuts, "Leading lepton #eta (OF, njet = 1)",  0)
   ### #makePlot(lumi, lumi_str, treeMC, treeDA, "Lep2_eta_Edge","Lep2_eta_OF_njet1", 20, -4, 4, cuts.AddList([cuts.OF, cuts.nj25eq1]),cuts, "Subleading lepton #eta (OF, njet = 1)",  0)          





   ### 
   ### makePlot(lumi, lumi_str, treeMC, treeDA, "met_Edge","met", 20, 0, 400, cuts.AddList([cuts.slep0jetNoMET, cuts.OF]),cuts, "E_{T}^{miss} [GeV]",  0)
   ### makePlot(lumi, lumi_str, treeMC, treeDA, "mt2_Edge","mt2", 20, 0, 400, cuts.AddList([cuts.slep0jetNoMT2, cuts.OF]),cuts, "mt2 [GeV]",  0)
   ### #makePlot(lumi, lumi_str, treeMC, sigs, "lepsMll_Edge","mll_HT", 40, 0, 400, cuts.AddList([cuts.slepHT,dPhiCut]),cuts, "mll HT > 0 ",  0, signames)
   ### #makePlot(lumi, lumi_str, treeMC, sigs, "metl1DPhi_Edge","metl1DPhi", 20, 0, 4, cuts.AddList([cuts.slepTOT,dPhiCut,"met_Edge > 100"]),cuts, "#Delta#Phi(E_{T}^{miss}, l_{1})  E_{T}^{miss} > 100 ",  0, signames)

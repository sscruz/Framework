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

    

def makeSummaryTable3l4l(plot_3l, plot_4l):
    wz3l, wz3l_e, zz4l, zz4l_e = 0., 0., 0., 0.,
    bkg3l, bkg3l_e, bkg4l, bkg4l_e = 0., 0., 0., 0.,
    obs3l, obs4l = 0, 0
    rat3l, rat3l_e, rat4l, rat4l_e = 0., 0., 0., 0.,
    
    for i, h_3l in enumerate(plot_3l.histos):
        if not i: continue
        tmp_double = r.Double()
        if 'WZ' in h_3l.GetName():
            wz3l = h_3l.IntegralAndError(1, h_3l.GetNbinsX()+1, tmp_double)
            wz3l_e = tmp_double
        elif not 'DATA' in h_3l.GetName(): 
            tmp_yield = h_3l.IntegralAndError(1, h_3l.GetNbinsX()+1, tmp_double)
            bkg3l  += tmp_yield
            bkg3l_e = math.sqrt(bkg3l_e**2 + tmp_double**2)
        elif 'DATA' in h_3l.GetName():
            obs3l = h_3l.Integral()
        
    for i, h_4l in enumerate(plot_4l.histos):
        if not i: continue
        tmp_double = r.Double()
        if 'ZZ' in h_4l.GetName():
            zz4l = h_4l.IntegralAndError(1, h_4l.GetNbinsX()+1, tmp_double)
            zz4l_e = tmp_double
        elif not 'DATA' in h_4l.GetName(): 
            tmp_yield = h_4l.IntegralAndError(1, h_4l.GetNbinsX()+1, tmp_double)
            bkg4l  += tmp_yield
            bkg4l_e = math.sqrt(bkg4l_e**2 + tmp_double**2)
        elif 'DATA' in h_4l.GetName():
            obs4l = h_4l.Integral()
        
    obsSub3l, obsSub3l_e, obsSub4l, obsSub4l_e = obs3l-bkg3l, math.sqrt(bkg3l_e**2 + obs3l), obs4l-bkg4l, math.sqrt(bkg4l_e**2 + obs4l)
    rat3l, rat3l_e = helper.ratioError(wz3l, wz3l_e, obsSub3l, obsSub3l_e)
    rat4l, rat4l_e = helper.ratioError(zz4l, zz4l_e, obsSub4l, obsSub4l_e)

    table3l4lComparisonString = '''\\documentclass[12pt,a4paper]{{article}}
\\begin{{document}}
\\begin{{table}}[hbtp] 
\\begin{{center}} 
\\bgroup 
\\def\\arraystretch{{1.2}} 
\\caption{{3-lepton and 4-lepton control regions. Signal MC is WZ-$>$3lnu in the 3-lepton region and ZZ-$>$4l in the 4-lepton region.}} 
\\label{{tab:tab3l4lcontrol}} 
\\begin{{tabular}}{{l| c c }} 
              & 3-lepton region & 4-lepton region \\\\ \\hline \\hline
signal MC        & {wz3l:.2f}     $\\pm$  {wz3l_e:.2f}       & {zz4l:.2f}     $\\pm$  {zz4l_e:.2f}       \\\\ \\hline
bkg. MC          & {bkg3l:.2f}  $\\pm$  {bkg3l_e:.2f}        & {bkg4l:.2f}  $\\pm$  {bkg4l_e:.2f}    \\\\ \\hline \\hline
\\textbf{{data}}             & \\textbf{{{obs3l}}}                               & \\textbf{{{obs4l}}}      \\\\ \\hline
data-bkg.        & {obsSub3l:.2f}   $\\pm$  {obsSub3l_e:.2f} & {obsSub4l:.2f}   $\\pm$  {obsSub4l_e:.2f}     \\\\ \\hline \\hline
sig./(data-bkg.) & {rat3l:.2f}   $\\pm$  {rat3l_e:.2f}       & {rat4l:.2f}   $\\pm$  {rat4l_e:.2f}     \\\\ \\hline

\\end{{tabular}} 
\\egroup 
\\end{{center}} 
\\end{{table}} 
\\end{{document}}'''.format(wz3l=wz3l, wz3l_e=wz3l_e, zz4l=zz4l, zz4l_e=zz4l_e,
                            bkg3l=bkg3l, bkg3l_e=bkg3l_e, bkg4l=bkg4l,bkg4l_e=bkg4l_e, 
                            obs3l=int(obs3l),obs4l=int(obs4l),
                            obsSub3l=obsSub3l, obsSub3l_e=obsSub3l_e, obsSub4l=obsSub4l, obsSub4l_e=obsSub4l_e,
                            rat3l=rat3l, rat3l_e=rat3l_e, rat4l=rat4l, rat4l_e=rat4l_e)
    compTableFile = open('plots/DataMC/4.0invfb/tables/controlRegion3l4l_comparison.tex','w')
    compTableFile.write(table3l4lComparisonString)
    compTableFile.close()

def makePlot(lumi, lumi_str, treeDA, treeMC, var, name, nbin, xmin, xmax, theCut, cuts,labelx, logx, returnplot = False):

    theCutDATA = cuts.AddList([theCut, cuts.trigger])
    MC = treeMC.getTH1F(lumi, "hMC_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx)
    MCS = treeMC.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, theCut, "", labelx)
    DATA = treeDA.getTH1F(lumi, "hDATA_%s"%(name), var, nbin, xmin, xmax, theCutDATA, '', labelx)

    maxValmc = MC.GetMaximum() #GetBinContent(MC.GetMaximumBin())
    maxValdata = DATA.GetMaximum() #GetBinContent(DATA.GetMaximumBin())
    maxVal = max(maxValmc, maxValdata)
    if not logx:
        MC.GetYaxis().SetRangeUser(0, 1.3*maxVal)
        MCS.SetMaximum(1.3*maxVal)
        DATA.GetYaxis().SetRangeUser(0, 1.3*maxVal)
    else:
        MC.GetYaxis().SetRangeUser(0.1, 2.0*maxVal)
        MCS.SetMaximum(2.0*maxVal)
        DATA.GetYaxis().SetRangeUser(0.1, 2.0*maxVal)
   
    print name
    SetOwnership(MC, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(MCS, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(DATA, 0 )   # 0 = release (not keep), 1 = keep
    plot = Canvas.Canvas('DataMC/%s/plot_%s'%(lumi_str,name), 'png,pdf,root', 0.65, 0.7, 0.85, 0.9)
    plot.addStack(MCS, "HIST", 1, 1)
    plot.addHisto(DATA, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plot.saveRatio(1, 1, logx, lumi, DATA, MC)
    
    if returnplot:
        return plot
    else:
        del plot
        return

def rebinEWino(histo):
    a= [0,50,100,150,225,300,325]
    _arr = array.array('d', a)
    h_ret = r.TH1F(histo.GetName(),histo.GetTitle(), 6, _arr)
    for ii,ib in enumerate(a):
        if ib == a[-1]: continue
        tmp_err = r.Double()
        if not ib == 300:
            tmp_cont = histo.IntegralAndError(histo.FindBin(ib), histo.FindBin(a[ii+1]-1.), tmp_err )
        else:
            tmp_cont = histo.IntegralAndError(histo.FindBin(ib), histo.GetNbinsX()+1, tmp_err )
        h_ret.SetBinContent(ii+1, tmp_cont)
        h_ret.SetBinError  (ii+1, tmp_err )
    return copy.deepcopy(h_ret)

def makeSummaryEWino(plot_EWino):
    f_ttbar = r.TFile('plots/ttbarMET/lumi800invpb/closure_highStats_withDataCorrected.root','READ')
    f_vince = r.TFile('plots/ttbarMET/lumi800invpb/data_SR_EWK_hists.root','READ')
    c_ttbar = f_ttbar.Get('c1_n2')
    h_ttPredMC = c_ttbar.FindObject("estimated ttbar")
    h_ttPredDA = c_ttbar.FindObject("data est. corrected")
    h_vince = copy.deepcopy(f_vince.Get('h_templ_met'))
    h_vince.Rebin(25)
    
    h_vincePred = h_ttPredMC.Clone('vincePredictionFinal')
    h_vincePred.Reset()
    for ib in range(1,h_vincePred.GetNbinsX()+2):
        h_vincePred.SetBinContent(ib, h_vince.GetBinContent(ib))
        h_vincePred.SetBinError  (ib, h_vince.GetBinError  (ib))
        if ib == h_vincePred.FindBin(301.):
            tmp_err = r.Double()
            tmp_cont = h_vince.IntegralAndError(ib, h_vince.GetNbinsX()+1, tmp_err)
            h_vincePred.SetBinContent(ib, tmp_cont)
            h_vincePred.SetBinError  (ib, tmp_err )
            break

    for i,h_ew in enumerate(plot_EWino.histos):
        if 'WZ'   in h_ew.GetName():
            h_wz = copy.deepcopy(h_ew)
        if 'ZZ'   in h_ew.GetName():
            h_zz = copy.deepcopy(h_ew)
        if 'rare' in h_ew.GetName():
            h_rare = copy.deepcopy(h_ew)
        if 'DATA' in h_ew.GetName():
            h_data = copy.deepcopy(h_ew)
    
    h_vincePred = rebinEWino(h_vincePred); h_vincePred.SetFillColor(r.kAzure+5)
    h_wz        = rebinEWino(h_wz       ); h_wz       .SetFillColor(r.kRed+2); h_wz.GetXaxis().SetTitle('MET (GeV)')
    h_zz        = rebinEWino(h_zz       ); h_wz       .SetFillColor(r.kGreen+1)
    h_rare      = rebinEWino(h_rare     ); h_rare     .SetFillColor(r.kGray)
    h_ttPredDA  = rebinEWino(h_ttPredDA ); h_ttPredDA .SetFillColor(r.kYellow+1)
    h_data      = rebinEWino(h_data     );
    bin75 = h_data.FindBin(75.)
    bkg = h_wz.GetBinContent(bin75)+h_zz.GetBinContent(bin75)+h_rare.GetBinContent(bin75)+h_ttPredDA.GetBinContent(bin75)
    data = h_data.GetBinContent(bin75)
    dy = h_vincePred.GetBinContent(bin75)

    h_vincePred.Scale((data-bkg)/dy)

    stack = r.THStack()
    stack.Add(h_wz       )
    stack.Add(h_zz       )
    stack.Add(h_rare     )
    stack.Add(h_ttPredDA )
    stack.Add(h_vincePred)

    totalBKG = h_vincePred.Clone('totalBkg')
    totalBKG.Add(h_wz)
    totalBKG.Add(h_zz)
    totalBKG.Add(h_rare)
    totalBKG.Add(h_ttPredDA)
    totalBKG.SetFillColor(1)
    totalBKG.SetMarkerSize(0.)
    totalBKG.SetFillStyle(3004)

    res = {}
    for i in range(1,totalBKG.GetNbinsX()+1):
        res['ra%d'%i] = h_rare     .GetBinContent(i)
        res['ra%d_e'%i] = h_rare     .GetBinError  (i)
        res['wz%d'%i] = h_wz       .GetBinContent(i)
        res['wz%d_e'%i] = h_wz       .GetBinError  (i)
        res['zz%d'%i] = h_zz       .GetBinContent(i)
        res['zz%d_e'%i] = h_zz       .GetBinError  (i)
        res['dy%d'%i] = h_vincePred.GetBinContent(i)
        res['dy%d_e'%i] = h_vincePred.GetBinError  (i)
        res['tt%d'%i] = h_ttPredDA .GetBinContent(i)
        res['tt%d_e'%i] = h_ttPredDA .GetBinError  (i)
        res['to%d'%i] = totalBKG .GetBinContent(i)
        res['to%d_e'%i] = totalBKG .GetBinError  (i)
        res['ob%d'%i] = h_data     .GetBinContent(i)
        res['ob%d_e'%i] = h_data     .GetBinError  (i)
        res['di%d'%i], res['di%d_e'%i] = helper.ratioError(h_data.GetBinContent(i), h_data.GetBinError(i), totalBKG.GetBinContent(i), totalBKG.GetBinError(i) )

    stack.SetMaximum(1.3*max(totalBKG.GetMaximum(),h_data.GetMaximum()) )
    plot = Canvas.Canvas('ttbarMET/%s/plot_resultFancyEWino'%('0.8invfb'), 'png,pdf,root', 0.7, 0.7, 0.95, 0.9)
    plot.addStack(stack, "HIST", 1, 1)
    plot.addHisto(totalBKG, "E2,SAME", "", "", r.kBlack, 1, -1)
    plot.addHisto(h_data, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plot.saveRatio(1, 1, 0, 0.8, h_data, totalBKG)

    
    tableResultsEWino = '''\\documentclass{{article}}
\\usepackage[a4paper,margin=1in,landscape]{{geometry}}
\\begin{{document}}
\\begin{{table}}[hbtp] 
\\begin{{center}} 
\\bgroup 
\\def\\arraystretch{{1.2}} 
\\caption{{Result for the EWKino search for 0.8 fb$^{{-1}}$.}} 
\\label{{tab:ewinoResult}} 
\\begin{{tabular}}{{l| c c c c c c}} 
 MET region      & 0-50 GeV                      & 50-100 GeV                   & 100-150 GeV                  & 150-225 GeV                  & 225-300 GeV                 & 300+ GeV \\\\ \\hline
rare             & {ra1:.2f} $\\pm$  {ra1_e:.2f} & {ra2:.2f} $\\pm$ {ra2_e:.2f} & {ra3:.2f} $\\pm$ {ra3_e:.2f} & {ra4:.2f} $\\pm$ {ra4_e:.2f} & {ra5:.2f} $\\pm$ {ra5_e:.2f}& {ra6:.2f} $\\pm$ {ra6_e:.2f}           \\\\
WZ3l             & {wz1:.2f} $\\pm$  {wz1_e:.2f} & {wz2:.2f} $\\pm$ {wz2_e:.2f} & {wz3:.2f} $\\pm$ {wz3_e:.2f} & {wz4:.2f} $\\pm$ {wz4_e:.2f} & {wz5:.2f} $\\pm$ {wz5_e:.2f}& {wz6:.2f} $\\pm$ {wz6_e:.2f}           \\\\
ZZ               & {zz1:.2f} $\\pm$  {zz1_e:.2f} & {zz2:.2f} $\\pm$ {zz2_e:.2f} & {zz3:.2f} $\\pm$ {zz3_e:.2f} & {zz4:.2f} $\\pm$ {zz4_e:.2f} & {zz5:.2f} $\\pm$ {zz5_e:.2f}& {zz6:.2f} $\\pm$ {zz6_e:.2f}           \\\\
DY prediction    & {dy1:.2f} $\\pm$  {dy1_e:.2f} & {dy2:.2f} $\\pm$ {dy2_e:.2f} & {dy3:.2f} $\\pm$ {dy3_e:.2f} & {dy4:.2f} $\\pm$ {dy4_e:.2f} & {dy5:.2f} $\\pm$ {dy5_e:.2f}& {dy6:.2f} $\\pm$ {dy6_e:.2f}           \\\\
tt prediction    & {tt1:.2f} $\\pm$  {tt1_e:.2f} & {tt2:.2f} $\\pm$ {tt2_e:.2f} & {tt3:.2f} $\\pm$ {tt3_e:.2f} & {tt4:.2f} $\\pm$ {tt4_e:.2f} & {tt5:.2f} $\\pm$ {tt5_e:.2f}& {tt6:.2f} $\\pm$ {tt6_e:.2f}           \\\\ \\hline
total bkg        & {to1:.2f} $\\pm$  {to1_e:.2f} & {to2:.2f} $\\pm$ {to2_e:.2f} & {to3:.2f} $\\pm$ {to3_e:.2f} & {to4:.2f} $\\pm$ {to4_e:.2f} & {to5:.2f} $\\pm$ {to5_e:.2f}& {to6:.2f} $\\pm$ {to6_e:.2f}           \\\\ \\hline
observed         & {ob1:.2f}                     & {ob2:.2f}                    & {ob3:.2f}                    & {ob4:.2f}                    & {ob5:.2f}                   & {ob6:.2f}                          \\\\
obs./pred.       & {di1:.2f} $\\pm$  {di1_e:.2f} & {di2:.2f} $\\pm$  {di2_e:.2f}& {di3:.2f} $\\pm$  {di3_e:.2f}& {di4:.2f} $\\pm$  {di4_e:.2f}& {di5:.2f} $\\pm$ {di5_e:.2f}& {di6:.2f} $\\pm$ {di6_e:.2f}      \\\\
\\end{{tabular}} 
\\egroup 
\\end{{center}} 
\\end{{table}} 
\\end{{document}}'''.format(
ra1=res['ra1'],  ra1_e=res['ra1_e'], ra2=res['ra2'], ra2_e=res['ra2_e'], ra3=res['ra3'], ra3_e=res['ra3_e'], ra4=res['ra4'], ra4_e=res['ra4_e'], ra5=res['ra5'], ra5_e=res['ra5_e'], ra6=res['ra6'], ra6_e=res['ra6_e'],
wz1=res['wz1'],  wz1_e=res['wz1_e'], wz2=res['wz2'], wz2_e=res['wz2_e'], wz3=res['wz3'], wz3_e=res['wz3_e'], wz4=res['wz4'], wz4_e=res['wz4_e'], wz5=res['wz5'], wz5_e=res['wz5_e'], wz6=res['wz6'], wz6_e=res['wz6_e'],
zz1=res['zz1'],  zz1_e=res['zz1_e'], zz2=res['zz2'], zz2_e=res['zz2_e'], zz3=res['zz3'], zz3_e=res['zz3_e'], zz4=res['zz4'], zz4_e=res['zz4_e'], zz5=res['zz5'], zz5_e=res['zz5_e'], zz6=res['zz6'], zz6_e=res['zz6_e'],
dy1=res['dy1'],  dy1_e=res['dy1_e'], dy2=res['dy2'], dy2_e=res['dy2_e'], dy3=res['dy3'], dy3_e=res['dy3_e'], dy4=res['dy4'], dy4_e=res['dy4_e'], dy5=res['dy5'], dy5_e=res['dy5_e'], dy6=res['dy6'], dy6_e=res['dy6_e'],
tt1=res['tt1'],  tt1_e=res['tt1_e'], tt2=res['tt2'], tt2_e=res['tt2_e'], tt3=res['tt3'], tt3_e=res['tt3_e'], tt4=res['tt4'], tt4_e=res['tt4_e'], tt5=res['tt5'], tt5_e=res['tt5_e'], tt6=res['tt6'], tt6_e=res['tt6_e'],
to1=res['to1'],  to1_e=res['to1_e'], to2=res['to2'], to2_e=res['to2_e'], to3=res['to3'], to3_e=res['to3_e'], to4=res['to4'], to4_e=res['to4_e'], to5=res['to5'], to5_e=res['to5_e'], to6=res['to6'], to6_e=res['to6_e'],
ob1=res['ob1'],                      ob2=res['ob2'],                     ob3=res['ob3'],                     ob4=res['ob4'],                     ob5=res['ob5'],                     ob6=res['ob6'],                   
di1=res['di1'],  di1_e=res['di1_e'], di2=res['di2'], di2_e=res['di2_e'], di3=res['di3'], di3_e=res['di3_e'], di4=res['di4'], di4_e=res['di4_e'], di5=res['di5'], di5_e=res['di5_e'], di6=res['di6'], di6_e=res['di6_e'])
    compTableFile = open('plots/ttbarMET/0.8invfb/tables/resultTable0p8invfb.tex','w')
    compTableFile.write(tableResultsEWino)
    compTableFile.close()

    return copy.deepcopy(h_vincePred)
    
    

if __name__ == "__main__":

    print bcolors.HEADER
    print '#######################################################################'
    print '                  Starting r_SFOF analysis...                          '
    print '#######################################################################' + bcolors.ENDC

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()

    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    #dyDatasets = ['DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    dyDatasets = ['DYJetsToLL_M50_HT100to200_ext', 'DYJetsToLL_M50_HT200to400_ext', 'DYJetsToLL_M50_HT400to600_ext', 'DYJetsToLL_M50_HT600toInf_ext']
    mcDatasets = ['ZZ', 'TTJets_DiLepton_total', 'WWTo2L2Nu', 'WZTo2L2Q', 'WZTo3LNu', 'TTZToLLNuNu', 'TTWToLNu']
    mcDatasets += dyDatasets
    daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2_runs_271036_275125', 'DoubleEG_Run2016B-PromptReco-v2_runs_271036_275125', 'MuonEG_Run2016B-PromptReco-v2_runs_271036_275125']



    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    #lumi = 2.6 ; maxrun = 999999; lumi_str = '2.6invfb'
    lumi = 4.0 ; maxrun = 275125; lumi_str = '4.0invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()
    labelmll = 'm_{ll} [GeV]'
    labelmet = 'MET [GeV]'
    labelnjet = "N. Jets"
    
    ## ewino_SR = makePlot(0.8, '0.8invfb', treeDA, treeMC, "met_Edge", "met_ewino_SR", 12, 0, 300, cuts.AddList([cuts.ewinoSR, 'run_Edge <= 274240']), cuts, labelmet, 0, True)
    ## makeSummaryEWino(ewino_SR)
    
    ## #3l and 4l plots
    plot_3l = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge"    , "met_3lregion_AF", 14, 60, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region3l]), cuts, labelmet, 0, True)
    plot_4l = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_4lregion_AF", 10,  0, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region4l]), cuts, labelmll, 0, True)
    makeSummaryTable3l4l(plot_3l, plot_4l)

    ## makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_met50_2jets_OF", 20,  0,  200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj2, 'met_Edge > 50']), cuts, 'M_{T2}', 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_met50_2jets_SF", 20,  0,  200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj2, 'met_Edge > 50']), cuts, 'M_{T2}', 1)

    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_mt2_140to160_2jets_SF", 10,  0,  300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj2, 'mt2_Edge > 140 && mt2_Edge < 160']), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_mt2_140to160_2jets_ee", 10,  0,  300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj2, 'mt2_Edge > 140 && mt2_Edge < 160']), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_mt2_140to160_2jets_mm", 10,  0,  300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj2, 'mt2_Edge > 140 && mt2_Edge < 160']), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nBJetMedium25_Edge", "nb_mt2_140to160_2jets_SF",  3,  -0.5,    2.5, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj2, 'mt2_Edge > 140 && mt2_Edge < 160']), cuts, labelnjet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_mt2_140to160_2jets_SF",  28, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj2, 'mt2_Edge > 140 && mt2_Edge < 160']), cuts, labelmll, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "nj_mt2_140to160_2jets_SF",  5,  1,    6, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj2, 'mt2_Edge > 140 && mt2_Edge < 160']), cuts, labelnjet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_met50_2jets_SF", 20,  0,  200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj2, 'met_Edge > 50']), cuts, 'M_{T2}', 1)

    ## plots for fmll and r0b1b
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge"      , "fmll_region", 20,    0,  200, cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge >= 1', cuts.nj2, 'j1MetDPhi_Edge > 1.', cuts.ZmassExtended, 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.']), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nBJetMedium25_Edge", "r0b1b_region", 2, -0.5,  1.5, cuts.AddList([cuts.goodLepton, cuts.OF, 'j1MetDPhi_Edge > 1.', cuts.nj2, 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.', cuts.ZmassExtended]), cuts, 'N_{b-jets}', 0)

    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nVert_Edge", "nvtx_inclusive_OF", 40,  0,  40, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj1]), cuts, 'n_{vertices}', 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nVert_Edge", "nvtx_inclusive_SF", 40,  0,  40, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj1]), cuts, 'n_{vertices}', 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nVert_Edge", "nvtx_inclusive_EE", 40,  0,  40, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj1]), cuts, 'n_{vertices}', 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nVert_Edge", "nvtx_inclusive_MM", 40,  0,  40, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj1]), cuts, 'n_{vertices}', 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_inclusive_SF", 40, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj1]), cuts, labelmll, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_inclusive_OF", 40, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj1]), cuts, labelmll, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_inclusive_ee", 40, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj1]), cuts, labelmll, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_inclusive_mm", 40, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj1]), cuts, labelmll, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_inclusive_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj1]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_inclusive_OF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj1]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_inclusive_ee", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj1]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_inclusive_mm", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj1]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.DYControlRegionNoMll]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.DYControlRegionNoMll]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.DYControlRegionNoMll]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.DYControlRegionNoMll]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DYJets_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.DYControlRegionNoMET]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DYJets_OF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.DYControlRegionNoMET]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DYJets_ee", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.DYControlRegionNoMET]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DYJets_mm", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.DYControlRegionNoMET]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbarcontrol_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbarcontrol_OF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbarcontrol_ee", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbarcontrol_mm", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.JetMETBaseline]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.JetMETBaseline]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.JetMETBaseline]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.JetMETBaseline]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_signal_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj2]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_signal_OF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj2]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_signal_ee", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj2]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_signal_mm", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj2]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_SF", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.SF]), cuts, labelnjet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_OF", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.OF]), cuts, labelnjet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_ee", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.ee]), cuts, labelnjet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_mm", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.mm]), cuts, labelnjet, 1)

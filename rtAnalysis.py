############################################################
############################################################
##          _          _            _              _      ##
##         /\ \       /\ \         /\ \           /\ \    ##
##        /  \ \     /  \ \____   /  \ \         /  \ \   ##
##       / /\ \ \   / /\ \_____\ / /\ \_\       / /\ \ \  ##
##      / / /\ \_\ / / /\/___  // / /\/_/      / / /\ \_\ ##
##     / /_/_ \/_// / /   / / // / / ______   / /_/_ \/_/ ##
##    / /____/\  / / /   / / // / / /\_____\ / /____/\    ##
##   / /\____\/ / / /   / / // / /  \/____ // /\____\/    ##
##  / / /______ \ \ \__/ / // / /_____/ / // / /______    ##
## / / /_______\ \ \___\/ // / /______\/ // / /_______\   ##
## \/__________/  \/_____/ \/___________/ \/__________/   ##
############################################################
############################################################
                                                       

import ROOT as r
from ROOT import gROOT, TCanvas, TFile, TF1, TGraphAsymmErrors
import math,sys,optparse, copy, subprocess
import math


import include.LeptonSF
import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Rounder    as rounder


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\032.4m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'




def makeTable(DATAnumee, DATAnummm, DATAnumOF, DATAdenee, DATAdenmm, DATAdenOF, effMlleevalue, effMllmmvalue , effMllOFvalue, eevale, mmvale, emvale,  rt, rtUnc, rtSyst, dataMC):
    line0 = '  \hline'
    line1 = '  ee &' 
    line2 = '  \mu\mu   &' 
    line3 = '  e\mu   &' 
    line4 = '  && RT   '
    if dataMC == 'DATA':
        line0 += ' Data & nominator  & den & $\epsilon_{\mathrm{trigger}}$ $\pm \sigma_{stat}}$ \\\\'
    else:
        line0 += ' MC & nominator  & den & $\epsilon_{\mathrm{trigger}}$ $\pm \sigma_{stat}}$ \\\\'
    line1 += ' %.f & %.f &    %.3f $\\pm$ %.3f       %s' %(DATAnumee.Integral(), DATAdenee.Integral(), effMlleevalue  , eevale  ,   '\\\\')
    line2 += ' %.f & %.f &    %.3f $\\pm$ %.3f       %s' %(DATAnummm.Integral(), DATAdenmm.Integral(), effMllmmvalue  , mmvale  ,   '\\\\')
    line3 += ' %.f & %.f &    %.3f $\\pm$ %.3f       %s' %(DATAnumOF.Integral(), DATAdenOF.Integral(), effMllOFvalue  , emvale  ,   '\\\\')
    line4 += ' %.3f $\\pm$ %.3f  $\\pm$ %.3f     %s' %(rt, rtUnc, rtSyst  ,   '\\\\')
    line0 += '\\hline'; line4 += '\\hline';
                                                                                                                                                                                     
    helper.ensureDirectory('plots/rt/%s/'%lumi_str)
    helper.ensureDirectory('plots/rt/%s/tables/'%lumi_str)
    compTableFile = open('plots/rt/%s/tables/resultTable_%s%s_%s.txt'%(lumi_str, str(lumi).replace('.','p'), "rt", dataMC),'w')
    compTableFile.write(line0+'\n')
    compTableFile.write(line1+'\n')
    compTableFile.write(line2+'\n')                                                                                             
    compTableFile.write(line3+'\n')                                                                                             
    compTableFile.write(line4+'\n')                                                                                             
    print line0
    print line1
    print line2                                                                                                                                                                      
    print line3                                                                                                                                
    print line4                                                                                                                                                                                        


############################################################
def getRT(eff_ee, unc_ee, eff_mm, unc_mm, eff_em, unc_em):
    if(eff_em == 0 ):return [0, 0, 0]
    if(eff_ee == 0 ):return [0, 0, 0]
    if(eff_mm == 0 ):return [0, 0, 0]
    RT = math.sqrt(eff_ee*eff_mm)/eff_em
    uncRTee = (math.sqrt(eff_mm)/eff_em)*0.5*(unc_ee/math.sqrt(eff_ee)) 
    uncRTmm = (math.sqrt(eff_ee)/eff_em)*0.5*(unc_mm/math.sqrt(eff_mm)) 
    uncRTem = math.sqrt(eff_ee*eff_mm)*((unc_em)/(eff_em*eff_em)) 
    uncRT = math.sqrt(uncRTee * uncRTee + uncRTmm * uncRTmm + uncRTem * uncRTem)
    uncsysRTee = (math.sqrt(eff_mm)/eff_em)*0.5*(0.05*eff_ee/math.sqrt(eff_ee)) 
    uncsysRTmm = (math.sqrt(eff_ee)/eff_em)*0.5*(0.05*eff_mm/math.sqrt(eff_mm)) 
    uncsysRTem = math.sqrt(eff_ee*eff_mm)*((0.05*eff_em)/(eff_em*eff_em)) 
    uncsysRT = math.sqrt(uncsysRTee * uncsysRTee + uncsysRTmm * uncsysRTmm + uncsysRTem * uncsysRTem)

    return [RT, uncRT, uncsysRT]


def RT(ref, eff_ee, eff_mm, eff_em):

    RTratio = copy.deepcopy(ref)
    RTratio.SetName("RT_" + ref.GetName())
    RTratioee = copy.deepcopy(ref)
    RTratioee.SetName("RT_ee" + ref.GetName())
    RTratiomm = copy.deepcopy(ref)
    RTratiomm.SetName("RT_mm" + ref.GetName())
    RTratioem = copy.deepcopy(ref)
    RTratioem.SetName("RT_em" + ref.GetName())

    RTratio.GetYaxis().SetTitle("RT")
    RTratio.GetXaxis().SetTitle(ref.GetXaxis().GetTitle())
    RTratioee.GetYaxis().SetTitle("RT")
    RTratioee.GetXaxis().SetTitle(ref.GetXaxis().GetTitle())
    RTratiomm.GetYaxis().SetTitle("RT")
    RTratiomm.GetXaxis().SetTitle(ref.GetXaxis().GetTitle())
    RTratioem.GetYaxis().SetTitle("RT")
    RTratioem.GetXaxis().SetTitle(ref.GetXaxis().GetTitle())
    RTratioee.Reset()
    RTratiomm.Reset()
    RTratioem.Reset()

    xee = eff_ee.GetX()
    yee = eff_ee.GetY()
    eyhee = eff_ee.GetEYhigh()
    eylee = eff_ee.GetEYlow()
    for i in range(0, eff_ee.GetN()):
        j = RTratioee.FindBin(xee[i])
        RTratioee.SetBinContent(j, yee[i])
        RTratioee.SetBinError(j, max(eyhee[i], eylee[i]))

    xmm = eff_mm.GetX()
    ymm = eff_mm.GetY()
    eyhmm = eff_mm.GetEYhigh()
    eylmm = eff_mm.GetEYlow()
    for i in range(0, eff_mm.GetN()):
        j = RTratiomm.FindBin(xmm[i])
        RTratiomm.SetBinContent(j, ymm[i])
        RTratiomm.SetBinError(j, max(eyhmm[i], eylmm[i]))

    xem = eff_em.GetX()
    yem = eff_em.GetY()
    eyhem = eff_em.GetEYhigh()
    eylem = eff_em.GetEYlow()
    for i in range(0, eff_em.GetN()):
        j = RTratioem.FindBin(xem[i])
        RTratioem.SetBinContent(j, yem[i])
        RTratioem.SetBinError(j, max(eyhem[i], eylem[i]))

    for i in range(1, RTratioee.GetNbinsX()+2):
        [rt, uncrt, uncsys] = getRT(RTratioee.GetBinContent(i), RTratioee.GetBinError(i), RTratiomm.GetBinContent(i), RTratiomm.GetBinError(i), RTratioem.GetBinContent(i), RTratioem.GetBinError(i)) 
        RTratio.SetBinContent(i, rt)
        RTratio.SetBinError(i, uncrt)
    
    #SetOwnership(RTratio, 0 )   # 0 = release (not keep), 1 = keep
    #SetOwnership(RTratioee, 0 )   # 0 = release (not keep), 1 = keep
    #SetOwnership(RTratiomm, 0 )   # 0 = release (not keep), 1 = keep
    #SetOwnership(RTratioem, 0 )   # 0 = release (not keep), 1 = keep
    #del RTratioee
    #del RTratiomm
    #del RTratioem
    return RTratio


def getTriggerEffs(num, den, dataMC):
    trig = den.Clone("trig_" + den.GetName())
    trig.Reset()
    print dataMC, "################################################################################################################################" 
    if dataMC == "data":
        option = 'v'
    else:
        option = 'pois'
            
    errs = TGraphAsymmErrors(num, den, option)
   # x = errs.GetX()
   # y = errs.GetY()
   # eyh = errs.GetEYhigh()
   # eyl = errs.GetEYlow()
   # for i in range(1, num.GetNbinsX()+1):
   #     if(x[i] != 0 or y[i] != 0):
   #         trig.SetBinContent(i, y[i])
   #         trig.SetBinError(i, max(eyh[i], eyl[i]))

    return errs


def saveInFile(theFile, measuredValueMC, measuredValueUncMC, measuredValueSystMC, measuredValueData, measuredValueUncData, measuredValueSystData):


    foutput = open(theFile + "_aux", 'w')
    for line in open(theFile).readlines():
        if line.find("rt") != -1:
            if line.find("DATA") != -1:
                foutput.write('rt          region          DATA        %.4f      %0.4f       %.4f\n'%(measuredValueData, measuredValueUncData, measuredValueSystData))
            if line.find("MC") != -1:
                foutput.write('rt          region          MC          %.4f      %0.4f       %.4f\n'%(measuredValueMC, measuredValueUncMC, measuredValueSystMC))
        else:
            foutput.write(line)

    foutput.close()
    subprocess.call(['mv ' + theFile + "_aux " + theFile], shell=True)



##Main body of the analysis
if __name__ == '__main__':

    print bcolors.HEADER
    print '#######################################################################'
    print '                  Starting RT analysis...                          '
    print '#######################################################################' + bcolors.ENDC

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples_forRt.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()

    theFile = opts.ingredientsFile
    print bcolors.HEADER + '[RTAnalysis] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    makeDependencyPlots = True 


    dyDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50HTskimmed', 'DYJetsToLL_M50_HT100to200','DYJetsToLL_M50_HT200to400', 'DYJetsToLL_M50_HT400to600', 'DYJetsToLL_M50_HT600to800', 'DYJetsToLL_M50_HT800to1200', 'DYJetsToLL_M50_HT1200to2500']
    ttDatasets = ['TTJets','TTJets_SingleLeptonFromT']
    stDatasets = ['TToLeptons_sch', 'T_tch_powheg', 'TBar_tch_powheg', 'T_tWch_noFullHad', 'TBar_tWch_noFullHad_ext', 'tZq_ll']
    ttzDatasets = ['TTZ_LO', 'TTLLJets_m1to10', 'TWZ', 'TTWToLNu', 'TTW_LO',  'TTWZ', 'TTZH', 'TTZZ', 'TTGJets']
    zz2lDatasets = ['ZZTo2L2Nu', 'GluGluToContinToZZTo2e2nu', 'GluGluToContinToZZTo2mu2nu']
    zz4lDatasets = ['ZZTo4L',  'GGHZZ4L_ext', 'VBF_HToZZTo4L']
    wwDatasets = ['WWTo2L2Nu', 'WWTo1L1Nu2Q', 'WJetsToLNu_LO']
    wzDatasets = ['WZTo3LNu_amcatnlo']
    raDatasets = ['WWW_4F', 'WZG', 'WZZ', 'ZZZ', 'TTHnobb_pow','TTTT' ]
    mcDatasets = zz4lDatasets + zz2lDatasets + ttzDatasets + raDatasets + wwDatasets +wzDatasets + stDatasets+  ttDatasets + dyDatasets                           



    daDatasets16B = ['DoubleEG_Run2016B_03Feb2017_ver2_v2_runs_273150_275376',
                   'DoubleMuon_Run2016B_03Feb2017_ver2_v2_runs_273150_275376', 
                   'MuonEG_Run2016B_03Feb2017_ver2_v2_runs_273150_275376'] 
 
    daDatasets16C = ['DoubleEG_Run2016C_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016C_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016C_03Feb2017_v1_runs_271036_284044']
    
    daDatasets16D = ['DoubleEG_Run2016D_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016D_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016D_03Feb2017_v1_runs_271036_284044']
 
    daDatasets16E = ['DoubleEG_Run2016E_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016E_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016E_03Feb2017_v1_runs_271036_284044']
 
    daDatasets16F = ['DoubleEG_Run2016F_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016F_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016F_03Feb2017_v1_runs_271036_284044'] 

    daDatasets16G = ['DoubleEG_Run2016G_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016G_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016G_03Feb2017_v1_runs_271036_284044']

    daDatasets16H = ['DoubleEG_Run2016H_03Feb2017_ver2_v1_runs_281085_284035',
                   'DoubleEG_Run2016H_03Feb2017_ver3_v1_runs_284036_284044',
                   'DoubleMuon_Run2016H_03Feb2017_ver2_v1_runs_281085_284035',
                   'DoubleMuon_Run2016H_03Feb2017_ver3_v1_runs_284036_284044',
                   'MuonEG_Run2016H_03Feb2017_ver2_v1_runs_281085_284035', 
                   'MuonEG_Run2016H_03Feb2017_ver3_v1_runs_284036_284044'] 

                                                                                  
    daDatasets17B = ['DoubleEG_Run2017B_17Nov2017_v1_runs_297046_299329',     
                   'DoubleMuon_Run2017B_17Nov2017_v1_runs_297046_299329',
                   'MuonEG_Run2017B_17Nov2017_v1_runs_297046_299329']    
                                                                              
    daDatasets17C = ['DoubleEG_Run2017C_17Nov2017_v1_runs_299368_302029',
                   'DoubleMuon_Run2017C_17Nov2017_v1_runs_299368_302029',
                   'MuonEG_Run2017C_17Nov2017_v1_runs_299368_302029']    
    
    daDatasets17D = ['DoubleEG_Run2017D_17Nov2017_v1_runs_302030_303434',
                   'DoubleMuon_Run2017D_17Nov2017_v1_runs_302030_303434',
                   'MuonEG_Run2017D_17Nov2017_v1_runs_302030_303434']    
                                                                              
    daDatasets17E = ['DoubleEG_Run2017E_17Nov2017_v1_runs_303824_304797',
                   'DoubleMuon_Run2017E_17Nov2017_v1_runs_303824_304797',
                   'MuonEG_Run2017E_17Nov2017_v1_runs_303824_304797']    
                                                                              
    daDatasets17F = ['DoubleEG_Run2017F_17Nov2017_v1_runs_305040_306462',
                  'DoubleMuon_Run2017F_17Nov2017_v1_runs_305040_306462',
                  'MuonEG_Run2017F_17Nov2017_v1_runs_305040_306462']          



    daDatasets16 = daDatasets16B + daDatasets16C + daDatasets16D +daDatasets16E + daDatasets16F + daDatasets16G + daDatasets16H     
    daDatasets17 = daDatasets17B + daDatasets17C + daDatasets17D +daDatasets17E + daDatasets17F       
    treeMC = Sample.Tree(helper.selectSamples("samples.dat", mcDatasets, 'MC'), 'MC'  , 0)
    treeDA16 = Sample.Tree(helper.selectSamples("samplesEdge.dat", daDatasets16, 'DA'), 'DATA', 1)
    treeDA17 = Sample.Tree(helper.selectSamples("samples.dat", daDatasets17, 'DA'), 'DATA', 1)

    print bcolors.HEADER + '[RTAnalysis] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    
    lumi = 41.9 ; maxrun = 999999; lumi_str = '41.9invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()
    labelx = 'm_{ll} [GeV]'
    labelmet = 'E_{T}^{miss} [GeV]'
    labelmt2 = 'mT2 [GeV]'
    labelpt1 = 'p_{T} leading [GeV]'
    labelpt2 = 'p_{T} subleading [GeV]'
    labeleta1 = '#eta leading'
    labeleta2 = '#eta subleading'
    mllbins = [20, 70, 81, 101, 111, 300]
    metbins = [0, 50, 100, 150]
    mt2bins = [0, 20, 40, 60, 80, 100, 160]
    ptbins = [ 25, 30, 40, 50, 70, 100, 150, 200]
    etabins = [-3.0, -2.4, -1.8, -1.2, -0.6, 0, 0.6, 1.2, 1.8, 2.4, 3.0]
    specialcut = ''
    kf = 'noKFactor'
    DATA16denMlleevalue =     treeDA16.getTH1F(lumi, "DATA16deneevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.ee]), '', labelx,"1", kf)
    DATA16numMlleevalue =     treeDA16.getTH1F(lumi, "DATA16numeevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.ee]), '', labelx,"1", kf)
    DATA16denMllmmvalue =     treeDA16.getTH1F(lumi, "DATA16denmmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.mm]), '', labelx,"1", kf)
    DATA16numMllmmvalue =     treeDA16.getTH1F(lumi, "DATA16nummmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.mm]), '', labelx,"1", kf)
    DATA16denMllSFvalue =     treeDA16.getTH1F(lumi, "DATA16denSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.SF]), '', labelx,"1", kf)
    DATA16numMllSFvalue =     treeDA16.getTH1F(lumi, "DATA16numSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.SF]), '', labelx,"1", kf)
    DATA16denMllOFvalue =     treeDA16.getTH1F(lumi, "DATA16denOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.OF]), '', labelx,"1", kf)
    DATA16numMllOFvalue =     treeDA16.getTH1F(lumi, "DATA16numOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.OF]), '', labelx,"1", kf)
   
    DATA17denMlleevalue =     treeDA17.getTH1F(lumi, "DATA17deneevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.ee]), '', labelx,"1", kf)
    DATA17numMlleevalue =     treeDA17.getTH1F(lumi, "DATA17numeevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.ee]), '', labelx,"1", kf)
    DATA17denMllmmvalue =     treeDA17.getTH1F(lumi, "DATA17denmmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.mm]), '', labelx,"1", kf)
    DATA17numMllmmvalue =     treeDA17.getTH1F(lumi, "DATA17nummmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.mm]), '', labelx,"1", kf)
    DATA17denMllSFvalue =     treeDA17.getTH1F(lumi, "DATA17denSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.SF]), '', labelx,"1", kf)
    DATA17numMllSFvalue =     treeDA17.getTH1F(lumi, "DATA17numSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.SF]), '', labelx,"1", kf)
    DATA17denMllOFvalue =     treeDA17.getTH1F(lumi, "DATA17denOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.OF]), '', labelx,"1", kf)
    DATA17numMllOFvalue =     treeDA17.getTH1F(lumi, "DATA17numOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.OF]), '', labelx,"1", kf)

    MCdenMlleevalue =     treeMC.getTH1F(lumi, "MCdeneevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.ee]), '', labelx,"1", kf)
    MCnumMlleevalue =     treeMC.getTH1F(lumi, "MCnumeevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.ee]), '', labelx,"1", kf)
    MCdenMllmmvalue =     treeMC.getTH1F(lumi, "MCdenmmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.mm]), '', labelx,"1", kf)
    MCnumMllmmvalue =     treeMC.getTH1F(lumi, "MCnummmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.mm]), '', labelx,"1", kf)
    MCdenMllSFvalue =     treeMC.getTH1F(lumi, "MCdenSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.SF]), '', labelx,"1", kf)
    MCnumMllSFvalue =     treeMC.getTH1F(lumi, "MCnumSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.SF]), '', labelx,"1", kf)
    MCdenMllOFvalue =     treeMC.getTH1F(lumi, "MCdenOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.OF]), '', labelx,"1", kf)
    MCnumMllOFvalue =     treeMC.getTH1F(lumi, "MCnumOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.OF]), '', labelx,"1", kf)
    
    effMlleevalueDATA16 = getTriggerEffs(DATA16numMlleevalue, DATA16denMlleevalue, "data")
    print "effMlleevalueDATA ", effMlleevalueDATA16.GetY()[0]
    effMllmmvalueDATA16 = getTriggerEffs(DATA16numMllmmvalue, DATA16denMllmmvalue, "data")
    print "effMllmmvalueDATA ", effMllmmvalueDATA16.GetY()[0]
    effMllSFvalueDATA16 = getTriggerEffs(DATA16numMllSFvalue, DATA16denMllSFvalue, "data")
    print "effMllSFvalueDATA ", effMllSFvalueDATA16.GetY()[0]
    effMllOFvalueDATA16 = getTriggerEffs(DATA16numMllOFvalue, DATA16denMllOFvalue, "data")
    print "effMllOFvalueDATA ", effMllOFvalueDATA16.GetY()[0]                                  
    effMlleevalueDATA17 = getTriggerEffs(DATA16numMlleevalue, DATA17denMlleevalue, "data")
    print "effMlleevalueDATA17 ", effMlleevalueDATA17.GetY()[0]
    effMllmmvalueDATA17 = getTriggerEffs(DATA17numMllmmvalue, DATA17denMllmmvalue, "data")
    print "effMllmmvalueDATA17 ", effMllmmvalueDATA17.GetY()[0]
    effMllSFvalueDATA17 = getTriggerEffs(DATA17numMllSFvalue, DATA17denMllSFvalue, "data")
    print "effMllSFvalueDATA17 ", effMllSFvalueDATA17.GetY()[0]
    effMllOFvalueDATA17 = getTriggerEffs(DATA17numMllOFvalue, DATA17denMllOFvalue, "data")
    print "effMllOFvalueDATA17 ", effMllOFvalueDATA17.GetY()[0]                                  
    
    
    
    print "MCnumMlleevalue ", MCnumMlleevalue
    print "MCdenMlleevalue ", MCdenMlleevalue
    print "MCnumMllmmvalue ", MCnumMllmmvalue
    print "MCdenMllmmvalue ", MCdenMllmmvalue
    print "MCnumMllSFvalue ", MCnumMllSFvalue
    print "MCdenMllSFvalue ", MCdenMllSFvalue
    print "MCnumMllOFvalue ", MCnumMllOFvalue
    print "MCdenMllOFvalue ", MCdenMllOFvalue
    
    
    effMlleevalueMC   = getTriggerEffs(MCnumMlleevalue, MCdenMlleevalue, "MC")
    print "effMlleevalueMC ", effMlleevalueMC.GetY()[0]
    effMllmmvalueMC   = getTriggerEffs(MCnumMllmmvalue, MCdenMllmmvalue, "MC")
    print "effMllmmvalueMC ", effMllmmvalueMC.GetY()[0]
    effMllSFvalueMC   = getTriggerEffs(MCnumMllSFvalue, MCdenMllSFvalue, "MC")
    print "effMllSFvalueMC ", effMllSFvalueMC.GetY()[0]
    effMllOFvalueMC   = getTriggerEffs(MCnumMllOFvalue, MCdenMllOFvalue, "MC")
    print "effMllOFvalueMC ", effMllOFvalueMC.GetY()[0]

    eevalda16_   = effMlleevalueDATA16.GetY()
    eevalda16eh_ = effMlleevalueDATA16.GetEYhigh()
    eevalda16el_ = effMlleevalueDATA16.GetEYlow()
    eevalda16    = eevalda16_[0]
    eevalda16e   = max(eevalda16eh_[0], eevalda16el_[0])
    mmvalda16_   = effMllmmvalueDATA16.GetY()
    mmvalda16eh_ = effMllmmvalueDATA16.GetEYhigh()
    mmvalda16el_ = effMllmmvalueDATA16.GetEYlow()
    mmvalda16    = mmvalda16_[0]
    mmvalda16e   = max(mmvalda16eh_[0], mmvalda16el_[0])
    emvalda16_   = effMllOFvalueDATA16.GetY()
    emvalda16eh_ = effMllOFvalueDATA16.GetEYhigh()
    emvalda16el_ = effMllOFvalueDATA16.GetEYlow()
    emvalda16    = emvalda16_[0]
    emvalda16e   = max(emvalda16eh_[0], emvalda16el_[0])

    eevalda17_   = effMlleevalueDATA17.GetY()
    eevalda17eh_ = effMlleevalueDATA17.GetEYhigh()
    eevalda17el_ = effMlleevalueDATA17.GetEYlow()
    eevalda17    = eevalda17_[0]
    eevalda17e   = max(eevalda17eh_[0], eevalda17el_[0])
    mmvalda17_   = effMllmmvalueDATA17.GetY()
    mmvalda17eh_ = effMllmmvalueDATA17.GetEYhigh()
    mmvalda17el_ = effMllmmvalueDATA17.GetEYlow()
    mmvalda17    = mmvalda17_[0]
    mmvalda17e   = max(mmvalda17eh_[0], mmvalda17el_[0])
    emvalda17_   = effMllOFvalueDATA17.GetY()
    emvalda17eh_ = effMllOFvalueDATA17.GetEYhigh()
    emvalda17el_ = effMllOFvalueDATA17.GetEYlow()
    emvalda17    = emvalda17_[0]
    emvalda17e   = max(emvalda17eh_[0], emvalda17el_[0])

    eevalmc_   = effMlleevalueMC.GetY()
    eevalmceh_ = effMlleevalueMC.GetEYhigh()
    eevalmcel_ = effMlleevalueMC.GetEYlow()
    eevalmc    = eevalmc_[0]
    eevalmce   = max(eevalmceh_[0], eevalmcel_[0])
    mmvalmc_   = effMllmmvalueMC.GetY()
    mmvalmceh_ = effMllmmvalueMC.GetEYhigh()
    mmvalmcel_ = effMllmmvalueMC.GetEYlow()
    mmvalmc    = mmvalmc_[0]
    mmvalmce   = max(mmvalmceh_[0], mmvalmcel_[0])
    emvalmc_   = effMllOFvalueMC.GetY()
    emvalmceh_ = effMllOFvalueMC.GetEYhigh()
    emvalmcel_ = effMllOFvalueMC.GetEYlow()
    emvalmc    = emvalmc_[0]
    emvalmce   = max(emvalmceh_[0], emvalmcel_[0])

    [da16rt, da16uncrt, da16systrt] = getRT(eevalda16, eevalda16e, mmvalda16, mmvalda16e, emvalda16, emvalda16e)
    [da17rt, da17uncrt, da17systrt] = getRT(eevalda17, eevalda17e, mmvalda17, mmvalda17e, emvalda17, emvalda17e)
    [mcrt, mcuncrt, mcsystrt] = getRT(eevalmc, eevalmce, mmvalmc, mmvalmce, emvalmc, emvalmce)
    print 'Measured RT value data 2016 ', da16rt, ' +/- ', da16uncrt, ' +/- ', da16systrt
    print 'Measured RT value data 2017 ', da17rt, ' +/- ', da17uncrt, ' +/- ', da17systrt
    print 'Measured RT value MC   ', mcrt, ' +/- ', mcuncrt, ' +/- ', mcsystrt
    saveInFile(theFile, mcrt, mcuncrt, mcsystrt, da17rt, da17uncrt, da17systrt)
    makeTable(DATA16numMlleevalue, DATA16numMllmmvalue, DATA16numMllOFvalue, DATA16denMlleevalue, DATA16denMllmmvalue, DATA16denMllOFvalue, eevalda16, mmvalda16, emvalda16, eevalda16e, mmvalda16e, emvalda16e, da16rt, da16uncrt, da16systrt, "DATA")
    makeTable(DATA17numMlleevalue, DATA17numMllmmvalue, DATA17numMllOFvalue, DATA17denMlleevalue, DATA17denMllmmvalue, DATA17denMllOFvalue, eevalda17, mmvalda17, emvalda17, eevalda17e, mmvalda17e, emvalda17e, da17rt, da17uncrt, da17systrt, "DATA")
    makeTable(MCnumMlleevalue, MCnumMllmmvalue, MCnumMllOFvalue, MCdenMlleevalue, MCdenMllmmvalue, MCdenMllOFvalue, eevalmc, mmvalmc, emvalmc, eevalmce, mmvalmce, emvalmce, mcrt, mcuncrt, mcsystrt, "MC")
    ######################## Calculation of total values #############################

    if makeDependencyPlots: 
        DATA16denMllee = treeDA16.getTH1F(lumi, "DATA16denee", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.ee]), '', labelx,"1", kf)
        DATA16numMllee = treeDA16.getTH1F(lumi, "DATA16numee", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.ee]), '', labelx,"1", kf)
        DATA16denMllmm = treeDA16.getTH1F(lumi, "DATA16denmm", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.mm]), '', labelx,"1", kf)
        DATA16numMllmm = treeDA16.getTH1F(lumi, "DATA16nummm", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.mm]), '', labelx,"1", kf)
        DATA16denMllSF = treeDA16.getTH1F(lumi, "DATA16denSF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.SF]), '', labelx,"1", kf)
        DATA16numMllSF = treeDA16.getTH1F(lumi, "DATA16numSF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.SF]), '', labelx,"1", kf)
        DATA16denMllOF = treeDA16.getTH1F(lumi, "DATA16denOF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.OF]), '', labelx,"1", kf)
        DATA16numMllOF = treeDA16.getTH1F(lumi, "DATA16numOF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.OF]), '', labelx,"1", kf)
        DATA16denpt1ee = treeDA16.getTH1F(lumi, "DATA16denpt1ee", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.ee]), '', labelpt1,"1", kf)
        DATA16numpt1ee = treeDA16.getTH1F(lumi, "DATA16numpt1ee", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.ee]), '', labelpt1,"1", kf)
        DATA16denpt1mm = treeDA16.getTH1F(lumi, "DATA16denpt1mm", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.mm]), '', labelpt1,"1", kf)
        DATA16numpt1mm = treeDA16.getTH1F(lumi, "DATA16numpt1mm", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.mm]), '', labelpt1,"1", kf)
        DATA16denpt1SF = treeDA16.getTH1F(lumi, "DATA16denpt1SF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.SF]), '', labelpt1,"1", kf)
        DATA16numpt1SF = treeDA16.getTH1F(lumi, "DATA16numpt1SF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.SF]), '', labelpt1,"1", kf)
        DATA16denpt1OF = treeDA16.getTH1F(lumi, "DATA16denpt1OF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.OF]), '', labelpt1,"1", kf)
        DATA16numpt1OF = treeDA16.getTH1F(lumi, "DATA16numpt1OF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.OF]), '', labelpt1,"1", kf)
        DATA16denpt2ee = treeDA16.getTH1F(lumi, "DATA16denpt2ee", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.ee]), '', labelpt2,"1", kf)
        DATA16numpt2ee = treeDA16.getTH1F(lumi, "DATA16numpt2ee", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.ee]), '', labelpt2,"1", kf)
        DATA16denpt2mm = treeDA16.getTH1F(lumi, "DATA16denpt2mm", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.mm]), '', labelpt2,"1", kf)
        DATA16numpt2mm = treeDA16.getTH1F(lumi, "DATA16numpt2mm", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.mm]), '', labelpt2,"1", kf)
        DATA16denpt2SF = treeDA16.getTH1F(lumi, "DATA16denpt2SF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.SF]), '', labelpt2,"1", kf)
        DATA16numpt2SF = treeDA16.getTH1F(lumi, "DATA16numpt2SF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.SF]), '', labelpt2,"1", kf)
        DATA16denpt2OF = treeDA16.getTH1F(lumi, "DATA16denpt2OF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.OF]), '', labelpt2,"1", kf)
        DATA16numpt2OF = treeDA16.getTH1F(lumi, "DATA16numpt2OF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.OF]), '', labelpt2,"1", kf)
        DATA16deneta1ee = treeDA16.getTH1F(lumi, "DATA16deneta1ee", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den16,  cuts.ee]), '', labeleta1,"1", kf)
        DATA16numeta1ee = treeDA16.getTH1F(lumi, "DATA16numeta1ee", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num16,  cuts.ee]), '', labeleta1,"1", kf)
        DATA16deneta1mm = treeDA16.getTH1F(lumi, "DATA16deneta1mm", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den16,  cuts.mm]), '', labeleta1,"1", kf)
        DATA16numeta1mm = treeDA16.getTH1F(lumi, "DATA16numeta1mm", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num16,  cuts.mm]), '', labeleta1,"1", kf)
        DATA16deneta1SF = treeDA16.getTH1F(lumi, "DATA16deneta1SF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den16,  cuts.SF]), '', labeleta1,"1", kf)
        DATA16numeta1SF = treeDA16.getTH1F(lumi, "DATA16numeta1SF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num16,  cuts.SF]), '', labeleta1,"1", kf)
        DATA16deneta1OF = treeDA16.getTH1F(lumi, "DATA16deneta1OF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den16,  cuts.OF]), '', labeleta1,"1", kf)
        DATA16numeta1OF = treeDA16.getTH1F(lumi, "DATA16numeta1OF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num16,  cuts.OF]), '', labeleta1,"1", kf)
        DATA16deneta2ee = treeDA16.getTH1F(lumi, "DATA16deneta2ee", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den16,  cuts.ee]), '', labeleta2,"1", kf)
        DATA16numeta2ee = treeDA16.getTH1F(lumi, "DATA16numeta2ee", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num16,  cuts.ee]), '', labeleta2,"1", kf)
        DATA16deneta2mm = treeDA16.getTH1F(lumi, "DATA16deneta2mm", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den16,  cuts.mm]), '', labeleta2,"1", kf)
        DATA16numeta2mm = treeDA16.getTH1F(lumi, "DATA16numeta2mm", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num16,  cuts.mm]), '', labeleta2,"1", kf)
        DATA16deneta2SF = treeDA16.getTH1F(lumi, "DATA16deneta2SF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den16,  cuts.SF]), '', labeleta2,"1", kf)
        DATA16numeta2SF = treeDA16.getTH1F(lumi, "DATA16numeta2SF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num16,  cuts.SF]), '', labeleta2,"1", kf)
        DATA16deneta2OF = treeDA16.getTH1F(lumi, "DATA16deneta2OF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den16,  cuts.OF]), '', labeleta2,"1", kf)
        DATA16numeta2OF = treeDA16.getTH1F(lumi, "DATA16numeta2OF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num16,  cuts.OF]), '', labeleta2,"1", kf)
        DATA16denMETee =  treeDA16.getTH1F(lumi, "DATA16deneevalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.ee]), '', labelmet,"1", kf)
        DATA16numMETee =  treeDA16.getTH1F(lumi, "DATA16numeevalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.ee]), '', labelmet,"1", kf)
        DATA16denMETmm =  treeDA16.getTH1F(lumi, "DATA16denmmvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.mm]), '', labelmet,"1", kf)
        DATA16numMETmm =  treeDA16.getTH1F(lumi, "DATA16nummmvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.mm]), '', labelmet,"1", kf)
        DATA16denMETSF =  treeDA16.getTH1F(lumi, "DATA16denSFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.SF]), '', labelmet,"1", kf)
        DATA16numMETSF =  treeDA16.getTH1F(lumi, "DATA16numSFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.SF]), '', labelmet,"1", kf)
        DATA16denMETOF =  treeDA16.getTH1F(lumi, "DATA16denOFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.OF]), '', labelmet,"1", kf)
        DATA16numMETOF =  treeDA16.getTH1F(lumi, "DATA16numOFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.OF]), '', labelmet,"1", kf)
        DATA16denmt2ee =  treeDA16.getTH1F(lumi, "DATA16deneevalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.ee]), '', labelmt2,"1", kf)
        DATA16nummt2ee =  treeDA16.getTH1F(lumi, "DATA16numeevalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.ee]), '', labelmt2,"1", kf)
        DATA16denmt2mm =  treeDA16.getTH1F(lumi, "DATA16denmmvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.mm]), '', labelmt2,"1", kf)
        DATA16nummt2mm =  treeDA16.getTH1F(lumi, "DATA16nummmvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.mm]), '', labelmt2,"1", kf)
        DATA16denmt2SF =  treeDA16.getTH1F(lumi, "DATA16denSFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.SF]), '', labelmt2,"1", kf)
        DATA16nummt2SF =  treeDA16.getTH1F(lumi, "DATA16numSFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.SF]), '', labelmt2,"1", kf)
        DATA16denmt2OF =  treeDA16.getTH1F(lumi, "DATA16denOFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den16,  cuts.OF]), '', labelmt2,"1", kf)
        DATA16nummt2OF =  treeDA16.getTH1F(lumi, "DATA16numOFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num16,  cuts.OF]), '', labelmt2,"1", kf)       

        DATA17denMllee = treeDA17.getTH1F(lumi, "DATA17denee", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.ee]), '', labelx,"1", kf)
        DATA17numMllee = treeDA17.getTH1F(lumi, "DATA17numee", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.ee]), '', labelx,"1", kf)
        DATA17denMllmm = treeDA17.getTH1F(lumi, "DATA17denmm", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.mm]), '', labelx,"1", kf)
        DATA17numMllmm = treeDA17.getTH1F(lumi, "DATA17nummm", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.mm]), '', labelx,"1", kf)
        DATA17denMllSF = treeDA17.getTH1F(lumi, "DATA17denSF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.SF]), '', labelx,"1", kf)
        DATA17numMllSF = treeDA17.getTH1F(lumi, "DATA17numSF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.SF]), '', labelx,"1", kf)
        DATA17denMllOF = treeDA17.getTH1F(lumi, "DATA17denOF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.OF]), '', labelx,"1", kf)
        DATA17numMllOF = treeDA17.getTH1F(lumi, "DATA17numOF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.OF]), '', labelx,"1", kf)
        DATA17denpt1ee = treeDA17.getTH1F(lumi, "DATA17denpt1ee", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.ee]), '', labelpt1,"1", kf)
        DATA17numpt1ee = treeDA17.getTH1F(lumi, "DATA17numpt1ee", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.ee]), '', labelpt1,"1", kf)
        DATA17denpt1mm = treeDA17.getTH1F(lumi, "DATA17denpt1mm", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.mm]), '', labelpt1,"1", kf)
        DATA17numpt1mm = treeDA17.getTH1F(lumi, "DATA17numpt1mm", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.mm]), '', labelpt1,"1", kf)
        DATA17denpt1SF = treeDA17.getTH1F(lumi, "DATA17denpt1SF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.SF]), '', labelpt1,"1", kf)
        DATA17numpt1SF = treeDA17.getTH1F(lumi, "DATA17numpt1SF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.SF]), '', labelpt1,"1", kf)
        DATA17denpt1OF = treeDA17.getTH1F(lumi, "DATA17denpt1OF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.OF]), '', labelpt1,"1", kf)
        DATA17numpt1OF = treeDA17.getTH1F(lumi, "DATA17numpt1OF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.OF]), '', labelpt1,"1", kf)
        DATA17denpt2ee = treeDA17.getTH1F(lumi, "DATA17denpt2ee", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.ee]), '', labelpt2,"1", kf)
        DATA17numpt2ee = treeDA17.getTH1F(lumi, "DATA17numpt2ee", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.ee]), '', labelpt2,"1", kf)
        DATA17denpt2mm = treeDA17.getTH1F(lumi, "DATA17denpt2mm", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.mm]), '', labelpt2,"1", kf)
        DATA17numpt2mm = treeDA17.getTH1F(lumi, "DATA17numpt2mm", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.mm]), '', labelpt2,"1", kf)
        DATA17denpt2SF = treeDA17.getTH1F(lumi, "DATA17denpt2SF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.SF]), '', labelpt2,"1", kf)
        DATA17numpt2SF = treeDA17.getTH1F(lumi, "DATA17numpt2SF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.SF]), '', labelpt2,"1", kf)
        DATA17denpt2OF = treeDA17.getTH1F(lumi, "DATA17denpt2OF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.OF]), '', labelpt2,"1", kf)
        DATA17numpt2OF = treeDA17.getTH1F(lumi, "DATA17numpt2OF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.OF]), '', labelpt2,"1", kf)
        DATA17deneta1ee = treeDA17.getTH1F(lumi, "DATA17deneta1ee", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den17,  cuts.ee]), '', labeleta1,"1", kf)
        DATA17numeta1ee = treeDA17.getTH1F(lumi, "DATA17numeta1ee", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num17,  cuts.ee]), '', labeleta1,"1", kf)
        DATA17deneta1mm = treeDA17.getTH1F(lumi, "DATA17deneta1mm", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den17,  cuts.mm]), '', labeleta1,"1", kf)
        DATA17numeta1mm = treeDA17.getTH1F(lumi, "DATA17numeta1mm", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num17,  cuts.mm]), '', labeleta1,"1", kf)
        DATA17deneta1SF = treeDA17.getTH1F(lumi, "DATA17deneta1SF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den17,  cuts.SF]), '', labeleta1,"1", kf)
        DATA17numeta1SF = treeDA17.getTH1F(lumi, "DATA17numeta1SF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num17,  cuts.SF]), '', labeleta1,"1", kf)
        DATA17deneta1OF = treeDA17.getTH1F(lumi, "DATA17deneta1OF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den17,  cuts.OF]), '', labeleta1,"1", kf)
        DATA17numeta1OF = treeDA17.getTH1F(lumi, "DATA17numeta1OF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num17,  cuts.OF]), '', labeleta1,"1", kf)
        DATA17deneta2ee = treeDA17.getTH1F(lumi, "DATA17deneta2ee", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den17,  cuts.ee]), '', labeleta2,"1", kf)
        DATA17numeta2ee = treeDA17.getTH1F(lumi, "DATA17numeta2ee", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num17,  cuts.ee]), '', labeleta2,"1", kf)
        DATA17deneta2mm = treeDA17.getTH1F(lumi, "DATA17deneta2mm", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den17,  cuts.mm]), '', labeleta2,"1", kf)
        DATA17numeta2mm = treeDA17.getTH1F(lumi, "DATA17numeta2mm", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num17,  cuts.mm]), '', labeleta2,"1", kf)
        DATA17deneta2SF = treeDA17.getTH1F(lumi, "DATA17deneta2SF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den17,  cuts.SF]), '', labeleta2,"1", kf)
        DATA17numeta2SF = treeDA17.getTH1F(lumi, "DATA17numeta2SF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num17,  cuts.SF]), '', labeleta2,"1", kf)
        DATA17deneta2OF = treeDA17.getTH1F(lumi, "DATA17deneta2OF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den17,  cuts.OF]), '', labeleta2,"1", kf)
        DATA17numeta2OF = treeDA17.getTH1F(lumi, "DATA17numeta2OF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num17,  cuts.OF]), '', labeleta2,"1", kf)
        DATA17denMETee =  treeDA17.getTH1F(lumi, "DATA17deneevalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.ee]), '', labelmet,"1", kf)
        DATA17numMETee =  treeDA17.getTH1F(lumi, "DATA17numeevalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.ee]), '', labelmet,"1", kf)
        DATA17denMETmm =  treeDA17.getTH1F(lumi, "DATA17denmmvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.mm]), '', labelmet,"1", kf)
        DATA17numMETmm =  treeDA17.getTH1F(lumi, "DATA17nummmvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.mm]), '', labelmet,"1", kf)
        DATA17denMETSF =  treeDA17.getTH1F(lumi, "DATA17denSFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.SF]), '', labelmet,"1", kf)
        DATA17numMETSF =  treeDA17.getTH1F(lumi, "DATA17numSFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.SF]), '', labelmet,"1", kf)
        DATA17denMETOF =  treeDA17.getTH1F(lumi, "DATA17denOFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.OF]), '', labelmet,"1", kf)
        DATA17numMETOF =  treeDA17.getTH1F(lumi, "DATA17numOFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.OF]), '', labelmet,"1", kf)
        DATA17denmt2ee =  treeDA17.getTH1F(lumi, "DATA17deneevalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.ee]), '', labelmt2,"1", kf)
        DATA17nummt2ee =  treeDA17.getTH1F(lumi, "DATA17numeevalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.ee]), '', labelmt2,"1", kf)
        DATA17denmt2mm =  treeDA17.getTH1F(lumi, "DATA17denmmvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.mm]), '', labelmt2,"1", kf)
        DATA17nummt2mm =  treeDA17.getTH1F(lumi, "DATA17nummmvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.mm]), '', labelmt2,"1", kf)
        DATA17denmt2SF =  treeDA17.getTH1F(lumi, "DATA17denSFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.SF]), '', labelmt2,"1", kf)
        DATA17nummt2SF =  treeDA17.getTH1F(lumi, "DATA17numSFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.SF]), '', labelmt2,"1", kf)
        DATA17denmt2OF =  treeDA17.getTH1F(lumi, "DATA17denOFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.OF]), '', labelmt2,"1", kf)
        DATA17nummt2OF =  treeDA17.getTH1F(lumi, "DATA17numOFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.OF]), '', labelmt2,"1", kf)       
                                                                                                                                                                                        
        MCdenMllee =   treeMC.getTH1F(lumi, "MCdenee", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.ee]), '', labelx,"1", kf)
        MCnumMllee =   treeMC.getTH1F(lumi, "MCnumee", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.ee]), '', labelx,"1", kf)
        MCdenMllmm =   treeMC.getTH1F(lumi, "MCdenmm", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.mm]), '', labelx,"1", kf)
        MCnumMllmm =   treeMC.getTH1F(lumi, "MCnummm", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.mm]), '', labelx,"1", kf)
        MCdenMllSF =   treeMC.getTH1F(lumi, "MCdenSF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.SF]), '', labelx,"1", kf)
        MCnumMllSF =   treeMC.getTH1F(lumi, "MCnumSF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.SF]), '', labelx,"1", kf)
        MCdenMllOF =   treeMC.getTH1F(lumi, "MCdenOF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.OF]), '', labelx,"1", kf)
        MCnumMllOF =   treeMC.getTH1F(lumi, "MCnumOF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.OF]), '', labelx,"1", kf)
                                                                                                                                                                                        
        MCdenpt1ee =   treeMC.getTH1F(lumi, "MCdenpt1ee", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.ee]), '', labelpt1,"1", kf)
        MCnumpt1ee =   treeMC.getTH1F(lumi, "MCnumpt1ee", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.ee]), '', labelpt1,"1", kf)
        MCdenpt1mm =   treeMC.getTH1F(lumi, "MCdenpt1mm", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.mm]), '', labelpt1,"1", kf)
        MCnumpt1mm =   treeMC.getTH1F(lumi, "MCnumpt1mm", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.mm]), '', labelpt1,"1", kf)
        MCdenpt1SF =   treeMC.getTH1F(lumi, "MCdenpt1SF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.SF]), '', labelpt1,"1", kf)
        MCnumpt1SF =   treeMC.getTH1F(lumi, "MCnumpt1SF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.SF]), '', labelpt1,"1", kf)
        MCdenpt1OF =   treeMC.getTH1F(lumi, "MCdenpt1OF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.OF]), '', labelpt1,"1", kf)
        MCnumpt1OF =   treeMC.getTH1F(lumi, "MCnumpt1OF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.OF]), '', labelpt1,"1", kf)
        
        MCdenpt2ee =   treeMC.getTH1F(lumi, "MCdenpt2ee", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.ee]), '', labelpt2,"1", kf)
        MCnumpt2ee =   treeMC.getTH1F(lumi, "MCnumpt2ee", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.ee]), '', labelpt2,"1", kf)
        MCdenpt2mm =   treeMC.getTH1F(lumi, "MCdenpt2mm", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.mm]), '', labelpt2,"1", kf)
        MCnumpt2mm =   treeMC.getTH1F(lumi, "MCnumpt2mm", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.mm]), '', labelpt2,"1", kf)
        MCdenpt2SF =   treeMC.getTH1F(lumi, "MCdenpt2SF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.SF]), '', labelpt2,"1", kf)
        MCnumpt2SF =   treeMC.getTH1F(lumi, "MCnumpt2SF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.SF]), '', labelpt2,"1", kf)
        MCdenpt2OF =   treeMC.getTH1F(lumi, "MCdenpt2OF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.OF]), '', labelpt2,"1", kf)
        MCnumpt2OF =   treeMC.getTH1F(lumi, "MCnumpt2OF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.OF]), '', labelpt2,"1", kf)
                                                                                                                                                                                        
        MCdeneta1ee =  treeMC.getTH1F(lumi, "MCdeneta1ee", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den17,  cuts.ee]), '', labeleta1,"1", kf)
        MCnumeta1ee =  treeMC.getTH1F(lumi, "MCnumeta1ee", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num17,  cuts.ee]), '', labeleta1,"1", kf)
        MCdeneta1mm =  treeMC.getTH1F(lumi, "MCdeneta1mm", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den17,  cuts.mm]), '', labeleta1,"1", kf)
        MCnumeta1mm =  treeMC.getTH1F(lumi, "MCnumeta1mm", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num17,  cuts.mm]), '', labeleta1,"1", kf)
        MCdeneta1SF =  treeMC.getTH1F(lumi, "MCdeneta1SF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den17,  cuts.SF]), '', labeleta1,"1", kf)
        MCnumeta1SF =  treeMC.getTH1F(lumi, "MCnumeta1SF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num17,  cuts.SF]), '', labeleta1,"1", kf)
        MCdeneta1OF =  treeMC.getTH1F(lumi, "MCdeneta1OF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den17,  cuts.OF]), '', labeleta1,"1", kf)
        MCnumeta1OF =  treeMC.getTH1F(lumi, "MCnumeta1OF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num17,  cuts.OF]), '', labeleta1,"1", kf)
        
        MCdeneta2ee =  treeMC.getTH1F(lumi, "MCdeneta2ee", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den17,  cuts.ee]), '', labeleta2,"1", kf)
        MCnumeta2ee =  treeMC.getTH1F(lumi, "MCnumeta2ee", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num17,  cuts.ee]), '', labeleta2,"1", kf)
        MCdeneta2mm =  treeMC.getTH1F(lumi, "MCdeneta2mm", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den17,  cuts.mm]), '', labeleta2,"1", kf)
        MCnumeta2mm =  treeMC.getTH1F(lumi, "MCnumeta2mm", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num17,  cuts.mm]), '', labeleta2,"1", kf)
        MCdeneta2SF =  treeMC.getTH1F(lumi, "MCdeneta2SF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den17,  cuts.SF]), '', labeleta2,"1", kf)
        MCnumeta2SF =  treeMC.getTH1F(lumi, "MCnumeta2SF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num17,  cuts.SF]), '', labeleta2,"1", kf)
        MCdeneta2OF =  treeMC.getTH1F(lumi, "MCdeneta2OF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den17,  cuts.OF]), '', labeleta2,"1", kf)
        MCnumeta2OF =  treeMC.getTH1F(lumi, "MCnumeta2OF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num17,  cuts.OF]), '', labeleta2,"1", kf)
                                                                                                                                                                                       
        MCdenMETee =   treeMC.getTH1F(lumi, "MCdeneevalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.ee]), '', labelmet,"1", kf)
        MCnumMETee =   treeMC.getTH1F(lumi, "MCnumeevalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.ee]), '', labelmet,"1", kf)
        MCdenMETmm =   treeMC.getTH1F(lumi, "MCdenmmvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.mm]), '', labelmet,"1", kf)
        MCnumMETmm =   treeMC.getTH1F(lumi, "MCnummmvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.mm]), '', labelmet,"1", kf)
        MCdenMETSF =   treeMC.getTH1F(lumi, "MCdenSFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.SF]), '', labelmet,"1", kf)
        MCnumMETSF =   treeMC.getTH1F(lumi, "MCnumSFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.SF]), '', labelmet,"1", kf)
        MCdenMETOF =   treeMC.getTH1F(lumi, "MCdenOFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.OF]), '', labelmet,"1", kf)
        MCnumMETOF =   treeMC.getTH1F(lumi, "MCnumOFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.OF]), '', labelmet,"1", kf)
        
        MCdenmt2ee =   treeMC.getTH1F(lumi, "MCdeneevalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.ee]), '', labelmt2,"1", kf)
        MCnummt2ee =   treeMC.getTH1F(lumi, "MCnumeevalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.ee]), '', labelmt2,"1", kf)
        MCdenmt2mm =   treeMC.getTH1F(lumi, "MCdenmmvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.mm]), '', labelmt2,"1", kf)
        MCnummt2mm =   treeMC.getTH1F(lumi, "MCnummmvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.mm]), '', labelmt2,"1", kf)
        MCdenmt2SF =   treeMC.getTH1F(lumi, "MCdenSFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.SF]), '', labelmt2,"1", kf)
        MCnummt2SF =   treeMC.getTH1F(lumi, "MCnumSFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.SF]), '', labelmt2,"1", kf)
        MCdenmt2OF =   treeMC.getTH1F(lumi, "MCdenOFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den17,  cuts.OF]), '', labelmt2,"1", kf)
        MCnummt2OF =   treeMC.getTH1F(lumi, "MCnumOFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num17,  cuts.OF]), '', labelmt2,"1", kf)  
        
        
        DATA16effMETee =  getTriggerEffs(DATA16numMETee,  DATA16denMETee,  "data")
        DATA16effMETmm =  getTriggerEffs(DATA16numMETmm,  DATA16denMETmm,  "data")    
        DATA16effMETSF =  getTriggerEffs(DATA16numMETSF,  DATA16denMETSF,  "data")    
        DATA16effMETOF =  getTriggerEffs(DATA16numMETOF,  DATA16denMETOF,  "data")
        DATA16effmt2ee =  getTriggerEffs(DATA16nummt2ee,  DATA16denmt2ee,  "data")
        DATA16effmt2mm =  getTriggerEffs(DATA16nummt2mm,  DATA16denmt2mm,  "data")
        DATA16effmt2SF =  getTriggerEffs(DATA16nummt2SF,  DATA16denmt2SF,  "data")
        DATA16effmt2OF =  getTriggerEffs(DATA16nummt2OF,  DATA16denmt2OF,  "data")
        DATA16effMllee =  getTriggerEffs(DATA16numMllee,  DATA16denMllee,  "data")
        DATA16effMllmm =  getTriggerEffs(DATA16numMllmm,  DATA16denMllmm,  "data")
        DATA16effMllSF =  getTriggerEffs(DATA16numMllSF,  DATA16denMllSF,  "data")
        DATA16effMllOF =  getTriggerEffs(DATA16numMllOF,  DATA16denMllOF,  "data")
        DATA16effpt1ee =  getTriggerEffs(DATA16numpt1ee,  DATA16denpt1ee,  "data")
        DATA16effpt1mm =  getTriggerEffs(DATA16numpt1mm,  DATA16denpt1mm,  "data")
        DATA16effpt1SF =  getTriggerEffs(DATA16numpt1SF,  DATA16denpt1SF,  "data")
        DATA16effpt1OF =  getTriggerEffs(DATA16numpt1OF,  DATA16denpt1OF,  "data")
        DATA16effpt2ee =  getTriggerEffs(DATA16numpt2ee,  DATA16denpt2ee,  "data")
        DATA16effpt2mm =  getTriggerEffs(DATA16numpt2mm,  DATA16denpt2mm,  "data")
        DATA16effpt2SF =  getTriggerEffs(DATA16numpt2SF,  DATA16denpt2SF,  "data")
        DATA16effpt2OF =  getTriggerEffs(DATA16numpt2OF,  DATA16denpt2OF,  "data")
        DATA16effeta1ee = getTriggerEffs(DATA16numeta1ee, DATA16deneta1ee, "data")
        DATA16effeta1mm = getTriggerEffs(DATA16numeta1mm, DATA16deneta1mm, "data")
        DATA16effeta1SF = getTriggerEffs(DATA16numeta1SF, DATA16deneta1SF, "data")
        DATA16effeta1OF = getTriggerEffs(DATA16numeta1OF, DATA16deneta1OF, "data")
        DATA16effeta2ee = getTriggerEffs(DATA16numeta2ee, DATA16deneta2ee, "data")
        DATA16effeta2mm = getTriggerEffs(DATA16numeta2mm, DATA16deneta2mm, "data")
        DATA16effeta2SF = getTriggerEffs(DATA16numeta2SF, DATA16deneta2SF, "data")
        DATA16effeta2OF = getTriggerEffs(DATA16numeta2OF, DATA16deneta2OF, "data")
       
        
        DATA17effMETee =  getTriggerEffs(DATA17numMETee,  DATA17denMETee,  "data")
        DATA17effMETmm =  getTriggerEffs(DATA17numMETmm,  DATA17denMETmm,  "data")
        DATA17effMETSF =  getTriggerEffs(DATA17numMETSF,  DATA17denMETSF,  "data")
        DATA17effMETOF =  getTriggerEffs(DATA17numMETOF,  DATA17denMETOF,  "data")
        DATA17effmt2ee =  getTriggerEffs(DATA17nummt2ee,  DATA17denmt2ee,  "data")
        DATA17effmt2mm =  getTriggerEffs(DATA17nummt2mm,  DATA17denmt2mm,  "data")
        DATA17effmt2SF =  getTriggerEffs(DATA17nummt2SF,  DATA17denmt2SF,  "data")
        DATA17effmt2OF =  getTriggerEffs(DATA17nummt2OF,  DATA17denmt2OF,  "data")
        DATA17effMllee =  getTriggerEffs(DATA17numMllee,  DATA17denMllee,  "data")
        DATA17effMllmm =  getTriggerEffs(DATA17numMllmm,  DATA17denMllmm,  "data")
        DATA17effMllSF =  getTriggerEffs(DATA17numMllSF,  DATA17denMllSF,  "data")
        DATA17effMllOF =  getTriggerEffs(DATA17numMllOF,  DATA17denMllOF,  "data")
        DATA17effpt1ee =  getTriggerEffs(DATA17numpt1ee,  DATA17denpt1ee,  "data")
        DATA17effpt1mm =  getTriggerEffs(DATA17numpt1mm,  DATA17denpt1mm,  "data")
        DATA17effpt1SF =  getTriggerEffs(DATA17numpt1SF,  DATA17denpt1SF,  "data")
        DATA17effpt1OF =  getTriggerEffs(DATA17numpt1OF,  DATA17denpt1OF,  "data")
        DATA17effpt2ee =  getTriggerEffs(DATA17numpt2ee,  DATA17denpt2ee,  "data")
        DATA17effpt2mm =  getTriggerEffs(DATA17numpt2mm,  DATA17denpt2mm,  "data")
        DATA17effpt2SF =  getTriggerEffs(DATA17numpt2SF,  DATA17denpt2SF,  "data")
        DATA17effpt2OF =  getTriggerEffs(DATA17numpt2OF,  DATA17denpt2OF,  "data")
        DATA17effeta1ee = getTriggerEffs(DATA17numeta1ee, DATA17deneta1ee, "data")
        DATA17effeta1mm = getTriggerEffs(DATA17numeta1mm, DATA17deneta1mm, "data")
        DATA17effeta1SF = getTriggerEffs(DATA17numeta1SF, DATA17deneta1SF, "data")
        DATA17effeta1OF = getTriggerEffs(DATA17numeta1OF, DATA17deneta1OF, "data")
        DATA17effeta2ee = getTriggerEffs(DATA17numeta2ee, DATA17deneta2ee, "data")
        DATA17effeta2mm = getTriggerEffs(DATA17numeta2mm, DATA17deneta2mm, "data")
        DATA17effeta2SF = getTriggerEffs(DATA17numeta2SF, DATA17deneta2SF, "data")
        DATA17effeta2OF = getTriggerEffs(DATA17numeta2OF, DATA17deneta2OF, "data")

        MCeffMETee =  getTriggerEffs(MCnumMETee, MCdenMETee, "MC")
        MCeffMETmm =  getTriggerEffs(MCnumMETmm, MCdenMETmm, "MC")  
        MCeffMETSF =  getTriggerEffs(MCnumMETSF, MCdenMETSF, "MC")  
        MCeffMETOF =  getTriggerEffs(MCnumMETOF, MCdenMETOF, "MC")
        MCeffmt2ee =  getTriggerEffs(MCnummt2ee, MCdenmt2ee, "MC")
        MCeffmt2mm =  getTriggerEffs(MCnummt2mm, MCdenmt2mm, "MC")
        MCeffmt2SF =  getTriggerEffs(MCnummt2SF, MCdenmt2SF, "MC")
        MCeffmt2OF =  getTriggerEffs(MCnummt2OF, MCdenmt2OF, "MC")
        MCeffMllee =  getTriggerEffs(MCnumMllee, MCdenMllee, "MC")
        MCeffMllmm =  getTriggerEffs(MCnumMllmm, MCdenMllmm, "MC")
        MCeffMllSF =  getTriggerEffs(MCnumMllSF, MCdenMllSF, "MC")
        MCeffMllOF =  getTriggerEffs(MCnumMllOF, MCdenMllOF, "MC")
        MCeffpt1ee =  getTriggerEffs(MCnumpt1ee, MCdenpt1ee, "MC")
        MCeffpt1mm =  getTriggerEffs(MCnumpt1mm, MCdenpt1mm, "MC")
        MCeffpt1SF =  getTriggerEffs(MCnumpt1SF, MCdenpt1SF, "MC")
        MCeffpt1OF =  getTriggerEffs(MCnumpt1OF, MCdenpt1OF, "MC")
        MCeffpt2ee =  getTriggerEffs(MCnumpt2ee, MCdenpt2ee, "MC")
        MCeffpt2mm =  getTriggerEffs(MCnumpt2mm, MCdenpt2mm, "MC")
        MCeffpt2SF =  getTriggerEffs(MCnumpt2SF, MCdenpt2SF, "MC")
        MCeffpt2OF =  getTriggerEffs(MCnumpt2OF, MCdenpt2OF, "MC")
        MCeffeta1ee = getTriggerEffs(MCnumeta1ee, MCdeneta1ee, "MC")
        MCeffeta1mm = getTriggerEffs(MCnumeta1mm, MCdeneta1mm, "MC")
        MCeffeta1SF = getTriggerEffs(MCnumeta1SF, MCdeneta1SF, "MC")
        MCeffeta1OF = getTriggerEffs(MCnumeta1OF, MCdeneta1OF, "MC")
        MCeffeta2ee = getTriggerEffs(MCnumeta2ee, MCdeneta2ee, "MC")
        MCeffeta2mm = getTriggerEffs(MCnumeta2mm, MCdeneta2mm, "MC")
        MCeffeta2SF = getTriggerEffs(MCnumeta2SF, MCdeneta2SF, "MC")
        MCeffeta2OF = getTriggerEffs(MCnumeta2OF, MCdeneta2OF, "MC")
        
        DATA16RTMET  = RT(DATA16numMETSF,  DATA16effMETee,  DATA16effMETmm,  DATA16effMETOF)
        DATA16RTmt2  = RT(DATA16nummt2SF,  DATA16effmt2ee,  DATA16effmt2mm,  DATA16effmt2OF)
        DATA16RTMll  = RT(DATA16numMllSF,  DATA16effMllee,  DATA16effMllmm,  DATA16effMllOF)
        DATA16RTpt1  = RT(DATA16numpt1SF,  DATA16effpt1ee,  DATA16effpt1mm,  DATA16effpt1OF)
        DATA16RTpt2  = RT(DATA16numpt2SF,  DATA16effpt2ee,  DATA16effpt2mm,  DATA16effpt2OF)
        DATA16RTeta1 = RT(DATA16numeta1SF, DATA16effeta1ee, DATA16effeta1mm, DATA16effeta1OF)
        DATA16RTeta2 = RT(DATA16numeta2SF, DATA16effeta2ee, DATA16effeta2mm, DATA16effeta2OF)
       
        DATA17RTMET  = RT(DATA17numMETSF,  DATA17effMETee,  DATA17effMETmm,  DATA17effMETOF)
        DATA17RTmt2  = RT(DATA17nummt2SF,  DATA17effmt2ee,  DATA17effmt2mm,  DATA17effmt2OF)
        DATA17RTMll  = RT(DATA17numMllSF,  DATA17effMllee,  DATA17effMllmm,  DATA17effMllOF)
        DATA17RTpt1  = RT(DATA17numpt1SF,  DATA17effpt1ee,  DATA17effpt1mm,  DATA17effpt1OF)
        DATA17RTpt2  = RT(DATA17numpt2SF,  DATA17effpt2ee,  DATA17effpt2mm,  DATA17effpt2OF)
        DATA17RTeta1 = RT(DATA17numeta1SF, DATA17effeta1ee, DATA17effeta1mm, DATA17effeta1OF)
        DATA17RTeta2 = RT(DATA17numeta2SF, DATA17effeta2ee, DATA17effeta2mm, DATA17effeta2OF)

        MCRTMET  = RT(MCnumMETSF,  MCeffMETee,  MCeffMETmm,  MCeffMETOF)
        MCRTmt2  = RT(MCnummt2SF,  MCeffmt2ee,  MCeffmt2mm,  MCeffmt2OF)
        MCRTMll  = RT(MCnumMllSF,  MCeffMllee,  MCeffMllmm,  MCeffMllOF)
        MCRTpt1  = RT(MCnumpt1SF,  MCeffpt1ee,  MCeffpt1mm,  MCeffpt1OF)
        MCRTpt2  = RT(MCnumpt2SF,  MCeffpt2ee,  MCeffpt2mm,  MCeffpt2OF)
        MCRTeta1 = RT(MCnumeta1SF, MCeffeta1ee, MCeffeta1mm, MCeffeta1OF)
        MCRTeta2 = RT(MCnumeta2SF, MCeffeta2ee, MCeffeta2mm, MCeffeta2OF)
        
        effRTMll = Canvas.Canvas('rt/%s/plot_rt_mll_'%(lumi_str), 'png,pdf',  0.6, 0.3, 0.8, 0.5)
        h_auxrtMll = r.TH1F("h_auxRTMll", "", 1, 0, 250)
        h_auxrtMll.GetYaxis().SetRangeUser(0.8, 1.2)
        h_auxrtMll.GetXaxis().SetRangeUser(0, 250)
        h_auxrtMll.GetXaxis().SetTitle(labelx)
        effRTMll.addLine(h_auxrtMll.GetXaxis().GetXmin(), da16rt, h_auxrtMll.GetXaxis().GetXmax(), da16rt ,r.kBlack)
        effRTMll.addLine(h_auxrtMll.GetXaxis().GetXmin(), da17rt, h_auxrtMll.GetXaxis().GetXmax(), da17rt ,r.kBlue)
        effRTMll.addHisto(h_auxrtMll, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
        effRTMll.addHisto(DATA16RTMll, 'PE,SAME', 'R_{T} data 2016', 'PL', r.kBlack , 1, 0)
        effRTMll.addHisto(DATA17RTMll, 'PE,SAME', 'R_{T} data 2017', 'PL', r.kBlue , 1, 0)
        effRTMll.addHisto(MCRTMll, 'PE,SAME', 'R_{T} MC', 'PL', r.kGreen , 1, 0)
        effRTMll.addBand(h_auxrtMll.GetXaxis().GetXmin(), da16rt-da16systrt, h_auxrtMll.GetXaxis().GetXmax(), da16rt+da16systrt, r.kOrange+6, 0.2)
        effRTMll.addBand(h_auxrtMll.GetXaxis().GetXmin(), da17rt-da17systrt, h_auxrtMll.GetXaxis().GetXmax(), da17rt+da17systrt, r.kPink+6, 0.2)
        effRTMll.addLatex (0.55, 0.3, 'Mean R_{T} data 2016: %.3f '%(da16rt))
        effRTMll.addLatex (0.55, 0.25, 'Mean R_{T} data 2017: %.3f '%(da17rt))
        effRTMll.addLatex (0.55, 0.2, 'Mean R_{T} MC: %.3f '%(mcrt))
        effRTMll.save(1, 1, 0, lumi, 0.8, 1.2)                                                            
        
        effRTMET = Canvas.Canvas('rt/%s/plot_rt_met_'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
        h_auxrtMET = r.TH1F("h_auxRTMET", "", 1, 0, 200)
        h_auxrtMET.GetYaxis().SetRangeUser(0.8, 1.2)
        h_auxrtMET.GetXaxis().SetRangeUser(0, 150)
        h_auxrtMET.GetXaxis().SetTitle(labelmet)
        effRTMET.addLine(h_auxrtMET.GetXaxis().GetXmin(), da16rt, h_auxrtMET.GetXaxis().GetXmax(), da16rt ,r.kBlack)
        effRTMET.addLine(h_auxrtMET.GetXaxis().GetXmin(), da17rt, h_auxrtMET.GetXaxis().GetXmax(), da17rt ,r.kBlue)
        effRTMET.addHisto(h_auxrtMET, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
        effRTMET.addHisto(DATA16RTMET, 'PE,SAME', 'R_{T} data 2016', 'PL', r.kBlack , 1, 0)
        effRTMET.addHisto(DATA17RTMET, 'PE,SAME', 'R_{T} data 2017', 'PL', r.kBlue , 1, 0)
        effRTMET.addHisto(MCRTMET, 'PE,SAME', 'R_{T} MC', 'PL', r.kGreen , 1, 0)
        effRTMET.addBand(h_auxrtMET.GetXaxis().GetXmin(), da16rt-da16systrt, h_auxrtMET.GetXaxis().GetXmax(), da16rt+da16systrt, r.kOrange+6, 0.2)
        effRTMET.addBand(h_auxrtMET.GetXaxis().GetXmin(), da17rt-da17systrt, h_auxrtMET.GetXaxis().GetXmax(), da17rt+da17systrt, r.kPink+6, 0.2)
        effRTMET.addLatex (0.6, 0.3, 'Mean R_{T} data 2016: %.3f '%(da16rt))
        effRTMET.addLatex (0.6, 0.25, 'Mean R_{T} data 2017: %.3f '%(da17rt))
        effRTMET.addLatex (0.6, 0.2, 'Mean R_{T} MC  : %.3f '%(mcrt))
        effRTMET.save(1, 1, 0, lumi, 0.8, 1.2)                                                                                                  
        
        
        effRTmt2 = Canvas.Canvas('rt/%s/plot_rt_mt2_'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
        h_auxrtmt2 = r.TH1F("h_auxRTMET", "", 1, 0, 160)
        h_auxrtmt2.GetYaxis().SetRangeUser(0.8, 1.2)
        h_auxrtmt2.GetXaxis().SetRangeUser(0, 160)
        h_auxrtmt2.GetXaxis().SetTitle(labelmt2)
        effRTmt2.addLine(h_auxrtmt2.GetXaxis().GetXmin(), da16rt, h_auxrtmt2.GetXaxis().GetXmax(), da16rt ,r.kBlack)
        effRTmt2.addLine(h_auxrtmt2.GetXaxis().GetXmin(), da17rt, h_auxrtmt2.GetXaxis().GetXmax(), da17rt ,r.kBlue)
        effRTmt2.addHisto(h_auxrtmt2, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
        effRTmt2.addHisto(DATA16RTmt2, 'PE,SAME', 'R_{T} data 2016', 'PL', r.kBlack , 1, 0)
        effRTmt2.addHisto(DATA17RTmt2, 'PE,SAME', 'R_{T} data 2017', 'PL', r.kBlue , 1, 0)
        effRTmt2.addHisto(MCRTmt2, 'PE,SAME', 'R_{T} MC', 'PL', r.kGreen , 1, 0)
        effRTmt2.addBand(h_auxrtmt2.GetXaxis().GetXmin(), da16rt-da16systrt, h_auxrtmt2.GetXaxis().GetXmax(), da16rt+da16systrt, r.kOrange+6, 0.2)
        effRTmt2.addBand(h_auxrtmt2.GetXaxis().GetXmin(), da17rt-da17systrt, h_auxrtmt2.GetXaxis().GetXmax(), da17rt+da17systrt, r.kPink+6, 0.2)
        effRTmt2.addLatex (0.6, 0.3, 'Mean R_{T} data 2016: %.3f '%(da16rt))
        effRTmt2.addLatex (0.6, 0.25, 'Mean R_{T} data 2017: %.3f '%(da17rt))
        effRTmt2.addLatex (0.6, 0.2, 'Mean R_{T} MC  : %.3f '%(mcrt))
        effRTmt2.save(1, 1, 0, lumi, 0.8, 1.2)                                                                                                 

        effRTpt1 = Canvas.Canvas('rt/%s/plot_rt_pt1_'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
        h_auxrtpt1 = r.TH1F("h_auxRTpt1", "", 1, 20, 150)
        h_auxrtpt1.GetYaxis().SetRangeUser(0.8, 1.2)
        h_auxrtpt1.GetXaxis().SetRangeUser(20, 150)
        h_auxrtpt1.GetXaxis().SetTitle(labelpt1)
        effRTpt1.addLine(h_auxrtpt1.GetXaxis().GetXmin(), da16rt, h_auxrtpt1.GetXaxis().GetXmax(), da16rt ,r.kBlack)
        effRTpt1.addLine(h_auxrtpt1.GetXaxis().GetXmin(), da17rt, h_auxrtpt1.GetXaxis().GetXmax(), da17rt ,r.kBlue)
        effRTpt1.addHisto(h_auxrtpt1, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
        effRTpt1.addHisto(DATA16RTpt1, 'PE,SAME', 'RT data 2016', 'PL', r.kBlack , 1, 0)
        effRTpt1.addHisto(DATA17RTpt1, 'PE,SAME', 'RT data 2017', 'PL', r.kBlue , 1, 0)
        effRTpt1.addHisto(MCRTpt1, 'PE,SAME', 'RT MC', 'PL', r.kGreen , 1, 0)
        effRTpt1.addBand(h_auxrtpt1.GetXaxis().GetXmin(), da16rt-da16systrt, h_auxrtpt1.GetXaxis().GetXmax(), da16rt+dasystrt, r.kOrange+6, 0.2)
        effRTpt1.addBand(h_auxrtpt1.GetXaxis().GetXmin(), da17rt-da17systrt, h_auxrtpt1.GetXaxis().GetXmax(), da17rt+dasystrt, r.kPink+6, 0.2)
        effRTpt1.addLatex (0.6, 0.3, 'Mean R_{T} data 2016: %.3f '%(da16rt))
        effRTpt1.addLatex (0.6, 0.25, 'Mean R_{T} data 2017: %.3f '%(da17rt))
        effRTpt1.addLatex (0.6, 0.2, 'Mean R_{T} MC  : %.3f '%(mcrt))
        effRTpt1.save(1, 1, 0, lumi, 0.8, 1.2)

        effRTpt2 = Canvas.Canvas('rt/%s/plot_rt_pt2_'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
        h_auxrtpt2 = r.TH1F("h_auxRTpt2", "", 1, 20, 150)
        h_auxrtpt2.GetYaxis().SetRangeUser(0.8, 1.2)
        h_auxrtpt2.GetXaxis().SetRangeUser(20, 150)
        h_auxrtpt2.GetXaxis().SetTitle(labelpt2)
        effRTpt2.addLine(h_auxrtpt2.GetXaxis().GetXmin(), da16rt, h_auxrtpt2.GetXaxis().GetXmax(), da16rt ,r.kBlack)
        effRTpt2.addLine(h_auxrtpt2.GetXaxis().GetXmin(), da17rt, h_auxrtpt2.GetXaxis().GetXmax(), da17rt ,r.kBlue)
        effRTpt2.addHisto(h_auxrtpt2, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
        effRTpt2.addHisto(DATA16RTpt2, 'PE,SAME', 'R_{T} data 2016', 'PL', r.kBlack , 1, 0)
        effRTpt2.addHisto(DATA17RTpt2, 'PE,SAME', 'R_{T} data 2017', 'PL', r.kBlue , 1, 0)
        effRTpt2.addHisto(MCRTpt2, 'PE,SAME', 'R_{T} MC', 'PL', r.kGreen , 1, 0)
        effRTpt2.addBand(h_auxrtpt2.GetXaxis().GetXmin(), da16rt-da16systrt, h_auxrtpt2.GetXaxis().GetXmax(), da16rt+da16systrt, r.kOrange+6, 0.2)
        effRTpt2.addBand(h_auxrtpt2.GetXaxis().GetXmin(), da17rt-da17systrt, h_auxrtpt2.GetXaxis().GetXmax(), da17rt+da17systrt, r.kPink+6, 0.2)
        effRTpt2.addLatex (0.6, 0.3, 'Mean R_{T} data 2016: %.3f '%(da16rt))
        effRTpt2.addLatex (0.6, 0.25, 'Mean R_{T} data 2017: %.3f '%(da17rt))
        effRTpt2.addLatex (0.6, 0.2, 'Mean R_{T} MC: %.3f '%(mcrt))
        effRTpt2.save(1, 1, 0, lumi, 0.8, 1.2)

        effRTeta1 = Canvas.Canvas('rt/%s/plot_rt_eta1_'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
        h_auxrteta1 = r.TH1F("h_auxRTeta1", "", 1, -2.4, 2.4)
        h_auxrteta1.GetYaxis().SetRangeUser(0.8, 1.2)
        h_auxrteta1.GetXaxis().SetRangeUser(-2.4, 2.4)
        h_auxrteta1.GetXaxis().SetTitle(labeleta1)
        effRTeta1.addLine(h_auxrteta1.GetXaxis().GetXmin(), da16rt, h_auxrteta1.GetXaxis().GetXmax(), da16rt ,r.kBlack)
        effRTeta1.addLine(h_auxrteta1.GetXaxis().GetXmin(), da17rt, h_auxrteta1.GetXaxis().GetXmax(), da17rt ,r.kBlue)
        effRTeta1.addHisto(h_auxrteta1, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
        effRTeta1.addHisto(DATA16RTeta1, 'PE,SAME', 'RT data 2016', 'PL', r.kBlack , 1, 0)
        effRTeta1.addHisto(DATA17RTeta1, 'PE,SAME', 'RT data 2017', 'PL', r.kBlue , 1, 0)
        effRTeta1.addHisto(MCRTeta1, 'PE,SAME', 'RT MC', 'PL', r.kGreen , 1, 0)
        effRTeta1.addBand(h_auxrteta1.GetXaxis().GetXmin(), da16rt-da16systrt, h_auxrteta1.GetXaxis().GetXmax(), da16rt+da16systrt, r.kOrange+6, 0.2)
        effRTeta1.addBand(h_auxrteta1.GetXaxis().GetXmin(), da17rt-da17systrt, h_auxrteta1.GetXaxis().GetXmax(), da17rt+da17systrt, r.kPink+6, 0.2)
        effRTeta1.addLatex (0.6, 0.3, 'Mean R_{T} data 2016: %.3f '%(da16rt))
        effRTeta1.addLatex (0.6, 0.25, 'Mean R_{T} data 2017: %.3f '%(da17rt))
        effRTeta1.addLatex (0.6, 0.2, 'Mean R_{T} MC  : %.3f '%(mcrt))
        effRTeta1.save(1, 1, 0, lumi, 0.8, 1.2)

        effRTeta2 = Canvas.Canvas('rt/%s/plot_rt_eta2_'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
        h_auxrteta2 = r.TH1F("h_auxRTeta2", "", 1, -2.4, 2.4)
        h_auxrteta2.GetYaxis().SetRangeUser(0.8, 1.2)
        h_auxrteta2.GetXaxis().SetRangeUser(-2.4, 2.4)
        h_auxrteta2.GetXaxis().SetTitle(labeleta2)
        effRTeta2.addLine(h_auxrteta2.GetXaxis().GetXmin(), da16rt, h_auxrteta2.GetXaxis().GetXmax(), da16rt ,r.kBlack)
        effRTeta2.addLine(h_auxrteta2.GetXaxis().GetXmin(), da17rt, h_auxrteta2.GetXaxis().GetXmax(), da17rt ,r.kBlue)
        effRTeta2.addHisto(h_auxrteta2, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
        effRTeta2.addHisto(DATA16RTeta2, 'PE,SAME', 'RT data 2016' , 'PL', r.kBlack , 1, 0)
        effRTeta2.addHisto(DATA17RTeta2, 'PE,SAME', 'RT data 2017' , 'PL', r.kBlue , 1, 0)
        effRTeta2.addHisto(MCRTeta2, 'PE,SAME', 'RT MC', 'PL', r.kGreen , 1, 0)
        effRTeta2.addBand(h_auxrteta2.GetXaxis().GetXmin(), da16rt-da16systrt, h_auxrteta2.GetXaxis().GetXmax(), da16rt+da16systrt, r.kOrange+6, 0.2)
        effRTeta2.addBand(h_auxrteta2.GetXaxis().GetXmin(), da17rt-da17systrt, h_auxrteta2.GetXaxis().GetXmax(), da17rt+da17systrt, r.kPink+6, 0.2)
        effRTeta2.addLatex (0.6, 0.3, 'Mean R_{T} data 2016: %.3f '%(da16rt))
        effRTeta2.addLatex (0.6, 0.25, 'Mean R_{T} data 2017: %.3f '%(da17rt))
        effRTeta2.addLatex (0.6, 0.2, 'Mean R_{T} MC  : %.3f '%(mcrt))
        effRTeta2.save(1, 1, 0, lumi, 0.8, 1.2)
   
        effMll = Canvas.Canvas('rt/%s/plot_eff_mll_'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
        h_auxMll = r.TH1F("h_auxMll", "", 1, 0, 250)
        h_auxMll.GetYaxis().SetRangeUser(0, 2)
        h_auxMll.GetXaxis().SetRangeUser(0, 250)
        h_auxMll.GetXaxis().SetTitle(labelx)
        effMll.addHisto(h_auxMll, 'h', '', '', r.kRed+1, 1, 0)
        DATA16effMllee.SetMarkerStyle(kOpenCircle)
        DATA16effMllmm.SetMarkerStyle(kOpenCircle)
        DATA16effMllSF.SetMarkerStyle(kOpenCircle)
        DATA16effMllOF.SetMarkerStyle(kOpenCircle)
        DATA17effMllee.SetMarkerStyle(kFullTriangleUp)
        DATA17effMllmm.SetMarkerStyle(kFullTriangleUp)
        DATA17effMllSF.SetMarkerStyle(kFullTriangleUp)
        DATA17effMllOF.SetMarkerStyle(kFullTriangleUp)
        effMll.addHisto(DATA16effMllee, 'PE,SAME', '2016 Double Electron', 'PL', r.kRed+1 , 1, 0)
        effMll.addHisto(DATA16effMllmm, 'PE,SAME', '2016 Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effMll.addHisto(DATA16effMllSF, 'PE,SAME', '2016 Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effMll.addHisto(DATA16effMllOF, 'PE,SAME', '2016 Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effMll.addHisto(DATA17effMllee, 'PE,SAME', '2017 Double Electron', 'PL', r.kRed+1 , 1, 0)
        effMll.addHisto(DATA17effMllmm, 'PE,SAME', '2017 Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effMll.addHisto(DATA17effMllSF, 'PE,SAME', '2017 Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effMll.addHisto(DATA17effMllOF, 'PE,SAME', '2017 Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        
        
        
        effMll.save(1, 1, 0, lumi, 0.2, 1.8)                                                                    

        effmt2 = Canvas.Canvas('rt/%s/plot_eff_mt2_'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
        h_auxmt2 = r.TH1F("h_auxmt2", "", 1, 0, 160)
        h_auxmt2.GetYaxis().SetRangeUser(0, 2)
        h_auxmt2.GetXaxis().SetRangeUser(0, 160)
        h_auxmt2.GetXaxis().SetTitle(labelx)
        effmt2.addHisto(h_auxmt2, 'h', '', '', r.kRed+1, 1, 0)
        DATA16effmt2ee.SetMarkerStyle(kOpenCircle)
        DATA16effmt2mm.SetMarkerStyle(kOpenCircle)
        DATA16effmt2SF.SetMarkerStyle(kOpenCircle)
        DATA16effmt2OF.SetMarkerStyle(kOpenCircle)
        DATA17effmt2ee.SetMarkerStyle(kFullTriangleUp)
        DATA17effmt2mm.SetMarkerStyle(kFullTriangleUp)
        DATA17effmt2SF.SetMarkerStyle(kFullTriangleUp)
        DATA17effmt2OF.SetMarkerStyle(kFullTriangleUp)
        effmt2.addHisto(DATA16effmt2ee, 'PE,SAME', '2016 Double Electron', 'PL', r.kRed+1 , 1, 0)
        effmt2.addHisto(DATA16effmt2mm, 'PE,SAME', '2016 Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effmt2.addHisto(DATA16effmt2SF, 'PE,SAME', '2016 Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effmt2.addHisto(DATA16effmt2OF, 'PE,SAME', '2016 Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effmt2.addHisto(DATA17effmt2ee, 'PE,SAME', '2017 Double Electron', 'PL', r.kRed+1 , 1, 0)
        effmt2.addHisto(DATA17effmt2mm, 'PE,SAME', '2017 Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effmt2.addHisto(DATA17effmt2SF, 'PE,SAME', '2017 Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effmt2.addHisto(DATA17effmt2OF, 'PE,SAME', '2017 Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effmt2.save(1, 1, 0, lumi, 0.2, 1.8)                                                                    

        effpt1 = Canvas.Canvas('rt/%s/plot_eff_pt1_'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
        h_auxpt1 = r.TH1F("h_auxpt1", "", 1, 0, 150)
        h_auxpt1.GetYaxis().SetRangeUser(0, 2)
        h_auxpt1.GetXaxis().SetRangeUser(0, 150)
        h_auxpt1.GetXaxis().SetTitle(labelpt1)
        effpt1.addHisto(h_auxpt1, 'h', '', '', r.kRed+1, 1, 0)
        DATA16effpt1ee.SetMarkerStyle(kOpenCircle)
        DATA16effpt1mm.SetMarkerStyle(kOpenCircle)
        DATA16effpt1SF.SetMarkerStyle(kOpenCircle)
        DATA16effpt1OF.SetMarkerStyle(kOpenCircle)
        DATA17effpt1ee.SetMarkerStyle(kFullTriangleUp)
        DATA17effpt1mm.SetMarkerStyle(kFullTriangleUp)
        DATA17effpt1SF.SetMarkerStyle(kFullTriangleUp)
        DATA17effpt1OF.SetMarkerStyle(kFullTriangleUp)
        effpt1.addHisto(DATA16effpt1ee, 'PE,SAME', '2016 Double Electron', 'PL', r.kRed+1 , 1, 0)
        effpt1.addHisto(DATA16effpt1mm, 'PE,SAME', '2016 Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effpt1.addHisto(DATA16effpt1SF, 'PE,SAME', '2016 Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effpt1.addHisto(DATA16effpt1OF, 'PE,SAME', '2016 Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effpt1.addHisto(DATA17effpt1ee, 'PE,SAME', '2017 Double Electron', 'PL', r.kRed+1 , 1, 0)
        effpt1.addHisto(DATA17effpt1mm, 'PE,SAME', '2017 Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effpt1.addHisto(DATA17effpt1SF, 'PE,SAME', '2017 Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effpt1.addHisto(DATA17effpt1OF, 'PE,SAME', '2017 Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effpt1.save(1, 1, 0, lumi, 0.2, 1.8)

        effpt2 = Canvas.Canvas('rt/%s/plot_eff_pt2_'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
        h_auxpt2 = r.TH1F("h_auxpt2", "", 1, 0, 150)
        h_auxpt2.GetYaxis().SetRangeUser(0, 2)
        h_auxpt2.GetXaxis().SetRangeUser(0, 150)
        h_auxpt2.GetXaxis().SetTitle(labelpt2)
        DATA16effpt2ee.SetMarkerStyle(kOpenCircle)
        DATA16effpt2mm.SetMarkerStyle(kOpenCircle)
        DATA16effpt2SF.SetMarkerStyle(kOpenCircle)
        DATA16effpt2OF.SetMarkerStyle(kOpenCircle)
        DATA17effpt2ee.SetMarkerStyle(kFullTriangleUp)
        DATA17effpt2mm.SetMarkerStyle(kFullTriangleUp)
        DATA17effpt2SF.SetMarkerStyle(kFullTriangleUp)
        DATA17effpt2OF.SetMarkerStyle(kFullTriangleUp)
        effpt2.addHisto(h_auxpt2, 'h', '', '', r.kRed+1, 1, 0)
        effpt2.addHisto(DATA16effpt2ee, 'PE,SAME', '2016 Double Electron', 'PL', r.kRed+1 , 1, 0)
        effpt2.addHisto(DATA16effpt2mm, 'PE,SAME', '2016 Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effpt2.addHisto(DATA16effpt2SF, 'PE,SAME', '2016 Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effpt2.addHisto(DATA16effpt2OF, 'PE,SAME', '2016 Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effpt2.addHisto(DATA17effpt2ee, 'PE,SAME', '2017 Double Electron', 'PL', r.kRed+1 , 1, 0)
        effpt2.addHisto(DATA17effpt2mm, 'PE,SAME', '2017 Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effpt2.addHisto(DATA17effpt2SF, 'PE,SAME', '2017 Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effpt2.addHisto(DATA17effpt2OF, 'PE,SAME', '2017 Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effpt2.save(1, 1, 0, lumi, 0.2, 1.8)
 
        effeta1 = Canvas.Canvas('rt/%s/plot_eff_eta1_'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
        h_auxeta1 = r.TH1F("h_auxeta1", "", 1, -2.4, 2.4)
        h_auxeta1.GetYaxis().SetRangeUser(0, 2)
        h_auxeta1.GetXaxis().SetRangeUser(-2.4, 2.4)
        h_auxeta1.GetXaxis().SetTitle(labeleta1)
        effeta1.addHisto(h_auxeta1, 'h', '', '', r.kRed+1, 1, 0)
        DATA16effeta1ee.SetMarkerStyle(kOpenCircle)
        DATA16effeta1mm.SetMarkerStyle(kOpenCircle)
        DATA16effeta1SF.SetMarkerStyle(kOpenCircle)
        DATA16effeta1OF.SetMarkerStyle(kOpenCircle)
        DATA17effeta1ee.SetMarkerStyle(kFullTriangleUp)
        DATA17effeta1mm.SetMarkerStyle(kFullTriangleUp)
        DATA17effeta1SF.SetMarkerStyle(kFullTriangleUp)
        DATA17effeta1OF.SetMarkerStyle(kFullTriangleUp)
        effeta1.addHisto(DATA16effeta1ee, 'PE,SAME', '2016 Double Electron', 'PL', r.kRed+1 , 1, 0)
        effeta1.addHisto(DATA16effeta1mm, 'PE,SAME', '2016 Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effeta1.addHisto(DATA16effeta1SF, 'PE,SAME', '2016 Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effeta1.addHisto(DATA16effeta1OF, 'PE,SAME', '2016 Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effeta1.addHisto(DATA17effeta1ee, 'PE,SAME', '2017 Double Electron', 'PL', r.kRed+1 , 1, 0)
        effeta1.addHisto(DATA17effeta1mm, 'PE,SAME', '2017 Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effeta1.addHisto(DATA17effeta1SF, 'PE,SAME', '2017 Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effeta1.addHisto(DATA17effeta1OF, 'PE,SAME', '2017 Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effeta1.save(1, 1, 0, lumi, 0.2, 1.8)

        effeta2 = Canvas.Canvas('rt/%s/plot_eff_eta2_'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
        h_auxeta2 = r.TH1F("h_auxeta2", "", 1, -2.4, 2.4)
        h_auxeta2.GetYaxis().SetRangeUser(0, 2)
        h_auxeta2.GetXaxis().SetRangeUser(-2.4, 2.4)
        h_auxeta2.GetXaxis().SetTitle(labeleta2)
        effeta2.addHisto(h_auxeta2, 'h', '', '', r.kRed+1, 1, 0)
        DATA16effeta2ee.SetMarkerStyle(kOpenCircle)
        DATA16effeta2mm.SetMarkerStyle(kOpenCircle)
        DATA16effeta2SF.SetMarkerStyle(kOpenCircle)
        DATA16effeta2OF.SetMarkerStyle(kOpenCircle)
        DATA17effeta2ee.SetMarkerStyle(kFullTriangleUp)
        DATA17effeta2mm.SetMarkerStyle(kFullTriangleUp)
        DATA17effeta2SF.SetMarkerStyle(kFullTriangleUp)
        DATA17effeta2OF.SetMarkerStyle(kFullTriangleUp)
        effeta2.addHisto(DATA16effeta2ee, 'PE,SAME', '2016 Double Electron', 'PL', r.kRed+1 , 1, 0)
        effeta2.addHisto(DATA16effeta2mm, 'PE,SAME', '2016 Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effeta2.addHisto(DATA16effeta2SF, 'PE,SAME', '2016 Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effeta2.addHisto(DATA16effeta2OF, 'PE,SAME', '2016 Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effeta2.addHisto(DATA17effeta2ee, 'PE,SAME', '2017 Double Electron', 'PL', r.kRed+1 , 1, 0)
        effeta2.addHisto(DATA17effeta2mm, 'PE,SAME', '2017 Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effeta2.addHisto(DATA17effeta2SF, 'PE,SAME', '2017 Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effeta2.addHisto(DATA17effeta2OF, 'PE,SAME', '2017 Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effeta2.save(1, 1, 0, lumi, 0.2, 1.8)
 
        print 'Measured RT value data 2016 ', da16rt, ' +/- ', da16uncrt, ' +/- ', da16systrt
        print 'Measured RT value data 2017 ', da17rt, ' +/- ', da17uncrt, ' +/- ', da17systrt
        print 'Measured RT value MC   ', mcrt, ' +/- ', mcuncrt, ' +/- ', mcsystrt


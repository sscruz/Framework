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
    if(eff_em == 0 ):
        print "Division by zero"
        return [0, 0, 0]
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
    mcDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO', 'TTJets_DiLepton', 'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT',  'T_tch_powheg', 'TBar_tch_powheg', 'WWTo2L2Nu',  'WZTo3LNu','WZTo2L2Q', 'ZZTo4L', 'ZZTo2L2Nu', 'ZZTo2L2Q', 'WWW', 'WWZ', 'WZZ', 'ZZZ', 'TWZ', 'tZq_ll', 'TTZToLLNuNu' , 'TTZToQQ', 'TTLLJets_m1to10', 'TTWToLNu','TTWToQQ',  'TTTT', 'TTHnobb_pow', 'VHToNonbb',  'GGHZZ4L',  'WJetsToLNu_LO']
    
    daDatasetsB = ['DoubleEG_Run2016B_03Feb2017_ver2_v2_runs_273150_275376',
                   'DoubleMuon_Run2016B_03Feb2017_ver2_v2_runs_273150_275376', 
                   'MuonEG_Run2016B_03Feb2017_ver2_v2_runs_273150_275376', 
                   'SingleMuon_Run2016B_03Feb2017_ver2_v2_runs_273150_275376',
                   'SingleElectron_Run2016B_03Feb2017_ver2_v2_runs_273150_275376',
                   'MET_Run2016B_03Feb2017_ver2_v2_runs_273150_275376', 
                   'JetHT_Run2016B_03Feb2017_ver2_v2_runs_273150_275376']             
 
    daDatasetsC = ['DoubleEG_Run2016C_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016C_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016C_03Feb2017_v1_runs_271036_284044', 
                   'SingleMuon_Run2016C_03Feb2017_v1_runs_271036_284044',
                   'SingleElectron_Run2016C_03Feb2017_v1_runs_271036_284044',
                   'MET_Run2016C_03Feb2017_v1_runs_271036_284044', 
                   'JetHT_Run2016C_03Feb2017_v1_runs_271036_284044']                    
    
    daDatasetsD = ['DoubleEG_Run2016D_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016D_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016D_03Feb2017_v1_runs_271036_284044',
                   'SingleMuon_Run2016D_03Feb2017_v1_runs_271036_284044',
                   'SingleElectron_Run2016D_03Feb2017_v1_runs_271036_284044',
                   'MET_Run2016D_03Feb2017_v1_runs_271036_284044', 
                   'JetHT_Run2016D_03Feb2017_v1_runs_271036_284044']          
 
    daDatasetsE = ['DoubleEG_Run2016E_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016E_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016E_03Feb2017_v1_runs_271036_284044',
                   'SingleMuon_Run2016E_03Feb2017_v1_runs_271036_284044',
                   'SingleElectron_Run2016E_03Feb2017_v1_runs_271036_284044',
                   'MET_Run2016E_03Feb2017_v1_runs_271036_284044', 
                   'JetHT_Run2016E_03Feb2017_v1_runs_271036_284044']          
 
    daDatasetsF = ['DoubleEG_Run2016F_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016F_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016F_03Feb2017_v1_runs_271036_284044', 
                   'SingleMuon_Run2016F_03Feb2017_v1_runs_271036_284044',        
                   'SingleElectron_Run2016F_03Feb2017_v1_runs_271036_284044',
                   'MET_Run2016F_03Feb2017_v1_runs_271036_284044', 
                   'JetHT_Run2016F_03Feb2017_v1_runs_271036_284044']          

    daDatasetsG = ['DoubleEG_Run2016G_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016G_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016G_03Feb2017_v1_runs_271036_284044',
                   'SingleMuon_Run2016G_03Feb2017_v1_runs_271036_284044',       
                   'SingleElectron_Run2016G_03Feb2017_v1_runs_271036_284044',
                   'MET_Run2016G_03Feb2017_v1_runs_271036_284044', 
                   'JetHT_Run2016G_03Feb2017_v1_runs_271036_284044']          

    daDatasetsH = ['DoubleEG_Run2016H_03Feb2017_ver2_v1_runs_281085_284035',
                   'DoubleEG_Run2016H_03Feb2017_ver3_v1_runs_284036_284044',
                   'DoubleMuon_Run2016H_03Feb2017_ver2_v1_runs_281085_284035',
                   'DoubleMuon_Run2016H_03Feb2017_ver3_v1_runs_284036_284044',
                   'MuonEG_Run2016H_03Feb2017_ver2_v1_runs_281085_284035', 
                   'MuonEG_Run2016H_03Feb2017_ver3_v1_runs_284036_284044', 
                   'SingleElectron_Run2016H_03Feb2017_ver2_v1_runs_281085_284035',
                   'SingleElectron_Run2016H_03Feb2017_ver3_v1_runs_284036_284044',
                   'SingleMuon_Run2016H_03Feb2017_ver2_v1_runs_281085_284035',
                   'SingleMuon_Run2016H_03Feb2017_ver3_v1_runs_284036_284044',
                   'MET_Run2016H_03Feb2017_ver2_v1_runs_281085_284035', 
                   'MET_Run2016H_03Feb2017_ver3_v1_runs_284036_284044',
                   'JetHT_Run2016H_03Feb2017_ver2_v1_runs_281085_284035',
                   'JetHT_Run2016H_03Feb2017_ver3_v1_runs_284036_284044']    
 
    daDatasets = daDatasetsB + daDatasetsC + daDatasetsD +daDatasetsE + daDatasetsF + daDatasetsG + daDatasetsH     
    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print bcolors.HEADER + '[RTAnalysis] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    
    lumi = 35.9 ; maxrun = 999999; lumi_str = '35.9invfb'
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
    DATAdenMlleevalue =   treeDA.getTH1F(lumi, "DATAdeneevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den,  cuts.ee]), '', labelx,"1", kf)
    DATAnumMlleevalue =     treeDA.getTH1F(lumi, "DATAnumeevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.ee]), '', labelx,"1", kf)
    DATAdenMllmmvalue =   treeDA.getTH1F(lumi, "DATAdenmmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den, cuts.mm]), '', labelx,"1", kf)
    DATAnumMllmmvalue =     treeDA.getTH1F(lumi, "DATAnummmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.mm]), '', labelx,"1", kf)
    DATAdenMllSFvalue =   treeDA.getTH1F(lumi, "DATAdenSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den, cuts.SF]), '', labelx,"1", kf)
    DATAnumMllSFvalue =     treeDA.getTH1F(lumi, "DATAnumSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.SF]), '', labelx,"1", kf)
    DATAdenMllOFvalue =   treeDA.getTH1F(lumi, "DATAdenOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den, cuts.OF]), '', labelx,"1", kf)
    DATAnumMllOFvalue =     treeDA.getTH1F(lumi, "DATAnumOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.OF]), '', labelx,"1", kf)
    
                                                                                                                                                                                                 
    MCdenMlleevalue =   treeMC.getTH1F(lumi, "MCdeneevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den, cuts.ee]), '', labelx,"1", kf)
    MCnumMlleevalue =     treeMC.getTH1F(lumi, "MCnumeevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.ee]), '', labelx,"1", kf)
    MCdenMllmmvalue =   treeMC.getTH1F(lumi, "MCdenmmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den, cuts.mm]), '', labelx,"1", kf)
    MCnumMllmmvalue =     treeMC.getTH1F(lumi, "MCnummmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.mm]), '', labelx,"1", kf)
    MCdenMllSFvalue =   treeMC.getTH1F(lumi, "MCdenSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den, cuts.SF]), '', labelx,"1", kf)
    MCnumMllSFvalue =     treeMC.getTH1F(lumi, "MCnumSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.SF]), '', labelx,"1", kf)
    MCdenMllOFvalue =   treeMC.getTH1F(lumi, "MCdenOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.den, cuts.OF]), '', labelx,"1", kf)
    MCnumMllOFvalue =     treeMC.getTH1F(lumi, "MCnumOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.OF]), '', labelx,"1", kf)
    
    effMlleevalueDATA = getTriggerEffs(DATAnumMlleevalue, DATAdenMlleevalue, "data")
    print "effMlleevalueDATA ", effMlleevalueDATA.GetY()[0]
    effMllmmvalueDATA = getTriggerEffs(DATAnumMllmmvalue, DATAdenMllmmvalue, "data")
    print "effMllmmvalueDATA ", effMllmmvalueDATA.GetY()[0]
    effMllSFvalueDATA = getTriggerEffs(DATAnumMllSFvalue, DATAdenMllSFvalue, "data")
    print "effMllSFvalueDATA ", effMllSFvalueDATA.GetY()[0]
    effMllOFvalueDATA = getTriggerEffs(DATAnumMllOFvalue, DATAdenMllOFvalue, "data")
    print "effMllOFvalueDATA ", effMllOFvalueDATA.GetY()[0]
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

    eevalda_ = effMlleevalueDATA.GetY()
    eevaldaeh_ = effMlleevalueDATA.GetEYhigh()
    eevaldael_ = effMlleevalueDATA.GetEYlow()
    eevalda = eevalda_[0]
    eevaldae = max(eevaldaeh_[0], eevaldael_[0])
    mmvalda_ = effMllmmvalueDATA.GetY()
    mmvaldaeh_ = effMllmmvalueDATA.GetEYhigh()
    mmvaldael_ = effMllmmvalueDATA.GetEYlow()
    mmvalda = mmvalda_[0]
    mmvaldae = max(mmvaldaeh_[0], mmvaldael_[0])
    emvalda_ = effMllOFvalueDATA.GetY()
    emvaldaeh_ = effMllOFvalueDATA.GetEYhigh()
    emvaldael_ = effMllOFvalueDATA.GetEYlow()
    emvalda = emvalda_[0]
    emvaldae = max(emvaldaeh_[0], emvaldael_[0])

    eevalmc_ = effMlleevalueMC.GetY()
    print "eevalmc_ ", eevalmc_
    eevalmceh_ = effMlleevalueMC.GetEYhigh()
    eevalmcel_ = effMlleevalueMC.GetEYlow()
    eevalmc = eevalmc_[0]
    eevalmce = max(eevalmceh_[0], eevalmcel_[0])
    print "eevalmc ", eevalmc
    mmvalmc_ = effMllmmvalueMC.GetY()
    print "mmvalmc_ ", mmvalmc_
    mmvalmceh_ = effMllmmvalueMC.GetEYhigh()
    mmvalmcel_ = effMllmmvalueMC.GetEYlow()
    mmvalmc = mmvalmc_[0]
    mmvalmce = max(mmvalmceh_[0], mmvalmcel_[0])
    print "mmvalmc ", mmvalmc
    emvalmc_ = effMllOFvalueMC.GetY()
    emvalmceh_ = effMllOFvalueMC.GetEYhigh()
    emvalmcel_ = effMllOFvalueMC.GetEYlow()
    emvalmc = emvalmc_[0]
    emvalmce = max(emvalmceh_[0], emvalmcel_[0])

    [dart, dauncrt, dasystrt] = getRT(eevalda, eevaldae, mmvalda, mmvaldae, emvalda, emvaldae)
    [mcrt, mcuncrt, mcsystrt] = getRT(eevalmc, eevalmce, mmvalmc, mmvalmce, emvalmc, emvalmce)
    print 'Measured RT value data ', dart, ' +/- ', dauncrt, ' +/- ', dasystrt
    print 'Measured RT value MC   ', mcrt, ' +/- ', mcuncrt, ' +/- ', mcsystrt
    saveInFile(theFile, mcrt, mcuncrt, mcsystrt, dart, dauncrt, dasystrt)
    makeTable(DATAnumMlleevalue, DATAnumMllmmvalue, DATAnumMllOFvalue, DATAdenMlleevalue, DATAdenMllmmvalue, DATAdenMllOFvalue, eevalda, mmvalda, emvalda, eevaldae, mmvaldae, emvaldae, dart, dauncrt, dasystrt, "DATA")
    makeTable(MCnumMlleevalue, MCnumMllmmvalue, MCnumMllOFvalue, MCdenMlleevalue, MCdenMllmmvalue, MCdenMllOFvalue, eevalmc, mmvalmc, emvalmc, eevalmce, mmvalmce, emvalmce, mcrt, mcuncrt, mcsystrt, "MC")
    ######################## Calculation of total values #############################

    if makeDependencyPlots: 
        DATAdenMllee = treeDA.getTH1F(lumi, "DATAdenee", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.ee]), '', labelx,"1", kf)
        DATAnumMllee =   treeDA.getTH1F(lumi, "DATAnumee", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.ee]), '', labelx,"1", kf)
        DATAdenMllmm = treeDA.getTH1F(lumi, "DATAdenmm", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.mm]), '', labelx,"1", kf)
        DATAnumMllmm =   treeDA.getTH1F(lumi, "DATAnummm", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.mm]), '', labelx,"1", kf)
        DATAdenMllSF = treeDA.getTH1F(lumi, "DATAdenSF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.SF]), '', labelx,"1", kf)
        DATAnumMllSF =   treeDA.getTH1F(lumi, "DATAnumSF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.SF]), '', labelx,"1", kf)
        DATAdenMllOF = treeDA.getTH1F(lumi, "DATAdenOF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.OF]), '', labelx,"1", kf)
        DATAnumMllOF =   treeDA.getTH1F(lumi, "DATAnumOF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.OF]), '', labelx,"1", kf)
        
        DATAdenpt1ee = treeDA.getTH1F(lumi, "DATAdenpt1ee", "Lep1_pt_Edge", ptbins , 1, 1, cuts.AddList([specialcut, cuts.den, cuts.ee]), '', labelpt1,"1", kf)
        DATAnumpt1ee =   treeDA.getTH1F(lumi, "DATAnumpt1ee", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.ee]), '', labelpt1,"1", kf)
        DATAdenpt1mm = treeDA.getTH1F(lumi, "DATAdenpt1mm", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.mm]), '', labelpt1,"1", kf)
        DATAnumpt1mm =   treeDA.getTH1F(lumi, "DATAnumpt1mm", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.mm]), '', labelpt1,"1", kf)
        DATAdenpt1SF = treeDA.getTH1F(lumi, "DATAdenpt1SF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.SF]), '', labelpt1,"1", kf)
        DATAnumpt1SF =   treeDA.getTH1F(lumi, "DATAnumpt1SF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.SF]), '', labelpt1,"1", kf)
        DATAdenpt1OF = treeDA.getTH1F(lumi, "DATAdenpt1OF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.OF]), '', labelpt1,"1", kf)
        DATAnumpt1OF =   treeDA.getTH1F(lumi, "DATAnumpt1OF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.OF]), '', labelpt1,"1", kf)
                                                                                                                                                                                        
        DATAdenpt2ee = treeDA.getTH1F(lumi, "DATAdenpt2ee", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.ee]), '', labelpt2,"1", kf)
        DATAnumpt2ee =   treeDA.getTH1F(lumi, "DATAnumpt2ee", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.ee]), '', labelpt2,"1", kf)
        DATAdenpt2mm = treeDA.getTH1F(lumi, "DATAdenpt2mm", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.mm]), '', labelpt2,"1", kf)
        DATAnumpt2mm =   treeDA.getTH1F(lumi, "DATAnumpt2mm", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.mm]), '', labelpt2,"1", kf)
        DATAdenpt2SF = treeDA.getTH1F(lumi, "DATAdenpt2SF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.SF]), '', labelpt2,"1", kf)
        DATAnumpt2SF =   treeDA.getTH1F(lumi, "DATAnumpt2SF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.SF]), '', labelpt2,"1", kf)
        DATAdenpt2OF = treeDA.getTH1F(lumi, "DATAdenpt2OF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.OF]), '', labelpt2,"1", kf)
        DATAnumpt2OF =   treeDA.getTH1F(lumi, "DATAnumpt2OF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.OF]), '', labelpt2,"1", kf)
                                                                                                                                                                                        
        DATAdeneta1ee = treeDA.getTH1F(lumi, "DATAdeneta1ee", "Lep1_eta_Edge", etabins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.ee]), '', labeleta1,"1", kf)
        DATAnumeta1ee =   treeDA.getTH1F(lumi, "DATAnumeta1ee", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num,  cuts.ee]), '', labeleta1,"1", kf)
        DATAdeneta1mm = treeDA.getTH1F(lumi, "DATAdeneta1mm", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den, cuts.mm]), '', labeleta1,"1", kf)
        DATAnumeta1mm =   treeDA.getTH1F(lumi, "DATAnumeta1mm", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num,  cuts.mm]), '', labeleta1,"1", kf)
        DATAdeneta1SF = treeDA.getTH1F(lumi, "DATAdeneta1SF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den, cuts.SF]), '', labeleta1,"1", kf)
        DATAnumeta1SF =   treeDA.getTH1F(lumi, "DATAnumeta1SF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num,  cuts.SF]), '', labeleta1,"1", kf)
        DATAdeneta1OF = treeDA.getTH1F(lumi, "DATAdeneta1OF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den, cuts.OF]), '', labeleta1,"1", kf)
        DATAnumeta1OF =   treeDA.getTH1F(lumi, "DATAnumeta1OF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num,  cuts.OF]), '', labeleta1,"1", kf)
                                                                                                                                                                                        
        DATAdeneta2ee = treeDA.getTH1F(lumi, "DATAdeneta2ee", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den, cuts.ee]), '', labeleta2,"1", kf)
        DATAnumeta2ee =   treeDA.getTH1F(lumi, "DATAnumeta2ee", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num,  cuts.ee]), '', labeleta2,"1", kf)
        DATAdeneta2mm = treeDA.getTH1F(lumi, "DATAdeneta2mm", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den, cuts.mm]), '', labeleta2,"1", kf)
        DATAnumeta2mm =   treeDA.getTH1F(lumi, "DATAnumeta2mm", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num,  cuts.mm]), '', labeleta2,"1", kf)
        DATAdeneta2SF = treeDA.getTH1F(lumi, "DATAdeneta2SF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den, cuts.SF]), '', labeleta2,"1", kf)
        DATAnumeta2SF =   treeDA.getTH1F(lumi, "DATAnumeta2SF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num,  cuts.SF]), '', labeleta2,"1", kf)
        DATAdeneta2OF = treeDA.getTH1F(lumi, "DATAdeneta2OF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den, cuts.OF]), '', labeleta2,"1", kf)
        DATAnumeta2OF =   treeDA.getTH1F(lumi, "DATAnumeta2OF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num,  cuts.OF]), '', labeleta2,"1", kf)
                                                                                                                                                                                        
        DATAdenMETee =   treeDA.getTH1F(lumi, "DATAdeneevalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.ee]), '', labelmet,"1", kf)
        DATAnumMETee =     treeDA.getTH1F(lumi, "DATAnumeevalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.ee]), '', labelmet,"1", kf)
        DATAdenMETmm =   treeDA.getTH1F(lumi, "DATAdenmmvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.mm]), '', labelmet,"1", kf)
        DATAnumMETmm =     treeDA.getTH1F(lumi, "DATAnummmvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.mm]), '', labelmet,"1", kf)
        DATAdenMETSF =   treeDA.getTH1F(lumi, "DATAdenSFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.SF]), '', labelmet,"1", kf)
        DATAnumMETSF =     treeDA.getTH1F(lumi, "DATAnumSFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.SF]), '', labelmet,"1", kf)
        DATAdenMETOF =   treeDA.getTH1F(lumi, "DATAdenOFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.OF]), '', labelmet,"1", kf)
        DATAnumMETOF =     treeDA.getTH1F(lumi, "DATAnumOFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.OF]), '', labelmet,"1", kf)
                                                                                                                                                                                        
        DATAdenmt2ee =   treeDA.getTH1F(lumi, "DATAdeneevalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.ee]), '', labelmt2,"1", kf)
        DATAnummt2ee =     treeDA.getTH1F(lumi, "DATAnumeevalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.ee]), '', labelmt2,"1", kf)
        DATAdenmt2mm =   treeDA.getTH1F(lumi, "DATAdenmmvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.mm]), '', labelmt2,"1", kf)
        DATAnummt2mm =     treeDA.getTH1F(lumi, "DATAnummmvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.mm]), '', labelmt2,"1", kf)
        DATAdenmt2SF =   treeDA.getTH1F(lumi, "DATAdenSFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.SF]), '', labelmt2,"1", kf)
        DATAnummt2SF =     treeDA.getTH1F(lumi, "DATAnumSFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.SF]), '', labelmt2,"1", kf)
        DATAdenmt2OF =   treeDA.getTH1F(lumi, "DATAdenOFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.OF]), '', labelmt2,"1", kf)
        DATAnummt2OF =     treeDA.getTH1F(lumi, "DATAnumOFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.OF]), '', labelmt2,"1", kf)            
                                                                                                                                                                                        
        MCdenMllee = treeMC.getTH1F(lumi, "MCdenee", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.ee]), '', labelx,"1", kf)
        MCnumMllee =   treeMC.getTH1F(lumi, "MCnumee", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.ee]), '', labelx,"1", kf)
        MCdenMllmm = treeMC.getTH1F(lumi, "MCdenmm", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.mm]), '', labelx,"1", kf)
        MCnumMllmm =   treeMC.getTH1F(lumi, "MCnummm", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.mm]), '', labelx,"1", kf)
        MCdenMllSF = treeMC.getTH1F(lumi, "MCdenSF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.SF]), '', labelx,"1", kf)
        MCnumMllSF =   treeMC.getTH1F(lumi, "MCnumSF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.SF]), '', labelx,"1", kf)
        MCdenMllOF = treeMC.getTH1F(lumi, "MCdenOF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.OF]), '', labelx,"1", kf)
        MCnumMllOF =   treeMC.getTH1F(lumi, "MCnumOF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.OF]), '', labelx,"1", kf)
                                                                                                                                                                                        
        MCdenpt1ee = treeMC.getTH1F(lumi, "MCdenpt1ee", "Lep1_pt_Edge", ptbins , 1, 1, cuts.AddList([specialcut, cuts.den, cuts.ee]), '', labelpt1,"1", kf)
        MCnumpt1ee =   treeMC.getTH1F(lumi, "MCnumpt1ee", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.ee]), '', labelpt1,"1", kf)
        MCdenpt1mm = treeMC.getTH1F(lumi, "MCdenpt1mm", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.mm]), '', labelpt1,"1", kf)
        MCnumpt1mm =   treeMC.getTH1F(lumi, "MCnumpt1mm", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.mm]), '', labelpt1,"1", kf)
        MCdenpt1SF = treeMC.getTH1F(lumi, "MCdenpt1SF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.SF]), '', labelpt1,"1", kf)
        MCnumpt1SF =   treeMC.getTH1F(lumi, "MCnumpt1SF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.SF]), '', labelpt1,"1", kf)
        MCdenpt1OF = treeMC.getTH1F(lumi, "MCdenpt1OF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.OF]), '', labelpt1,"1", kf)
        MCnumpt1OF =   treeMC.getTH1F(lumi, "MCnumpt1OF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.OF]), '', labelpt1,"1", kf)
        
        MCdenpt2ee = treeMC.getTH1F(lumi, "MCdenpt2ee", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.ee]), '', labelpt2,"1", kf)
        MCnumpt2ee =   treeMC.getTH1F(lumi, "MCnumpt2ee", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.ee]), '', labelpt2,"1", kf)
        MCdenpt2mm = treeMC.getTH1F(lumi, "MCdenpt2mm", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.mm]), '', labelpt2,"1", kf)
        MCnumpt2mm =   treeMC.getTH1F(lumi, "MCnumpt2mm", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.mm]), '', labelpt2,"1", kf)
        MCdenpt2SF = treeMC.getTH1F(lumi, "MCdenpt2SF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.SF]), '', labelpt2,"1", kf)
        MCnumpt2SF =   treeMC.getTH1F(lumi, "MCnumpt2SF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.SF]), '', labelpt2,"1", kf)
        MCdenpt2OF = treeMC.getTH1F(lumi, "MCdenpt2OF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.OF]), '', labelpt2,"1", kf)
        MCnumpt2OF =   treeMC.getTH1F(lumi, "MCnumpt2OF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.OF]), '', labelpt2,"1", kf)
                                                                                                                                                                                        
        MCdeneta1ee = treeMC.getTH1F(lumi, "MCdeneta1ee", "Lep1_eta_Edge", etabins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.ee]), '', labeleta1,"1", kf)
        MCnumeta1ee =   treeMC.getTH1F(lumi, "MCnumeta1ee", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num,  cuts.ee]), '', labeleta1,"1", kf)
        MCdeneta1mm = treeMC.getTH1F(lumi, "MCdeneta1mm", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den, cuts.mm]), '', labeleta1,"1", kf)
        MCnumeta1mm =   treeMC.getTH1F(lumi, "MCnumeta1mm", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num,  cuts.mm]), '', labeleta1,"1", kf)
        MCdeneta1SF = treeMC.getTH1F(lumi, "MCdeneta1SF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den, cuts.SF]), '', labeleta1,"1", kf)
        MCnumeta1SF =   treeMC.getTH1F(lumi, "MCnumeta1SF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num,  cuts.SF]), '', labeleta1,"1", kf)
        MCdeneta1OF = treeMC.getTH1F(lumi, "MCdeneta1OF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den, cuts.OF]), '', labeleta1,"1", kf)
        MCnumeta1OF =   treeMC.getTH1F(lumi, "MCnumeta1OF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num,  cuts.OF]), '', labeleta1,"1", kf)
        
        MCdeneta2ee = treeMC.getTH1F(lumi, "MCdeneta2ee", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den, cuts.ee]), '', labeleta2,"1", kf)
        MCnumeta2ee =   treeMC.getTH1F(lumi, "MCnumeta2ee", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num,  cuts.ee]), '', labeleta2,"1", kf)
        MCdeneta2mm = treeMC.getTH1F(lumi, "MCdeneta2mm", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den, cuts.mm]), '', labeleta2,"1", kf)
        MCnumeta2mm =   treeMC.getTH1F(lumi, "MCnumeta2mm", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num,  cuts.mm]), '', labeleta2,"1", kf)
        MCdeneta2SF = treeMC.getTH1F(lumi, "MCdeneta2SF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den, cuts.SF]), '', labeleta2,"1", kf)
        MCnumeta2SF =   treeMC.getTH1F(lumi, "MCnumeta2SF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num,  cuts.SF]), '', labeleta2,"1", kf)
        MCdeneta2OF = treeMC.getTH1F(lumi, "MCdeneta2OF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.den, cuts.OF]), '', labeleta2,"1", kf)
        MCnumeta2OF =   treeMC.getTH1F(lumi, "MCnumeta2OF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.num,  cuts.OF]), '', labeleta2,"1", kf)
                                                                                                                                                                                        
        MCdenMETee =   treeMC.getTH1F(lumi, "MCdeneevalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.ee]), '', labelmet,"1", kf)
        MCnumMETee =     treeMC.getTH1F(lumi, "MCnumeevalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.ee]), '', labelmet,"1", kf)
        MCdenMETmm =   treeMC.getTH1F(lumi, "MCdenmmvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.mm]), '', labelmet,"1", kf)
        MCnumMETmm =     treeMC.getTH1F(lumi, "MCnummmvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.mm]), '', labelmet,"1", kf)
        MCdenMETSF =   treeMC.getTH1F(lumi, "MCdenSFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.SF]), '', labelmet,"1", kf)
        MCnumMETSF =     treeMC.getTH1F(lumi, "MCnumSFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.SF]), '', labelmet,"1", kf)
        MCdenMETOF =   treeMC.getTH1F(lumi, "MCdenOFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.OF]), '', labelmet,"1", kf)
        MCnumMETOF =     treeMC.getTH1F(lumi, "MCnumOFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.OF]), '', labelmet,"1", kf)
        
        MCdenmt2ee =   treeMC.getTH1F(lumi, "MCdeneevalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.ee]), '', labelmt2,"1", kf)
        MCnummt2ee =     treeMC.getTH1F(lumi, "MCnumeevalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.ee]), '', labelmt2,"1", kf)
        MCdenmt2mm =   treeMC.getTH1F(lumi, "MCdenmmvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.mm]), '', labelmt2,"1", kf)
        MCnummt2mm =     treeMC.getTH1F(lumi, "MCnummmvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.mm]), '', labelmt2,"1", kf)
        MCdenmt2SF =   treeMC.getTH1F(lumi, "MCdenSFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.SF]), '', labelmt2,"1", kf)
        MCnummt2SF =     treeMC.getTH1F(lumi, "MCnumSFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.SF]), '', labelmt2,"1", kf)
        MCdenmt2OF =   treeMC.getTH1F(lumi, "MCdenOFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.den, cuts.OF]), '', labelmt2,"1", kf)
        MCnummt2OF =     treeMC.getTH1F(lumi, "MCnumOFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.num,  cuts.OF]), '', labelmt2,"1", kf )  
        
        
        DATAeffMETee =  getTriggerEffs(DATAnumMETee, DATAdenMETee, "data")
        DATAeffMETmm =  getTriggerEffs(DATAnumMETmm, DATAdenMETmm, "data")    
        DATAeffMETSF =  getTriggerEffs(DATAnumMETSF, DATAdenMETSF, "data")    
        DATAeffMETOF =  getTriggerEffs(DATAnumMETOF, DATAdenMETOF, "data")
        DATAeffmt2ee =  getTriggerEffs(DATAnummt2ee, DATAdenmt2ee, "data")
        DATAeffmt2mm =  getTriggerEffs(DATAnummt2mm, DATAdenmt2mm, "data")
        DATAeffmt2SF =  getTriggerEffs(DATAnummt2SF, DATAdenmt2SF, "data")
        DATAeffmt2OF =  getTriggerEffs(DATAnummt2OF, DATAdenmt2OF, "data")
        DATAeffMllee =  getTriggerEffs(DATAnumMllee, DATAdenMllee, "data")
        DATAeffMllmm =  getTriggerEffs(DATAnumMllmm, DATAdenMllmm, "data")
        DATAeffMllSF =  getTriggerEffs(DATAnumMllSF, DATAdenMllSF, "data")
        DATAeffMllOF =  getTriggerEffs(DATAnumMllOF, DATAdenMllOF, "data")
        DATAeffpt1ee =  getTriggerEffs(DATAnumpt1ee, DATAdenpt1ee, "data")
        DATAeffpt1mm =  getTriggerEffs(DATAnumpt1mm, DATAdenpt1mm, "data")
        DATAeffpt1SF =  getTriggerEffs(DATAnumpt1SF, DATAdenpt1SF, "data")
        DATAeffpt1OF =  getTriggerEffs(DATAnumpt1OF, DATAdenpt1OF, "data")
        DATAeffpt2ee =  getTriggerEffs(DATAnumpt2ee, DATAdenpt2ee, "data")
        DATAeffpt2mm =  getTriggerEffs(DATAnumpt2mm, DATAdenpt2mm, "data")
        DATAeffpt2SF =  getTriggerEffs(DATAnumpt2SF, DATAdenpt2SF, "data")
        DATAeffpt2OF =  getTriggerEffs(DATAnumpt2OF, DATAdenpt2OF, "data")
        DATAeffeta1ee = getTriggerEffs(DATAnumeta1ee, DATAdeneta1ee, "data")
        DATAeffeta1mm = getTriggerEffs(DATAnumeta1mm, DATAdeneta1mm, "data")
        DATAeffeta1SF = getTriggerEffs(DATAnumeta1SF, DATAdeneta1SF, "data")
        DATAeffeta1OF = getTriggerEffs(DATAnumeta1OF, DATAdeneta1OF, "data")
        DATAeffeta2ee = getTriggerEffs(DATAnumeta2ee, DATAdeneta2ee, "data")
        DATAeffeta2mm = getTriggerEffs(DATAnumeta2mm, DATAdeneta2mm, "data")
        DATAeffeta2SF = getTriggerEffs(DATAnumeta2SF, DATAdeneta2SF, "data")
        DATAeffeta2OF = getTriggerEffs(DATAnumeta2OF, DATAdeneta2OF, "data")
        
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
        
        DATARTMET  = RT(DATAnumMETSF, DATAeffMETee, DATAeffMETmm, DATAeffMETOF)
        DATARTmt2  = RT(DATAnummt2SF, DATAeffmt2ee, DATAeffmt2mm, DATAeffmt2OF)
        DATARTMll  = RT(DATAnumMllSF, DATAeffMllee, DATAeffMllmm, DATAeffMllOF)
        DATARTpt1  = RT(DATAnumpt1SF, DATAeffpt1ee, DATAeffpt1mm, DATAeffpt1OF)
        DATARTpt2  = RT(DATAnumpt2SF, DATAeffpt2ee, DATAeffpt2mm, DATAeffpt2OF)
        DATARTeta1 = RT(DATAnumeta1SF, DATAeffeta1ee, DATAeffeta1mm, DATAeffeta1OF)
        DATARTeta2 = RT(DATAnumeta2SF, DATAeffeta2ee, DATAeffeta2mm, DATAeffeta2OF)
        
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
        effRTMll.addLine(h_auxrtMll.GetXaxis().GetXmin(), dart, h_auxrtMll.GetXaxis().GetXmax(), dart ,r.kBlue)
        effRTMll.addHisto(h_auxrtMll, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
        effRTMll.addHisto(DATARTMll, 'PE,SAME', 'R_{T} data', 'PL', r.kBlack , 1, 0)
        effRTMll.addHisto(MCRTMll, 'PE,SAME', 'R_{T} MC', 'PL', r.kGreen , 1, 0)
        effRTMll.addBand(h_auxrtMll.GetXaxis().GetXmin(), dart-dasystrt, h_auxrtMll.GetXaxis().GetXmax(), dart+dasystrt, r.kOrange+6, 0.2)
        effRTMll.addLatex (0.55, 0.25, 'Mean R_{T} data: %.3f '%(dart))
        effRTMll.addLatex (0.55, 0.2, 'Mean R_{T} MC: %.3f '%(mcrt))
        effRTMll.save(1, 1, 0, lumi, 0.8, 1.2)                                                            
        
        effRTMET = Canvas.Canvas('rt/%s/plot_rt_met_'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
        h_auxrtMET = r.TH1F("h_auxRTMET", "", 1, 0, 200)
        h_auxrtMET.GetYaxis().SetRangeUser(0.8, 1.2)
        h_auxrtMET.GetXaxis().SetRangeUser(0, 150)
        h_auxrtMET.GetXaxis().SetTitle(labelmet)
        effRTMET.addLine(h_auxrtMET.GetXaxis().GetXmin(), dart, h_auxrtMET.GetXaxis().GetXmax(), dart ,r.kBlue)
        effRTMET.addHisto(h_auxrtMET, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
        effRTMET.addHisto(DATARTMET, 'PE,SAME', 'R_{T} data', 'PL', r.kBlack , 1, 0)
        effRTMET.addHisto(MCRTMET, 'PE,SAME', 'R_{T} MC', 'PL', r.kGreen , 1, 0)
        effRTMET.addBand(h_auxrtMET.GetXaxis().GetXmin(), dart-dasystrt, h_auxrtMET.GetXaxis().GetXmax(), dart+dasystrt, r.kOrange+6, 0.2)
        effRTMET.addLatex (0.6, 0.25, 'Mean R_{T} data: %.3f '%(dart))
        effRTMET.addLatex (0.6, 0.2, 'Mean R_{T} MC  : %.3f '%(mcrt))
        effRTMET.save(1, 1, 0, lumi, 0.8, 1.2)                                                                                                  
        
        
        effRTmt2 = Canvas.Canvas('rt/%s/plot_rt_mt2_'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
        h_auxrtmt2 = r.TH1F("h_auxRTMET", "", 1, 0, 160)
        h_auxrtmt2.GetYaxis().SetRangeUser(0.8, 1.2)
        h_auxrtmt2.GetXaxis().SetRangeUser(0, 160)
        h_auxrtmt2.GetXaxis().SetTitle(labelmt2)
        effRTmt2.addLine(h_auxrtmt2.GetXaxis().GetXmin(), dart, h_auxrtmt2.GetXaxis().GetXmax(), dart ,r.kBlue)
        effRTmt2.addHisto(h_auxrtmt2, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
        effRTmt2.addHisto(DATARTmt2, 'PE,SAME', 'R_{T} data', 'PL', r.kBlack , 1, 0)
        effRTmt2.addHisto(MCRTmt2, 'PE,SAME', 'R_{T} MC', 'PL', r.kGreen , 1, 0)
        effRTmt2.addBand(h_auxrtmt2.GetXaxis().GetXmin(), dart-dasystrt, h_auxrtmt2.GetXaxis().GetXmax(), dart+dasystrt, r.kOrange+6, 0.2)
        effRTmt2.addLatex (0.6, 0.25, 'Mean R_{T} data: %.3f '%(dart))
        effRTmt2.addLatex (0.6, 0.2, 'Mean R_{T} MC  : %.3f '%(mcrt))
        effRTmt2.save(1, 1, 0, lumi, 0.8, 1.2)                                                                                                 

        effRTpt1 = Canvas.Canvas('rt/%s/plot_rt_pt1_'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
        h_auxrtpt1 = r.TH1F("h_auxRTpt1", "", 1, 20, 150)
        h_auxrtpt1.GetYaxis().SetRangeUser(0.8, 1.2)
        h_auxrtpt1.GetXaxis().SetRangeUser(20, 150)
        h_auxrtpt1.GetXaxis().SetTitle(labelpt1)
        effRTpt1.addLine(h_auxrtpt1.GetXaxis().GetXmin(), dart, h_auxrtpt1.GetXaxis().GetXmax(), dart ,r.kBlue)
        effRTpt1.addHisto(h_auxrtpt1, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
        effRTpt1.addHisto(DATARTpt1, 'PE,SAME', 'RT data', 'PL', r.kBlack , 1, 0)
        effRTpt1.addHisto(MCRTpt1, 'PE,SAME', 'RT MC', 'PL', r.kGreen , 1, 0)
        effRTpt1.addBand(h_auxrtpt1.GetXaxis().GetXmin(), dart-dasystrt, h_auxrtpt1.GetXaxis().GetXmax(), dart+dasystrt, r.kOrange+6, 0.2)
        effRTpt1.addLatex (0.6, 0.25, 'Mean R_{T} data: %.3f '%(dart))
        effRTpt1.addLatex (0.6, 0.2, 'Mean R_{T} MC  : %.3f '%(mcrt))
        effRTpt1.save(1, 1, 0, lumi, 0.8, 1.2)

        effRTpt2 = Canvas.Canvas('rt/%s/plot_rt_pt2_'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
        h_auxrtpt2 = r.TH1F("h_auxRTpt2", "", 1, 20, 150)
        h_auxrtpt2.GetYaxis().SetRangeUser(0.8, 1.2)
        h_auxrtpt2.GetXaxis().SetRangeUser(20, 150)
        h_auxrtpt2.GetXaxis().SetTitle(labelpt2)
        effRTpt2.addLine(h_auxrtpt2.GetXaxis().GetXmin(), dart, h_auxrtpt2.GetXaxis().GetXmax(), dart ,r.kBlue)
        effRTpt2.addHisto(h_auxrtpt2, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
        effRTpt2.addHisto(DATARTpt2, 'PE,SAME', 'R_{T} data', 'PL', r.kBlack , 1, 0)
        effRTpt2.addHisto(MCRTpt2, 'PE,SAME', 'R_{T} MC', 'PL', r.kGreen , 1, 0)
        effRTpt2.addBand(h_auxrtpt2.GetXaxis().GetXmin(), dart-dasystrt, h_auxrtpt2.GetXaxis().GetXmax(), dart+dasystrt, r.kOrange+6, 0.2)
        effRTpt2.addLatex (0.6, 0.25, 'Mean R_{T} data: %.3f '%(dart))
        effRTpt2.addLatex (0.6, 0.2, 'Mean R_{T} MC: %.3f '%(mcrt))
        effRTpt2.save(1, 1, 0, lumi, 0.8, 1.2)

        effRTeta1 = Canvas.Canvas('rt/%s/plot_rt_eta1_'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
        h_auxrteta1 = r.TH1F("h_auxRTeta1", "", 1, -2.4, 2.4)
        h_auxrteta1.GetYaxis().SetRangeUser(0.8, 1.2)
        h_auxrteta1.GetXaxis().SetRangeUser(-2.4, 2.4)
        h_auxrteta1.GetXaxis().SetTitle(labeleta1)
        effRTeta1.addLine(h_auxrteta1.GetXaxis().GetXmin(), dart, h_auxrteta1.GetXaxis().GetXmax(), dart ,r.kBlue)
        effRTeta1.addHisto(h_auxrteta1, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
        effRTeta1.addHisto(DATARTeta1, 'PE,SAME', 'RT data', 'PL', r.kBlack , 1, 0)
        effRTeta1.addHisto(MCRTeta1, 'PE,SAME', 'RT MC', 'PL', r.kGreen , 1, 0)
        effRTeta1.addBand(h_auxrteta1.GetXaxis().GetXmin(), dart-dasystrt, h_auxrteta1.GetXaxis().GetXmax(), dart+dasystrt, r.kOrange+6, 0.2)
        effRTeta1.addLatex (0.6, 0.25, 'Mean R_{T} data: %.3f '%(dart))
        effRTeta1.addLatex (0.6, 0.2, 'Mean R_{T} MC  : %.3f '%(mcrt))
        effRTeta1.save(1, 1, 0, lumi, 0.8, 1.2)

        effRTeta2 = Canvas.Canvas('rt/%s/plot_rt_eta2_'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
        h_auxrteta2 = r.TH1F("h_auxRTeta2", "", 1, -2.4, 2.4)
        h_auxrteta2.GetYaxis().SetRangeUser(0.8, 1.2)
        h_auxrteta2.GetXaxis().SetRangeUser(-2.4, 2.4)
        h_auxrteta2.GetXaxis().SetTitle(labeleta2)
        effRTeta2.addLine(h_auxrteta2.GetXaxis().GetXmin(), dart, h_auxrteta2.GetXaxis().GetXmax(), dart ,r.kBlue)
        effRTeta2.addHisto(h_auxrteta2, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
        effRTeta2.addHisto(DATARTeta2, 'PE,SAME', 'RT data' , 'PL', r.kBlack , 1, 0)
        effRTeta2.addHisto(MCRTeta2, 'PE,SAME', 'RT MC', 'PL', r.kGreen , 1, 0)
        effRTeta2.addBand(h_auxrteta2.GetXaxis().GetXmin(), dart-dasystrt, h_auxrteta2.GetXaxis().GetXmax(), dart+dasystrt, r.kOrange+6, 0.2)
        effRTeta2.addLatex (0.6, 0.25, 'Mean R_{T} data: %.3f '%(dart))
        effRTeta2.addLatex (0.6, 0.2, 'Mean R_{T} MC  : %.3f '%(mcrt))
        effRTeta2.save(1, 1, 0, lumi, 0.8, 1.2)
   
        effMll = Canvas.Canvas('rt/%s/plot_eff_mll_'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
        h_auxMll = r.TH1F("h_auxMll", "", 1, 0, 250)
        h_auxMll.GetYaxis().SetRangeUser(0, 2)
        h_auxMll.GetXaxis().SetRangeUser(0, 250)
        h_auxMll.GetXaxis().SetTitle(labelx)
        effMll.addHisto(h_auxMll, 'h', '', '', r.kRed+1, 1, 0)
        effMll.addHisto(DATAeffMllee, 'PE,SAME', 'Double Electron', 'PL', r.kRed+1 , 1, 0)
        effMll.addHisto(DATAeffMllmm, 'PE,SAME', 'Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effMll.addHisto(DATAeffMllSF, 'PE,SAME', 'Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effMll.addHisto(DATAeffMllOF, 'PE,SAME', 'Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effMll.save(1, 1, 0, lumi, 0.2, 1.8)                                                                    

        effmt2 = Canvas.Canvas('rt/%s/plot_eff_mt2_'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
        h_auxmt2 = r.TH1F("h_auxmt2", "", 1, 0, 160)
        h_auxmt2.GetYaxis().SetRangeUser(0, 2)
        h_auxmt2.GetXaxis().SetRangeUser(0, 160)
        h_auxmt2.GetXaxis().SetTitle(labelx)
        effmt2.addHisto(h_auxmt2, 'h', '', '', r.kRed+1, 1, 0)
        effmt2.addHisto(DATAeffmt2ee, 'PE,SAME', 'Double Electron', 'PL', r.kRed+1 , 1, 0)
        effmt2.addHisto(DATAeffmt2mm, 'PE,SAME', 'Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effmt2.addHisto(DATAeffmt2SF, 'PE,SAME', 'Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effmt2.addHisto(DATAeffmt2OF, 'PE,SAME', 'Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effmt2.save(1, 1, 0, lumi, 0.2, 1.8)                                                                    

        effpt1 = Canvas.Canvas('rt/%s/plot_eff_pt1_'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
        h_auxpt1 = r.TH1F("h_auxpt1", "", 1, 0, 150)
        h_auxpt1.GetYaxis().SetRangeUser(0, 2)
        h_auxpt1.GetXaxis().SetRangeUser(0, 150)
        h_auxpt1.GetXaxis().SetTitle(labelpt1)
        effpt1.addHisto(h_auxpt1, 'h', '', '', r.kRed+1, 1, 0)
        effpt1.addHisto(DATAeffpt1ee, 'PE,SAME', 'Double Electron', 'PL', r.kRed+1 , 1, 0)
        effpt1.addHisto(DATAeffpt1mm, 'PE,SAME', 'Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effpt1.addHisto(DATAeffpt1SF, 'PE,SAME', 'Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effpt1.addHisto(DATAeffpt1OF, 'PE,SAME', 'Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effpt1.save(1, 1, 0, lumi, 0.2, 1.8)

        effpt2 = Canvas.Canvas('rt/%s/plot_eff_pt2_'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
        h_auxpt2 = r.TH1F("h_auxpt2", "", 1, 0, 150)
        h_auxpt2.GetYaxis().SetRangeUser(0, 2)
        h_auxpt2.GetXaxis().SetRangeUser(0, 150)
        h_auxpt2.GetXaxis().SetTitle(labelpt2)
        effpt2.addHisto(h_auxpt2, 'h', '', '', r.kRed+1, 1, 0)
        effpt2.addHisto(DATAeffpt2ee, 'PE,SAME', 'Double Electron', 'PL', r.kRed+1 , 1, 0)
        effpt2.addHisto(DATAeffpt2mm, 'PE,SAME', 'Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effpt2.addHisto(DATAeffpt2SF, 'PE,SAME', 'Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effpt2.addHisto(DATAeffpt2OF, 'PE,SAME', 'Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effpt2.save(1, 1, 0, lumi, 0.2, 1.8)
 
        effeta1 = Canvas.Canvas('rt/%s/plot_eff_eta1_'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
        h_auxeta1 = r.TH1F("h_auxeta1", "", 1, -2.4, 2.4)
        h_auxeta1.GetYaxis().SetRangeUser(0, 2)
        h_auxeta1.GetXaxis().SetRangeUser(-2.4, 2.4)
        h_auxeta1.GetXaxis().SetTitle(labeleta1)
        effeta1.addHisto(h_auxeta1, 'h', '', '', r.kRed+1, 1, 0)
        effeta1.addHisto(DATAeffeta1ee, 'PE,SAME', 'Double Electron', 'PL', r.kRed+1 , 1, 0)
        effeta1.addHisto(DATAeffeta1mm, 'PE,SAME', 'Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effeta1.addHisto(DATAeffeta1SF, 'PE,SAME', 'Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effeta1.addHisto(DATAeffeta1OF, 'PE,SAME', 'Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effeta1.save(1, 1, 0, lumi, 0.2, 1.8)

        effeta2 = Canvas.Canvas('rt/%s/plot_eff_eta2_'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
        h_auxeta2 = r.TH1F("h_auxeta2", "", 1, -2.4, 2.4)
        h_auxeta2.GetYaxis().SetRangeUser(0, 2)
        h_auxeta2.GetXaxis().SetRangeUser(-2.4, 2.4)
        h_auxeta2.GetXaxis().SetTitle(labeleta2)
        effeta2.addHisto(h_auxeta2, 'h', '', '', r.kRed+1, 1, 0)
        effeta2.addHisto(DATAeffeta2ee, 'PE,SAME', 'Double Electron', 'PL', r.kRed+1 , 1, 0)
        effeta2.addHisto(DATAeffeta2mm, 'PE,SAME', 'Double Muon', 'PL', r.kBlue+1 , 1, 0)
        effeta2.addHisto(DATAeffeta2SF, 'PE,SAME', 'Same Flavor', 'PL', r.kGreen+1 , 1, 0)
        effeta2.addHisto(DATAeffeta2OF, 'PE,SAME', 'Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
        effeta2.save(1, 1, 0, lumi, 0.2, 1.8)
 
        print 'Measured RT value data ', dart, ' +/- ', dauncrt, ' +/- ', dasystrt
        print 'Measured RT value MC   ', mcrt, ' +/- ', mcuncrt, ' +/- ', mcsystrt


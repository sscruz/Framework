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




def makeTable(DATAnumeratoree, DATAnumeratormm, DATAnumeratorOF, DATAdenominatoree, DATAdenominatormm, DATAdenominatorOF, effMlleevalue, effMllmmvalue , effMllOFvalue, eevale, mmvale, emvale, dataMC):
    line0 = '  \hline'
    line1 = '  ee &' 
    line2 = '  \mu\mu   &' 
    line3 = '  e\mu   &' 
    line0 += ' & nominator  & denominator & $\epsilon_{\mathrm{trigger}}$ $\pm \sigma_{stat}}$ \\\\'
    line1 += ' %.f & %.f &    %.3f $\\pm$ %.3f       %s' %(DATAnumeratoree.Integral(), DATAdenominatoree.Integral(), effMlleevalue  , eevale  ,   '\\\\')
    line2 += ' %.f & %.f &    %.3f $\\pm$ %.3f       %s' %(DATAnumeratormm.Integral(), DATAdenominatormm.Integral(), effMllmmvalue  , mmvale  ,   '\\\\')
    line3 += ' %.f & %.f &    %.3f $\\pm$ %.3f       %s' %(DATAnumeratorOF.Integral(), DATAdenominatorOF.Integral(), effMllOFvalue  , emvale  ,   '\\\\')
    line0 += '\\hline'; line3 += '\\hline';
                                                                                                                                                                                     
    helper.ensureDirectory('plots/rt/%s/'%lumi_str)
    helper.ensureDirectory('plots/rt/%s/tables/'%lumi_str)
    compTableFile = open('plots/rt/%s/tables/resultTable_%s%s.txt'%(lumi_str, str(lumi).replace('.','p'), "rt"),'w')
    compTableFile.write(line0+'\n')
    compTableFile.write(line1+'\n')
    compTableFile.write(line2+'\n')                                                                                             
    compTableFile.write(line3+'\n')                                                                                             
    if dataMC == 'DATA':
        print "RT for Data"
    else:
        print "RT for MC"
    print line0
    print line1
    print line2                                                                                                                                                                      
    print line3                                                                                                                                                                                         


############################################################
def getRT(eff_ee, unc_ee, eff_mm, unc_mm, eff_em, unc_em):
    if(eff_em == 0 ):
        print "Division by zero"
        return [0, 0, 0]
    if (eff_mm == 0):
        eff_mm = 0.0001
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


def getTriggerEffs(num, den):

    trig = den.Clone("trig_" + den.GetName())
    trig.Reset()



    errs = TGraphAsymmErrors(num, den, 'v')
    #x = errs.GetX()
    #y = errs.GetY()
    #eyh = errs.GetEYhigh()
    #eyl = errs.GetEYlow()
    #for i in range(1, num.GetNbinsX()+1):
    #    if(x[i] != 0 or y[i] != 0):
    #        trig.SetBinContent(i, y[i])
    #        trig.SetBinError(i, max(eyh[i], eyl[i]))

    return errs


def saveInFile(theFile, measuredValueMC, measuredValueUncMC, measuredValueSystMC, measuredValueData, measuredValueUncData, measuredValueSystData, dataMC):


    foutput = open(theFile + "_aux", 'w')
    for line in open(theFile).readlines():
        if line.find("rt") != -1:
            if line.find(dataMC) != -1:
                foutput.write('rt          region          %s        %.4f      %0.4f       %.4f\n'%(dataMC, measuredValueData, measuredValueUncData, measuredValueSystData))
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
    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    
    mcDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO', 'TTJets_DiLepton', 'ZZTo4L','GGHZZ4L',  'WZTo3LNu', 'WWW', 'WWZ','ZZZ', 'tZq_ll','WWTo2L2Nu', 'ZZTo2L2Nu', 'WZTo2L2Q','TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT', 'TTTT',  'TTWToQQ', 'TTZToLLNuNu' ,'TTWToLNu', 'WJetsToLNu_LO']
    
    daDatasets = ['DoubleEG_Run2016H-PromptReco-v2_runs_281207_284035_part2',             
                'DoubleEG_Run2016H-PromptReco-v3_runs_284036_284044',        
                'JetHT_Run2016B_23Sep2016_v3_runs_273150_275376',      
                'JetHT_Run2016C_23Sep2016_v1_runs_271036_284044',                
                'JetHT_Run2016D_23Sep2016_v1_runs_271036_284044',               
                'JetHT_Run2016E_23Sep2016_v1_runs_271036_284044',             
                'DoubleEG_Run2016D_23Sep2016_v1_runs_271036_284044',                 
                'JetHT_Run2016F_23Sep2016_v1_runs_271036_284044',                    
                'JetHT_Run2016G_23Sep2016_v1_runs_271036_284044',                    
                'JetHT_Run2016H-PromptReco-v2_runs_281207_284035',                   
                'JetHT_Run2016H-PromptReco-v3_runs_284036_284044',                   
                'MET_Run2016B_23Sep2016_v3_runs_273150_275376',       
                'MET_Run2016C_23Sep2016_v1_runs_271036_284044',       
                'MET_Run2016D_23Sep2016_v1_runs_271036_284044',       
                'MET_Run2016E_23Sep2016_v1_runs_271036_284044',       
                'MET_Run2016F_23Sep2016_v1_runs_271036_284044',       
                'DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part1',       
                'MET_Run2016G_23Sep2016_v1_runs_271036_284044',       
                'MET_Run2016H-PromptReco-v3_runs_284036_284044',       
                'MET_Run2016H-PromptReco-v2_runs_281207_284035',       
                'DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376_part1',       
                'DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376_part2',       
                'DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part2',       
                'DoubleEG_Run2016E_23Sep2016_v1_runs_271036_284044',       
                'DoubleEG_Run2016C_23Sep2016_v1_runs_271036_284044',       
                'SingleElectron_Run2016B_23Sep2016_v3_runs_273150_275376',       
                'DoubleEG_Run2016F_23Sep2016_v1_runs_271036_284044',       
                'SingleElectron_Run2016E_23Sep2016_v1_runs_271036_284044',       
                'SingleElectron_Run2016C_23Sep2016_v1_runs_271036_284044',       
                'SingleElectron_Run2016F_23Sep2016_v1_runs_271036_284044',       
                'SingleElectron_Run2016H-PromptReco-v2_runs_281207_284035',       
                'SingleElectron_Run2016H-PromptReco-v3_runs_284036_284044',       
                'DoubleEG_Run2016H-PromptReco-v2_runs_281207_284035_part1',       
                'SingleElectron_Run2016D_23Sep2016_v1_runs_271036_284044',       
                'SingleElectron_Run2016G_23Sep2016_v1_runs_271036_284044',       
                'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part2',       
                'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part2',       
                'DoubleMuon_Run2016C_23Sep2016_v1_runs_271036_284044',       
                'DoubleMuon_Run2016F_23Sep2016_v1_runs_271036_284044',       
                'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part1',       
                'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part1',       
                'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part3',       
                'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part2',       
                'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part1',       
                'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part1',       
                'DoubleMuon_Run2016H-PromptReco-v3_runs_284036_284044',       
                'MuonEG_Run2016B_23Sep2016_v3_runs_273150_275376',       
                'DoubleMuon_Run2016H-PromptReco-v2_runs_281207_284035_part1',       
                'MuonEG_Run2016C_23Sep2016_v1_runs_271036_284044',       
                'DoubleMuon_Run2016H-PromptReco-v2_runs_281207_284035_part2',       
                'MuonEG_Run2016E_23Sep2016_v1_runs_271036_284044',       
                'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part3',       
                'MuonEG_Run2016D_23Sep2016_v1_runs_271036_284044',       
                'MuonEG_Run2016H-PromptReco-v3_runs_284036_284044',       
                'MuonEG_Run2016F_23Sep2016_v1_runs_271036_284044',       
                'DoubleMuon_Run2016H-PromptReco-v2_runs_281207_284035_part3',       
                'MuonEG_Run2016G_23Sep2016_v1_runs_271036_284044',       
                'SingleMuon_Run2016C_23Sep2016_v1_runs_271036_284044',       
                'SingleMuon_Run2016E_23Sep2016_v1_runs_271036_284044',       
                'SingleMuon_Run2016B_23Sep2016_v3_runs_273150_275376',       
                'SingleMuon_Run2016D_23Sep2016_v1_runs_271036_284044',       
                'SingleMuon_Run2016F_23Sep2016_v1_runs_271036_284044',       
                'SingleMuon_Run2016H-PromptReco-v3_runs_284036_284044',       
                'MuonEG_Run2016H-PromptReco-v2_runs_281207_284035',       
                'SingleMuon_Run2016G_23Sep2016_v1_runs_271036_284044',       
                'SingleMuon_Run2016H-PromptReco-v2_runs_281207_284035']               

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    
    lumi = 36.2 ; maxrun = 999999; lumi_str = '36.2invfb'
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

    DATAdenominatorMllee = treeDA.getTH1F(lumi, "DATAdenominatoree", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelx)
    print "DATAdenominatorMllee", DATAdenominatorMllee.Integral()
    DATAnumeratorMllee =   treeDA.getTH1F(lumi, "DATAnumeratoree", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelx)
    print "DATAnumeratorMllee", DATAnumeratorMllee.Integral()
    DATAdenominatorMllmm = treeDA.getTH1F(lumi, "DATAdenominatormm", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelx)
    DATAnumeratorMllmm =   treeDA.getTH1F(lumi, "DATAnumeratormm", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelx)
    DATAdenominatorMllSF = treeDA.getTH1F(lumi, "DATAdenominatorSF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelx)
    DATAnumeratorMllSF =   treeDA.getTH1F(lumi, "DATAnumeratorSF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelx)
    DATAdenominatorMllOF = treeDA.getTH1F(lumi, "DATAdenominatorOF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelx)
    DATAnumeratorMllOF =   treeDA.getTH1F(lumi, "DATAnumeratorOF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelx)
    
    DATAdenominatorpt1ee = treeDA.getTH1F(lumi, "DATAdenominatorpt1ee", "Lep1_pt_Edge", ptbins , 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelpt1)
    DATAnumeratorpt1ee =   treeDA.getTH1F(lumi, "DATAnumeratorpt1ee", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelpt1)
    DATAdenominatorpt1mm = treeDA.getTH1F(lumi, "DATAdenominatorpt1mm", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelpt1)
    DATAnumeratorpt1mm =   treeDA.getTH1F(lumi, "DATAnumeratorpt1mm", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelpt1)
    DATAdenominatorpt1SF = treeDA.getTH1F(lumi, "DATAdenominatorpt1SF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelpt1)
    DATAnumeratorpt1SF =   treeDA.getTH1F(lumi, "DATAnumeratorpt1SF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelpt1)
    DATAdenominatorpt1OF = treeDA.getTH1F(lumi, "DATAdenominatorpt1OF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelpt1)
    DATAnumeratorpt1OF =   treeDA.getTH1F(lumi, "DATAnumeratorpt1OF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelpt1)
 
    DATAdenominatorpt2ee = treeDA.getTH1F(lumi, "DATAdenominatorpt2ee", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelpt2)
    DATAnumeratorpt2ee =   treeDA.getTH1F(lumi, "DATAnumeratorpt2ee", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelpt2)
    DATAdenominatorpt2mm = treeDA.getTH1F(lumi, "DATAdenominatorpt2mm", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelpt2)
    DATAnumeratorpt2mm =   treeDA.getTH1F(lumi, "DATAnumeratorpt2mm", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelpt2)
    DATAdenominatorpt2SF = treeDA.getTH1F(lumi, "DATAdenominatorpt2SF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelpt2)
    DATAnumeratorpt2SF =   treeDA.getTH1F(lumi, "DATAnumeratorpt2SF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelpt2)
    DATAdenominatorpt2OF = treeDA.getTH1F(lumi, "DATAdenominatorpt2OF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelpt2)
    DATAnumeratorpt2OF =   treeDA.getTH1F(lumi, "DATAnumeratorpt2OF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelpt2)
 
    DATAdenominatoreta1ee = treeDA.getTH1F(lumi, "DATAdenominatoreta1ee", "Lep1_eta_Edge", etabins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labeleta1)
    DATAnumeratoreta1ee =   treeDA.getTH1F(lumi, "DATAnumeratoreta1ee", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labeleta1)
    DATAdenominatoreta1mm = treeDA.getTH1F(lumi, "DATAdenominatoreta1mm", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labeleta1)
    DATAnumeratoreta1mm =   treeDA.getTH1F(lumi, "DATAnumeratoreta1mm", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labeleta1)
    DATAdenominatoreta1SF = treeDA.getTH1F(lumi, "DATAdenominatoreta1SF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labeleta1)
    DATAnumeratoreta1SF =   treeDA.getTH1F(lumi, "DATAnumeratoreta1SF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labeleta1)
    DATAdenominatoreta1OF = treeDA.getTH1F(lumi, "DATAdenominatoreta1OF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labeleta1)
    DATAnumeratoreta1OF =   treeDA.getTH1F(lumi, "DATAnumeratoreta1OF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labeleta1)
 
    DATAdenominatoreta2ee = treeDA.getTH1F(lumi, "DATAdenominatoreta2ee", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labeleta2)
    DATAnumeratoreta2ee =   treeDA.getTH1F(lumi, "DATAnumeratoreta2ee", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labeleta2)
    DATAdenominatoreta2mm = treeDA.getTH1F(lumi, "DATAdenominatoreta2mm", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labeleta2)
    DATAnumeratoreta2mm =   treeDA.getTH1F(lumi, "DATAnumeratoreta2mm", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labeleta2)
    DATAdenominatoreta2SF = treeDA.getTH1F(lumi, "DATAdenominatoreta2SF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labeleta2)
    DATAnumeratoreta2SF =   treeDA.getTH1F(lumi, "DATAnumeratoreta2SF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labeleta2)
    DATAdenominatoreta2OF = treeDA.getTH1F(lumi, "DATAdenominatoreta2OF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labeleta2)
    DATAnumeratoreta2OF =   treeDA.getTH1F(lumi, "DATAnumeratoreta2OF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labeleta2)
 
    DATAdenominatorMlleevalue =   treeDA.getTH1F(lumi, "DATAdenominatoreevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelx)
    DATAnumeratorMlleevalue =     treeDA.getTH1F(lumi, "DATAnumeratoreevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelx)
    DATAdenominatorMllmmvalue =   treeDA.getTH1F(lumi, "DATAdenominatormmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelx)
    DATAnumeratorMllmmvalue =     treeDA.getTH1F(lumi, "DATAnumeratormmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelx)
    DATAdenominatorMllSFvalue =   treeDA.getTH1F(lumi, "DATAdenominatorSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelx)
    DATAnumeratorMllSFvalue =     treeDA.getTH1F(lumi, "DATAnumeratorSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelx)
    DATAdenominatorMllOFvalue =   treeDA.getTH1F(lumi, "DATAdenominatorOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelx)
    DATAnumeratorMllOFvalue =     treeDA.getTH1F(lumi, "DATAnumeratorOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelx)
    
    DATAdenominatorMETee =   treeDA.getTH1F(lumi, "DATAdenominatoreevalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelmet)
    DATAnumeratorMETee =     treeDA.getTH1F(lumi, "DATAnumeratoreevalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelmet)
    DATAdenominatorMETmm =   treeDA.getTH1F(lumi, "DATAdenominatormmvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelmet)
    DATAnumeratorMETmm =     treeDA.getTH1F(lumi, "DATAnumeratormmvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelmet)
    DATAdenominatorMETSF =   treeDA.getTH1F(lumi, "DATAdenominatorSFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelmet)
    DATAnumeratorMETSF =     treeDA.getTH1F(lumi, "DATAnumeratorSFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelmet)
    DATAdenominatorMETOF =   treeDA.getTH1F(lumi, "DATAdenominatorOFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelmet)
    DATAnumeratorMETOF =     treeDA.getTH1F(lumi, "DATAnumeratorOFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelmet)

    DATAdenominatormt2ee =   treeDA.getTH1F(lumi, "DATAdenominatoreevalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelmt2)
    DATAnumeratormt2ee =     treeDA.getTH1F(lumi, "DATAnumeratoreevalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelmt2)
    DATAdenominatormt2mm =   treeDA.getTH1F(lumi, "DATAdenominatormmvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelmt2)
    DATAnumeratormt2mm =     treeDA.getTH1F(lumi, "DATAnumeratormmvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelmt2)
    DATAdenominatormt2SF =   treeDA.getTH1F(lumi, "DATAdenominatorSFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelmt2)
    DATAnumeratormt2SF =     treeDA.getTH1F(lumi, "DATAnumeratorSFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelmt2)
    DATAdenominatormt2OF =   treeDA.getTH1F(lumi, "DATAdenominatorOFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelmt2)
    DATAnumeratormt2OF =     treeDA.getTH1F(lumi, "DATAnumeratorOFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelmt2)                

    MCdenominatorMllee = treeMC.getTH1F(lumi, "MCdenominatoree", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelx)
    MCnumeratorMllee =   treeMC.getTH1F(lumi, "MCnumeratoree", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelx)
    MCdenominatorMllmm = treeMC.getTH1F(lumi, "MCdenominatormm", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelx)
    MCnumeratorMllmm =   treeMC.getTH1F(lumi, "MCnumeratormm", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelx)
    MCdenominatorMllSF = treeMC.getTH1F(lumi, "MCdenominatorSF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelx)
    MCnumeratorMllSF =   treeMC.getTH1F(lumi, "MCnumeratorSF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelx)
    MCdenominatorMllOF = treeMC.getTH1F(lumi, "MCdenominatorOF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelx)
    MCnumeratorMllOF =   treeMC.getTH1F(lumi, "MCnumeratorOF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelx)
    
    MCdenominatorpt1ee = treeMC.getTH1F(lumi, "MCdenominatorpt1ee", "Lep1_pt_Edge", ptbins , 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelpt1)
    MCnumeratorpt1ee =   treeMC.getTH1F(lumi, "MCnumeratorpt1ee", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelpt1)
    MCdenominatorpt1mm = treeMC.getTH1F(lumi, "MCdenominatorpt1mm", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelpt1)
    MCnumeratorpt1mm =   treeMC.getTH1F(lumi, "MCnumeratorpt1mm", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelpt1)
    MCdenominatorpt1SF = treeMC.getTH1F(lumi, "MCdenominatorpt1SF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelpt1)
    MCnumeratorpt1SF =   treeMC.getTH1F(lumi, "MCnumeratorpt1SF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelpt1)
    MCdenominatorpt1OF = treeMC.getTH1F(lumi, "MCdenominatorpt1OF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelpt1)
    MCnumeratorpt1OF =   treeMC.getTH1F(lumi, "MCnumeratorpt1OF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelpt1)
    
    MCdenominatorpt2ee = treeMC.getTH1F(lumi, "MCdenominatorpt2ee", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelpt2)
    MCnumeratorpt2ee =   treeMC.getTH1F(lumi, "MCnumeratorpt2ee", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelpt2)
    MCdenominatorpt2mm = treeMC.getTH1F(lumi, "MCdenominatorpt2mm", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelpt2)
    MCnumeratorpt2mm =   treeMC.getTH1F(lumi, "MCnumeratorpt2mm", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelpt2)
    MCdenominatorpt2SF = treeMC.getTH1F(lumi, "MCdenominatorpt2SF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelpt2)
    MCnumeratorpt2SF =   treeMC.getTH1F(lumi, "MCnumeratorpt2SF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelpt2)
    MCdenominatorpt2OF = treeMC.getTH1F(lumi, "MCdenominatorpt2OF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelpt2)
    MCnumeratorpt2OF =   treeMC.getTH1F(lumi, "MCnumeratorpt2OF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelpt2)
   
    MCdenominatoreta1ee = treeMC.getTH1F(lumi, "MCdenominatoreta1ee", "Lep1_eta_Edge", etabins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labeleta1)
    MCnumeratoreta1ee =   treeMC.getTH1F(lumi, "MCnumeratoreta1ee", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labeleta1)
    MCdenominatoreta1mm = treeMC.getTH1F(lumi, "MCdenominatoreta1mm", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labeleta1)
    MCnumeratoreta1mm =   treeMC.getTH1F(lumi, "MCnumeratoreta1mm", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labeleta1)
    MCdenominatoreta1SF = treeMC.getTH1F(lumi, "MCdenominatoreta1SF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labeleta1)
    MCnumeratoreta1SF =   treeMC.getTH1F(lumi, "MCnumeratoreta1SF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labeleta1)
    MCdenominatoreta1OF = treeMC.getTH1F(lumi, "MCdenominatoreta1OF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labeleta1)
    MCnumeratoreta1OF =   treeMC.getTH1F(lumi, "MCnumeratoreta1OF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labeleta1)
    
    MCdenominatoreta2ee = treeMC.getTH1F(lumi, "MCdenominatoreta2ee", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labeleta2)
    MCnumeratoreta2ee =   treeMC.getTH1F(lumi, "MCnumeratoreta2ee", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labeleta2)
    MCdenominatoreta2mm = treeMC.getTH1F(lumi, "MCdenominatoreta2mm", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labeleta2)
    MCnumeratoreta2mm =   treeMC.getTH1F(lumi, "MCnumeratoreta2mm", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labeleta2)
    MCdenominatoreta2SF = treeMC.getTH1F(lumi, "MCdenominatoreta2SF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labeleta2)
    MCnumeratoreta2SF =   treeMC.getTH1F(lumi, "MCnumeratoreta2SF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labeleta2)
    MCdenominatoreta2OF = treeMC.getTH1F(lumi, "MCdenominatoreta2OF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labeleta2)
    MCnumeratoreta2OF =   treeMC.getTH1F(lumi, "MCnumeratoreta2OF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labeleta2)
                                                                                                                                                                                                 
    MCdenominatorMlleevalue =   treeMC.getTH1F(lumi, "MCdenominatoreevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelx)
    MCnumeratorMlleevalue =     treeMC.getTH1F(lumi, "MCnumeratoreevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelx)
    MCdenominatorMllmmvalue =   treeMC.getTH1F(lumi, "MCdenominatormmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelx)
    MCnumeratorMllmmvalue =     treeMC.getTH1F(lumi, "MCnumeratormmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelx)
    MCdenominatorMllSFvalue =   treeMC.getTH1F(lumi, "MCdenominatorSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelx)
    MCnumeratorMllSFvalue =     treeMC.getTH1F(lumi, "MCnumeratorSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelx)
    MCdenominatorMllOFvalue =   treeMC.getTH1F(lumi, "MCdenominatorOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelx)
    MCnumeratorMllOFvalue =     treeMC.getTH1F(lumi, "MCnumeratorOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelx)
    
    MCdenominatorMETee =   treeMC.getTH1F(lumi, "MCdenominatoreevalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelmet)
    MCnumeratorMETee =     treeMC.getTH1F(lumi, "MCnumeratoreevalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelmet)
    MCdenominatorMETmm =   treeMC.getTH1F(lumi, "MCdenominatormmvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelmet)
    MCnumeratorMETmm =     treeMC.getTH1F(lumi, "MCnumeratormmvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelmet)
    MCdenominatorMETSF =   treeMC.getTH1F(lumi, "MCdenominatorSFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelmet)
    MCnumeratorMETSF =     treeMC.getTH1F(lumi, "MCnumeratorSFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelmet)
    MCdenominatorMETOF =   treeMC.getTH1F(lumi, "MCdenominatorOFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelmet)
    MCnumeratorMETOF =     treeMC.getTH1F(lumi, "MCnumeratorOFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelmet)
    
    MCdenominatormt2ee =   treeMC.getTH1F(lumi, "MCdenominatoreevalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelmt2)
    MCnumeratormt2ee =     treeMC.getTH1F(lumi, "MCnumeratoreevalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelmt2)
    MCdenominatormt2mm =   treeMC.getTH1F(lumi, "MCdenominatormmvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelmt2)
    MCnumeratormt2mm =     treeMC.getTH1F(lumi, "MCnumeratormmvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelmt2)
    MCdenominatormt2SF =   treeMC.getTH1F(lumi, "MCdenominatorSFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelmt2)
    MCnumeratormt2SF =     treeMC.getTH1F(lumi, "MCnumeratorSFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelmt2)
    MCdenominatormt2OF =   treeMC.getTH1F(lumi, "MCdenominatorOFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelmt2)
    MCnumeratormt2OF =     treeMC.getTH1F(lumi, "MCnumeratorOFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelmt2)                ######################## Calculation of total values #############################
    effMlleevalueDATA = getTriggerEffs(DATAnumeratorMlleevalue, DATAdenominatorMlleevalue)
    effMllmmvalueDATA = getTriggerEffs(DATAnumeratorMllmmvalue, DATAdenominatorMllmmvalue)
    effMllSFvalueDATA = getTriggerEffs(DATAnumeratorMllSFvalue, DATAdenominatorMllSFvalue)
    effMllOFvalueDATA = getTriggerEffs(DATAnumeratorMllOFvalue, DATAdenominatorMllOFvalue)
    effMlleevalueMC   = getTriggerEffs(MCnumeratorMlleevalue, MCdenominatorMlleevalue)
    effMllmmvalueMC   = getTriggerEffs(MCnumeratorMllmmvalue, MCdenominatorMllmmvalue)
    effMllSFvalueMC   = getTriggerEffs(MCnumeratorMllSFvalue, MCdenominatorMllSFvalue)
    effMllOFvalueMC   = getTriggerEffs(MCnumeratorMllOFvalue, MCdenominatorMllOFvalue)

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
    eevalmceh_ = effMlleevalueMC.GetEYhigh()
    eevalmcel_ = effMlleevalueMC.GetEYlow()
    eevalmc = eevalmc_[0]
    eevalmce = max(eevalmceh_[0], eevalmcel_[0])
    mmvalmc_ = effMllmmvalueMC.GetY()
    mmvalmceh_ = effMllmmvalueMC.GetEYhigh()
    mmvalmcel_ = effMllmmvalueMC.GetEYlow()
    mmvalmc = mmvalmc_[0]
    mmvalmce = max(mmvalmceh_[0], mmvalmcel_[0])
    emvalmc_ = effMllOFvalueMC.GetY()
    emvalmceh_ = effMllOFvalueMC.GetEYhigh()
    emvalmcel_ = effMllOFvalueMC.GetEYlow()
    emvalmc = emvalmc_[0]
    emvalmce = max(emvalmceh_[0], emvalmcel_[0])

    [dart, dauncrt, dasystrt] = getRT(eevalda, eevaldae, mmvalda, mmvaldae, emvalda, emvaldae)
    [mcrt, mcuncrt, mcsystrt] = getRT(eevalmc, eevalmce, mmvalmc, mmvalmce, emvalmc, emvalmce)

    makeTable(DATAnumeratorMlleevalue, DATAnumeratorMllmmvalue, DATAnumeratorMllOFvalue, DATAdenominatorMlleevalue, DATAdenominatorMllmmvalue, DATAdenominatorMllOFvalue, eevalda, mmvalda, emvalda, eevaldae, mmvaldae, emvaldae, "DATA")
    makeTable(MCnumeratorMlleevalue, MCnumeratorMllmmvalue, MCnumeratorMllOFvalue, MCdenominatorMlleevalue, MCdenominatorMllmmvalue, MCdenominatorMllOFvalue, eevalmc, mmvalmc, emvalmc, eevalmce, mmvalmce, emvalmce, "MC")
    ######################## Calculation of total values #############################

    DATAeffMETee =  getTriggerEffs(DATAnumeratorMETee, DATAdenominatorMETee)
    DATAeffMETmm =  getTriggerEffs(DATAnumeratorMETmm, DATAdenominatorMETmm)    
    DATAeffMETSF =  getTriggerEffs(DATAnumeratorMETSF, DATAdenominatorMETSF)    
    DATAeffMETOF =  getTriggerEffs(DATAnumeratorMETOF, DATAdenominatorMETOF)
   
    DATAeffmt2ee =  getTriggerEffs(DATAnumeratormt2ee, DATAdenominatormt2ee)
    DATAeffmt2mm =  getTriggerEffs(DATAnumeratormt2mm, DATAdenominatormt2mm)
    DATAeffmt2SF =  getTriggerEffs(DATAnumeratormt2SF, DATAdenominatormt2SF)
    DATAeffmt2OF =  getTriggerEffs(DATAnumeratormt2OF, DATAdenominatormt2OF)

    DATAeffMllee =  getTriggerEffs(DATAnumeratorMllee, DATAdenominatorMllee)
    DATAeffMllmm =  getTriggerEffs(DATAnumeratorMllmm, DATAdenominatorMllmm)
    DATAeffMllSF =  getTriggerEffs(DATAnumeratorMllSF, DATAdenominatorMllSF)
    DATAeffMllOF =  getTriggerEffs(DATAnumeratorMllOF, DATAdenominatorMllOF)
    
    DATAeffpt1ee =  getTriggerEffs(DATAnumeratorpt1ee, DATAdenominatorpt1ee)
    DATAeffpt1mm =  getTriggerEffs(DATAnumeratorpt1mm, DATAdenominatorpt1mm)
    DATAeffpt1SF =  getTriggerEffs(DATAnumeratorpt1SF, DATAdenominatorpt1SF)
    DATAeffpt1OF =  getTriggerEffs(DATAnumeratorpt1OF, DATAdenominatorpt1OF)
    
    DATAeffpt2ee =  getTriggerEffs(DATAnumeratorpt2ee, DATAdenominatorpt2ee)
    DATAeffpt2mm =  getTriggerEffs(DATAnumeratorpt2mm, DATAdenominatorpt2mm)
    DATAeffpt2SF =  getTriggerEffs(DATAnumeratorpt2SF, DATAdenominatorpt2SF)
    DATAeffpt2OF =  getTriggerEffs(DATAnumeratorpt2OF, DATAdenominatorpt2OF)
    
    DATAeffeta1ee = getTriggerEffs(DATAnumeratoreta1ee, DATAdenominatoreta1ee)
    DATAeffeta1mm = getTriggerEffs(DATAnumeratoreta1mm, DATAdenominatoreta1mm)
    DATAeffeta1SF = getTriggerEffs(DATAnumeratoreta1SF, DATAdenominatoreta1SF)
    DATAeffeta1OF = getTriggerEffs(DATAnumeratoreta1OF, DATAdenominatoreta1OF)
    
    DATAeffeta2ee = getTriggerEffs(DATAnumeratoreta2ee, DATAdenominatoreta2ee)
    DATAeffeta2mm = getTriggerEffs(DATAnumeratoreta2mm, DATAdenominatoreta2mm)
    DATAeffeta2SF = getTriggerEffs(DATAnumeratoreta2SF, DATAdenominatoreta2SF)
    DATAeffeta2OF = getTriggerEffs(DATAnumeratoreta2OF, DATAdenominatoreta2OF)

    MCeffMETee =  getTriggerEffs(MCnumeratorMETee, MCdenominatorMETee)
    MCeffMETmm =  getTriggerEffs(MCnumeratorMETmm, MCdenominatorMETmm)  
    MCeffMETSF =  getTriggerEffs(MCnumeratorMETSF, MCdenominatorMETSF)  
    MCeffMETOF =  getTriggerEffs(MCnumeratorMETOF, MCdenominatorMETOF)
    MCeffmt2ee =  getTriggerEffs(MCnumeratormt2ee, MCdenominatormt2ee)
    MCeffmt2mm =  getTriggerEffs(MCnumeratormt2mm, MCdenominatormt2mm)
    MCeffmt2SF =  getTriggerEffs(MCnumeratormt2SF, MCdenominatormt2SF)
    MCeffmt2OF =  getTriggerEffs(MCnumeratormt2OF, MCdenominatormt2OF)
    MCeffMllee =  getTriggerEffs(MCnumeratorMllee, MCdenominatorMllee)
    MCeffMllmm =  getTriggerEffs(MCnumeratorMllmm, MCdenominatorMllmm)
    MCeffMllSF =  getTriggerEffs(MCnumeratorMllSF, MCdenominatorMllSF)
    MCeffMllOF =  getTriggerEffs(MCnumeratorMllOF, MCdenominatorMllOF)
    MCeffpt1ee =  getTriggerEffs(MCnumeratorpt1ee, MCdenominatorpt1ee)
    MCeffpt1mm =  getTriggerEffs(MCnumeratorpt1mm, MCdenominatorpt1mm)
    MCeffpt1SF =  getTriggerEffs(MCnumeratorpt1SF, MCdenominatorpt1SF)
    MCeffpt1OF =  getTriggerEffs(MCnumeratorpt1OF, MCdenominatorpt1OF)
    MCeffpt2ee =  getTriggerEffs(MCnumeratorpt2ee, MCdenominatorpt2ee)
    MCeffpt2mm =  getTriggerEffs(MCnumeratorpt2mm, MCdenominatorpt2mm)
    MCeffpt2SF =  getTriggerEffs(MCnumeratorpt2SF, MCdenominatorpt2SF)
    MCeffpt2OF =  getTriggerEffs(MCnumeratorpt2OF, MCdenominatorpt2OF)
    MCeffeta1ee = getTriggerEffs(MCnumeratoreta1ee, MCdenominatoreta1ee)
    MCeffeta1mm = getTriggerEffs(MCnumeratoreta1mm, MCdenominatoreta1mm)
    MCeffeta1SF = getTriggerEffs(MCnumeratoreta1SF, MCdenominatoreta1SF)
    MCeffeta1OF = getTriggerEffs(MCnumeratoreta1OF, MCdenominatoreta1OF)
    MCeffeta2ee = getTriggerEffs(MCnumeratoreta2ee, MCdenominatoreta2ee)
    MCeffeta2mm = getTriggerEffs(MCnumeratoreta2mm, MCdenominatoreta2mm)
    MCeffeta2SF = getTriggerEffs(MCnumeratoreta2SF, MCdenominatoreta2SF)
    MCeffeta2OF = getTriggerEffs(MCnumeratoreta2OF, MCdenominatoreta2OF)

    DATARTMET  = RT(DATAnumeratorMETSF, DATAeffMETee, DATAeffMETmm, DATAeffMETOF)
    DATARTmt2  = RT(DATAnumeratormt2SF, DATAeffmt2ee, DATAeffmt2mm, DATAeffmt2OF)
    DATARTMll  = RT(DATAnumeratorMllSF, DATAeffMllee, DATAeffMllmm, DATAeffMllOF)
    DATARTpt1  = RT(DATAnumeratorpt1SF, DATAeffpt1ee, DATAeffpt1mm, DATAeffpt1OF)
    DATARTpt2  = RT(DATAnumeratorpt2SF, DATAeffpt2ee, DATAeffpt2mm, DATAeffpt2OF)
    DATARTeta1 = RT(DATAnumeratoreta1SF, DATAeffeta1ee, DATAeffeta1mm, DATAeffeta1OF)
    DATARTeta2 = RT(DATAnumeratoreta2SF, DATAeffeta2ee, DATAeffeta2mm, DATAeffeta2OF)
    
    MCRTMET  = RT(MCnumeratorMETSF,  MCeffMETee,  MCeffMETmm,  MCeffMETOF)
    MCRTmt2  = RT(MCnumeratormt2SF,  MCeffmt2ee,  MCeffmt2mm,  MCeffmt2OF)
    MCRTMll  = RT(MCnumeratorMllSF,  MCeffMllee,  MCeffMllmm,  MCeffMllOF)
    MCRTpt1  = RT(MCnumeratorpt1SF,  MCeffpt1ee,  MCeffpt1mm,  MCeffpt1OF)
    MCRTpt2  = RT(MCnumeratorpt2SF,  MCeffpt2ee,  MCeffpt2mm,  MCeffpt2OF)
    MCRTeta1 = RT(MCnumeratoreta1SF, MCeffeta1ee, MCeffeta1mm, MCeffeta1OF)
    MCRTeta2 = RT(MCnumeratoreta2SF, MCeffeta2ee, MCeffeta2mm, MCeffeta2OF)

    effRTMll = Canvas.Canvas('rt/%s/plot_rt_mll'%(lumi_str), 'png,pdf',  0.6, 0.3, 0.8, 0.5)
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
    
    effRTMET = Canvas.Canvas('rt/%s/plot_rt_met'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
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

    
    effRTmt2 = Canvas.Canvas('rt/%s/plot_rt_mt2'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
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

    effRTpt1 = Canvas.Canvas('rt/%s/plot_rt_pt1'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
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

    effRTpt2 = Canvas.Canvas('rt/%s/plot_rt_pt2'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
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

    effRTeta1 = Canvas.Canvas('rt/%s/plot_rt_eta1'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
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

    effRTeta2 = Canvas.Canvas('rt/%s/plot_rt_eta2'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
    h_auxrteta2 = r.TH1F("h_auxRTeta2", "", 1, -2.4, 2.4)
    h_auxrteta2.GetYaxis().SetRangeUser(0.8, 1.2)
    h_auxrteta2.GetXaxis().SetRangeUser(-2.4, 2.4)
    h_auxrteta2.GetXaxis().SetTitle(labeleta2)
    effRTeta2.addLine(h_auxrteta2.GetXaxis().GetXmin(), dart, h_auxrteta2.GetXaxis().GetXmax(), dart ,r.kBlue)
    effRTeta2.addHisto(h_auxrteta2, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
    effRTeta2.addHisto(DATARTeta2, 'PE,SAME', 'RT', 'PL', r.kBlack , 1, 0)
    effRTeta2.addBand(h_auxrteta2.GetXaxis().GetXmin(), dart-dasystrt, h_auxrteta2.GetXaxis().GetXmax(), dart+dasystrt, r.kOrange+6, 0.2)
    effRTeta2.addLatex (0.6, 0.25, 'Mean R_{T} data: %.3f '%(dart))
    effRTeta2.addLatex (0.6, 0.2, 'Mean R_{T} MC  : %.3f '%(mcrt))
    effRTeta2.save(1, 1, 0, lumi, 0.8, 1.2)
   
    effMll = Canvas.Canvas('rt/%s/plot_eff_mll'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
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

    effmt2 = Canvas.Canvas('rt/%s/plot_eff_mt2'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
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

    effpt1 = Canvas.Canvas('rt/%s/plot_eff_pt1'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
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

    effpt2 = Canvas.Canvas('rt/%s/plot_eff_pt2'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
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
 
    effeta1 = Canvas.Canvas('rt/%s/plot_eff_eta1'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
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

    effeta2 = Canvas.Canvas('rt/%s/plot_eff_eta2'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
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
    saveInFile(theFile, 0, 0, 0, dart, dauncrt, dasystrt, "DATA")
    saveInFile(theFile, 0, 0, 0, mcrt, mcuncrt, mcsystrt, "MC")


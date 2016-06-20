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
import math,sys,optparse
import math



import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Rounder    as rounder


############################################################
def getRT(eff_ee, unc_ee, eff_mm, unc_mm, eff_em, unc_em):
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


def RT(eff_ee, eff_mm, eff_em):

    RTratio = eff_ee.Clone("RT_" + eff_ee.GetName())

    RTratio.GetYaxis().SetTitle("RT")
    RTratio.GetXaxis().SetTitle(eff_mm.GetXaxis().GetTitle())

    for i in range(0, eff_mm.GetNbinsX()+1):

        if(eff_em.GetBinContent(i) != 0 and eff_ee.GetBinContent(i) != 0 and eff_mm.GetBinContent(i) != 0):
            [rt, uncrt, uncsys] = getRT(eff_ee.GetBinContent(i), eff_ee.GetBinError(i), eff_mm.GetBinContent(i), eff_mm.GetBinError(i), eff_em.GetBinContent(i), eff_em.GetBinError(i)) 
            RTratio.SetBinContent(i, rt)
            RTratio.SetBinError(i, math.sqrt(uncrt*uncrt+uncsys*uncsys))

    return RTratio


def getTriggerEffs(num, den):

    trig = den.Clone("trig_" + den.GetName())
    trig.Reset()



    errs = TGraphAsymmErrors(num, den, 'v')
    x = errs.GetX()
    y = errs.GetY()
    eyh = errs.GetEYhigh()
    eyl = errs.GetEYlow()
    for i in range(0, num.GetNbinsX()+1):
        if(x[i] != 0 and y[i] != 0):
            trig.SetBinContent(i, y[i])
            trig.SetBinError(i, max(eyh[i], eyl[i]))

    return trig


if __name__ == '__main__':



##Main body of the analysis
if __name__ == '__main__':



    print bcolors.HEADER
    print '#######################################################################'
    print '                  Starting RT analysis...                          '
    print '#######################################################################' + bcolors.ENDC

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()

    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    mcDatasets = ['TTLep_pow']
    daDatasets = ['DoubleMuon_Run2016B_PromptReco_v2_runs_273150_273730', 'DoubleEG_Run2016B_PromptReco_v2_runs_273150_273730', 'MuonEG_Run2016B_PromptReco_v2_runs_273150_273730']

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    lumi = 2.1 ; maxrun = 999999
    lumi_str = '2.1invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()



    DATAdenominatorMllee =   treeMC.getTH1F(lumi, "DATAdenominatoree", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelx)
    DATAnumeratorMllee =     treeMC.getTH1F(lumi, "DATAnumeratoree", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.trigger, cuts.ee]), '', labelx)
    DATAdenominatorMllmm =   treeMC.getTH1F(lumi, "DATAdenominatormm", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelx)
    DATAnumeratorMllmm =     treeMC.getTH1F(lumi, "DATAnumeratormm", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.trigger, cuts.mm]), '', labelx)
    DATAdenominatorMllSF =   treeMC.getTH1F(lumi, "DATAdenominatorSF", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelx)
    DATAnumeratorMllSF =     treeMC.getTH1F(lumi, "DATAnumeratorSF", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.trigger, cuts.SF]), '', labelx)
    DATAdenominatorMllOF =   treeMC.getTH1F(lumi, "DATAdenominatorOF", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelx)
    DATAnumeratorMllOF =     treeMC.getTH1F(lumi, "DATAnumeratorOF", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.trigger, cuts.OF]), '', labelx)

    DATAdenominatorMlleevalue =   treeMC.getTH1F(lumi, "DATAdenominatoreevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelx)
    DATAnumeratorMlleevalue =     treeMC.getTH1F(lumi, "DATAnumeratoreevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.trigger, cuts.ee]), '', labelx)
    DATAdenominatorMllmmvalue =   treeMC.getTH1F(lumi, "DATAdenominatormmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelx)
    DATAnumeratorMllmmvalue =     treeMC.getTH1F(lumi, "DATAnumeratormmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.trigger, cuts.mm]), '', labelx)
    DATAdenominatorMllSFvalue =   treeMC.getTH1F(lumi, "DATAdenominatorSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelx)
    DATAnumeratorMllSFvalue =     treeMC.getTH1F(lumi, "DATAnumeratorSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.trigger, cuts.SF]), '', labelx)
    DATAdenominatorMllOFvalue =   treeMC.getTH1F(lumi, "DATAdenominatorOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelx)
    DATAnumeratorMllOFvalue =     treeMC.getTH1F(lumi, "DATAnumeratorOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.trigger, cuts.OF]), '', labelx)

    DATAdenominatorptee =   treeMC.getTH1F(lumi, "DATAdenominatoree", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelx)
    DATAnumeratorptee =     treeMC.getTH1F(lumi, "DATAnumeratoree", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.trigger, cuts.ee]), '', labelx)
    DATAdenominatorptmm =   treeMC.getTH1F(lumi, "DATAdenominatormm", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelx)
    DATAnumeratorptmm =     treeMC.getTH1F(lumi, "DATAnumeratormm", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.trigger, cuts.mm]), '', labelx)
    DATAdenominatorptSF =   treeMC.getTH1F(lumi, "DATAdenominatorSF", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelx)
    DATAnumeratorptSF =     treeMC.getTH1F(lumi, "DATAnumeratorSF", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.trigger, cuts.SF]), '', labelx)
    DATAdenominatorOF =   treeMC.getTH1F(lumi, "DATAdenominatorOF", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelx)
    DATAnumeratorMllOF =     treeMC.getTH1F(lumi, "DATAnumeratorOF", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.trigger, cuts.OF]), '', labelx)

 





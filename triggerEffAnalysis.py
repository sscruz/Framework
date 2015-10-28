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
import Rounder as rounder
import math

import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample




def getRT(self, eff_ee, unc_ee, eff_mm, unc_mm, eff_em, unc_em):
    RT = math.sqrt(eff_ee*eff_mm)/eff_em
    uncRTee = (math.sqrt(eff_mm)/eff_em)*0.5*(unc_ee/math.sqrt(eff_ee)) 
    uncRTmm = (math.sqrt(eff_ee)/eff_em)*0.5*(unc_mm/math.sqrt(eff_mm)) 
    uncRTem = math.sqrt(eff_ee*eff_mm)*((unc_em)/(eff_em*eff_em])) 
    uncRT = math.sqrt(uncRTee * uncRTee + uncRTmm * uncRTmm + uncRTem * uncRTem)
    uncsysRTee = (math.sqrt(eff_mm)/eff_em)*0.5*(0.05*eff_ee/math.sqrt(eff_ee)) 
    uncsysRTmm = (math.sqrt(eff_ee)/eff_em)*0.5*(0.05*eff_mm/math.sqrt(eff_mm)) 
    uncsysRTem = math.sqrt(eff_ee*eff_mm)*((0.05*eff_em)/(eff_em*eff_em])) 
    uncsysRT = math.sqrt(uncsysRTee * uncsysRTee + uncsysRTmm * uncsysRTmm + uncsysRTem * uncsysRTem)

    return [RT, uncRT, uncsysRT]


def getTriggerEffs(tree, numcut, dencut, var, varname, binning, lumi):


    passHisto = tree.getTH1F(lumi, 'passHisto', var, binning[0], binning[1], binning[2], numcut, '', varname)
    allHisto  = tree.getTH1F(lumi, 'allHisto',  var, binning[0], binning[1], binning[2], dencut, '', varname)
    passYields = tree.getYields(lumi, 't.lepsMll_Edge', 0, 1000, numcut)
    allYields = tree.getYields(lumi, 't.lepsMll_Edge', 0, 1000, dencut)
    eff = passYields[0]/allYields[0]
    unc = math.sqrt((passYields[1]*passYields[1])/(allYields[0]*allYields[0]) + (passYields[0]*passYields[0]*allYields[1]*allYields[1])/(allYields[0]*allYields[0]*allYields[0]*allYields[0]))    

    errs = TGraphAsymmErrors(passHisto, allHisto, 'a')

    errs.GetHistogram().GetXaxis().SetTitle(varname)

    return [errs, eff, unc]


if __name__ == '__main__':

    parser = optparse.OptionParser(usage='usage: %prog [options] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    (options, args) = parser.parse_args()

    if len(args) != 1:
      parser.error('wrong number of arguments')

    inputFileName = args[0]
    print 'Going to load DATA and MC trees...'
    mcDatasets = ['TTLep_pow', 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    daDatasets = ['HTMHT_Run2015D_05Oct_v1_runs_246908_258751', 'HTMHT_Run2015D_v4_runs_246908_258751',
                  'JetHT_Run2015D_05Oct_v1_runs_246908_258751', 'JetHT_Run2015D_v4_runs_246908_258751']
    treeMC = Sample.Tree(helper.selectSamples(inputFileName, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(inputFileName, daDatasets, 'DA'), 'DATA', 1)
    #tree = treeMC
    print 'Trees successfully loaded...'

    lumi = 1.28
    lumi_str = 'lumi'+str(lumi).replace('.', 'p')
   
    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()
 
    vetoSignalCR = "(!(" + cuts.METJetsSignalRegion + ")&& !(" + cuts.METJetsControlRegion + "))"
    denominator_EE_c = cuts.AddList([cuts.GoodLeptonNoTriggeree(), cuts.central, cuts.triggerHT, cuts.HT, vetoSignalCR])
    denominator_MM_c = cuts.AddList([cuts.GoodLeptonNoTriggermm(), cuts.central, cuts.triggerHT, cuts.HT, vetoSignalCR])
    denominator_EM_c = cuts.AddList([cuts.GoodLeptonNoTriggerOF(), cuts.central, cuts.triggerHT, cuts.HT, vetoSignalCR])
    numerator_EE_c = cuts.AddList([denominator_EE, cuts.trigEEc])
    numerator_MM_c = cuts.AddList([denominator_MM, cuts.trigMMc])
    numerator_EM_c = cuts.AddList([denominator_EM, cuts.trigEMc])
    denominator_EE_f = cuts.AddList([cuts.GoodLeptonNoTriggeree(), cuts.forward, cuts.triggerHT, cuts.HT, vetoSignalCR])
    denominator_MM_f = cuts.AddList([cuts.GoodLeptonNoTriggermm(), cuts.forward, cuts.triggerHT, cuts.HT, vetoSignalCR])
    denominator_EM_f = cuts.AddList([cuts.GoodLeptonNoTriggerOF(), cuts.forward, cuts.triggerHT, cuts.HT, vetoSignalCR])
    numerator_EE_f = cuts.AddList([denominator_EE, cuts.trigEEc])
    numerator_MM_f = cuts.AddList([denominator_MM, cuts.trigMMc])
    numerator_EM_f = cuts.AddList([denominator_EM, cuts.trigEMc])


    for dType in ['MC', 'DATA']:
        [eff_mll_ee_ce, eff_ee_ce, unc_ee_ce] = getTriggerEffs(treeDA, numerator_EE_c, denominator_EE_c, 't.lepsMll_Edge', 'm_{ee} (GeV)', [20,  0, 200], lumi)
        [eff_pt1_ee_ce, eff_ee_ce, unc_ee_ce] = getTriggerEffs(treeDA, numerator_EE_c, denominator_EE_c, 't.Lep1_pt_Edge', 'p_{T}^{e,lead} (GeV)', [10,  20, 120], lumi)
        [eff_pt2_ee_ce, eff_ee_ce, unc_ee_ce] = getTriggerEffs(treeDA, numerator_EE_c, denominator_EE_c, 't.Lep2_pt_Edge', 'p_{T}^{e,trail} (GeV)', [10,  20, 120], lumi)
        [eff_mll_ee_fo, eff_ee_ce, unc_ee_fo] = getTriggerEffs(treeDA, numerator_EE_f, denominator_EE_f, 't.lepsMll_Edge', 'm_{ee} (GeV)', [20,  0, 200], lumi)
        [eff_pt1_ee_fo, eff_ee_ce, unc_ee_fo] = getTriggerEffs(treeDA, numerator_EE_f, denominator_EE_f, 't.Lep1_pt_Edge', 'p_{T}^{e,lead} (GeV)', [10,  20, 120], lumi)
        [eff_pt2_ee_fo, eff_ee_ce, unc_ee_fo] = getTriggerEffs(treeDA, numerator_EE_f, denominator_EE_f, 't.Lep2_pt_Edge', 'p_{T}^{e,trail} (GeV)', [10,  20, 120], lumi)
    
        [eff_mll_mm_ce, eff_mm_ce, unc_mm_ce] = getTriggerEffs(trmmDA, numerator_MM_c, denominator_MM_c, 't.lepsMll_Edge', 'm_{mm} (GeV)', [20,  0, 200], lumi)
        [eff_pt1_mm_ce, eff_mm_ce, unc_mm_ce] = getTriggerEffs(trmmDA, numerator_MM_c, denominator_MM_c, 't.Lep1_pt_Edge', 'p_{T}^{m,lead} (GeV)', [10,  20, 120], lumi)
        [eff_pt2_mm_ce, eff_mm_ce, unc_mm_ce] = getTriggerEffs(trmmDA, numerator_MM_c, denominator_MM_c, 't.Lep2_pt_Edge', 'p_{T}^{m,trail} (GeV)', [10,  20, 120], lumi)
        [eff_mll_mm_fo, eff_mm_ce, unc_mm_fo] = getTriggerEffs(trmmDA, numerator_MM_f, denominator_MM_f, 't.lepsMll_Edge', 'm_{mm} (GeV)', [20,  0, 200], lumi)
        [eff_pt1_mm_fo, eff_mm_ce, unc_mm_fo] = getTriggerEffs(trmmDA, numerator_MM_f, denominator_MM_f, 't.Lep1_pt_Edge', 'p_{T}^{m,lead} (GeV)', [10,  20, 120], lumi)
        [eff_pt2_mm_fo, eff_mm_ce, unc_mm_fo] = getTriggerEffs(trmmDA, numerator_MM_f, denominator_MM_f, 't.Lep2_pt_Edge', 'p_{T}^{m,trail} (GeV)', [10,  20, 120], lumi)

        [eff_mll_em_ce, eff_em_ce, unc_em_ce] = getTriggerEffs(tremDA, numerator_EM_c, denominator_EM_c, 't.lepsMll_Edge', 'm_{em} (GeV)', [20,  0, 200], lumi)
        [eff_pt1_em_ce, eff_em_ce, unc_em_ce] = getTriggerEffs(tremDA, numerator_EM_c, denominator_EM_c, 't.Lep1_pt_Edge', 'p_{T}^{lead} (GeV)', [10,  20, 120], lumi)
        [eff_pt2_em_ce, eff_em_ce, unc_em_ce] = getTriggerEffs(tremDA, numerator_EM_c, denominator_EM_c, 't.Lep2_pt_Edge', 'p_{T}^{trail} (GeV)', [10,  20, 120], lumi)
        [eff_mll_em_fo, eff_em_ce, unc_em_fo] = getTriggerEffs(tremDA, numerator_EM_f, denominator_EM_f, 't.lepsMll_Edge', 'm_{em} (GeV)', [20,  0, 200], lumi)
        [eff_pt1_em_fo, eff_em_ce, unc_em_fo] = getTriggerEffs(tremDA, numerator_EM_f, denominator_EM_f, 't.Lep1_pt_Edge', 'p_{T}^{lead} (GeV)', [10,  20, 120], lumi)
        [eff_pt2_em_fo, eff_em_ce, unc_em_fo] = getTriggerEffs(tremDA, numerator_EM_f, denominator_EM_f, 't.Lep2_pt_Edge', 'p_{T}^{trail} (GeV)', [10,  20, 120], lumi)
    
        [RT_c, uncRT_c, sysRT_c] = getRT(eff_ee_ce, unc_ee_ce, eff_mm_ce, unc_mm_ce, eff_em_ce, unc_em_ce)
        [RT_f, uncRT_f, sysRT_f] = getRT(eff_ee_fo, unc_ee_fo, eff_mm_fo, unc_mm_fo, eff_em_fo, unc_em_fo)


        a = rounder.Rounder()
        print effs 
        print "------------ Summary of values for", dType, "---------------"
        print "EE efficiency central", a.toStringB(eff_ee_ce, unc_ee_ce)
        print "MM efficiency central", a.toStringB(eff_mm_ce, unc_mm_ce)
        print "EM efficiency central", a.toStringB(eff_em_ce, unc_em_ce)
        print "EE efficiency forward", a.toStringB(eff_ee_fo, unc_ee_fo)
        print "MM efficiency forward", a.toStringB(eff_mm_fo, unc_mm_fo)
        print "EM efficiency forward", a.toStringB(eff_em_fo, unc_em_fo)
        print "RT central", a.toStringB(RT_c, uncRT_c), "+/-", a.toString(sysRT_c)
        print "RT central", a.toStringB(RT_f, uncRT_f), "+/-", a.toString(sysRT_f)

        ing = ingredients("ingredients.dat", dType)
        ing.UpdateRTValues([RT_c, uncRT_c, sysRT_c], [RT_f, uncRT_f, sysRT_f])







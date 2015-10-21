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


def getTriggerEffs(tree, numcut, dencut, var, varname, binning, lumi):


    passHisto = tree.getTH1F(lumi, 'passHisto', var, binning[0], binning[1], binning[2], numcut, '', varname)
    allHisto  = tree.getTH1F(lumi, 'allHisto',  var, binning[0], binning[1], binning[2], dencut, '', varname)
    passYields = tree.getYields(lumi, 't.lepsMll_Edge', 0, 1000, numcut)
    allYields = tree.getYields(lumi, 't.lepsMll_Edge', 0, 1000, dencut)
    eff = passYields[0]/allYields[0]
    unc = math.sqrt((passYields[1]*passYields[1])/(allYields[0]*allYields[0]) + (passYields[0]*passYields[0]*allYields[1]*allYields[1])/(allYields[0]*allYields[0]*allYields[0]*allYields[0]))    


    for i in range(1,passHisto.GetNbinsX()+1):
        print 'at variable %s events passing/total %.2f  of  %.2f' %(varname, passHisto.GetBinContent(i), allHisto.GetBinContent(i) )


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

    denominator_EE = cuts.AddList([cuts.GoodLeptonNoTriggeree(), cuts.triggerHT, cuts.HT])
    denominator_MM = cuts.AddList([cuts.GoodLeptonNoTriggermm(), cuts.triggerHT, cuts.HT])
    denominator_EM = cuts.AddList([cuts.GoodLeptonNoTriggerOF(), cuts.triggerHT, cuts.HT])
    numerator_EE = cuts.AddList([denominator_EE, cuts.trigEEc])
    numerator_MM = cuts.AddList([denominator_MM, cuts.trigMMc])
    numerator_EM = cuts.AddList([denominator_EM, cuts.trigEMc])

    
    ##sys.exit(0)

    effs = []
    uncs = []
    for flavor in ['ee', 'mm', 'em']:

        if   flavor is 'ee': 
            trigger = cuts.trigEEc
            mllvar  = 'm_{ee} (GeV)'
            pt1var  = 'p_{T}^{e,lead} (GeV)'
            pt2var  = 'p_{T}^{e,trail} (GeV)'
            numcut = numerator_EE
            dencut = denominator_EE
            tag    = "EE"
        elif flavor is 'mm': 
            trigger = cuts.trigMMc
            mllvar  = 'm_{#mu#mu} (GeV)'
            pt1var  = 'p_{T}^{#mu,lead} (GeV)'
            pt2var  = 'p_{T}^{#mu,trail} (GeV)'
            numcut = numerator_MM
            dencut = denominator_MM
            tag    = "MM"
        elif flavor is 'em': 
            trigger = cuts.trigEMc
            mllvar  = 'm_{e#mu} (GeV)'
            pt1var  = 'p_{T}^{lep,lead} (GeV)'
            pt2var  = 'p_{T}^{lep,trail} (GeV)'
            numcut = numerator_EM
            dencut = denominator_EM
            tag    = "EM"
        else: 
            print 'something is wrong...'

        [eff_mll, eff, unc]  = getTriggerEffs(treeDA, numcut, dencut, 't.lepsMll_Edge', mllvar, [20,  0, 200], lumi)
        effs.append(eff)
        uncs.append(unc)
        [eff_l1pt, eff, unc] = getTriggerEffs(treeDA, numcut, dencut, 't.Lep1_pt_Edge', pt1var, [10, 20, 120], lumi)
        effs.append(eff)
        uncs.append(unc)
        [eff_l2pt, eff, unc] = getTriggerEffs(treeDA, numcut, dencut, 't.Lep2_pt_Edge', pt2var, [10, 20, 120], lumi)
        effs.append(eff)
        uncs.append(unc)

        plot_mll = Canvas.Canvas('trigger/%s/plot_eff_%s_mll'%(lumi_str, flavor), 'png', 0.6, 0.6, 0.8, 0.8)
        eff_mll.GetHistogram().Draw()
        eff_mll.GetYaxis().SetRangeUser(-0.05, 1.05)
        eff_mll.Draw('apz')
        plot_mll.save(0, 0, 0, lumi)

        plot_l1pt = Canvas.Canvas('trigger/%s/plot_eff_%s_l1pt'%(lumi_str, flavor), 'png', 0.6, 0.6, 0.8, 0.8)
        eff_l1pt.GetHistogram().Draw()
        eff_l1pt.GetYaxis().SetRangeUser(-0.05, 1.05)
        eff_l1pt.Draw('apz')
        plot_l1pt.save(0, 0, 0, lumi)

        plot_l2pt = Canvas.Canvas('trigger/%s/plot_eff_%s_l2pt'%(lumi_str, flavor), 'png', 0.6, 0.6, 0.8, 0.8)
        eff_l2pt.GetHistogram().Draw()
        eff_l2pt.GetYaxis().SetRangeUser(-0.05, 1.05)
        eff_l2pt.Draw('apz')
        plot_l2pt.save(0, 0, 0, lumi)

        del plot_mll, plot_l1pt, plot_l2pt


    a = rounder.Rounder()
 
    print "Summary of values"
    print "EE efficiency ", a.toStringB(effs[0][0], effs[0][1])
    print "MM efficiency ", a.toStringB(effs[1][0], effs[1][1])
    print "EM efficiency ", a.toStringB(effs[2][0], effs[2][1])

    RT = math.sqrt(effs[0][0]*effs[1][0])/effs[2][0]
    uncRTee = (math.sqrt(effs[1][0])/effs[2][0])*0.5*(effs[0][1]/math.sqrt(effs[0][0])) 
    uncRTmm = (math.sqrt(effs[0][0])/effs[2][0])*0.5*(effs[1][1]/math.sqrt(effs[1][0])) 
    uncRTem = math.sqrt(effs[0][0]*effs[1][0])*(math.sqrt(effs[2][1])/(effs[2][0]*effs[2][0])) 
    uncRT = math.sqrt(uncRTee * uncRTee + uncRTmm * uncRTmm + uncRTem * uncRTem)










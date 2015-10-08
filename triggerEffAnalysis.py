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

import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample


def getTriggerEffs(tree, cut, addCut, var, varname, binning, lumi):


    passHisto = tree.getTH1F(lumi, 'passHisto', var, binning[0], binning[1], binning[2], cut[:-1]+'&&'+addCut+')', '', varname)
    allHisto  = tree.getTH1F(lumi, 'allHisto',  var, binning[0], binning[1], binning[2], cut                     , '', varname)

    for i in range(1,passHisto.GetNbinsX()+1):
        print 'at variable %s events passing/total %.2f  of  %.2f' %(varname, passHisto.GetBinContent(i), allHisto.GetBinContent(i) )

    ##ratio = passHisto.Clone('eff_' + passHisto.GetName())
    ##ratio.Divide(allHisto)
    ##ratio.GetYaxis().SetTitle('trigger eff.')

    errs = TGraphAsymmErrors(passHisto, allHisto, 'a')

    errs.GetHistogram().GetXaxis().SetTitle(varname)

    return errs


if __name__ == '__main__':

    parser = optparse.OptionParser(usage='usage: %prog [options] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    (options, args) = parser.parse_args()

    if len(args) != 1:
      parser.error('wrong number of arguments')

    inputFileName = args[0]
    print 'Going to load DATA and MC trees...'
    mcDatasets = ['TTJets', 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    daDatasets = ['HTMHT_Run2015C', 'JetHT_Run2015C',
                  'HTMHT_Run2015D', 'JetHT_Run2015D']
    #daDatasets = ['DoubleMuon_Run2015C', 'DoubleEG_Run2015C', 'MuonEG_Run2015C',
    #              'DoubleMuon_Run2015D', 'DoubleEG_Run2015D', 'MuonEG_Run2015D']
    treeMC = Sample.Tree(helper.selectSamples(inputFileName, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(inputFileName, daDatasets, 'DA'), 'DATA', 1)
    #tree = treeMC
    print 'Trees successfully loaded...'

    lumi = 0.225
   
    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()
    
    ht_ht = getTriggerEffs(treeDA, cuts.AddList([cuts.goodLepton,cuts.Central()]), cuts.triggerHT, 't.htJet35j_Edge', 'H_{T} (GeV)', [20, 0, 800], lumi)
   
    plot_ht_ht = Canvas.Canvas('trigger/plot_eff_ht_ht', 'png', 0.6, 0.6, 0.8, 0.8)
    ht_ht.GetHistogram().Draw()
    ht_ht.GetYaxis().SetRangeUser(-0.05, 1.05)
    ht_ht.Draw('apz')
    plot_ht_ht.save(0, 0, 0, lumi)
    
    ##sys.exit(0)

    for flavor in ['ee', 'mm', 'em']:

        if   flavor is 'ee': 
            trigger = cuts.trigEEc
            mllvar  = 'm_{ee} (GeV)'
            pt1var  = 'p_{T}^{e,lead} (GeV)'
            pt2var  = 'p_{T}^{e,trail} (GeV)'
            lepcut = cuts.GoodLeptonee()
        elif flavor is 'mm': 
            trigger = cuts.trigMMc
            mllvar  = 'm_{#mu#mu} (GeV)'
            pt1var  = 'p_{T}^{#mu,lead} (GeV)'
            pt2var  = 'p_{T}^{#mu,trail} (GeV)'
            lepcut = cuts.GoodLeptonmm()
        elif flavor is 'em': 
            trigger = cuts.trigEMc
            mllvar  = 'm_{e#mu} (GeV)'
            pt1var  = 'p_{T}^{lep,lead} (GeV)'
            pt2var  = 'p_{T}^{lep,trail} (GeV)'
            lepcut = cuts.GoodLeptonOF()
        else: 
            print 'something is wrong...'

        eff_mll  = getTriggerEffs(treeDA, cuts.AddList([lepcut,cuts.Central(),cuts.triggerHT, 't.htJet35j_Edge > 400']), trigger, 't.lepsMll_Edge', mllvar, [20,  0, 200], lumi)
        eff_l1pt = getTriggerEffs(treeDA, cuts.AddList([lepcut,cuts.Central(),cuts.triggerHT, 't.htJet35j_Edge > 400']), trigger, 't.Lep1_pt_Edge', pt1var, [10, 20, 120], lumi)
        eff_l2pt = getTriggerEffs(treeDA, cuts.AddList([lepcut,cuts.Central(),cuts.triggerHT, 't.htJet35j_Edge > 400']), trigger, 't.Lep2_pt_Edge', pt2var, [10, 20, 120], lumi)

        plot_mll = Canvas.Canvas('trigger/plot_eff_'+flavor+'_mll', 'png', 0.6, 0.6, 0.8, 0.8)
        eff_mll.GetHistogram().Draw()
        eff_mll.GetYaxis().SetRangeUser(-0.05, 1.05)
        eff_mll.Draw('apz')
        plot_mll.save(0, 0, 0, lumi)

        plot_l1pt = Canvas.Canvas('trigger/plot_eff_'+flavor+'_l1pt', 'png', 0.6, 0.6, 0.8, 0.8)
        eff_l1pt.GetHistogram().Draw()
        eff_l1pt.GetYaxis().SetRangeUser(-0.05, 1.05)
        eff_l1pt.Draw('apz')
        plot_l1pt.save(0, 0, 0, lumi)

        plot_l2pt = Canvas.Canvas('trigger/plot_eff_'+flavor+'_l2pt', 'png', 0.6, 0.6, 0.8, 0.8)
        eff_l2pt.GetHistogram().Draw()
        eff_l2pt.GetYaxis().SetRangeUser(-0.05, 1.05)
        eff_l2pt.Draw('apz')
        plot_l2pt.save(0, 0, 0, lumi)

        del plot_mll, plot_l1pt, plot_l2pt


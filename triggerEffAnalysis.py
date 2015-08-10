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
import math as math
import sys

from optparse import OptionParser
from ROOT import gROOT, TCanvas, TFile, TF1, TGraphAsymmErrors
from Sample import Sample, Block, Tree
from CutManager import CutManager
from Canvas import Canvas


def getTriggerEffs(tree, cut, addCut, var, varname, binning, lumi):


    passHisto = tree.getTH1F(lumi, 'passHisto', var, binning[0], binning[1], binning[2], cut[:-1]+'&&'+addCut+')', '', varname)
    allHisto  = tree.getTH1F(lumi, 'allHisto',  var, binning[0], binning[1], binning[2], cut                     , '', varname)

    for i in range(1,passHisto.GetNbinsX()+1):
        print 'at variable %s events passing/total %.2f  of  %.2f' %(varname, passHisto.GetBinContent(i), allHisto.GetBinContent(i) )

    ratio = passHisto.Clone('eff_' + passHisto.GetName())
    ratio.Divide(allHisto)
    ratio.GetYaxis().SetTitle('trigger eff.')

    errs = TGraphAsymmErrors(passHisto, allHisto, 'a')

    errs.GetHistogram().GetXaxis().SetTitle(varname)

    return errs


if __name__ == '__main__':

    parser = OptionParser(usage='usage: %prog [options] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    (options, args) = parser.parse_args()

    if len(args) != 1:
      parser.error('wrong number of arguments')

    inputFileName = args[0]
    tree = Tree(inputFileName, 'MC', 0)
   
    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager()
    
    ht_ht = getTriggerEffs(tree, cuts.AddList([cuts.goodLepton,cuts.Central()]), 'HLT_HT > 0', 't.htJet35j_Edge', 'H_{T} (GeV)', [20,   0,  800], 4.) ## ,'(HLT_DoubleMu > 0 || HLT_DoubleEl > 0|| HLT_MuEG > 0)'
   
    plot_ht_ht = Canvas('plot_eff_ht_ht', 'png', 0.6, 0.6, 0.8, 0.8)
    ht_ht.GetHistogram().Draw()
    ht_ht.Draw('apz')
    plot_ht_ht.save(0, 0, 0, 4.0)
    
    #sys.exit(0)

    for flavor in ['ee', 'mm', 'em']:

        if   flavor is 'ee': 
            trigger = 'HLT_DoubleEl'
            mllvar  = 'm_{ee} (GeV)'
            pt1var  = 'p_{T}^{e,lead} (GeV)'
            pt2var  = 'p_{T}^{e,trail} (GeV)'
            lepcut = cuts.GoodLeptonee()
        elif flavor is 'mm': 
            trigger = 'HLT_DoubleMu'
            mllvar  = 'm_{#mu#mu} (GeV)'
            pt1var  = 'p_{T}^{#mu,lead} (GeV)'
            pt2var  = 'p_{T}^{#mu,trail} (GeV)'
            lepcut = cuts.GoodLeptonmm()
        elif flavor is 'em': 
            trigger = 'HLT_MuEG'
            mllvar  = 'm_{e#mu} (GeV)'
            pt1var  = 'p_{T}^{lep,lead} (GeV)'
            pt2var  = 'p_{T}^{lep,trail} (GeV)'
            lepcut = cuts.GoodLeptonOF()
        else: 
            print 'something is wrong...'

        eff_mll  = getTriggerEffs(tree, cuts.AddList([lepcut,cuts.Central(),'HLT_HT > 0', 't.htJet35j_Edge > 400']), trigger, 't.lepsMll_Edge'  , mllvar, [20,  0, 200], 4.)
        eff_l1pt = getTriggerEffs(tree, cuts.AddList([lepcut,cuts.Central(),'HLT_HT > 0', 't.htJet35j_Edge > 400']), trigger, 't.Lep_Edge_pt[0]', pt1var, [10, 20, 120], 4.)
        eff_l2pt = getTriggerEffs(tree, cuts.AddList([lepcut,cuts.Central(),'HLT_HT > 0', 't.htJet35j_Edge > 400']), trigger, 't.Lep_Edge_pt[1]', pt2var, [10, 20, 120], 4.)

        plot_mll = Canvas('plot_eff_'+flavor+'_mll', 'png', 0.6, 0.6, 0.8, 0.8)
        eff_mll.GetHistogram().Draw()
        eff_mll.Draw('apz')
        plot_mll.save(0, 0, 0, 4.0)

        plot_l1pt = Canvas('plot_eff_'+flavor+'_l1pt', 'png', 0.6, 0.6, 0.8, 0.8)
        eff_l1pt.GetHistogram().Draw()
        eff_l1pt.Draw('apz')
        plot_l1pt.save(0, 0, 0, 4.0)

        plot_l2pt = Canvas('plot_eff_'+flavor+'_l2pt', 'png', 0.6, 0.6, 0.8, 0.8)
        eff_l2pt.GetHistogram().Draw()
        eff_l2pt.Draw('apz')
        plot_l2pt.save(0, 0, 0, 4.0)

        del plot_mll, plot_l1pt, plot_l2pt


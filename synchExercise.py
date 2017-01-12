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
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats
import math, sys, optparse, os, subprocess

import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.nll



if __name__ == '__main__':

    parser = optparse.OptionParser(usage='usage: %prog [options] FilenameWithSamples', version='%prog 1.0')
    ##parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    (options, args) = parser.parse_args()


    if len(args) != 1:
      parser.error('wrong number of arguments')

    inputFileName = args[0]

    print 'Going to load DATA and MC trees...'

    synchDatasets = ['DoubleEG_Run2016*_23Sep2016_v*',
                     'DoubleMuon_Run2016*_23Sep2016_v*',
                     'MuonEG_Run2016*_23Sep2016_v*']

    treeSynch = Sample.Tree(helper.selectSamples(inputFileName, synchDatasets, 'SYNCH'), 'SYNCH'  , 1)

    print 'Trees successfully loaded...'

    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()

    outputName = 'synch_baseline_nonttbarlike_allMll_20170112'

    ## define the cut you want to have applied
    cutString = cuts.AddList([cuts.protection, cuts.Baseline, cuts.Zveto, cuts.OF, cuts.EdgeBaseline, cuts.ttBarLike])
    print 'producing event list for the following cutstring:\r==================================================='
    print cutString
    print '==================================================='


    for ind,sam in enumerate(treeSynch.blocks[0].samples):
        actualTree = treeSynch.blocks[0].samples[ind].ttree

        ## define the variables you want to have scanned
        scanString  = "run_Edge:lumi_Edge:evt_Edge:"
        scanString += "Lep1_pt_Edge:Lep1_eta_Edge:Lep1_phi_Edge:Lep1_pdgId_Edge:Lep1_miniRelIso_Edge:"
        scanString += "Lep2_pt_Edge:Lep2_eta_Edge:Lep2_phi_Edge:Lep2_pdgId_Edge:Lep2_miniRelIso_Edge:"
        scanString += "lepsMll_Edge:met_Edge:mt2_Edge:nJetSel_Edge:nBJetMedium35_Edge:abs(j1MetDPhi_Edge):abs(j2MetDPhi_Edge)"

        lenScan = len(scanString.split(':'))


        actualTree.SetScanField(-1)
        save = os.dup( sys.stdout.fileno() )
        if not ind:
            of_name = 'evtLists/{on}.txt'.format(on=outputName)
            newout = file( of_name, 'w' )
        os.dup2( newout.fileno(), sys.stdout.fileno() )
        actualTree.Scan(scanString, cutString, 'colsize=12')#colsize=8 precision=8 col=8d:12d');
        os.dup2( save, sys.stdout.fileno() )
        if ind == len(treeSynch.blocks[0].samples)-1:
            newout.close()

subprocess.call(['python']+['evtLists/overlap.py', '{of}'.format(of=of_name)] )#+['/dev/null'],stderr=subprocess.PIPE)

   

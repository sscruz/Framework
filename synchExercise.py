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
import math, sys, optparse, os

import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample



if __name__ == '__main__':

    parser = optparse.OptionParser(usage='usage: %prog [options] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    parser.add_option('-t', '--trigger', action='store', type='int', dest='triggerFlag', default='1', help='Trigger cut. Set to 0 if you want to run without trigger')
    (options, args) = parser.parse_args()


    if len(args) != 1:
      parser.error('wrong number of arguments')

    inputFileName = args[0]

    print 'Going to load DATA and MC trees...'
    synchDatasets = ['TTLep_pow']
    treeSynch = Sample.Tree(helper.selectSamples(inputFileName, synchDatasets, 'MC'), 'MC'  , 0)

    print 'Trees successfully loaded...'

    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()

    actualTree = treeSynch.blocks[0].samples[0].ttree

    ## define the variables you want to have scanned
    scanString  = "evt:"
    scanString += "Lep1_pt_Edge:Lep1_eta_Edge:Lep1_phi_Edge:Lep1_pdgId_Edge:Lep1_miniRelIso_Edge:Lep1_mvaIdSpring15_Edge:t.Lep1_minTauDR_Edge:"
    scanString += "Lep2_pt_Edge:Lep2_eta_Edge:Lep2_phi_Edge:Lep2_pdgId_Edge:Lep2_miniRelIso_Edge:Lep2_mvaIdSpring15_Edge:t.Lep2_minTauDR_Edge:"
    scanString += "lepsMll_Edge:met_pt:met_rawPt:nJetSel_Edge:nBJetMedium35_Edge:rhoCN:(Lep1_pdgId_Edge*Lep2_pdgId_Edge)"

    ## define the cut you want to have applied
    cutString = cuts.goodLepton

    actualTree.SetScanField(-1)
    save = os.dup( sys.stdout.fileno() )
    newout = file( 'synch_cerneth.txt', 'w' )
    os.dup2( newout.fileno(), sys.stdout.fileno() )
    actualTree.Scan(scanString, cutString, '')#colsize=8 precision=8 col=8d:12d');
    os.dup2( save, sys.stdout.fileno() )
    newout.close()

    

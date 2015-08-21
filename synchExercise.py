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
import sys, os

from optparse import OptionParser
from ROOT import gROOT, TCanvas, TFile, TF1, TGraphAsymmErrors
from Sample import Sample, Block, Tree
from CutManager import CutManager
from Canvas import Canvas


if __name__ == '__main__':

    parser = OptionParser(usage='usage: %prog [options] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    (options, args) = parser.parse_args()

    if len(args) != 1:
      parser.error('wrong number of arguments')

    inputFileName = args[0]
    tree = Tree(inputFileName, 'MC', 0)
   
    tree.printTree()

    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager()

    actualTree = tree.blocks[0].samples[0].ttree

    ## define the variables you want to have scanned
    scanString  = "run:evt:"
    scanString += "Lep1_pt_Edge:Lep1_eta_Edge:Lep1_phi_Edge:Lep1_pdgId_Edge:Lep1_mvaIdPhys14_Edge:"
    scanString += "Lep2_pt_Edge:Lep2_eta_Edge:Lep2_phi_Edge:Lep2_pdgId_Edge:Lep2_mvaIdPhys14_Edge:"
    scanString += "lepsMll_Edge:met_pt:nJetSel_Edge:(Lep1_pdgId_Edge*Lep2_pdgId_Edge)"

    ## define the cut you want to have applied
    cutString = cuts.goodLepton

    actualTree.SetScanField(-1)
    save = os.dup( sys.stdout.fileno() )
    newout = file( 'synch_cerneth.txt', 'w' )
    os.dup2( newout.fileno(), sys.stdout.fileno() )
    actualTree.Scan(scanString, cutString, 'colsize=6 precision=6 col=6d:12d');
    os.dup2( save, sys.stdout.fileno() )
    newout.close()

    

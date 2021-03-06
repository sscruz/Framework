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

    print "caca"
    parser = optparse.OptionParser(usage='usage: %prog [options] FilenameWithSamples', version='%prog 1.0')
    ##parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    (options, args) = parser.parse_args()


    if len(args) != 1:
      parser.error('wrong number of arguments')

    inputFileName = args[0]

    print 'Going to load DATA and MC trees...'
    #synchDatasets = ['TTLep_pow']
    ## 1p3 fb-1 synchDatasets = ['DoubleMuon_Run2015C_25ns_05Oct_v1_runs_246908_258714' , 'DoubleEG_Run2015C_25ns_05Oct_v1_runs_246908_258714' , 'MuonEG_Run2015C_25ns_05Oct_v1_runs_246908_258714' ,
    ## 1p3 fb-1                  'DoubleMuon_Run2015D_05Oct_v1_runs_246908_258751'      , 'DoubleEG_Run2015D_05Oct_v1_runs_246908_258751'      , 'MuonEG_Run2015D_05Oct_v2_runs_246908_258751'      ,
    ## 1p3 fb-1                  'DoubleMuon_Run2015D_v4_runs_246908_258751'            , 'DoubleEG_Run2015D_v4_runs_246908_258751'            , 'MuonEG_Run2015D_v4_runs_246908_258751'            ]
    #mcDatasets = ['TTLep_pow', 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    synchDatasets = ['HTMHT_Run2015C_25ns-05Oct_v1_runs_246908_260628', 'HTMHT_Run2015D-05Oct_v1_runs_246908_260628', 'HTMHT_Run2015D_v4_runs_246908_260628', 'JetHT_Run2015C_25ns-05Oct_v1_runs_246908_260628', 'JetHT_Run2015D-05Oct_v1_runs_246908_260628', 'JetHT_Run2015D_v4_runs_246908_260628']

    #synchDatasets = ['DoubleMuon_Run2015C_25ns-05Oct_v1_runs_246908_260628' , 'DoubleEG_Run2015C_25ns-05Oct_v1_runs_246908_260628' , 'MuonEG_Run2015C_25ns-05Oct_v1_runs_246908_260628' ,                     'DoubleMuon_Run2015D-05Oct_v1_runs_246908_260628'      , 'DoubleEG_Run2015D-05Oct_v1_runs_246908_260628'      , 'MuonEG_Run2015D-05Oct_v2_runs_246908_260628'      ,                      'DoubleMuon_Run2015D_v4_runs_246908_260628'            , 'DoubleEG_Run2015D_v4_runs_246908_260628'            , 'MuonEG_Run2015D_v4_runs_246908_260628'            ]
    treeSynch = Sample.Tree(helper.selectSamples(inputFileName, synchDatasets, 'SYNCH'), 'SYNCH'  , 1)

    print 'Trees successfully loaded...'

    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()

    for i in synchDatasets:
        ind = synchDatasets.index(i)
        actualTree = treeSynch.blocks[0].samples[ind].ttree

        ## define the variables you want to have scanned
        scanString  = "run:lumi:evt"
        #scanString += "Lep1_pt_Edge:Lep1_eta_Edge:Lep1_phi_Edge:Lep1_pdgId_Edge:Lep1_miniRelIso_Edge:Lep1_mvaIdSpring15_Edge:t.Lep1_minTauDR_Edge:"
        #scanString += "Lep2_pt_Edge:Lep2_eta_Edge:Lep2_phi_Edge:Lep2_pdgId_Edge:Lep2_miniRelIso_Edge:Lep2_mvaIdSpring15_Edge:t.Lep2_minTauDR_Edge:"
        #scanString += "lepsMll_Edge:met_pt:met_rawPt:nJetSel_Edge:nBJetMedium35_Edge:rhoCN:(Lep1_pdgId_Edge*Lep2_pdgId_Edge)"

        ## define the cut you want to have applied
        #cutString = cuts.AddList([cuts.DYControlRegion, cuts.DYmass, cuts.GoodLeptonee(), cuts.Central()])
        cutString = cuts.AddList([cuts.triggerCalculation, cuts.ee, cuts.Central()]) 

        #print cutString

        actualTree.SetScanField(-1)
        save = os.dup( sys.stdout.fileno() )
        newout = file( 'evtLists/synchOnFullDataset_%d.txt'%(ind), 'w' )
        os.dup2( newout.fileno(), sys.stdout.fileno() )
        actualTree.Scan(scanString, cutString, 'colsize=12')#colsize=8 precision=8 col=8d:12d');
        os.dup2( save, sys.stdout.fileno() )
        newout.close()

   

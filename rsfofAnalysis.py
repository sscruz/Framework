#####################################################################
######                                                              #
###### 8=========D         OO                        8==========D   #  
###### OO                  ||                                ,88   #
###### ||                  ||                              ,88'    #  
###### ||O---\     ,adPPYb,|| ,adPPYb,d8  ,adPPYba,      ,88'      #
###### ||O---/    a8'    `Y||a8'    `Y88 a8P_____88    ,88'        #
###### ||         8b       ||8b       88 8PP'''''''  ,88'          #
###### \/         '8a,   ,d||'8a,   ,d88 '8b,   ,aa 88'            #
###### 8=========D `'8bbdP'\/  `'YbbdP'Y8  `'Ybbd8'' 8==========D   #
######                       aa,    ,88                             #
######                         'Y8bbdP'                             #
######                                                              #
#####################################################################

import ROOT as r
import math as math

from optparse import OptionParser
from ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats
from Sample import Sample, Block, Tree
from CutManager import CutManager
from Canvas import Canvas


def make_rsfof(histo_sf, histo_of):

    ratio = histo_sf.Clone('rsfof_' + histo_sf.GetName())
    ratio.Divide(histo_of)
    ratio.GetYaxis().SetTitle('r_{SFOF}')

    fit = TF1('myfit','pol0', ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax())
    
    ratio.Fit('myfit')

    ratio.GetYaxis().SetRangeUser(0.5,1.5)

    ## ratio.Draw()
    ## fit.Draw('same')
    ## ps = ratio.GetListOfFunctions().FindObject('stats')
    ## print 'this is ps:', ps
    ## ps.SetX1NDC(0.15)
    ## ps.SetX2NDC(0.55)

    ## ps.SetY1NDC(0.15)
    ## ps.SetY2NDC(0.25)

    f = open('txts/'+ratio.GetName()+'_values.txt', 'w')
    for i in range(1, ratio.GetNbinsX()+1):
        min, max = ratio.GetBinLowEdge(i), ratio.GetBinLowEdge(i)+ratio.GetBinWidth(i)
        print    'R_SFOF in [%.2f, %.2f] GeV:\t%.3f +- %.3f'    %(min, max, ratio.GetBinContent(i), ratio.GetBinError(i) )
        f.write( 'R_SFOF in [%.2f, %.2f] GeV:\t%.3f +- %.3f \n' %(min, max, ratio.GetBinContent(i), ratio.GetBinError(i) ) )
    f.close()


    return ratio



if __name__ == '__main__':

    parser = OptionParser(usage='usage: %prog [options] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    parser.add_option('-r', '--region', action='store', type='string', dest='region', default='inclusive', help='region for which to produce Rsfof plots. \n choices are \'inclusive\', \'ttjets\', \'signal\'. default is \'inclusive\'')
    parser.add_option('-t', '--trigger', action='store', type='int', dest='triggerFlag', default='1', help='Trigger cut. Set to 0 if you want to run without trigger')
    (options, args) = parser.parse_args()


    if len(args) != 1:
      parser.error('wrong number of arguments')

    inputFileName = args[0]
    tree = Tree(inputFileName, 'MC', 0)
   
    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager()

    if   options.region == 'inclusive':
       r_cut = cuts.nj2
    elif options.region == 'ttjets':
       r_cut = cuts.METJetsControlRegion
    elif options.region == 'signal':
       r_cut = cuts.METJetsSignalRegion


    if not options.triggerFlag:
        print 'running without trigger...'

    #for var in ['mll', 'met', 'nj', 'nvx']:
    for var in ['mll']:

        print '======================='
        print 'producing rsfof for', var


        if   var =='met':
            t_var  = 'met_pt'
            bin    = [20, 20, 250]
            v_name = 'E_{T}^{miss} (GeV)'

        elif var =='mll':
            t_var = 't.lepsMll_Edge'
            bin    = [[20,70,81,101,120,300], 1, 1]
            v_name = 'm_{ll} (GeV)'

        elif var =='nj':
            t_var = 't.nJetSel_Edge'
            bin    = [ 7,  0.5,   7.5]
            v_name = 'N_{jets}'

        elif var =='nvx':
            t_var = 'nVert'
            bin    = [ 10,  1,  30]
            v_name = 'N_{vertices}'


        trigsf = '(HLT_DoubleMu >0 || HLT_DoubleEl > 0)'
        trigof = '(HLT_MuEG > 0)'

        trig_suffix = ''
        if not options.triggerFlag:
            trigsf = '(HLT_DoubleMu > -1 || HLT_DoubleEl > -1)'
            trigof = '(HLT_MuEG > -1)'
            trig_suffix = '_notrig'

        print trigsf, trigof
        
        for eta in ['central', 'forward']:

            print '...in', eta

            if   eta == 'central': etacut = cuts.Central()
            elif eta == 'forward': etacut = cuts.Forward()

            sf = tree.getTH1F(4, var+'sf_'+eta+'_'+options.region, t_var, bin[0], bin[1], bin[2], cuts.AddList([r_cut, cuts.GoodLeptonSF(), etacut, trigsf]), '', v_name)
            of = tree.getTH1F(4, var+'of_'+eta+'_'+options.region, t_var, bin[0], bin[1], bin[2], cuts.AddList([r_cut, cuts.GoodLeptonOF(), etacut, trigof]), '', v_name)

   
            rsfof= make_rsfof(sf, of)
            plot_rsfof= Canvas('plot_rsfof_'+var+'_'+eta+'_'+options.region+trig_suffix, 'png,pdf', 0.6, 0.6, 0.8, 0.8)
            plot_rsfof.addHisto(rsfof, 'PE', 'OF', 'L', r.kBlack, 1, 0)
            plot_rsfof.addLine(rsfof.GetXaxis().GetXmin(), 1., rsfof.GetXaxis().GetXmax(),1., 3)
            plot_rsfof.save(0, 0, 0, 4.0)
            

            del sf, of, rsfof, plot_rsfof



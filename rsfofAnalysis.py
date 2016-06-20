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
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats
import math, sys, optparse, copy, subprocess

import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample



class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'




def saveInFile(theFile, measuredValueMC, measuredValueUncMC, measuredValueData, measuredValueUncData):

    foutput = open(theFile + "_aux", 'w')
    for line in open(theFile).readlines():
        if line.find("rsfof") != -1 and line.find("direct") != -1:
            if line.find("DATA") != -1:
                foutput.write('rsfof       direct          DATA        %.4f      %0.4f       %.4f\n'%(measuredValueData, measuredValueUncData, measuredValueUncMC))
            else:
                foutput.write('rsfof       direct          MC          %.4f      %0.4f       %.4f\n'%(measuredValueMC, measuredValueUncMC, 0))
        else:
            foutput.write(line)

    foutput.close()
    subprocess.call(['mv ' + theFile + "_aux " + theFile], shell=True)
   


def make_rsfof(histo_sf, histo_of, dataMC):

    ratio = histo_sf.Clone('rsfof_' + histo_sf.GetName())
    ratio.Divide(histo_of)
    ratio.GetYaxis().SetTitle('r_{SFOF}')

    doFit = 0
    if doFit:
        fit = TF1('myfit','pol0', ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax())
        fit.SetLineColor(r.kBlack if dataMC  == 'DATA' else r.kRed+1)
        fit.SetLineStyle(2)
        
        ratio.Fit('myfit','0')

    ratio.GetYaxis().SetRangeUser(0.5,1.5)

    return ratio




def runAnalysis(treeDA, treeMC, cuts, specialcut, tag, save, ingredientsFile):


    labelx = "m_{ll} [GeV]"
    labelmet = "MET [GeV]"
    labelnjet = "N. Jets"


    #####Main mll control plot
    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Starting mll plot' + bcolors.ENDC
    
    MCControlSF =        treeMC.getTH1F(lumi, "MCControlSF", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectControlRegion, cuts.SF]), '', labelx)
    MCControlOF =        treeMC.getTH1F(lumi, "MCControlOF", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectControlRegion, cuts.OF]), '', labelx)
    MCSignalSF =         treeMC.getTH1F(lumi, "MCSignalSF", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectSignalRegion, cuts.SF]), '', labelx)
    MCSignalOF =         treeMC.getTH1F(lumi, "MCSignalOF", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectSignalRegion, cuts.OF]), '', labelx)
    MCControlSFvalue =   treeMC.getTH1F(lumi, "MCControlSFvalue", "lepsMll_Edge", [20, 300], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectControlRegionRestrictive, cuts.SF]), '', labelx)
    MCControlOFvalue =   treeMC.getTH1F(lumi, "MCControlOFvalue", "lepsMll_Edge", [20, 300], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectControlRegionRestrictive, cuts.OF]), '', labelx)
    MCSignalSFvalue =    treeMC.getTH1F(lumi, "MCSignalSFvalue", "lepsMll_Edge", [20, 300], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectSignalRegionRestrictive, cuts.SF]), '', labelx)
    MCSignalOFvalue =    treeMC.getTH1F(lumi, "MCSignalOFvalue", "lepsMll_Edge", [20, 300], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectSignalRegionRestrictive, cuts.OF]), '', labelx)
    DataControlSF =      treeDA.getTH1F(lumi, "DataControlSF", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut,cuts.goodLepton,cuts.RSFOFDirectControlRegion, cuts.SF]), '', labelx)
    DataControlOF =      treeDA.getTH1F(lumi, "DataControlOF", "lepsMll_Edge", [20, 70, 81, 101, 111, 300], 1, 1, cuts.AddList([specialcut,cuts.goodLepton,cuts.RSFOFDirectControlRegion, cuts.OF]), '', labelx)
    DataControlSFvalue = treeDA.getTH1F(lumi, "DataControlSFvalue", "lepsMll_Edge", [20, 300], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectControlRegionRestrictive, cuts.SF]), '', labelx)
    DataControlOFvalue = treeDA.getTH1F(lumi, "DataControlOFvalue", "lepsMll_Edge", [20, 300], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectControlRegionRestrictive, cuts.OF]), '', labelx)
    
    MCControl =          make_rsfof(MCControlSF, MCControlOF, "MC")
    MCSignal =           make_rsfof(MCSignalSF, MCSignalOF, "MC")
    MCControlvalue =     make_rsfof(MCControlSFvalue, MCControlOFvalue, "MC")
    MCSignalvalue =      make_rsfof(MCSignalSFvalue, MCSignalOFvalue, "MC")
    DataControl =        make_rsfof(MCControlSF, MCControlOF, "DATA")
    DataControlvalue =   make_rsfof(MCControlSFvalue, MCControlOFvalue, "DATA")

    measuredValueMC =         MCControlvalue.GetBinContent(1)
    measuredValueUncMC =      MCControlvalue.GetBinError(1)
    measuredValueData =       DataControlvalue.GetBinContent(1)
    measuredValueUncData =    DataControlvalue.GetBinError(1)
    measuredValueUncTotData = math.sqrt(measuredValueUncData**2 + measuredValueUncMC**2)

    plot_rsfof = Canvas.Canvas('rsfof/%s_%s/plot_rsfof_mll'%(lumi_str,tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
    plot_rsfof.addHisto(MCControl, 'PE', 'ttjets region - MC', 'PL', r.kRed+1 , 1, 0)
    plot_rsfof.addHisto(MCControlvalue, 'PE,SAME', 'ttjets measurement - MC', 'PL', r.kBlue , 1, 0)
    plot_rsfof.addHisto(MCSignal, 'PE,SAME', 'Signal region - MC', 'PL', r.kRed+1 , 1, 0)
    plot_rsfof.addHisto(DataControl, 'PE,SAME', 'ttjets region - DATA', 'PL', r.kBlack , 1, 0)
    plot_rsfof.addHisto(DataControlvalue, 'PE,SAME', 'ttjets measurement - DATA', 'PL', r.kGreen , 1, 0)
    plot_rsfof.addBand (MCControl.GetXaxis().GetXmin(), measuredValueData-measuredValueUncTotData, MCControl.GetXaxis().GetXmax(), measuredValueData+measuredValueUncTotData, r.kGreen, 0.2)
    plot_rsfof.addLine (MCControl.GetXaxis().GetXmin(), measuredValueData, MCControl.GetXaxis().GetXmax(), measuredValueData, r.kGreen)
    plot_rsfof.save(1, 1, 0, lumi, 0.2, 1.8)
    
    if save==True:  
        saveInFile(ingredientsFile, measuredValueMC, measuredValueUncMC, measuredValueData, measuredValueUncData)    


    ######Met and JET plots 
    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Starting jet and met plot' + bcolors.ENDC
    MCSignalSFMET =      treeMC.getTH1F(lumi, "MCSignalSFMET", "met_Edge", [20, 40, 60, 80, 100], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFDirectSignalRegionRestrictiveNoMET, cuts.SF]), '', labelmet)
    MCSignalOFMET =      treeMC.getTH1F(lumi, "MCSignalOFMET", "met_Edge", [20, 40, 60, 80, 100], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFDirectSignalRegionRestrictiveNoMET, cuts.OF]), '', labelmet)
    MCSignalSFJet =      treeMC.getTH1F(lumi, "MCSignalSFJets", "nJetSel_Edge", 10, 0, 10, cuts.AddList([cuts.goodLepton, cuts.RSFOFDirectSignalRegionRestrictiveNoJet, cuts.SF]), '', labelnjet)
    MCSignalOFJet =      treeMC.getTH1F(lumi, "MCSignalOFJets", "nJetSel_Edge", 10, 0, 10, cuts.AddList([cuts.goodLepton, cuts.RSFOFDirectSignalRegionRestrictiveNoJet, cuts.OF]), '', labelnjet)
    
    MCSignalMET   =      make_rsfof(MCSignalSFMET, MCSignalOFMET, "MC")
    MCSignalJet   =      make_rsfof(MCSignalSFJet, MCSignalOFJet, "MC")

    plot_rsfofmet = Canvas.Canvas('rsfof/%s_%s/plot_rsfof_met'%(lumi_str,tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
    plot_rsfofmet.addHisto(MCSignalMET, 'PE', 'ttjets region - MC', 'PL', r.kRed+1 , 1, 0)
    plot_rsfofmet.addBand(MCSignalMET.GetXaxis().GetXmin(), measuredValueMC-measuredValueUncMC, MCSignalMET.GetXaxis().GetXmax(), measuredValueMC+measuredValueUncMC, r.kGreen, 0.2)
    plot_rsfofmet.addLine(MCSignalMET.GetXaxis().GetXmin(), measuredValueMC, MCSignalMET.GetXaxis().GetXmax(), measuredValueMC, r.kGreen)
    plot_rsfofmet.save(1, 1, 0, lumi, 0.2, 1.8)

    plot_rsfofjet = Canvas.Canvas('rsfof/%s_%s/plot_rsfof_jet'%(lumi_str,tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
    plot_rsfofjet.addHisto(MCSignalJet, 'PE', 'ttjets region - MC', 'PL', r.kRed+1 , 1, 0)
    plot_rsfofjet.addBand(MCSignalJet.GetXaxis().GetXmin(), measuredValueMC-measuredValueUncMC, MCSignalJet.GetXaxis().GetXmax(), measuredValueMC+measuredValueUncMC, r.kGreen, 0.2)
    plot_rsfofjet.addLine(MCSignalJet.GetXaxis().GetXmin(), measuredValueMC, MCSignalJet.GetXaxis().GetXmax(), measuredValueMC, r.kGreen)
    plot_rsfofjet.save(1, 1, 0, lumi, 0.2, 1.8)
     



##Main body of the analysis
if __name__ == '__main__':

    print bcolors.HEADER 
    print '#######################################################################'
    print '                  Starting r_SFOF analysis...                          ' 
    print '#######################################################################' + bcolors.ENDC

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()

    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    mcDatasets = ['TTJets_DiLepton', 'DYJetsToLL_M10to50', 'DYJetsToLL_M50', 'WW', 'WZ', 'ZZ']
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
    
    runAnalysis(treeDA, treeMC, cuts, '', 'nocut', True, opts.ingredientsFile)
 

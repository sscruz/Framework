#####################################################################
######                                                              #
###### 88888888888         88                        888888888888   #  
###### 88                  88                                 ,88   #
###### 88                  88                               ,88"    #  
###### 88aaaaa     ,adPPYb,88  ,adPPYb,d8  ,adPPYba,      ,88"      #
###### 88"""""    a8"    `Y88 a8"    `Y88 a8P_____88    ,88"        #
###### 88         8b       88 8b       88 8PP"""""""  ,88"          #
###### 88         "8a,   ,d88 "8a,   ,d88 "8b,   ,aa 88"            #
###### 88888888888 `"8bbdP"Y8  `"YbbdP"Y8  `"Ybbd8"' 888888888888   #
######                       aa,    ,88                             #
######                         "Y8bbdP"                             #
######                                                              #
#####################################################################

import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors
import math, sys, optparse, array, copy, subprocess

import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Rounder    as Rounder

class bcolors: 
    HEADER = '\033[95m' 
    OKBLUE = '\033[94m' 
    OKGREEN = '\033[92m' 
    WARNING = '\033[93m' 
    FAIL = '\033[91m' 
    ENDC = '\033[0m' 
    BOLD = '\033[1m' 
    UNDERLINE = '\033[4m' 



def saveInFile(theFile, measuredValueMC, measuredValueUncMC, measuredValueUncSystMC, measuredValueData, measuredValueUncData, measuredValueUncSystData, tag):

    foutput = open(theFile + "_aux", 'w')
    for line in open(theFile).readlines():
        if line.find("rmue") != -1 and line.find(tag) != -1:
            if line.find("DATA") != -1:
                foutput.write('rmue        alone           DATA        %.4f      %0.4f       %.4f\n'%(measuredValueData, measuredValueUncData, measuredValueUncSystData))
            else:
                foutput.write('rmue        alone           MC          %.4f      %0.4f       %.4f\n'%(measuredValueMC, measuredValueUncMC, measuredValueUncSystMC))
        else:
            foutput.write(line)

    foutput.close()
    subprocess.call(['mv ' + theFile + "_aux " + theFile], shell=True)



def calc_rmue(Nmm, Nee, Emm, Eee):

    val = [0, 0]
    print Nmm, Nee
    if(Nmm >= 0 and Nee > 0):
        val[0] = math.sqrt(Nmm/Nee)
        val[1] = math.sqrt(0.25 * Emm * Emm / (Nmm * Nee) + 0.25 * Eee * Eee * Nmm / (Nee * Nee * Nee))
    else:
        print bcolors.HEADER + "[rmueAnalysis] " + bcolors.FAIL + " The yields in Nee and Nmm are <=0 " + bcolors.ENDC
    return val


def make_rmue(histo_mm, histo_ee):

    ratio = histo_mm.Clone("rmue_" + histo_mm.GetName())
    ratio.GetYaxis().SetTitle("r_{#mu e}")
    ratio.GetXaxis().SetTitle(histo_mm.GetXaxis().GetTitle())
  
    for i in range(0, histo_mm.GetNbinsX()+1):
        Nmm = histo_mm.GetBinContent(i)
        Nee = histo_ee.GetBinContent(i)
        Emm = histo_mm.GetBinError(i)
        Eee = histo_ee.GetBinError(i)
        if(Nee != 0):
            if(Nee == 0 or Nmm == 0):
              continue
            val = calc_rmue(Nmm, Nee, Emm, Eee)
            ratio.SetBinContent(i, val[0])
            ratio.SetBinError(i, val[1])
    return ratio



def convertToFactor(histo, getGraph=False):
    tmp_histo = copy.deepcopy(histo)
    #tmp_histo.Sumw2()
    for i in range(tmp_histo.GetNbinsX()+1):
        rmue     = tmp_histo.GetBinContent(i)
        rmue_err = tmp_histo.GetBinError  (i)
        if rmue:
            fac = 0.5*(rmue + 1./rmue)
            err = 0.5*((1. - 1./(rmue**2))*rmue_err)
        else:
            fac = 0.
            err = 0.
        tmp_histo.SetBinContent(i, fac)
        tmp_histo.SetBinError  (i, err)
    if getGraph:
        return TGraphErrors(tmp_histo)
    else:
        return tmp_histo


def makeAnalysis(treeDA, treeMC, cuts, specialcut, tag, save, ingredientsFile):

    print bcolors.HEADER + '[rmueAnalysis] ' + bcolors.OKBLUE + 'Producing histograms...' + bcolors.ENDC
    MCDYControlMllee =         treeMC.getTH1F(lumi, "MCDYControlMllee", "lepsMll_Edge", 20, 20, 200, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoMll, cuts.ee]), '', labelx)
    MCDYControlMllmm =         treeMC.getTH1F(lumi, "MCDYControlMllmm", "lepsMll_Edge", 20, 20, 200, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoMll, cuts.mm]), '', labelx)
    DATADYControlMllee =       treeDA.getTH1F(lumi, "DATADYControlMllee", "lepsMll_Edge", 20, 20, 200, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoMll, cuts.ee]), '', labelx)
    DATADYControlMllmm =       treeDA.getTH1F(lumi, "DATADYControlMllmm", "lepsMll_Edge", 20, 20, 200, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoMll, cuts.mm]), '', labelx)
    MCDYControlMlleevalue =    treeMC.getTH1F(lumi, "MCDYControlMlleevalue", "lepsMll_Edge", 1, 20, 300, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegion, cuts.ee]), '', labelx)
    MCDYControlMllmmvalue =    treeMC.getTH1F(lumi, "MCDYControlMllmmvalue", "lepsMll_Edge", 1, 20, 300, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegion, cuts.mm]), '', labelx)
    DATADYControlMlleevalue =  treeDA.getTH1F(lumi, "DATADYControlMlleevalue", "lepsMll_Edge", 1, 20, 300, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegion, cuts.ee]), '', labelx)
    DATADYControlMllmmvalue =  treeDA.getTH1F(lumi, "DATADYControlMllmmvalue", "lepsMll_Edge", 1, 20, 300, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegion, cuts.mm]), '', labelx)
    MCDYControlMETee =         treeMC.getTH1F(lumi, "MCDYControlMETee", "met_Edge", 5, 20, 200, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoMET, cuts.ee]), '', labelmet)
    MCDYControlMETmm =         treeMC.getTH1F(lumi, "MCDYControlMETmm", "met_Edge", 5, 20, 200, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoMET, cuts.mm]), '', labelmet)
    DATADYControlMETee =       treeDA.getTH1F(lumi, "DATADYControlMETee", "met_Edge", 5, 20, 200, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoMET, cuts.ee]), '', labelmet)
    DATADYControlMETmm =       treeDA.getTH1F(lumi, "DATADYControlMETmm", "met_Edge", 5, 20, 200, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoMET, cuts.mm]), '', labelmet)
    MCDYControlJetee =         treeMC.getTH1F(lumi, "MCDYControlJetee", "met_Edge", 5, 20, 200, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoJet, cuts.ee]), '', labelnjet)
    MCDYControlJetmm =         treeMC.getTH1F(lumi, "MCDYControlJetmm", "met_Edge", 5, 20, 200, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoJet, cuts.mm]), '', labelnjet)
    DATADYControlJetee =       treeDA.getTH1F(lumi, "DATADYControlJetee", "met_Edge", 5, 20, 200, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoJet, cuts.ee]), '', labelnjet)
    DATADYControlJetmm =       treeDA.getTH1F(lumi, "DATADYControlJetmm", "met_Edge", 5, 20, 200, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoJet, cuts.mm]), '', labelnjet)
 
    MCDYControlMll =           make_rmue(MCDYControlMllmm, MCDYControlMllee)
    DATADYControlMll =         make_rmue(DATADYControlMllmm, DATADYControlMllee)
    MCDYControlMllvalue =      make_rmue(MCDYControlMllmmvalue, MCDYControlMlleevalue)
    DATADYControlMllvalue =    make_rmue(DATADYControlMllmmvalue, DATADYControlMlleevalue)
    MCDYControlMET =           make_rmue(MCDYControlMETmm, MCDYControlMETee)
    DATADYControlMET =         make_rmue(DATADYControlMETmm, DATADYControlMETee)
    MCDYControlJet =           make_rmue(MCDYControlMETmm, MCDYControlMETee)
    DATADYControlJet =         make_rmue(DATADYControlJetmm, DATADYControlJetee)
   
    factorMCDYControlMll =          convertToFactor(MCDYControlMll)
    factorDATADYControlMll =        convertToFactor(DATADYControlMll)
    factorMCDYControlMllvalue =     convertToFactor(MCDYControlMllvalue)
    factorDATADYControlMllvalue =   convertToFactor(DATADYControlMllvalue)
    factorMCDYControlMET =          convertToFactor(MCDYControlMET)
    factorDATADYControlMET =        convertToFactor(DATADYControlMET)
    factorMCDYControlJet =          convertToFactor(MCDYControlJet)
    factorDATADYControlJet =        convertToFactor(DATADYControlJet)
    
    systematicForrmue = 0.1 
    MCrmuemeasured =               MCDYControlMllvalue.GetBinContent(1)
    MCrmuemeasuredUnc =            MCDYControlMllvalue.GetBinError(1)
    MCrmuemeasuredUncSyst =        MCrmuemeasured * systematicForrmue
    DATArmuemeasured =             DATADYControlMllvalue.GetBinContent(1)
    DATArmuemeasuredUnc =          DATADYControlMllvalue.GetBinError(1)
    DATArmuemeasuredUncSyst =      DATArmuemeasured * systematicForrmue
    DATArmuemeasuredUncTot  =      math.sqrt(DATArmuemeasuredUnc**2 + DATArmuemeasuredUncSyst**2)
    factorMCrmuemeasured =             factorMCDYControlMllvalue.GetBinContent(1)
    factorMCrmuemeasuredUnc =          factorMCDYControlMllvalue.GetBinError(1)
    factorMCrmuemeasuredUncSyst =      (1.0-1.0/MCrmuemeasured**2) * MCrmuemeasuredUncSyst if MCrmuemeasured != 0 else 1000.00
    factorMCrmuemeasuredUncTot  =      math.sqrt(MCrmuemeasuredUnc**2 + MCrmuemeasuredUncSyst**2)
    factorDATArmuemeasured =             factorDATADYControlMllvalue.GetBinContent(1)
    factorDATArmuemeasuredUnc =          factorDATADYControlMllvalue.GetBinError(1)
    factorDATArmuemeasuredUncSyst =      (1.0-1.0/DATArmuemeasured**2) * DATArmuemeasuredUncSyst if DATArmuemeasured != 0 else 1000.00
    factorDATArmuemeasuredUncTot  =      math.sqrt(DATArmuemeasuredUnc**2 + DATArmuemeasuredUncSyst**2)

    if save==True:
        saveInFile(ingredientsFile, MCrmuemeasured, MCrmuemeasuredUnc, MCrmuemeasuredUncSyst, DATArmuemeasured, DATArmuemeasuredUnc, DATArmuemeasuredUncSyst, 'alone')
        saveInFile(ingredientsFile, factorMCrmuemeasured, factorMCrmuemeasuredUnc, factorMCrmuemeasuredUncSyst, factorDATArmuemeasured, factorDATArmuemeasuredUnc, factorDATArmuemeasuredUncSyst, 'factor')


    print bcolors.HEADER + '[rmueAnalysis] ' + bcolors.OKBLUE + 'Producing plots...' + bcolors.ENDC
    plot_rmue_mll = Canvas.Canvas('rmue/%s_%s/plot_rmue_mll'%(lumi_str, tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
    plot_rmue_mll.addHisto(MCDYControlMll, 'PE', 'MC', 'PL', r.kRed+1 , 1, 0)
    plot_rmue_mll.addHisto(DATADYControlMll, 'PE,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_rmue_mll.addBand(MCDYControlMll.GetXaxis().GetXmin(), DATArmuemeasured-DATArmuemeasuredUncTot, MCDYControlMll.GetXaxis().GetXmax(), DATArmuemeasured+DATArmuemeasuredUncTot, r.kGreen, 0.2)
    plot_rmue_mll.addLine(MCDYControlMll.GetXaxis().GetXmin(), DATArmuemeasured, MCDYControlMll.GetXaxis().GetXmax(), DATArmuemeasured,r.kGreen)
    plot_rmue_mll.save(1, 1, 0, lumi, 0.2, 1.8)
    
    plot_rmue_met = Canvas.Canvas('rmue/%s_%s/plot_rmue_met'%(lumi_str, tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
    plot_rmue_met.addHisto(MCDYControlMET, 'PE', 'MC', 'PL', r.kRed+1 , 1, 0)
    plot_rmue_met.addHisto(DATADYControlMET, 'PE,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_rmue_met.addBand(MCDYControlMET.GetXaxis().GetXmin(), DATArmuemeasured-DATArmuemeasuredUncTot, MCDYControlMET.GetXaxis().GetXmax(), DATArmuemeasured+DATArmuemeasuredUncTot, r.kGreen, 0.2)
    plot_rmue_met.addLine(MCDYControlMET.GetXaxis().GetXmin(), DATArmuemeasured, MCDYControlMET.GetXaxis().GetXmax(), DATArmuemeasured,r.kGreen)
    plot_rmue_met.save(1, 1, 0, lumi, 0.2, 1.8)
  
    plot_rmue_jet = Canvas.Canvas('rmue/%s_%s/plot_rmue_jet'%(lumi_str, tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
    plot_rmue_jet.addHisto(MCDYControlJet, 'PE', 'MC', 'PL', r.kRed+1 , 1, 0)
    plot_rmue_jet.addHisto(DATADYControlJet, 'PE,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_rmue_jet.addBand(MCDYControlJet.GetXaxis().GetXmin(), DATArmuemeasured-DATArmuemeasuredUncTot, MCDYControlJet.GetXaxis().GetXmax(), DATArmuemeasured+DATArmuemeasuredUncTot, r.kGreen, 0.2)
    plot_rmue_jet.addLine(MCDYControlJet.GetXaxis().GetXmin(), DATArmuemeasured, MCDYControlJet.GetXaxis().GetXmax(), DATArmuemeasured,r.kGreen)
    plot_rmue_jet.save(1, 1, 0, lumi, 0.2, 1.8)
    
    plot_factor_mll = Canvas.Canvas('rmue/%s_%s/plot_factor_mll'%(lumi_str, tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
    plot_factor_mll.addHisto(factorMCDYControlMll, 'PE', 'MC', 'PL', r.kRed+1 , 1, 0)
    plot_factor_mll.addHisto(factorDATADYControlMll, 'PE,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_factor_mll.addBand(factorMCDYControlMll.GetXaxis().GetXmin(), factorDATArmuemeasured-factorDATArmuemeasuredUncTot, factorMCDYControlMll.GetXaxis().GetXmax(), factorDATArmuemeasured+factorDATArmuemeasuredUncTot, r.kGreen, 0.2)
    plot_factor_mll.addLine(factorMCDYControlMll.GetXaxis().GetXmin(), factorDATArmuemeasured, factorMCDYControlMll.GetXaxis().GetXmax(), factorDATArmuemeasured,r.kGreen)
    plot_factor_mll.save(1, 1, 0, lumi, 0.2, 1.8)
    
 
 
if __name__ == "__main__":

    print bcolors.HEADER
    print '#######################################################################'
    print '                  Starting r_mue analysis...                          '
    print '#######################################################################' + bcolors.ENDC

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()

    print bcolors.HEADER + '[rmueAnalysis] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    mcDatasets = ['TTJets_DiLepton', 'TTJets_DiLepton_ext', 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2', 'DoubleEG_Run2016B-PromptReco-v2', 'MuonEG_Run2016B-PromptReco-v2']

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print bcolors.HEADER + '[rmueAnalysis] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    #lumi = 0.864 ; maxrun = 274240
    lumi = 2.66  ; maxrun = 999999
    lumi_str = str(lumi).replace('.','p')+'invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()
    cuts = CutManager.CutManager()

    labelx = "m_{ll} [GeV]"
    labelmet = "MET [GeV]"
    labelnjet = "N. Jets"

    ####Cuts needed by rmue
    cuts = CutManager.CutManager()

    makeAnalysis(treeDA, treeMC, cuts, '', 'nocut', True, opts.ingredientsFile)


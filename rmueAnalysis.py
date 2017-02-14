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

import include.LeptonSF
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


def makeTable(MCmm, MCee, DATAmm, DATAee, MCrmue, MCrmueUnc, MCrmueUncSyst, DATArmue, DATArmueUnc, DATArmueUncSyst):
    line0 = '  \hline'
    line1 = '  Data &' 
    line2 = '  MC   &' 
    line0 += ' & $\mathrm{N_{\mu\mu}}$  & $\mathrm{N_{ee}}$ & $\mathrm{r_{\mu e } \pm \sigma_{stat} \pm \sigma_{syst}}$ \\\\'
    line1 += ' %.f & %.f &    %.2f $\\pm$ %.2f $\\pm$ %.2f      %s' %(DATAmm.Integral(), DATAee.Integral(), DATArmue  , DATArmueUnc  , DATArmueUncSyst,  '\\\\')
    line2 += ' %.f & %.f &    %.2f $\\pm$ %.2f $\\pm$ %.2f      %s' %(MCmm.Integral(), MCee.Integral(), MCrmue  , MCrmueUnc  , MCrmueUncSyst,  '\\\\')
    line0 += '\\hline'; line2 += '\\hline';

    helper.ensureDirectory('plots/rmue/%s/'%lumi_str)
    helper.ensureDirectory('plots/rmue/%s/tables/'%lumi_str)
    compTableFile = open('plots/rmue/%s/tables/resultTable_%s%s.txt'%(lumi_str, str(lumi).replace('.','p'), "rmue"),'w')
    compTableFile.write(line0+'\n')
    compTableFile.write(line1+'\n')
    compTableFile.write(line2+'\n')                                                                                             
    print line0
    print line1
    print line2                                                                                                                                                                      


def saveInFile(theFile, measuredValueMC, measuredValueUncMC, measuredValueUncSystMC, measuredValueData, measuredValueUncData, measuredValueUncSystData, tag):

    foutput = open(theFile + "_aux", 'w')
    for line in open(theFile).readlines():
        if line.find("rmue") != -1 and line.find(tag) != -1 and  line.find("alone") != -1 :
            if line.find("DATA") != -1:
                foutput.write('rmue        alone           DATA        %.4f      %0.4f       %.4f\n'%(measuredValueData, measuredValueUncData, measuredValueUncSystData))
            else:
                foutput.write('rmue        alone           MC          %.4f      %0.4f       %.4f\n'%(measuredValueMC, measuredValueUncMC, measuredValueUncSystMC))           
        elif line.find("rmue") != -1 and line.find(tag) != -1 and  line.find("factor") != -1 :
            if line.find("DATA") != -1:
                foutput.write('rmue        factor           DATA        %.4f      %0.4f       \n'%(getFactor(measuredValueData, measuredValueUncData, measuredValueUncSystData)))
            else:
                foutput.write('rmue        factor           MC          %.4f      %0.4f       \n'%(getFactor(measuredValueMC, measuredValueUncMC, measuredValueUncSystMC)))           
        
        elif line.find("rmue") != -1 and line.find(tag) != -1 and  line.find("coeffA") != -1 :
            if line.find("DATA") != -1:
                foutput.write('rmue        coeffA           DATA        %.4f      %0.4f       \n'%(measuredValueData, measuredValueUncSystData))
            else:
                foutput.write('rmue        coeffA           MC        %.4f      %0.4f       \n'%(measuredValueMC, measuredValueUncSystMC))

        elif line.find("rmue") != -1 and line.find(tag) != -1 and  line.find("coeffB") != -1 :
            if line.find("DATA") != -1:
                foutput.write('rmue        coeffB           DATA        %.4f      %0.4f       \n'%(measuredValueData, measuredValueUncSystData))
            else:
                foutput.write('rmue        coeffB           MC        %.4f      %0.4f       \n'%(measuredValueMC, measuredValueUncSystMC))
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
        tmp_histo.GetYaxis().SetTitle("0.5(r_{#mu e} + 1/r_{#mu e})")
    if getGraph:
        return TGraphErrors(tmp_histo)
    else:
        return tmp_histo                                                         


def getFactor(rmue, rmue_err, rmue_err_syst):
    fac = 0.5*(rmue + 1./rmue)
    err = 0.5*((1. - 1./(rmue**2))*math.sqrt(rmue_err**2 + rmue_err_syst**2) )
    return fac, err                                                     


def makeAnalysis(treeDA, treeMC, cuts, specialcut, tag, save, ingredientsFile):
    print bcolors.HEADER + '[rmueAnalysis] ' + bcolors.OKBLUE + 'Producing histograms...' + bcolors.ENDC
    MCDYControlMllee =         treeMC.getTH1F(lumi, "MCDYControlMllee", "lepsMll_Edge", 20, 20, 200, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoMll, cuts.ee]), '', labelx)
    MCDYControlMllmm =         treeMC.getTH1F(lumi, "MCDYControlMllmm", "lepsMll_Edge", 20, 20, 200, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoMll, cuts.mm]), '', labelx)
    DATADYControlMllee =       treeDA.getTH1F(lumi, "DATADYControlMllee", "lepsMll_Edge", 20, 20, 200, cuts.AddList([specialcut, cuts.trigger, cuts.goodLepton, cuts.DYControlRegionNoMll, cuts.ee]), '', labelx)
    DATADYControlMllmm =       treeDA.getTH1F(lumi, "DATADYControlMllmm", "lepsMll_Edge", 20, 20, 200, cuts.AddList([specialcut, cuts.trigger, cuts.goodLepton, cuts.DYControlRegionNoMll, cuts.mm]), '', labelx)
    MCDYControlMlleevalue =    treeMC.getTH1F(lumi, "MCDYControlMlleevalue", "lepsMll_Edge", 1, 20, 300, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegion, cuts.ee]), '', labelx)
    MCDYControlMllmmvalue =    treeMC.getTH1F(lumi, "MCDYControlMllmmvalue", "lepsMll_Edge", 1, 20, 300, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegion, cuts.mm]), '', labelx)
    DATADYControlMlleevalue =  treeDA.getTH1F(lumi, "DATADYControlMlleevalue", "lepsMll_Edge", 1, 20, 300, cuts.AddList([specialcut, cuts.trigger, cuts.goodLepton, cuts.DYControlRegion, cuts.ee]), '', labelx)
    DATADYControlMllmmvalue =  treeDA.getTH1F(lumi, "DATADYControlMllmmvalue", "lepsMll_Edge", 1, 20, 300, cuts.AddList([specialcut, cuts.trigger, cuts.goodLepton, cuts.DYControlRegion, cuts.mm]), '', labelx)
    MCDYControlMETee =         treeMC.getTH1F(lumi, "MCDYControlMETee", "met_Edge", 5, 20, 200, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoMET, cuts.ee]), '', labelmet)
    MCDYControlMETmm =         treeMC.getTH1F(lumi, "MCDYControlMETmm", "met_Edge", 5, 20, 200, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoMET, cuts.mm]), '', labelmet)
    DATADYControlMETee =       treeDA.getTH1F(lumi, "DATADYControlMETee", "met_Edge", 5, 20, 200, cuts.AddList([specialcut, cuts.trigger, cuts.goodLepton, cuts.DYControlRegionNoMET, cuts.ee]), '', labelmet)
    DATADYControlMETmm =       treeDA.getTH1F(lumi, "DATADYControlMETmm", "met_Edge", 5, 20, 200, cuts.AddList([specialcut, cuts.trigger, cuts.goodLepton, cuts.DYControlRegionNoMET, cuts.mm]), '', labelmet)
    MCDYControlJetee =         treeMC.getTH1F(lumi, "MCDYControlJetee", "nJetSel_Edge", 10, 0, 10, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoJet, cuts.ee]), '', labelnjet)
    MCDYControlJetmm =         treeMC.getTH1F(lumi, "MCDYControlJetmm", "nJetSel_Edge", 10, 0, 10, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoJet, cuts.mm]), '', labelnjet)
    DATADYControlJetee =       treeDA.getTH1F(lumi, "DATADYControlJetee", "nJetSel_Edge", 10, 0, 10, cuts.AddList([specialcut, cuts.trigger, cuts.goodLepton, cuts.DYControlRegionNoJet, cuts.ee]), '', labelnjet)
    DATADYControlJetmm =       treeDA.getTH1F(lumi, "DATADYControlJetmm", "nJetSel_Edge", 10, 0, 10, cuts.AddList([specialcut, cuts.trigger, cuts.goodLepton, cuts.DYControlRegionNoJet, cuts.mm]), '', labelnjet)
    MCDYControlmt2ee =         treeMC.getTH1F(lumi, "MCDYControlmt2ee", "mt2_Edge", 8, 0, 160, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoJet, cuts.ee]), '', labelmt2)
    MCDYControlmt2mm =         treeMC.getTH1F(lumi, "MCDYControlmt2mm", "mt2_Edge", 8, 0, 160, cuts.AddList([specialcut, cuts.goodLepton, cuts.DYControlRegionNoJet, cuts.mm]), '', labelmt2)
    DATADYControlmt2ee =       treeDA.getTH1F(lumi, "DATADYControlmt2ee", "mt2_Edge", 8, 0, 160, cuts.AddList([specialcut, cuts.trigger, cuts.goodLepton, cuts.DYControlRegionNoJet, cuts.ee]), '', labelmt2)
    DATADYControlmt2mm =       treeDA.getTH1F(lumi, "DATADYControlmt2mm", "mt2_Edge", 8, 0, 160, cuts.AddList([specialcut, cuts.trigger, cuts.goodLepton, cuts.DYControlRegionNoJet, cuts.mm]), '', labelmt2)

    MCDYControlPT2ee =    treeMC.getTH1F(lumi, "MCDYControlPT2ee", "Lep2_pt_Edge", 8, 20, 200, cuts.AddList([specialcut,  cuts.DYControlRegion, cuts.ee]), '', labelpt2)
    MCDYControlPT2mm =    treeMC.getTH1F(lumi, "MCDYControlPT2mm", "Lep2_pt_Edge", 8, 20, 200, cuts.AddList([specialcut,  cuts.DYControlRegion, cuts.mm]), '', labelpt2)
    DATADYControlPT2ee =  treeDA.getTH1F(lumi, "DATADYControlPT2ee", "Lep2_pt_Edge", 8, 20, 200, cuts.AddList([specialcut, cuts.trigger, cuts.DYControlRegion, cuts.ee]), '', labelpt2)
    DATADYControlPT2mm =  treeDA.getTH1F(lumi, "DATADYControlPT2mm", "Lep2_pt_Edge", 8, 20, 200, cuts.AddList([specialcut, cuts.trigger, cuts.DYControlRegion, cuts.mm]), '', labelpt2)

    
    print "got all TH1Fs" 
    MCDYControlMll =           make_rmue(MCDYControlMllmm, MCDYControlMllee)
    DATADYControlMll =         make_rmue(DATADYControlMllmm, DATADYControlMllee)
    MCDYControlMllvalue =      make_rmue(MCDYControlMllmmvalue, MCDYControlMlleevalue)
    DATADYControlMllvalue =    make_rmue(DATADYControlMllmmvalue, DATADYControlMlleevalue)
    MCDYControlMET =           make_rmue(MCDYControlMETmm, MCDYControlMETee)
    DATADYControlMET =         make_rmue(DATADYControlMETmm, DATADYControlMETee)
    MCDYControlJet =           make_rmue(MCDYControlJetmm, MCDYControlJetee)
    DATADYControlJet =         make_rmue(DATADYControlJetmm, DATADYControlJetee)
    MCDYControlmt2 =           make_rmue(MCDYControlmt2mm, MCDYControlmt2ee)
    DATADYControlmt2 =         make_rmue(DATADYControlmt2mm, DATADYControlmt2ee)

    DATADYControlPT2   =         make_rmue(DATADYControlPT2mm, DATADYControlPT2ee)
    MCDYControlPT2 =         make_rmue(MCDYControlPT2mm, MCDYControlPT2ee)


    print 'fits to pt2'
    fun_da = r.TF1('fun_da','[0]+[1]/x', 20, 400)
    fun_mc = r.TF1('fun_mc','[0]+[1]/x', 20, 400)
    fun_da.SetLineWidth(2); fun_mc.SetLineWidth(2)
    DATADYControlPT2.Fit('fun_da','rQ'); 
    data_a = fun_da.GetParameter(0); data_a_e = fun_da.GetParError(0)
    data_b = fun_da.GetParameter(1); data_b_e = fun_da.GetParError(1)
    MCDYControlPT2.Fit('fun_mc','rQ')
    mc_a = fun_mc.GetParameter(0); mc_a_e = fun_mc.GetParError(0)
    mc_b = fun_mc.GetParameter(1); mc_b_e = fun_mc.GetParError(1)
    da_string = 'Fit on data: (%4.2f #pm %4.2f) + (%4.2f #pm %4.2f) p_{T}^{-1}'%(data_a,data_a_e, data_b, data_b_e)
    mc_string = 'Fit on MC:   (%4.2f #pm %4.2f) + (%4.2f #pm %4.2f) p_{T}^{-1}'%(mc_a,mc_a_e, mc_b, mc_b_e)

    plot_rmue_pt2 = Canvas.Canvas('rmue/%s_%s/plot_rmue_pT2'%(lumi_str, tag), 'png,pdf', 0.6, 0.2, 0.75, 0.35)
    plot_rmue_pt2.addHisto(MCDYControlPT2, 'PE', 'MC', 'PL', r.kRed+1, 1, 0)
    plot_rmue_pt2.addHisto(DATADYControlPT2, 'PE,SAME', 'DATA', 'PL', r.kBlack , 1, 1)
    plot_rmue_pt2.addHisto(fun_da, 'L,SAME', 'Fit - Data', 'L', r.kBlack, 1,2)
    plot_rmue_pt2.addHisto(fun_mc, 'L,SAME', 'Fit - MC', 'L', r.kRed+1, 1,3)
    plot_rmue_pt2.addLatex(0.2,0.2, da_string, 42, 0.02)
    plot_rmue_pt2.addLatex(0.2,0.25, mc_string, 42, 0.02)
    plot_rmue_pt2.save(1, 1, 0, lumi, 0.2, 1.8)




    print "made all the rmues" 
    factorMCDYControlMll =          convertToFactor(MCDYControlMll)
    factorDATADYControlMll =        convertToFactor(DATADYControlMll)
    factorMCDYControlMllvalue =     convertToFactor(MCDYControlMllvalue)
    factorDATADYControlMllvalue =   convertToFactor(DATADYControlMllvalue)
    factorMCDYControlMET =          convertToFactor(MCDYControlMET)
    factorDATADYControlMET =        convertToFactor(DATADYControlMET)
    factorMCDYControlJet =          convertToFactor(MCDYControlJet)
    factorDATADYControlJet =        convertToFactor(DATADYControlJet)
    factorMCDYControlmt2 =          convertToFactor(MCDYControlmt2)
    factorDATADYControlmt2 =        convertToFactor(DATADYControlmt2)

    print "convert all factors" 
    systematicForrmue = 0.1 
    MCrmuemeasured =               MCDYControlMllvalue.GetBinContent(1)
    print "MC rmue ", MCrmuemeasured 
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
    print "ready to save"
    if save==True:
        saveInFile(ingredientsFile, MCrmuemeasured, MCrmuemeasuredUnc, MCrmuemeasuredUncSyst, DATArmuemeasured, DATArmuemeasuredUnc, DATArmuemeasuredUncSyst, 'alone')
        saveInFile(ingredientsFile, factorMCrmuemeasured, factorMCrmuemeasuredUnc, factorMCrmuemeasuredUncSyst, factorDATArmuemeasured, factorDATArmuemeasuredUnc, factorDATArmuemeasuredUncSyst, 'factor')
        saveInFile(ingredientsFile, mc_a, 0., mc_a_e, data_a, 0., data_a_e,  'coeffA')
        saveInFile(ingredientsFile, mc_b, 0., mc_b_e, data_b, 0., data_b_e,  'coeffB')
    MCDYControlPT2.Fit('fun_mc','rQ')
    mc_a = fun_mc.GetParameter(0); mc_a_e = fun_mc.GetParError(0)
    mc_b = fun_mc.GetParameter(1); mc_b_e = fun_mc.GetParError(1)
      
    makeTable(MCDYControlMllmm, MCDYControlMllee, DATADYControlMllmm, DATADYControlMllee, MCrmuemeasured, MCrmuemeasuredUnc, MCrmuemeasuredUncSyst, DATArmuemeasured, DATArmuemeasuredUnc, DATArmuemeasuredUncSyst)
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

    plot_factor_met = Canvas.Canvas('rmue/%s_%s/plot_factor_met'%(lumi_str, tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
    plot_factor_met.addHisto(factorMCDYControlMET, 'PE', 'MC', 'PL', r.kRed+1 , 1, 0)
    plot_factor_met.addHisto(factorDATADYControlMET, 'PE,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_factor_met.addBand(factorMCDYControlMET.GetXaxis().GetXmin(), factorDATArmuemeasured-factorDATArmuemeasuredUncTot, factorMCDYControlMET.GetXaxis().GetXmax(), factorDATArmuemeasured+factorDATArmuemeasuredUncTot, r.kGreen, 0.2)
    plot_factor_met.addLine(factorMCDYControlMET.GetXaxis().GetXmin(), factorDATArmuemeasured, factorMCDYControlMET.GetXaxis().GetXmax(), factorDATArmuemeasured,r.kGreen)
    plot_factor_met.save(1, 1, 0, lumi, 0.2, 1.8)

    plot_factor_jet = Canvas.Canvas('rmue/%s_%s/plot_factor_jet'%(lumi_str, tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
    plot_factor_jet.addHisto(factorMCDYControlJet, 'PE', 'MC', 'PL', r.kRed+1 , 1, 0)
    plot_factor_jet.addHisto(factorDATADYControlJet, 'PE,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_factor_jet.addBand(factorMCDYControlJet.GetXaxis().GetXmin(), factorDATArmuemeasured-factorDATArmuemeasuredUncTot, factorMCDYControlJet.GetXaxis().GetXmax(), factorDATArmuemeasured+factorDATArmuemeasuredUncTot, r.kGreen, 0.2)
    plot_factor_jet.addLine(factorMCDYControlJet.GetXaxis().GetXmin(), factorDATArmuemeasured, factorMCDYControlJet.GetXaxis().GetXmax(), factorDATArmuemeasured,r.kGreen)
    plot_factor_jet.save(1, 1, 0, lumi, 0.2, 1.8)                                                                                                                                                                                                 
    plot_factor_mt2 = Canvas.Canvas('rmue/%s_%s/plot_factor_mt2'%(lumi_str, tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
    plot_factor_mt2.addHisto(factorMCDYControlmt2, 'PE', 'MC', 'PL', r.kRed+1 , 1, 0)
    plot_factor_mt2.addHisto(factorDATADYControlmt2, 'PE,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_factor_mt2.addBand(factorMCDYControlmt2.GetXaxis().GetXmin(), factorDATArmuemeasured-factorDATArmuemeasuredUncTot, factorMCDYControlmt2.GetXaxis().GetXmax(), factorDATArmuemeasured+factorDATArmuemeasuredUncTot, r.kGreen, 0.2)
    plot_factor_mt2.addLine(factorMCDYControlmt2.GetXaxis().GetXmin(), factorDATArmuemeasured, factorMCDYControlmt2.GetXaxis().GetXmax(), factorDATArmuemeasured,r.kGreen)
    plot_factor_mt2.save(1, 1, 0, lumi, 0.2, 1.8)                                                                                                                                                                                                 

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

    mcDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO', 'TTJets_DiLepton', 'ZZTo4L','GGHZZ4L',  'WZTo3LNu', 'WWW', 'WWZ','ZZZ', 'tZq_ll','WWTo2L2Nu', 'ZZTo2L2Nu', 'WZTo2L2Q','TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT', 'TTTT',  'TTWToQQ',  'TTZToLLNuNu' ,'TTWToLNu', 'WJetsToLNu_LO']
 

    daDatasetsB = ['DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376_part1',
                   'DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376_part2',
                   'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part1',
                   'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part2',
                   'MuonEG_Run2016B_23Sep2016_v3_runs_273150_275376', 
                   'DoubleEG_Run2016B_23Sep2016_v3_runs_recovery', 
                   'MuonEG_Run2016B_23Sep2016_v3_runs_recovery']

    daDatasetsC = ['DoubleEG_Run2016C_23Sep2016_v1_runs_271036_284044',
                   'DoubleMuon_Run2016C_23Sep2016_v1_runs_271036_284044',
                   'MuonEG_Run2016C_23Sep2016_v1_runs_271036_284044']

    daDatasetsD = ['DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part2',
                   'DoubleEG_Run2016D_23Sep2016_v1_runs_271036_284044',
                   'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part1',
                   'MuonEG_Run2016D_23Sep2016_v1_runs_271036_284044']

    daDatasetsE = ['DoubleEG_Run2016E_23Sep2016_v1_runs_271036_284044',
                   'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part2',
                   'MuonEG_Run2016E_23Sep2016_v1_runs_271036_284044',
                   'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part1']

    daDatasetsF = ['DoubleEG_Run2016F_23Sep2016_v1_runs_271036_284044',
                   'DoubleMuon_Run2016F_23Sep2016_v1_runs_271036_284044',
                   'MuonEG_Run2016F_23Sep2016_v1_runs_271036_284044']

    daDatasetsG = ['DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part2',
                   'DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part1',
                   'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part3',
                   'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part2',
                   'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part1',
                   'MuonEG_Run2016G_23Sep2016_v1_runs_271036_284044']

    daDatasetsH = ['DoubleEG_Run2016H-PromptReco-v2_runs_281207_284035_part1',
                   'DoubleEG_Run2016H-PromptReco-v2_runs_281207_284035_part2',
                   'DoubleEG_Run2016H-PromptReco-v3_runs_284036_284044',
                   'DoubleMuon_Run2016H-PromptReco-v2_runs_281207_284035_part1',
                   'DoubleMuon_Run2016H-PromptReco-v3_runs_284036_284044',
                   'DoubleMuon_Run2016H-PromptReco-v2_runs_281207_284035_part2',
                   'MuonEG_Run2016H-PromptReco-v3_runs_284036_284044',
                   'DoubleMuon_Run2016H-PromptReco-v2_runs_281207_284035_part3',
                   'MuonEG_Run2016H-PromptReco-v2_runs_281207_284035']                      



    daDatasets = daDatasetsB + daDatasetsC + daDatasetsD +daDatasetsE + daDatasetsF + daDatasetsG + daDatasetsH       


                 

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print bcolors.HEADER + '[rmueAnalysis] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    lumi = 36.4 ; maxrun = 276811; lumi_str = '36.4invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()
    cuts = CutManager.CutManager()

    labelx = "m_{ll} [GeV]"
    labelmet = "E_{T}^{miss} [GeV]"
    labelnjet = "N. Jets"
    labelmt2 = "mt2"
    labelpt2 = "p_{T}^{lep2}"
    #Cuts needed by rmue
    cuts = CutManager.CutManager()

    makeAnalysis(treeDA, treeMC, cuts, '', 'nocut', True, opts.ingredientsFile)


    


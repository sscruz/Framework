#####################################################################
######                                                              #
###### ==========                                   ===========     #  
###### ||                 ||                                ,88     #
###### ||                 ||                              ,88'      #  
###### ||---\     ,adPPYb,|| ,adPPYb,d8  ,adPPYba,      ,88'        #
###### ||---/    a8'    `Y||a8'    `Y88 a8P_____88    ,88'          #
###### ||        8b       ||8b       88 8PP'''''''  ,88'            #
###### ||        '8a,   ,d||'8a,   ,d88 '8b,   ,aa 88'              #
###### ========== `'8bbdP'    `'YbbdP'Y8  `'Ybbd8'' ===========     #
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

def makeResultsTable(rinout_1_mc, rinout_2_mc, rinout_3_mc, rinout_4_mc, rinout_5_mc, rinout_6_mc, rinout_1_da, rinout_2_da, rinout_3_da, rinout_4_da, rinout_5_da, rinout_6_da,lumi, lumi_str):
    
    line0  = '\\begin{tabular}{ccc}'                                                                 
    line1  = '\\hline'                                                                 
    line2  = '\\hline'                                                                 
    line3 = '   $\mathrm{r}_{\mathrm{inout}}$  & '
    line35 = '\\hline'
    line4 = '   MC   & '
    line5 = '   DATA & ' 
    line6 =  '\\hline'
    line7 =  '\\hline'
    line8 = '\\end{{tabular}}' 
    line3 += '\\& $\mathrm{m}_{ll}$ in [20,86] $\mathrm{m}_{ll}$ in [96, 150] $\mathrm{m}_{ll}$ in [150, 200] $\mathrm{m}_{ll}$ in [200,300] $\mathrm{m}_{ll}$ in [300, 400] $\mathrm{m}_{ll}$ above 400    \\hline  \\\\'
    line4 += '  %.3f $\\pm$ %.3f &  %.3f $\\pm$ %.3f &  %.3f $\\pm$ %.3f &  %.3f $\\pm$ %.3f &  %.3f $\\pm$ %.3f &  %.3f $\\pm$ %.3f' %(rinout_1_mc[0], rinout_1_mc[1],rinout_2_mc[0], rinout_2_mc[1], rinout_3_mc[0], rinout_3_mc[1], rinout_4_mc[0], rinout_4_mc[1], rinout_5_mc[0], rinout_5_mc[1], rinout_6_mc[0], rinout_6_mc[1] )
    line5 += '  %.3f $\\pm$ %.3f &  %.3f $\\pm$ %.3f &  %.3f $\\pm$ %.3f &  %.3f $\\pm$ %.3f &  %.3f $\\pm$ %.3f &  %.3f $\\pm$ %.3f' %(rinout_1_da[0], rinout_1_da[1],rinout_2_da[0], rinout_2_da[1], rinout_3_da[0], rinout_3_da[1], rinout_4_da[0], rinout_4_da[1], rinout_5_da[0], rinout_5_da[1], rinout_6_da[0], rinout_6_da[1] )

    helper.ensureDirectory('plots/rinout/%s/'%lumi_str)
    helper.ensureDirectory('plots/rionut/%s/tables/'%lumi_str)
    compTableFile = open('plots/rinout/%s/tables/resultTable_%s%s.tex'%(lumi_str, str(lumi).replace('.','p'), "rinout"),'w')
    compTableFile.write(line0+'\n')
    compTableFile.write(line1+'\n')
    compTableFile.write(line2+'\n')
    compTableFile.write(line3+'\n')
    compTableFile.write(line35+'\n')
    compTableFile.write(line4+'\n')
    compTableFile.write(line5+'\n')
    compTableFile.write(line6+'\n')
    compTableFile.write(line7+'\n')
    compTableFile.write(line8+'\n')
    compTableFile.close()
                                                                                                                                                                        
                                                                                                                                                                         
def readFromFile(theFile, dataMC): 

    for line in open(theFile).readlines():
        if line.find("rsfof") != -1 and line.find("final") != -1:
            if line.find(dataMC) != -1:
                splitline = line.split(" ")
                x = 0
                if dataMC == 'MC': x = 2
                rsfof = splitline[26+x]
                stat = splitline[32+x]
                syst = splitline[38+x]
                arr = []
                arr.append(rsfof)
                arr.append(stat)
                arr.append(syst)
                if arr[0] == '':
                    print "warning probably not right rsfof!!!", arr
                return map(float, arr)
                                                                                    
                                                                                                
def saveInFile(theFile, region, dataMC, rinout, stat, syst):
    foutput = open(theFile + "_aux", 'w')
    for line in open(theFile).readlines():
        if line.find("rinout") != -1 and line.find(dataMC) != -1  and line.find(region) != -1:
            foutput.write('rinout      %s           %s        %.4f      %0.4f       %.4f\n'%(region, dataMC, rinout, stat, syst))
        else:
            foutput.write(line)

    foutput.close()
    subprocess.call(['mv ' + theFile + "_aux " + theFile], shell=True)                                                                                                          
   
                                                                                                                                                                                                                                                          
def calc_frac(inSF, outSF, Ein, Eout):
    val = [0, 0]
    if(inSF != 0):
        val[0] = outSF/inSF
        val[1] = math.sqrt((Eout * Eout) / (inSF * inSF) + (Ein * Ein) * (outSF*outSF) / (inSF * inSF * inSF * inSF))

    return val


def calc_rinout(inSF,outSF, EinSF, EoutSF):

    Ein  = math.sqrt(EinSF*EinSF)
    Eout = math.sqrt(EoutSF*EoutSF)

    return calc_frac(inSF, outSF, Ein, Eout)                                                                                         

def subtractTT(SF, OF, rsfof, dataMC):
    
      corr = SF.Clone("rinout_" + SF.GetName())
      corr.Add(OF , -rsfof)
      print ' i have %.3f events in the SF sample' %(SF.Integral())
      print ' i have %.3f events in the OF sample' %(OF.Integral())
      print ' i have %.3f events in the subtracted sample' %(corr.Integral())
      return corr                                                           


def make_rinout_histo(histo_out ,  histo_in, ymin, ymax):
    print "calculating some rinout histos!"

    ratio = histo_out.Clone("rinout_" + histo_out.GetName())

    # set every bin to the value at the Z
   # for _bin in range(1,finalHisto.GetNbinsX()+1):
   #     finalHisto.SetBinContent(_bin, finalHisto.GetBinContent(finalHisto.FindBin(2.)))
   #     finalHisto.SetBinError  (_bin, finalHisto.GetBinError  (finalHisto.FindBin(2.)))

    ratio.GetYaxis().SetTitle("r_{out/in}")
    ratio.Divide(histo_in)
    ratio.GetYaxis().SetRangeUser(ymin, ymax)
    return ratio                                                                                            



def runAnalysis(lumi, treeDA, treeMC, cuts, specialcut, tag, save, ingredientsFile):

    makeNicePlot = False
    saveInFile = False
    labelx = "m_{ll} [GeV]"
    labelmet = "MET [GeV]"
    labelnjet = "N. Jets"
    inclusivebins  = range(20, 450, 5)
    metbins  = [0, 10, 20, 30, 40, 50, 60, 70, 90]
    njbins  = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.5]
    coarsebins  = [20, 86, 96, 150, 200, 300, 400]

    #####Main mll control plot
    print bcolors.HEADER + '[RinoutAnalysis] ' + bcolors.OKBLUE + 'Starting mll plot' + bcolors.ENDC
    logx = 1
    if makeNicePlot:
       #make the nice finely binned plots for the presentation
        mll_SF_inclusive_mc = treeMC.getTH1F(lumi,"MC_inclusive_mll_SF","lepsMll_Edge",inclusivebins,1,1,cuts.AddList([cuts.goodLepton, cuts.SF, cuts.DYControlRegionNoMll, cuts.MET100]), '', labelx)
        mll_OF_inclusive_mc = treeMC.getTH1F(lumi,"MC_inclusive_mll_OF","lepsMll_Edge",inclusivebins,1,1,cuts.AddList([cuts.goodLepton, cuts.OF, cuts.DYControlRegionNoMll, cuts.MET100]), '', labelx)
        mll_SF_inclusive_da = treeDA.getTH1F(lumi,"DA_inclusive_mll_SF","lepsMll_Edge",inclusivebins,1,1,cuts.AddList([cuts.goodLepton, cuts.SF, cuts.DYControlRegionNoMll, cuts.MET100]), '', labelx)
        mll_OF_inclusive_da = treeDA.getTH1F(lumi,"DA_inclusive_mll_OF","lepsMll_Edge",inclusivebins,1,1,cuts.AddList([cuts.goodLepton, cuts.OF, cuts.DYControlRegionNoMll, cuts.MET100]), '', labelx)
    
        mll_corrected_inclusive_mc = subtractTT(mll_SF_inclusive_mc, mll_OF_inclusive_mc, readFromFile(ingredientsFile, "MC")[0],  "MC")
        mll_corrected_inclusive_da = subtractTT(mll_SF_inclusive_da, mll_OF_inclusive_da, readFromFile(ingredientsFile, "DATA")[0], "DATA")
    
        plot_mll_inclusive_SF = Canvas.Canvas('rinout/%s_%s/plot_mll_inclusive_SF'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)
        plot_mll_inclusive_SF.addHisto(mll_SF_inclusive_mc, 'HIST', 'SF MC', 'L', r.kGreen+1 , 1, 0)
        plot_mll_inclusive_SF.addHisto(mll_SF_inclusive_da, 'E1,SAME', 'SF DATA', 'PL', r.kBlack , 1, 0)
        plot_mll_inclusive_SF.saveRatio(1, 1, logx, lumi, mll_SF_inclusive_da, mll_SF_inclusive_mc)                                                                                
       
        plot_mll_inclusive_OF = Canvas.Canvas('rinout/%s_%s/plot_mll_inclusive_OF'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)
        plot_mll_inclusive_OF.addHisto(mll_OF_inclusive_mc, 'HIST', 'OF MC', 'L', r.kGreen+1 , 1, 0)
        plot_mll_inclusive_OF.addHisto(mll_OF_inclusive_da, 'E1,SAME', 'OF DATA', 'PL', r.kBlack , 1, 0)
        plot_mll_inclusive_OF.saveRatio(1, 1, logx, lumi, mll_SF_inclusive_da, mll_SF_inclusive_mc)                                                                                
    
        plot_mll_inclusive_corrected = Canvas.Canvas('rinout/%s_%s/plot_mll_inclusive_corrected'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)
        plot_mll_inclusive_corrected.addHisto(mll_corrected_inclusive_mc, 'HIST', 'SF-OF MC', 'L', r.kGreen+1 , 1, 0)
        plot_mll_inclusive_corrected.addHisto(mll_corrected_inclusive_da, 'E1,SAME', 'SF-OF DATA', 'PL', r.kBlack , 1, 0)
        plot_mll_inclusive_corrected.saveRatio(1, 1, logx, lumi, mll_corrected_inclusive_da, mll_corrected_inclusive_mc)                                                                                

    #make the histos for the rinout calculationi
    mll_SF_coarse_mc = treeMC.getTH1F(lumi, "MC_inclusive_mll_SF", "lepsMll_Edge", coarsebins, 1, 1, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.DYControlRegionNoMll,cuts.MET100]), '', labelx)
    mll_OF_coarse_mc = treeMC.getTH1F(lumi, "MC_inclusive_mll_OF", "lepsMll_Edge", coarsebins, 1, 1, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.DYControlRegionNoMll,cuts.MET100]), '', labelx)
    mll_SF_coarse_da = treeDA.getTH1F(lumi, "DA_inclusive_mll_SF", "lepsMll_Edge", coarsebins, 1, 1, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.DYControlRegionNoMll,cuts.MET100]), '', labelx)
    mll_OF_coarse_da = treeDA.getTH1F(lumi, "DA_inclusive_mll_OF", "lepsMll_Edge", coarsebins, 1, 1, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.DYControlRegionNoMll,cuts.MET100]), '', labelx)
    print "mll_SF_coarse_mc" , mll_SF_coarse_mc.Integral()
    print "mll_SF_coarse_da" , mll_SF_coarse_da.Integral()
    mll_corrected_mc = subtractTT(mll_SF_coarse_mc, mll_OF_coarse_mc, readFromFile(ingredientsFile, "MC")[0],  "MC")
    mll_corrected_da = subtractTT(mll_SF_coarse_da, mll_OF_coarse_da, readFromFile(ingredientsFile, "DATA")[0], "DATA")
    print "mll_corrected_mc", mll_corrected_mc.Integral()                                                                                                                     
    print "mll_corrected_da", mll_corrected_da.Integral()                                                                                                                     
    plot_mll_coarse_SF = Canvas.Canvas('rinout/%s_%s/plot_mll_SF'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)
    plot_mll_coarse_SF.addHisto(mll_SF_coarse_mc, 'HIST', 'SF MC', 'L', r.kGreen+1 , 1, 0)
    plot_mll_coarse_SF.addHisto(mll_SF_coarse_da, 'E1,SAME', 'SF DATA', 'PL', r.kBlack , 1, 0)
    plot_mll_coarse_SF.saveRatio(1, 1, logx, lumi, mll_SF_coarse_da, mll_SF_coarse_mc)                                                                                
    plot_mll_coarse_OF = Canvas.Canvas('rinout/%s_%s/plot_mll_OF'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)
    plot_mll_coarse_OF.addHisto(mll_OF_coarse_mc, 'HIST', 'OF MC', 'L', r.kGreen+1 , 1, 0)
    plot_mll_coarse_OF.addHisto(mll_OF_coarse_da, 'E1,SAME', 'OF DATA', 'PL', r.kBlack , 1, 0)
    plot_mll_coarse_OF.saveRatio(1, 1, logx, lumi, mll_SF_coarse_da, mll_SF_coarse_mc)                                                                                
    plot_mll_coarse_corrected = Canvas.Canvas('rinout/%s_%s/plot_mll_corrected'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)
    plot_mll_coarse_corrected.addHisto(mll_corrected_mc, 'HIST', 'SF-OF MC', 'L', r.kGreen+1 , 1, 0)
    plot_mll_coarse_corrected.addHisto(mll_corrected_da, 'E1,SAME', 'SF-OF DATA', 'PL', r.kBlack , 1, 0)
    plot_mll_coarse_corrected.saveRatio(1, 1, logx, lumi, mll_corrected_da, mll_corrected_mc)                                                                                

    rinout_1_mc  = calc_rinout(mll_corrected_mc.GetBinContent(2), mll_corrected_mc.GetBinContent(1) , mll_corrected_mc.GetBinError(2), mll_corrected_mc.GetBinError(1))
    rinout_z_mc  = calc_rinout(mll_corrected_mc.GetBinContent(2), mll_corrected_mc.GetBinContent(2) , mll_corrected_mc.GetBinError(2), mll_corrected_mc.GetBinError(2))
    rinout_2_mc  = calc_rinout(mll_corrected_mc.GetBinContent(2), mll_corrected_mc.GetBinContent(3) , mll_corrected_mc.GetBinError(2), mll_corrected_mc.GetBinError(3))
    rinout_3_mc  = calc_rinout(mll_corrected_mc.GetBinContent(2), mll_corrected_mc.GetBinContent(4) , mll_corrected_mc.GetBinError(2), mll_corrected_mc.GetBinError(4))
    rinout_4_mc  = calc_rinout(mll_corrected_mc.GetBinContent(2), mll_corrected_mc.GetBinContent(5) , mll_corrected_mc.GetBinError(2), mll_corrected_mc.GetBinError(5))
    rinout_5_mc  = calc_rinout(mll_corrected_mc.GetBinContent(2), mll_corrected_mc.GetBinContent(6) , mll_corrected_mc.GetBinError(2), mll_corrected_mc.GetBinError(6))
    rinout_6_mc  = calc_rinout(mll_corrected_mc.GetBinContent(2), mll_corrected_mc.GetBinContent(7) , mll_corrected_mc.GetBinError(2), mll_corrected_mc.GetBinError(7))      
    rinout_1_da  = calc_rinout(mll_corrected_da.GetBinContent(2), mll_corrected_da.GetBinContent(1) , mll_corrected_da.GetBinError(2), mll_corrected_da.GetBinError(1))
    rinout_z_da  = calc_rinout(mll_corrected_da.GetBinContent(2), mll_corrected_da.GetBinContent(2) , mll_corrected_da.GetBinError(2), mll_corrected_da.GetBinError(2))
    rinout_2_da  = calc_rinout(mll_corrected_da.GetBinContent(2), mll_corrected_da.GetBinContent(3) , mll_corrected_da.GetBinError(2), mll_corrected_da.GetBinError(3))
    rinout_3_da  = calc_rinout(mll_corrected_da.GetBinContent(2), mll_corrected_da.GetBinContent(4) , mll_corrected_da.GetBinError(2), mll_corrected_da.GetBinError(4))
    rinout_4_da  = calc_rinout(mll_corrected_da.GetBinContent(2), mll_corrected_da.GetBinContent(5) , mll_corrected_da.GetBinError(2), mll_corrected_da.GetBinError(5))
    rinout_5_da  = calc_rinout(mll_corrected_da.GetBinContent(2), mll_corrected_da.GetBinContent(6) , mll_corrected_da.GetBinError(2), mll_corrected_da.GetBinError(6))
    rinout_6_da  = calc_rinout(mll_corrected_da.GetBinContent(2), mll_corrected_da.GetBinContent(7) , mll_corrected_da.GetBinError(2), mll_corrected_da.GetBinError(7))      
    
    print 'rinout_1_da ', rinout_1_da
    print 'rinout_1_mc ', rinout_1_mc
    print 'rinout_2_da ', rinout_2_da
    print 'rinout_2_mc ', rinout_2_mc
    print 'rinout_3_da ', rinout_3_da
    print 'rinout_3_mc ', rinout_3_mc
    print 'rinout_4_da ', rinout_4_da
    print 'rinout_4_mc ', rinout_4_mc
    print 'rinout_5_da ', rinout_5_da
    print 'rinout_5_mc ', rinout_5_mc
    print 'rinout_6_da ', rinout_6_da
    print 'rinout_6_mc ', rinout_6_mc    
    if saveInFile:
        makeResultsTable(rinout_1_mc, rinout_2_mc , rinout_3_mc, rinout_4_mc, rinout_5_mc, rinout_6_mc,rinout_1_da, rinout_2_da , rinout_3_da, rinout_4_da, rinout_5_da, rinout_6_da, lumi, str(lumi))
        saveInFile(ingredientsFile, 'dy_m1', 'DATA', rinout_1_da[0], rinout_1_da[1], 0.2)
        saveInFile(ingredientsFile, 'dy_oz', 'DATA', rinout_z_da[0], rinout_z_da[1], 0.2)
        saveInFile(ingredientsFile, 'dy_m2', 'DATA', rinout_2_da[0], rinout_2_da[1], 0.2)
        saveInFile(ingredientsFile, 'dy_m3', 'DATA', rinout_3_da[0], rinout_3_da[1], 0.2)
        saveInFile(ingredientsFile, 'dy_m4', 'DATA', rinout_4_da[0], rinout_4_da[1], 0.2)
        saveInFile(ingredientsFile, 'dy_m5', 'DATA', rinout_5_da[0], rinout_5_da[1], 0.2)
        saveInFile(ingredientsFile, 'dy_m6', 'DATA', rinout_6_da[0], rinout_6_da[1], 0.2)
        saveInFile(ingredientsFile, 'dy_m1', 'MC  ', rinout_1_mc[0], rinout_1_mc[1], 0.2)
        saveInFile(ingredientsFile, 'dy_oz', 'MC  ', rinout_z_mc[0], rinout_z_mc[1], 0.2)
        saveInFile(ingredientsFile, 'dy_m2', 'MC  ', rinout_2_mc[0], rinout_2_mc[1], 0.2)
        saveInFile(ingredientsFile, 'dy_m3', 'MC  ', rinout_3_mc[0], rinout_3_mc[1], 0.2)
        saveInFile(ingredientsFile, 'dy_m4', 'MC  ', rinout_4_mc[0], rinout_4_mc[1], 0.2)
        saveInFile(ingredientsFile, 'dy_m5', 'MC  ', rinout_5_mc[0], rinout_5_mc[1], 0.2)
        saveInFile(ingredientsFile, 'dy_m6', 'MC  ', rinout_6_mc[0], rinout_6_mc[1], 0.2)
        print "done writing in the ingredients.dat"

    met_1_SF_mc= treeMC.getTH1F(lumi, "MC_met_mll_20-86_SF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll20_86]),'', labelmet)
    met_1_SF_da= treeDA.getTH1F(lumi, "DA_met_mll_20-86_SF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll20_86]),'', labelmet)
    met_1_OF_mc= treeMC.getTH1F(lumi, "MC_met_mll_20-86_OF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll20_86]),'', labelmet)
    met_1_OF_da= treeDA.getTH1F(lumi, "DA_met_mll_20-86_OF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll20_86]),'', labelmet)
    met_z_SF_mc= treeMC.getTH1F(lumi, "MC_met_mll_86-96_SF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll86_96]),'', labelmet)
    met_z_SF_da= treeDA.getTH1F(lumi, "DA_met_mll_86-96_SF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll86_96]),'', labelmet)
    met_z_OF_mc= treeMC.getTH1F(lumi, "MC_met_mll_86-96_OF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll86_96]),'', labelmet)
    met_z_OF_da= treeDA.getTH1F(lumi, "DA_met_mll_86-96_OF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll86_96]),'', labelmet)
    met_2_SF_mc= treeMC.getTH1F(lumi, "MC_met_mll_96-150_SF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll96_150]),'', labelmet)
    met_2_SF_da= treeDA.getTH1F(lumi, "DA_met_mll_96-150_SF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll96_150]),'', labelmet)
    met_2_OF_mc= treeMC.getTH1F(lumi, "MC_met_mll_96-150_OF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll96_150]),'', labelmet)
    met_2_OF_da= treeDA.getTH1F(lumi, "DA_met_mll_96-150_OF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll96_150]),'', labelmet)
    met_3_SF_mc= treeMC.getTH1F(lumi, "MC_met_mll_150-200_SF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll150_200]),'', labelmet)
    met_3_SF_da= treeDA.getTH1F(lumi, "DA_met_mll_150-200_SF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll150_200]),'', labelmet)
    met_3_OF_mc= treeMC.getTH1F(lumi, "MC_met_mll_150-200_OF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll150_200]),'', labelmet)
    met_3_OF_da= treeDA.getTH1F(lumi, "DA_met_mll_150-200_OF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll150_200]),'', labelmet)
    met_4_SF_mc= treeMC.getTH1F(lumi, "MC_met_mll_200-300_SF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll200_300]),'', labelmet)
    met_4_SF_da= treeDA.getTH1F(lumi, "DA_met_mll_200-300_SF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll200_300]),'', labelmet)
    met_4_OF_mc= treeMC.getTH1F(lumi, "MC_met_mll_200-300_OF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll200_300]),'', labelmet)
    met_4_OF_da= treeDA.getTH1F(lumi, "DA_met_mll_200-300_OF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll200_300]),'', labelmet)
    met_5_SF_mc= treeMC.getTH1F(lumi, "MC_met_mll_300-400_SF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll300_400]),'', labelmet)
    met_5_SF_da= treeDA.getTH1F(lumi, "DA_met_mll_300-400_SF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll300_400]),'', labelmet)
    met_5_OF_mc= treeMC.getTH1F(lumi, "MC_met_mll_300-400_OF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll300_400]),'', labelmet)
    met_5_OF_da= treeDA.getTH1F(lumi, "DA_met_mll_300-400_OF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll300_400]),'', labelmet)
    met_6_SF_mc= treeMC.getTH1F(lumi, "MC_met_mll_400_SF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll400]),'', labelmet)
    met_6_SF_da= treeDA.getTH1F(lumi, "DA_met_mll_400_SF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll400]),'', labelmet)
    met_6_OF_mc= treeMC.getTH1F(lumi, "MC_met_mll_400_OF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll400]),'', labelmet)
    met_6_OF_da= treeDA.getTH1F(lumi, "DA_met_mll_400_OF","met_Edge",metbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.MET100, cuts.DYControlRegionNoMETNoMll,cuts.mll400]),'', labelmet)
    print "met_1_SF_mc", met_1_SF_mc.Integral()
    print "met_1_SF_da", met_1_SF_da.Integral()
    #subtract the OF from the total
    met_1_corrected_mc  = subtractTT(met_1_SF_mc, met_1_OF_mc, readFromFile(ingredientsFile, "MC")[0],  "MC")
    met_1_corrected_da  = subtractTT(met_1_SF_da, met_1_OF_da, readFromFile(ingredientsFile, "DATA")[0],  "DATA")
    met_z_corrected_mc  = subtractTT(met_z_SF_mc, met_z_OF_mc, readFromFile(ingredientsFile, "MC")[0],  "MC")
    met_z_corrected_da  = subtractTT(met_z_SF_da, met_z_OF_da, readFromFile(ingredientsFile, "DATA")[0],  "DATA")
    met_2_corrected_mc  = subtractTT(met_2_SF_mc, met_2_OF_mc, readFromFile(ingredientsFile, "MC")[0],  "MC")
    met_2_corrected_da  = subtractTT(met_2_SF_da, met_2_OF_da, readFromFile(ingredientsFile, "DATA")[0],  "DATA")
    met_3_corrected_mc  = subtractTT(met_3_SF_mc, met_3_OF_mc, readFromFile(ingredientsFile, "MC")[0],  "MC")
    met_3_corrected_da  = subtractTT(met_3_SF_da, met_3_OF_da, readFromFile(ingredientsFile, "DATA")[0],  "DATA")
    met_4_corrected_mc  = subtractTT(met_4_SF_mc, met_4_OF_mc, readFromFile(ingredientsFile, "MC")[0],  "MC")
    met_4_corrected_da  = subtractTT(met_4_SF_da, met_4_OF_da, readFromFile(ingredientsFile, "DATA")[0],  "DATA")
    met_5_corrected_mc  = subtractTT(met_5_SF_mc, met_5_OF_mc, readFromFile(ingredientsFile, "MC")[0],  "MC")
    met_5_corrected_da  = subtractTT(met_5_SF_da, met_5_OF_da, readFromFile(ingredientsFile, "DATA")[0],  "DATA")
    met_6_corrected_mc  = subtractTT(met_6_SF_mc, met_6_OF_mc, readFromFile(ingredientsFile, "MC")[0],  "MC")
    met_6_corrected_da  = subtractTT(met_6_SF_da, met_6_OF_da, readFromFile(ingredientsFile, "DATA")[0],  "DATA")
   
    print "met_1_corrected_mc", met_1_corrected_mc.Integral()
    print "met_1_corrected_da", met_1_corrected_da.Integral()
    #calculate rinout with the corrected (OF subtracted) samples
    met_rinout_1_mc = make_rinout_histo(met_1_corrected_mc, met_z_corrected_mc, 0., 1. )
    met_rinout_1_da = make_rinout_histo(met_1_corrected_da, met_z_corrected_da, 0., 1. )  
    met_rinout_2_mc = make_rinout_histo(met_2_corrected_mc, met_z_corrected_mc, 0., 0.5 )
    met_rinout_2_da = make_rinout_histo(met_2_corrected_da, met_z_corrected_da, 0., 0.5 )
    met_rinout_3_mc = make_rinout_histo(met_3_corrected_mc, met_z_corrected_mc, 0., 0.05) 
    met_rinout_3_da = make_rinout_histo(met_3_corrected_da, met_z_corrected_da, 0., 0.05)
    met_rinout_4_mc = make_rinout_histo(met_4_corrected_mc, met_z_corrected_mc, 0., 0.05)  
    met_rinout_4_da = make_rinout_histo(met_4_corrected_da, met_z_corrected_da, 0., 0.05)
    met_rinout_5_mc = make_rinout_histo(met_5_corrected_mc, met_z_corrected_mc, 0., 0.05) 
    met_rinout_5_da = make_rinout_histo(met_5_corrected_da, met_z_corrected_da, 0., 0.05)
    met_rinout_6_mc = make_rinout_histo(met_6_corrected_mc, met_z_corrected_mc, 0., 0.05) 
    met_rinout_6_da = make_rinout_histo(met_6_corrected_da, met_z_corrected_da, 0., 0.05)
    print "met_rinout_1_mc", met_rinout_1_mc.Integral()
    print "met_rinout_1_da", met_rinout_1_da.Integral()

    plot_met_rinout_1 = Canvas.Canvas('rinout/%s_%s/plot_met_rinout_mll_20-86'%(lumi_str,tag), 'png,pdf', 0.6, 0.7, 0.85, 0.9)
    plot_met_rinout_1.addHisto(met_rinout_1_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_met_rinout_1.addHisto(met_rinout_1_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_met_rinout_1.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 20-86"%(rinout_1_da[0]) )
    plot_met_rinout_1.addBand(met_rinout_1_mc.GetXaxis().GetXmin(), rinout_1_da[0]*0.8, met_rinout_1_mc.GetXaxis().GetXmax(), rinout_1_da[0]*1.2, r.kOrange+6, 0.2)
    plot_met_rinout_1.addLine(met_rinout_1_mc.GetXaxis().GetXmin(), rinout_1_da[0], met_rinout_1_mc.GetXaxis().GetXmax(), rinout_1_da[0] ,r.kBlue)
    plot_met_rinout_1.saveRatio(1, 1, 0, lumi, met_rinout_1_da, met_rinout_1_mc)                                                                                                  
    plot_met_rinout_2 = Canvas.Canvas('rinout/%s_%s/plot_met_rinout_mll_86-150'%(lumi_str,tag), 'png,pdf', 0.6, 0.7, 0.85, 0.9)
    plot_met_rinout_2.addHisto(met_rinout_2_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_met_rinout_2.addHisto(met_rinout_2_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_met_rinout_2.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 86-150"%(rinout_2_da[0]) )
    plot_met_rinout_2.addBand(met_rinout_2_mc.GetXaxis().GetXmin(), rinout_2_da[0]*0.8, met_rinout_2_mc.GetXaxis().GetXmax(), rinout_2_da[0]*1.2, r.kOrange+6, 0.2)
    plot_met_rinout_2.addLine(met_rinout_2_mc.GetXaxis().GetXmin(), rinout_2_da[0], met_rinout_2_mc.GetXaxis().GetXmax(), rinout_2_da[0] ,r.kBlue)
    plot_met_rinout_2.saveRatio(1, 1, 0, lumi, met_rinout_2_da, met_rinout_2_mc)                                                        
    plot_met_rinout_3 = Canvas.Canvas('rinout/%s_%s/plot_met_rinout_mll_150-200'%(lumi_str,tag), 'png,pdf', 0.6, 0.7, 0.85, 0.9)           
    plot_met_rinout_3.addHisto(met_rinout_3_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_met_rinout_3.addHisto(met_rinout_3_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_met_rinout_3.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 150-200"%(rinout_3_da[0]) )
    plot_met_rinout_3.addBand(met_rinout_3_mc.GetXaxis().GetXmin(), rinout_3_da[0]*0.8, met_rinout_3_mc.GetXaxis().GetXmax(), rinout_3_da[0]*1.2, r.kOrange+6, 0.2)
    plot_met_rinout_3.addLine(met_rinout_3_mc.GetXaxis().GetXmin(), rinout_3_da[0], met_rinout_3_mc.GetXaxis().GetXmax(), rinout_3_da[0] ,r.kBlue)
    plot_met_rinout_3.saveRatio(1, 1, 0, lumi, met_rinout_3_da, met_rinout_3_mc)                                                        
    plot_met_rinout_4 = Canvas.Canvas('rinout/%s_%s/plot_met_rinout_mll_200-300'%(lumi_str,tag), 'png,pdf', 0.6, 0.7, 0.85, 0.9)           
    plot_met_rinout_4.addHisto(met_rinout_4_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_met_rinout_4.addHisto(met_rinout_4_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_met_rinout_4.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 200-300"%(rinout_4_da[0]) )
    plot_met_rinout_4.addBand(met_rinout_4_mc.GetXaxis().GetXmin(), rinout_4_da[0]*0.8, met_rinout_4_mc.GetXaxis().GetXmax(), rinout_4_da[0]*1.2, r.kOrange+6, 0.2)
    plot_met_rinout_4.addLine(met_rinout_4_mc.GetXaxis().GetXmin(), rinout_4_da[0], met_rinout_4_mc.GetXaxis().GetXmax(), rinout_4_da[0] ,r.kBlue)
    plot_met_rinout_4.saveRatio(1, 1, 0, lumi, met_rinout_4_da, met_rinout_4_mc)                                                        
    plot_met_rinout_5 = Canvas.Canvas('rinout/%s_%s/plot_met_rinout_mll_300-400'%(lumi_str,tag), 'png,pdf', 0.6, 0.7, 0.85, 0.9)           
    plot_met_rinout_5.addHisto(met_rinout_5_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_met_rinout_5.addHisto(met_rinout_5_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_met_rinout_5.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 300-400"%(rinout_5_da[0]) )
    plot_met_rinout_5.addBand(met_rinout_5_mc.GetXaxis().GetXmin(), rinout_5_da[0]*0.8, met_rinout_5_mc.GetXaxis().GetXmax(), rinout_5_da[0]*1.2, r.kOrange+6, 0.2)
    plot_met_rinout_5.addLine(met_rinout_5_mc.GetXaxis().GetXmin(), rinout_5_da[0], met_rinout_5_mc.GetXaxis().GetXmax(), rinout_5_da[0] ,r.kBlue)
    plot_met_rinout_5.saveRatio(1, 1, 0, lumi, met_rinout_5_da, met_rinout_5_mc)                                                        
    plot_met_rinout_6 = Canvas.Canvas('rinout/%s_%s/plot_met_rinout_mll_400'%(lumi_str,tag), 'png,pdf', 0.6, 0.7, 0.85, 0.9)       
    plot_met_rinout_6.addHisto(met_rinout_6_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_met_rinout_6.addHisto(met_rinout_6_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_met_rinout_6.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} > 400"%(rinout_6_da[0]) )
    plot_met_rinout_6.addBand(met_rinout_6_mc.GetXaxis().GetXmin(), rinout_6_da[0]*0.8, met_rinout_6_mc.GetXaxis().GetXmax(), rinout_6_da[0]*1.2, r.kOrange+6, 0.2)
    plot_met_rinout_6.addLine(met_rinout_6_mc.GetXaxis().GetXmin(), rinout_6_da[0], met_rinout_6_mc.GetXaxis().GetXmax(), rinout_6_da[0] ,r.kBlue)
    plot_met_rinout_6.saveRatio(1, 1, 0, lumi, met_rinout_6_da, met_rinout_6_mc)                                                                                 
    
    #make njet plot##############################################################################################################################################################################
    nj_1_SF_mc=treeMC.getTH1F(lumi,"M_njet_mll_20-86_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF, cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll20_86]),'', labelnjet) 
    nj_1_SF_da=treeDA.getTH1F(lumi,"D_njet_mll_20-86_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF, cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll20_86]),'', labelnjet)
    nj_1_OF_mc=treeMC.getTH1F(lumi,"M_njet_mll_20-86_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF, cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll20_86]),'', labelnjet)
    nj_1_OF_da=treeDA.getTH1F(lumi,"D_njet_mll_20-86_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF, cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll20_86]),'', labelnjet)
    nj_z_SF_mc=treeMC.getTH1F(lumi,"M_njet_mll_86-96_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF, cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll86_96]),'', labelnjet)
    nj_z_SF_da=treeDA.getTH1F(lumi,"D_njet_mll_86-96_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF, cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll86_96]),'', labelnjet)
    nj_z_OF_mc=treeMC.getTH1F(lumi,"M_njet_mll_86-96_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF, cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll86_96]),'', labelnjet)
    nj_z_OF_da=treeDA.getTH1F(lumi,"D_njet_mll_86-96_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF, cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll86_96]),'', labelnjet)
    nj_2_SF_mc=treeMC.getTH1F(lumi,"M_njet_mll_96-150_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF, cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll96_150]),'',labelnjet)
    nj_2_SF_da=treeDA.getTH1F(lumi,"D_njet_mll_96-150_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF, cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll96_150]),'',labelnjet)
    nj_2_OF_mc=treeMC.getTH1F(lumi,"M_njet_mll_96-150_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF, cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll96_150]),'',labelnjet)
    nj_2_OF_da=treeDA.getTH1F(lumi,"D_njet_mll_96-150_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF, cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll96_150]),'',labelnjet)
    nj_3_SF_mc=treeMC.getTH1F(lumi,"M_njet_mll_150-200_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll150_200]),'',labelnjet)
    nj_3_SF_da=treeDA.getTH1F(lumi,"D_njet_mll_150-200_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll150_200]),'',labelnjet)
    nj_3_OF_mc=treeMC.getTH1F(lumi,"M_njet_mll_150-200_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll150_200]),'',labelnjet)
    nj_3_OF_da=treeDA.getTH1F(lumi,"D_njet_mll_150-200_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll150_200]),'',labelnjet)
    nj_4_SF_mc=treeMC.getTH1F(lumi,"M_njet_mll_200-300_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll200_300]),'',labelnjet)
    nj_4_SF_da=treeDA.getTH1F(lumi,"D_njet_mll_200-300_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll200_300]),'',labelnjet)
    nj_4_OF_mc=treeMC.getTH1F(lumi,"M_njet_mll_200-300_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll200_300]),'',labelnjet)
    nj_4_OF_da=treeDA.getTH1F(lumi,"D_njet_mll_200-300_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll200_300]),'',labelnjet)
    nj_5_SF_mc=treeMC.getTH1F(lumi,"M_njet_mll_300-400_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll300_400]),'',labelnjet)
    nj_5_SF_da=treeDA.getTH1F(lumi,"D_njet_mll_300-400_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF,cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll300_400]),'',labelnjet)
    nj_5_OF_mc=treeMC.getTH1F(lumi,"M_njet_mll_300-400_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll300_400]),'',labelnjet)
    nj_5_OF_da=treeDA.getTH1F(lumi,"D_njet_mll_300-400_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF,cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll300_400]),'',labelnjet)
    nj_6_SF_mc=treeMC.getTH1F(lumi,"M_njet_mll_400_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF, cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll400]),'',labelnjet)
    nj_6_SF_da=treeDA.getTH1F(lumi,"D_njet_mll_400_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.SF, cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll400]),'',labelnjet)
    nj_6_OF_mc=treeMC.getTH1F(lumi,"M_njet_mll_400_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF, cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll400]),'',labelnjet)
    nj_6_OF_da=treeDA.getTH1F(lumi,"D_njet_mll_400_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.goodLepton,cuts.OF, cuts.DYControlRegionNoNjetNoMll,cuts.MET100,cuts.mll400]),'',labelnjet)
                                                                                                                             
    nj_1_corrected_mc  = subtractTT(nj_1_SF_mc, nj_1_OF_mc, readFromFile(ingredientsFile, "MC")[0],  "MC")
    nj_1_corrected_da  = subtractTT(nj_1_SF_da, nj_1_OF_da, readFromFile(ingredientsFile, "DATA")[0],  "DATA")
    nj_z_corrected_mc  = subtractTT(nj_z_SF_mc, nj_z_OF_mc, readFromFile(ingredientsFile, "MC")[0],  "MC")
    nj_z_corrected_da  = subtractTT(nj_z_SF_da, nj_z_OF_da, readFromFile(ingredientsFile, "DATA")[0],  "DATA")
    nj_2_corrected_mc  = subtractTT(nj_2_SF_mc, nj_2_OF_mc, readFromFile(ingredientsFile, "MC")[0],  "MC")
    nj_2_corrected_da  = subtractTT(nj_2_SF_da, nj_2_OF_da, readFromFile(ingredientsFile, "DATA")[0],  "DATA")
    nj_3_corrected_mc  = subtractTT(nj_3_SF_mc, nj_3_OF_mc, readFromFile(ingredientsFile, "MC")[0],  "MC")
    nj_3_corrected_da  = subtractTT(nj_3_SF_da, nj_3_OF_da, readFromFile(ingredientsFile, "DATA")[0],  "DATA")
    nj_4_corrected_mc  = subtractTT(nj_4_SF_mc, nj_4_OF_mc, readFromFile(ingredientsFile, "MC")[0],  "MC")
    nj_4_corrected_da  = subtractTT(nj_4_SF_da, nj_4_OF_da, readFromFile(ingredientsFile, "DATA")[0],  "DATA")
    nj_5_corrected_mc  = subtractTT(nj_5_SF_mc, nj_5_OF_mc, readFromFile(ingredientsFile, "MC")[0],  "MC")
    nj_5_corrected_da  = subtractTT(nj_5_SF_da, nj_5_OF_da, readFromFile(ingredientsFile, "DATA")[0],  "DATA")
    nj_6_corrected_mc  = subtractTT(nj_6_SF_mc, nj_6_OF_mc, readFromFile(ingredientsFile, "MC")[0],  "MC")
    nj_6_corrected_da  = subtractTT(nj_6_SF_da, nj_6_OF_da, readFromFile(ingredientsFile, "DATA")[0],  "DATA")
    
    nj_rinout_1_mc = make_rinout_histo(nj_1_corrected_mc, nj_z_corrected_mc, 0., 1. )
    nj_rinout_1_da = make_rinout_histo(nj_1_corrected_da, nj_z_corrected_da, 0., 1. )  
    nj_rinout_2_mc = make_rinout_histo(nj_2_corrected_mc, nj_z_corrected_mc, 0., 0.5 )
    nj_rinout_2_da = make_rinout_histo(nj_2_corrected_da, nj_z_corrected_da, 0., 0.5 )
    nj_rinout_3_mc = make_rinout_histo(nj_3_corrected_mc, nj_z_corrected_mc, 0., 0.05 ) 
    nj_rinout_3_da = make_rinout_histo(nj_3_corrected_da, nj_z_corrected_da, 0., 0.05 )
    nj_rinout_4_mc = make_rinout_histo(nj_4_corrected_mc, nj_z_corrected_mc, 0., 0.05 )  
    nj_rinout_4_da = make_rinout_histo(nj_4_corrected_da, nj_z_corrected_da, 0., 0.05 )
    nj_rinout_5_mc = make_rinout_histo(nj_5_corrected_mc, nj_z_corrected_mc, 0., 0.05 ) 
    nj_rinout_5_da = make_rinout_histo(nj_5_corrected_da, nj_z_corrected_da, 0., 0.05 )
    nj_rinout_6_mc = make_rinout_histo(nj_6_corrected_mc, nj_z_corrected_mc, 0., 0.05 ) 
    nj_rinout_6_da = make_rinout_histo(nj_6_corrected_da, nj_z_corrected_da, 0., 0.05 )
                                                                                                                                                                                       
    plot_nj_rinout_1 = Canvas.Canvas('rinout/%s_%s/plot_njet_rinout_mll_20-86'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)
    plot_nj_rinout_1.addHisto(nj_rinout_1_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_nj_rinout_1.addHisto(nj_rinout_1_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_nj_rinout_1.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 20-86"%(rinout_1_da[0]) )
    plot_nj_rinout_1.addBand(nj_rinout_1_mc.GetXaxis().GetXmin(), rinout_1_da[0]*0.8, nj_rinout_1_mc.GetXaxis().GetXmax(), rinout_1_da[0]*1.2, r.kOrange+6, 0.2)
    plot_nj_rinout_1.addLine(nj_rinout_1_mc.GetXaxis().GetXmin(), rinout_1_da[0], nj_rinout_1_mc.GetXaxis().GetXmax(), rinout_1_da[0] ,r.kBlue)
    plot_nj_rinout_1.saveRatio(1, 1, 0, lumi, nj_rinout_1_da, nj_rinout_1_mc)                                                                                                  
    plot_nj_rinout_2 = Canvas.Canvas('rinout/%s_%s/plot_njet_rinout_mll_86-150'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)
    plot_nj_rinout_2.addHisto(nj_rinout_2_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_nj_rinout_2.addHisto(nj_rinout_2_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_nj_rinout_2.addBand(nj_rinout_2_mc.GetXaxis().GetXmin(), rinout_2_da[0]*0.8, nj_rinout_2_mc.GetXaxis().GetXmax(), rinout_2_da[0]*1.2, r.kOrange+6, 0.2)
    plot_nj_rinout_2.addLine(nj_rinout_2_mc.GetXaxis().GetXmin(), rinout_2_da[0], nj_rinout_2_mc.GetXaxis().GetXmax(), rinout_2_da[0] ,r.kBlue)
    plot_nj_rinout_2.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 86-150"%(rinout_2_da[0]) )
    plot_nj_rinout_2.saveRatio(1, 1, 0, lumi, nj_rinout_2_da, nj_rinout_2_mc)                                                        
    plot_nj_rinout_3 = Canvas.Canvas('rinout/%s_%s/plot_njet_rinout_mll_150-200'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)           
    plot_nj_rinout_3.addHisto(nj_rinout_3_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_nj_rinout_3.addHisto(nj_rinout_3_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_nj_rinout_3.addBand(nj_rinout_3_mc.GetXaxis().GetXmin(), rinout_3_da[0]*0.8, nj_rinout_3_mc.GetXaxis().GetXmax(), rinout_3_da[0]*1.2, r.kOrange+6, 0.2)
    plot_nj_rinout_3.addLine(nj_rinout_3_mc.GetXaxis().GetXmin(), rinout_3_da[0], nj_rinout_3_mc.GetXaxis().GetXmax(), rinout_3_da[0] ,r.kBlue) 
    plot_nj_rinout_3.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 150-200"%(rinout_3_da[0]) )
    plot_nj_rinout_3.saveRatio(1, 1, 0, lumi, nj_rinout_3_da, nj_rinout_3_mc)                                                        
    plot_nj_rinout_4 = Canvas.Canvas('rinout/%s_%s/plot_njet_rinout_mll_200-300'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)           
    plot_nj_rinout_4.addHisto(nj_rinout_4_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_nj_rinout_4.addHisto(nj_rinout_4_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_nj_rinout_4.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 200-300"%(rinout_4_da[0]) )
    plot_nj_rinout_4.addBand(nj_rinout_4_mc.GetXaxis().GetXmin(), rinout_4_da[0]*0.8, nj_rinout_4_mc.GetXaxis().GetXmax(), rinout_4_da[0]*1.2, r.kOrange+6, 0.2)
    plot_nj_rinout_4.addLine(nj_rinout_4_mc.GetXaxis().GetXmin(), rinout_4_da[0], nj_rinout_4_mc.GetXaxis().GetXmax(), rinout_4_da[0] ,r.kBlue) 
    plot_nj_rinout_4.saveRatio(1, 1, 0, lumi, nj_rinout_4_da, nj_rinout_4_mc)                                                        
    plot_nj_rinout_5 = Canvas.Canvas('rinout/%s_%s/plot_njet_rinout_mll_300-400'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)           
    plot_nj_rinout_5.addHisto(nj_rinout_5_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_nj_rinout_5.addHisto(nj_rinout_5_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_nj_rinout_5.addBand(nj_rinout_5_mc.GetXaxis().GetXmin(), rinout_5_da[0]*0.8, nj_rinout_5_mc.GetXaxis().GetXmax(), rinout_5_da[0]*1.2, r.kOrange+6, 0.2)
    plot_nj_rinout_5.addLine(nj_rinout_5_mc.GetXaxis().GetXmin(), rinout_5_da[0], nj_rinout_5_mc.GetXaxis().GetXmax(), rinout_5_da[0] ,r.kBlue) 
    plot_nj_rinout_5.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 300-400"%(rinout_5_da[0]) )
    plot_nj_rinout_5.saveRatio(1, 1, 0, lumi, nj_rinout_5_da, nj_rinout_5_mc)                                                        
    plot_nj_rinout_6 = Canvas.Canvas('rinout/%s_%s/plot_njet_rinout_mll_400'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)       
    plot_nj_rinout_6.addHisto(nj_rinout_6_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_nj_rinout_6.addHisto(nj_rinout_6_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_nj_rinout_6.addBand(nj_rinout_6_mc.GetXaxis().GetXmin(), rinout_6_da[0]*0.8, nj_rinout_6_mc.GetXaxis().GetXmax(), rinout_6_da[0]*1.2, r.kOrange+6, 0.2)
    plot_nj_rinout_6.addLine(nj_rinout_6_mc.GetXaxis().GetXmin(), rinout_6_da[0], nj_rinout_6_mc.GetXaxis().GetXmax(), rinout_6_da[0] ,r.kBlue) 
    plot_nj_rinout_6.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} > 400"%(rinout_6_da[0]) )
    plot_nj_rinout_6.saveRatio(1, 1, 0, lumi, nj_rinout_6_da, nj_rinout_6_mc)                                                                                              
    

if __name__ == '__main__':

    print bcolors.HEADER 
    print '#######################################################################'
    print '                  Starting r_inout analysis...                          ' 
    print '#######################################################################' + bcolors.ENDC

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()

    print bcolors.HEADER + '[RinoutAnalysis] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    mcDatasets = ['TTJets_DiLepton', 'TTJets_DiLepton_ext',   'DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO','ZZTo4L', 'WZTo3LNu', 'WWW', 'WWZ','WZZ', 'ZZZ',  'TTZToLLNuNu' ,'TTWToLNu', 'T_tWch', 'TBar_tWch' , 'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT', 'TToLeptons_sch', 'TToLeptons_tch_powheg', 'TBarToLeptons_tch_powheg', 'TTHnobb_pow', 'VHToNonbb', 'WJetsToLNu_LO']
    daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2_runs_273150_275376', 'DoubleEG_Run2016B-PromptReco-v2_runs_273150_275376', 'MuonEG_Run2016B-PromptReco-v2_runs_273150_275376',
                  'DoubleMuon_Run2016C-PromptReco-v2_runs_275420_276283', 'DoubleEG_Run2016C-PromptReco-v2_runs_275420_276283', 'MuonEG_Run2016C-PromptReco-v2_runs_275420_276283',
                  'DoubleMuon_Run2016D-PromptReco-v2_runs_276315_276811', 'DoubleEG_Run2016D-PromptReco-v2_runs_276315_276811', 'MuonEG_Run2016D-PromptReco-v2_runs_276315_276811']

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)
  
    print bcolors.HEADER + '[RinoutAnalysis] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    maxrun = 999999
    lumi = 12.9 ; maxrun = 276811; lumi_str = '12.9invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 

    cuts = CutManager.CutManager()
    
    runAnalysis(lumi, treeDA, treeMC, cuts, '', 'nocut', True, opts.ingredientsFile)
 

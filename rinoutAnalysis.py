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

import include.LeptonSF
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

def makeResultsTable(mc, da, rinout_20_60_mc, rinout_60_86_mc,rinout_96_150_mc, rinout_150_200_mc, rinout_200_300_mc, rinout_300_400_mc, rinout_400_mc, rinout_20_60_da,  rinout_60_86_da, rinout_96_150_da, rinout_150_200_da, rinout_200_300_da, rinout_300_400_da, rinout_400_da,lumi, lumi_str):
    line00 = '&Data N_{in}:  %.0f $\\pm$ %.0f  & & MC N_{in}:  %.0f $\\pm$ %.0f \\\\'%(da.GetBinContent(3), da.GetBinError(3),mc.GetBinContent(3), mc.GetBinError(3))
    line0  = 'm_{ll} & N_{out} & R_{inout} & N_{out} & R_{inout} \\\\\hline'                   
    line1  = '20-60& %.0f $\\pm$ %.0f &  %.3f $\\pm$ %.3f & %.0f $\\pm$ %.0f & %.3f $\\pm$ %.3f  \\\\'%(da.GetBinContent(1), da.GetBinError(1), rinout_20_60_da[0], rinout_20_60_da[1], mc.GetBinContent(1), mc.GetBinError(1), rinout_20_60_mc[0], rinout_20_60_mc[1]) 
    line2  = '60-86& %.0f $\\pm$ %.0f &  %.3f $\\pm$ %.3f & %.0f $\\pm$ %.0f & %.3f $\\pm$ %.3f   \\\\'%(da.GetBinContent(2), da.GetBinError(2), rinout_60_86_da[0], rinout_60_86_da[1], mc.GetBinContent(2), mc.GetBinError(2), rinout_60_86_mc[0], rinout_60_86_mc[1]) 
    line3  = '96-150& %.0f $\\pm$ %.0f &  %.3f $\\pm$ %.3f & %.0f $\\pm$ %.0f & %.3f $\\pm$ %.3f   \\\\'%(da.GetBinContent(4), da.GetBinError(4), rinout_96_150_da[0], rinout_96_150_da[1], mc.GetBinContent(4), mc.GetBinError(4), rinout_96_150_mc[0], rinout_96_150_mc[1]) 
    line4  = '150-200& %.0f $\\pm$ %.0f &  %.3f $\\pm$ %.3f & %.0f $\\pm$ %.0f & %.3f $\\pm$ %.3f   \\\\'%(da.GetBinContent(5), da.GetBinError(5), rinout_150_200_da[0], rinout_150_200_da[1], mc.GetBinContent(5), mc.GetBinError(5), rinout_150_200_mc[0], rinout_150_200_mc[1]) 
    line5  = '200-300& %.0f $\\pm$ %.0f &  %.3f $\\pm$ %.3f & %.0f $\\pm$ %.0f & %.3f $\\pm$ %.3f   \\\\'%(da.GetBinContent(6), da.GetBinError(6), rinout_200_300_da[0], rinout_200_300_da[1], mc.GetBinContent(6), mc.GetBinError(6), rinout_200_300_mc[0], rinout_200_300_mc[1]) 
    line6  = '300-400& %.0f $\\pm$ %.0f &  %.3f $\\pm$ %.3f & %.0f $\\pm$ %.0f & %.3f $\\pm$ %.3f   \\\\'%(da.GetBinContent(7), da.GetBinError(7), rinout_300_400_da[0], rinout_300_400_da[1], mc.GetBinContent(7), mc.GetBinError(7), rinout_300_400_mc[0], rinout_300_400_mc[1]) 
    line7  = '+400& %.0f $\\pm$ %.0f &  %.3f $\\pm$ %.3f & %.0f $\\pm$ %.0f & %.3f $\\pm$ %.3f   \\\\'%(da.GetBinContent(8), da.GetBinError(8), rinout_400_da[0], rinout_400_da[1], mc.GetBinContent(8), mc.GetBinError(8), rinout_400_mc[0], rinout_400_mc[1]) 
    
    helper.ensurePath('plots/rinout/%s/tables/resultTable_%s%s_OLD.txt'%(lumi_str, str(lumi).replace('.','p'), "rinout"))
    compTableFile = open('plots/rinout/%s/tables/resultTable_%s%s_OLD.txt'%(lumi_str, str(lumi).replace('.','p'), "rinout"),'w')
    compTableFile.write(line00+'\n')
    compTableFile.write(line0+'\n')
    compTableFile.write(line1+'\n')
    compTableFile.write(line2+'\n')
    compTableFile.write(line3+'\n')
    compTableFile.write(line4+'\n')
    compTableFile.write(line5+'\n')
    compTableFile.write(line6+'\n')
    compTableFile.write(line7+'\n')
    compTableFile.close()
                                                                                                                                                                        
                                                                                                                                                                         
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


def calc_rinout(inSF,outSF, EinSF, EoutSF, syst):
    
    Ein  = math.sqrt(EinSF*EinSF +(syst*inSF)**2)
    Eout = math.sqrt(EoutSF*EoutSF +(syst*outSF)**2)
    print "Ein " ,Ein
    print "Eout " ,Eout
    return calc_frac(inSF, outSF, Ein, Eout)                                                                                         

def subtractTT(SF, OF, rsfof, dataMC):
    
      corr = SF.Clone("rinout_" + SF.GetName())
      corr.Add(OF , -rsfof)
      print '#################################################################'
      print ' ', dataMC
      print ' i have %.3f events in the SF sample' %(SF.Integral())
      print ' i have %.3f events in the OF sample' %(OF.Integral())
      print ' i have %.3f events in the subtracted sample' %(corr.Integral())
      print '#################################################################'
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

    makeNicePlot = True
    saveThisInFile = True
    labelx = "m_{ll} [GeV]"
    labelmet = "MET [GeV]"
    labelnjet = "N. Jets"
    inclusivebins  = range(20, 450, 5)
    metbins  = [0, 10, 20, 30, 40, 50, 60, 70, 90]
    njbins  = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.5]
    coarsebins  = [20, 60, 86, 96, 150, 200, 300, 400]

    #####Main mll control plot
    print bcolors.HEADER + '[RinoutAnalysis] ' + bcolors.OKBLUE + 'Starting mll plot' + bcolors.ENDC
    logx = 1
    if makeNicePlot:
       #make the nice finely binned plots for the presentation
        mll_SF_inclusive_mc = treeMC.getTH1F(lumi,"MC_inclusive_mll_SF","lepsMll_Edge",inclusivebins,1,1,cuts.AddList([cuts.goodLepton, cuts.SF, cuts.DYControlRegionNoMll]), '', labelx)
        mll_OF_inclusive_mc = treeMC.getTH1F(lumi,"MC_inclusive_mll_OF","lepsMll_Edge",inclusivebins,1,1,cuts.AddList([cuts.goodLepton, cuts.OF, cuts.DYControlRegionNoMll]), '', labelx)
        mll_SF_inclusive_da = treeDA.getTH1F(lumi,"DA_inclusive_mll_SF","lepsMll_Edge",inclusivebins,1,1,cuts.AddList([cuts.goodLepton, cuts.SF, cuts.DYControlRegionNoMll, cuts.trigger]),'', labelx)
        mll_OF_inclusive_da = treeDA.getTH1F(lumi,"DA_inclusive_mll_OF","lepsMll_Edge",inclusivebins,1,1,cuts.AddList([cuts.goodLepton, cuts.OF, cuts.DYControlRegionNoMll, cuts.trigger]),'', labelx)
    
        mll_corrected_inclusive_mc = subtractTT(mll_SF_inclusive_mc, mll_OF_inclusive_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
        mll_corrected_inclusive_da = subtractTT(mll_SF_inclusive_da, mll_OF_inclusive_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0], "DATA")
    
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
    mll_SF_coarse_mc = treeMC.getTH1F(lumi, "MC_inclusive_mll_SF", "lepsMll_Edge", coarsebins, 1, 1, cuts.AddList([cuts.mT2_80, cuts.SF, cuts.DYControlRegionNoMll]), '', labelx)
    mll_OF_coarse_mc = treeMC.getTH1F(lumi, "MC_inclusive_mll_OF", "lepsMll_Edge", coarsebins, 1, 1, cuts.AddList([cuts.mT2_80, cuts.OF, cuts.DYControlRegionNoMll]), '', labelx)
    mll_SF_coarse_da = treeDA.getTH1F(lumi, "DA_inclusive_mll_SF", "lepsMll_Edge", coarsebins, 1, 1, cuts.AddList([cuts.mT2_80, cuts.SF, cuts.DYControlRegionNoMll]), '', labelx)
    mll_OF_coarse_da = treeDA.getTH1F(lumi, "DA_inclusive_mll_OF", "lepsMll_Edge", coarsebins, 1, 1, cuts.AddList([cuts.mT2_80, cuts.OF, cuts.DYControlRegionNoMll]), '', labelx)
    mll_corrected_mc = subtractTT(mll_SF_coarse_mc, mll_OF_coarse_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    mll_corrected_da = subtractTT(mll_SF_coarse_da, mll_OF_coarse_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0], "DATA")
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

    rinout_20_60_mc  = calc_rinout(mll_corrected_mc.GetBinContent(3), mll_corrected_mc.GetBinContent(1) , mll_corrected_mc.GetBinError(3), mll_corrected_mc.GetBinError(1), 0.5)
    rinout_60_86_mc  = calc_rinout(mll_corrected_mc.GetBinContent(3), mll_corrected_mc.GetBinContent(2) , mll_corrected_mc.GetBinError(3), mll_corrected_mc.GetBinError(2), 0.5)
    rinout_86_96_mc  = calc_rinout(mll_corrected_mc.GetBinContent(3), mll_corrected_mc.GetBinContent(3) , mll_corrected_mc.GetBinError(3), mll_corrected_mc.GetBinError(3), 0.5)
    rinout_96_150_mc  = calc_rinout(mll_corrected_mc.GetBinContent(3), mll_corrected_mc.GetBinContent(4) , mll_corrected_mc.GetBinError(3), mll_corrected_mc.GetBinError(4), 0.5)
    rinout_150_200_mc  = calc_rinout(mll_corrected_mc.GetBinContent(3), mll_corrected_mc.GetBinContent(5) , mll_corrected_mc.GetBinError(3), mll_corrected_mc.GetBinError(5), 0.5)
    rinout_200_300_mc  = calc_rinout(mll_corrected_mc.GetBinContent(3), mll_corrected_mc.GetBinContent(6) , mll_corrected_mc.GetBinError(3), mll_corrected_mc.GetBinError(6), 1.)
    rinout_300_400_mc  = calc_rinout(mll_corrected_mc.GetBinContent(3), mll_corrected_mc.GetBinContent(7) , mll_corrected_mc.GetBinError(3), mll_corrected_mc.GetBinError(7), 1.)
    rinout_400_mc  = calc_rinout(mll_corrected_mc.GetBinContent(3), mll_corrected_mc.GetBinContent(8) , mll_corrected_mc.GetBinError(3), mll_corrected_mc.GetBinError(8), 1.)      
    rinout_20_60_da  = calc_rinout(mll_corrected_da.GetBinContent(3), mll_corrected_da.GetBinContent(1) , mll_corrected_da.GetBinError(3), mll_corrected_da.GetBinError(1), 0.5)
    rinout_60_86_da  = calc_rinout(mll_corrected_da.GetBinContent(3), mll_corrected_da.GetBinContent(2) , mll_corrected_da.GetBinError(3), mll_corrected_da.GetBinError(2), 0.5)
    rinout_86_96_da  = calc_rinout(mll_corrected_da.GetBinContent(3), mll_corrected_da.GetBinContent(3) , mll_corrected_da.GetBinError(3), mll_corrected_da.GetBinError(3), 0.5)
    rinout_96_150_da  = calc_rinout(mll_corrected_da.GetBinContent(3), mll_corrected_da.GetBinContent(4) , mll_corrected_da.GetBinError(3), mll_corrected_da.GetBinError(4), 0.5)
    rinout_150_200_da  = calc_rinout(mll_corrected_da.GetBinContent(3), mll_corrected_da.GetBinContent(5) , mll_corrected_da.GetBinError(3), mll_corrected_da.GetBinError(5), 0.5)
    rinout_200_300_da  = calc_rinout(mll_corrected_da.GetBinContent(3), mll_corrected_da.GetBinContent(6) , mll_corrected_da.GetBinError(3), mll_corrected_da.GetBinError(6),1.)
    rinout_300_400_da  = calc_rinout(mll_corrected_da.GetBinContent(3), mll_corrected_da.GetBinContent(7) , mll_corrected_da.GetBinError(3), mll_corrected_da.GetBinError(7), 1.)
    rinout_400_da  = calc_rinout(mll_corrected_da.GetBinContent(3), mll_corrected_da.GetBinContent(8) , mll_corrected_da.GetBinError(3), mll_corrected_da.GetBinError(8), 1.)      
    
    print 'rinout_20_60_da ', rinout_20_60_da[0] , rinout_20_60_da[1] 
    print 'rinout_60_86_da ', rinout_60_86_da[0] , rinout_60_86_da[1] 
    print 'rinout_86_96_da ', rinout_86_96_da[0] , rinout_86_96_da[1] 
    print 'rinout_96_150_da ', rinout_96_150_da[0] , rinout_96_150_da[1] 
    print 'rinout_150_200_da ', rinout_150_200_da[0] , rinout_150_200_da[1] 
    print 'rinout_200_300_da ', rinout_200_300_da[0] , rinout_200_300_da[1] 
    print 'rinout_300_400_da ', rinout_300_400_da[0] , rinout_300_400_da[1] 
    print 'rinout_400_da ', rinout_400_da[0] , rinout_400_da[1] 
    if saveThisInFile:
        makeResultsTable(mll_corrected_mc, mll_corrected_da, rinout_20_60_mc, rinout_60_86_mc, rinout_96_150_mc, rinout_150_200_mc, rinout_200_300_mc, rinout_300_400_mc, rinout_400_mc, rinout_20_60_da, rinout_60_86_da, rinout_96_150_da, rinout_150_200_da, rinout_200_300_da, rinout_300_400_da, rinout_400_da,  lumi, str(lumi))
        print "ingredientsFile", ingredientsFile
        saveInFile(ingredientsFile, 'dy_m20_60__', 'DATA', rinout_20_60_da[0], rinout_20_60_da[1],rinout_20_60_da[0]* 0.5)
        saveInFile(ingredientsFile, 'dy_m60_86__', 'DATA', rinout_60_86_da[0], rinout_60_86_da[1], rinout_60_86_da[0]*0.5)
        saveInFile(ingredientsFile, 'dy_m86_96__', 'DATA', rinout_86_96_da[0], rinout_86_96_da[1], rinout_86_96_da[0]*0.5)
        saveInFile(ingredientsFile, 'dy_m96_150_', 'DATA', rinout_96_150_da[0], rinout_96_150_da[1], rinout_96_150_da[0]*0.5)
        saveInFile(ingredientsFile, 'dy_m150_200', 'DATA', rinout_150_200_da[0], rinout_150_200_da[1], rinout_150_200_da[0]*1.)
        saveInFile(ingredientsFile, 'dy_m200_300', 'DATA', rinout_200_300_da[0], rinout_200_300_da[1], rinout_200_300_da[0]*1.)
        saveInFile(ingredientsFile, 'dy_m300_400', 'DATA', rinout_300_400_da[0], rinout_300_400_da[1], rinout_300_400_da[0]*1.)
        saveInFile(ingredientsFile, 'dy_m400____', 'DATA', rinout_400_da[0]      , rinout_400_da[1], rinout_400_da[0]*1.)
        saveInFile(ingredientsFile, 'dy_m20_60__', 'MC  ', rinout_20_60_mc[0], rinout_20_60_mc[1], rinout_20_60_mc[0]*0.5)
        saveInFile(ingredientsFile, 'dy_m60_86__', 'MC  ', rinout_60_86_mc[0], rinout_60_86_mc[1], rinout_60_86_mc[0]*0.5)
        saveInFile(ingredientsFile, 'dy_m86_96__', 'MC  ', rinout_86_96_mc[0], rinout_86_96_mc[1], rinout_86_96_mc[0]*0.5)
        saveInFile(ingredientsFile, 'dy_m96_150_', 'MC  ', rinout_96_150_mc[0], rinout_96_150_mc[1], rinout_96_150_mc[0]*0.5)
        saveInFile(ingredientsFile, 'dy_m150_200', 'MC  ', rinout_150_200_mc[0], rinout_150_200_mc[1], rinout_150_200_mc[0]*1.)
        saveInFile(ingredientsFile, 'dy_m200_300', 'MC  ', rinout_200_300_mc[0], rinout_200_300_mc[1], rinout_200_300_mc[0]*1.)
        saveInFile(ingredientsFile, 'dy_m300_400', 'MC  ', rinout_300_400_mc[0], rinout_300_400_mc[1], rinout_300_400_mc[0]*1.)
        saveInFile(ingredientsFile, 'dy_m400____', 'MC  ', rinout_400_mc[0], rinout_400_mc[1], rinout_400_mc[0]*1.)
        print "done writing in the ingredients.dat"

    met_1_SF_mc= treeMC.getTH1F(lumi,"MC_20-60SF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll20_60]),'', labelmet)
    met_1_SF_da= treeDA.getTH1F(lumi,"DA_20-60SF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll20_60, cuts.trigger]),'', labelmet)
    met_1_OF_mc= treeMC.getTH1F(lumi,"MC_20-60OF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll20_60]),'', labelmet)
    met_1_OF_da= treeDA.getTH1F(lumi,"DA_20-60OF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll20_60, cuts.trigger]),'', labelmet)
    met_2_SF_mc= treeMC.getTH1F(lumi,"MC_60-86SF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll60_86]),'', labelmet)
    met_2_SF_da= treeDA.getTH1F(lumi,"DA_60-86SF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll60_86, cuts.trigger]),'', labelmet)
    met_2_OF_mc= treeMC.getTH1F(lumi,"MC_60-86OF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll60_86]),'', labelmet)
    met_2_OF_da= treeDA.getTH1F(lumi,"DA_60-86OF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll60_86, cuts.trigger]),'', labelmet)
    met_z_SF_mc= treeMC.getTH1F(lumi,"MC_86-96SF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll86_96]),'', labelmet)
    met_z_SF_da= treeDA.getTH1F(lumi,"DA_86-96SF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll86_96, cuts.trigger]),'', labelmet)
    met_z_OF_mc= treeMC.getTH1F(lumi,"MC_86-96OF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll86_96]),'', labelmet)
    met_z_OF_da= treeDA.getTH1F(lumi,"DA_86-96OF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll86_96, cuts.trigger]),'', labelmet)
    met_3_SF_mc= treeMC.getTH1F(lumi,"MC_96-150SF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll96_150]),'', labelmet)
    met_3_SF_da= treeDA.getTH1F(lumi,"DA_96-150SF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll96_150, cuts.trigger]),'', labelmet)
    met_3_OF_mc= treeMC.getTH1F(lumi,"MC_96-150OF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll96_150]),'', labelmet)
    met_3_OF_da= treeDA.getTH1F(lumi,"DA_96-150OF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll96_150, cuts.trigger]),'', labelmet)
    met_4_SF_mc= treeMC.getTH1F(lumi,"MC_150-200SF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll150_200]),'', labelmet)
    met_4_SF_da= treeDA.getTH1F(lumi,"DA_150-200SF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll150_200, cuts.trigger]),'',labelmet)
    met_4_OF_mc= treeMC.getTH1F(lumi,"MC_150-200OF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll150_200]),'', labelmet)
    met_4_OF_da= treeDA.getTH1F(lumi,"DA_150-200OF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll150_200, cuts.trigger]),'',labelmet)
    met_5_SF_mc= treeMC.getTH1F(lumi,"MC_200-300SF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll200_300]),'', labelmet)
    met_5_SF_da= treeDA.getTH1F(lumi,"DA_200-300SF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll200_300, cuts.trigger]),'',labelmet)
    met_5_OF_mc= treeMC.getTH1F(lumi,"MC_200-300OF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll200_300]),'', labelmet)
    met_5_OF_da= treeDA.getTH1F(lumi,"DA_200-300OF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll200_300, cuts.trigger]),'',labelmet)
    met_6_SF_mc= treeMC.getTH1F(lumi,"MC_300-400SF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll300_400]),'', labelmet)
    met_6_SF_da= treeDA.getTH1F(lumi,"DA_300-400SF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll300_400, cuts.trigger]),'',labelmet)
    met_6_OF_mc= treeMC.getTH1F(lumi,"MC_300-400OF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll300_400]),'', labelmet)
    met_6_OF_da= treeDA.getTH1F(lumi,"DA_300-400OF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll300_400, cuts.trigger]),'',labelmet)
    met_7_SF_mc= treeMC.getTH1F(lumi,"MC_400_SF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll400]),'', labelmet)
    met_7_SF_da= treeDA.getTH1F(lumi,"DA_400_SF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll400, cuts.trigger]),'', labelmet)
    met_7_OF_mc= treeMC.getTH1F(lumi,"MC_400_OF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll400]),'', labelmet)
    met_7_OF_da= treeDA.getTH1F(lumi,"DA_400_OF","met_Edge",metbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.METl100, cuts.DYControlRegionNoMllNoMET,cuts.mll400, cuts.trigger]),'', labelmet)
    #subtract the OF from the total
    met_1_corrected_mc  = subtractTT(met_1_SF_mc, met_1_OF_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    met_1_corrected_da  = subtractTT(met_1_SF_da, met_1_OF_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0],  "DATA")
    met_2_corrected_mc  = subtractTT(met_2_SF_mc, met_2_OF_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    met_2_corrected_da  = subtractTT(met_2_SF_da, met_2_OF_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0],  "DATA")
    met_z_corrected_mc  = subtractTT(met_z_SF_mc, met_z_OF_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    met_z_corrected_da  = subtractTT(met_z_SF_da, met_z_OF_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0],  "DATA")
    met_3_corrected_mc  = subtractTT(met_3_SF_mc, met_3_OF_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    met_3_corrected_da  = subtractTT(met_3_SF_da, met_3_OF_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0],  "DATA")
    met_4_corrected_mc  = subtractTT(met_4_SF_mc, met_4_OF_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    met_4_corrected_da  = subtractTT(met_4_SF_da, met_4_OF_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0],  "DATA")
    met_5_corrected_mc  = subtractTT(met_5_SF_mc, met_5_OF_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    met_5_corrected_da  = subtractTT(met_5_SF_da, met_5_OF_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0],  "DATA")
    met_6_corrected_mc  = subtractTT(met_6_SF_mc, met_6_OF_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    met_6_corrected_da  = subtractTT(met_6_SF_da, met_6_OF_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0],  "DATA")
    met_7_corrected_mc  = subtractTT(met_7_SF_mc, met_7_OF_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    met_7_corrected_da  = subtractTT(met_7_SF_da, met_7_OF_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0],  "DATA")
   
    #calculate rinout with the corrected (OF subtracted) samples
    met_rinout_20_60_mc = make_rinout_histo(met_1_corrected_mc, met_z_corrected_mc, 0., 0.5 )
    met_rinout_20_60_da = make_rinout_histo(met_1_corrected_da, met_z_corrected_da, 0., 0.5 )  
    met_rinout_60_86_mc = make_rinout_histo(met_2_corrected_mc, met_z_corrected_mc, 0., 1. )
    met_rinout_60_86_da = make_rinout_histo(met_2_corrected_da, met_z_corrected_da, 0., 1. )  
    met_rinout_96_150_mc = make_rinout_histo(met_3_corrected_mc, met_z_corrected_mc, 0., 0.2) 
    met_rinout_96_150_da = make_rinout_histo(met_3_corrected_da, met_z_corrected_da, 0., 0.2)
    met_rinout_150_200_mc = make_rinout_histo(met_4_corrected_mc, met_z_corrected_mc, 0., 0.03)  
    met_rinout_150_200_da = make_rinout_histo(met_4_corrected_da, met_z_corrected_da, 0., 0.03)
    met_rinout_200_300_mc = make_rinout_histo(met_5_corrected_mc, met_z_corrected_mc, 0., 0.03) 
    met_rinout_200_300_da = make_rinout_histo(met_5_corrected_da, met_z_corrected_da, 0., 0.03)
    met_rinout_300_400_mc = make_rinout_histo(met_6_corrected_mc, met_z_corrected_mc, 0., 0.03) 
    met_rinout_300_400_da = make_rinout_histo(met_6_corrected_da, met_z_corrected_da, 0., 0.03)
    met_rinout_400_mc     = make_rinout_histo(met_7_corrected_mc, met_z_corrected_mc, 0., 0.03) 
    met_rinout_400_da     = make_rinout_histo(met_7_corrected_da, met_z_corrected_da, 0., 0.03) 

    plot_met_rinout_20_60 = Canvas.Canvas('rinout/%s_%s/plot_met_rinout_mll_20-60'%(lumi_str,tag), 'png,pdf', 0.6, 0.7, 0.85, 0.9)
    plot_met_rinout_20_60.addHisto(met_rinout_20_60_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_met_rinout_20_60.addHisto(met_rinout_20_60_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_met_rinout_20_60.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 20-60"%(rinout_20_60_da[0]) )
    plot_met_rinout_20_60.addBand(met_rinout_20_60_mc.GetXaxis().GetXmin(), rinout_20_60_da[0]*0.8, met_rinout_20_60_mc.GetXaxis().GetXmax(), rinout_20_60_da[0]*1.2, r.kOrange+6, 0.2)
    plot_met_rinout_20_60.addLine(met_rinout_20_60_mc.GetXaxis().GetXmin(), rinout_20_60_da[0], met_rinout_20_60_mc.GetXaxis().GetXmax(), rinout_20_60_da[0] ,r.kBlue)
    plot_met_rinout_20_60.saveRatio(1, 1, 0, lumi, met_rinout_20_60_da, met_rinout_20_60_mc)                                                                                                  
    
    plot_met_rinout_60_86 = Canvas.Canvas('rinout/%s_%s/plot_met_rinout_mll_60-86'%(lumi_str,tag), 'png,pdf', 0.6, 0.7, 0.85, 0.9)
    plot_met_rinout_60_86.addHisto(met_rinout_60_86_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)                                                                                         
    plot_met_rinout_60_86.addHisto(met_rinout_60_86_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)                                                                                   
    plot_met_rinout_60_86.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 60-86"%(rinout_20_60_da[0]) )                                                                             
    plot_met_rinout_60_86.addBand(met_rinout_60_86_mc.GetXaxis().GetXmin(), rinout_60_86_da[0]*0.8, met_rinout_60_86_mc.GetXaxis().GetXmax(), rinout_60_86_da[0]*1.2, r.kOrange+6, 0.2)         
    plot_met_rinout_60_86.addLine(met_rinout_60_86_mc.GetXaxis().GetXmin(), rinout_60_86_da[0], met_rinout_60_86_mc.GetXaxis().GetXmax(), rinout_60_86_da[0] ,r.kBlue)                          
    plot_met_rinout_60_86.saveRatio(1, 1, 0, lumi, met_rinout_60_86_da, met_rinout_60_86_mc)                                                                                            
    
    plot_met_rinout_96_150 = Canvas.Canvas('rinout/%s_%s/plot_met_rinout_mll_96-150'%(lumi_str,tag), 'png,pdf', 0.6, 0.7, 0.85, 0.9)           
    plot_met_rinout_96_150.addHisto(met_rinout_96_150_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_met_rinout_96_150.addHisto(met_rinout_96_150_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_met_rinout_96_150.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 96-150"%(rinout_96_150_da[0]) )
    plot_met_rinout_96_150.addBand(met_rinout_96_150_mc.GetXaxis().GetXmin(), rinout_96_150_da[0]*0.8, met_rinout_96_150_mc.GetXaxis().GetXmax(), rinout_96_150_da[0]*1.2, r.kOrange+6, 0.2)
    plot_met_rinout_96_150.addLine(met_rinout_96_150_mc.GetXaxis().GetXmin(), rinout_96_150_da[0], met_rinout_96_150_mc.GetXaxis().GetXmax(), rinout_96_150_da[0] ,r.kBlue)
    plot_met_rinout_96_150.saveRatio(1, 1, 0, lumi, met_rinout_96_150_da, met_rinout_96_150_mc)                                                        
    
    plot_met_rinout_150_200 = Canvas.Canvas('rinout/%s_%s/plot_met_rinout_mll_150-200'%(lumi_str,tag), 'png,pdf', 0.6, 0.7, 0.85, 0.9)           
    plot_met_rinout_150_200.addHisto(met_rinout_150_200_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_met_rinout_150_200.addHisto(met_rinout_150_200_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_met_rinout_150_200.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 150-200"%(rinout_150_200_da[0]) )
    plot_met_rinout_150_200.addBand(met_rinout_150_200_mc.GetXaxis().GetXmin(), rinout_150_200_da[0]*0.8, met_rinout_150_200_mc.GetXaxis().GetXmax(), rinout_150_200_da[0]*1.2, r.kOrange+6, 0.2)
    plot_met_rinout_150_200.addLine(met_rinout_150_200_mc.GetXaxis().GetXmin(), rinout_150_200_da[0], met_rinout_150_200_mc.GetXaxis().GetXmax(), rinout_150_200_da[0] ,r.kBlue)
    plot_met_rinout_150_200.saveRatio(1, 1, 0, lumi, met_rinout_150_200_da, met_rinout_150_200_mc)                                                        
    
    plot_met_rinout_200_300 = Canvas.Canvas('rinout/%s_%s/plot_met_rinout_mll_200-300'%(lumi_str,tag), 'png,pdf', 0.6, 0.7, 0.85, 0.9)           
    plot_met_rinout_200_300.addHisto(met_rinout_200_300_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_met_rinout_200_300.addHisto(met_rinout_200_300_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_met_rinout_200_300.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 200-300"%(rinout_200_300_da[0]) )
    plot_met_rinout_200_300.addBand(met_rinout_200_300_mc.GetXaxis().GetXmin(), rinout_200_300_da[0]*0.8, met_rinout_200_300_mc.GetXaxis().GetXmax(), rinout_200_300_da[0]*1.2, r.kOrange+6, 0.2)
    plot_met_rinout_200_300.addLine(met_rinout_200_300_mc.GetXaxis().GetXmin(), rinout_200_300_da[0], met_rinout_200_300_mc.GetXaxis().GetXmax(), rinout_200_300_da[0] ,r.kBlue)
    plot_met_rinout_200_300.saveRatio(1, 1, 0, lumi, met_rinout_200_300_da, met_rinout_200_300_mc)                                                        
    
    plot_met_rinout_300_400 = Canvas.Canvas('rinout/%s_%s/plot_met_rinout_mll_300-400'%(lumi_str,tag), 'png,pdf', 0.6, 0.7, 0.85, 0.9)       
    plot_met_rinout_300_400.addHisto(met_rinout_300_400_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_met_rinout_300_400.addHisto(met_rinout_300_400_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_met_rinout_300_400.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 300- 400"%(rinout_300_400_da[0]) )
    plot_met_rinout_300_400.addBand(met_rinout_300_400_mc.GetXaxis().GetXmin(), rinout_300_400_da[0]*0.8, met_rinout_300_400_mc.GetXaxis().GetXmax(), rinout_300_400_da[0]*1.2, r.kOrange+6, 0.2)
    plot_met_rinout_300_400.addLine(met_rinout_300_400_mc.GetXaxis().GetXmin(), rinout_300_400_da[0], met_rinout_300_400_mc.GetXaxis().GetXmax(), rinout_300_400_da[0] ,r.kBlue)
    plot_met_rinout_300_400.saveRatio(1, 1, 0, lumi, met_rinout_300_400_da, met_rinout_300_400_mc)                                                                                                  
    plot_met_rinout_400 = Canvas.Canvas('rinout/%s_%s/plot_met_rinout_mll_400'%(lumi_str,tag), 'png,pdf', 0.6, 0.7, 0.85, 0.9)       
    plot_met_rinout_400.addHisto(met_rinout_400_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_met_rinout_400.addHisto(met_rinout_400_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_met_rinout_400.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} + 400"%(rinout_300_400_da[0]) )
    plot_met_rinout_400.addBand(met_rinout_400_mc.GetXaxis().GetXmin(), rinout_400_da[0]*0.8, met_rinout_400_mc.GetXaxis().GetXmax(), rinout_400_da[0]*1.2, r.kOrange+6, 0.2)
    plot_met_rinout_400.addLine(met_rinout_400_mc.GetXaxis().GetXmin(), rinout_400_da[0], met_rinout_400_mc.GetXaxis().GetXmax(), rinout_400_da[0] ,r.kBlue)
    plot_met_rinout_400.saveRatio(1, 1, 0, lumi, met_rinout_400_da, met_rinout_400_mc)                                                                                                 

    #make njet plot##############################################################################################################################################################################
    nj_1_SF_mc=treeMC.getTH1F(lumi,"M_n20-60_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll20_60]),'', labelnjet) 
    nj_1_SF_da=treeDA.getTH1F(lumi,"D_n20-60_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll20_60]),'', labelnjet)
    nj_1_OF_mc=treeMC.getTH1F(lumi,"M_n20-60_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll20_60]),'', labelnjet)
    nj_1_OF_da=treeDA.getTH1F(lumi,"D_n20-60_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll20_60]),'', labelnjet)
    nj_2_SF_mc=treeMC.getTH1F(lumi,"M_n60-86_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll60_86]),'', labelnjet)
    nj_2_SF_da=treeDA.getTH1F(lumi,"D_n60-86_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll60_86]),'', labelnjet)
    nj_2_OF_mc=treeMC.getTH1F(lumi,"M_n60-86_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll60_86]),'', labelnjet)
    nj_2_OF_da=treeDA.getTH1F(lumi,"D_n60-86_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll60_86]),'', labelnjet)
    nj_z_SF_mc=treeMC.getTH1F(lumi,"M_n86-96_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll86_96]),'', labelnjet)
    nj_z_SF_da=treeDA.getTH1F(lumi,"D_n86-96_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll86_96]),'', labelnjet)
    nj_z_OF_mc=treeMC.getTH1F(lumi,"M_n86-96_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll86_96]),'', labelnjet)
    nj_z_OF_da=treeDA.getTH1F(lumi,"D_n86-96_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll86_96]),'', labelnjet)
    nj_3_SF_mc=treeMC.getTH1F(lumi,"M_n96-150_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll96_150]),'',labelnjet)
    nj_3_SF_da=treeDA.getTH1F(lumi,"D_n96-150_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll96_150]),'',labelnjet)
    nj_3_OF_mc=treeMC.getTH1F(lumi,"M_n96-150_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll96_150]),'',labelnjet)
    nj_3_OF_da=treeDA.getTH1F(lumi,"D_n96-150_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll96_150]),'',labelnjet)
    nj_4_SF_mc=treeMC.getTH1F(lumi,"M_n150-200_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll150_200]),'',labelnjet)
    nj_4_SF_da=treeDA.getTH1F(lumi,"D_n150-200_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll150_200]),'',labelnjet)
    nj_4_OF_mc=treeMC.getTH1F(lumi,"M_n150-200_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll150_200]),'',labelnjet)
    nj_4_OF_da=treeDA.getTH1F(lumi,"D_n150-200_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll150_200]),'',labelnjet)
    nj_5_SF_mc=treeMC.getTH1F(lumi,"M_n200-300_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll200_300]),'',labelnjet)
    nj_5_SF_da=treeDA.getTH1F(lumi,"D_n200-300_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll200_300]),'',labelnjet)
    nj_5_OF_mc=treeMC.getTH1F(lumi,"M_n200-300_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll200_300]),'',labelnjet)
    nj_5_OF_da=treeDA.getTH1F(lumi,"D_n200-300_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll200_300]),'',labelnjet)
    nj_6_SF_mc=treeMC.getTH1F(lumi,"M_n300-400_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll300_400]),'',labelnjet)
    nj_6_SF_da=treeDA.getTH1F(lumi,"D_n300-400_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF,cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll300_400]),'',labelnjet)
    nj_6_OF_mc=treeMC.getTH1F(lumi,"M_n300-400_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll300_400]),'',labelnjet)
    nj_6_OF_da=treeDA.getTH1F(lumi,"D_n300-400_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF,cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll300_400]),'',labelnjet)
    nj_7_SF_mc=treeMC.getTH1F(lumi,"M_n400_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll400]),'',labelnjet)
    nj_7_SF_da=treeDA.getTH1F(lumi,"D_n400_SF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.SF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll400]),'',labelnjet)
    nj_7_OF_mc=treeMC.getTH1F(lumi,"M_n400_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll400]),'',labelnjet)
    nj_7_OF_da=treeDA.getTH1F(lumi,"D_n400_OF","nJetSel_Edge",njbins,1,1,cuts.AddList([cuts.mT2_80,cuts.OF, cuts.DYControlRegionNoMllNoNJet,cuts.METl100,cuts.mll400]),'',labelnjet)
                                                                                                                             
    nj_1_corrected_mc  = subtractTT(nj_1_SF_mc, nj_1_OF_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    nj_1_corrected_da  = subtractTT(nj_1_SF_da, nj_1_OF_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0],  "DATA")
    nj_2_corrected_mc  = subtractTT(nj_2_SF_mc, nj_2_OF_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    nj_2_corrected_da  = subtractTT(nj_2_SF_da, nj_2_OF_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0],  "DATA")
    nj_z_corrected_mc  = subtractTT(nj_z_SF_mc, nj_z_OF_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    nj_z_corrected_da  = subtractTT(nj_z_SF_da, nj_z_OF_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0],  "DATA")
    nj_3_corrected_mc  = subtractTT(nj_3_SF_mc, nj_3_OF_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    nj_3_corrected_da  = subtractTT(nj_3_SF_da, nj_3_OF_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0],  "DATA")
    nj_4_corrected_mc  = subtractTT(nj_4_SF_mc, nj_4_OF_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    nj_4_corrected_da  = subtractTT(nj_4_SF_da, nj_4_OF_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0],  "DATA")
    nj_5_corrected_mc  = subtractTT(nj_5_SF_mc, nj_5_OF_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    nj_5_corrected_da  = subtractTT(nj_5_SF_da, nj_5_OF_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0],  "DATA")
    nj_6_corrected_mc  = subtractTT(nj_6_SF_mc, nj_6_OF_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    nj_6_corrected_da  = subtractTT(nj_6_SF_da, nj_6_OF_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0],  "DATA")
    nj_7_corrected_mc  = subtractTT(nj_7_SF_mc, nj_7_OF_mc, helper.readFromFileRsfofD(ingredientsFile, "MC")[0],  "MC")
    nj_7_corrected_da  = subtractTT(nj_7_SF_da, nj_7_OF_da, helper.readFromFileRsfofD(ingredientsFile, "DATA")[0],  "DATA")
    
    nj_rinout_20_60_mc = make_rinout_histo(nj_1_corrected_mc, nj_z_corrected_mc, 0., 0.5 )
    nj_rinout_20_60_da = make_rinout_histo(nj_1_corrected_da, nj_z_corrected_da, 0., 0.5 )  
    nj_rinout_60_86_mc = make_rinout_histo(nj_2_corrected_mc, nj_z_corrected_mc, 0., 1. )
    nj_rinout_60_86_da = make_rinout_histo(nj_2_corrected_da, nj_z_corrected_da, 0., 1. )  
    nj_rinout_96_150_mc = make_rinout_histo(nj_3_corrected_mc, nj_z_corrected_mc, 0., 0.2 ) 
    nj_rinout_96_150_da = make_rinout_histo(nj_3_corrected_da, nj_z_corrected_da, 0., 0.2 )
    nj_rinout_150_200_mc = make_rinout_histo(nj_4_corrected_mc, nj_z_corrected_mc, 0., 0.05 )  
    nj_rinout_150_200_da = make_rinout_histo(nj_4_corrected_da, nj_z_corrected_da, 0., 0.05 )
    nj_rinout_200_300_mc = make_rinout_histo(nj_5_corrected_mc, nj_z_corrected_mc, 0., 0.02 ) 
    nj_rinout_200_300_da = make_rinout_histo(nj_5_corrected_da, nj_z_corrected_da, 0., 0.02 )
    nj_rinout_300_400_mc = make_rinout_histo(nj_6_corrected_mc, nj_z_corrected_mc, 0., 0.02 ) 
    nj_rinout_300_400_da = make_rinout_histo(nj_6_corrected_da, nj_z_corrected_da, 0., 0.02 )
    nj_rinout_400_mc = make_rinout_histo(nj_7_corrected_mc, nj_z_corrected_mc, 0., 0.02 ) 
    nj_rinout_400_da = make_rinout_histo(nj_7_corrected_da, nj_z_corrected_da, 0., 0.02 )
                                                                                                                                                                                       
    plot_nj_rinout_20_60 = Canvas.Canvas('rinout/%s_%s/plot_njet_rinout_mll_20-60'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)
    plot_nj_rinout_20_60.addHisto(nj_rinout_20_60_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_nj_rinout_20_60.addHisto(nj_rinout_20_60_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_nj_rinout_20_60.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 20-60"%(rinout_20_60_da[0]) )
    plot_nj_rinout_20_60.addBand(nj_rinout_20_60_mc.GetXaxis().GetXmin(), rinout_20_60_da[0]*0.8, nj_rinout_20_60_mc.GetXaxis().GetXmax(), rinout_20_60_da[0]*1.2, r.kOrange+6, 0.2)
    plot_nj_rinout_20_60.addLine(nj_rinout_20_60_mc.GetXaxis().GetXmin(), rinout_20_60_da[0], nj_rinout_20_60_mc.GetXaxis().GetXmax(), rinout_20_60_da[0] ,r.kBlue)
    plot_nj_rinout_20_60.saveRatio(1, 1, 0, lumi, nj_rinout_20_60_da, nj_rinout_20_60_mc)                                                                                                  
    plot_nj_rinout_60_86 = Canvas.Canvas('rinout/%s_%s/plot_njet_rinout_mll_60-86'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)
    plot_nj_rinout_60_86.addHisto(nj_rinout_60_86_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_nj_rinout_60_86.addHisto(nj_rinout_60_86_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_nj_rinout_60_86.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 60-86"%(rinout_60_86_da[0]) )
    plot_nj_rinout_60_86.addBand(nj_rinout_60_86_mc.GetXaxis().GetXmin(), rinout_60_86_da[0]*0.8, nj_rinout_60_86_mc.GetXaxis().GetXmax(), rinout_60_86_da[0]*1.2, r.kOrange+6, 0.2)
    plot_nj_rinout_60_86.addLine(nj_rinout_60_86_mc.GetXaxis().GetXmin(), rinout_60_86_da[0], nj_rinout_60_86_mc.GetXaxis().GetXmax(), rinout_60_86_da[0] ,r.kBlue)
    plot_nj_rinout_60_86.saveRatio(1, 1, 0, lumi, nj_rinout_60_86_da, nj_rinout_60_86_mc)                                                                                         
    plot_nj_rinout_96_150 = Canvas.Canvas('rinout/%s_%s/plot_njet_rinout_mll_96-150'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)           
    plot_nj_rinout_96_150.addHisto(nj_rinout_96_150_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_nj_rinout_96_150.addHisto(nj_rinout_96_150_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_nj_rinout_96_150.addBand(nj_rinout_96_150_mc.GetXaxis().GetXmin(), rinout_96_150_da[0]*0.8, nj_rinout_96_150_mc.GetXaxis().GetXmax(), rinout_96_150_da[0]*1.2, r.kOrange+6, 0.2)
    plot_nj_rinout_96_150.addLine(nj_rinout_96_150_mc.GetXaxis().GetXmin(), rinout_96_150_da[0], nj_rinout_96_150_mc.GetXaxis().GetXmax(), rinout_96_150_da[0] ,r.kBlue) 
    plot_nj_rinout_96_150.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 96-150"%(rinout_96_150_da[0]) )
    plot_nj_rinout_96_150.saveRatio(1, 1, 0, lumi, nj_rinout_96_150_da, nj_rinout_96_150_mc)                                                        
    plot_nj_rinout_150_200 = Canvas.Canvas('rinout/%s_%s/plot_njet_rinout_mll_150-200'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)           
    plot_nj_rinout_150_200.addHisto(nj_rinout_150_200_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_nj_rinout_150_200.addHisto(nj_rinout_150_200_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_nj_rinout_150_200.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 150-200"%(rinout_150_200_da[0]) )
    plot_nj_rinout_150_200.addBand(nj_rinout_150_200_mc.GetXaxis().GetXmin(), rinout_150_200_da[0]*0.8, nj_rinout_150_200_mc.GetXaxis().GetXmax(), rinout_150_200_da[0]*1.2, r.kOrange+6, 0.2)
    plot_nj_rinout_150_200.addLine(nj_rinout_150_200_mc.GetXaxis().GetXmin(), rinout_150_200_da[0], nj_rinout_150_200_mc.GetXaxis().GetXmax(), rinout_150_200_da[0] ,r.kBlue) 
    plot_nj_rinout_150_200.saveRatio(1, 1, 0, lumi, nj_rinout_150_200_da, nj_rinout_150_200_mc)                                                        
    plot_nj_rinout_200_300 = Canvas.Canvas('rinout/%s_%s/plot_njet_rinout_mll_200-300'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)           
    plot_nj_rinout_200_300.addHisto(nj_rinout_200_300_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_nj_rinout_200_300.addHisto(nj_rinout_200_300_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_nj_rinout_200_300.addBand(nj_rinout_200_300_mc.GetXaxis().GetXmin(), rinout_200_300_da[0]*0.8, nj_rinout_200_300_mc.GetXaxis().GetXmax(), rinout_200_300_da[0]*1.2, r.kOrange+6, 0.2)
    plot_nj_rinout_200_300.addLine(nj_rinout_200_300_mc.GetXaxis().GetXmin(), rinout_200_300_da[0], nj_rinout_200_300_mc.GetXaxis().GetXmax(), rinout_200_300_da[0] ,r.kBlue) 
    plot_nj_rinout_200_300.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 200-300"%(rinout_200_300_da[0]) )
    plot_nj_rinout_200_300.saveRatio(1, 1, 0, lumi, nj_rinout_200_300_da, nj_rinout_200_300_mc)                                                        
    plot_nj_rinout_300_400 = Canvas.Canvas('rinout/%s_%s/plot_njet_rinout_mll_300-400'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)       
    plot_nj_rinout_300_400.addHisto(nj_rinout_300_400_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_nj_rinout_300_400.addHisto(nj_rinout_300_400_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_nj_rinout_300_400.addBand(nj_rinout_300_400_mc.GetXaxis().GetXmin(), rinout_300_400_da[0]*0.8, nj_rinout_300_400_mc.GetXaxis().GetXmax(), rinout_300_400_da[0]*1.2, r.kOrange+6, 0.2)
    plot_nj_rinout_300_400.addLine(nj_rinout_300_400_mc.GetXaxis().GetXmin(), rinout_300_400_da[0], nj_rinout_300_400_mc.GetXaxis().GetXmax(), rinout_300_400_da[0] ,r.kBlue) 
    plot_nj_rinout_300_400.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 300-400"%(rinout_300_400_da[0]) )
    plot_nj_rinout_300_400.saveRatio(1, 1, 0, lumi, nj_rinout_300_400_da, nj_rinout_300_400_mc)                                                                                                  
    plot_nj_rinout_400 = Canvas.Canvas('rinout/%s_%s/plot_njet_rinout_mll_400'%(lumi_str,tag), 'png,pdf', 0.7, 0.7, 0.85, 0.9)       
    plot_nj_rinout_400.addHisto(nj_rinout_400_mc, 'E1', 'MC', 'L', r.kGreen+1 , 1, 0)
    plot_nj_rinout_400.addHisto(nj_rinout_400_da, 'E1,SAME', 'DATA', 'PL', r.kBlack , 1, 0)
    plot_nj_rinout_400.addBand(nj_rinout_400_mc.GetXaxis().GetXmin(), rinout_400_da[0]*0.8, nj_rinout_400_mc.GetXaxis().GetXmax(), rinout_400_da[0]*1.2, r.kOrange+6, 0.2)
    plot_nj_rinout_400.addLine(nj_rinout_400_mc.GetXaxis().GetXmin(), rinout_400_da[0], nj_rinout_400_mc.GetXaxis().GetXmax(), rinout_400_da[0] ,r.kBlue) 
    plot_nj_rinout_400.addLatex (0.6, 0.6, "<r_{inout}> %.3f in m_{ll} 300-400"%(rinout_400_da[0]) )
    plot_nj_rinout_400.saveRatio(1, 1, 0, lumi, nj_rinout_400_da, nj_rinout_400_mc)                                                                                                  

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

    DYDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    ttDatasets = ['TTJets_DiLepton']
    mcDatasets = ['TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT',  'T_tch_powheg', 'TBar_tch_powheg', 'WWTo2L2Nu',  'WZTo3LNu','WZTo2L2Q', 'ZZTo4L', 'ZZTo2L2Nu', 'ZZTo2L2Q', 'WWW', 'WWZ', 'WZZ', 'ZZZ', 'TWZ', 'tZq_ll', 'TTZToLLNuNu' , 'TTZToQQ', 'TTLLJets_m1to10', 'TTWToLNu','TTWToQQ',  'TTTT', 'TTHnobb_pow', 'VHToNonbb',  'GGHZZ4L',  'WJetsToLNu_LO']
    mcDatasets += ttDatasets
    mcDatasets += DYDatasets

    daDatasetsB = ['DoubleEG_Run2016B_03Feb2017_ver2_v2_runs_273150_275376',
                   'DoubleMuon_Run2016B_03Feb2017_ver2_v2_runs_273150_275376',
                   'MuonEG_Run2016B_03Feb2017_ver2_v2_runs_273150_275376']    
                                                                              
                                                                              
    daDatasetsC = ['DoubleEG_Run2016C_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016C_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016C_03Feb2017_v1_runs_271036_284044']    
    
    daDatasetsD = ['DoubleEG_Run2016D_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016D_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016D_03Feb2017_v1_runs_271036_284044']    
                                                                              
    daDatasetsE = ['DoubleEG_Run2016E_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016E_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016E_03Feb2017_v1_runs_271036_284044']    
                                                                              
    daDatasetsF = ['DoubleEG_Run2016F_03Feb2017_v1_runs_271036_284044',
                  'DoubleMuon_Run2016F_03Feb2017_v1_runs_271036_284044',
                  'MuonEG_Run2016F_03Feb2017_v1_runs_271036_284044']  
                                                                              
    daDatasetsG = ['DoubleEG_Run2016G_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016G_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016G_03Feb2017_v1_runs_271036_284044']    
                                                                              
    daDatasetsH = ['DoubleEG_Run2016H_03Feb2017_ver2_v1_runs_281085_284035',
                   'DoubleEG_Run2016H_03Feb2017_ver3_v1_runs_284036_284044',
                   'DoubleMuon_Run2016H_03Feb2017_ver2_v1_runs_281085_284035',
                   'DoubleMuon_Run2016H_03Feb2017_ver3_v1_runs_284036_284044',
                   'MuonEG_Run2016H_03Feb2017_ver2_v1_runs_281085_284035', 
                   'MuonEG_Run2016H_03Feb2017_ver3_v1_runs_284036_284044']    

    daDatasets = daDatasetsB + daDatasetsC + daDatasetsD +daDatasetsE + daDatasetsF + daDatasetsG + daDatasetsH       


    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)
  
    print bcolors.HEADER + '[RinoutAnalysis] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    maxrun = 999999
    lumi = 36.8 ; maxrun = 999999; lumi_str = '36.8invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 

    cuts = CutManager.CutManager()
    
    runAnalysis(lumi, treeDA, treeMC, cuts, '', 'nocut', True, opts.ingredientsFile)
 

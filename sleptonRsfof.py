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
  

def makeTable(DATAee, DATAmm, DATASF, DATAOF, MCee, MCmm, MCSF, MCOF,DATArsfof, DATArsfof_e, MCrsfof, MCrsfof_e):
    line0 = '  \hline'
    line1 = '  ee &' 
    line2 = '  \mu\mu   &' 
    line3 = '  OF   &' 
    line4 = '  SF   &' 
    line5 = '  R_{\mathrm{SF/OF}}   &' 
    line0 += ' & Data  & MC   \\\\'
    line1 += ' %.f & %.f   \\\ ' %(DATAee.Integral(), MCee.Integral())
    line2 += ' %.f & %.f   \\\ ' %(DATAmm.Integral(), MCmm.Integral())
    line3 += ' %.f & %.f   \\\ ' %(DATAOF.Integral(), MCOF.Integral())
    line4 += ' %.f & %.f   \\\ ' %(DATASF.Integral(), MCSF.Integral())
    line5 += ' %.3f \pm %.3f & %.3f \pm %.3f  \\\ ' %(DATArsfof,DATArsfof_e, MCrsfof, MCrsfof_e)
    line0 += '\\hline'; line4 += '\\hline'; line5 += '\\hline';
                                                                                                                                                                                     
    helper.ensureDirectory('plots/rsfof/%s/'%lumi_str)
    helper.ensureDirectory('plots/rsfof/%s/tables/'%lumi_str)
    compTableFile = open('plots/rsfof/%s/tables/resultTable_%s%s.txt'%(lumi_str, str(lumi).replace('.','p'), "rt"),'w')
    compTableFile.write(line0+'\n')
    compTableFile.write(line1+'\n')
    compTableFile.write(line2+'\n')                                                                                             
    compTableFile.write(line3+'\n')                                                                                             
    compTableFile.write(line4+'\n')                                                                                             
    compTableFile.write(line5+'\n')                                                                                             
    print line0
    print line1
    print line2                                                                                                                                                                      
    print line3                                                                                                                                             
    print line4                                                                                                                                                      
    print line5                                                                                                                                                     


                                                                                                                                                                              
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

    ratio.GetYaxis().SetRangeUser(0,3)

    return ratio


def runAnalysis(lumi, treeDAEdge, treeMCEdge, treeDASlep, treeMCSlep, treeMCSlepFS, cuts, specialcut, tag, save, ingredientsFile):

    labelx = "m_{ll} [GeV]"
    labelmet = "p_{T}^{miss} [GeV]"
    labelmt2 = "M_{T2} [GeV]"
    labelnj = "N. Jets"
    labelpt = "p_{T}^{leading} [GeV]"
    labelHT = "HT [GeV]"

    #####Main mll control plot
    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Starting mll plot' + bcolors.ENDC
    bins = [20,40, 60, 76,106, 120, 150, 180, 210, 250, 300, 350, 400] 
    MCCRSFvalue  = treeMCEdge.getTH1F(lumi, "MCCRSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFCR, cuts.SF]), '', labelx, "1", kf)
    MCCROFvalue  = treeMCEdge.getTH1F(lumi, "MCCROFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFCR, cuts.OF]), '', labelx, "1", kf)
    DataCRSFvalue= treeDAEdge.getTH1F(lumi,"DataCRSFvalue","lepsMll_Edge",[20,1000], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFCR, cuts.SF]), '', labelx, "1", kf)
    DataCROFvalue= treeDAEdge.getTH1F(lumi,"DataCROFvalue","lepsMll_Edge",[20,1000], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFCR, cuts.OF]), '', labelx, "1", kf)    
    MCCRvalue   = make_rsfof(MCCRSFvalue, MCCROFvalue, "MC")
    DataCRvalue = make_rsfof(DataCRSFvalue, DataCROFvalue, "DA")

    measuredValueMC =         MCCRvalue.GetBinContent(1)
    measuredValueUncMC =      MCCRvalue.GetBinError(1)
    measuredValueData =       DataCRvalue.GetBinContent(1)
    if(measuredValueUncMC < measuredValueData*0.04):
        measuredValueUncMC = measuredValueData*0.04
    measuredValueUncData =    DataCRvalue.GetBinError(1)
    measuredValueUncTotData = math.sqrt(measuredValueUncData**2 + measuredValueUncMC**2)
    #print measuredValueUncTotData
    
    print "MC:   " ,measuredValueMC
    print "Data: " ,measuredValueData, " +- ",measuredValueUncTotData
   # makeTable(DataCRee, DataCRmm, DataCRSF, DataCROF, MCCRee, MCCRmm, MCCRSF, MCCROF,measuredValueData, measuredValueUncData, measuredValueMC, measuredValueUncMC)
    doMll = False
    if doMll:
        print "doing mll"
        MCCRSF       = treeMCEdge.getTH1F(lumi, "MCCSF", "lepsMll_Edge", bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFCR, cuts.SF]), '', labelx, "1", kf)
        MCCROF       = treeMCEdge.getTH1F(lumi, "MCCOF", "lepsMll_Edge", bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFCR, cuts.OF]), '', labelx, "1", kf)
        DataCRSF     = treeDAEdge.getTH1F(lumi,"DataCRSF","lepsMll_Edge",bins,1,1,cuts.AddList([specialcut,cuts.goodLepton,cuts.RSFOFCR, cuts.SF]), '', labelx, "1", kf)
        DataCROF     = treeDAEdge.getTH1F(lumi,"DataCROF","lepsMll_Edge",bins,1,1,cuts.AddList([specialcut,cuts.goodLepton,cuts.RSFOFCR, cuts.OF]), '', labelx, "1", kf)
        DataCR      = make_rsfof(DataCRSF, DataCROF, "DA")
        MCCR        = make_rsfof(MCCRSF, MCCROF, "MC")
        
        plot_rsfof = Canvas.Canvas('rsfof/%s_%s/plot_rsfof_mll'%(lumi_str,tag), 'png,pdf',  0.37, 0.6, 0.62, 0.8)
        plot_rsfof.addBand (MCCR.GetXaxis().GetXmin(), measuredValueData-measuredValueUncTotData, MCCR.GetXaxis().GetXmax(), measuredValueData+measuredValueUncTotData, r.kOrange-5, 0.4)
        plot_rsfof.addLine (MCCR.GetXaxis().GetXmin(), measuredValueData, MCCR.GetXaxis().GetXmax(), measuredValueData, r.kGreen-2)
        plot_rsfof.addLine(MCCR.GetXaxis().GetXmin(), 1, MCCR.GetXaxis().GetXmax(), 1, r.kBlack)
        plot_rsfof.addHisto(MCCR, 'E1', 'R_{SFOF} control region - MC', 'PL', r.kGreen-2 , 1, 0)
        plot_rsfof.addHisto(DataCR, 'E1,SAME', 'R_{SFOF} control region - Data', 'PL', r.kBlack , 1, 0)
        plot_rsfof.save(1, 1, 0, lumi, 0.8, 1.8)
    
    ######Met and JET plots 
    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Starting jet and met plot' + bcolors.ENDC
    doNJets = True
    if doNJets:    
     #   print "Doing NJets FS Only"
     #   MCCRSFnjFS =treeMCSlepFS.getTH1F(lumi,"MCCRSFnjFS","nJet25_Edge",[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFCRNoJet, cuts.SF]), '', labelnj)
     #   MCCROFnjFS =treeMCSlepFS.getTH1F(lumi,"MCCROFnjFS","nJet25_Edge",[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFCRNoJet, cuts.OF]), '', labelnj)
     #   DACRSFnjFS =treeDASlep.getTH1F(lumi,"DACRSFnjFS","nJet25_Edge",[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFCRNoJet, cuts.SF]), '', labelnj)
     #   DACROFnjFS =treeDASlep.getTH1F(lumi,"DACROFnjFS","nJet25_Edge",[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFCRNoJet, cuts.OF]), '', labelnj)
     #   MCCRnjFS   =make_rsfof(MCCRSFnjFS, MCCROFnjFS, "MC")
     #   DACRnjFS   =make_rsfof(DACRSFnjFS, DACROFnjFS, "DA")
     #                                                                                                                                                                        
     #   plot_rsfofnjFS = Canvas.Canvas('rsfof/%s_%s/plot_rsfof_njetFSonly'%(lumi_str,tag), 'png,pdf', 0.35, 0.6, 0.6, 0.8)
     #   plot_rsfofnjFS.addBand(MCCRnjFS.GetXaxis().GetXmin(), measuredValueData-measuredValueUncTotData, MCCRnjFS.GetXaxis().GetXmax(), measuredValueData+measuredValueUncTotData, r.kOrange-5, 0.4)
     #   plot_rsfofnjFS.addLine(MCCRnjFS.GetXaxis().GetXmin(), measuredValueData, MCCRnjFS.GetXaxis().GetXmax(), measuredValueData, r.kGreen-2)
     #   plot_rsfofnjFS.addLine(MCCRnjFS.GetXaxis().GetXmin(), 1, MCCRnjFS.GetXaxis().GetXmax(), 1, r.kBlack)
     #   plot_rsfofnjFS.addHisto(MCCRnjFS, 'E1', 'R_{SFOF} control region - FS MC', 'PL', r.kGreen-2 , 1, 0)
     #   plot_rsfofnjFS.addHisto(DACRnjFS, 'E1, SAME', 'R_{SFOF} control region - Data', 'PL', r.kBlack , 1, 0)
     #   plot_rsfofnjFS.save(1, 1, 0, lumi, 0.8, 1.8)                                                                                                                                                     
#  
    
        print "Doing NJets "
        MCCRSFnj =treeMCSlep.getTH1F(lumi,"MCCRSleptonSFnj","nJet25_Edge",[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFCRNoJet,cuts.MT2MET, cuts.SF]), '', "N. Jets w/ m_{ll} > 170 GeV", "1", kf)
        MCCROFnj =treeMCSlep.getTH1F(lumi,"MCCRSleptonOFnj","nJet25_Edge",[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFCRNoJet,cuts.MT2MET, cuts.OF]), '', "N. Jets w/ m_{ll} > 170 GeV", "1", kf)
        DACRSFnj =treeDASlep.getTH1F(lumi,"DACRSleptonSFnj","nJet25_Edge",[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFCRNoJet,cuts.MT2MET, cuts.SF]), '', "N. Jets w/ m_{ll} > 170 GeV", "1", kf)
        DACROFnj =treeDASlep.getTH1F(lumi,"DACRSleptonOFnj","nJet25_Edge",[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFCRNoJet,cuts.MT2MET, cuts.OF]), '', "N. Jets w/ m_{ll} > 170 GeV", "1", kf)
        MCCRnj   =make_rsfof(MCCRSFnj, MCCROFnj, "MC")
        DACRnj   =make_rsfof(DACRSFnj, DACROFnj, "DA")
                                                                                                                                                                             
        plot_rsfofnj = Canvas.Canvas('rsfof/%s_%s/plot_rsfof_njet'%(lumi_str,tag), 'png,pdf', 0.37, 0.6, 0.62, 0.8)
        plot_rsfofnj.addBand(MCCRnj.GetXaxis().GetXmin(), measuredValueData-measuredValueUncTotData, MCCRnj.GetXaxis().GetXmax(), measuredValueData+measuredValueUncTotData, r.kOrange-5, 0.4)
        plot_rsfofnj.addLine(MCCRnj.GetXaxis().GetXmin(), measuredValueData, MCCRnj.GetXaxis().GetXmax(), measuredValueData, r.kGreen-2)
        plot_rsfofnj.addLine(MCCRnj.GetXaxis().GetXmin(), 1, MCCRnj.GetXaxis().GetXmax(), 1, r.kBlack)
        plot_rsfofnj.addHisto(MCCRnj, 'E1', 'R_{SFOF} control region - MC', 'PL', r.kGreen-2 , 1, 0)
        plot_rsfofnj.addHisto(DACRnj, 'E1, SAME', 'R_{SFOF} control region - Data', 'PL', r.kBlack , 1, 0)
        plot_rsfofnj.save(1, 1, 0, lumi, 0.8, 1.8)                                                                                                                                                    
#    
    print "Doing MET"
    MCCRSFMET =treeMCEdge.getTH1F(lumi,"MCCRSFMET", "met_Edge", [100, 150, 200, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFCRNoMET, cuts.SF]),'', labelmet, "1", kf)
    MCCROFMET =treeMCEdge.getTH1F(lumi,"MCCROFMET", "met_Edge", [100, 150, 200, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFCRNoMET, cuts.OF]),'', labelmet, "1", kf)
    DACRSFMET =treeDAEdge.getTH1F(lumi,"DACRSFMET", "met_Edge", [100, 150, 200, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFCRNoMET, cuts.SF]),'', labelmet, "1", kf)
    DACROFMET =treeDAEdge.getTH1F(lumi,"DACROFMET", "met_Edge", [100, 150, 200, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFCRNoMET, cuts.OF]),'', labelmet, "1", kf)
    MCCRMET   =make_rsfof(MCCRSFMET, MCCROFMET, "MC")
    DACRMET   =make_rsfof(DACRSFMET, DACROFMET, "DA")
    plot_rsfofmet = Canvas.Canvas('rsfof/%s_%s/plot_rsfof_met'%(lumi_str,tag), 'png,pdf', 0.37, 0.6, 0.62, 0.8)
    plot_rsfofmet.addBand(MCCRMET.GetXaxis().GetXmin(), measuredValueData-measuredValueUncTotData, MCCRMET.GetXaxis().GetXmax(), measuredValueData+measuredValueUncTotData, r.kOrange-5, 0.4)
    plot_rsfofmet.addLine(MCCRMET.GetXaxis().GetXmin(), 1, MCCRMET.GetXaxis().GetXmax(), 1, r.kBlack)
    plot_rsfofmet.addLine(MCCRMET.GetXaxis().GetXmin(), measuredValueData, MCCRMET.GetXaxis().GetXmax(), measuredValueData, r.kGreen-2)
    plot_rsfofmet.addHisto(MCCRMET, 'E1', 'R_{SFOF} control region - MC', 'PL', r.kGreen-2 , 1, 0)
    plot_rsfofmet.addHisto(DACRMET, 'E1, SAME', 'R_{SFOF} control region - Data', 'PL', r.kBlack , 1, 0)
    plot_rsfofmet.save(1, 1, 0, lumi, 0.8, 1.8)                                                                                                                                                 

    
    print "Doing MT2"
    MCCRSFmt2 =treeMCEdge.getTH1F(lumi,"MCCRSFmt2","mt2_Edge",[0,20,40,60,80,100], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFCR, cuts.SF]),'', labelmt2, "1", kf)
    MCCROFmt2 =treeMCEdge.getTH1F(lumi,"MCCROFmt2","mt2_Edge",[0,20,40,60,80,100], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFCR, cuts.OF]),'', labelmt2, "1", kf)
    DACRSFmt2 =treeDAEdge.getTH1F(lumi,"DACRSFmt2","mt2_Edge",[0,20,40,60,80,100], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFCR, cuts.SF]), '', labelmt2, "1", kf)
    DACROFmt2 =treeDAEdge.getTH1F(lumi,"DACROFmt2","mt2_Edge",[0,20,40,60,80,100], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFCR, cuts.OF]), '', labelmt2, "1", kf)
    MCCRmt2   =make_rsfof(MCCRSFmt2, MCCROFmt2, "MC")
    DACRmt2   =make_rsfof(DACRSFmt2, DACROFmt2, "DA")

    plot_rsfofmt2 = Canvas.Canvas('rsfof/%s_%s/plot_rsfof_mt2'%(lumi_str,tag), 'png,pdf', 0.37, 0.6, 0.62, 0.8)
    plot_rsfofmt2.addBand(MCCRmt2.GetXaxis().GetXmin(), measuredValueData-measuredValueUncTotData, MCCRmt2.GetXaxis().GetXmax(), measuredValueData+measuredValueUncTotData, r.kOrange-5, 0.4)
    plot_rsfofmt2.addLine(MCCRmt2.GetXaxis().GetXmin(), measuredValueData, MCCRmt2.GetXaxis().GetXmax(), measuredValueData, r.kGreen-2)
    plot_rsfofmt2.addLine(MCCRmt2.GetXaxis().GetXmin(), 1, MCCRmt2.GetXaxis().GetXmax(), 1, r.kBlack)
    plot_rsfofmt2.addHisto(MCCRmt2, 'E1', 'R_{SFOF} control region - MC', 'PL', r.kGreen-2 , 1, 0)
    plot_rsfofmt2.addHisto(DACRmt2, 'E1, SAME', 'R_{SFOF} control region - Data', 'PL', r.kBlack , 1, 0)
    plot_rsfofmt2.save(1, 1, 0, lumi, 0.8, 1.8)                                                                                                                                                       



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

    #mcDatasets = ['ZZTo4L', 'WZTo3LNu', 'WWW', 'WWZ',  'WZZ', 'ZZZ', 'DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO', 'WWTo2L2Nu', 'TTTo2L2Nu', 'ZZTo2L2Nu', 'TTHnobb_pow', 'VHToNonbb', 'TWZ', 'WZTo2L2Q',  'TBar_tch_powheg', 'T_tch_powheg', 'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT', 'TTTT',  'TTZToQQ', 'TTWToQQ', 'TTLLJets_m1to10', 'TTWToLNu_ext2']

    dyDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    zzDatasets = ['ZZTo2L2Nu']
    wzDatasets = ['WZTo3LNu']
    ttzDatasets = ['TTZToLLNuNu_ext2']
    othersDatasets = ['WWZ', 'WZZ', 'ZZZ', 'TWZ', 'tZq_ll']
    fsDatasets = ['TTTT','TTTo2L2Nu', 'TBar_tch_powheg', 'T_tch_powheg', 'WWTo2L2Nu','GGWWTo2L2Nu', 'WWW','TTWToLNu_ext2',  'TTWToQQ', 'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT']                       
    kf = "noKFactor" 
    #mcDatasets = fsDatasets 
    mcDatasets = fsDatasets+dyDatasets + othersDatasets + zzDatasets + wzDatasets + ttzDatasets
    
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


    #daDatasets = daDatasetsE
    daDatasets = daDatasetsB + daDatasetsC + daDatasetsD +daDatasetsE + daDatasetsF + daDatasetsG + daDatasetsH
    
    treeMCEdge = Sample.Tree(helper.selectSamples('samplesEdge.dat', mcDatasets, 'MC'), 'MC'  , 0, isScan = 0, isOnEOS = 0)
    treeMCSlepton = Sample.Tree(helper.selectSamples('samplesUnskimmedEOS.dat', mcDatasets, 'MC'), 'MC'  , 0, isScan = 0, isOnEOS = 1)
    treeMCSleptonFS = Sample.Tree(helper.selectSamples('samplesUnskimmedEOS.dat', fsDatasets, 'MC'), 'MC'  , 0, isScan = 0, isOnEOS = 1)
    treeDAEdge = Sample.Tree(helper.selectSamples('samplesEdge.dat', daDatasets, 'DA'), 'DATA', 1, isScan = 0, isOnEOS = 0)
    treeDASlepton = Sample.Tree(helper.selectSamples('samplesUnskimmedEOS.dat', daDatasets, 'DA'), 'DATA', 1, isScan = 0, isOnEOS = 1)
  
    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    maxrun = 999999
    #lumi = 4.1 ; maxrun = 999999; lumi_str = '4.1'
    lumi = 35.9 ; maxrun = 999999; lumi_str = '35.9'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 

    cuts = CutManager.CutManager()
    
    runAnalysis(lumi, treeDAEdge, treeMCEdge, treeDASlepton, treeMCSlepton, treeMCSleptonFS,cuts, '', 'nocut', True, opts.ingredientsFile)
 

import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TH1F, TGraph
import math, sys, optparse, copy, re, array, itertools, os


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables




if __name__ == "__main__":
   

    gROOT.SetBatch(True)
    gROOT.ProcessLine('.L include/tdrstyle.C')
    r.setTDRStyle() 
    r.gStyle.SetPaintTextFormat("4.2f")
    
    
    global lumi
    lumi= 12.9
    mcDatasets = ['TT_pow_ext34']
    siDatasets = ['SMS_T6bbllslepton_mSbottom600to775_mLSP150to725',
                  'SMS_T6bbllslepton_mSbottom400to575_mLSP150to550',
                  'SMS_T6bbllslepton_mSbottom800to950_mLSP150to900']
    treeMC = Sample.Tree(helper.selectSamples('samples.dat', mcDatasets, 'MC'), 'MC', 0, isScan = 0)
    treeSI = Sample.Tree(helper.selectSamples('samples.dat', siDatasets, 'SI'), 'SI', 0, isScan = 1)
    cuts = CutManager.CutManager()

    theCutsNum = [cuts.SignalRegionBaseLineNoTrigger, cuts.SF]
    theCutsDen = [cuts.goodLepton,cuts.SF]

    for pu in ['greater','smaller']:
        puCut = ['nVert_Edge > 20'] if 'greater' in pu else ['nVert_Edge < 20']
        denominador  = treeSI.getTH2F(lumi,"dem","GenSusyMScan2_Edge:GenSusyMScan1_Edge",20,400.-25,950+25,32,125.,925.,cuts.AddList(theCutsDen+puCut),"","m_{#tilde{b}} [GeV]","m_{#tilde{#chi}_{2}^{0}} [GeV]")
        numerador    = treeSI.getTH2F(lumi,"num","GenSusyMScan2_Edge:GenSusyMScan1_Edge",20,400.-25,950+25,32,125.,925.,cuts.AddList(theCutsNum+puCut),"","m_{#tilde{b}} [GeV]","m_{#tilde{#chi}_{2}^{0}} [GeV]")
        numerador.SetTitle('n_{Vert} > 20' if 'greater' in pu else 'n_{Vert} > 20')
        numerador.Divide(denominador)
        numerador.GetZaxis().SetRangeUser(0.,1.)
        c = r.TCanvas()
        numerador.Draw('colz text')
        c.SaveAs('plots/EffVsPU/%s.pdf'%pu)
        c.SaveAs('plots/EffVsPU/%s.png'%pu)

    for mp in [[900,200],[900,700],[500,150],[500,400]]:
        high    = treeSI.getTH1F(lumi,"high","lepsMll_Edge",50,0.,1000.,cuts.AddList(theCutsNum+['nVert_Edge > 20']+['GenSusyMScan1_Edge == %.0f && GenSusyMScan2_Edge == %.0f'%(mp[0], mp[1])]),"","m_{ll} [GeV]",ylabel="A.U.")
        small   = treeSI.getTH1F(lumi,"smal","lepsMll_Edge",50,0.,1000.,cuts.AddList(theCutsNum+['nVert_Edge < 20']+['GenSusyMScan1_Edge == %.0f && GenSusyMScan2_Edge == %.0f'%(mp[0], mp[1])]),"","m_{ll} [GeV]",ylabel="A.U.")
        high.Scale(1/high.Integral())
        small.Scale(1/small.Integral())
        high.SetLineWidth(2); small.SetLineWidth(2)
        plot = Canvas.Canvas('EffVsPU/ttbar_sig_%d_%d'%(mp[0],mp[1]),'png,pdf',0.61,0.6,0.71,0.72)
        plot.addHisto(high, 'HIST',"n_{Vertices} > 20","L",r.kRed, 1,0)        
        plot.addHisto(small,'HIST,SAME',"n_{Vertices} < 20","L",r.kBlue,1,0)
        plot.addLatex(0.61,0.80,'m_{#tilde{b}} = %d GeV'%(mp[0]))
        plot.addLatex(0.61,0.75,'m_{LSP} = %d GeV'%(mp[1]))

        plot.save(1,True,False,lumi)

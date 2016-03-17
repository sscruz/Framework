import os
import sys
import ROOT as r

from include.setTDRStyle import setTDRStyle
from include.helper import createMyColors
from include.helper import myColors
from ROOT import TCanvas, TH1F, TPad, THStack, TGraphErrors, TLatex, TLine
import math
	
def makeOverviewPlotWithOnZ():

	colors = createMyColors()	

	
	
	histObs = TH1F("histObs","histObs",17,0,17)
	
	histObs.SetMarkerColor(r.kBlack)
	histObs.SetLineColor(r.kBlack)
	histObs.SetMarkerStyle(20)
	
	histTotal = TH1F("histPred","histPred",17,0,17)
	histFlavSym = TH1F("histFlavSym","histFlavSym",17,0,17)
	histDY = TH1F("histDY","histDY",17,0,17)
	histMC = TH1F("histMC","histMC",17,0,17)
	
	hCanvas = TCanvas("hCanvas", "Distribution", 800,800)
	
	plotPad = TPad("plotPad","plotPad",0,0,1,1)
	style=setTDRStyle()
	style.SetPadBottomMargin(0.3)
	plotPad.UseCurrentStyle()
	plotPad.Draw()	
	plotPad.cd()	
	



	### on-Z numbers from PAS draft

	histObs.SetBinContent(1,28)
	histObs.SetBinContent(2,7)
	histObs.SetBinContent(3,6)		
	histObs.SetBinContent(4,6)

	histObs.SetBinContent(5,21)
	histObs.SetBinContent(6,6)
	histObs.SetBinContent(7,1)		
	histObs.SetBinContent(8,3)

	histObs.SetBinContent(9,20)
	histObs.SetBinContent(10,10)
	histObs.SetBinContent(11,2)		
	histObs.SetBinContent(12,0)

	histObs.SetBinContent(13,45)
	histObs.SetBinContent(14,23)
	histObs.SetBinContent(15,4)		
	histObs.SetBinContent(16,3)
	
	histObs.SetBinContent(17,14)

	
	#~ names = ["low-Mass central","below-Z central","on-Z central","above-Z central","high-Mass central","low-Mass forward","below-Z forward","on-Z forward","above-Z forward","high-Mass forward","low-Mass central","below-Z central","on-Z central","above-Z central","high-Mass central","low-Mass forward","below-Z forward","on-Z forward","above-Z forward","high-Mass forward","low-Mass central","below-Z central","on-Z central","above-Z central","high-Mass central","low-Mass forward","below-Z forward","on-Z forward","above-Z forward","high-Mass forward"]
	#~ names = ["#geq 0 b-tags c","= 0 b-tags c","#geq 1 b-tags c","#geq 0 b-tags f","= 0 b-tags f","#geq 1 b-tags f","#geq 0 b-tags c","= 0 b-tags c","#geq 1 b-tags c","#geq 0 b-tags f","= 0 b-tags f","#geq 1 b-tags f","#geq 0 b-tags c","= 0 b-tags c","#geq 1 b-tags c","#geq 0 b-tags f","= 0 b-tags f","#geq 1 b-tags f","#geq 0 b-tags c","= 0 b-tags c","#geq 1 b-tags c","#geq 0 b-tags f","= 0 b-tags f","#geq 1 b-tags f","#geq 0 b-tags c","#geq 1 b-tags c","= 0 b-tags c","#geq 0 b-tags f","= 0 b-tags f","#geq 1 b-tags f"]
	names = ["E_{T}^{miss} 100-150 GeV","E_{T}^{miss} 150-225 GeV","E_{T}^{miss} 225-300 GeV","E_{T}^{miss} > 300 GeV","E_{T}^{miss} 100-150 GeV","E_{T}^{miss} 150-225 GeV","E_{T}^{miss} 225-300 GeV","E_{T}^{miss} > 300 GeV","E_{T}^{miss} 100-150 GeV","E_{T}^{miss} 150-225 GeV","E_{T}^{miss} 225-300 GeV","E_{T}^{miss} > 300 GeV","E_{T}^{miss} 100-150 GeV","E_{T}^{miss} 150-225 GeV","E_{T}^{miss} 225-300 GeV","E_{T}^{miss} > 300 GeV","ATLAS SR"]
	
	for index, name in enumerate(names):
	
		histObs.GetXaxis().SetBinLabel(index+1,name)
	
#~ 
	histTotal.SetBinContent(1, 29.1)
	histTotal.SetBinContent(2, 9.1)
	histTotal.SetBinContent(3, 3.4)
	histTotal.SetBinContent(4, 2.1)

	histTotal.SetBinContent(5, 14.3)
	histTotal.SetBinContent(6, 6.9)
	histTotal.SetBinContent(7, 6.1)
	histTotal.SetBinContent(8, 1.5)
	
	histTotal.SetBinContent(9, 23.6)
	histTotal.SetBinContent(10, 8.2)
	histTotal.SetBinContent(11, 0.8)
	histTotal.SetBinContent(12, 1.5)

	histTotal.SetBinContent(13, 44.7)
	histTotal.SetBinContent(14, 16.8)
	histTotal.SetBinContent(15, 0.6)
	histTotal.SetBinContent(16, 1.5)
	
	histTotal.SetBinContent(17, 12.3)	

	histDY.SetBinContent(1, 24.3)
	histDY.SetBinContent(2, 4.6)
	histDY.SetBinContent(3, 1.5)
	histDY.SetBinContent(4, 1.1)

	histDY.SetBinContent(5, 4.5)
	histDY.SetBinContent(6, 1.4)
	histDY.SetBinContent(7, 0.7)
	histDY.SetBinContent(8, 0.2)
	
	histDY.SetBinContent(9, 10.0)
	histDY.SetBinContent(10, 3.2)
	histDY.SetBinContent(11, 0.3)
	histDY.SetBinContent(12, 0.1)

	histDY.SetBinContent(13, 5.0)
	histDY.SetBinContent(14, 1.6)
	histDY.SetBinContent(15, 0.4)
	histDY.SetBinContent(16, 0.3)
	
	histDY.SetBinContent(17, 3.9)
	
	histFlavSym.SetBinContent(1,3.2)
	histFlavSym.SetBinContent(2,3.2)
	histFlavSym.SetBinContent(3,1.1)
	histFlavSym.SetBinContent(4,0)

	histFlavSym.SetBinContent(5,9.5)
	histFlavSym.SetBinContent(6,5.3)
	histFlavSym.SetBinContent(7,5.3)
	histFlavSym.SetBinContent(8,1.1)
	
	histFlavSym.SetBinContent(9,12.6)
	histFlavSym.SetBinContent(10,4.2)
	histFlavSym.SetBinContent(11,0)
	histFlavSym.SetBinContent(12,1.1)

	histFlavSym.SetBinContent(13,38.9)
	histFlavSym.SetBinContent(14,14.7)
	histFlavSym.SetBinContent(15,0.0)
	histFlavSym.SetBinContent(16,1.1)
	
	histFlavSym.SetBinContent(17,6.13)	


	histMC.SetBinContent(1,1.6)
	histMC.SetBinContent(2,1.3)
	histMC.SetBinContent(3,0.8)
	histMC.SetBinContent(4,1)

	histMC.SetBinContent(5,0.3)
	histMC.SetBinContent(6,0.2)
	histMC.SetBinContent(7,0.1)
	histMC.SetBinContent(8,0.2)
	
	histMC.SetBinContent(9,1)
	histMC.SetBinContent(10,0.8)
	histMC.SetBinContent(11,0.5)
	histMC.SetBinContent(12,0.3)

	histMC.SetBinContent(13,0.8)
	histMC.SetBinContent(14,0.5)
	histMC.SetBinContent(15,0.2)
	histMC.SetBinContent(16,0.1)
	
	histMC.SetBinContent(17,2.1)	

	
	
	
	errGraph = r.TGraphAsymmErrors()
	
	for i in range(1,histFlavSym.GetNbinsX()+1):
		errGraph.SetPoint(i,i-0.5,histTotal.GetBinContent(i))
		

	errGraph.SetPointError(1,0.5,0.5,4.7,5.3)
	errGraph.SetPointError(2,0.5,0.5,1.9,3.2)
	errGraph.SetPointError(3,0.5,0.5,1.0,2.5)
	errGraph.SetPointError(4,0.5,0.5,0.7,1.4)

	errGraph.SetPointError(5,0.5,0.5,3.2,4.4)
	errGraph.SetPointError(6,0.5,0.5,2.3,3.6)
	errGraph.SetPointError(7,0.5,0.5,2.3,3.6)
	errGraph.SetPointError(8,0.5,0.5,0.9,2.4)

	errGraph.SetPointError(9,0.5,0.5,3.7,4.9)
	errGraph.SetPointError(10,0.5,0.5,2.1,3.4)
	errGraph.SetPointError(11,0.5,0.5,0.2,1.2)
	errGraph.SetPointError(12,0.5,0.5,0.9,2.4)

	errGraph.SetPointError(13,0.5,0.5,6.6,7.7)
	errGraph.SetPointError(14,0.5,0.5,3.9,5.1)
	errGraph.SetPointError(15,0.5,0.5,0.3,1.2)
	errGraph.SetPointError(16,0.5,0.5,0.9,2.4)

	errGraph.SetPointError(17,0.5,0.5,2.8,4.0)	

	errGraph.SetFillColor(myColors["MyBlue"])
	errGraph.SetFillStyle(3001)	

	histFlavSym.SetLineColor(r.kBlue+3)
	histFlavSym.SetLineWidth(2)
	
	histDY.SetLineColor(r.kGreen+3)
	histDY.SetFillColor(r.kGreen+3)
	histDY.SetFillStyle(3002)
	
	histMC.SetLineColor(r.kViolet+3)
	histMC.SetFillColor(r.kViolet+3)
	histMC.SetFillStyle(3002)

	
	
	#~ histDY.Add(histMC)
	stack = THStack()
	stack.Add(histMC)	
	stack.Add(histDY)	
	stack.Add(histFlavSym)

	
	
	histObs.GetYaxis().SetRangeUser(0,100)
	histObs.GetYaxis().SetTitle("Events")
	histObs.LabelsOption("v")

	histObs.UseCurrentStyle()
	histObs.Draw("pe")

	
	
	#~ hCanvas.DrawFrame(-0.5,0,30.5,65,"; %s ; %s" %("","Events"))
	
	latex = TLatex()
	latex.SetTextFont(42)
	latex.SetTextAlign(31)
	latex.SetTextSize(0.04)
	latex.SetNDC(True)
	latexCMS = TLatex()
	latexCMS.SetTextFont(61)
	#latexCMS.SetTextAlign(31)
	latexCMS.SetTextSize(0.06)
	latexCMS.SetNDC(True)
	latexCMSExtra = TLatex()
	latexCMSExtra.SetTextFont(52)
	#latexCMSExtra.SetTextAlign(31)
	latexCMSExtra.SetTextSize(0.045)
	latexCMSExtra.SetNDC(True)		
	


	intlumi = TLatex()
	intlumi.SetTextAlign(12)
	intlumi.SetTextSize(0.03)
	intlumi.SetNDC(True)		

	latex.DrawLatex(0.95, 0.96, "%s fb^{-1} (13 TeV)"%"2.3")
	
	cmsExtra = "Preliminary"
	latexCMS.DrawLatex(0.19,0.88,"CMS")
	if "Simulation" in cmsExtra:
		yLabelPos = 0.81	
	else:
		yLabelPos = 0.84	

	latexCMSExtra.DrawLatex(0.19,yLabelPos,"%s"%(cmsExtra))

	
	leg1 = r.TLegend(0.4, 0.85, 0.925, 0.95,"","brNDC")
	leg1.SetNColumns(2)
	leg1.SetFillColor(10)
	leg1.SetLineColor(10)
	leg1.SetShadowColor(0)
	leg1.SetBorderSize(1)
	
	
	leg1.AddEntry(histObs,"Data    ","pe")
	leg1.AddEntry(histTotal, "Total backgrounds","l")
	
	
	leg2 = r.TLegend(0.41, 0.75, 0.925, 0.85,"","brNDC")
	leg2.SetNColumns(3)
	leg2.SetFillColor(10)
	leg2.SetLineColor(10)
	leg2.SetShadowColor(0)
	leg2.SetBorderSize(1)
	
	leg2.AddEntry(errGraph,"Total uncert. ", "f")	
	leg2.AddEntry(histDY,"Z-jets  ", "f")
	leg2.AddEntry(histMC,"Rare", "f")
	
	

	stack.Draw("samehist")
	errGraph.Draw("same02")
	
	histObs.Draw("pesame")
	
	leg1.Draw("same")
	leg2.Draw("same")

	
	
	line1 = TLine(8,0,8,45)
	line2 = TLine(16,0,16,45)
	line3 = TLine(18,0,18,45)
	line4 = TLine(24,0,24,45)

	line1.SetLineColor(r.kBlack)
	line2.SetLineColor(r.kBlack)
	line3.SetLineColor(r.kBlack)
	line4.SetLineColor(r.kBlack)

	line1.SetLineWidth(2)
	line2.SetLineWidth(2)
	line3.SetLineWidth(4)
	line4.SetLineWidth(2)

	line1.Draw("same")
	line2.Draw("same")
	
	line5 = TLine(4,0,4,35)
	line6 = TLine(12,0,12,35)
	line7 = TLine(15,0,15,180)
	line8 = TLine(21,0,21,180)
	line9 = TLine(27,0,27,180)

	line5.SetLineColor(r.kBlack)
	line6.SetLineColor(r.kBlack)
	line7.SetLineColor(r.kBlack)
	line8.SetLineColor(r.kBlack)
	line9.SetLineColor(r.kBlack)

	line5.SetLineWidth(2)
	line6.SetLineWidth(2)
	line7.SetLineWidth(2)
	line8.SetLineWidth(2)
	line9.SetLineWidth(2)
	line5.SetLineStyle(r.kDashed)
	line6.SetLineStyle(r.kDashed)
	line7.SetLineStyle(r.kDashed)
	line8.SetLineStyle(r.kDashed)
	line9.SetLineStyle(r.kDashed)

	line5.Draw("same")
	line6.Draw("same")


	label = TLatex()
	label.SetTextAlign(12)
	label.SetTextSize(0.04)
	label.SetTextColor(r.kBlack)	
	#~ label.SetTextAngle(45)	
	
	label.DrawLatex(2.,60,"#splitline{N_{jets} = 2-3}{H_{T} > 400 GeV}")
	label.DrawLatex(10.5,60,"N_{jets} #geq 4")
	
	label = TLatex()
	label.SetTextAlign(12)
	label.SetTextSize(0.04)
	label.SetTextColor(r.kBlack)	
	label.SetTextAngle(45)	
	
	label.DrawLatex(1.5,30,"N_{b} = 0")
	label.DrawLatex(5.5,30,"N_{b} #geq 1")
	label.DrawLatex(9.5,30,"N_{b} = 0")
	label.DrawLatex(13.5,30,"N_{b} #geq 1")



	plotPad.RedrawAxis()
	
	hCanvas.Print("onZOverview.pdf")
	hCanvas.Print("onZOverview.root")


def makeOverviewPlotEdge():

        obs_cen = []
        obs_fwd = [] 
        tot_cen = []
        tot_fwd = []
        dy_cen = []
        dy_fwd = []
        other_cen = []
        other_fwd = []
        lines = open("resultsEdge.txt").readlines()
        for l in lines:
          s = l.split()
          obs_cen.append(float(s[0]))
          obs_fwd.append(float(s[1]))
          tot_cen.append(float(s[2]))
          tot_fwd.append(float(s[3]))
          dy_cen.append(float(s[4]))
          dy_fwd.append(float(s[5]))
          other_cen.append(float(s[6]))
          other_fwd.append(float(s[7]))
  

	colors = createMyColors()	
	
	histObs = TH1F("histObs","histObs",30,0,30)
	
	histObs.SetMarkerColor(r.kBlack)
	histObs.SetLineColor(r.kBlack)
	histObs.SetMarkerStyle(20)
	
	histPred = TH1F("histPred","histPred",30,0,30)
	histFlavSym = TH1F("histFlavSym","histFlavSym",30,0,30)
	histDY = TH1F("histDY","histDY",30,0,30)
	histMC = TH1F("histMC","histMC",30,0,30)
	
	hCanvas = TCanvas("hCanvas", "Distribution", 800,800)
	
	plotPad = TPad("plotPad","plotPad",0,0,1,1)
	style=setTDRStyle()
	style.SetPadBottomMargin(0.3)
	plotPad.UseCurrentStyle()
	plotPad.Draw()	
	plotPad.cd()	
		### on-Z numbers from PAS draft
        histObs.SetBinContent(1,obs_cen[0])
	histObs.SetBinContent(2,obs_cen[2])
	histObs.SetBinContent(3,obs_cen[4])		
	histObs.SetBinContent(4,obs_fwd[0])
	histObs.SetBinContent(5,obs_fwd[2])
	histObs.SetBinContent(6,obs_fwd[4])
	histObs.SetBinContent(7,obs_cen[6])
	histObs.SetBinContent(8,obs_cen[8])
	histObs.SetBinContent(9,obs_cen[10])		
	histObs.SetBinContent(0,obs_fwd[6])
	histObs.SetBinContent(11,obs_fwd[8])
	histObs.SetBinContent(12,obs_fwd[10])
	histObs.SetBinContent(13,obs_cen[12])
	histObs.SetBinContent(14,obs_cen[14])
	histObs.SetBinContent(15,obs_cen[16])		
	histObs.SetBinContent(16,obs_fwd[12])
	histObs.SetBinContent(17,obs_fwd[14])
	histObs.SetBinContent(18,obs_fwd[16])
	histObs.SetBinContent(19,obs_cen[18])
	histObs.SetBinContent(20,obs_cen[20])
	histObs.SetBinContent(21,obs_cen[22])		
	histObs.SetBinContent(22,obs_fwd[18])
	histObs.SetBinContent(23,obs_fwd[20])
	histObs.SetBinContent(24,obs_fwd[22])
	histObs.SetBinContent(25,obs_cen[24])
	histObs.SetBinContent(26,obs_cen[26])
	histObs.SetBinContent(27,obs_cen[28])		
	histObs.SetBinContent(28,obs_fwd[24])
	histObs.SetBinContent(29,obs_fwd[26])
	histObs.SetBinContent(30,obs_fwd[28])

        histPred.SetBinContent(1,tot_cen[0])
	histPred.SetBinContent(2,tot_cen[2])
	histPred.SetBinContent(3,tot_cen[4])		
	histPred.SetBinContent(4,tot_fwd[0])
	histPred.SetBinContent(5,tot_fwd[2])
	histPred.SetBinContent(6,tot_fwd[4])
	histPred.SetBinContent(7,tot_cen[6])
	histPred.SetBinContent(8,tot_cen[8])
	histPred.SetBinContent(9,tot_cen[10])		
	histPred.SetBinContent(0,tot_fwd[6])
	histPred.SetBinContent(11,tot_fwd[8])
	histPred.SetBinContent(12,tot_fwd[10])
	histPred.SetBinContent(13,tot_cen[12])
	histPred.SetBinContent(14,tot_cen[14])
	histPred.SetBinContent(15,tot_cen[16])		
	histPred.SetBinContent(16,tot_fwd[12])
	histPred.SetBinContent(17,tot_fwd[14])
	histPred.SetBinContent(18,tot_fwd[16])
	histPred.SetBinContent(19,tot_cen[18])
	histPred.SetBinContent(20,tot_cen[20])
	histPred.SetBinContent(21,tot_cen[22])		
	histPred.SetBinContent(22,tot_fwd[18])
	histPred.SetBinContent(23,tot_fwd[20])
	histPred.SetBinContent(24,tot_fwd[22])
	histPred.SetBinContent(25,tot_cen[24])
	histPred.SetBinContent(26,tot_cen[26])
	histPred.SetBinContent(27,tot_cen[28])		
	histPred.SetBinContent(28,tot_fwd[24])
	histPred.SetBinContent(29,tot_fwd[26])
	histPred.SetBinContent(30,tot_fwd[28])

        histPred.SetBinError(1,tot_cen[1])
	histPred.SetBinError(2,tot_cen[3])
	histPred.SetBinError(3,tot_cen[5])		
	histPred.SetBinError(4,tot_fwd[1])
	histPred.SetBinError(5,tot_fwd[3])
	histPred.SetBinError(6,tot_fwd[5])
	histPred.SetBinError(7,tot_cen[7])
	histPred.SetBinError(8,tot_cen[9])
	histPred.SetBinError(9,tot_cen[11])		
	histPred.SetBinError(0,tot_fwd[7])
	histPred.SetBinError(11,tot_fwd[9])
	histPred.SetBinError(12,tot_fwd[11])
	histPred.SetBinError(13,tot_cen[13])
	histPred.SetBinError(14,tot_cen[15])
	histPred.SetBinError(15,tot_cen[17])		
	histPred.SetBinError(16,tot_fwd[13])
	histPred.SetBinError(17,tot_fwd[15])
	histPred.SetBinError(18,tot_fwd[17])
	histPred.SetBinError(19,tot_cen[19])
	histPred.SetBinError(20,tot_cen[21])
	histPred.SetBinError(21,tot_cen[23])		
	histPred.SetBinError(22,tot_fwd[19])
	histPred.SetBinError(23,tot_fwd[21])
	histPred.SetBinError(24,tot_fwd[23])
	histPred.SetBinError(25,tot_cen[25])
	histPred.SetBinError(26,tot_cen[27])
	histPred.SetBinError(27,tot_cen[29])		
	histPred.SetBinError(28,tot_fwd[25])
	histPred.SetBinError(29,tot_fwd[27])
	histPred.SetBinError(30,tot_fwd[29])
 
        histDY.SetBinContent(1,dy_cen[0])
	histDY.SetBinContent(2,dy_cen[2])
	histDY.SetBinContent(3,dy_cen[4])		
	histDY.SetBinContent(4,dy_fwd[0])
	histDY.SetBinContent(5,dy_fwd[2])
	histDY.SetBinContent(6,dy_fwd[4])
	histDY.SetBinContent(7,dy_cen[6])
	histDY.SetBinContent(8,dy_cen[8])
	histDY.SetBinContent(9,dy_cen[10])		
	histDY.SetBinContent(0,dy_fwd[6])
	histDY.SetBinContent(11,dy_fwd[8])
	histDY.SetBinContent(12,dy_fwd[10])
	histDY.SetBinContent(13,dy_cen[12])
	histDY.SetBinContent(14,dy_cen[14])
	histDY.SetBinContent(15,dy_cen[16])		
	histDY.SetBinContent(16,dy_fwd[12])
	histDY.SetBinContent(17,dy_fwd[14])
	histDY.SetBinContent(18,dy_fwd[16])
	histDY.SetBinContent(19,dy_cen[18])
	histDY.SetBinContent(20,dy_cen[20])
	histDY.SetBinContent(21,dy_cen[22])		
	histDY.SetBinContent(22,dy_fwd[18])
	histDY.SetBinContent(23,dy_fwd[20])
	histDY.SetBinContent(24,dy_fwd[22])
	histDY.SetBinContent(25,dy_cen[24])
	histDY.SetBinContent(26,dy_cen[26])
	histDY.SetBinContent(27,dy_cen[28])		
	histDY.SetBinContent(28,dy_fwd[24])
	histDY.SetBinContent(29,dy_fwd[26])
	histDY.SetBinContent(30,dy_fwd[28])

        histDY.SetBinError(1,dy_cen[1])
	histDY.SetBinError(2,dy_cen[3])
	histDY.SetBinError(3,dy_cen[5])		
	histDY.SetBinError(4,dy_fwd[1])
	histDY.SetBinError(5,dy_fwd[3])
	histDY.SetBinError(6,dy_fwd[5])
	histDY.SetBinError(7,dy_cen[7])
	histDY.SetBinError(8,dy_cen[9])
	histDY.SetBinError(9,dy_cen[11])		
	histDY.SetBinError(0,dy_fwd[7])
	histDY.SetBinError(11,dy_fwd[9])
	histDY.SetBinError(12,dy_fwd[11])
	histDY.SetBinError(13,dy_cen[13])
	histDY.SetBinError(14,dy_cen[15])
	histDY.SetBinError(15,dy_cen[17])		
	histDY.SetBinError(16,dy_fwd[13])
	histDY.SetBinError(17,dy_fwd[15])
	histDY.SetBinError(18,dy_fwd[17])
	histDY.SetBinError(19,dy_cen[19])
	histDY.SetBinError(20,dy_cen[21])
	histDY.SetBinError(21,dy_cen[23])		
	histDY.SetBinError(22,dy_fwd[19])
	histDY.SetBinError(23,dy_fwd[21])
	histDY.SetBinError(24,dy_fwd[23])
	histDY.SetBinError(25,dy_cen[25])
	histDY.SetBinError(26,dy_cen[27])
	histDY.SetBinError(27,dy_cen[29])		
	histDY.SetBinError(28,dy_fwd[25])
	histDY.SetBinError(29,dy_fwd[27])
	histDY.SetBinError(30,dy_fwd[29])
	
        histMC.SetBinContent(1,other_cen[0])
	histMC.SetBinContent(2,other_cen[2])
	histMC.SetBinContent(3,other_cen[4])		
	histMC.SetBinContent(4,other_fwd[0])
	histMC.SetBinContent(5,other_fwd[2])
	histMC.SetBinContent(6,other_fwd[4])
	histMC.SetBinContent(7,other_cen[6])
	histMC.SetBinContent(8,other_cen[8])
	histMC.SetBinContent(9,other_cen[10])		
	histMC.SetBinContent(0,other_fwd[6])
	histMC.SetBinContent(11,other_fwd[8])
	histMC.SetBinContent(12,other_fwd[10])
	histMC.SetBinContent(13,other_cen[12])
	histMC.SetBinContent(14,other_cen[14])
	histMC.SetBinContent(15,other_cen[16])		
	histMC.SetBinContent(16,other_fwd[12])
	histMC.SetBinContent(17,other_fwd[14])
	histMC.SetBinContent(18,other_fwd[16])
	histMC.SetBinContent(19,other_cen[18])
	histMC.SetBinContent(20,other_cen[20])
	histMC.SetBinContent(21,other_cen[22])		
	histMC.SetBinContent(22,other_fwd[18])
	histMC.SetBinContent(23,other_fwd[20])
	histMC.SetBinContent(24,other_fwd[22])
	histMC.SetBinContent(25,other_cen[24])
	histMC.SetBinContent(26,other_cen[26])
	histMC.SetBinContent(27,other_cen[28])		
	histMC.SetBinContent(28,other_fwd[24])
	histMC.SetBinContent(29,other_fwd[26])
	histMC.SetBinContent(30,other_fwd[28])
	
        histMC.SetBinError(1,other_cen[1])
	histMC.SetBinError(2,other_cen[3])
	histMC.SetBinError(3,other_cen[5])		
	histMC.SetBinError(4,other_fwd[1])
	histMC.SetBinError(5,other_fwd[3])
	histMC.SetBinError(6,other_fwd[5])
	histMC.SetBinError(7,other_cen[7])
	histMC.SetBinError(8,other_cen[9])
	histMC.SetBinError(9,other_cen[11])		
	histMC.SetBinError(0,other_fwd[7])
	histMC.SetBinError(11,other_fwd[9])
	histMC.SetBinError(12,other_fwd[11])
	histMC.SetBinError(13,other_cen[13])
	histMC.SetBinError(14,other_cen[15])
	histMC.SetBinError(15,other_cen[17])		
	histMC.SetBinError(16,other_fwd[13])
	histMC.SetBinError(17,other_fwd[15])
	histMC.SetBinError(18,other_fwd[17])
	histMC.SetBinError(19,other_cen[19])
	histMC.SetBinError(20,other_cen[21])
	histMC.SetBinError(21,other_cen[23])		
	histMC.SetBinError(22,other_fwd[19])
	histMC.SetBinError(23,other_fwd[21])
	histMC.SetBinError(24,other_fwd[23])
	histMC.SetBinError(25,other_cen[25])
	histMC.SetBinError(26,other_cen[27])
	histMC.SetBinError(27,other_cen[29])		
	histMC.SetBinError(28,other_fwd[25])
	histMC.SetBinError(29,other_fwd[27])
	histMC.SetBinError(30,other_fwd[29])
	

  	
	names = ["inclusive c", "b-Veto c", "b-Tagged c", "inclusive f", "b-Veto f", "b-Tagged f",
                "inclusive c", "b-Veto c", "b-Tagged c", "inclusive f", "b-Veto f", "b-Tagged f",
	        "inclusive c", "b-Veto c", "b-Tagged c", "inclusive f", "b-Veto f", "b-Tagged f",
          	"inclusive c", "b-Veto c", "b-Tagged c", "inclusive f", "b-Veto f", "b-Tagged f",
	        "inclusive c", "b-Veto c", "b-Tagged c", "inclusive f", "b-Veto f", "b-Tagged f"]

	
	for index, name in enumerate(names):
		histObs.GetXaxis().SetBinLabel(index+1,name)
	

        for i in range(1,histMC.GetNbinsX()+1):
            pred = histPred.GetBinContent(i)
            pred_e = histPred.GetBinError(i)
            dy = histDY.GetBinContent(i)
            dy_e = histDY.GetBinError(i)
            mc = histMC.GetBinContent(i)
            mc_e = histMC.GetBinError(i)
 
            histFlavSym.SetBinContent(i, pred-dy-mc)
            histFlavSym.SetBinError(i, math.sqrt(pred_e*pred_e-dy_e*dy_e-mc_e*mc_e))

	
	errGraph = r.TGraphAsymmErrors()
	for i in range(1,histFlavSym.GetNbinsX()+1):
		errGraph.SetPoint(i,i-0.5,histPred.GetBinContent(i))
		errGraph.SetPointError(i, 0.5, 0.5,histPred.GetBinError(i), histPred.GetBinError(i))
		

	errGraph.SetFillColor(myColors["MyBlue"])
	errGraph.SetFillStyle(3001)	

	histFlavSym.SetLineColor(r.kBlue+3)
	histFlavSym.SetLineWidth(2)
	
	histDY.SetLineColor(r.kGreen+3)
	histDY.SetFillColor(r.kGreen+3)
	histDY.SetFillStyle(3002)
	
	histMC.SetLineColor(r.kViolet+3)
	histMC.SetFillColor(r.kViolet+3)
	histMC.SetFillStyle(3002)

	
	
	#~ histDY.Add(histMC)
	stack = THStack()
	stack.Add(histMC)	
	stack.Add(histDY)	
	stack.Add(histFlavSym)

	
	
	histObs.GetYaxis().SetRangeUser(0,900)
	histObs.GetYaxis().SetTitle("Events")
	histObs.LabelsOption("v")

	histObs.UseCurrentStyle()
	histObs.Draw("pe")

	
	
	#~ hCanvas.DrawFrame(-0.5,0,30.5,65,"; %s ; %s" %("","Events"))
	
	latex = TLatex()
	latex.SetTextFont(42)
	latex.SetTextAlign(31)
	latex.SetTextSize(0.04)
	latex.SetNDC(True)
	latexCMS = TLatex()
	latexCMS.SetTextFont(61)
	#latexCMS.SetTextAlign(31)
	latexCMS.SetTextSize(0.06)
	latexCMS.SetNDC(True)
	latexCMSExtra = TLatex()
	latexCMSExtra.SetTextFont(52)
	#latexCMSExtra.SetTextAlign(31)
	latexCMSExtra.SetTextSize(0.045)
	latexCMSExtra.SetNDC(True)		
        print "Holaaaaaaaaaaaaaaa" 
	


	intlumi = TLatex()
	intlumi.SetTextAlign(12)
	intlumi.SetTextSize(0.03)
	intlumi.SetNDC(True)		

	latex.DrawLatex(0.95, 0.96, "%s fb^{-1} (13 TeV)"%"2.3")
	
	cmsExtra = "Preliminary"
	latexCMS.DrawLatex(0.19,0.88,"CMS")
	if "Simulation" in cmsExtra:
		yLabelPos = 0.81	
	else:
		yLabelPos = 0.84	

	latexCMSExtra.DrawLatex(0.19,yLabelPos,"%s"%(cmsExtra))

	
	leg1 = r.TLegend(0.4, 0.85, 0.925, 0.95,"","brNDC")
	leg1.SetNColumns(2)
	leg1.SetFillColor(10)
	leg1.SetLineColor(10)
	leg1.SetShadowColor(0)
	leg1.SetBorderSize(1)
	
	
	leg1.AddEntry(histObs,"Data    ","pe")
	leg1.AddEntry(histPred, "Total backgrounds","l")
	
	
	leg2 = r.TLegend(0.41, 0.75, 0.925, 0.85,"","brNDC")
	leg2.SetNColumns(3)
	leg2.SetFillColor(10)
	leg2.SetLineColor(10)
	leg2.SetShadowColor(0)
	leg2.SetBorderSize(1)
	
	leg2.AddEntry(errGraph,"Total uncert. ", "f")	
	leg2.AddEntry(histDY,"Z-jets  ", "f")
	leg2.AddEntry(histMC,"Rare", "f")
	
	

	stack.Draw("samehist")
	errGraph.Draw("same02")
	
	histObs.Draw("pesame")
	
	leg1.Draw("same")
	leg2.Draw("same")

	

	line1 = r.TLine(6,0,6,350)
	line2 = r.TLine(12,0,12,350)
	line3 = r.TLine(18,0,18,350)
	line4 = r.TLine(24,0,24,350)

	line1.SetLineColor(r.kBlack)
	line2.SetLineColor(r.kBlack)
	line3.SetLineColor(r.kBlack)
	line4.SetLineColor(r.kBlack)

	line1.SetLineWidth(2)
	line2.SetLineWidth(2)
	line3.SetLineWidth(4)
	line4.SetLineWidth(2)

	line1.Draw("same")
	line2.Draw("same")
	line3.Draw("same")
	line4.Draw("same")
	
	line5 = r.TLine(3,0,3,320)
	line6 = r.TLine(9,0,9,320)
	line7 = r.TLine(15,0,15,320)
	line8 = r.TLine(21,0,21,320)
	line9 = r.TLine(27,0,27,320)

	line5.SetLineColor(r.kBlack)
	line6.SetLineColor(r.kBlack)
	line7.SetLineColor(r.kBlack)
	line8.SetLineColor(r.kBlack)
	line9.SetLineColor(r.kBlack)

	line5.SetLineWidth(2)
	line6.SetLineWidth(2)
	line7.SetLineWidth(2)
	line8.SetLineWidth(2)
	line9.SetLineWidth(2)
	line5.SetLineStyle(r.kDashed)
	line6.SetLineStyle(r.kDashed)
	line7.SetLineStyle(r.kDashed)
	line8.SetLineStyle(r.kDashed)
	line9.SetLineStyle(r.kDashed)

	line5.Draw("same")
	line6.Draw("same")
	line7.Draw("same")
	line8.Draw("same")
	line9.Draw("same")


	label = r.TLatex()
	label.SetTextAlign(12)
	label.SetTextSize(0.04)
	label.SetTextColor(r.kBlack)	
	label.SetTextAngle(45)	
	
	label.DrawLatex(1,380,"low-mass")
	label.DrawLatex(8,380,"on-Z")
	label.DrawLatex(13,380,"high-mass")
	label.DrawLatex(19,380,"below-Z")
	label.DrawLatex(25,380,"above-Z")


	plotPad.RedrawAxis()
	
	hCanvas.Print("edgeOverview.pdf")
	hCanvas.Print("edgeOverview.root")

	


	
def main():
	
        makeOverviewPlotEdge()	
	#makeOverviewPlotWithOnZ()	


main()

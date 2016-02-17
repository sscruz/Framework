import os
import sys
import ROOT as r

from include.setTDRStyle import setTDRStyle
from include.helper import createMyColors
from include.helper import myColors
from ROOT import TCanvas, TH1F, TPad, THStack, TGraphErrors, TLatex, TLine

	
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

	
def main():
	
	
	makeOverviewPlotWithOnZ()	


main()

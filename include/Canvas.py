from ROOT import TCanvas, TLegend, TPad, TLine, TLatex, TH1F, THStack, TGraphErrors, TLine, TPaveStats, TGraph, TArrow
import ROOT as r
import os

class Canvas:
   'Common base class for all Samples'

   def __init__(self, name, format, x1, y1, x2, y2):
      self.name = name
      self.format = format
      self.plotNames = [name + "." + i for i in format.split(',')]
      self.myCanvas = TCanvas(name, name)
      self.ToDraw = []
      self.orderForLegend = []
      self.histos = []
      self.lines = []
      self.arrows= []
      self.latexs= []
      self.bands = []
      self.options = []
      self.labels = []      
      self.labelsOption = []
      self.myLegend = TLegend(x1, y1, x2, y2)
      self.myLegend.SetFillColor(r.kWhite)
      self.myLegend.SetTextFont(42)
      self.myLegend.SetTextSize(0.04)
      self.myLegend.SetLineWidth(0)
      self.myLegend.SetBorderSize(0)


   def banner(self, isData, lumi):
    
      latex = TLatex()
      latex.SetNDC();
      latex.SetTextAngle(0);
      latex.SetTextColor(r.kBlack);
      latex.SetTextFont(42);
      latex.SetTextAlign(31);
      latex.SetTextSize(0.05);
      latex.DrawLatex(0.21, 0.93, "CMS")

      latexb = TLatex()
      latexb.SetNDC();
      latexb.SetTextAngle(0);
      latexb.SetTextColor(r.kBlack);
      latexb.SetTextFont(42);
      latexb.SetTextAlign(31);
      latexb.SetTextSize(0.03);
 
      if(isData):
        latexb.DrawLatex(0.33, 0.93, "Preliminary")
      else:
        latexb.DrawLatex(0.33, 0.93, "Simulation")

      text_lumi = str(lumi) + " fb^{-1} (13 TeV)"
      latexc = TLatex()
      latexc.SetNDC();
      latexc.SetTextAngle(0);
      latexc.SetTextColor(r.kBlack);
      latexc.SetTextFont(42);
      latexc.SetTextAlign(31);
      latexc.SetTextSize(0.05);
      latexc.DrawLatex(0.90, 0.93, text_lumi)

   def banner2(self, isData, lumi):
    
      latex = TLatex()
      latex.SetNDC();
      latex.SetTextAngle(0);
      latex.SetTextColor(r.kBlack);
      latex.SetTextFont(42);
      latex.SetTextAlign(31);
      latex.SetTextSize(0.05);
      latex.DrawLatex(0.23, 0.93, "CMS")

      latexb = TLatex()
      latexb.SetNDC();
      latexb.SetTextAngle(0);
      latexb.SetTextColor(r.kBlack);
      latexb.SetTextFont(42);
      latexb.SetTextAlign(31);
      latexb.SetTextSize(0.03);
 
      if(isData):
        latexb.DrawLatex(0.37, 0.93, "Preliminary")
      else:
        latexb.DrawLatex(0.37, 0.93, "Simulation")

      text_lumi = str(lumi) + " fb^{-1} (13 TeV)"
      latexc = TLatex()
      latexc.SetNDC();
      latexc.SetTextAngle(0);
      latexc.SetTextColor(r.kBlack);
      latexc.SetTextFont(42);
      latexc.SetTextAlign(31);
      latexc.SetTextSize(0.05);
      latexc.DrawLatex(0.90, 0.93, text_lumi)

   def addBand(self, x1, y1, x2, y2, color, opacity):

      grshade = TGraph(4)
      grshade.SetPoint(0,x1,y1)
      grshade.SetPoint(1,x2,y1)
      grshade.SetPoint(2,x2,y2)
      grshade.SetPoint(3,x1,y2)
      #grshade.SetFillStyle(3001)
      grshade.SetFillColorAlpha(color, opacity)
      self.bands.append(grshade)

   def addLine(self, x1, y1, x2, y2, color, thickness = 0.):
      line = TLine(x1,y1,x2,y2)
      line.SetLineColor(color)
      if thickness:
          line.SetLineWidth(thickness)
      self.lines.append(line)

   def addArrow(self, x1, y1, x2, y2, color, option, thickness = 0.):
      arrow = TArrow(x1,y1,x2,y2, 0.05, option)
      arrow.SetLineColor(color)
      if thickness:
          arrow.SetLineWidth(thickness)
      self.arrows.append(arrow)

   def addLatex(self, x1, y1, text, font=42, size = 0.04):
      lat = [x1, y1, text, font, size]
      self.latexs.append(lat)

 
   def addHisto(self, h, option, label, labelOption, color, ToDraw, orderForLegend):

      if(color != ""):
          h.SetLineColor(color)
          h.SetMarkerColor(color)
      if(label == ""):
          label = h.GetTitle()

      self.histos.append(h)
      self.options.append(option)
      self.labels.append(label)
      self.labelsOption.append(labelOption)
      self.ToDraw.append(ToDraw)
      self.orderForLegend.append(orderForLegend)

   def addGraph(self, h, option, label, labelOption, color, ToDraw, orderForLegend):

      if(color != ""):
          h.SetLineColor(color)
          h.SetMarkerColor(color)
      if(label == ""):
          label = h.GetTitle()

      self.histos.append(h)
      self.options.append(option)
      self.labels.append(label)
      self.labelsOption.append(labelOption)
      self.ToDraw.append(ToDraw)
      self.orderForLegend.append(orderForLegend)


   def addStack(self, h, option, ToDraw, orderForLegend):

      legendCounter = orderForLegend
      if(orderForLegend < len(self.orderForLegend)):
          legendCounter = len(self.orderForLegend)

      self.addHisto(h, option, "", "", "", ToDraw, -1)  
      for h_c in h.GetHists():
          self.addHisto(h_c, "H", h_c.GetTitle(), "F", "", 0, legendCounter)
          legendCounter = legendCounter + 1
       

 
   def makeLegend(self):

      for i in range(0, len(self.histos)):
          for j in range(0, len(self.orderForLegend)):
              if(self.orderForLegend[j] != -1 and self.orderForLegend[j] == i):
                  self.myLegend.AddEntry(self.histos[j], self.labels[j], self.labelsOption[j])
          

   def ensurePath(self, _path):
      d = os.path.dirname(_path)
      if not os.path.exists(d):
         os.makedirs(d)

   def saveRatio(self, legend, isData, log, lumi, hdata, hMC, r_ymin=0, r_ymax=2):

      self.myCanvas.cd()

      pad1 = TPad("pad1", "pad1", 0, 0.2, 1, 1.0)
      pad1.SetBottomMargin(0.1)
      pad1.Draw()
      pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.2)
      pad2.SetTopMargin(0.1);
      pad2.SetBottomMargin(0.3);
      pad2.Draw();

      pad1.cd()
      if(log):
          pad1.SetLogy(1)

      for i in range(0, len(self.histos)):
          if(self.ToDraw[i] != 0):
              self.histos[i].Draw(self.options[i])

      if(legend):
          self.makeLegend()
          self.myLegend.Draw()

      for band in self.bands:
          band.Draw('f')
  
      for line in self.lines:
          line.Draw()
  
      for arrow in self.arrows:
          arrow.Draw()
  
      for latex in self.latexs:
          lat = TLatex()
          lat.SetNDC()
          lat.SetTextSize(latex[-1])
          lat.SetTextFont(latex[-2])
          lat.DrawLatex(latex[0], latex[1], latex[2])
  
      
      ratio = hdata.Clone("ratio")
      ratio.Divide(hMC)

      ratio.SetTitle("")
      ratio.GetYaxis().SetRangeUser(r_ymin, r_ymax);
      ratio.GetYaxis().SetTitle("Ratio");
      ratio.GetYaxis().CenterTitle();
      ratio.GetYaxis().SetLabelSize(0.12);
      ratio.GetXaxis().SetLabelSize(0.12);
      ratio.GetYaxis().SetTitleOffset(0.3);
      ratio.GetYaxis().SetNdivisions(4);
      ratio.GetYaxis().SetTitleSize(0.14);
      ratio.GetXaxis().SetTitleSize(0.14);
      ratio.GetXaxis().SetTitle('');
      ratio.SetMarkerStyle(20);
      ratio.SetMarkerSize(0.8*ratio.GetMarkerSize());
      ratio.SetMarkerColor(r.kGray+3);
      ratio.SetLineColor(r.kGray+3);


      pad2.cd();  
      line = TLine(ratio.GetBinLowEdge(1), 1, ratio.GetBinLowEdge(ratio.GetNbinsX()+1), 1)
      line.SetLineColor(r.kRed);
      ratio.Draw();
      line.Draw();
      ratio.Draw("E,SAME");

      pad1.cd()
      self.banner(isData, lumi)
      for plotName in self.plotNames:
          path = 'plots/'+plotName
          self.ensurePath(path)
          self.myCanvas.SaveAs(path)


      #del self.myCanvas


   def save(self, legend, isData, log, lumi):

      self.myCanvas.cd()
      
      if(log):
          self.myCanvas.GetPad(0).SetLogy(1)
     
      for i in range(0, len(self.histos)):
          if(self.ToDraw[i] != 0):        
              self.histos[i].Draw(self.options[i])

      for band in self.bands:
          band.Draw('f')
  
      for line in self.lines:
          line.Draw()
  
      for arrow in self.arrows:
          arrow.Draw()
  
      for latex in self.latexs:
          lat = TLatex()
          lat.SetNDC()
          lat.SetTextSize(latex[-1])
          lat.SetTextFont(latex[-2])
          lat.DrawLatex(latex[0], latex[1], latex[2])
  
      ## ps = self.histos[0].GetListOfFunctions().FindObject('stat')
      ## if ps:
      ##   ps.SetX1NDC(0.15)
      ##   ps.SetX2NDC(0.55)

      ##   ps.SetY1NDC(0.15)
      ##   ps.SetY2NDC(0.25)
            

      if(legend):
          self.makeLegend()
          self.myLegend.Draw()

      self.banner2(isData, lumi)
      for plotName in self.plotNames:
          path = 'plots/'+plotName
          self.ensurePath(path)
          self.myCanvas.SaveAs(path)

      #del self.myCanvas




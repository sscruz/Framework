from ROOT import TCanvas, TLegend, TPad, TLine, TLatex, TH1F, THStack, TGraphErrors, TLine, TPaveStats, TGraph, TArrow
import ROOT as r
import os, copy

class Canvas:
   'Common base class for all Samples'

   def __init__(self, name, _format, x1, y1, x2, y2, ww=0, hh=0):
      self.name = name
      self.format = _format
      self.plotNames = [name + "." + i for i in _format.split(',')]
      self.myCanvas = TCanvas(name, name) if not ww else TCanvas(name, name, ww, hh)
      self.ToDraw = []
      self.orderForLegend = []
      self.histos = []
      self.lines = []
      self.arrows= []
      self.latexs= []
      self.bands = []
      self.ratiobands = []
      self.options = []
      self.labels = []      
      self.labelsOption = []
      self.myLegend = TLegend(x1, y1, x2, y2)
      self.myLegend.SetFillColor(r.kWhite)
      self.myLegend.SetTextFont(42)
      self.myLegend.SetTextSize(0.04)
      self.myLegend.SetLineWidth(0)
      self.myLegend.SetBorderSize(0)
      r.gStyle.SetPadRightMargin(0.05)


   def banner(self, isData, lumi):





      latex = TLatex()

      latex.SetNDC()

      latex.SetTextColor(r.kBlack)

      latex.SetTextFont(61)

      latex.SetLineWidth(2)

      latex.SetTextAlign(31)

      latex.SetTextSize(0.04)

      latex.DrawLatex(0.21, 0.93, "CMS")



      latexb = TLatex()

      latexb.SetNDC();

      latexb.SetTextColor(r.kBlack);

      latexb.SetLineWidth(2)

      latexb.SetTextFont(52);

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

      latexc.SetTextSize(0.04);

      latexc.DrawLatex(0.90, 0.93, text_lumi)

   def banner2(self, isData, lumi):





      latex = TLatex()

      latex.SetNDC()

      latex.SetTextColor(r.kBlack)

      latex.SetTextFont(61)

      latex.SetLineWidth(2)

      latex.SetTextAlign(31)

      latex.SetTextSize(0.04)

      latex.DrawLatex(0.21, 0.93, "CMS")



      latexb = TLatex()

      latexb.SetNDC();

      latexb.SetTextColor(r.kBlack);

      latexb.SetLineWidth(2)

      latexb.SetTextFont(52);

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

      latexc.SetTextSize(0.04);

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

   def addBandRatio(self, x1, y1, x2, y2, color, opacity):
      grshade = TGraph(4)
      grshade.SetPoint(0,x1,y1)
      grshade.SetPoint(1,x2,y1)
      grshade.SetPoint(2,x2,y2)
      grshade.SetPoint(3,x1,y2)
      #grshade.SetFillStyle(3001)
      grshade.SetFillColorAlpha(color, opacity)
      self.ratiobands.append(grshade)

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

   def makeOFHisto(self, histo):
      nbinsx = histo.GetNbinsX()
      xmin = histo.GetXaxis().GetXmin(); xmax = histo.GetXaxis().GetXmax();
      newhisto = r.TH1F(histo.GetName() +'_withOFBin', histo.GetTitle()+'_withOFBin', nbinsx+1, xmin, xmax+(xmax-xmin)/nbinsx)
      newhisto.Sumw2()
      newhisto.SetMarkerColor(histo.GetMarkerColor())
      newhisto.SetMarkerStyle(histo.GetMarkerStyle())
      newhisto.SetMarkerSize (histo.GetMarkerSize ())
      newhisto.SetLineColor(histo.GetLineColor())
      newhisto.SetLineStyle(histo.GetLineStyle())
      newhisto.SetMaximum(histo.GetMaximum())
      for i in range(1,nbinsx+2):
         newhisto.SetBinContent(i,histo.GetBinContent(i))
         newhisto.SetBinError  (i,histo.GetBinError  (i))
      return newhisto
        
 
   def addHisto(self, h, option, label, labelOption, color, ToDraw, orderForLegend, doOF = False):

      if(color != ""):
          h.SetLineColor(color)
          h.SetMarkerColor(color)
      if(label == ""):
          label = h.GetTitle()

      self.histos.append(h if not doOF else self.makeOFHisto(h))
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


   def addStack(self, h, option, ToDraw, orderForLegend, labels = []):

      legendCounter = orderForLegend
      if(orderForLegend < len(self.orderForLegend)):
          legendCounter = len(self.orderForLegend)

      self.addHisto(h, option, "", "", "", ToDraw, -1)  
      for ind, h_c in enumerate(h.GetHists()):
          self.addHisto(h_c, "H", "%s (%4.1f)"%(h_c.GetTitle(),h_c.Integral()) if len(labels) == 0 else labels[ind], "F", "", 0, legendCounter)
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

      pad1 = TPad("pad1", "pad1", 0, 0.25, 1, 1.0)
      pad1.SetBottomMargin(0.1)
      pad1.Draw()
      pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.25)
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
  
      
      if type(hMC) != list:
          hMClist = [hMC]
      else: hMClist = hMC

      ratios = []

      for tmp_hMC in hMClist:
          #print 'making ratio for', tmp_hMC.GetName()
          ind = hMClist.index(tmp_hMC)
          tmp_ratio = hdata.Clone(tmp_hMC.GetName()+'_ratio')
          tmp_ratio.Divide(tmp_hMC)
          ## print 'at histogram', hdata.GetName()
          ## print 'ranges', hdata.GetXaxis().GetXmin(), hdata.GetXaxis().GetXmax(), hdata.GetNbinsX()
          ## print 'dividing by ', tmp_hMC.GetName()
          ## print 'ranges', tmp_hMC.GetXaxis().GetXmin(), tmp_hMC.GetXaxis().GetXmax(), tmp_hMC.GetNbinsX()

          tmp_ratio.SetTitle("")
          tmp_ratio.GetYaxis().SetRangeUser(r_ymin, r_ymax);
          tmp_ratio.GetYaxis().SetTitle("Data/Pred");
          tmp_ratio.GetYaxis().CenterTitle();
          tmp_ratio.GetYaxis().SetLabelSize(0.12);
          tmp_ratio.GetXaxis().SetLabelSize(0.12);
          tmp_ratio.GetYaxis().SetTitleOffset(0.3);
          tmp_ratio.GetYaxis().SetNdivisions(4);
          tmp_ratio.GetYaxis().SetTitleSize(0.14);
          tmp_ratio.GetXaxis().SetTitleSize(0.14);
          tmp_ratio.GetXaxis().SetTitle('');
          tmp_ratio.SetMarkerStyle(tmp_hMC.GetMarkerStyle());
          #tmp_ratio.SetMarkerSize (tmp_hMC.GetMarkerSize());
          tmp_ratio.SetMarkerSize(0)
          tmp_ratio.SetFillColor(r.kBlue-5)
          tmp_ratio.SetMarkerColor(tmp_hMC.GetMarkerColor());
          #tmp_ratio.SetMarkerSize(0.6*tmp_ratio.GetMarkerSize());
          #tmp_ratio.SetMarkerColor(r.kBlack if len(hMClist) == 1 else tmp_hMC.GetMarkerColor());
          #tmp_ratio.SetLineColor  (r.kBlack if len(hMClist) == 1 else tmp_hMC.GetLineColor  ());
          tmp_ratio.SetLineColor  (tmp_hMC.GetLineColor());
          tmp_ratio.SetLineStyle(tmp_hMC.GetLineStyle())
          ratios.append(tmp_ratio)
          xmin = tmp_ratio.GetBinLowEdge(1)
          xmax = tmp_ratio.GetBinLowEdge(tmp_ratio.GetNbinsX()+1)

      #tmp_ratio.Draw("E,SAME");
      pad2.cd();  
      for rat in ratios:
          rat.Draw('P,E2,same');
          points = rat.Clone()
          points.SetMarkerStyle(r.kFullCircle)
          points.Draw("P,same")

      for band in self.ratiobands:
          band.Draw('f')

      line = TLine(xmin, 1, xmax, 1)
      line.SetLineColor(r.kGray+2);
      line.Draw('');

      pad1.cd()
      self.banner(isData, lumi)
      for plotName in self.plotNames:
          path = 'plots/'+plotName
          self.ensurePath(path)
          self.myCanvas.SaveAs(path)

      pad1.IsA().Destructor(pad1) 
      pad2.IsA().Destructor(pad2) 
      self.myLegend.IsA().Destructor(self.myLegend)
      self.myCanvas.IsA().Destructor(self.myCanvas)

 
   def save(self, legend, isData, log, lumi, ymin=0, ymax=0):

      self.myCanvas.cd()
      
      if(log):
          self.myCanvas.GetPad(0).SetLogy(1)
     
      for i in range(0, len(self.histos)):
          if(self.ToDraw[i] != 0):        
              if ymin and ymax:
                  self.histos[i].GetYaxis().SetRangeUser(ymin, ymax)
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

      self.myLegend.IsA().Destructor(self.myLegend)
      self.myCanvas.IsA().Destructor(self.myCanvas)

 #del self.myCanvas




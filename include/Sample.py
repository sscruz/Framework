import ROOT as r
from array import array
from ROOT import TTree, TFile, TCut, TH1F, TH2F, THStack, TCanvas


class Sample:
   'Common base class for all Samples'

   def __init__(self, name, location, friendlocation, xsection, isdata):
      self.name = name
      self.location = location
      self.xSection = xsection
      self.isData = isdata
      self.tfile = TFile(self.location+self.name+'/treeProducerSusyEdge/tree.root')
      self.ftfile = TFile(self.location+'/'+friendlocation+'/evVarFriend_'+self.name+'.root')
      self.ttree = self.tfile.Get('tree')
      self.ttree.AddFriend('sf/t',self.ftfile)
      if not self.isData:
        gw = 0.
        for i in self.ttree:
            gw = abs(i.genWeight)
            if gw: break
        self.count = self.tfile.Get('SumGenWeights').GetBinContent(1)/abs(gw)
      else:
        self.count = self.tfile.Get('Count').GetEntries()
      self.lumWeight = 1.0
      if(self.isData == 0):
        self.lumWeight = self.xSection / self.count
      self.puWeight = "1.0"
      if(self.isData == 0):
        self.puWeight = "t.PileupW_Edge"
#        self.puWeight = "1.0"
   def printSample(self):
      print "#################################"
      print "Sample Name: ", self.name
      print "Sample Location: ", self.location
      print "Sample XSection: ", self.xSection
      print "Sample IsData: ", self.isData
      print "Sample LumWeight: ", self.lumWeight
      print "#################################"


   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel):
 
      if(xmin == xmax):
        h = TH1F(name, "", len(nbin)-1, array('d', nbin))
        ylabel = "# events"
      else:
        h = TH1F(name, "", nbin, xmin, xmax)
        bw = int((xmax-xmin)/nbin)
        ylabel = "Events / " + str(bw) + " GeV"
      h.Sumw2()
      h.GetXaxis().SetTitle(xlabel)
      h.GetYaxis().SetTitle(ylabel)

      addCut = ""
      if self.isData:
        if(name.find("DoubleMuon") != -1):
          addCut = "(!((Lep1_pdgId_Edge * Lep2_pdgId_Edge == -121) || (Lep1_pdgId_Edge * Lep2_pdgId_Edge == -143)))"
          cut = cut + "* ( " + addCut + " )"
        if(name.find("DoubleEG") != -1):
          addCut = "(!((Lep1_pdgId_Edge * Lep2_pdgId_Edge == -169) || (Lep1_pdgId_Edge * Lep2_pdgId_Edge == -143)))"
          cut = cut + "* ( " + addCut + " )"
        if(name.find("MuonEG") != -1):
          addCut = "(!((Lep1_pdgId_Edge * Lep2_pdgId_Edge == -121) || (Lep1_pdgId_Edge * Lep2_pdgId_Edge == -169)))"
          cut = cut + "* ( " + addCut + " )"
           
      if(self.isData == 0):
        cut = cut + "* ( " + str(self.lumWeight*lumi) + " * genWeight/abs(genWeight) * " + self.puWeight + " )" 
      
      self.ttree.Project(name, var, cut, options) 
      return h

   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel):
   
     h = TH2F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax)
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     if(self.isData == 0):
        cut = cut + "* ( " + str(self.lumWeight*lumi) + " * genWeight/abs(genWeight) * " + self.puWeight + " )" 
     
     self.ttree.Project(name, var, cut, options) 
     return h

class Block:
   'Common base class for all Sample Blocks'

   def __init__(self, name, label, color, isdata):
      self.name  = name
      self.color = color
      self.isData = isdata
      self.label = label
      self.samples = []

   def printBlock(self):

      print "####################"
      print "Block Name: ", self.name
      print "Block Color: ", self.color
      print "Block IsData: ", self.isData
      print "####################"
      print "This block contains the following Samples"

      for l in self.samples:
        l.printSample()
     

   def addSample(self, s):
      self.samples.append(s)

   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel):

     if(xmin == xmax):
       h = TH1F(name, "", len(nbin)-1, array('d', nbin))
       ylabel = "# events"
     else:
       h = TH1F(name, "", nbin, xmin, xmax)
       bw = int((xmax-xmin)/nbin)
       ylabel = "Events / " + str(bw) + " GeV"
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)

     for s in self.samples:
       AuxName = "auxT1_sample" + s.name
       haux = s.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel)
       h.Add(haux)
       del haux


     h.SetLineColor(self.color)
     h.SetMarkerColor(self.color)
     h.SetTitle(self.label)

     return h

   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel):
   
     h = TH2F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax)
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     for s in self.samples:
     
       AuxName = "auxT2_block" + s.name
       haux = s.getTH2F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel)
       h.Add(haux)
       del haux

     return h   

       

class Tree:
   'Common base class for a physics meaningful tree'

   def __init__(self, fileName, name, isdata):
      self.name  = name
      self.isData = isdata
      self.blocks = []
      self.parseFileName(fileName)

   def parseFileName(self, fileName):

      f = open(fileName)

      for l in f.readlines():

        if (l[0] == "#" or len(l) < 2):
          continue

        splitedLine = str.split(l)
        block       = splitedLine[0]
        theColor    = splitedLine[1]
        name        = splitedLine[2]
        label       = splitedLine[3]
        location    = splitedLine[4]
        flocation   = splitedLine[5]
        xsection    = float(splitedLine[6])
        isdata      = int(splitedLine[7])

        color = 0
        plusposition = theColor.find("+")
        if(plusposition == -1):
          color = eval(theColor)
        else:
          color = eval(theColor[0:plusposition])
          color = color + int(theColor[plusposition+1:len(theColor)])

        sample = Sample(name, location, flocation, xsection, isdata)
        coincidentBlock = [l for l in self.blocks if l.name == block]

        if(coincidentBlock == []):

          newBlock = Block(block, label, color, isdata)
          newBlock.addSample(sample)
          self.addBlock(newBlock)

        else:

          coincidentBlock[0].addSample(sample)





   def printTree(self):

      print "######"
      print "Tree Name: ", self.name
      print "Tree IsData: ", self.isData
      print "######"
      print "This Tree contains the following Blocks"

      for l in self.blocks:
        l.printBlock()
     

   def addBlock(self, b):
      self.blocks.append(b)



   def getYields(self, lumi, var, xmin, xmax, cut):
  
      h = self.getTH1F(lumi, "yields", var, 1, xmin, xmax, cut, "", "")
      nbinmin = h.FindBin(xmin)
      nbinmax = h.FindBin(xmax)
      error = r.Double()
      value = h.IntegralAndError(nbinmin, nbinmax, error)
      y = [value, error]
      
      del h
      return y

   def getStack(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel):
   
     hs = THStack(name, "")
     for b in self.blocks:
     
       AuxName = "auxStack_block_" + name + "_" + b.name
       haux = b.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel)
       haux.SetFillColor(b.color)
       hs.Add(haux)
       del haux


     can_aux = TCanvas("can_aux")
     can_aux.cd()
     hs.Draw()
     del can_aux

     hs.GetXaxis().SetTitle(xlabel)
     b = int((xmax-xmin)/nbin)
     ylabel = "Events / " + str(b) + " GeV"
     hs.GetYaxis().SetTitle(ylabel)

     return hs   


   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel):
   
     if(xmin == xmax):
       h = TH1F(name, "", len(nbin)-1, array('d', nbin))
       ylabel = "# events"
     else:
       h = TH1F(name, "", nbin, xmin, xmax)
       bw = int((xmax-xmin)/nbin)
       ylabel = "Events / " + str(bw) + " GeV"
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     for b in self.blocks:
       AuxName = "auxh1_block_" + name + "_" + b.name
       haux = b.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel)
       h.Add(haux)
       del haux

     return h   

   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel):
   
     h = TH2F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax)
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     for b in self.blocks:
     
       AuxName = "aux_block" + name + "_" + b.name
       haux = b.getTH2F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel)
       h.Add(haux)
       del haux

     return h   











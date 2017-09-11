import ROOT as r
from array import array
from ROOT import TTree, TFile, TCut, TH1F, TH2F, TH3F, THStack, TCanvas
import include.LeptonSF
import include.FastSimSF

class Sample:
   'Common base class for all Samples'

   def __init__(self, name, friendlocation, xsection, isdata, doKfactor, isScan, isOnEOS):
      self.name = name
      self.location = friendlocation
      self.xSection = xsection
      self.doKfactor = doKfactor
      self.isData = isdata
      if isOnEOS:
          self.ftfile = TFile(friendlocation)
      else:
          ftfileloc = friendlocation+'/evVarFriend_'+self.name+'.root' 
          self.ftfile = TFile(ftfileloc)                                    
      self.ttree = self.ftfile.Get('sf/t')
      self.isScan = isScan
      if not self.isData and not self.isScan:
        gw = 0.
        for i in self.ttree:
            gw = abs(i.genWeight_Edge)
            if gw: break
        self.count = self.ftfile.Get('SumGenWeights').GetBinContent(1)/abs(gw)
      else:
        #self.count = self.ftfile.Get('Count').GetBinContent(1)
        self.count = self.ftfile.Get('sf/t').GetEntries()
      self.puWeight   = '1.0'
      self.SFWeight   = '1.0'
      self.btagWeight = '1.0'
      self.triggWeight = '1.0'
      self.ISRWeight  = '1.0'

      #print self.name
      #print isdata
      if not self.isData and not self.isScan  > 0:
        self.lumWeight = self.xSection / self.count
        self.puWeight    = "PileupW_Edge"
        self.btagWeight  = "weight_btagsf_Edge"
        self.SFWeight = "LepSF(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSF(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge)"
        #print "self.name ", self.name
        print self.name
        print "weight ", self.lumWeight*35.9
        #self.triggWeight = "weight_trigger_Edge"

      if self.isScan > 0:
        self.lumWeight  =  1.0
        self.xSection = self.isScan
        self.puWeight    = "1.0"
        self.btagWeight  = "weight_btagsf_Edge"
        self.SFWeight = "LepSF(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSF(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge)*LepSFFastSim(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSFFastSim(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge)"
        self.ISRWeight = 'ISRweight_Edge'
        if self.isScan == 1:
            self.smsCount =   self.ftfile.Get('CountSMS')
        else:
            self.smsCount =  self.ftfile.Get('sf/t').GetEntries()
            self.lumWeight = self.xSection / self.smsCount
   def printSample(self):
      print "#################################"
      print "Sample Name: ", self.name
      print "Sample Location: ", self.location
      print "Sample XSection: ", self.xSection
      print "Sample IsData: ", self.isData
      print "Sample LumWeight: ", self.lumWeight
      print "#################################"


   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel,  extraWeight, doKfactorGENVar):
   #def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel, ofBin=True, extraWeight="1", ylabel = 'Events', doKfactorGENVar = "noKFactor"):
#      if doKfactorGENVar == 'doTGraph':
#          fileD =  TFile("kfactors.root");
#          hScale = fileD.Get("kfactorsPt");
#          hScale.Fit("pol6")
#          fit = hScale.GetFunction("pol6");
#          fileD.Close()    
      ofBin = True
      ylabel = 'Events'
      if(xmin == xmax):                                          
          _nbins = len(nbin)-1
          _arr = array('d', nbin)
          h = TH1F(name+'_noOF', "", _nbins, _arr)
          _newarr = _arr + array('d', [ 2*_arr[-1]-_arr[-2] ])
          h_of = TH1F(name, "", _nbins+1, _newarr)                                     
      else:
           h = TH1F(name+'_noOF', "", nbin, xmin, xmax)
           bw = int((xmax-xmin)/nbin)
           #if ylabel == "au":
           #    ylabel = "A.U"
           #else:
           ylabel = ylabel +"/ " + str(bw) + " GeV"
           h_of = TH1F(name, '', nbin+1, xmin, xmax+bw)
      h.Sumw2()
      h.GetXaxis().SetTitle(xlabel)
      h.GetYaxis().SetTitle(ylabel)

      h_of.Sumw2()
      h_of.GetXaxis().SetTitle(xlabel)
      h_of.GetYaxis().SetTitle(ylabel)

      addCut = ""
      #if self.isData:
      #  if(name.find("DoubleMuon") != -1):
      #    addCut = "(!((Lep1_pdgId_Edge * Lep2_pdgId_Edge == -121) || (Lep1_pdgId_Edge * Lep2_pdgId_Edge == -143)))"
      #    cut = cut + "* ( " + addCut + " ) "
      #  if(name.find("DoubleEG") != -1):
      #    addCut = "(!((Lep1_pdgId_Edge * Lep2_pdgId_Edge == -169) || (Lep1_pdgId_Edge * Lep2_pdgId_Edge == -143)))"
      #    cut = cut + "* ( " + addCut + " ) "
      #  if(name.find("MuonEG") != -1):
      #    addCut = "(!((Lep1_pdgId_Edge * Lep2_pdgId_Edge == -121) || (Lep1_pdgId_Edge * Lep2_pdgId_Edge == -169)))"
      #    cut = cut + "* ( " + addCut + " ) "
       
      #fileD =  TFile("kfactors.root");
      kf = "1"  
      if(self.isData == 0):
          if (self.doKfactor == 1 ): #this is the kfactor for ZZto4l
#              if doKfactorGENVar == 'doTGraph':
#                  kf = "("
#                  scalebins = [[0, 5], [5,10],[10, 15],[15, 20], [20,25], [25,30], [30,35], [35,40], [40,45], [45,50], [50,55], [55,60], [60,65], [65,70],[70, 75], [75,80], [80,85], [85,90], [90,95],[95, 100], [100, 1000]]
#                  for q in scalebins:   
#                      mini = float(min(scalebins[scalebins.index(q)]))
#                      maxi = float(max(scalebins[scalebins.index(q)]))
#                      mid = float((mini+maxi)/2)
#                      p0 = fit.GetParameter(0) 
#                      p1 = fit.GetParameter(1) 
#                      p2 = fit.GetParameter(2) 
#                      p3 = fit.GetParameter(3) 
#                      p4 = fit.GetParameter(4) 
#                      p5 = fit.GetParameter(5) 
#                      p6 = fit.GetParameter(6) 
#                      #x = hScale.GetY() 
#                      #kf = "( "+str(p0)+"+ ("+str(p1)+"*abs(GENptZZ_Edge)))"
#                      kf = "( "+str(p0)+"+ (("+str(p1)+")*GENptZZ_Edge) + (("+str(p2)+")*GENptZZ_Edge**2)+(("+str(p3)+")*GENptZZ_Edge**3)+(("+str(p4)+")*GENptZZ_Edge**4)+(("+str(p5)+")*GENptZZ_Edge**5)+(("+str(p6)+")*GENptZZ_Edge**6))*(GENptZZ_Edge < 100)"
#                      #kf = "( "+str(p0)+"+ (("+str(p1)+")*GENptZZ_Edge) + (("+str(p2)+")*GENptZZ_Edge**2)+(("+str(p3)+")*GENptZZ_Edge**3)+(("+str(p4)+")*GENptZZ_Edge**4)+(("+str(p5)+")*GENptZZ_Edge**5)+(("+str(p6)+")*GENptZZ_Edge**6))*(GENptZZ_Edge < 100) + 1.584*(GENptZZ_Edge > 100)"
#                      #kf = "( "+str(p0)+"+ ("+str(p1)+"*abs(GENptZZ_Edge)) + ("+str(p2)+"*abs(GENptZZ_Edge)^2)+("+str(p3)+"*abs(GENptZZ_Edge)^3)+("+str(p4)+"*abs(GENptZZ_Edge)^4)+("+str(p5)+"*abs(GENptZZ_Edge)^5)+("+str(p6)+"*abs(GENptZZ_Edge)^6))"
#                          #kf = kf +  str(x[scalebins.index(q)]) +"*(abs(GENptZZ_Edge)>" +str(mini)+"&&abs(GENptZZ_Edge)<="+ str(maxi) +")"
#                      #else:
#                      #    kf = kf + "+"+ str(x[scalebins.index(q)]) +"*(abs(GENptZZ_Edge)>" +str(mini)+" &&abs(GENptZZ_Edge)<="+ str(maxi) +")"
#                  #kf = kf+ ")"       
#                  print kf
              if doKfactorGENVar == 'ZZmass':
                  kf = "(1.23613*(abs(GENmassZZ_Edge)>0.0&&abs(GENmassZZ_Edge)<=25.0)+1.1755*(abs(GENmassZZ_Edge)>25.0 &&abs(GENmassZZ_Edge)<=50.0 )+ 1.1704*(abs(GENmassZZ_Edge)>50.0 &&abs(GENmassZZ_Edge)<=75.0 )+ 1.0314*(abs(GENmassZZ_Edge)>75.0 &&abs(GENmassZZ_Edge)<=100.0)+1.0528*(abs(GENmassZZ_Edge)>100.0&&abs(GENmassZZ_Edge)<=125.0)+1.1128*(abs(GENmassZZ_Edge)>125.0&&abs(GENmassZZ_Edge)<=150.0)+1.1336*(abs(GENmassZZ_Edge)>150.0&&abs(GENmassZZ_Edge)<=175.0)+1.1035*(abs(GENmassZZ_Edge)>175.0&&abs(GENmassZZ_Edge)<=200.0)+1.1005*(abs(GENmassZZ_Edge)>200.0&&abs(GENmassZZ_Edge)<=225.0)+1.1097*(abs(GENmassZZ_Edge)>225.0&&abs(GENmassZZ_Edge)<=250.0)+1.1206*(abs(GENmassZZ_Edge)>250.0&&abs(GENmassZZ_Edge)<=275.0)+1.1158*(abs(GENmassZZ_Edge)>275.0&&abs(GENmassZZ_Edge)<=300.0)+1.1390*(abs(GENmassZZ_Edge)>300.0&&abs(GENmassZZ_Edge)<=325.0)+1.1485*(abs(GENmassZZ_Edge)>325.0&&abs(GENmassZZ_Edge)<=350.0)+1.1461*(abs(GENmassZZ_Edge)>350.0&&abs(GENmassZZ_Edge)<=375.0)+1.1457*(abs(GENmassZZ_Edge)>375.0&&abs(GENmassZZ_Edge)<=400.0)+1.1382*(abs(GENmassZZ_Edge)>400.0&&abs(GENmassZZ_Edge)<=425.0)+1.1552*(abs(GENmassZZ_Edge)>425.0&&abs(GENmassZZ_Edge)<=450.0)+1.1367*(abs(GENmassZZ_Edge)>450.0&&abs(GENmassZZ_Edge)<=475.0)+1.1322*(abs(GENmassZZ_Edge)>475.0))"
              if doKfactorGENVar == 'ZZpt':
                  kf =  "(0.64155*(abs(GENptZZ_Edge)>0.0&&abs(GENptZZ_Edge)<=5.0)+1.0998*(abs(GENptZZ_Edge)>5.0 &&abs(GENptZZ_Edge)<=10.0) +1.2939*(abs(GENptZZ_Edge)>10.0&&abs(GENptZZ_Edge)<=15.0)+1.3785*(abs(GENptZZ_Edge)>15.0&&abs(GENptZZ_Edge)<=20.0)+1.4243*(abs(GENptZZ_Edge)>20.0&&abs(GENptZZ_Edge)<=25.0)+1.4503*(abs(GENptZZ_Edge)>25.0&&abs(GENptZZ_Edge)<=30.0)+1.4701*(abs(GENptZZ_Edge)>30.0&&abs(GENptZZ_Edge)<=35.0)+1.4882*(abs(GENptZZ_Edge)>35.0&&abs(GENptZZ_Edge)<=40.0)+1.5057*(abs(GENptZZ_Edge)>40.0&&abs(GENptZZ_Edge)<=45.0)+1.5021*(abs(GENptZZ_Edge)>45.0&&abs(GENptZZ_Edge)<=50.0)+1.5091*(abs(GENptZZ_Edge)>50.0&&abs(GENptZZ_Edge)<=55.0)+1.5246*(abs(GENptZZ_Edge)>55.0&&abs(GENptZZ_Edge)<=60.0)+1.5240*(abs(GENptZZ_Edge)>60.0&&abs(GENptZZ_Edge)<=65.0)+1.5241*(abs(GENptZZ_Edge)>65.0&&abs(GENptZZ_Edge)<=70.0)+1.5542*(abs(GENptZZ_Edge)>70.0&&abs(GENptZZ_Edge)<=75.0)+1.5254*(abs(GENptZZ_Edge)>75.0&&abs(GENptZZ_Edge)<=80.0)+1.5789*(abs(GENptZZ_Edge)>80.0&&abs(GENptZZ_Edge)<=85.0)+1.5303*(abs(GENptZZ_Edge)>85.0&&abs(GENptZZ_Edge)<=90.0)+1.5614*(abs(GENptZZ_Edge)>90.0&&abs(GENptZZ_Edge)<=95.0)+1.5446*(abs(GENptZZ_Edge)>95.0&&abs(GENptZZ_Edge)<=100.0)+1.5722*(abs(GENptZZ_Edge)>100.0))"
              if doKfactorGENVar == 'ZZdPhi':
                  kf = "(1.5158*(abs(GENphiZZ_Edge)>0.0&&abs(GENphiZZ_Edge)<=0.1)+1.4962*(abs(GENphiZZ_Edge)>0.1&&abs(GENphiZZ_Edge)<=0.2)+1.4955*(abs(GENphiZZ_Edge)>0.2&&abs(GENphiZZ_Edge)<=0.3)+1.4832*(abs(GENphiZZ_Edge)>0.3&&abs(GENphiZZ_Edge)<=0.4)+1.4655*(abs(GENphiZZ_Edge)>0.4&&abs(GENphiZZ_Edge)<=0.5)+1.4915*(abs(GENphiZZ_Edge)>0.5&&abs(GENphiZZ_Edge)<=0.6)+1.4411*(abs(GENphiZZ_Edge)>0.6&&abs(GENphiZZ_Edge)<=0.7)+1.4408*(abs(GENphiZZ_Edge)>0.7&&abs(GENphiZZ_Edge)<=0.8)+1.4143*(abs(GENphiZZ_Edge)>0.8&&abs(GENphiZZ_Edge)<=0.9)+1.4225*(abs(GENphiZZ_Edge)>0.9&&abs(GENphiZZ_Edge)<=1.0)+1.4010*(abs(GENphiZZ_Edge)>1.0&&abs(GENphiZZ_Edge)<=1.1)+1.4085*(abs(GENphiZZ_Edge)>1.1&&abs(GENphiZZ_Edge)<=1.2)+1.3812*(abs(GENphiZZ_Edge)>1.2&&abs(GENphiZZ_Edge)<=1.3)+1.3705*(abs(GENphiZZ_Edge)>1.3&&abs(GENphiZZ_Edge)<=1.4)+1.3473*(abs(GENphiZZ_Edge)>1.4&&abs(GENphiZZ_Edge)<=1.5)+1.3401*(abs(GENphiZZ_Edge)>1.5&&abs(GENphiZZ_Edge)<=1.6)+1.3126*(abs(GENphiZZ_Edge)>1.6&&abs(GENphiZZ_Edge)<=1.7)+1.2900*(abs(GENphiZZ_Edge)>1.7&&abs(GENphiZZ_Edge)<=1.8)+1.2553*(abs(GENphiZZ_Edge)>1.8&&abs(GENphiZZ_Edge)<=1.9)+1.2544*(abs(GENphiZZ_Edge)>1.9&&abs(GENphiZZ_Edge)<=2.0)+1.2241*(abs(GENphiZZ_Edge)>2.0&&abs(GENphiZZ_Edge)<=2.1)+1.1788*(abs(GENphiZZ_Edge)>2.1&&abs(GENphiZZ_Edge)<=2.2)+1.1626*(abs(GENphiZZ_Edge)>2.2&&abs(GENphiZZ_Edge)<=2.3)+1.1054*(abs(GENphiZZ_Edge)>2.3&&abs(GENphiZZ_Edge)<=2.4)+1.0747*(abs(GENphiZZ_Edge)>2.4&&abs(GENphiZZ_Edge)<=2.5)+1.0218*(abs(GENphiZZ_Edge)>2.5&&abs(GENphiZZ_Edge)<=2.6)+0.9463*(abs(GENphiZZ_Edge)>2.6&&abs(GENphiZZ_Edge)<=2.7)+0.8574*(abs(GENphiZZ_Edge)>2.7&&abs(GENphiZZ_Edge)<=2.8)+0.7166*(abs(GENphiZZ_Edge)>2.8&&abs(GENphiZZ_Edge)<=2.9)+1.1328*(abs(GENphiZZ_Edge)>2.9&&abs(GENphiZZ_Edge)<=3.1416))"
              if doKfactorGENVar == 'noKFactor':
                  kf = "1"
          if (self.doKfactor == 2): #this is the kfactor for ZZto2l2nu
              if doKfactorGENVar == 'ZZmass':
                  kf = "(1.25094*(abs(GENmassZZ_Edge)>0.0&&abs(GENmassZZ_Edge)<=25.0)+1.2245*(abs(GENmassZZ_Edge)>25.0 &&abs(GENmassZZ_Edge)<=50.0 )+1.1928*(abs(GENmassZZ_Edge)>50.0 &&abs(GENmassZZ_Edge)<=75.0 )+1.0459*(abs(GENmassZZ_Edge)>75.0 &&abs(GENmassZZ_Edge)<=100.0)+1.0832*(abs(GENmassZZ_Edge)>100.0&&abs(GENmassZZ_Edge)<=125.0)+1.0999*(abs(GENmassZZ_Edge)>125.0&&abs(GENmassZZ_Edge)<=150.0)+1.1669*(abs(GENmassZZ_Edge)>150.0&&abs(GENmassZZ_Edge)<=175.0)+1.1039*(abs(GENmassZZ_Edge)>175.0&&abs(GENmassZZ_Edge)<=200.0)+1.1059*(abs(GENmassZZ_Edge)>200.0&&abs(GENmassZZ_Edge)<=225.0)+1.1069*(abs(GENmassZZ_Edge)>225.0&&abs(GENmassZZ_Edge)<=250.0)+1.1119*(abs(GENmassZZ_Edge)>250.0&&abs(GENmassZZ_Edge)<=275.0)+1.1352*(abs(GENmassZZ_Edge)>275.0&&abs(GENmassZZ_Edge)<=300.0)+1.1189*(abs(GENmassZZ_Edge)>300.0&&abs(GENmassZZ_Edge)<=325.0)+1.1389*(abs(GENmassZZ_Edge)>325.0&&abs(GENmassZZ_Edge)<=350.0)+1.1546*(abs(GENmassZZ_Edge)>350.0&&abs(GENmassZZ_Edge)<=375.0)+1.1734*(abs(GENmassZZ_Edge)>375.0&&abs(GENmassZZ_Edge)<=400.0)+1.2009*(abs(GENmassZZ_Edge)>400.0&&abs(GENmassZZ_Edge)<=425.0)+1.1891*(abs(GENmassZZ_Edge)>425.0&&abs(GENmassZZ_Edge)<=450.0)+1.1854*(abs(GENmassZZ_Edge)>450.0&&abs(GENmassZZ_Edge)<=475.0)+1.12864*(abs(GENmassZZ_Edge)>475.0))"
              if doKfactorGENVar == 'ZZpt':
                  kf  = "(0.7436*(abs(GENptZZ_Edge)>0.0&&abs(GENptZZ_Edge)<=5.0)+1.14789*(abs(GENptZZ_Edge)>5.0&&abs(GENptZZ_Edge)<=10.0)+1.33815*(abs(GENptZZ_Edge)>10.0&&abs(GENptZZ_Edge)<=15.0)+1.41420*(abs(GENptZZ_Edge)>15.0&&abs(GENptZZ_Edge)<=20.0)+1.45511*(abs(GENptZZ_Edge)>20.0&&abs(GENptZZ_Edge)<=25.0)+1.47569*(abs(GENptZZ_Edge)>25.0&&abs(GENptZZ_Edge)<=30.0)+1.49053*(abs(GENptZZ_Edge)>30.0&&abs(GENptZZ_Edge)<=35.0)+1.50622*(abs(GENptZZ_Edge)>35.0&&abs(GENptZZ_Edge)<=40.0)+1.50328*(abs(GENptZZ_Edge)>40.0&&abs(GENptZZ_Edge)<=45.0)+1.52186*(abs(GENptZZ_Edge)>45.0&&abs(GENptZZ_Edge)<=50.0)+1.52043*(abs(GENptZZ_Edge)>50.0&&abs(GENptZZ_Edge)<=55.0)+1.53977*(abs(GENptZZ_Edge)>55.0&&abs(GENptZZ_Edge)<=60.0)+1.53491*(abs(GENptZZ_Edge)>60.0&&abs(GENptZZ_Edge)<=65.0)+1.51772*(abs(GENptZZ_Edge)>65.0&&abs(GENptZZ_Edge)<=70.0)+1.54494*(abs(GENptZZ_Edge)>70.0&&abs(GENptZZ_Edge)<=75.0)+1.57762*(abs(GENptZZ_Edge)>75.0&&abs(GENptZZ_Edge)<=80.0)+1.55078*(abs(GENptZZ_Edge)>80.0&&abs(GENptZZ_Edge)<=85.0)+1.57078*(abs(GENptZZ_Edge)>85.0&&abs(GENptZZ_Edge)<=90.0)+1.56162*(abs(GENptZZ_Edge)>90.0&&abs(GENptZZ_Edge)<=95.0)+1.54183*(abs(GENptZZ_Edge)>95.0&&abs(GENptZZ_Edge)<=100.0)+1.58485*(abs(GENptZZ_Edge)>100.0))"
              if doKfactorGENVar == 'ZZdPhi':
                  kf = "(1.513834489150*(abs(GENphiZZ_Edge)>0.0&&abs(GENphiZZ_Edge)<=0.1)+1.54173*(abs(GENphiZZ_Edge)>0.1&&abs(GENphiZZ_Edge)<=0.2)+1.49782*(abs(GENphiZZ_Edge)>0.2&&abs(GENphiZZ_Edge)<=0.3)+1.53495*(abs(GENphiZZ_Edge)>0.3&&abs(GENphiZZ_Edge)<=0.4)+1.47821*(abs(GENphiZZ_Edge)>0.4&&abs(GENphiZZ_Edge)<=0.5)+1.50433*(abs(GENphiZZ_Edge)>0.5&&abs(GENphiZZ_Edge)<=0.6)+1.52062*(abs(GENphiZZ_Edge)>0.6&&abs(GENphiZZ_Edge)<=0.7)+1.50701*(abs(GENphiZZ_Edge)>0.7&&abs(GENphiZZ_Edge)<=0.8)+1.49424*(abs(GENphiZZ_Edge)>0.8&&abs(GENphiZZ_Edge)<=0.9)+1.45053*(abs(GENphiZZ_Edge)>0.9&&abs(GENphiZZ_Edge)<=1.0)+1.46081*(abs(GENphiZZ_Edge)>1.0&&abs(GENphiZZ_Edge)<=1.1)+1.47160*(abs(GENphiZZ_Edge)>1.1&&abs(GENphiZZ_Edge)<=1.2)+1.46770*(abs(GENphiZZ_Edge)>1.2&&abs(GENphiZZ_Edge)<=1.3)+1.42240*(abs(GENphiZZ_Edge)>1.3&&abs(GENphiZZ_Edge)<=1.4)+1.39718*(abs(GENphiZZ_Edge)>1.4&&abs(GENphiZZ_Edge)<=1.5)+1.37559*(abs(GENphiZZ_Edge)>1.5&&abs(GENphiZZ_Edge)<=1.6)+1.39190*(abs(GENphiZZ_Edge)>1.6&&abs(GENphiZZ_Edge)<=1.7)+1.36856*(abs(GENphiZZ_Edge)>1.7&&abs(GENphiZZ_Edge)<=1.8)+1.31788*(abs(GENphiZZ_Edge)>1.8&&abs(GENphiZZ_Edge)<=1.9)+1.31401*(abs(GENphiZZ_Edge)>1.9&&abs(GENphiZZ_Edge)<=2.0)+1.27464*(abs(GENphiZZ_Edge)>2.0&&abs(GENphiZZ_Edge)<=2.1)+1.24234*(abs(GENphiZZ_Edge)>2.1&&abs(GENphiZZ_Edge)<=2.2)+1.24472*(abs(GENphiZZ_Edge)>2.2&&abs(GENphiZZ_Edge)<=2.3)+1.14625*(abs(GENphiZZ_Edge)>2.3&&abs(GENphiZZ_Edge)<=2.4)+1.10780*(abs(GENphiZZ_Edge)>2.4&&abs(GENphiZZ_Edge)<=2.5)+1.04205*(abs(GENphiZZ_Edge)>2.5&&abs(GENphiZZ_Edge)<=2.6)+0.97360*(abs(GENphiZZ_Edge)>2.6&&abs(GENphiZZ_Edge)<=2.7)+0.87216*(abs(GENphiZZ_Edge)>2.7&&abs(GENphiZZ_Edge)<=2.8)+0.73450*(abs(GENphiZZ_Edge)>2.8&&abs(GENphiZZ_Edge)<=2.9)+1.16315*(abs(GENphiZZ_Edge)>2.9&&abs(GENphiZZ_Edge)<=3.1416))"
              if doKfactorGENVar == 'noKFactor':
                  kf = "1"
          
          if (self.doKfactor == 0): #all other MC has no NNLO/NLO reweighting
              kf = "1"  
          cut = cut + "* ( " + str(self.lumWeight*lumi) + " * genWeight_Edge/abs(genWeight_Edge) * " + self.puWeight +" * "+self.SFWeight+" * "+self.btagWeight+" * "+extraWeight +" * "+kf + " )" 
      else: 
         addDataFilters = "&&(  (Flag_eeBadScFilter_Edge == 1  ))"
         cut = "("+ cut + addDataFilters+ ")" + "* (" + extraWeight +")"
      
      #fileD.Close()                                                                                                                
      if (self.doKfactor == 1): print "doing ", doKfactorGENVar, "for ZZ4l kfactor!"
      #print "kfactor is ", kf
      if (self.doKfactor == 2): print "doing ", doKfactorGENVar, "for  ZZ2l kfactor!"
      self.ttree.Project(h.GetName(), var, cut, options)
      for _bin in range(1, h.GetNbinsX()+2):
          h_of.SetBinContent(_bin, h.GetBinContent(_bin))
          h_of.SetBinError  (_bin, h.GetBinError  (_bin))
      #fileD.Close()                                                                                                                
      return (h_of if ofBin else h)

    
   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel, extraWeight):
   
     if(xmin == xmax) and (ymax == ymin):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),len(nbiny)-1, array('d', nbiny))
     elif (xmin == xmax):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),nbiny,ymin,ymax)
     elif (ymin == ymax):
        h = TH2F(name, "", nbinx,xmin,xmax,len(nbiny)-1, array('d', nbiny))
     else: 
        h = TH2F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax)
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     if(self.isData == 0):
        cut = cut + "* ( " + str(self.lumWeight*lumi) + " * genWeight_Edge/abs(genWeight_Edge) * " + self.puWeight + " * " + self.SFWeight + " * " + self.btagWeight + " * " +  self.triggWeight  + "*" + extraWeight + " )" 
     else: 
        cut = cut + "* ( " + extraWeight + ")"
     self.ttree.Project(name, var, cut, options) 
     return h

   def getTH3F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax, cut, options, xlabel, ylabel, zlabel, extraWeight):
   
     if(xmin == xmax) and (ymax == ymin) and (zmax == zmin):
        h = TH3F(name, "", len(nbinx)-1, array('d', nbinx), len(nbiny)-1, array('d', nbiny), len(nbinz)-1, array('d', nbinz))
     else: 
        h = TH3F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax)
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     h.GetZaxis().SetTitle(zlabel)
     
     if(self.isData == 0):
        cut = cut + "* ( " + str(self.lumWeight*lumi) + " * genWeight_Edge/abs(genWeight_Edge) * " + self.puWeight + " * " + self.SFWeight + " * " + self.btagWeight + " * " +  self.triggWeight  + "*" + extraWeight + " )" 
     else: 
        cut = cut + "* ( " + extraWeight + ")"
     self.ttree.Project(name, var, cut, options) 
     return h

class Block:
   'Common base class for all Sample Blocks'

   def __init__(self, name, label, color, isdata, doKfactor):
      self.name  = name
      self.color = color
      self.isData = isdata
      self.doKfactor = doKfactor
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

   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel, extraWeight,doKFactorGENVar):
   #def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel, ofBin = True, extraWeight='1',ylabel = "Events" , doKFactorGENVar = 'noKFactor'):
     for _is,s in enumerate(self.samples):
       
       AuxName = "auxT1_sample" + s.name
       haux = s.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel, extraWeight, doKFactorGENVar)
       if not _is:
          h = haux.Clone(name+'_blockHisto')
       else:
          h.Add(haux)
       del haux

     h.SetLineColor(self.color)
     h.SetMarkerColor(self.color)
     h.SetTitle(self.label)

     return h

   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel, extraWeight):
     if(xmin == xmax) and (ymax == ymin):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),len(nbiny)-1, array('d', nbiny))
     elif (xmin == xmax):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),nbiny,ymin,ymax)
     elif (ymin == ymax):
        h = TH2F(name, "", nbinx,xmin,xmax,len(nbiny)-1, array('d', nbiny))
     else: 
        h = TH2F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax)

     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     for s in self.samples:
     
       AuxName = "auxT2_block" + s.name
       haux = s.getTH2F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel, extraWeight)
       h.Add(haux)
       del haux

     return h   

   def getTH3F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax, cut, options, xlabel, ylabel, zlabel, extraWeight):
     if(xmin == xmax) and (ymax == ymin) and (zmax == zmin):
        h = TH3F(name, "", len(nbinx)-1, array('d', nbinx),len(nbiny)-1, array('d', nbiny),len(nbinz)-1, array('d', nbinz))
     else: 
        h = TH3F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax)

     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     h.GetZaxis().SetTitle(zlabel)
     
     for s in self.samples:
     
       AuxName = "auxT3_block" + s.name
       haux = s.getTH3F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax, cut, options, xlabel, ylabel, zlabel,extraWeight)
       h.Add(haux)
       del haux

     return h   

       

class Tree:
   'Common base class for a physics meaningful tree'

   def __init__(self, fileName, name, isdata, isScan = False, isOnEOS = 0):
      print fileName
      self.name  = name
      self.isData = isdata
      self.blocks = []
      self.isScan = isScan
      self.isOnEOS = isOnEOS
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
        flocation   = splitedLine[4]
        xsection    = float(splitedLine[5])
        isdata      = int(splitedLine[6])
        doKfactor   = int(splitedLine[7])

        color = 0
        plusposition = theColor.find("+")
        if(plusposition == -1):
          color = eval(theColor)
        else:
          color = eval(theColor[0:plusposition])
          color = color + int(theColor[plusposition+1:len(theColor)])

        sample = Sample(name, flocation, xsection, isdata, doKfactor, self.isScan, self.isOnEOS)
        coincidentBlock = [l for l in self.blocks if l.name == block]
        if(coincidentBlock == []):

          newBlock = Block(block, label, color, isdata, doKfactor)
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
     print "xlabel ", xlabel
     for b in self.blocks:
     
       AuxName = "auxStack_block_" + name + "_" + b.name
       haux = b.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel, "1", 'noKFactor')
       haux.SetFillColor(b.color)
       hs.Add(haux)
       del haux


     can_aux = TCanvas("can_%s_%s"%(name, b.name))
     can_aux.cd()
     hs.Draw()

     del can_aux

     if xmax != xmin:
       hs.GetXaxis().SetTitle(xlabel)
       b = int((xmax-xmin)/nbin)
       ylabel = "Events / " + str(b) + " GeV"
     else:     
       ylabel = "# events"
   
     hs.GetYaxis().SetTitle(ylabel)
     return hs   


   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel, extraWeight, doKFactorGENVar):
     
     for ib,b in enumerate(self.blocks):
       AuxName = "auxh1_block_" + name + "_" + b.name
       haux = b.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel, extraWeight, doKFactorGENVar)
       if not ib:
          h = haux.Clone(name+'_treeHisto')
       else:
          h.Add(haux)
       del haux
       
       return h

   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, extraWeight='1'):
     if(xmin == xmax) and (ymax == ymin):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),len(nbiny)-1, array('d', nbiny))
     elif (xmin == xmax):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),nbiny,ymin,ymax)
     elif (ymin == ymax):
        h = TH2F(name, "", nbinx,xmin,xmax,len(nbiny)-1, array('d', nbiny))
     else: 
        h = TH2F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax)
        
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     for b in self.blocks:
     
       AuxName = "aux_block" + name + "_" + b.name
       haux = b.getTH2F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel, extraWeight)
       h.Add(haux)
       del haux

     return h   

   def getTH3F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax, cut, options, xlabel, ylabel, zlabel, extraWeight='1'):
     if(xmin == xmax) and (ymax == ymin) and (zmax == zmin):
        h = TH3F(name, "", len(nbinx)-1, array('d', nbinx),len(nbiny)-1, array('d', nbiny), len(nbinz)-1, array('d', nbinz))
     else: 
        h = TH3F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax)
        
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     h.GetZaxis().SetTitle(zlabel)
     
     for b in self.blocks:
     
       AuxName = "aux_block" + name + "_" + b.name
       haux = b.getTH3F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax, cut, options, xlabel, ylabel, zlabel, extraWeight)
       h.Add(haux)
       del haux

     return h   


from ROOT import TGraphErrors
import sys


class collection:
    def __init__(self, bins, vname):
        self.bins  = bins
        self.vname = vname
        self.cen_mc, self.cen_mc_gr = 0, 0
        self.fwd_mc, self.fwd_mc_gr = 0, 0
        self.cen_da, self.cen_da_gr = 0, 0
        self.fwd_da, self.fwd_da_gr = 0, 0

    def getHisto(self, dataMC, eta):
        if   (dataMC, eta) == ('MC'  , 'central'): 
            return self.cen_mc
        elif (dataMC, eta) == ('MC'  , 'forward'): 
            return self.fwd_mc
        elif (dataMC, eta) == ('DATA', 'central'): 
            return self.cen_da
        elif (dataMC, eta) == ('DATA', 'forward'): 
            return self.fwd_da

    def getGraph(self, dataMC, eta):
        if   (dataMC, eta) == ('MC'  , 'central'): 
            return self.cen_mc_gr
        elif (dataMC, eta) == ('MC'  , 'forward'): 
            return self.fwd_mc_gr
        elif (dataMC, eta) == ('DATA', 'central'): 
            return self.cen_da_gr
        elif (dataMC, eta) == ('DATA', 'forward'): 
            return self.fwd_da_gr

    def setHisto(self, histo, dataMC, eta):
        if   (dataMC, eta) == ('MC'  , 'central'): 
            self.cen_mc = histo
            self.cen_mc_gr = TGraphErrors(histo)
        elif (dataMC, eta) == ('MC'  , 'forward'): 
            self.fwd_mc = histo
            self.fwd_mc_gr = TGraphErrors(histo)
        elif (dataMC, eta) == ('DATA', 'central'): 
            self.cen_da = histo
            self.cen_da_gr = TGraphErrors(histo)
        elif (dataMC, eta) == ('DATA', 'forward'): 
            self.fwd_da = histo
            self.fwd_da_gr = TGraphErrors(histo)
        else:
            print 'you are not calling setHisto correctly'

    def printValues(self):
        for _histo in [self.cen_mc, self.cen_da, self.fwd_mc, self.fwd_da]:
            if _histo == 0: continue
            if _histo in [self.cen_mc, self.cen_da]:
                eta = 'central'
            else:
                eta = 'forward'
            for _bin in range(1,_histo.GetNbinsX()+1):
                print '%-6s in [%.0f, %.0f] in %s: %.3f +- %.3f' %( self.vname,
                      _histo.GetXaxis().GetBinLowEdge(_bin),
                      _histo.GetXaxis().GetBinUpEdge (_bin),
                      eta,
                      _histo.GetBinContent(_bin),
                      _histo.GetBinError(_bin))

    def saveInFile(self, pattern, systErr, findBin = 0):
        print 'writing calculated values into file...'
        filename = 'ingredients.dat'
        f = open(filename, 'r')
        lines = f.readlines()
        newlines = []
        for line in lines:
            appended = False
            for t in ['MC', 'DATA'] if self.cen_da else ['MC']:
                if all(s in line for s in pattern+[t]):
                    newlines.append('%-6s \t %-15s %-6s \t %.4f \t %.4f \t %-6s \t %.4f \t %.4f \t %-6s\n' %(
                            line.split()[0], line.split()[1], t,
                            self.getHisto(t, 'central').GetBinContent(1 if not findBin else self.getHisto(t, 'central').FindBin(findBin)),
                            self.getHisto(t, 'central').GetBinError  (1 if not findBin else self.getHisto(t, 'central').FindBin(findBin)),
                            str(systErr),
                            self.getHisto(t, 'forward').GetBinContent(1 if not findBin else self.getHisto(t, 'forward').FindBin(findBin)),
                            self.getHisto(t, 'forward').GetBinError  (1 if not findBin else self.getHisto(t, 'forward').FindBin(findBin)),
                            str(systErr) ))
                    appended = True
    
            if not appended:
                newlines.append(line)
        f.close()
        g = open(filename, 'w')
        g.writelines(newlines)

class region():
    def __init__(self, name, cuts, rvars, bins, doData):
        self.name      = name
        self.cuts      = cuts
        self.rvars     = rvars
        self.bins      = bins
        self.doData    = doData
        if len(rvars) is not len(bins):
            print 'length of variables and bins has to be equal'
            sys.exit('exiting!')
        self.setVariables()

    def setVariables(self):
        for v in self.rvars:
            if v == 'mll':
                self.mll      = collection(self.bins[self.rvars.index(v)], v)
                self.mll_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.mll_pred = collection(self.bins[self.rvars.index(v)], v)
            if v == 'met':
                self.met      = collection(self.bins[self.rvars.index(v)], v)
                self.met_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.met_pred = collection(self.bins[self.rvars.index(v)], v)
            if v == 'nj' :
                self.nj       = collection(self.bins[self.rvars.index(v)], v)
                self.nj_dy       = collection(self.bins[self.rvars.index(v)], v)
                self.nj_pred  = collection(self.bins[self.rvars.index(v)], v)
            if v == 'nb' :
                self.nb       = collection(self.bins[self.rvars.index(v)], v)
                self.nb_dy       = collection(self.bins[self.rvars.index(v)], v)
                self.nb_pred  = collection(self.bins[self.rvars.index(v)], v)
            if v == 'jzb' :
                self.jzb      = collection(self.bins[self.rvars.index(v)], v)
                self.jzb_pred = collection(self.bins[self.rvars.index(v)], v)

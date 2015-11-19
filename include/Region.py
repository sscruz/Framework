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
        s_eta = 'cen' if eta == 'central' else 'fwd'
        dm = dataMC[:2].lower()
        return getattr(self, '%s_%s'%(s_eta, dm))

    def getGraph(self, dataMC, eta):
        s_eta = 'cen' if eta == 'central' else 'fwd'
        dm = dataMC[:2].lower()
        return getattr(self, '%s_%s_gr'%(s_eta, dm))

    def setHisto(self, histo, dataMC, eta):
        s_eta = 'cen' if eta == 'central' else 'fwd'
        dm = dataMC[:2].lower()
        setattr(self, '%s_%s'   %(s_eta, dm), histo)
        setattr(self, '%s_%s_gr'%(s_eta, dm), TGraphErrors(histo))

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
        if type(systErr) != list:
            systErr = [systErr, systErr]
        for line in lines:
            appended = False
            for t in ['MC', 'DATA'] if self.cen_da else ['MC']:
                if all(s in line for s in pattern+[t]):
                    newlines.append('%-6s      %-15s %-6s      %.4f      %.4f      %.4f      %.4f      %.4f      %.4f\n' %(
                            line.split()[0], line.split()[1], t,
                            self.getHisto(t, 'central').GetBinContent(1 if not findBin else self.getHisto(t, 'central').FindBin(findBin)),
                            self.getHisto(t, 'central').GetBinError  (1 if not findBin else self.getHisto(t, 'central').FindBin(findBin)),
                            systErr[0],
                            self.getHisto(t, 'forward').GetBinContent(1 if not findBin else self.getHisto(t, 'forward').FindBin(findBin)),
                            self.getHisto(t, 'forward').GetBinError  (1 if not findBin else self.getHisto(t, 'forward').FindBin(findBin)),
                            systErr[1] ))
                    appended = True
    
            if not appended:
                newlines.append(line)
        f.close()
        g = open(filename, 'w')
        g.writelines(newlines)
        g.close()

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
                self.mll_ra   = collection(self.bins[self.rvars.index(v)], v)
                self.mll_tt   = collection(self.bins[self.rvars.index(v)], v)
                self.mll_pred = collection(self.bins[self.rvars.index(v)], v)
            if v == 'met':
                self.met      = collection(self.bins[self.rvars.index(v)], v)
                self.met_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.met_ra   = collection(self.bins[self.rvars.index(v)], v)
                self.met_tt   = collection(self.bins[self.rvars.index(v)], v)
                self.met_pred = collection(self.bins[self.rvars.index(v)], v)
            if v == 'nj' :
                self.nj       = collection(self.bins[self.rvars.index(v)], v)
                self.nj_dy    = collection(self.bins[self.rvars.index(v)], v)
                self.nj_ra    = collection(self.bins[self.rvars.index(v)], v)
                self.nj_tt    = collection(self.bins[self.rvars.index(v)], v)
                self.nj_pred  = collection(self.bins[self.rvars.index(v)], v)
            if v == 'nb' :
                self.nb       = collection(self.bins[self.rvars.index(v)], v)
                self.nb_dy    = collection(self.bins[self.rvars.index(v)], v)
                self.nb_ra    = collection(self.bins[self.rvars.index(v)], v)
                self.nb_tt    = collection(self.bins[self.rvars.index(v)], v)
                self.nb_pred  = collection(self.bins[self.rvars.index(v)], v)
            if v == 'nvtx' :
                self.nvtx     = collection(self.bins[self.rvars.index(v)], v)
                self.nvtx_dy  = collection(self.bins[self.rvars.index(v)], v)
                self.nvtx_ra  = collection(self.bins[self.rvars.index(v)], v)
                self.nvtx_tt  = collection(self.bins[self.rvars.index(v)], v)
                self.nvtx_pred = collection(self.bins[self.rvars.index(v)], v)
            if v == 'jzb' :
                self.jzb      = collection(self.bins[self.rvars.index(v)], v)
                self.jzb_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.jzb_ra   = collection(self.bins[self.rvars.index(v)], v)
                self.jzb_tt   = collection(self.bins[self.rvars.index(v)], v)
                self.jzb_pred = collection(self.bins[self.rvars.index(v)], v)
            if v == 'mlb' :
                self.mlb      = collection(self.bins[self.rvars.index(v)], v)
                self.mlb_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.mlb_ra   = collection(self.bins[self.rvars.index(v)], v)
                self.mlb_tt   = collection(self.bins[self.rvars.index(v)], v)
                self.mlb_pred = collection(self.bins[self.rvars.index(v)], v)
            if v == 'dr' :
                self.dr      = collection(self.bins[self.rvars.index(v)], v)
                self.dr_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.dr_ra   = collection(self.bins[self.rvars.index(v)], v)
                self.dr_tt   = collection(self.bins[self.rvars.index(v)], v)
                self.dr_pred = collection(self.bins[self.rvars.index(v)], v)
            if v == 'l1iso' :
                self.l1iso      = collection(self.bins[self.rvars.index(v)], v)
                self.l1iso_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.l1iso_ra   = collection(self.bins[self.rvars.index(v)], v)
                self.l1iso_tt   = collection(self.bins[self.rvars.index(v)], v)
                self.l1iso_pred = collection(self.bins[self.rvars.index(v)], v)
            if v == 'l2iso' :
                self.l2iso      = collection(self.bins[self.rvars.index(v)], v)
                self.l2iso_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.l2iso_ra   = collection(self.bins[self.rvars.index(v)], v)
                self.l2iso_tt   = collection(self.bins[self.rvars.index(v)], v)
                self.l2iso_pred = collection(self.bins[self.rvars.index(v)], v)
            if v == 'l1pt' :
                self.l1pt      = collection(self.bins[self.rvars.index(v)], v)
                self.l1pt_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.l1pt_ra   = collection(self.bins[self.rvars.index(v)], v)
                self.l1pt_tt   = collection(self.bins[self.rvars.index(v)], v)
                self.l1pt_pred = collection(self.bins[self.rvars.index(v)], v)
            if v == 'l2pt' :
                self.l2pt      = collection(self.bins[self.rvars.index(v)], v)
                self.l2pt_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.l2pt_ra   = collection(self.bins[self.rvars.index(v)], v)
                self.l2pt_tt   = collection(self.bins[self.rvars.index(v)], v)
                self.l2pt_pred = collection(self.bins[self.rvars.index(v)], v)
            if v == 'ht' :
                self.ht      = collection(self.bins[self.rvars.index(v)], v)
                self.ht_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.ht_ra   = collection(self.bins[self.rvars.index(v)], v)
                self.ht_tt   = collection(self.bins[self.rvars.index(v)], v)
                self.ht_pred = collection(self.bins[self.rvars.index(v)], v)

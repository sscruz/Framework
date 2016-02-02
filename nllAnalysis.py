import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats
import math, sys, optparse, copy, re, array


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables
import nllAnalysis


def getLikelihood(listofhistos, inverse = False):
    lh = 1.
    for histo in listofhistos:
        r = histo.GetRandom(); 
        lh = lh*histo.GetBinContent(histo.FindBin(r)) 
    return -1.*math.log(lh if not inverse else 1./lh)

def comparePDFs(st=0):
    pdf_file = r.TFile('pdfs_version6.root')
    r.gStyle.SetOptStat(0)
    if   st==1: metvar = 'st' 
    elif st==2: metvar = 'zpt' 
    else: metvar = 'met'
    for varname in ['nll']:#'ldr', 'met', 'mlb', 'nll']:
        if   varname == 'ldr':
            pdf_histo    = copy.deepcopy(pdf_file.Get('emu_data_pdf_histo_ldr_ds_cuts_of_sr_met150'))
            pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_pdf_histo_ldr_ds_cuts_of_sr_met150'   ))
            ymax = 0.006
        elif varname == 'met':
            pdf_histo    = copy.deepcopy(pdf_file.Get('emu_data_pdf_histo_met_ds_cuts_of_sr_met150'))
            pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_pdf_histo_met_ds_cuts_of_sr_met150'   ))
            ymax = 0.02
        elif varname == 'mlb':
            pdf_histo    = copy.deepcopy(pdf_file.Get('emu_data_pdf_histo_mlb_ds_cuts_of_sr_met150'))
            pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_pdf_histo_mlb_ds_cuts_of_sr_met150'   ))
            ymax = 0.008
        elif varname == 'nll':
            ldr_pdf_histo    = copy.deepcopy(pdf_file.Get('emu_data_pdf_histo_ldr_ds_cuts_of_sr_met150')); ldr_pdf_histo   .Scale(1./ldr_pdf_histo   .Integral())
            ldr_pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_pdf_histo_ldr_ds_cuts_of_sr_met150'   )); ldr_pdf_histo_mc.Scale(1./ldr_pdf_histo_mc.Integral())
            met_pdf_histo    = copy.deepcopy(pdf_file.Get('emu_data_pdf_histo_%s_ds_cuts_of_sr_met150'%metvar)); met_pdf_histo   .Scale(1./met_pdf_histo   .Integral())
            met_pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_pdf_histo_%s_ds_cuts_of_sr_met150'   %metvar)); met_pdf_histo_mc.Scale(1./met_pdf_histo_mc.Integral())
            mlb_pdf_histo    = copy.deepcopy(pdf_file.Get('emu_data_pdf_histo_mlb_ds_cuts_of_sr_met150')); mlb_pdf_histo   .Scale(1./mlb_pdf_histo   .Integral())
            mlb_pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_pdf_histo_mlb_ds_cuts_of_sr_met150'   )); mlb_pdf_histo_mc.Scale(1./mlb_pdf_histo_mc.Integral())
            #pdf_histo    = r.TH1F('lh_histo_da', 'lh_histo_da', 51, 12., 36.48)
            #pdf_histo_mc = r.TH1F('lh_histo_mc', 'lh_histo_mc', 51, 12., 36.48)
            pdf_histo    = r.TH1F('lh_histo_da', 'lh_histo_da', 400, 12., 36.48)
            pdf_histo_mc = r.TH1F('lh_histo_mc', 'lh_histo_mc', 400, 12., 36.48)
            for i in range(500000):
                lh    = getLikelihood(ldr_pdf_histo   , met_pdf_histo   , mlb_pdf_histo   )
                lh_mc = getLikelihood(ldr_pdf_histo_mc, met_pdf_histo_mc, mlb_pdf_histo_mc)
                pdf_histo   .Fill(lh)
                pdf_histo_mc.Fill(lh_mc)
            ymax = 0.0300

        

        for i in range(pdf_histo.GetNbinsX()+1):
            pdf_histo.SetBinError(i, 0.); pdf_histo_mc.SetBinError(i, 0.)
        pdf_histo   .Scale(1./pdf_histo.Integral())   ; pdf_histo   .SetTitle(''); pdf_histo   .GetYaxis().SetRangeUser(1e-6, ymax); pdf_histo   .SetLineWidth(2);
        pdf_histo_mc.Scale(1./pdf_histo_mc.Integral()); pdf_histo_mc.SetTitle(''); pdf_histo_mc.GetYaxis().SetRangeUser(1e-6, ymax); pdf_histo_mc.SetLineWidth(2);

        pdf_histo.GetXaxis().SetTitle(varname)

        plot = Canvas.Canvas('signalTestsWithMet/pdf_comparison_%s%s'%(varname, '_'+metvar if st else ''), 'png,pdf', 0.7, 0.55, 0.95, 0.89, 600, 600)
        li = 0
        plot.addHisto(pdf_histo   , 'hist'     , 'pdf - data', 'L' , r.kBlack, 1,  li); li+=1
        plot.addHisto(pdf_histo_mc, 'hist,same', 'pdf - mc'  , 'L' , r.kRed  , 1,  li); li+=1
        plot.saveRatio(1, 0, 0, 1., pdf_histo, pdf_histo_mc, 0.5, 1.5)
    return plot, copy.deepcopy(pdf_histo)

def makePDFChecks(treeTT, cuts, treeSIG = 0, sigMasses = 0):
    pdf_file = r.TFile('pdfs_version12.root')
    r.gStyle.SetOptStat(0)
    plots = []; rets = [];
    lhvar = 'zpt'
    for varname in ['nll']:#'met', 'lepdr', 'sum_mlb', 'nll']:
        nbins = 50
        if   varname == 'lepdr':
            var = 'lepsDR_Edge'
            pdf_histo    = copy.deepcopy(pdf_file.Get('emu_data_pdf_histo_ldr_ds_cuts_of_sr_met150'))
            pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_pdf_histo_ldr_ds_cuts_of_sr_met150'   ))
            pdf_histo.Rebin(20); pdf_histo_mc.Rebin(20)
            xmin = pdf_histo.GetXaxis().GetXmin(); xmax = pdf_histo.GetXaxis().GetXmax(); xint = (xmax-xmin)/nbins; nbins = 49
            ymax = 0.10
        elif varname == 'sum_mlb':
            var = 'sum_mlb_Edge'
            pdf_histo    = copy.deepcopy(pdf_file.Get('emu_data_pdf_histo_mlb_ds_cuts_of_sr_met150'))
            pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_pdf_histo_mlb_ds_cuts_of_sr_met150'   ))
            pdf_histo.Rebin(20); pdf_histo_mc.Rebin(20)
            xmin = 0.; xmax = 750-15.0; xint = 15.0; nbins = 49
            ymax = 0.14
        elif varname == 'st':
            var = 'st_Edge'
            pdf_histo    = copy.deepcopy(pdf_file.Get('emu_data_pdf_histo_st_ds_cuts_of_sr_met150'))
            pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_pdf_histo_st_ds_cuts_of_sr_met150'   ))
            pdf_histo.Rebin(20); pdf_histo_mc.Rebin(20)
            xmin = 100.; xmax = 1000-18.0; xint = 18.0; nbins = 49
            ymax = 0.14
        elif varname == 'zpt':
            var = 'lepsZPt_Edge'
            pdf_histo    = copy.deepcopy(pdf_file.Get('emu_data_pdf_histo_zpt_ds_cuts_of_sr_met150'))
            pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_pdf_histo_zpt_ds_cuts_of_sr_met150'   ))
            pdf_histo.Rebin(20); pdf_histo_mc.Rebin(20)
            xmin = 0.; xmax = 600.-12.0; xint = 12.0; nbins = 49
            ymax = 0.25
        elif varname == 'met':
            var = 'met_pt'
            pdf_histo    = copy.deepcopy(pdf_file.Get('emu_data_pdf_histo_met_ds_cuts_of_sr_met150'))
            pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_pdf_histo_met_ds_cuts_of_sr_met150'   ))
            pdf_histo.Rebin(20); pdf_histo_mc.Rebin(20)
            xmin = 0.; xmax = 750-15.0; xint = 15.0; nbins = 49
            ymax = 0.32
        elif varname == 'invzpt':
            var = 'cum_ana_zpt_data_Edge'
            pdf_histo    = copy.deepcopy(pdf_file.Get('em_data_pdf_histo_met_ds_cuts_of_sr_met150'))
            pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_pdf_histo_met_ds_cuts_of_sr_met150'  ))
            pdf_histo.Rebin(20); pdf_histo_mc.Rebin(20)
            xmin = 0.; xmax = 1.; xint = xmax/50; xmax = xmax - xint; nbins = 49;
            ymax = 0.25
        elif varname == 'invall':
            var = 'TMath::Power(cum_ana_mlb_data_Edge*cum_ana_zpt_data_Edge*cum_ana_ldr_data_Edge*cum_ana_met_data_Edge,0.25)'
            pdf_histo    = copy.deepcopy(pdf_file.Get('em_data_pdf_histo_mlb_ds_cuts_of_sr_met150'))
            pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_pdf_histo_mlb_ds_cuts_of_sr_met150'  ))
            pdf_histo.Rebin(20); pdf_histo_mc.Rebin(20)
            xmin = 0.; xmax = 4.; xint = xmax/50; xmax = xmax - xint; nbins = 49;#-0.08; xint = 0.08; nbins = 49
            ymax = 0.25
            print 'numbers', xmin, xmax, xint, nbins
        elif varname == 'nll':
            pdfOrFit = 'pdf'
            var = '-1.*TMath::Log(lh_{foo}met_data_Edge*lh_{foo}mlb_data_Edge*lh_{foo}ldr_data_Edge*lh_{foo}zpt_data_Edge)'.format(foo='' if pdfOrFit == 'pdf' else 'ana_')
            ldr_pdf_histo    = copy.deepcopy(pdf_file.Get('em_data_%s_histo_ldr_ds_cuts_of_sr_met150'%pdfOrFit)); ldr_pdf_histo   .Scale(1./ldr_pdf_histo   .Integral())
            ldr_pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_%s_histo_ldr_ds_cuts_of_sr_met150'  %pdfOrFit)); ldr_pdf_histo_mc.Scale(1./ldr_pdf_histo_mc.Integral())
            zpt_pdf_histo    = copy.deepcopy(pdf_file.Get('em_data_%s_histo_zpt_ds_cuts_of_sr_met150'%pdfOrFit)); zpt_pdf_histo   .Scale(1./zpt_pdf_histo   .Integral())
            zpt_pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_%s_histo_zpt_ds_cuts_of_sr_met150'  %pdfOrFit)); zpt_pdf_histo_mc.Scale(1./zpt_pdf_histo_mc.Integral())
            mlb_pdf_histo    = copy.deepcopy(pdf_file.Get('em_data_%s_histo_mlb_ds_cuts_of_sr_met150'%pdfOrFit)); mlb_pdf_histo   .Scale(1./mlb_pdf_histo   .Integral())
            mlb_pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_%s_histo_mlb_ds_cuts_of_sr_met150'  %pdfOrFit)); mlb_pdf_histo_mc.Scale(1./mlb_pdf_histo_mc.Integral())
            met_pdf_histo    = copy.deepcopy(pdf_file.Get('em_data_%s_histo_met_ds_cuts_of_sr_met150'%pdfOrFit)); met_pdf_histo   .Scale(1./met_pdf_histo   .Integral())
            met_pdf_histo_mc = copy.deepcopy(pdf_file.Get('tt_mc_%s_histo_met_ds_cuts_of_sr_met150'  %pdfOrFit)); met_pdf_histo_mc.Scale(1./met_pdf_histo_mc.Integral())
            pdf_histo    = r.TH1F('lh_histo_da', 'lh_histo_da', 50, 15., 39.5)
            pdf_histo_mc = r.TH1F('lh_histo_mc', 'lh_histo_mc', 50, 15., 39.5)
            for i in range(500000):
                lh    = getLikelihood([ldr_pdf_histo   , zpt_pdf_histo   , mlb_pdf_histo   , met_pdf_histo   ], inverse = False)
                lh_mc = getLikelihood([ldr_pdf_histo_mc, zpt_pdf_histo_mc, mlb_pdf_histo_mc, met_pdf_histo_mc], inverse = False)
                #print lh_mc
                pdf_histo   .Fill(lh)
                pdf_histo_mc.Fill(lh_mc)
            xmin = 15.; xmax = 40.; nbins = 49; xint = (xmax-xmin)/(nbins+1); xmax = xmax - xint;
            ymax = 0.25
            print 'before everything', xmin, xmax, nbins, xint
            #xmin = -30.; xmax = 36.; xint = (xmax-xmin)/nbins; nbins = 50
            #ymax = 0.30

        for i in range(pdf_histo.GetNbinsX()+1):
            pdf_histo.SetBinError(i, 0.); pdf_histo_mc.SetBinError(i, 0.)
        pdf_histo   .Scale(1./pdf_histo.Integral())   ; pdf_histo   .SetTitle(''); pdf_histo   .GetYaxis().SetRangeUser(1e-6, ymax); pdf_histo   .SetLineWidth(2);
        pdf_histo_mc.Scale(1./pdf_histo_mc.Integral()); pdf_histo_mc.SetTitle(''); pdf_histo_mc.GetYaxis().SetRangeUser(1e-6, ymax); pdf_histo_mc.SetLineWidth(2);

        cuts_sr_sf = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonNoTriggerSF()])
        cuts_sr_of = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonNoTriggerOF()])

        r.gStyle.SetOptStat(0)
        h_tt_sr_sf  = treeTT.getTH1F(1., '%s_sr_sf_tt' %varname, var, nbins , xmin, xmax , cuts_sr_sf, '', varname)
        h_tt_sr_of  = treeTT.getTH1F(1., '%s_sr_of_tt' %varname, var, nbins , xmin, xmax , cuts_sr_of, '', varname)
        print 'nbins, mxin, xmax ', nbins, xmin, xmax
        print 'nbins, mxin, xmax histo ', h_tt_sr_sf.GetNbinsX(), h_tt_sr_sf.GetXaxis().GetXmin(), h_tt_sr_sf.GetXaxis().GetXmax()
        print 'nbins, mxin, xmax pdf histo ', pdf_histo.GetNbinsX(), pdf_histo.GetXaxis().GetXmin(), pdf_histo.GetXaxis().GetXmax()
        rets.append(copy.deepcopy(h_tt_sr_sf))
        rets.append(copy.deepcopy(h_tt_sr_of))

        sigHistos = []
        if treeSIG:
            for mass in sigMasses:
                newCut = cuts_sr_sf.replace(cuts.twoLeptons, 't.nPairLep_Edge > 0')
                addStr = '&& GenSusyMScan1 == %d && GenSusyMScan2 == %d'%(mass[0], mass[1])
                tmp_histo = treeSIG.getTH1F(1., '%s_sr_sf_sig' %varname, var, nbins , xmin, xmax , '('+newCut+addStr+')', '', varname)
                if tmp_histo.Integral():
                    tmp_histo.Scale(1./tmp_histo.Integral())
                    sigHistos.append(tmp_histo)
        print sigHistos[0].GetXaxis().GetXmin(), sigHistos[0].GetXaxis().GetXmax(),sigHistos[0].GetNbinsX()
                
        h_tt_sr_sf.SetMarkerSize(0.7); h_tt_sr_sf.SetMarkerStyle(20); h_tt_sr_sf.SetMarkerColor(r.kGreen); h_tt_sr_sf.Scale(1./h_tt_sr_sf.Integral())
        h_tt_sr_of.SetMarkerSize(0.7); h_tt_sr_of.SetMarkerStyle(20); h_tt_sr_of.SetMarkerColor(r.kBlue ); h_tt_sr_of.Scale(1./h_tt_sr_of.Integral())

        h_tt_sr_sf.GetYaxis().SetRangeUser(1e-6, ymax)

        plot = Canvas.Canvas('inversePDFTests/pdf_check_%s%s'%(varname, '_'+lhvar if varname == 'nll' else ''), 'png,pdf', 0.7, 0.55, 0.95, 0.89, 600, 600)
        li = 0
        if not 'inv' in varname: plot.addHisto(pdf_histo   , 'hist'     , 'pdf - data', 'L' , r.kBlack, 1,  li); li+=1
        if not 'inv' in varname: plot.addHisto(pdf_histo_mc, 'hist,same', 'pdf - mc'  , 'L' , r.kRed  , 1,  li); li+=1
        plot.addHisto(h_tt_sr_sf  , 'pe, same' , 'tt - SR SF', 'PL', r.kGreen, 1,  li); li+=1
        plot.addHisto(h_tt_sr_of  , 'pe, same' , 'tt - SR OF', 'PL', r.kBlue , 1,  li); li+=1
        for ind,histo in enumerate(sigHistos):
            print 'adding signal histo for mass %d/%d'%(sigMasses[ind][0], sigMasses[ind][1])
            histo.SetMarkerSize(0.7); histo.SetMarkerStyle(20); histo.SetMarkerColor(ind+4);
            plot.addHisto(histo, 'pe, same', '%d/%d'%(sigMasses[ind][0], sigMasses[ind][1])  , 'PL' , ind+5  , 1,  li); li+=1
            for c in [19., 19.5, 20., 20.5, 24., 24.5, 25., 25.5, 26.]:
                print 'cut %.1f: 1.-sig.eff %.2f bkg.eff %.2f' %(c, 1.-histo.Integral(1, histo.FindBin(c)), pdf_histo.Integral(1, pdf_histo.FindBin(c)))
        plot.saveRatio(1, 0, 0 if varname!='met' else 0 , 1., pdf_histo, [h_tt_sr_of, h_tt_sr_sf, pdf_histo_mc])
    return plot, rets

def makeLikelihoodControl(treeSIG, treeTT, cuts, cuts_norm, massx, massy):
    pdf_file = r.TFile('pdfs_version0.root')
    r.gStyle.SetOptStat(0)
    plots = []
    nbins =  50
    for varname in ['sumptmet']:
    #for varname in ['lepdr', 'sum_mlb', 'l1metdphi']:
        if   varname == 'lepdr':
            var = 'lepsDR_Edge'
            pdf_histo = copy.deepcopy(pdf_file.Get('pdf_lepsdr'))
            pdf_histo.Rebin(20)
            xmin = pdf_histo.GetXaxis().GetXmin(); xmax = pdf_histo.GetXaxis().GetXmax()
            ymax = 0.12
        elif varname == 'sum_mlb':
            var = 'sum_mlb_Edge'
            pdf_histo = copy.deepcopy(pdf_file.Get('pdf_summlb'))
            pdf_histo.Rebin(20)
            xmin = pdf_histo.GetXaxis().GetXmin(); xmax = pdf_histo.GetXaxis().GetXmax()
            ymax = 0.20
        elif varname == 'l1metdphi':
            var = 'metl1DPhi_Edge'
            pdf_histo = copy.deepcopy(pdf_file.Get('pdf_mldphi'))
            pdf_histo.Rebin(20)
            xmin = pdf_histo.GetXaxis().GetXmin(); xmax = pdf_histo.GetXaxis().GetXmax()
            ymax = 0.10
        elif varname == 'sumptmet':
            var = 'Lep1_pt_Edge+Lep2_pt_Edge+met_pt'
            pdf_histo = copy.deepcopy(pdf_file.Get('pdf_mldphi'))
            pdf_histo.Rebin(20)
            xmin = 0.; xmax =500.
            ymax = 0.10
        for i in range(pdf_histo.GetNbinsX()+1):
            pdf_histo.SetBinError(i, 0.)
            
        xint = (xmax-xmin)/nbins
        print 'xint', xint
        test_tt  = treeTT .getTH1F(1., '%s_tt' %varname, var, nbins-1 , xmin, (nbins)*xint if varname != 'sum_mlb' else 980, cuts.AddList([cuts_norm ]), '', varname)
        test_tt .Scale(1./test_tt .Integral())
        test_tt .SetLineWidth(1)
        test_tt.GetYaxis().SetRangeUser(0., ymax)
        pdf_histo.SetMarkerStyle(20)
        pdf_histo.SetMarkerSize (0.8)
        pdf_histo.SetMarkerColor(r.kGray+1)
        pdf_histo.Scale(1./pdf_histo.Integral())
        print 'pdf histo: xmin %.1f xmax %.2f nbins %.2f' %(xmin, xmax, pdf_histo.GetNbinsX())
        print xmax

        plot = Canvas.Canvas('signalTests/%s'%(varname), 'png,pdf', 0.6, 0.60, 0.90, 0.89)
        plot.addHisto(test_tt  , 'hist'  , 'TTBAR - MC'  , 'L' , r.kBlack , 1,  0)
        if varname != 'sumptmet':
            plot.addHisto(pdf_histo, 'p,same', 'pdf - data'  , 'PL' , r.kGray-2, 1,  1)
        ind = 0
        tmp_histos = []
        for m in massx:
            tmp_massx = m
            tmp_massy = massy[ind]
            print 'getting signal histo for mx %.0f and my %.0f' %(tmp_massx, tmp_massy)
            tmp_histo = treeSIG.getTH1F(1., '%s_sms'%varname, var, nbins-1 , xmin, (nbins)*xint if varname != 'sum_mlb' else 980, cuts.AddList([cuts_norm, 'GenSusyMScan1 == %.0f && GenSusyMScan2 == %.0f'%(tmp_massx, tmp_massy)]), '', varname)
            print tmp_histo.GetNbinsX(), tmp_histo.GetXaxis().GetXmin(), tmp_histo.GetXaxis().GetXmax()
            tmp_histo.Scale(1./tmp_histo.Integral())
            tmp_histo.SetLineWidth(2)
            plot.addHisto(tmp_histo, 'hist,same', 'SIGNAL %.0f/%.0f'%(tmp_massx, tmp_massy), 'L' , ind+2 , 1,  ind+2)
            tmp_histos.append(tmp_histo)
            ind += 1
        if varname != 'sumptmet':
            plot.saveRatio(1, 0, 0 , 1., test_tt, pdf_histo, 0., 2.0)
        else:
            plot.saveRatio(1, 0, 0 , 1., test_tt, tmp_histos, 0., 2.0)
        plots.append(plot)
        print '============================================================== \n Chi2:'
        print pdf_histo.Chi2Test(test_tt,'UW,P')
        print '=============================================================='
    return plots, pdf_histo


def makeLikelihoodPlot(varname, var, nbins, xmin, xmax, treeSIG, treeTT, cuts, cuts_norm, massx, massy):
    r.gStyle.SetOptStat(0)
    test_tt  = treeTT .getTH1F(1., '%s_tt' %varname, var, nbins , xmin, xmax, cuts.AddList([cuts_norm ]), '', varname)
    test_tt .Scale(1./test_tt .Integral())
    test_tt .SetLineWidth(2)
    test_tt.GetYaxis().SetRangeUser(0., 0.12)

    plot = Canvas.Canvas('signalTests/%s'%(varname), 'png,pdf', 0.6, 0.50, 0.85, 0.82)
    plot.addHisto(test_tt , 'hist'     , 'TTBAR'                          , 'L' , r.kBlack, 1,  0)
    ind = 0
    for m in massx:
        tmp_massx = m
        tmp_massy = massy[ind]
        print 'getting signal histo for mx %.0f and my %.0f' %(tmp_massx, tmp_massy)
        tmp_histo = treeSIG.getTH1F(1., '%s_sms'%varname, var, nbins , xmin, xmax, cuts.AddList([cuts_norm, 'GenSusyMScan1 == %.0f && GenSusyMScan2 == %.0f'%(tmp_massx, tmp_massy)]), '', varname)
        tmp_histo.Scale(1./tmp_histo.Integral())
        tmp_histo.SetLineWidth(2)
        plot.addHisto(tmp_histo, 'hist,same', 'SIGNAL %.0f/%.0f'%(tmp_massx, tmp_massy), 'L' , ind+2 , 1,  ind+1)
        ind += 1
    plot.save(1, 0, 0 , 1.)
    return plot

 
if __name__ == "__main__":


    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    (opts, args) = parser.parse_args()

    print 'running with these options'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)
    print ' \n\n'
    print 'Going to load the tree(s)...'

    global treeTT, treeSIG, treeDY, lumi
    ttDatasets = ['TTLep_pow']
    treeTT = Sample.Tree(helper.selectSamples(opts.sampleFile, ttDatasets, 'TT'), 'TT'  , 0, isScan = False)

    sigDatasets = ['SMS_T6bbllslepton_mSbottom-600To900_mLSP-200To800', 
                   'SMS_T6bbllslepton_mSbottom-400To550_mLSP-200To500' ]
    treeSIG = Sample.Tree(helper.selectSamples(opts.sampleFile, sigDatasets, 'SIG'), 'SIG'  , 0, isScan = True)

    cuts = CutManager.CutManager()

    lumi = 2.1
    ## this should be then sr-ID:m_slepton:m_sbottom for the final scan


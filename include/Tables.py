import math

def makeConciseTable(binnedSRincb, binnedSR0b, binnedSR1b, ingDA, ingMC, onZ):
    ret = []
    for eta in ['central', 'forward']:
        _eta = 'cen' if eta == 'central' else 'fwd'
        obshistoincb     = binnedSRincb.mll     .getHisto('DATA', eta)
        predhistoincb    = binnedSRincb.mll_pred.getHisto('DATA', eta)
        obshisto0b       = binnedSR0b.mll       .getHisto('DATA', eta)
        predhisto0b      = binnedSR0b.mll_pred  .getHisto('DATA', eta)
        obshisto1b       = binnedSR1b.mll       .getHisto('DATA', eta)
        predhisto1b      = binnedSR1b.mll_pred  .getHisto('DATA', eta)
        #tmp_histo_dy0b   = dyShapes['%db_%s_%s_binned'%(0, 'mc', eta)]
        #tmp_histo_dy1b   = dyShapes['%db_%s_%s_binned'%(1, 'mc', eta)]
        for i in range(1,obshisto0b.GetNbinsX()+1):
            if i == obshisto0b.GetNbinsX(): continue
            mr = 'lm' if i == 1 else 'bz' if i == 2 else 'oz' if i == 3 else 'az' if i == 4 else 'hm'
            tmp_rinout   = getattr(getattr(ingDA, 'rinout_dy_%s'%mr), '%s_val'%_eta)
            tmp_rinout_e = getattr(getattr(ingDA, 'rinout_dy_%s'%mr), '%s_err'%_eta) if mr != 'oz' else 0.
            onz_incb   = tmp_rinout*getattr(onZ, '%s_incb'%_eta)
            onz_0b     = tmp_rinout*getattr(onZ, '%s_0b'  %_eta)
            onz_1b     = tmp_rinout*getattr(onZ, '%s_1b'  %_eta)
            onz_incb_e = math.sqrt( (getattr(onZ, '%s_incb'%_eta)*tmp_rinout_e)**2 + (tmp_rinout*getattr(onZ, '%s_incb_e'%_eta))**2 )
            onz_0b_e   = math.sqrt( (getattr(onZ, '%s_0b'  %_eta)*tmp_rinout_e)**2 + (tmp_rinout*getattr(onZ, '%s_0b_e'  %_eta))**2 )
            onz_1b_e   = math.sqrt( (getattr(onZ, '%s_1b'  %_eta)*tmp_rinout_e)**2 + (tmp_rinout*getattr(onZ, '%s_1b_e'  %_eta))**2 )
            fs_incb = predhistoincb.GetBinContent(i); fs_incb_e = predhistoincb.GetBinError(i)
            fs_0b   = predhisto0b.GetBinContent(i)  ; fs_0b_e   = predhisto0b.GetBinError(i)
            fs_1b   = predhisto1b.GetBinContent(i)  ; fs_1b_e   = predhisto1b.GetBinError(i)
            tot_incb = fs_incb + onz_incb
            tot_0b   = fs_0b   + onz_0b
            tot_1b   = fs_1b   + onz_1b
            tot_incb_e = math.sqrt(fs_incb_e**2 + onz_incb_e**2)
            tot_0b_e   = math.sqrt(fs_0b_e  **2 + onz_0b_e  **2)
            tot_1b_e   = math.sqrt(fs_1b_e  **2 + onz_1b_e  **2)
            #print i, tmp_rinout, onz_0b, onz_0b_e

            ret.append('&  %6.1f $\\pm$ %6.1f  &  \\multirow{2}{*}{\\textbf{ %3d }} &  %6.1f $\\pm$ %6.1f  &  \\multirow{2}{*}{\\textbf{ %3d }} &   %6.1f $\\pm$ %6.1f  &  \\multirow{2}{*}{\\textbf{ %d }} \\\\ \r'%(
            tot_incb, tot_incb_e, obshistoincb.GetBinContent(i),
            tot_0b  , tot_0b_e  , obshisto0b  .GetBinContent(i),
            tot_1b  , tot_1b_e  , obshisto1b  .GetBinContent(i)))
            ret.append('& (%6.1f $\\pm$ %6.1f) &                                    & (%6.1f $\\pm$ %6.1f) &                                    &  (%6.1f $\\pm$ %6.1f) &                                   \\\\ \\cline{2-8} \r'%(
            onz_incb, onz_incb_e, onz_0b, onz_0b_e, onz_1b, onz_1b_e ))
    return ret

def makeConciseTableWith2b(binnedSRincb, binnedSR0b, binnedSR1b, binnedSR2b, ingDA, ingMC, onZ, dataMC='DATA'):
    ret = []
    isdata = (dataMC == 'DATA')
    for eta in ['central', 'forward']:
        _eta = 'cen' if eta == 'central' else 'fwd'
        obshistoincb     = binnedSRincb.mll     .getHisto(dataMC, eta)
        predhistoincb    = binnedSRincb.mll_pred.getHisto(dataMC, eta)
        obshisto0b       = binnedSR0b.mll       .getHisto(dataMC, eta)
        predhisto0b      = binnedSR0b.mll_pred  .getHisto(dataMC, eta)
        obshisto1b       = binnedSR1b.mll       .getHisto(dataMC, eta)
        predhisto1b      = binnedSR1b.mll_pred  .getHisto(dataMC, eta)
        obshisto2b       = binnedSR2b.mll       .getHisto(dataMC, eta)
        predhisto2b      = binnedSR2b.mll_pred  .getHisto(dataMC, eta)
        #tmp_histo_dy0b   = dyShapes['%db_%s_%s_binned'%(0, 'mc', eta)]
        #tmp_histo_dy1b   = dyShapes['%db_%s_%s_binned'%(1, 'mc', eta)]
        for i in range(1,obshisto0b.GetNbinsX()+1):
            if i == obshisto0b.GetNbinsX(): continue
            mr = 'lm' if i == 1 else 'bz' if i == 2 else 'oz' if i == 3 else 'az' if i == 4 else 'hm'
            tmp_rinout   = getattr(getattr(ingDA, 'rinout_dy_%s'%mr), '%s_val'%_eta)
            tmp_rinout_e = getattr(getattr(ingDA, 'rinout_dy_%s'%mr), '%s_err'%_eta) if mr != 'oz' else 0.
            onz_incb   = tmp_rinout*getattr(onZ, '%s_incb'%_eta)
            onz_0b     = tmp_rinout*getattr(onZ, '%s_0b'  %_eta)
            onz_1b     = tmp_rinout*getattr(onZ, '%s_1b'  %_eta)
            onz_2b     = tmp_rinout*getattr(onZ, '%s_2b'  %_eta)
            onz_incb_e = math.sqrt( (getattr(onZ, '%s_incb'%_eta)*tmp_rinout_e)**2 + (tmp_rinout*getattr(onZ, '%s_incb_e'%_eta))**2 )
            onz_0b_e   = math.sqrt( (getattr(onZ, '%s_0b'  %_eta)*tmp_rinout_e)**2 + (tmp_rinout*getattr(onZ, '%s_0b_e'  %_eta))**2 )
            onz_1b_e   = math.sqrt( (getattr(onZ, '%s_1b'  %_eta)*tmp_rinout_e)**2 + (tmp_rinout*getattr(onZ, '%s_1b_e'  %_eta))**2 )
            onz_2b_e   = math.sqrt( (getattr(onZ, '%s_2b'  %_eta)*tmp_rinout_e)**2 + (tmp_rinout*getattr(onZ, '%s_2b_e'  %_eta))**2 )
            fs_incb = predhistoincb.GetBinContent(i); fs_incb_e = predhistoincb.GetBinError(i)
            fs_0b   = predhisto0b.GetBinContent(i)  ; fs_0b_e   = predhisto0b.GetBinError(i)
            fs_1b   = predhisto1b.GetBinContent(i)  ; fs_1b_e   = predhisto1b.GetBinError(i)
            fs_2b   = predhisto2b.GetBinContent(i)  ; fs_2b_e   = predhisto2b.GetBinError(i)
            tot_incb = fs_incb + onz_incb
            tot_0b   = fs_0b   + onz_0b
            tot_1b   = fs_1b   + onz_1b
            tot_2b   = fs_2b   + onz_2b
            tot_incb_e = math.sqrt(fs_incb_e**2 + onz_incb_e**2)
            tot_0b_e   = math.sqrt(fs_0b_e  **2 + onz_0b_e  **2)
            tot_1b_e   = math.sqrt(fs_1b_e  **2 + onz_1b_e  **2)
            tot_2b_e   = math.sqrt(fs_2b_e  **2 + onz_2b_e  **2)
            #print i, tmp_rinout, onz_0b, onz_0b_e

            ret.append('&  %6.1f $\\pm$ %6.1f  &  \\multirow{2}{*}{\\textbf{ %3d }} &  %6.1f $\\pm$ %6.1f  &  \\multirow{2}{*}{\\textbf{ %3d }} &   %6.1f $\\pm$ %6.1f  &  \\multirow{2}{*}{\\textbf{ %d }} &   %6.1f $\\pm$ %6.1f  &  \\multirow{2}{*}{\\textbf{ %d }} \\\\ \r'%(
            tot_incb, tot_incb_e, obshistoincb.GetBinContent(i),
            tot_0b  , tot_0b_e  , obshisto0b  .GetBinContent(i),
            tot_1b  , tot_1b_e  , obshisto1b  .GetBinContent(i),
            tot_2b  , tot_2b_e  , obshisto2b  .GetBinContent(i)))
            ret.append('& (%6.1f $\\pm$ %6.1f) &                                    & (%6.1f $\\pm$ %6.1f) &                                    &  (%6.1f $\\pm$ %6.1f) &                                   &  (%6.1f $\\pm$ %6.1f) &                                   \\\\ \\cline{2-10} \r'%(
            onz_incb, onz_incb_e, onz_0b, onz_0b_e, onz_1b, onz_1b_e, onz_2b, onz_2b_e ))
    return ret


def makeRSFOFTable(ingDA, ingMC):
    ret = []
    ret.append('\\begin{table}[hbtp]')
    ret.append('\\begin{center}')
    ret.append('\\bgroup')
    ret.append('\\def\\arraystretch{1.2}')
    #ret.append('\\small')

    ret.append('\\caption{Calulation of all RSFOF numbers}')
    ret.append('\\label{tab:combinedRSFOF}')
    ret.append('\\begin{tabular}{l| c c| c c }')
    ret.append('\\multicolumn{1}{c}{} & \\multicolumn{2}{c}{\\textbf{Central}} & \\multicolumn{2}{c}{\\textbf{Forward}} \\\\ \\cline{2-5}')
    ret.append('& Data & MC & Data & MC \\\\ \\hline')
    ret.append('$\\frac{1}{2}$ $( r_{\\mu/e} + r_{\\mu/e}^{-1} )$   &  %.3f $\\pm$ %.3f  &   %.3f$\\pm$ %.3f   &   %.3f$\\pm$ %.3f &    %.3f$\\pm$ %.3f    \\\\'%(ingDA.rmue_factor.cen_val, ingDA.rmue_factor.cen_err, ingMC.rmue_factor.cen_val, ingMC.rmue_factor.cen_err, 
                                                                                                                                                                               ingDA.rmue_factor.fwd_val, ingDA.rmue_factor.fwd_err, ingMC.rmue_factor.fwd_val, ingMC.rmue_factor.fwd_err))
    ret.append('$R_{T}$                          &  %.3f $\\pm$ %.3f  &   %.3f$\\pm$ %.3f   &   %.3f$\\pm$ %.3f &    %.3f$\\pm$ %.3f    \\\\ \\hline'%(ingDA.rt_region.cen_val  , ingDA.rt_region.cen_err  , ingMC.rt_region.cen_val  , ingMC.rt_region.cen_err  , 
                                                                                                                                                       ingDA.rt_region.fwd_val  , ingDA.rt_region.fwd_err  , ingMC.rt_region.fwd_val  , ingMC.rt_region.fwd_err  ))
    ret.append('& \\multicolumn{4}{c}{\\Rsfof}  \\\\ \\hline')
    ret.append('from factorization        &  %.3f $\\pm$ %.3f   &  %.3f $\\pm$ %.3f       &  %.3f $\\pm$ %.3f  &   %.3f $\\pm$ %.3f     \\\\'%(ingDA.rsfof_fac_cen, ingDA.rsfof_fac_cen_e, ingMC.rsfof_fac_cen, ingMC.rsfof_fac_cen_e, 
                                                                                                                                               ingDA.rsfof_fac_fwd, ingDA.rsfof_fac_fwd_e, ingMC.rsfof_fac_fwd, ingMC.rsfof_fac_fwd_e))
    ret.append('     direct measurement   &  %.3f $\\pm$ %.3f   &  %.3f $\\pm$ %.3f       &  %.3f $\\pm$ %.3f  &   %.3f $\\pm$ %.3f     \\\\ \\hline'%(ingDA.rsfof_dir_cen, ingDA.rsfof_dir_cen_e, ingMC.rsfof_dir_cen, ingMC.rsfof_dir_cen_e, 
                                                                                                                                               ingDA.rsfof_dir_fwd, ingDA.rsfof_dir_fwd_e, ingMC.rsfof_dir_fwd, ingMC.rsfof_dir_fwd_e))
    ret.append('weighted average    &  \\textbf{%.3f $\\pm$ %.3f}   &  \\textbf{%.3f $\\pm$ %.3f}       &  \\textbf{%.3f $\\pm$ %.3f}  &   \\textbf{%.3f $\\pm$ %.3f}     \\\\'%(ingDA.rsfof_final_cen_val, ingDA.rsfof_final_cen_err, ingMC.rsfof_final_cen_val, ingMC.rsfof_final_cen_err, 
                                                                                                                                         ingDA.rsfof_final_fwd_val, ingDA.rsfof_final_fwd_err, ingMC.rsfof_final_fwd_val, ingMC.rsfof_final_fwd_err))
    ret.append('\\end{tabular}')
    ret.append('\\egroup')
    ret.append('\\end{center}')
    ret.append('\\end{table}')
    for i in ret: 
        print i
    return ret


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

            ret.append('&  %6.2f $\\pm$ %6.2f  &  \\multirow{2}{*}{\\textbf{ %3d }} &  %6.2f $\\pm$ %6.2f  &  \\multirow{2}{*}{\\textbf{ %3d }} &   %6.2f $\\pm$ %6.2f  &  \\multirow{2}{*}{\\textbf{ %d }} \\\\ \r'%(
            tot_incb, tot_incb_e, obshistoincb.GetBinContent(i),
            tot_0b  , tot_0b_e  , obshisto0b  .GetBinContent(i),
            tot_1b  , tot_1b_e  , obshisto1b  .GetBinContent(i)))
            ret.append('& (%6.2f $\\pm$ %6.2f) &                                    & (%6.2f $\\pm$ %6.2f) &                                    &  (%6.2f $\\pm$ %6.2f) &                                   \\\\ \\cline{2-8} \r'%(
            onz_incb, onz_incb_e, onz_0b, onz_0b_e, onz_1b, onz_1b_e ))
    return ret


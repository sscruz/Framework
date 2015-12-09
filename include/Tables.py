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

def makeDataCards(binnedSRList, dc_name, ingDA, onZ):
    print 'preparing datacards'

    for sr in binnedSRList:
        print '   ... for region %s'%(sr.name)
        b_string = 'incb' if 'inc' in sr.name else '0b' if '0b' in sr.name else '1b' if '1b' in sr.name else '2b' if '2b' in sr.name else 'nimps'
        ret = []
        for eta in ['central', 'forward']:
            print '      ... in %s'%eta
            _eta = 'cen' if eta == 'central' else 'fwd'
            obshisto = sr.mll     .getHisto('DATA', eta)
            predhisto= sr.mll_pred.getHisto('DATA', eta)
            for i in range(1,obshisto.GetNbinsX()+1):
                if i == obshisto.GetNbinsX(): continue
                mr   = 'lm'      if i == 1 else 'bz'     if i == 2 else 'oz'  if i == 3 else 'az'     if i == 4 else 'hm'
                mass = 'lowMass' if i == 1 else 'belowZ' if i == 2 else 'onZ' if i == 3 else 'aboveZ' if i == 4 else 'highMass' if i == 5 else 'tooHigh'
                print '         ... for mass region %s'%mass
                # get rinout and onZ prediction
                tmp_rinout   = getattr(getattr(ingDA, 'rinout_dy_%s'%mr), '%s_val'%_eta)
                tmp_rinout_e = getattr(getattr(ingDA, 'rinout_dy_%s'%mr), '%s_err'%_eta) if mr != 'oz' else 0.
                onz          = tmp_rinout*getattr(onZ, '%s_%s'%(_eta, b_string))
                onz_e        = math.sqrt((getattr(onZ, '%s_%s'%(_eta, b_string))*tmp_rinout_e)**2 + (tmp_rinout*getattr(onZ, '%s_%s_e'%(_eta, b_string)))**2 )
                # get FS background
                fs           = predhisto.GetBinContent(i)
                fs_e         = predhisto.GetBinError(i)
                # get the observed
                obs          = obshisto.GetBinContent(i)
                # get opposite flavor yield
                #fs_bkg       = predhisto.GetBinContent(i) ## this is wrong obviously!
                of_histo = sr.mll_of_central if eta == 'central' else sr.mll_of_forward
                of_yield     = of_histo.GetBinContent(i)
                # make a meaningful bin name
                bin_name = '%s_%s_%s'%(eta, mass, b_string)
                # get the final RSFOF
                rsfof   = getattr(ingDA, 'rsfof_final_%s_val'%(_eta))
                rsfof_e = getattr(ingDA, 'rsfof_final_%s_err'%(_eta))
                dc = '''# this is the datacard for bin {bin_name}
imax 1  number of channels
jmax 2  number of backgrounds
kmax *  number of nuisance parameters (sources of systematical uncertainties)
------------
# only one channel here with observation {obs}
bin            {bin_name}
observation    {obs}
------------
bin        {bin_name}     {bin_name}     {bin_name}
process    sig            FS             DY    
process    0              1              2
rate       XXRATEXX         {fs_bkg:.2f}       {dy_bkg:.2f}
------------
deltaS       lnN              1.15      -           -       
lumi         lnN              1.12      -           -       
FS_unc       lnN              -         {fs_unc:.2f}    -       
{bin_name}_fs_stat      gmN {of_yield}   -         {rsfof:.3f}        -       
DY_unc       lnN              -         -           {dy_unc:.2f}'''.format(bin_name=bin_name, obs=obs, fs_bkg=fs, fs_unc=1+rsfof_e, dy_bkg=onz, dy_unc=1+onz_e/onz, of_yield=int(of_yield), rsfof=rsfof)
                tmp_file = open('datacards/%s.txt'%(bin_name),'w')
                tmp_file.write(dc)
                tmp_file.close()

    print 'done with datacards'


def makeRinoutTable(region):
    mc_cen = region.mll.getHisto('MC'  , 'central')
    mc_fwd = region.mll.getHisto('MC'  , 'forward')
    da_cen = region.mll.getHisto('DATA', 'central')
    da_fwd = region.mll.getHisto('DATA', 'forward')

    table = '''
\\begin{{table}}[ht!]
\\bgroup
\\def\\arraystretch{{1.2}}
\\caption{{Measured values for \\rinout for data and MC in the different signal regions of the off-Z analysis.}}
\\label{{tab:rinoutvalues}}
\\begin{{center}}
\\begin{{tabular}}{{ l c c c c c}}
\\hline \\hline
      eta region         &       &  \\mll in [20, 70]  &   \\mll in [70, 81]   & \\mll in [81, 120] & \\mll in [120, 300]  \\\\ \\hline
\\multirow{{2}}{{*}}{{central}} &  MC   &  {mc_lm_cen:.3f} $\\pm$ {mc_lm_cen_e:.3f} &   {mc_bz_cen:.3f} $\\pm$ {mc_bz_cen_e:.3f}  & {mc_az_cen:.3f} $\\pm$ {mc_az_cen_e:.3f} & {mc_hm_cen:.3f} $\\pm$ {mc_hm_cen_e:.3f} \\\\
                               &  data &  {da_lm_cen:.3f} $\\pm$ {da_lm_cen_e:.3f} &   {da_bz_cen:.3f} $\\pm$ {da_bz_cen_e:.3f}  & {da_az_cen:.3f} $\\pm$ {da_az_cen_e:.3f} & {da_hm_cen:.3f} $\\pm$ {da_hm_cen_e:.3f} \\\\
\\multirow{{2}}{{*}}{{forward}} &  MC   &  {mc_lm_fwd:.3f} $\\pm$ {mc_lm_fwd_e:.3f} &   {mc_bz_fwd:.3f} $\\pm$ {mc_bz_fwd_e:.3f}  & {mc_az_fwd:.3f} $\\pm$ {mc_az_fwd_e:.3f} & {mc_hm_fwd:.3f} $\\pm$ {mc_hm_fwd_e:.3f} \\\\
                               &  data &  {da_lm_fwd:.3f} $\\pm$ {da_lm_fwd_e:.3f} &   {da_bz_fwd:.3f} $\\pm$ {da_bz_fwd_e:.3f}  & {da_az_fwd:.3f} $\\pm$ {da_az_fwd_e:.3f} & {da_hm_fwd:.3f} $\\pm$ {da_hm_fwd_e:.3f} \\\\
\\hline\\hline
\\end{{tabular}}
\\end{{center}}
\\egroup
\\end{{table}}
    '''.format(mc_lm_cen=mc_cen.GetBinContent(1), mc_lm_cen_e=mc_cen.GetBinError(1), mc_lm_fwd=mc_fwd.GetBinContent(1), mc_lm_fwd_e=mc_fwd.GetBinError(1),
               mc_bz_cen=mc_cen.GetBinContent(2), mc_bz_cen_e=mc_cen.GetBinError(2), mc_bz_fwd=mc_fwd.GetBinContent(2), mc_bz_fwd_e=mc_fwd.GetBinError(2),
               mc_az_cen=mc_cen.GetBinContent(4), mc_az_cen_e=mc_cen.GetBinError(4), mc_az_fwd=mc_fwd.GetBinContent(4), mc_az_fwd_e=mc_fwd.GetBinError(4),
               mc_hm_cen=mc_cen.GetBinContent(5), mc_hm_cen_e=mc_cen.GetBinError(5), mc_hm_fwd=mc_fwd.GetBinContent(5), mc_hm_fwd_e=mc_fwd.GetBinError(5),
               da_lm_cen=da_cen.GetBinContent(1), da_lm_cen_e=da_cen.GetBinError(1), da_lm_fwd=da_fwd.GetBinContent(1), da_lm_fwd_e=da_fwd.GetBinError(1),
               da_bz_cen=da_cen.GetBinContent(2), da_bz_cen_e=da_cen.GetBinError(2), da_bz_fwd=da_fwd.GetBinContent(2), da_bz_fwd_e=da_fwd.GetBinError(2),
               da_az_cen=da_cen.GetBinContent(4), da_az_cen_e=da_cen.GetBinError(4), da_az_fwd=da_fwd.GetBinContent(4), da_az_fwd_e=da_fwd.GetBinError(4),
               da_hm_cen=da_cen.GetBinContent(5), da_hm_cen_e=da_cen.GetBinError(5), da_hm_fwd=da_fwd.GetBinContent(5), da_hm_fwd_e=da_fwd.GetBinError(5))
    print table
    return table




import numpy as np
import VP4Dhyb_GYSutils
import scipy as sp
GYSut = VP4Dhyb_GYSutils.impall('VP4Dhyb_GYSutils')




#-----------------------------------------------------
# Computation of the growth rate and the frequency
#-----------------------------------------------------
def GR_Freq_Comp(f, fsquare, itime_init_interp, itime_end_interp, Nk, time, w_min, square):
    dt              = time[1] - time[0]
    time_tmp        = time[itime_init_interp:itime_end_interp+1].T[0]
    f_tmp           = f[itime_init_interp:itime_end_interp+1]
    Lx              = time_tmp[-1] - time_tmp[0]
    # Computation of the first approximation of omega: if square != 0, we must divide by 2
    [TFFPhi, Ffreq] = GYSut.Fourier1D( f_tmp, time_tmp )
    omegaapp        = np.abs(Ffreq[np.max(np.nonzero(np.abs(TFFPhi)*(Ffreq>w_min)==np.max(np.abs(TFFPhi)*(Ffreq>w_min))))])/(2 - (square == 0))
    # Computation of the phase of the curve
    nend            = int(2 * sp.pi / (omegaapp * dt))
    fsquare_ct      = np.sqrt(fsquare[itime_init_interp + range(0,nend+1)].T[0])
    iphase          = np.nonzero(fsquare_ct==np.min(fsquare_ct))[0]
    phase           = time_tmp[iphase][0]
    fsquare_tmp     = np.sqrt(fsquare[range(itime_init_interp,itime_init_interp+iphase,1) + range(itime_init_interp+iphase+1,itime_end_interp+1,1)].T[0])
    time_tmp_2      = time[range(itime_init_interp,itime_init_interp+iphase,1) + range(itime_init_interp+iphase+1,itime_end_interp+1,1)].T[0]
    # Array in which results will be placed
    Results         = np.zeros((Nk+1, 6), dtype=float, order='C')
    # Results will be computed for different omegas
    for d_omega in range(Nk+1):
        omega           = omegaapp + 2 * sp.pi / Lx * (2 * float(d_omega) / float(Nk) - 1)
        # Computation of the growth rate
        cos_omega_t     = np.abs(np.sin(omega*(time_tmp_2-phase)))
        fsquare_tmp_2   = fsquare_tmp / cos_omega_t
        poly_resu       = np.polyfit( time_tmp_2, np.log(fsquare_tmp_2), 1, full='True' )
        [p,s]           = poly_resu[0]
        residuals       = poly_resu[1]
        rank            = poly_resu[2]
        singular_values = poly_resu[3]
        rcond           = poly_resu[4]
        # Results of the loop are placed in Results
        Results[d_omega,0]    = np.linalg.norm(np.exp(s)*np.exp(p*time_tmp_2)*np.abs(np.sin(omega*(time_tmp_2-phase))) - fsquare_tmp)/np.sqrt(Lx)
        Results[d_omega,1:6]  = [omega, p, s, phase, residuals]
    #end for

    # The best result is taken
    L2_diff         = np.amin(Results[:,0])
    Ind_best        = np.nonzero(Results[:,0] == L2_diff)
    ib              = Ind_best[0][0]
    omega           = Results[ib,1]
    p               = Results[ib,2]
    s               = Results[ib,3]
    phase           = Results[ib,4]

    return [omega, p, s, L2_diff, phase, residuals]

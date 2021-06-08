'''
More diagnostics for FBPIC
'''
from openpmd_viewer import OpenPMDTimeSeries
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, pi, m_e, e
from scipy.optimize import curve_fit

sqrt_2pi=np.sqrt(2.*pi)

def gaussian_profile_square( r, peak, w ):
    '''
        returns the square of a gaussian profile
    '''
    return peak*np.exp(-np.square(r)/(w*w/2.))

def get_w_EL_vs_t(path, laser_pol='x', if_save_txt=False, if_save_plot=False, save_plot_type='pdf', peak_method=1):
    '''
        Get w and laser E-field amplitude vs t using the projection of laser electric field squared.
        laser_pol: laser polarization direction, can be 'x' or 'y'.
        peak_method: the method for getting the peak E-field
                     0 means find the maximum axial value of the absolute value of E-field
                     1 means find the amplitude of fit and normalize to the amplitude of E-field
        Returns:
        t: array of time
        w: array of laser radius
        EL: array of laser E-field amplitude
    '''
    ts = OpenPMDTimeSeries(path+"/hdf5", check_all_files=False)
    len_iter = len(ts.iterations)
    w_array = np.zeros(len_iter)
    peak_array = np.zeros(len_iter)
    for i in range(len_iter):
        F, info = ts.get_field(field='E', coord=laser_pol, iteration=ts.iterations[i], m='all')
        F_square_avg = np.average(np.square(F), axis=1)
        # only take half of F
        F_square_avg=F_square_avg[int(F_square_avg.size/2):]
        # only take half of r
        rspread=info.r[int(info.r.size/2):]
        if i<1:
            # Save initial peak E field
            EL0 = np.max(np.abs(F))
            # initialize w_guess
            #calculate RMS for the initial guess of fit
            w_guess = np.average(rspread, weights=F_square_avg)*sqrt_2pi
            w_guess = np.abs(w_guess)
            # in some cases the RMS can be very large, then manually reduce it
            half_rmax=rspread[-1]/2.
            w_guess = min(w_guess, half_rmax)
        popt, pcov = curve_fit( gaussian_profile_square, rspread, F_square_avg, p0=[F_square_avg[0], w_guess], bounds=((F_square_avg[0]*0.5, w_guess*0.1), (F_square_avg[0]*2., w_guess*2.)))
        w_array[i] = popt[1]
        if 0==peak_method:
            peak_array[i] = np.max(np.abs(F[F.shape[0]//2,:]))
        else:
            peak_array[i] = popt[0]
        # set w_guess for the next iteration
        w_guess = popt[1]

    if 1==peak_method:
        peak_array = np.sqrt(peak_array)
        # normalize laser electric field according to initial value of peak electric field
        peak_array = peak_array*(EL0/peak_array[0])

    if if_save_txt:
        data_save_name=path+"/w_vs_t.data"
        data=np.transpose(np.reshape(np.concatenate((ts.t, w_array, peak_array)), (3, len_iter)))
        np.savetxt(data_save_name, data, delimiter='    ', fmt='%1.5e')
        print('Data saved at '+data_save_name)
    if if_save_plot:
        # matplotlib.use('Agg') is to avoid calling X display while using plt.plot()
        # In a non-interactive mode such as job submission,
        # this is required to avoid error while calling X display.
        import matplotlib
        matplotlib.use('Agg')
        h_fig, h_ax1 = plt.subplots()
        color1 = 'k'
        h_ax1.plot(ts.t, w_array, color=color1)
        h_ax1.set_xlabel("$t$ [s]")
        h_ax1.set_ylabel("$w$ [m]", color=color1)
        h_ax1.tick_params(axis='y', labelcolor=color1)
        h_ax2 = h_ax1.twinx()
        color2 = 'r'
        h_ax2.plot(ts.t, peak_array, color=color2)
        h_ax2.set_ylabel("$E_L$ [m]", color=color2)
        h_ax2.tick_params(axis='y', labelcolor=color2)
        plt.tight_layout()
        plot_save_name=path+"/w_vs_t."+save_plot_type
        plt.savefig(plot_save_name)
        print('Plot saved at '+plot_save_name)

    return ts.t, w_array, peak_array

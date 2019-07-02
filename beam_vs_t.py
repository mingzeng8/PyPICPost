import outfile
import LorentzianFit as lf
import numpy as np
import matplotlib.pyplot as plt

# plot beam parameters vs time

def get_beam_parameters(os_file, n0_per_cc, gamma_threshold, gamma_spread_method="lfit"):
    ''' calculate beam parameters according to the raw output data
        return q_pC (charge in pico-coulomb), emittances (both transverse directions), popt[0] (center gamma), popt[2]*2. (absolute Full-Width-Half-Maximum spread)
    '''
    os_file.read_raw_x1()
    os_file.read_raw_x2()
    os_file.read_raw_x3()
    os_file.read_raw_p1()
    os_file.read_raw_p2()
    os_file.read_raw_p3()
    os_file.read_raw_q()
    os_file.read_raw_ene(ene_key_warning=False)
    os_file.select_raw_data(ene_low=gamma_threshold-1.)
    '''condlist = [os_file._raw_ene>(gamma_threshold-1.)]
    os_file._raw_q = np.select(condlist, [os_file._raw_q])'''
    if 0==len(os_file._raw_select_index):
        print("Warning: no particle is selected.")
        return os_file.time, 0.0, [0.0,0.0], 1.0, 00.5, False
    q_pC = os_file.calculate_q_pC(n0_per_cc, if_select = True)
    emittances = os_file.calculate_norm_rms_emittance_um(n0_per_cc, (2,3), if_select = True)
    if "rms"==gamma_spread_method:
        gamma, gamma_spread = os_file.raw_mean_rms_ene(if_select = True)
        #raw_mean_rms_ene() returns the ene average. +1 to obtain gamma
        gamma+=gamma
    elif "lfit"==gamma_spread_method:
        x, y = os_file.raw_hist_gamma(range_min=gamma_threshold, if_select = True)
        try:
            popt, _ = lf.FitLorentzian(x, y)
            gamma = popt[0]
            gamma_spread = abs(popt[2])*2.
        except:
            print("Warning: fit not found at t = {}.".format(os_file.time))
            gamma = 1.
            gamma_spread = 1.
    else: raise NotImplementedError("Method of {} for gamma spread of the beam has not been implemented yet!".format(gamma_spread_method))
    return os_file.time, q_pC, emittances, gamma, gamma_spread

def plot_beam_parameters_vs_t(path, species_name, n0_per_cc, code_name='osiris', gamma_threshold=1., start=0, count=0, stride=1, max_missing_file=0, gamma_spread_method="lfit"):
    t_array = list()
    q_pC_array = list()
    emittance_x_array = list()
    emittance_y_array = list()
    center_gamma_array = list()
    spread_gamma_array = list()
    missing_file = 0
    os_file = outfile.OutFile(code_name=code_name, path=path, field_name='raw',spec_name=species_name)
    for i in range(start, start+count*stride, stride):
        os_file.out_num=i
        try:
            os_file.open()
            #set missing_file = 0 if success
            missing_file = 0
        except IOError as err:
            missing_file = missing_file+1
            if missing_file>max_missing_file:
                print('Warning: unable to open file \'{0}\'. Iteration breaks. Exception message:\n{1}'.format(os_file.path_filename, err))
                break
            else:
                print('Warning! File {0} missing. Exception message:\n{1}'.format(os_file.path_filename, err))
                continue
        t, q_pC, emittances, center_gamma, spread_gamma = get_beam_parameters(os_file, n0_per_cc, gamma_threshold, gamma_spread_method=gamma_spread_method)
        os_file.close()
        t_array.append(t)
        q_pC_array.append(q_pC)
        emittance_x_array.append(emittances[0])
        emittance_y_array.append(emittances[1])
        center_gamma_array.append(center_gamma)
        spread_gamma_array.append(spread_gamma)
    h_f = plt.figure()
    h_ax = h_f.add_subplot(231)
    h_ax.plot(t_array, q_pC_array)
    plt.ylim()
    plt.xlabel('t [$\\omega_p$]')
    plt.ylabel('q [pC]')

    h_ax = h_f.add_subplot(232)
    h_ax.plot(t_array, emittance_x_array, 'k-', label='$\\epsilon_x$')
    h_ax.plot(t_array, emittance_y_array, 'r-', label='$\\epsilon_y$')
    plt.xlabel('t [$\\omega_p$]')
    plt.ylabel('$\epsilon$ [$\mu$m]')
    plt.legend()

    h_ax = h_f.add_subplot(233)
    h_ax.plot(t_array, center_gamma_array)
    plt.xlabel('t [$\\omega_p$]')
    plt.ylabel('$\\gamma$')

    h_ax = h_f.add_subplot(234)
    h_ax.plot(t_array, spread_gamma_array)
    plt.xlabel('t [$\\omega_p$]')
    plt.ylabel('$\\Delta \\gamma$')

    h_ax = h_f.add_subplot(235)
    h_ax.plot(t_array, np.divide(spread_gamma_array, center_gamma_array))
    plt.xlabel('t [$\\omega_p$]')
    plt.ylabel('$\\Delta \\gamma / \\gamma$')

    try:
        h_ax = h_f.add_subplot(236)
        os_file.plot_raw_hist_gamma(h_f, h_ax, if_select=True)
    except Exception as err:
        print("Exception message: {}".format(err))

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    plot_beam_parameters_vs_t(code_name='hipace', path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/X1_Shared_Pardis_Ming/50um800pC1.1e16', species_name='trailer', n0_per_cc=1.e16, gamma_threshold=1., start=0, count=99999, stride=20, max_missing_file=1, gamma_spread_method="rms")
    #plot_beam_parameters_vs_t(code_name='osiris', path='/home/zming/simulations/os2D/os_DRI3D19', species_name='plasma', n0_per_cc=4.0e16, gamma_threshold=10., start=33, count=999, stride=1)
    #plot_beam_parameters_vs_t(code_name='osiris', path='/home/zming/mnt/os_PT3D20', species_name='Ne_eK', n0_per_cc=2.43e18, gamma_threshold=50., start=10, count=999, stride=1)
    #plot_beam_parameters_vs_t(code_name='osiris', path='/home/zming/mnt/os_PT3D22', species_name='e', n0_per_cc=5.e17, gamma_threshold=50., start=8, count=999, stride=1)
    #plot_beam_parameters_vs_t(code_name='osiris', path='/home/zming/mnt/os_beam3D136', species_name='He_e', n0_per_cc=4.9e16, gamma_threshold=1., start=20, count=999, stride=1)
    #plot_beam_parameters_vs_t(code_name='hipace', path='/home/zming/simulations/os2D/hi_beam3D136', species_name='trailer', n0_per_cc=4.9e16, gamma_threshold=1., start=100, count=999, stride=10)

import outfile
import LorentzianFit as lf
import numpy as np
import matplotlib.pyplot as plt

# plot beam parameters vs time

def get_beam_parameters(outfile_object, n0_per_cc, gamma_threshold, gamma_spread_method="lfit"):
    ''' calculate beam parameters according to the raw output data
        return q_pC (charge in pico-coulomb), emittances (both transverse directions), popt[0] (center gamma), popt[2]*2. (absolute Full-Width-Half-Maximum spread)
    '''
    outfile_object.read_raw_x1()
    outfile_object.read_raw_x2()
    outfile_object.read_raw_x3()
    outfile_object.read_raw_p1()
    outfile_object.read_raw_p2()
    outfile_object.read_raw_p3()
    outfile_object.read_raw_q()
    outfile_object.read_raw_ene(ene_key_warning=False)
    outfile_object.select_raw_data(ene_low=gamma_threshold-1.)
    '''condlist = [outfile_object._raw_ene>(gamma_threshold-1.)]
    outfile_object._raw_q = np.select(condlist, [outfile_object._raw_q])'''
    if 0==len(outfile_object._raw_select_index):
        print("Warning: no particle is selected.")
        return outfile_object.time, 0.0, [0.0,0.0], 1.0, 00.5
    q_pC = outfile_object.calculate_q_pC(n0_per_cc, if_select = True)
    emittances = outfile_object.calculate_norm_rms_emittance_um(n0_per_cc, (2,3), if_select = True)
    if "rms"==gamma_spread_method:
        gamma, gamma_spread = outfile_object.raw_mean_rms_ene(if_select = True)
        #raw_mean_rms_ene() returns the ene average. +1 to obtain gamma
        gamma+=1
    elif "lfit"==gamma_spread_method:
        x, y = outfile_object.raw_hist_gamma(range_min=gamma_threshold, if_select = True)
        try:
            popt, _ = lf.FitLorentzian(x, y)
            gamma = popt[0]
            gamma_spread = abs(popt[2])*2.
        except:
            print("Warning: fit not found at t = {}.".format(outfile_object.time))
            gamma = 1.
            gamma_spread = 1.
    else: raise NotImplementedError("Method of {} for gamma spread of the beam has not been implemented yet!".format(gamma_spread_method))
    return outfile_object.time, q_pC, emittances, gamma, gamma_spread

def plot_beam_parameters_vs_t(path, species_name, n0_per_cc, code_name='osiris', gamma_threshold=1., use_num_list=False, start=0, count=0, stride=1, t_offset = 0., max_missing_file=0, gamma_spread_method="lfit", h_f=None, charge_abs=False, linestyle='-', label = None):
    '''
    max_missing_file is the maximum continue files missing allowed. HiPACE sometimes fails in dumping output files.
    '''
    t_array = list()
    q_pC_array = list()
    emittance_x_array = list()
    emittance_y_array = list()
    center_gamma_array = list()
    spread_gamma_array = list()
    missing_file = 0
    outfile_object = outfile.OutFile(code_name=code_name, path=path, field_name='raw',spec_name=species_name, use_num_list=use_num_list)
    for i in range(start, start+count*stride, stride):
        try: outfile_object.out_num=i
        except KeyError:
            print('Reaching the end of number list. Finishing...')
            break
        print('File number {}'.format(outfile_object.actual_num))
        try:
            outfile_object.open()
            #set missing_file = 0 if success
            missing_file = 0
        except IOError as err:
            missing_file += 1
            if missing_file>max_missing_file:
                print('Warning: Cannot open file number {}. Iteration breaks. Exception message:\n{}'.format(i, err))
                break
            else:
                print('Warning: File missing. Exception message:\n{}'.format(err))
                continue
        t, q_pC, emittances, center_gamma, spread_gamma = get_beam_parameters(outfile_object, n0_per_cc, gamma_threshold, gamma_spread_method=gamma_spread_method)
        outfile_object.close()
        t_array.append(t+t_offset)
        q_pC_array.append(q_pC)
        emittance_x_array.append(emittances[0])
        emittance_y_array.append(emittances[1])
        center_gamma_array.append(center_gamma)
        spread_gamma_array.append(spread_gamma)
    if h_f is None: h_f = plt.figure()
    h_ax = h_f.add_subplot(231)
    if charge_abs: q_pC_array = np.abs(q_pC_array)
    h_ax.plot(t_array, q_pC_array, linestyle = linestyle, label=label)
    plt.ylim()
    plt.xlabel('t [$\\omega_p$]')
    plt.ylabel('q [pC]')
    if label is not None: plt.legend()

    h_ax = h_f.add_subplot(232)
    h_ax.plot(t_array, emittance_x_array, color='k', linestyle = linestyle, label='$\\epsilon_x$')
    h_ax.plot(t_array, emittance_y_array, color='r', linestyle = linestyle, label='$\\epsilon_y$')
    plt.xlabel('t [$\\omega_p$]')
    plt.ylabel('$\epsilon$ [$\mu$m]')
    plt.legend()

    h_ax = h_f.add_subplot(233)
    h_ax.plot(t_array, center_gamma_array, linestyle = linestyle)
    plt.xlabel('t [$\\omega_p$]')
    plt.ylabel('$\\gamma$')

    h_ax = h_f.add_subplot(234)
    h_ax.plot(t_array, spread_gamma_array, linestyle = linestyle)
    plt.xlabel('t [$\\omega_p$]')
    plt.ylabel('$\\Delta \\gamma$')

    h_ax = h_f.add_subplot(235)
    h_ax.plot(t_array, np.divide(spread_gamma_array, center_gamma_array), linestyle = linestyle)
    plt.xlabel('t [$\\omega_p$]')
    plt.ylabel('$\\Delta \\gamma / \\gamma$')

    try:
        h_ax = h_f.add_subplot(236)
        outfile_object.plot_raw_hist_gamma(h_f, h_ax, if_select=True, linestyle = linestyle)
        if t_offset != 0.: h_ax.set_title('$t={0:.2f}$'.format(outfile_object.time+t_offset))
    except Exception as err:
        print("Exception message: {}".format(err))


if __name__ == '__main__':
    h_f = plt.figure()
    plot_beam_parameters_vs_t(code_name='osiris', path='/home/zming/mnt/JSCRATCH/50um600pC1.1e16thh', species_name='ramp', n0_per_cc=1.e16, gamma_threshold=1., start=3, count=775, stride=1, max_missing_file=0, gamma_spread_method="rms", charge_abs=True, h_f=h_f, label='OSIRIS')
    plot_beam_parameters_vs_t(code_name='hipace', path='/home/zming/mnt/JSCRATCH/50um600pC1.1e16thh', species_name='trailer', n0_per_cc=1.e16, gamma_threshold=1., use_num_list=True, start=0, count=119, stride=1, t_offset=0., max_missing_file=1, gamma_spread_method="rms", h_f=h_f, linestyle='--', label='HiPACE')
    #plot_beam_parameters_vs_t(code_name='osiris', path='/home/zming/simulations/os2D/os_DRI3D19', species_name='plasma', n0_per_cc=4.0e16, gamma_threshold=10., start=33, count=999, stride=1)
    #plot_beam_parameters_vs_t(code_name='osiris', path='/home/zming/mnt/os_PT3D20', species_name='Ne_eK', n0_per_cc=2.43e18, gamma_threshold=50., start=10, count=999, stride=1)
    #plot_beam_parameters_vs_t(code_name='osiris', path='/home/zming/mnt/os_PT3D22', species_name='e', n0_per_cc=5.e17, gamma_threshold=50., start=8, count=999, stride=1)
    #plot_beam_parameters_vs_t(code_name='osiris', path='/home/zming/mnt/JSCRATCH/os_beam3D142', species_name='He_e', n0_per_cc=4.9e16, gamma_threshold=1., start=15, count=999, stride=1)
    #plot_beam_parameters_vs_t(code_name='hipace', path='/home/zming/simulations/os2D/X1_Shared_Pardis_Ming/NegBox/50um600pC1.1e16', species_name='trailer', n0_per_cc=1.0e16, gamma_threshold=1., start=0, count=999, stride=20)
    plt.tight_layout()
    plt.show()

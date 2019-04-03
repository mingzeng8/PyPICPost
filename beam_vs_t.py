import outfile
import LorentzianFit as lf
import numpy as np
import matplotlib.pyplot as plt

# plot beam parameters vs time

def get_beam_parameters(os_file, n0_per_cc, dx, dy, dz, gamma_threshold):
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
    condlist = [os_file._raw_ene>(gamma_threshold-1.)]
    os_file._raw_q = np.select(condlist, [os_file._raw_q])
    q_pC = os_file.calculate_q_pC(n0_per_cc, dx, dy, dz)
    emittances = os_file.calculate_norm_rms_emittance_um(n0_per_cc, (2,3))
    x, y, _ = os_file.raw_hist_gamma(range_min=gamma_threshold)
    popt, _ = lf.FitLorentzian(x, y)
    return os_file.time, q_pC, emittances, popt[0], popt[2]*2.

def plot_beam_parameters_vs_t(path, species_name, n0_per_cc, dx, dy, dz, gamma_threshold=1., start=0, count=0, stride=1):
    t_array = list()
    q_pC_array = list()
    emittance_x_array = list()
    emittance_y_array = list()
    center_gamma_array = list()
    FWHM_gamma_array = list()
    os_file = outfile.OutFile(path=path, field_name='raw',spec_name=species_name)
    for i in range(start, start+count*stride, stride):
        os_file.out_num=i
        try:
            os_file.open()
        except IOError:
            print('Warning: unable to open file \'{0}\'. Iteration breaks.'.format(os_file.path_filename))
            break
        t, q_pC, emittances, center_gamma, FWHM_gamma = get_beam_parameters(os_file, n0_per_cc, dx, dy, dz, gamma_threshold)
        t_array.append(t)
        q_pC_array.append(q_pC)
        emittance_x_array.append(emittances[0])
        emittance_y_array.append(emittances[1])
        center_gamma_array.append(center_gamma)
        FWHM_gamma_array.append(FWHM_gamma)
    h_f = plt.figure()
    h_ax = h_f.add_subplot(221)
    h_ax.plot(t_array, q_pC_array)
    plt.ylim(top=0.)
    plt.xlabel('t [$\\omega_p$]')
    plt.ylabel('q [pC]')

    h_ax = h_f.add_subplot(222)
    h_ax.plot(t_array, emittance_x_array, 'k-', label='$\\epsilon_x$')
    h_ax.plot(t_array, emittance_y_array, 'r-', label='$\\epsilon_y$')
    plt.xlabel('t [$\\omega_p$]')
    plt.ylabel('$\epsilon$ [$\mu$m]')
    plt.legend()

    h_ax = h_f.add_subplot(223)
    h_ax.plot(t_array, center_gamma_array)
    plt.xlabel('t [$\\omega_p$]')
    plt.ylabel('$\\gamma$')

    h_ax = h_f.add_subplot(224)
    h_ax.plot(t_array, FWHM_gamma_array)
    plt.xlabel('t [$\\omega_p$]')
    plt.ylabel('$\\Delta \\gamma$')

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    plot_beam_parameters_vs_t('/home/zming/simulations/os2D/os_beam3D52', 'He_e', 4.9e16, 10./512, 8./2048, 8./256, 2., start=36, count=400, stride=1)

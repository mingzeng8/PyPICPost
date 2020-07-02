from outfile import *
from matplotlib.patches import Arc

def os_beam3D_plot(h_fig=None, h_ax=None, out_num=16, laser_field_name=None, trail_spec_name=None, if_psi=False, t_offset = None, **kwargs):
    if h_fig is None:
        h_fig = plt.figure()
    if h_ax is None:
        h_ax = h_fig.add_subplot(111)
    file1 = OutFile(path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/os_beam3D/os_beam3D146', field_name='charge', spec_name='e', out_num=out_num)
    file1.open(t_offset = t_offset)
    file1.read_data_slice()
    file1.close()
    file1.plot_data(h_fig, h_ax, cmap='gray', vmin=-5, if_colorbar=False, **kwargs)
    file1.spec_name='driver'
    file1.open(t_offset = t_offset)
    file1.read_data_slice()
    file1.close()
    file1.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_higher_range_transparent(plt.cm.hot), vmin=-18, if_colorbar=False, **kwargs)
    if if_psi:
        file1.field_name='psi'
        file1.open(t_offset = t_offset)
        file1.read_data_lineout()
        file1.close()
        file1.plot_data(h_fig, h_ax, ls='-',c='b', **kwargs)
    if trail_spec_name is not None:
        file1.field_name='charge'
        file1.spec_name=trail_spec_name
        file1.open(t_offset = t_offset)
        file1.read_data_slice()
        file1.close()
        file1.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_higher_range_transparent(), if_colorbar=False, **kwargs)
    if laser_field_name is not None:
        file1.field_name=laser_field_name
        file1.open(t_offset = t_offset)
        file1.read_data_slice()
        file1.close()
        #file1.data_profile2d()
        #file1.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_lower_range_transparent(plt.cm.Reds, transparency_transition_region=[0.15,0.4]), **kwargs)
        file1.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_middle_range_transparent(), if_colorbar=False, **kwargs)
    h_ax.set_xlim([0.,9.5])
    h_ax.set_ylim([-2.,2.])
    h_ax.set_aspect('equal','box')
    return h_fig, h_ax, file1.time

if __name__ == '__main__':
    #h_fig = plt.figure(figsize=[4., 4.])
    h_fig, h_axs = plt.subplots(4, 1, figsize=[4.5, 5.5], sharex=True)
    h_fig.subplots_adjust(hspace=0)
    #h_ax = h_fig.add_subplot(211)
    _, _, t = os_beam3D_plot(h_fig, h_axs[0], out_num=16, laser_field_name='e3', if_psi=True, if_z2zeta=True, t_offset = -10.)
    h_axs[0].text(6.6,1.4,'$\\omega_p t = {:.1f}$'.format(t+10.))

    delta = 0.7
    theta = 75.
    laser_start0 = [23.5-t, -4.]
    laser_start = [24.-t, -2.]
    laser_start[0] = laser_start0[0]+(laser_start[1]-laser_start0[1])/np.tan(theta/180.*np.pi)
    #print(laser_start)
    arrow_dz = .8
    h_axs[0].arrow(laser_start[0], laser_start[1], arrow_dz, arrow_dz*np.tan(theta/180.*np.pi), width=.05, color='g')
    arrow_middle_dz = 0.3
    arrow_r_dz = 1.
    h_axs[0].arrow(laser_start[0]+arrow_middle_dz, laser_start[1]+arrow_middle_dz*np.tan(theta/180.*np.pi), arrow_r_dz, -arrow_r_dz/np.tan(theta/180.*np.pi), width=.05, color='g')
    h_axs[0].text(25.-t, 1.5, '$z\'$')
    h_axs[0].text(25.7-t, -1.4, '$r\'$')
    h_axs[0].add_patch(Arc(laser_start, 1., 1., 0., 0., theta, color='g'))
    h_axs[0].text(24.6-t, -1.8, '$\\theta$')
    h_axs[0].set_title(None)
    h_axs[0].set_xlabel(None)
    h_axs[0].set_ylabel('$k_p x$, $\\psi$')

    # Draw the defocusing phase    
    h_axs[0].add_patch(Arc((21.35-t, 0.), .9, .5, 0., 0., 360., color='y', ls='--', lw=2.))
    #h_axs[0].text(21.0, .35, 'Defocusing phase')

    #h_ax = h_fig.add_subplot(212)
    _, _, t = os_beam3D_plot(h_fig, h_axs[1], out_num=17, laser_field_name='e3', trail_spec_name='He_e', if_z2zeta=True, t_offset = -10.)
    h_axs[1].set_title(None)
    h_axs[1].set_xlabel('$k_p \\zeta$')
    h_axs[1].set_ylabel('$k_p x$')
    h_axs[1].text(6.6,1.4,'$\\omega_p t = {:.1f}$'.format(t+10.))

    _, _, t = os_beam3D_plot(h_fig, h_axs[2], out_num=18, laser_field_name='e3', trail_spec_name='He_e', if_z2zeta=True, t_offset = -10.)
    h_axs[2].set_title(None)
    h_axs[2].set_xlabel('$k_p \\zeta$')
    h_axs[2].set_ylabel('$k_p x$')
    h_axs[2].text(6.6,1.4,'$\\omega_p t = {:.1f}$'.format(t+10.))

    _, _, t = os_beam3D_plot(h_fig, h_axs[3], out_num=30, trail_spec_name='He_e', if_z2zeta=True, t_offset = -10.)
    h_axs[3].set_title(None)
    h_axs[3].set_xlabel('$k_p \\zeta$')
    h_axs[3].set_ylabel('$k_p x$')
    h_axs[3].text(6.6,1.4,'$\\omega_p t = {:.1f}$'.format(t+10.))
    #plt.tight_layout()
    plt.show()

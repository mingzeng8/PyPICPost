from outfile import *
from matplotlib.patches import Arc

def os_beam3D_plot(h_fig=None, h_ax=None, out_num=16, laser_field_name=None, trail_spec_name=None, **kwargs):
    if h_fig is None:
        h_fig = plt.figure()
    if h_ax is None:
        h_ax = h_fig.add_subplot(111)
    file1 = OutFile(path='/home/zming/mnt/JSCRATCH/os_beam3D/os_beam3D146', field_name='charge', spec_name='e', out_num=out_num)
    file1.open()
    file1.read_data_slice()
    file1.close()
    file1.plot_data(h_fig, h_ax, cmap='gray', vmin=-5, if_colorbar=False, **kwargs)
    file1.spec_name='driver'
    file1.open()
    file1.read_data_slice()
    file1.close()
    file1.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_higher_range_transparent(plt.cm.hot), vmin=-18, if_colorbar=False, **kwargs)
    if trail_spec_name is not None:
        file1.field_name='charge'
        file1.spec_name=trail_spec_name
        file1.open()
        file1.read_data_slice()
        file1.close()
        file1.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_higher_range_transparent(), if_colorbar=False, **kwargs)
    if laser_field_name is not None:
        file1.field_name=laser_field_name
        file1.open()
        file1.read_data_slice()
        file1.close()
        #file1.data_profile2d()
        #file1.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_lower_range_transparent(plt.cm.Reds, transparency_transition_region=[0.15,0.4]), **kwargs)
        file1.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_middle_range_transparent(), if_colorbar=False, **kwargs)
    h_ax.set_ylim([-2.,2.])
    h_ax.set_aspect('equal','box')
    return h_fig, h_ax

if __name__ == '__main__':
    h_fig = plt.figure(figsize=[4., 4.])
    h_ax = h_fig.add_subplot(211)
    os_beam3D_plot(h_fig, h_ax, out_num=16, laser_field_name='e3')

    delta = 0.7
    theta = 75.
    laser_start0 = [23.5, -4.]
    laser_start = [24., -2.]
    laser_start[0] = laser_start0[0]+(laser_start[1]-laser_start0[1])/np.tan(theta/180.*np.pi)
    #print(laser_start)
    arrow_dz = .8
    h_ax.arrow(laser_start[0], laser_start[1], arrow_dz, arrow_dz*np.tan(theta/180.*np.pi), width=.05, color='g')
    arrow_middle_dz = 0.3
    arrow_r_dz = 1.
    h_ax.arrow(laser_start[0]+arrow_middle_dz, laser_start[1]+arrow_middle_dz*np.tan(theta/180.*np.pi), arrow_r_dz, -arrow_r_dz/np.tan(theta/180.*np.pi), width=.05, color='g')
    h_ax.text(25., 1.5, '$z\'$')
    h_ax.text(25.7, -1.4, '$r\'$')
    h_ax.add_patch(Arc(laser_start, 1., 1., 0., 0., theta, color='g'))
    h_ax.text(24.6, -1.8, '$\\theta$')
    h_ax.set_title(None)
    h_ax.set_xlabel(None)
    h_ax.set_ylabel('$k_p x$')

    h_ax = h_fig.add_subplot(212)
    os_beam3D_plot(h_fig, h_ax, out_num=18, laser_field_name='e3', trail_spec_name='He_e')
    h_ax.set_title(None)
    h_ax.set_xlabel('$k_p z$')
    h_ax.set_ylabel('$k_p x$')
    plt.tight_layout()
    plt.show()

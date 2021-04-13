# Examples and tests
if __name__ == '__main__':
    out_num=1
    def tmp_3D():
        h_fig = plt.figure(figsize=(5,5))
        file1 = OutFile(path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/os_beam3D/os_beam3D146', field_name='e1', out_num=52)
        h_ax = h_fig.add_subplot(111)
        file1.open()
        file1.read_data_slice()
        file1.plot_data(h_fig, h_ax)#, cmap=my_cmap.cmap_higher_range_transparent())
        h_ax.set_aspect('equal','box')
        file1.read_data()
        file1.close()
        E=file1._data
        print('E_max = {}, E_min = {}'.format(E.max(), E.min()))
        plt.tight_layout()
        plt.show()
    def X1plot():
        h_fig = plt.figure(figsize=(5,5))
        file1 = OutFile(code_name='osiris',path='/home/zming/mnt/JSCRATCH/X1/scan_2019_12_04/He_Ar_3.0/1.0', field_name='charge', spec_name='plasma', out_num=18)
        h_ax = h_fig.add_subplot(111)
        file1.open()
        file1.read_data_slice(dir=1)
        file1.plot_data(h_fig, h_ax, cmap='gray', vmin=-0.001, if_colorbar=False)
        h_ax.set_aspect('equal','box')
        file1.close()
        file1.spec_name = 'driver'
        file1.open()
        file1.read_data_slice(dir=1)
        file1.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_higher_range_transparent(plt.cm.hot), vmin=-0.002, if_colorbar=False)
        file1.close()
        h_ax.set_title(None)
        h_ax.set_xlabel('z [$\mu$m]')
        h_ax.set_ylabel('y [$\mu$m]')
        plt.tight_layout()

        h_fig = plt.figure(figsize=(5,5))
        file1.out_num = 64
        file1.spec_name = 'ramp'
        h_ax = h_fig.add_subplot(111)
        file1.open()
        file1.read_data_slice(dir=1)
        ramp_data = file1._data
        file1.close()
        file1.spec_name = 'plasma'
        file1.open()
        file1.read_data_slice(dir=1)
        plasma_data = file1._data
        #file1._data = np.select([ramp_data>=(-0.0003541), ramp_data<(-0.0003541)],[plasma_data+ramp_data, -0.0003541])
        file1._data += ramp_data
        file1.plot_data(h_fig, h_ax, cmap=plt.cm.gray, vmin=-0.001, if_colorbar=False)
        file1._data = np.select([ramp_data>=(-0.0006), ramp_data<(-0.0006)],[0., ramp_data])
        file1.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_higher_range_transparent(plt.cm.jet), vmin=-0.002, if_colorbar=False)
        h_ax.set_aspect('equal','box')
        file1.close()
        file1.spec_name = 'driver'
        file1.open()
        file1.read_data_slice(dir=1)
        file1.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_higher_range_transparent(plt.cm.hot, transparency_transition_region=[0.8,0.9]), vmin=-0.002, if_colorbar=False)
        file1.close()
        h_ax.set_title(None)
        h_ax.set_xlabel('z [$\mu$m]')
        h_ax.set_ylabel('y [$\mu$m]')
        plt.tight_layout()

        h_fig = plt.figure(figsize=(5,5))
        file1.out_num = 133
        file1.spec_name = 'plasma'
        h_ax = h_fig.add_subplot(111)
        file1.open()
        file1.read_data_slice(dir=1)
        file1.plot_data(h_fig, h_ax, cmap=plt.cm.gray, vmin=-0.001, if_colorbar=False)
        h_ax.set_aspect('equal','box')
        file1.close()
        file1.spec_name = 'ramp'
        file1.open()
        file1.read_data_slice(dir=1)
        file1.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_higher_range_transparent(plt.cm.jet, transparency_transition_region=[0.5,0.6]), vmin=-0.002, if_colorbar=False)
        h_ax.set_aspect('equal','box')
        file1.close()
        file1.spec_name = 'driver'
        file1.open()
        file1.read_data_slice(dir=1)
        file1.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_higher_range_transparent(plt.cm.hot), vmin=-0.002, if_colorbar=False)
        file1.close()
        h_ax.set_title(None)
        h_ax.set_xlabel('z [$\mu$m]')
        h_ax.set_ylabel('y [$\mu$m]')
        plt.tight_layout()

        plt.show()
    def cyl_m_2D():
        h_fig = plt.figure(figsize=(16.5,11))
        file1 = OutFile(path='/home/zming/simulations/os2D/os_laser3DQ3',field_name='charge',spec_name='plasma',out_num=out_num,cyl_m_num=0,cyl_m_re_im='re')
        h_ax = h_fig.add_subplot(221)
        file1.open()
        file1.read_data()
        file1.plot_data(h_fig, h_ax)#, if_log_colorbar=True)
        file1.close()
        file1.field_name='e3_cyl_m'
        h_ax = h_fig.add_subplot(222)
        file1.open()
        file1.read_data()
        file1.plot_data(h_fig, h_ax)#, if_log_colorbar=True)#, vmax=1.)
        file1.close()
        file1.cyl_m_num=1
        h_ax = h_fig.add_subplot(223)
        file1.open()
        file1.read_data()
        file1.plot_data(h_fig, h_ax)#, vmax=3., vmin=-3.)
        file1.close()
        #file1.cyl_m_re_im='im'
        h_ax = h_fig.add_subplot(224)
        file1.open()
        file1.read_data_project2d(dir = 0, if_abs = True)
        #file1.plot_data(h_fig, h_ax)
        file1.fit_for_W(h_fig, h_ax)
        print(file1._W)
        file1.close()
        plt.show()
    def tmp_hipace3D():
        file1 = OutFile(code_name='hipace', path='/home/zming/simulations/os2D/hi_beam3D105',field_name='raw',spec_name='trailer',out_num=990)
        file1.open(filename='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/X1_Shared_Pardis_Ming/NegBox/50um600pC1.1e16/RAW-beam-driver-000100.h5')
        file1.read_raw_x1()
        file1.read_raw_p1()
        file1.close()
        h_fig = plt.figure()
        h_ax = h_fig.add_subplot(111)
        h_ax.scatter(file1._raw_x1, file1._raw_p1, s=1, marker='.', c='r')
        file1.open(filename='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/X1_Shared_Pardis_Ming/NegBox/50um600pC1.1e16/RAW-plasma-000100.h5')
        file1.read_raw_x1()
        file1.read_raw_p1()
        file1.close()
        h_ax.scatter(file1._raw_x1, file1._raw_p1, s=1, marker='.', c='b')
        plt.show(block=True)
    def scatter():
        h_fig = plt.figure()
        h_ax1 = h_fig.add_subplot(211)
        h_ax2 = h_fig.add_subplot(212)
        file1 = OutFile(path='/home/zming/mnt/JSCRATCH/os_beam3D/os_beam3D143',field_name='raw',spec_name='driver',out_num=40)
        file1.open()
        file1.read_raw_x1()
        file1.read_raw_x2()
        file1.read_raw_x3()
        file1.close()
        h_ax1.scatter(file1._raw_x1, file1._raw_x2, s=1, marker='.', c='r')
        h_ax1.set_xlabel('$x_1$')
        h_ax1.set_ylabel('$x_2$')
        h_ax1.set_title('$t={}$'.format(file1.time))
        h_ax2.scatter(file1._raw_x1, file1._raw_x3, s=1, marker='.', c='r')
        h_ax2.set_xlabel('$x_1$')
        h_ax2.set_ylabel('$x_3$')
        h_ax2.set_title('$t={}$'.format(file1.time))

        file1.spec_name='He_e'        
        file1.open()
        file1.read_raw_x1()
        file1.read_raw_x2()
        file1.read_raw_x3()
        file1.close()
        h_ax1.scatter(file1._raw_x1, file1._raw_x2, s=1, marker='.', c='b')
        h_ax1.set_xlabel('$x_1$')
        h_ax1.set_ylabel('$x_2$')
        h_ax1.set_title('$t={}$'.format(file1.time))
        h_ax2.scatter(file1._raw_x1, file1._raw_x3, s=1, marker='.', c='b')
        h_ax2.set_xlabel('$x_1$')
        h_ax2.set_ylabel('$x_3$')
        h_ax2.set_title('$t={}$'.format(file1.time))
        plt.tight_layout()
        plt.show(block=True)
    def select_tag():
        file1 = OutFile(path='/home/zming/mnt/X1_Shared_Pardis_Ming/NegBox/50um600pC1.1e16h',field_name='raw',spec_name='ramp',out_num=20)
        file1.open()
        #print(file1.fileid.keys())
        file1.read_raw_q()
        file1.read_raw_x1()
        file1.read_raw_x2()
        file1.read_raw_x3()
        file1.read_raw_p1()
        file1.read_raw_p2()
        file1.read_raw_p3()
        h_fig = plt.figure()
        #file1.select_raw_data(x2_up=0.05, r_up=0.17)
        h_ax = h_fig.add_subplot(221)
        file1.plot_raw_hist2D(h_fig, h_ax, dims='p1x1', num_bins=128, range=None, cmap=None, if_select = True, if_colorbar=True, colorbar_orientation="vertical", if_log_colorbar=False)
        h_ax = h_fig.add_subplot(222)
        file1.plot_raw_hist2D(h_fig, h_ax, dims='p2x2', num_bins=128, range=None, cmap=None, if_select = True, if_colorbar=True, colorbar_orientation="vertical", if_log_colorbar=False)
        print('{} pC, {} um'.format(file1.calculate_q_pC(n0_per_cc = 4.9e16, if_select = True), file1.calculate_norm_rms_emittance_um(n0_per_cc = 4.9e16, directions=(2,3), if_select = True)))
        h_ax = h_fig.add_subplot(223)
        #h_ax.scatter(file1._raw_x1[file1._raw_select_index], file1._raw_x2[file1._raw_select_index], s=1, marker='.', c='r')
        h_ax = h_fig.add_subplot(224)
        #h_ax.scatter(file1._raw_x1[file1._raw_select_index], file1._raw_x3[file1._raw_select_index], s=1, marker='.', c='r')
        #file1.read_raw_tag()
        #file1.save_tag_file(if_select = True)
        file1.close()
        plt.show(block=True)
    def tracking():
        file1 = OutFile(code_name='osiris',path='/home/zming/mnt/JDATA/os_beam3D/os_beam3D110',field_name='tracks',average='',spec_name='He_e')
        file1.open()
        file1.plot_tracks(x_quant = b't', y_quant = b'x1')
        file1.plot_tracks(x_quant = b't', y_quant = b'x2')
        file1.plot_tracks(x_quant = b't', y_quant = b'x3')
        #file1.plot_tracks(mwin_v = 1., mwin_t_offset = 32.3)
        file1.close()
        plt.show()
    def oblique_laser_2D():
        h_fig = plt.figure(figsize=(5,5))
        file1 = OutFile(code_name='osiris',path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/os_oblique_laser_test2D3', field_name='e3', spec_name='plasma', out_num=12)
        h_ax = h_fig.add_subplot(111)
        file1.open()
        file1.read_data()
        #e1=file1._data
        file1.close()
        '''file1.field_name='e2'
        file1.open()
        file1.read_data()
        # get total e field
        file1._data = np.sqrt(np.square(file1._data)+np.square(e1))*np.sign(e1)
        file1.close()'''
        print('Maximum of e3 is {}'.format(file1._data.max()))
        file1.plot_data(h_fig, h_ax)#, vmin=-0.001)#, cmap=my_cmap.cmap_higher_range_transparent())
        angle = 30.
        z = np.array([0., 0.])
        x = np.array([-5., 4.9])
        z[1] = z[0]+(x[1]-x[0])*np.tan(angle/180.*np.pi)
        plt.plot(z,x,'r--')
        h_ax.set_aspect('equal','box')

        plt.tight_layout()
        plt.show()
    def oblique_laser_3D():
        h_fig = plt.figure(figsize=(5,5))
        file1 = OutFile(code_name='osiris',path='/home/zming/mnt/JSCRATCH/os_beam3D/os_beam3D148', field_name='e3', out_num=15)
        h_ax = h_fig.add_subplot(111)
        file1.open()
        file1.read_data_slice()
        #e1=file1._data
        file1.close()
        '''file1.field_name='e2'
        file1.open()
        file1.read_data()
        # get total e field
        file1._data = np.sqrt(np.square(file1._data)+np.square(e1))*np.sign(e1)
        file1.close()'''
        print('Maximum of e3 is {}'.format(file1._data.max()))
        file1.plot_data(h_fig, h_ax)#, vmin=-0.001)#, cmap=my_cmap.cmap_higher_range_transparent())
        angle = 15.
        z = np.array([23.5, 0.])
        x = np.array([-4., 3.9])
        z[1] = z[0]+(x[1]-x[0])*np.tan(angle/180.*np.pi)
        plt.plot(z,x,'r--')
        h_ax.set_aspect('equal','box')

        plt.tight_layout()
        plt.show()
    def beam_measure_beam3D():
        file1 = OutFile(code_name='osiris',path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/os_beam3D/os_beam3D264',field_name='raw',spec_name='He_e',use_num_list=True,out_num=50)
        file1.open()
        file1.read_raw_q()
        file1.read_raw_p1()
        file1.read_raw_x2()
        file1.read_raw_p2()
        file1.read_raw_x3()
        file1.read_raw_p3()
        file1.read_raw_ene()
        #file1.select_raw_data(ene_low=10.)
        print('gamma_min = {}'.format(file1._raw_ene.min()+1.))
        print('gamma_max = {}'.format(file1._raw_ene.max()+1.))
        ene_avg, ene_rms_spread = file1.raw_mean_rms_ene(if_select=True)
        emittance, twiss_parameters =file1.calculate_norm_rms_emittance_um(4.9e16, directions=(2,3), if_select=True)
        h_fig, h_ax = file1.plot_raw_hist2D(dims='p1x1', if_select=True, if_log_colorbar=True)
        #file1.data_center_of_mass2d(weigh_threshold=1e1)
        #h_ax2 = h_ax.twinx()
        #file1.data2D_slice_spread(weigh_threshold=1e-1)
        #file1.plot_data(h_fig, h_ax2, if_title=False, c='r')
        #h_axt = h_ax.twinx()
        #z, I = file1.raw_hist_current(n0_per_cc=4.9e16, num_bins=32, if_select=True)
        #h_axt.plot(z,I, c='b')
        print('t = {}'.format(file1.time))
        print('Q = {} pC'.format(file1.calculate_q_pC(4.9e16, if_select=True)))
        print('E = {} +- {} MeV'.format(ene_avg*0.511, ene_rms_spread*0.511))
        print('epsilon_x_norm = {} um, epsilon_y_norm = {} um'.format(emittance[0], emittance[1]))
        file1.close()
        plt.show()
    def beam_measure_hi():
        file1 = OutFile(code_name='hipace',path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/hi_beam3D/newhi_beam3D247rr',field_name='raw',spec_name='trailer',use_num_list=True,out_num=200)
        file1.open()
        file1.read_raw_q()
        file1.read_raw_p1()
        file1.read_raw_x2()
        file1.read_raw_p2()
        file1.read_raw_x3()
        file1.read_raw_p3()
        file1.read_raw_ene()
        file1.select_raw_data(r_up=1.2)
        ene_avg, ene_rms_spread = file1.raw_mean_rms_ene(if_select=True)
        emittance, twiss_parameters =file1.calculate_norm_rms_emittance_um(4.9e16, directions=(2,3), if_select=True)
        h_fig, h_ax = file1.plot_raw_hist2D(dims='p2x2', num_bins=[128,128], if_select=True, if_log_colorbar=False)
        h_ax2 = h_ax.twinx()
        file1.data2D_slice_spread(weigh_threshold=1e2, relative_spread=True)
        file1.plot_data(h_fig, h_ax2, if_title=False, c='r')
        plt.figure()
        z, I = file1.raw_hist_current(n0_per_cc=5.03e15, num_bins=32)
        plt.plot(z,I)
        plt.xlabel('$k_p z$')
        plt.ylabel('$I$ [A]')
        print('t = {}'.format(file1.time))
        print('Q = {} pC'.format(file1.calculate_q_pC(4.9e16, if_select=True)))
        print('E = {} +- {} MeV'.format(ene_avg*0.511, ene_rms_spread*0.511))
        print('epsilon_x_norm = {} um, epsilon_y_norm = {} um'.format(emittance[0], emittance[1]))
        file1.close()
        plt.show()
    def hi_beam3D():
        h_fig, h_axs = plt.subplots(1, 2, figsize=[5., 2.5])
        file1 = OutFile(code_name='hipace',path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/hi_beam3D/newhi_beam3D161',field_name='raw',spec_name='trailer',use_num_list=True,out_num=1735)
        file1.open()
        file1.read_raw_q()
        file1.read_raw_p1()
        file1.read_raw_x2()
        file1.read_raw_p2()
        file1.read_raw_x3()
        file1.read_raw_p3()
        file1.read_raw_ene()
        file1.select_raw_data(p1_low=10.)
        ene_avg, ene_rms_spread = file1.raw_mean_rms_ene(if_select=True)
        emittance, twiss_parameters =file1.calculate_norm_rms_emittance_um(4.9e16, directions=(2,3), if_select=True)
        file1.plot_raw_hist2D(h_fig, h_axs[0], dims='p1x1', num_bins=[256,128], if_select=True, if_colorbar=False, cmap=my_cmap.cmap_lower_range_transparent(original_cmap=plt.get_cmap('Reds')))
        h_axs[0].set_xlabel('$k_p\zeta$')
        h_axs[0].set_title(None)
        h_axt = h_axs[0].twinx()
        #file1.data2D_slice_spread(weigh_threshold=1e1)
        #file1.plot_data(h_fig, h_ax2, if_title=False, c='r')
        #plt.figure()
        z, I = file1.raw_hist_current(n0_per_cc=4.9e16, num_bins=32, if_select=True)
        h_axt.plot(z,I, c='b')
        #plt.xlabel('$k_p z$')
        #h_axt.set_ylabel('$I$ [A]')
        h_axt.minorticks_on()
        print('Time = {}'.format(file1.time))
        print('Q = {} pC'.format(file1.calculate_q_pC(4.9e16, if_select=True)))
        #print('E = {} +- {} MeV'.format(ene_avg*0.511, ene_rms_spread*0.511))
        print('epsilon_x_norm = {} um, epsilon_y_norm = {} um'.format(emittance[0], emittance[1]))
        file1.close()

        file1 = OutFile(code_name='hipace',path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/hi_beam3D/newhi_beam3D259',field_name='raw',spec_name='trailer',use_num_list=True,out_num=120)
        file1.open()
        file1.read_raw_q()
        file1.read_raw_p1()
        file1.read_raw_x2()
        file1.read_raw_p2()
        file1.read_raw_x3()
        file1.read_raw_p3()
        file1.read_raw_ene()
        file1.select_raw_data(p1_low=10.,r_up=0.4)
        ene_avg, ene_rms_spread = file1.raw_mean_rms_ene(if_select=True)
        emittance, twiss_parameters =file1.calculate_norm_rms_emittance_um(4.9e16, directions=(2,3), if_select=True)
        file1.plot_raw_hist2D(h_fig, h_axs[1], dims='p1x1', num_bins=[256,128], if_select=True, if_colorbar=False, cmap=my_cmap.cmap_lower_range_transparent(original_cmap=plt.get_cmap('Reds')))
        h_axs[1].set_xlabel('$k_p\zeta$')
        h_axs[1].set_ylabel(None)
        h_axs[1].set_title(None)
        h_axt = h_axs[1].twinx()
        #file1.data2D_slice_spread(weigh_threshold=1e1)
        #file1.plot_data(h_fig, h_ax2, if_title=False, c='r')
        #plt.figure()
        z, I = file1.raw_hist_current(n0_per_cc=4.9e16, num_bins=32, if_select=True)
        h_axt.plot(z,I, c='b')
        #plt.xlabel('$k_p z$')
        h_axt.set_ylabel('$I$ [A]')
        h_axt.minorticks_on()
        print('Time = {}'.format(file1.time))
        print('Q = {} pC'.format(file1.calculate_q_pC(4.9e16, if_select=True)))
        #print('E = {} +- {} MeV'.format(ene_avg*0.511, ene_rms_spread*0.511))
        print('epsilon_x_norm = {} um, epsilon_y_norm = {} um'.format(emittance[0], emittance[1]))
        file1.close()
        plt.tight_layout()
        plt.show()
    def beam_measure():
        file1 = OutFile(code_name='osiris',path='/home/zming/mnt/JSCRATCH/X1/scan_2020_6_10/tailoring7',field_name='raw',spec_name='driver',out_num=150)
        file1.open()
        file1.read_raw_q()
        file1.read_raw_p1()
        file1.read_raw_x2()
        file1.read_raw_p2()
        file1.read_raw_x3()
        file1.read_raw_p3()
        file1.read_raw_ene()
        print('gamma_min = {}'.format(file1._raw_ene.min()+1.))
        print('gamma_max = {}'.format(file1._raw_ene.max()+1.))
        #file1.select_raw_data(ene_up=2191.)
        ene_avg, ene_rms_spread = file1.raw_mean_rms_ene(if_select=True)
        emittance, twiss_parameters =file1.calculate_norm_rms_emittance_um(2.824e19, directions=(2,3), if_select=True)
        print('Q = {} pC'.format(file1.calculate_q_pC(2.824e19, if_select=True)))
        print('E = {} +- {} MeV'.format(ene_avg*0.511, ene_rms_spread*0.511))
        print('epsilon_x_norm = {} um, epsilon_y_norm = {} um'.format(emittance[0], emittance[1]))
        h_fig, h_ax = file1.plot_raw_hist2D(dims='p1x1', if_select=True, if_log_colorbar=True, if_colorbar=False)
        '''h_ax2 = h_ax.twinx()
        file1.data2D_slice_spread(weigh_threshold=1.e-1)
        file1.plot_data(h_fig, h_ax2, if_title=False)
        file1.plot_data(h_fig, h_ax2, c='r', if_title=False)'''
        file1.close()
        plt.show()
    def test_os():
        file1 = OutFile(code_name='osiris',path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/os_PDMU/os_PDMU2',field_name='raw',spec_name='mu',out_num=17)
        file1.open(filename = '/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/qp_hi_compare/qp_cinj/driver.h5')
        file1.read_raw_q()
        file1.read_raw_x1()
        file1.read_raw_x2()
        file1.read_raw_x3()
        #file1.read_raw_p1()
        #file1.read_raw_p2()
        #file1.read_raw_p3()
        file1.read_raw_gamma()
        file1.read_raw_ene()
        print('gamma_min = {}'.format(file1._raw_gamma.min()))
        print('gamma_max = {}'.format(file1._raw_gamma.max()))
        print('gamma_minux_pz_min = {}'.format(np.min(file1._raw_gamma-file1._raw_p1)))
        #file1.select_raw_data(ene_low=2.)
        ene_avg, ene_rms_spread = file1.raw_mean_rms_ene(if_select=True)
        emittance, twiss_parameters =file1.calculate_norm_rms_emittance_um(5.0334e15, directions=(2,3), if_select=True)
        print('Q = {} pC'.format(file1.calculate_q_pC(5.0334e15, if_select=True)))
        print('E = {} +- {} MeV'.format(ene_avg*0.511, ene_rms_spread*0.511))
        print('epsilon_x_norm = {} um, epsilon_y_norm = {} um'.format(emittance[0], emittance[1]))
        h_fig, h_ax = file1.plot_raw_hist2D(dims='p1x1', num_bins=[128,256], if_reread=False, if_select=True, if_log_colorbar=True, if_colorbar=True)#,range=[[19569.4, 19569.6], [5.,11.]])
        print(file1._raw_x1.max())
        file1.close()
        plt.show()
    def test_hi():
        file1 = OutFile(code_name='hipace',path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/hi_qp_compare/hi3',field_name='raw',spec_name='trailer', use_num_list=True, out_num=584)
        file1.open()
        print(file1._cell_size)
        file1.read_raw_q()
        file1.read_raw_p1()
        file1.read_raw_x2()
        file1.read_raw_p2()
        file1.read_raw_x3()
        file1.read_raw_p3()
        file1.read_raw_ene()
        print('gamma_min = {}'.format(file1._raw_ene.min()+1.))
        print('gamma_max = {}'.format(file1._raw_ene.max()+1.))
        #file1.select_raw_data(ene_up=2191.)
        ene_avg, ene_rms_spread = file1.raw_mean_rms_ene(if_select=True)
        emittance, twiss_parameters =file1.calculate_norm_rms_emittance_um(5.03e15, directions=(2,3), if_select=True)
        print('Q = {} pC'.format(file1.calculate_q_pC(5.03e15, if_select=True)))
        print('E = {} +- {} MeV'.format(ene_avg*0.511, ene_rms_spread*0.511))
        print('epsilon_x_norm = {} um, epsilon_y_norm = {} um'.format(emittance[0], emittance[1]))
        h_fig, h_ax = file1.plot_raw_hist2D(dims='p1x1', if_select=True, if_log_colorbar=True, if_colorbar=True, vmin=.05, vmax=3.e2)
        file1.close()
    def test_qp():
        file1 = OutFile(code_name='quickpic',path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/qp_hi_compare/qp_cinj',field_name='raw',spec_name='Beam0001', out_num=0)
        file1.open(cell_size_qp_raw = np.array([12/512, 12/512, 13/512]))
        file1.read_raw_q()
        file1.read_raw_p1()
        file1.read_raw_x2()
        file1.read_raw_p2()
        file1.read_raw_x3()
        file1.read_raw_p3()
        file1.read_raw_ene()
        print('gamma_min = {}'.format(file1._raw_ene.min()+1.))
        print('gamma_max = {}'.format(file1._raw_ene.max()+1.))
        #file1.select_raw_data(ene_up=2191.)
        ene_avg, ene_rms_spread = file1.raw_mean_rms_ene(if_select=True)
        emittance, twiss_parameters =file1.calculate_norm_rms_emittance_um(5.0334e15, directions=(2,3), if_select=True)
        print('Q = {} pC'.format(file1.calculate_q_pC(5.0334e15, if_select=True)))
        print('E = {} +- {} MeV'.format(ene_avg*0.511, ene_rms_spread*0.511))
        print('epsilon_x_norm = {} um, epsilon_y_norm = {} um'.format(emittance[0], emittance[1]))
        h_fig, h_ax = file1.plot_raw_hist2D(dims='p1x1', if_select=True, if_log_colorbar=True, if_colorbar=True)
        print(file1._raw_x1.min())
        file1.close()
    def test_profile():
        dir=2
        file1 = OutFile()
        file1.open(filename='/home/zming/PyHome/scripts_X1/Profile_plotter/myprofile.h5')
        file1.read_data_slice(dir=dir)
        file1.plot_data(vmin=0.5, vmax=1.5)
        file1.close()
        file1 = OutFile(path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/os_prof_file_test',field_name='charge',spec_name='e',out_num=7)
        file1.open()
        file1.read_data_slice(dir=dir)
        file1._data=np.abs(file1._data)
        file1.plot_data(vmin=0.5, vmax=1.5)
        file1.close()
        plt.show()

    def test_fb_lineout():
        from scipy.fft import fft, ifft
        file1 = OutFile(code_name = 'fbpic', path='/home/ming/mnt/CNG12/FB/DCII/diag_k=79_a0=3.3_w0=3.7_dphase=0.0_2mmdope', field_name='e2', use_num_list=True, out_num=0)
        file1.open()
        file1.read_data_slice()
        #file1.plot_data()
        file1.data_lineout2d(dir = 0, pos = None)
        h_f, h_ax = file1.plot_data(linewidth=0.2, c='k')
        file1._data = fft(file1._data)
        l = len(file1._data)
        r = 0.02
        ind1 = int(r*l)
        ind2 = l-ind1
        f2 = np.zeros(l)
        # f2 has the higher frequency
        f2[ind1:ind2] = file1._data[ind1:ind2]
        # file1._data has the lower frequency
        file1._data[ind1:ind2] = np.zeros(ind2-ind1)
        #plt.plot(file1._data)
        #plt.plot(f2)
        file1._data = ifft(file1._data)
        file1.plot_data(h_fig=h_f, h_ax=h_ax, c='r')
        file1._data = ifft(f2)
        file1.plot_data(h_fig=h_f, h_ax=h_ax, c='b')
        file1.spec_name = 'Ne_e'
        file1.read_raw('q')
        print('Q = {} pC'.format(file1.calculate_q_pC()))
        file1.close()
        plt.show()
        
    test_fb_lineout()

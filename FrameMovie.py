import outfile
import my_cmap
import os
import matplotlib.pyplot as plt
import numpy as np

class Frames:
    def __init__(self, code_name = 'osiris', simulation_path = None, frame_folder = 'Frames', plot_type = 0, use_num_list = False, start_num = 0, stride_num=1, count_num=1, if_scan_folder=False, driver_spec_name='driver', driver_vmin=None, driver_vmax=None, background_spec_name='e', background_vmin=None, background_vmax=None, trail_spec_name=None, trail_vmin=None, trail_vmax=None, if_e1=False, if_psi=False, dir=2, save_type='png', max_missing_file=0):
        self.code_name = code_name
        self.simulation_path = simulation_path
        self.frame_path = simulation_path+'/'+frame_folder
        self.use_num_list = use_num_list
        self.start_num = start_num
        self.stride_num = stride_num
        self.count_num = count_num
        # when if_scan_folder is True, scan the folder and generate a list of all available fils.
        self.if_scan_folder = if_scan_folder
        self.plot_type = plot_type
        self.driver_spec_name = driver_spec_name
        self.driver_vmin = driver_vmin
        self.driver_vmax = driver_vmax
        self.background_spec_name = background_spec_name
        self.background_vmin = background_vmin
        self.background_vmax = background_vmax
        self.trail_spec_name = trail_spec_name
        self.trail_vmin = trail_vmin
        self.trail_vmax = trail_vmax
        self.if_e1 = if_e1
        self.if_psi = if_psi
        self.dir = dir
        self.save_type = save_type
        #allowed number of missing files when doing plot loop. In HiPACE sometimes there are missing output files.
        self.max_missing_file = max_missing_file

################################property simulation_path################################
    def get_simulation_path(self):
        return self._simulation_path

    def set_simulation_path(self, value):
        if os.path.isdir(value):
            self._simulation_path = value
        else:
            raise IOError('\'{0}\' is not a folder or does not exist!'.format(value))

    simulation_path = property(get_simulation_path, set_simulation_path)

################################property frame_path################################
#frame_path is the folder for saving frames. Create the folder if it does not exist and its parent folder exist. Raise an exception if the path is a existing file or the folder is not empty.
    def get_frame_path(self):
        return self._frame_path

    def set_frame_path(self, value):
        tmp_path = os.path.abspath(value)
        head, tail = os.path.split(tmp_path)
        if os.path.isdir(head):
            if os.path.isfile(tmp_path):
                raise IOError('\'{0}\' is a file! Use a empty folder instead.'.format(tmp_path))
            elif os.path.isdir(tmp_path):
                if [] != os.listdir(tmp_path):
                    print('Warning: \'{0}\' is not empty! Existing files will not be overwritten.'.format(tmp_path))
            else:
                print('Warning: \'{0}\' does not exist! But we are creating it for use.'.format(tmp_path))
                os.makedirs(tmp_path)
            self._frame_path = tmp_path
        else:
            raise IOError('The parent dir \'{0}\' does not exist!'.format(head))

    frame_path = property(get_frame_path, set_frame_path)

################################property start_num################################
    def get_start_num(self):
        return self._start_num

    def set_start_num(self, value):
        if not isinstance(value, int):
            raise ValueError('start_num should be an interger!')
        if value<0:
            raise ValueError('start_num should not be negative!')
        self._start_num = value

    start_num = property(get_start_num, set_start_num)

################################property stride_num################################
    def get_stride_num(self):
        return self._stride_num

    def set_stride_num(self, value):
        if not isinstance(value, int):
            raise ValueError('stride_num should be an interger!')
        if value<1:
            raise ValueError('stride_num should not be at least 1!')
        self._stride_num = value

    stride_num = property(get_stride_num, set_stride_num)

################################property count_num################################
    def get_count_num(self):
        return self._count_num

    def set_count_num(self, value):
        if not isinstance(value, int):
            raise ValueError('count_num should be an interger!')
        if value<1:
            raise ValueError('count_num should not be at least 1!')
        self._count_num = value

    count_num = property(get_count_num, set_count_num)

################################property plot_type################################
    def get_plot_type(self):
        return self._plot_type

    def set_plot_type(self, value):
        self._plot_type = value

    plot_type = property(get_plot_type, set_plot_type)

################################method plot_beam_driven################################
#plot one frame for beam driven cases
    def plot_beam_driven(self, out_num):
        h_fig = plt.figure(figsize=(8,4.5))
        file1 = outfile.OutFile(code_name = self.code_name, path=self.simulation_path, field_name='charge', average='', spec_name=self.background_spec_name, use_num_list = self.use_num_list, out_num=out_num)
        h_ax = h_fig.add_subplot(111)
        h_ax.set_aspect('equal','box')
        file1.open()
        file1.read_data_slice(dir=self.dir)
        file1.plot_data(h_fig, h_ax, vmin=self.background_vmin, vmax=self.background_vmax, cmap='gray')
        file1._color_bar.set_label('$\\rho_e$')
        file1.close()

        file1.spec_name=self.driver_spec_name
        file1.open()
        file1.read_data_slice(dir=self.dir)
        file1.plot_data(h_fig, h_ax, vmin=self.driver_vmin, vmax=self.driver_vmax, cmap=my_cmap.cmap_higher_range_transparent(plt.cm.hot))
        file1._color_bar.set_label('$\\rho_d$')
        file1.close()

        if self.trail_spec_name is not None:
            file1.spec_name=self.trail_spec_name
            file1.open()
            file1.read_data_slice(dir=self.dir)
            file1.plot_data(h_fig, h_ax, vmin=self.trail_vmin, vmax=self.trail_vmax, cmap=my_cmap.cmap_higher_range_transparent())
            file1._color_bar.set_label('$\\rho_t$')
            file1.close()

        if self.if_e1:
            file1.field_name='e1'
            file1.open()
            file1.read_data_lineout()
            file1.plot_data(h_fig, h_ax, linestyle='-r', if_ylabel=False, multiple=1)
            file1.close()

        if self.if_psi:
            file1.field_name='psi'
            file1.open()
            file1.read_data_lineout()
            file1.plot_data(h_fig, h_ax, linestyle='-b', if_ylabel=False, multiple=1)
            file1.close()

        plt.tight_layout()
        return h_fig

################################method plot_laser_driven################################
#plot one frame for laser driven cases
    def plot_laser_driven(self, out_num):
        h_fig = plt.figure(figsize=(6.5,5))
        file1 = outfile.OutFile(path=self.simulation_path, field_name='charge', average='', spec_name=self.background_spec_name, use_num_list = self.use_num_list, out_num=out_num)
        h_ax = h_fig.add_subplot(111)
        #h_ax.set_aspect('equal', 'box')
        #plt.ylim(-10,10)
        file1.open()
        file1.read_data_slice(dir=self.dir)
        file1.plot_data(h_fig, h_ax, cmap='gray', vmin=-5.)
        file1._color_bar.set_label('$\\rho_e$')
        file1.close()

        if self.trail_spec_name is not None:
            file1.spec_name=self.trail_spec_name
            file1.open()
            file1.read_data_slice(dir=self.dir)
            file1.plot_data(h_fig, h_ax, vmin=self.trail_vmin, vmax=self.trail_vmax, cmap=my_cmap.cmap_higher_range_transparent())
            file1._color_bar.set_label('$\\rho_t$')
            file1.close()

        file1.field_name='e3'
        file1.open()
        file1.read_data_slice(dir=self.dir)
        file1.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_middle_range_transparent())
        file1._color_bar.set_label('$E_y$')
        file1.close()
        plt.tight_layout()
        return h_fig

################################method plot_p1x1################################
#plot one frame of p1x1 for species self.trail_spec_name
    def plot_p1x1(self, out_num):
        h_fig = plt.figure(figsize=(6.5,5))
        file1 = outfile.OutFile(path=self.simulation_path, field_name='p1x1', average='', spec_name=self.trail_spec_name, use_num_list = self.use_num_list, out_num=out_num)
        h_ax = h_fig.add_subplot(111)
        file1.open()
        file1.read_data()
        file1.plot_data(h_fig, h_ax, if_log_colorbar=True)
        #file1._color_bar.set_label('$\\rho_e$')
        file1.close()
        plt.tight_layout()
        return h_fig

################################method plot_x1x2################################
#plot one frame of phasespace x1x2 for species self.trail_spec_name. Phasespace is constructed from raw data
    def plot_x1x2_raw(self, out_num):
        h_fig = plt.figure(figsize=(6.5,5))
        file1 = outfile.OutFile(path=self.simulation_path, field_name='raw', spec_name=self.trail_spec_name, use_num_list = self.use_num_list, out_num=out_num)
        h_ax = h_fig.add_subplot(111)
        file1.open()
        file1.read_raw_q()
        file1.read_raw_x1()
        file1.read_raw_x2()
        file1.plot_raw_hist2D(h_fig, h_ax, dim=3, cmap=my_cmap.cmap_lower_range_transparent(), if_log_colorbar=False)
        file1.close()
        plt.tight_layout()
        return h_fig

################################method plot_x1x3################################
#plot one frame of phasespace x1x3 for species self.trail_spec_name. Phasespace is constructed from raw data
    def plot_x1x3_raw(self, out_num):
        h_fig = plt.figure(figsize=(6.5,5))
        file1 = outfile.OutFile(path=self.simulation_path, field_name='raw', spec_name=self.trail_spec_name, use_num_list = self.use_num_list, out_num=out_num)
        h_ax = h_fig.add_subplot(111)
        file1.open()
        file1.read_raw_q()
        file1.read_raw_x1()
        file1.read_raw_x3()
        file1.plot_raw_hist2D(h_fig, h_ax, dim=4, cmap=my_cmap.cmap_lower_range_transparent(), if_log_colorbar=False)
        file1.close()
        plt.tight_layout()
        return h_fig

################################method plot_x2p2################################
#plot one frame of phasespace x2p2 for species self.trail_spec_name. Phasespace is constructed from raw data
    def plot_x2p2_raw(self, out_num):
        h_fig = plt.figure(figsize=(6.5,5))
        file1 = outfile.OutFile(path=self.simulation_path, field_name='raw', spec_name=self.trail_spec_name, use_num_list = self.use_num_list, out_num=out_num)
        h_ax = h_fig.add_subplot(111)
        file1.open()
        file1.read_raw_q()
        file1.read_raw_x2()
        file1.read_raw_p2()
        file1.plot_raw_hist2D(h_fig, h_ax, dim=1, cmap=my_cmap.cmap_lower_range_transparent(), if_log_colorbar=False)
        file1.close()
        plt.tight_layout()
        return h_fig

################################method plot_x3p3################################
#plot one frame of phasespace x3p3 for species self.trail_spec_name. Phasespace is constructed from raw data
    def plot_x3p3_raw(self, out_num):
        h_fig = plt.figure(figsize=(6.5,5))
        file1 = outfile.OutFile(path=self.simulation_path, field_name='raw', spec_name=self.trail_spec_name, use_num_list = self.use_num_list, out_num=out_num)
        h_ax = h_fig.add_subplot(111)
        file1.open()
        file1.read_raw_q()
        file1.read_raw_x3()
        file1.read_raw_p3()
        file1.plot_raw_hist2D(h_fig, h_ax, dim=2, cmap=my_cmap.cmap_lower_range_transparent(), if_log_colorbar=False)
        file1.close()
        plt.tight_layout()
        return h_fig

################################method plot_save################################
#plot one frame and save using the method either plot_save_beam_driven or plot_save_laser_driven, depends on the "plot_type" property (0 for beam driver side view, 1 for laser driver side view, 2 for phase space "p1x1" of self.trail_spec_name)
    def plot_save(self, *args, **kwargs):
        if 'out_num' in kwargs:
            out_num = kwargs['out_num']
        else:
            out_num = args[0]
        save_file_name = '{0}/{1}.{2}'.format(self.frame_path, out_num, self.save_type)
        if os.path.isfile(save_file_name):
            print('Skipping existing number {}.'.format(out_num))
        else:
            print('Working on number {}.'.format(out_num))
            methods = (self.plot_beam_driven, self.plot_laser_driven, self.plot_p1x1, self.plot_x1x2_raw, self.plot_x1x3_raw, self.plot_x2p2_raw, self.plot_x3p3_raw)
            h_fig = methods[self.plot_type](*args, **kwargs)
            plt.savefig(save_file_name, format=self.save_type)
            plt.close(h_fig)

################################method save_frames################################
#save all frames
    def save_frames(self):
        print('Working on simulation \'{0}\' and saving frames at \'{1}\'.'.format(self.simulation_path, self.frame_path))
        missing_file = 0
        for i in range(self.start_num, self.start_num+self.count_num*self.stride_num, self.stride_num):
            try:
                self.plot_save(i)
                #set missing_file = 0 if success
                missing_file = 0
            except FileNotFoundError as err:
                missing_file = missing_file+1
                if missing_file>self.max_missing_file:
                    print('Iteration stops at frame number {0}. Exception message:\n{1}'.format(i, err))
                    break
                else: print('Warning! File No. {0} missing. Exception message:\n{1}'.format(i, err))
            except ValueError as err:
                print('Number {} does not exist. Seems all files are processed. Finishing...'.format(i))
                break

################################method make_movie################################
#make movie based on the saved frames
#    def make_movie(self):
#        print('Making movie.mpg under \'{0}\' - this may take a while.'.format(self.frame_path))
#        subprocess.call('mencoder \'{0}/*.png\' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o {0}/movie.mpg'.format(self.frame_path), shell=True)

if __name__ == '__main__':
    #frame1 = Frames(code_name = 'hipace', simulation_path = '/home/zming/mnt/X1_Shared_Pardis_Ming/NegBox/50um600pC1.1e16hh', frame_folder='Frames2', plot_type = 0, use_num_list = True, start_num = 0, stride_num=1, count_num=99999, background_spec_name='plasma_electrons', background_vmin=-5, driver_spec_name='driver', driver_vmin=-5, trail_spec_name='trailer', trail_vmax=0., trail_vmin=-11., if_e1=False, if_psi=False, max_missing_file=2, dir=2)
    frame1 = Frames(simulation_path = '/home/zming/mnt/X1_Shared_Pardis_Ming/NegBox/50um600pC1.1e16hh', frame_folder='os_Frames/Frames2', plot_type = 0, start_num = 0, stride_num=1, count_num=99999, background_spec_name='plasma', background_vmin=-5, driver_spec_name='driver', driver_vmin=-5, trail_spec_name='ramp', trail_vmax=0., trail_vmin=-11., if_e1=False, if_psi=False, dir=2)
    #frame1 = Frames(code_name = 'osiris', simulation_path = '/home/zming/simulations/os2D/os_DRI3D19', frame_folder='x3p3', plot_type = 6, start_num = 37, stride_num=1, count_num=99999, trail_spec_name='plasma')
    #frame1 = Frames(code_name = 'osiris', simulation_path = '/home/zming/simulations/os2D/os_DRI3D19', frame_folder='Frames1', plot_type = 0, start_num = 0, stride_num=1, count_num=99999, background_spec_name='plasma', background_vmin=-10, driver_spec_name='beam-driver', driver_vmin=-5, if_e1=True, if_psi=False, dir=1)
    #frame1 = Frames(code_name = 'osiris', simulation_path = '/home/zming/simulations/os2D/os_DRI3D19', frame_folder='p1x1', plot_type = 2, start_num = 16, stride_num=1, count_num=99999, trail_spec_name='plasma')
    #frame1 = Frames(code_name = 'osiris', simulation_path = '/home/zming/mnt/os_PT3D23', frame_folder='p1x1', plot_type = 2, start_num = 1, stride_num=1, count_num=99999, trail_spec_name='e')
    #frame1 = Frames(code_name = 'osiris', simulation_path = '/home/zming/mnt/os_PT3D23', frame_folder='Frames2', plot_type = 1, start_num = 0, stride_num=1, count_num=99999, background_spec_name='e', trail_spec_name=None, if_e1=False, if_psi=False, dir=2)
    #frame1 = Frames(code_name = 'osiris', simulation_path = '/home/zming/mnt/os_beam3D110', frame_folder='Frames2', plot_type = 0, start_num = 14, stride_num=1, count_num=99999, background_spec_name='e', background_vmin=-5, driver_spec_name='driver', driver_vmin=-5, trail_spec_name='He_e', trail_vmax=0, trail_vmin=-1, if_e1=True, if_psi=True, dir=2)
    #frame1 = Frames(code_name = 'hipace', simulation_path = '/home/zming/simulations/os2D/X1_Shared_Pardis_Ming/NegBox/50um600pC1.1e16', frame_folder='Frames2', plot_type = 0, start_num = 0, stride_num=20, count_num=99999, background_spec_name='plasma', background_vmin=-10, driver_spec_name='beam', driver_vmin=-50, if_e1=False, if_psi=False, max_missing_file=1, dir=2)
    frame1.save_frames()
    #frame1.make_movie()

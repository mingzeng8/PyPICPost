import outfile
import my_cmap
import os
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings(action='once')

# Do not have to start X-server. Save figures instead of show them.
import matplotlib as mpl
mpl.use('Agg')

class Frames:
    def __init__(self, code_name = 'osiris', simulation_path = None, frame_folder = None, plot_type = 'laser_driven', plot_spec_name=None, average='', use_num_list = False, start_num = 0, stride_num=1, count_num=1, laser_field_name='e3', if_laser_profile=False, laser_cmap=None, density_type='charge', driver_spec_name='driver', driver_cmap=None, driver_vmin=None, driver_vmax=None, background_spec_name='e', background_cmap=None, background_vmin=None, background_vmax=None, trail_spec_name=None, trail_cmap=None, trail_vmin=None, trail_vmax=None, if_e1=False, e1_multiple=1., if_psi=False, psi_multiple=1., if_driver_cm=False, if_trail_cm=False, dir=2, project_dir=None, save_type='png', max_missing_file=0, plot_range=None):
        self.code_name = code_name
        self.simulation_path = simulation_path
        if frame_folder is None:
            if code_name == 'osiris': frame_pref = 'os_Frames'
            elif code_name == 'hipace': frame_pref = 'hi_Frames'
            elif code_name == 'quickpic': frame_pref = 'qp_Frames'
            elif code_name == 'fbpic': frame_pref = 'fb_Frames'
            else: raise NotImplementedError('code_name {} not implemented.'.format(code_name))
            if plot_type not in {'laser_driven', 'beam_driven'}: frame_surf = '{}_{}'.format(plot_spec_name, plot_type)
            elif dir==0: frame_surf = 'slice_x-y'
            elif dir==1:
                if code_name == 'quickpic': frame_surf = 'slice_x-z'
                else: frame_surf = 'slice_y-z'
            elif dir==2:
                if code_name == 'quickpic': frame_surf = 'slice_y-z'
                else: frame_surf = 'slice_x-z'
            else: raise ValueError('dir should be 0, 1, or 2.')
            frame_folder = frame_pref+'/'+frame_surf
        self.frame_path = simulation_path+'/'+frame_folder
        # when use_num_list is True, scan the folder and generate a list of all available fils.
        self.use_num_list = use_num_list
        self.start_num = start_num
        self.stride_num = stride_num
        self.count_num = count_num
        self.plot_type = plot_type
        # For phasespace plots, i.e. plot_type = 'p1x1' or so on, one has to give plot_spec_name
        self.plot_spec_name = plot_spec_name
        self.laser_field_name = laser_field_name
        self.if_laser_profile = if_laser_profile
        self.laser_cmap = laser_cmap
        self.density_type = density_type
        self.driver_spec_name = driver_spec_name
        self.driver_cmap = driver_cmap
        self.driver_vmin = driver_vmin
        self.driver_vmax = driver_vmax
        self.background_spec_name = background_spec_name
        self.background_cmap = background_cmap
        self.background_vmin = background_vmin
        self.background_vmax = background_vmax
        self.trail_spec_name = trail_spec_name
        self.trail_cmap = trail_cmap
        self.trail_vmin = trail_vmin
        self.trail_vmax = trail_vmax
        self.if_e1 = if_e1
        self.e1_multiple = e1_multiple
        self.if_psi = if_psi
        self.psi_multiple = psi_multiple
        # if_driver_cm = True means plot the lineout of driver beam center of mass in the beam driven case
        self.if_driver_cm = if_driver_cm
        # if_trail_cm = True means plot the lineout of trail beam center of mass
        self.if_trail_cm = if_trail_cm
        self.dir = dir
        # If projection direction for the phasespace plots
        self.project_dir = project_dir
        self.save_type = save_type
        #allowed number of missing files when doing plot loop. In HiPACE sometimes there are missing output files.
        self.max_missing_file = max_missing_file
        self.plot_range = plot_range
        # Initialize outfile object, especially initialize avail_num_list if use_num_list.
        # We will look up for full grid dump first.
        if self.plot_type in {'beam_driven', 'laser_driven'}:
            try:
                # If use_num_list is True, outfile.OutFile() will try to setup a file list while initialize. But if slice dump is not found, we try to find full grid dump file and do slicing.
                # For QuickPIC and OSIRIS, fields may saved as slices
                # fld_slice in {1, 2, 3} while dir in (0, 1, 2)
                self.outfile = outfile.OutFile(code_name = code_name, path=simulation_path, field_name=self.density_type, average=average, spec_name=background_spec_name, use_num_list = use_num_list, out_num=start_num, fld_slice=dir+1)
                if not os.path.isfile(self.outfile.path_filename):
                    raise FileNotFoundError
            except FileNotFoundError:
                # For full grid dumps
                self.outfile = outfile.OutFile(code_name = code_name, path=simulation_path, field_name=self.density_type, average=average, spec_name=background_spec_name, use_num_list = use_num_list, out_num=start_num)
        else:
            self.outfile = outfile.OutFile(code_name = code_name, path=simulation_path, field_name='raw', spec_name=plot_spec_name, use_num_list = use_num_list, out_num=start_num)

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
        head, _ = os.path.split(tmp_path)
        head_of_head, _ = os.path.split(head)
        if os.path.isdir(head_of_head):
            if os.path.isfile(tmp_path):
                raise IOError('\'{0}\' is a file! Use a empty folder instead.'.format(tmp_path))
            elif os.path.isdir(tmp_path):
                if [] != os.listdir(tmp_path):
                    print('Warning: \'{0}\' is not empty!'.format(tmp_path))
            else:
                print('Warning: \'{0}\' does not exist! But we are creating it for use.'.format(tmp_path))
                os.makedirs(tmp_path)
            self._frame_path = tmp_path
        else:
            raise IOError('The parent dir \'{0}\' does not exist!'.format(head_of_head))

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

################################property laser_field_name################################
    def get_laser_field_name(self):
        return self._laser_field_name

    def set_laser_field_name(self, value):
        laser_fields = {'e1','e2','e3',None}
        if value not in laser_fields: raise ValueError('laser_field_name = {}, but allowed values are {}.'.format(value, laser_fields))
        self._laser_field_name = value

    laser_field_name = property(get_laser_field_name, set_laser_field_name)

################################property density_type################################
    def get_density_type(self):
        return self._density_type

    def set_density_type(self, value):
        allowed_density_types = {'charge','m'}
        if value not in allowed_density_types: raise ValueError('density_type = {}, but allowed values are {}.'.format(value, allowed_density_types))
        self._density_type = value

    density_type = property(get_density_type, set_density_type)

################################property laser_cmap################################
    def get_laser_cmap(self):
        return self._laser_cmap

    def set_laser_cmap(self, value):
        if value is None:
            # set to default
            if self.if_laser_profile: self._laser_cmap = plt.cm.Reds
            else: self._laser_cmap = plt.cm.bwr
        else: self._laser_cmap = value

    laser_cmap = property(get_laser_cmap, set_laser_cmap)

################################property driver_cmap################################
    def get_driver_cmap(self):
        return self._driver_cmap

    def set_driver_cmap(self, value):
        if value is None:
            # set to default
            if self.density_type=='charge': self._driver_cmap = plt.cm.hot
            else: self._driver_cmap = plt.cm.hot.reversed()
        else: self._driver_cmap = value

    driver_cmap = property(get_driver_cmap, set_driver_cmap)

################################property background_cmap################################
    def get_background_cmap(self):
        return self._background_cmap

    def set_background_cmap(self, value):
        if value is None:
            # set to default
            if self.density_type=='charge': self._background_cmap = 'gray'
            else: self._background_cmap = 'gray_r'
        else: self._background_cmap = value

    background_cmap = property(get_background_cmap, set_background_cmap)

################################property trail_cmap################################
    def get_trail_cmap(self):
        return self._trail_cmap

    def set_trail_cmap(self, value):
        if value is None:
            # set to default
            if self.density_type=='charge': self._trail_cmap = plt.cm.jet
            else: self._trail_cmap = plt.cm.jet.reversed()
        else: self._trail_cmap = value

    trail_cmap = property(get_trail_cmap, set_trail_cmap)

################################method plot_background################################
#plot background plasma
    def plot_background(self, h_fig, h_ax):
        self.outfile.field_name=self.density_type
        self.outfile.spec_name=self.background_spec_name
        self.outfile.open()
        if self.outfile.num_dimensions>=3: self.outfile.read_data_slice(dir=self.dir)
        else: self.outfile.read_data()
        self.outfile.plot_data(h_fig, h_ax, vmin=self.background_vmin, vmax=self.background_vmax, cmap=self.background_cmap)
        self.outfile._color_bar.set_label('$\\rho_e$' if self.density_type=='charge' else '$n_e$')
        self.outfile.close()

################################method plot_trail################################
#plot trailing beam
    def plot_trail(self, h_fig, h_ax):
        self.outfile.field_name=self.density_type
        self.outfile.spec_name=self.trail_spec_name
        try:
            self.outfile.open()
            if self.outfile.num_dimensions>=3: self.outfile.read_data_slice(dir=self.dir)
            else: self.outfile.read_data()
            self.outfile.plot_data(h_fig, h_ax, vmin=self.trail_vmin, vmax=self.trail_vmax, cmap=(my_cmap.cmap_higher_range_transparent(self.trail_cmap) if self.density_type=='charge' else my_cmap.cmap_lower_range_transparent(self.trail_cmap)))
            self.outfile._color_bar.set_label('$\\rho_t$' if self.density_type=='charge' else '$n_t$')
            if self.if_trail_cm:
                try:
                    self.outfile.data_center_of_mass2d()
                    self.outfile.plot_data(h_fig, h_ax, c='m', ls='--')
                # RuntimeError occurs when there is no particle
                except RuntimeError: pass
            self.outfile.close()
        except FileNotFoundError: warnings.warn("Trailer file '{}' not found. Trailer is not ploted.".format(self.outfile.path_filename))
        except KeyError: warnings.warn("'{}' does not contain '{}'. Trailer is not ploted.".format(self.outfile.path_filename, self.outfile.spec_name))

################################method plot_e1_line################################
#plot e1 lineout
    def plot_e1_line(self, h_fig, h_ax):
        self.outfile.field_name='e1'
        self.outfile.open()
        self.outfile.read_data_lineout()
        self.outfile.plot_data(h_fig, h_ax, if_ylabel=False, multiple=self.e1_multiple, c='r', ls='-')
        self.outfile.close()
        
################################method plot_psi_line################################
#plot psi lineout
    def plot_psi_line(self, h_fig, h_ax):
        self.outfile.field_name='psi'
        self.outfile.open()
        self.outfile.read_data_lineout()
        self.outfile.plot_data(h_fig, h_ax, if_ylabel=False, multiple=self.psi_multiple, c='b', ls='-')
        self.outfile.close()

################################method plot_e1_psi_line################################
#plot e1 and psi lineout
    def plot_e1_psi_line(self, h_fig, h_ax):
        if self.if_e1: self.plot_e1_line(h_fig, h_ax)
        if self.if_psi:
            try: self.plot_psi_line(h_fig, h_ax)
            except (FileNotFoundError, KeyError):
            # Obtain psi by integrating e1
                if not self.if_e1:
                    # Have not read e1, read e1 first
                    self.outfile.field_name='e1'
                    self.outfile.open()
                    self.outfile.read_data_lineout()
                    self.outfile.close()
                e1 = self.outfile._data
                if self.code_name != 'quickpic':
                    # For code except quickpic, psi is e1 integration from right to left
                    e1 = np.flip(e1)
                np.cumsum(e1*self.outfile._axis_slices[0].step, out=self.outfile._data)
                if self.code_name != 'quickpic':
                    # For code except quickpic, flip the data
                    self.outfile._data = np.flip(self.outfile._data)
                self.outfile.plot_data(h_fig, h_ax, if_ylabel=False, multiple=self.psi_multiple, c='b', ls='-')

################################method plot_beam_driven################################
#plot one frame for beam driven cases
    def plot_beam_driven(self, out_num, **kwargs):
        h_fig = plt.figure(**kwargs)
        self.outfile.out_num=out_num
        h_ax = h_fig.add_subplot(111)
        self.plot_background(h_fig, h_ax)

        if self.driver_spec_name is not None:
            self.outfile.field_name=self.density_type
            self.outfile.spec_name=self.driver_spec_name
            try:
                self.outfile.open()
                if self.outfile.num_dimensions>=3: self.outfile.read_data_slice(dir=self.dir)
                else: self.outfile.read_data()
                self.outfile.plot_data(h_fig, h_ax, vmin=self.driver_vmin, vmax=self.driver_vmax, cmap=(my_cmap.cmap_higher_range_transparent(self.driver_cmap) if self.density_type=='charge' else my_cmap.cmap_lower_range_transparent(self.driver_cmap)))
                self.outfile._color_bar.set_label('$\\rho_d$' if self.density_type=='charge' else '$n_d$')
                if self.if_driver_cm:
                    try:
                        self.outfile.data_center_of_mass2d()
                        self.outfile.plot_data(h_fig, h_ax, c='c', ls='--')
                    # RuntimeError occurs when there is no particle
                    except RuntimeError: pass
                self.outfile.close()
            except FileNotFoundError: warnings.warn('Driver file {} not found. Driver is not ploted.'.format(self.outfile.path_filename))
            except KeyError: warnings.warn("'{}' does not contain '{}'. Driver is not ploted.".format(self.outfile.path_filename, self.outfile.spec_name))

        if self.trail_spec_name is not None:
            self.plot_trail(h_fig, h_ax)

        self.plot_e1_psi_line(h_fig, h_ax)
        h_ax.set_aspect('equal','box')
        plt.tight_layout()
        return h_fig

################################method plot_laser_driven################################
#plot one frame for laser driven cases
    def plot_laser_driven(self, out_num, **kwargs):
        if 'figsize' in kwargs: figsize = kwargs['figsize']
        else: figsize = None
        h_fig = plt.figure(figsize=figsize)
        self.outfile.out_num=out_num
        h_ax = h_fig.add_subplot(111)
        self.plot_background(h_fig, h_ax)

        if self.laser_field_name is not None:
            self.outfile.field_name=self.laser_field_name
            try:
                self.outfile.open()
                if self.outfile.num_dimensions>=3: self.outfile.read_data_slice(dir=self.dir)
                else: self.outfile.read_data()
                # When if_laser_profile is Ture, plot the laser profile instead of original E-field
                if self.if_laser_profile:
                    self.outfile.data_profile2d()
                    self.outfile.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_lower_range_transparent(self.laser_cmap, transparency_transition_region=[0.15,0.4]))
                else: self.outfile.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_middle_range_transparent(self.laser_cmap))
                self.outfile._color_bar.set_label('$E_L$')
                self.outfile.close()
            except FileNotFoundError: warnings.warn('Laser file {} not found. Laser is not ploted.'.format(self.outfile.path_filename))
            except KeyError: warnings.warn("'{}' does not contain '{}'. Trailer is not ploted.".format(self.outfile.path_filename, self.field_name))

        if self.trail_spec_name is not None:
            self.plot_trail(h_fig, h_ax)

        self.plot_e1_psi_line(h_fig, h_ax)
        h_ax.set_aspect('equal', 'box')
        if 'ylim' in kwargs:
            plt.ylim(kwargs['ylim'])
        plt.tight_layout()
        return h_fig

################################method plot_p1x1################################
#plot one frame of p1x1 for species self.trail_spec_name
# This method is deprecated. It's function can be replaced by plot_phasespace_raw
#    def plot_p1x1(self, out_num):
#        h_fig = plt.figure(figsize=(6.5,5))
#        file1 = outfile.OutFile(path=self.simulation_path, field_name='p1x1', average='', spec_name=self.trail_spec_name, use_num_list = self.use_num_list, out_num=out_num)
#        h_ax = h_fig.add_subplot(111)
#        file1.open()
#        file1.read_data()
#        file1.plot_data(h_fig, h_ax, if_log_colorbar=True)
#        #file1._color_bar.set_label('$\\rho_e$')
#        file1.close()
#        plt.tight_layout()
#        return h_fig

################################method plot_phasespace_raw################################
#plot one frame of phasespace for species self.trail_spec_name. Phasespace is constructed from raw data. Phasespace type is determined by self.plot_type.
    def plot_phasespace_raw(self, out_num):
        h_fig = plt.figure(figsize=(6.5,5))
        self.outfile.spec_name=self.plot_spec_name
        self.outfile.out_num=out_num
        h_ax = h_fig.add_subplot(111)
        self.outfile.open()
        try:
            self.outfile.plot_raw_hist2D(h_fig, h_ax, dims=self.plot_type, cmap=my_cmap.cmap_lower_range_transparent(), if_log_colorbar=False, select_range=self.plot_range)
            self.outfile.close()
            if self.project_dir is not None:
                # project_dir == 0 or 1. For 0, flip x and y
                self.outfile.data_project2d(dir = self.project_dir)
                if_flip_xy=(not bool(self.project_dir))
                get_lim = (h_ax.get_xlim, h_ax.get_ylim)
                lim = get_lim[self.project_dir]()
                self.outfile.plot_data(h_fig, h_ax, if_flip_xy=(not bool(self.project_dir)), multiple=(lim[1]-lim[0])/2/self.outfile.data.max(), offset=lim[0], color='k')
        except RuntimeError:
            # RuntimeError is raised if too few particle is found in the raw file.
            print('Cannot plot 2D histogram. Try next...')
            self.outfile.close()
            plt.close(h_fig)
            return None
        self.outfile.close()
        plt.tight_layout()
        return h_fig

################################method plot_save################################
#plot one frame and save using the method either plot_save_beam_driven or plot_save_laser_driven, depends on the "plot_type" property (0 for beam driver side view, 1 for laser driver side view, 2 for phase space "p1x1" of self.trail_spec_name)
    def plot_save(self, out_num, overwrite = False, **kwargs):
        '''if 'out_num' in kwargs:
            out_num = kwargs['out_num']
        else:
            out_num = args[0]
        if 'overwrite' in kwargs:
            overwrite = kwargs['overwrite']
            kwargs.pop('overwrite')
        else:
            overwrite = False'''
        save_file_name = '{0}/{1}.{2}'.format(self.frame_path, out_num, self.save_type)
        if os.path.isfile(save_file_name):
            if overwrite:
                print('Overwritting number {}.'.format(out_num))
            else:
                print('Skipping existing number {}.'.format(out_num))
                return
        else:
            print('Working on number {}.'.format(out_num))
        if 'beam_driven' == self.plot_type: h_fig = self.plot_beam_driven(out_num=out_num, **kwargs)
        elif 'laser_driven' == self.plot_type: h_fig = self.plot_laser_driven(out_num=out_num, **kwargs)
        else: h_fig = self.plot_phasespace_raw(out_num=out_num, **kwargs)
        if h_fig is not None:
            plt.savefig(save_file_name, format=self.save_type)
            plt.close(h_fig)

################################method save_frames################################
#save all frames
    def save_frames(self, **kwargs):
        print('Working on simulation \'{0}\' and saving frames at \'{1}\'.'.format(self.simulation_path, self.frame_path))
        missing_file = 0
        for i in range(self.start_num, self.start_num+self.count_num*self.stride_num, self.stride_num):
            try:
                self.plot_save(i, **kwargs)
                #set missing_file = 0 if success
                missing_file = 0
            except FileNotFoundError as err:
                missing_file = missing_file+1
                if missing_file>self.max_missing_file:
                    print('Iteration stops at frame number {0}. Exception message:\n{1}'.format(i, err))
                    break
                else: print('Warning! File No. {0} missing. Exception message:\n{1}'.format(i, err))
            except Exception as err:
                print('Number {} does not exist. Seems all files are processed. Finishing...'.format(i))
                break

################################method make_movie################################
#make movie based on the saved frames
#    def make_movie(self):
#        print('Making movie.mpg under \'{0}\' - this may take a while.'.format(self.frame_path))
#        subprocess.call('mencoder \'{0}/*.png\' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o {0}/movie.mpg'.format(self.frame_path), shell=True)

if __name__ == '__main__':
    # An example
    frame1 = Frames(simulation_path = '/home/ming/mnt/BSCC_HOME/jobs/os_DCLBII/2', plot_type = 'laser_driven', use_num_list = True, start_num = 0, stride_num=1, count_num=99999, laser_field_name = 'e3', background_spec_name='e', background_vmin=-5, trail_spec_name='O_e', trail_vmin=-20, if_e1=True, if_psi=True, dir=2)
    #frame1.save_frames()

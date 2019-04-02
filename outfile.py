import numpy as np
import h5py
import matplotlib.pyplot as plt
import TwoDGaussianFit as tdgf
from scipy import signal, pi
from scipy import e as const_e
from scipy.constants import e as e_charge
from matplotlib.colors import LogNorm

float_type=np.float64

class OutFile:
    def __init__(self, code_name = 'osiris', path = '.', out_type = None, field_name = 'e3', spec_name = '', out_num = 0, average='', cyl_m_num = 0, cyl_m_re_im='re'):
##value digit_num##
        self.digit_num = 6
        self.code_name = code_name
        self.path = path
        self.spec_name = spec_name
        self._field_name_to_out_type = {'psi':'FLD', 'e1':'FLD', 'e2':'FLD', 'e3':'FLD', 'e3_cyl_m':'FLD_CYL_M', 'b1':'FLD', 'b2':'FLD', 'b3':'FLD', 'j1':'DENSITY', 'ene':'DENSITY', 'charge':'DENSITY', 'ion_charge':'ION', 'p1x1':'PHA', 'p2x2':'PHA', 'raw':'RAW'}
        #self.out_type = out_type
        self.field_name = field_name
        self.out_num = out_num
        self.average = average
        #2D cylindrical modes for fields
        self.cyl_m_num = cyl_m_num
        self.cyl_m_re_im = cyl_m_re_im

        self._num_dimensions = 0
        self._axis_labels_original = ('z', 'x', 'y', '$p_z$', '$p_x$', '$p_y$')
        self._axis_units_original = ('$c / \\omega_p$',)*3 + ('$m_ec$',)*3
        #self._axis_units_original = ('$k_p^{-1}$',)*3 + ('$m_ec$',)*3
        self._field_names = {'psi':'$\psi$', 'e1':'$E_z$', 'e2':'$E_x$', 'e3':'$E_y$', 'e3_cyl_m':'$E_y$', 'b1':'$B_z$', 'b2':'$B_x$', 'b3':'$B_y$', 'j1':'$J_z$', 'ene':'$E_k / n_p m_e c^2$', 'charge':'$\\rho$', 'ion_charge':'$\\rho$', 'p1x1':'$p_1x_1$ [arb. units]', 'p2x2':'$p_2x_2$ [arb. units]', 'beam_charge':'$\\rho_b$', 'plasma_charge':'$\\rho_e$'}
        self._accepted_out_type = {'DENSITY', 'FLD', 'FLD_CYL_M', 'ION', 'PHA', 'RAW'}
            

################################property code_name################################
    def get_code_name(self):
        return self._code_name

    def set_code_name(self, value):
        if value not in {'osiris',}:
            raise ValueError('code_name \'{0}\' not implemented!'.format(value))
        self._code_name = value

    code_name = property(get_code_name, set_code_name)

################################property path################################
    def get_path(self):
        return self._path

    def set_path(self, value):
        if not isinstance(value, str):
            raise ValueError('path should be a string!')
        self._path = value

    path = property(get_path, set_path)

################################property out_type################################
    def get_out_type(self):
        return self._out_type

    def set_out_type(self, value):
        #if value not in self._accepted_out_type:
        #    raise ValueError('out_type \'{0}\' not implemented!'.format(value))
        #self._out_type = value
        raise RuntimeError('Please do not set out_type directly!. It will be automatically set when setting field_name.')

    out_type = property(get_out_type, set_out_type)

################################property field_name################################
    def get_field_name(self):
        return self._field_name

    def set_field_name(self, value):
        if value not in self._field_name_to_out_type:
            raise ValueError('field_name \'{0}\' not yet implemented!'.format(value))
        self._field_name = value
        self._out_type = self._field_name_to_out_type[self._field_name]

    field_name = property(get_field_name, set_field_name)

################################property spec_name################################
    def get_spec_name(self):
        return self._spec_name

    def set_spec_name(self, value):
        if not isinstance(value, str):
            raise TypeError('spec_name should be a string!'.format(value))
        self._spec_name = value

    spec_name = property(get_spec_name, set_spec_name)

################################property out_num################################
    def get_out_num(self):
        return self._out_num

    def set_out_num(self, value):
        if value not in range(int(10**self.digit_num)):
            raise ValueError('out_num \'{0}\' not in range!'.format(value))
        self._out_num = value
##value num_of_zeros##
        if 0 == value:
            self.num_of_zeros = self.digit_num - 1
        else:
            self.num_of_zeros = self.digit_num - 1 - int(np.log10(value))

    out_num = property(get_out_num, set_out_num)

################################property average################################
    def get_average(self):
        return self._average

    def set_average(self, value):
        if value not in ('','-savg'):
            raise NotImplementedError('average = \'{0}\' not implemented!'.format(value))
        self._average = value

    average = property(get_average, set_average)

################################property path_filename################################
    def get_path_filename(self):
        if 'DENSITY' == self._out_type:
            return '{0}/MS/DENSITY/{2}/{1}{5}/{1}{5}-{2}-{3}{4}.h5'.format(self.path, self.field_name, self.spec_name, '0'*self.num_of_zeros, self.out_num, self.average)
        elif 'ION' == self._out_type:
            return '{0}/MS/ION/{2}/{1}/{1}-{2}-{3}{4}.h5'.format(self.path, self.field_name, self.spec_name, '0'*self.num_of_zeros, self.out_num)
        elif 'FLD' == self._out_type:
            return '{0}/MS/FLD/{1}{4}/{1}{4}-{2}{3}.h5'.format(self.path, self.field_name, '0'*self.num_of_zeros, self.out_num, self.average)
        elif 'FLD_CYL_M' == self._out_type:
            return '{0}/MS/FLD/MODE-{1}-{2}/{3}/{3}-{1}-{4}-{5}{6}.h5'.format(self.path, self.cyl_m_num, self.cyl_m_re_im.upper(), self.field_name, self.cyl_m_re_im.lower(), '0'*self.num_of_zeros, self.out_num)
        elif 'PHA' == self._out_type:
            return '{0}/MS/PHA/{1}/{2}/{1}-{2}-{3}{4}.h5'.format(self.path, self.field_name, self.spec_name, '0'*self.num_of_zeros, self.out_num)
        elif 'RAW' == self._out_type:
            return '{0}/MS/RAW/{1}/RAW-{1}-{2}{3}.h5'.format(self.path, self.spec_name, '0'*self.num_of_zeros, self.out_num)

    def set_path_filename(self):
        raise RuntimeError('You cannot set path_filename directly!')

    path_filename = property(get_path_filename, set_path_filename)

################################property num_dimensions################################
    def get_num_dimensions(self):
        return self._num_dimensions

    def set_num_dimensions(self, value):
        raise RuntimeError('You cannot set num_dimensions directly! The OutFile.open() method sets it automatically.')

    num_dimensions = property(get_num_dimensions, set_num_dimensions)

################################property axis_ranges################################
    def get_axis_range(self):
        return self._axis_range

    def set_axis_range(self, value):
        raise RuntimeError('You cannot set axis_range directly! The OutFile.open() method sets it automatically.')

    axis_range = property(get_axis_range, set_axis_range)

################################property data################################
    def get_data(self):
        return self._data

    def set_data(self, value):
        raise RuntimeError('You cannot set data directly! The OutFile.read_data() method sets it.')

    data = property(get_data, set_data)

################################property time################################
    def get_time(self):
        return self.fileid.attrs['TIME'][0]

    def set_time(self, value):
        raise RuntimeError('You cannot set time directly!')

    time = property(get_time, set_time)

################################method open################################
    def open(self):
        self.fileid = h5py.File(self.path_filename,'r')
        if self._out_type in self._accepted_out_type:
            if 'RAW' == self._out_type:
                return
##value _num_dimensions##
            try:
                #try to read the dimension from 'AXIS'
                self._num_dimensions = len(self.fileid['AXIS'].keys())
##value _axis_range##
                self._axis_range = np.zeros((2,self._num_dimensions), dtype=float_type)
                for i in range(self._num_dimensions):
                    self.fileid['AXIS/AXIS{0}'.format(i+1)].read_direct(self._axis_range, np.s_[:], np.s_[:,i])
            except KeyError:
                #when 'AXIS' does not exist, eg. output from HiPACE, read from the attribute
                xmax = self.fileid.attrs.get('XMAX')
                xmin = self.fileid.attrs.get('XMIN')
                self._num_dimensions = len(xmax)
##value _axis_range##
                self._axis_range = np.zeros((2,self._num_dimensions), dtype=float_type)
                for i in range(self._num_dimensions):
                    self._axis_range[0,i] = xmin[i]
                    self._axis_range[1,i] = xmax[i]
            for i in self.fileid.keys():
                if i!='AXIS':
                    break
##value _data_name_in_file##
            self._data_name_in_file = i
            #In osiris, the data 3 dimensions correspond to dir = 2, 1, 0, respectively. So we need to flip the shape tuple here.
##value _cell_size##
            self._cell_size = (self._axis_range[1, :] - self._axis_range[0, :]) / np.flipud(self.fileid[self._data_name_in_file].shape)

################################method close################################
    def close(self):
        self.fileid.close()

################################method read_raw_q################################
    def read_raw_q(self):
        self._raw_q = np.zeros(self.fileid['q'].shape, dtype=float_type)
        self.fileid['q'].read_direct(self._raw_q)

################################method read_raw_ene################################
    def read_raw_ene(self):
        try:
            self._raw_ene = np.zeros(self.fileid['ene'].shape, dtype=float_type)
            self.fileid['ene'].read_direct(self._raw_ene)
        except KeyError:
            print('Warning! Key \'ene\' does not exist in particle raw data! Reading p1, p2, p3 and doing ene=sqrt(p1^2+p2^2+p3^2)-1 instead. Please make sure p1, p2, p3 are read before this!')
            self._raw_ene = np.sqrt(np.square(self._raw_p1)+np.square(self._raw_p2)+np.square(self._raw_p3))-1.

################################method read_raw_x1################################
    def read_raw_x1(self):
        self._raw_x1 = np.zeros(self.fileid['x1'].shape, dtype=float_type)
        self.fileid['x1'].read_direct(self._raw_x1)

################################method read_raw_x2################################
    def read_raw_x2(self):
        self._raw_x2 = np.zeros(self.fileid['x2'].shape, dtype=float_type)
        self.fileid['x2'].read_direct(self._raw_x2)

################################method read_raw_x3################################
    def read_raw_x3(self):
        self._raw_x3 = np.zeros(self.fileid['x3'].shape, dtype=float_type)
        self.fileid['x3'].read_direct(self._raw_x3)

################################method read_raw_p1################################
    def read_raw_p1(self):
        self._raw_p1 = np.zeros(self.fileid['p1'].shape, dtype=float_type)
        self.fileid['p1'].read_direct(self._raw_p1)

################################method read_raw_p2################################
    def read_raw_p2(self):
        self._raw_p2 = np.zeros(self.fileid['p2'].shape, dtype=float_type)
        self.fileid['p2'].read_direct(self._raw_p2)

################################method read_raw_p3################################
    def read_raw_p3(self):
        self._raw_p3 = np.zeros(self.fileid['p3'].shape, dtype=float_type)
        self.fileid['p3'].read_direct(self._raw_p3)

################################method calculate_q_pC################################
# calculate the charge in unit of pico-coulumb
# please call read_raw_q() before this function
    def calculate_q_pC(self, n0_per_cc, dx_norm, dy_norm, dz_norm):
        '''
        n0_per_cc: reference density in simulation, in unit of per centimeter cube
        dx_norm, dy_norm, dz_norm: normalized dx, dy, dz
        one_over_k0_um: one over k0, in unit of micrometer
        one_over_k0_um = (3.541e-20*n0_per_cc)^-0.5
        '''
        return np.sum(self._raw_q)*e_charge*dx_norm*dy_norm*dz_norm/np.sqrt(n0_per_cc)*1.5e29

################################method calculate_norm_rms_emittance_um################################
# calculate the normalized rms emittance in unit of micrometer radian
# please call read_raw_x2(), read_raw_p2() and/or read_raw_x3(), read_raw_p3()
# and/or read_raw_x1(), read_raw_p1()
# and read_raw_q() before this function
    def calculate_norm_rms_emittance_um(self, n0_per_cc, directions=(2,)):
        '''
        n0_per_cc: reference density in simulation, in unit of per centimeter cube
        one_over_k0_um: one over k0, in unit of micrometer
        one_over_k0_um = (3.541e-20*n0_per_cc)^-0.5
        '''
        one_over_k0_um = (3.541e-20*n0_per_cc)**-0.5
        emittances=np.zeros(len(directions))
        weights=np.absolute(self._raw_q)
        sum_weight=np.sum(weights)
        for i in range(len(directions)):
            if 3==directions[i]:
                x=self._raw_x3
                p=self._raw_p3
            elif 1==directions[i]:
                x=self._raw_x1
                p=self._raw_p1
            else:
                x=self._raw_x2
                p=self._raw_p2
            term1=np.sum(np.multiply(np.square(x), weights))/sum_weight
            term2=np.sum(np.multiply(np.square(p), weights))/sum_weight
            term3=np.sum(np.multiply(np.multiply(x, p), weights))/sum_weight
            emittances[i]=np.sqrt(term1*term2-term3*term3)*one_over_k0_um
        return emittances

################################method read_data################################
    def read_data(self):
        self._data = np.zeros(self.fileid[self._data_name_in_file].shape, dtype=float_type)
        self.fileid[self._data_name_in_file].read_direct(self._data)
        self._axis_slices = [slice(self._axis_range[0, i], self._axis_range[1, i], self._cell_size[i]) for i in range(self._num_dimensions)]
        if 'p1x1' == self.field_name:
            self._axis_labels = [self._axis_labels_original[0], self._axis_labels_original[3]]
            self._axis_units = [self._axis_units_original[0], self._axis_units_original[3]]
        else:
            self._axis_labels = [self._axis_labels_original[i] for i in range(self._num_dimensions)]
            self._axis_units = [self._axis_units_original[i] for i in range(self._num_dimensions)]
        self._fig_title = 't = {0:.2f}'.format(self.fileid.attrs['TIME'][0])

################################method read_data_slice################################
    def read_data_slice(self, dir = 2, pos = None):
        '''dir is the direction perpendicular to the slice plane'''
        if 3 > self._num_dimensions:
            raise RuntimeError('Method OutFile.read_data_slice() cannot work on data with dimension number < 3!')
        if dir not in range (3):
            raise ValueError('Slice direction should be 0, 1 or 2!')
        tmp_len = self.fileid[self._data_name_in_file].shape[2-dir]
        if pos is None:
            #set pos at the middle of the box
            pos = (self._axis_range[0, dir] + self._axis_range[1, dir])/2.
            pos_index = int(tmp_len/2.)
        else:
            #get the index of the nearest grid. +0.5 here has a similar effect of rounding.
            pos_index = int((pos - self._axis_range[0, dir]) / self._cell_size[dir] + 0.5)
            if tmp_len<=pos_index:
                print('Warning: pos is larger than the upper bundary! Force slicing at the upper bundary.')
                pos_index = tmp_len-1
                pos = self._axis_range[1, dir]
            elif 0>pos_index:
                print('Warning: pos is smaller than the lower bundary! Force slicing at the lower bundary.')
                pos_index = 0
                pos = self._axis_range[0, dir]
        new_shape = [self.fileid[self._data_name_in_file].shape[i] for i in range(3) if 2-dir != i]
        self._data = np.zeros(new_shape, dtype=float_type)
        slice_tuple = tuple([slice(None, None, None) if 2-dir != i else pos_index for i in range(3)])
        self.fileid[self._data_name_in_file].read_direct(self._data, source_sel=slice_tuple)
        self._axis_slices = [slice(self._axis_range[0, i], self._axis_range[1, i], self._cell_size[i]) for i in range(3) if i!=dir]#this can be simplified to self._axis_slices = [self._axis_slices[i] for i in range(3) if i!=dir] but bug test is required
        self._axis_labels = [self._axis_labels_original[i] for i in range(self._num_dimensions) if i!=dir]
        self._axis_units = [self._axis_units_original[i] for i in range(self._num_dimensions) if i!=dir]
        self._fig_title = 't = {0:.2f}, slice at {1} = {2}'.format(self.fileid.attrs['TIME'][0], self._axis_labels_original[dir], pos)

################################method read_data_project################################
    def read_data_project(self, *args, **kwargs):
        '''Automatically detect the dimension of data and choose the proper project method'''
        self.read_data()
        project_methods = (self.data_project2d, self.data_project3d)
        return project_methods[self._data.ndim-2](*args, **kwargs)

################################method data_project3d################################
    def data_project3d(self, dir = 0, if_abs = False):
        '''dir is the direction the summation is taking through'''
        if dir not in range (3):
            raise ValueError('Project direction should be 0, 1 or 2!')
        if if_abs:
            self._data = np.absolute(self._data)
        self._data = np.sum(self._data, axis = 2-dir)/self.fileid[self._data_name_in_file].shape[2-dir]
        self._axis_slices = [self._axis_slices[i] for i in range(3) if i!=dir]
        self._axis_labels = [self._axis_labels_original[i] for i in range(self._num_dimensions) if i!=dir]
        self._axis_units = [self._axis_units_original[i] for i in range(self._num_dimensions) if i!=dir]
        self._fig_title = 't = {0:.2f}, project {1}along {2} direction'.format(self.fileid.attrs['TIME'][0], 'absolute value ' if if_abs else '', self._axis_labels_original[dir])

################################method data_project2d################################
    def data_project2d(self, dir = 0, if_abs = False):
        '''dir is the direction the summation is taking through'''
        if dir not in range (2):
            raise ValueError('Project direction should be 0 or 1!')
        if if_abs:
            self._data = np.absolute(self._data)
        self._data = np.sum(self._data, axis = 1-dir)/self.fileid[self._data_name_in_file].shape[1-dir]
        self._axis_slices = [self._axis_slices[1-dir]]
        self._axis_labels = [self._axis_labels_original[1-dir]]
        self._axis_units = [self._axis_units_original[1-dir]]
        self._fig_title = 't = {0:.2f}, project {1}along {2} direction'.format(self.fileid.attrs['TIME'][0], 'absolute value ' if if_abs else '', self._axis_labels_original[dir])

################################method read_data_lineout################################
    def read_data_lineout(self, dir = 0, pos =(0., 0.)):
        '''dir is the lineout direction'''
        if 2 > self._num_dimensions:
            raise RuntimeError('Method OutFile.read_data_lineout() cannot work on data with dimension number < 2!')
        if dir not in range(self._num_dimensions):
            raise ValueError('Lineout direction should be in range({0})!'.format(self._num_dimensions))
        #get the index of the nearest grid. +0.5 here has a similar effect of rounding.
        tmp_array = np.transpose(np.array([[self._axis_range[0, i], self._cell_size[i]] for i in range(self._num_dimensions) if i!=dir]))
        #bug occurs here if without [:self._num_dimensions-1]: pos_index has two elements, but it is expected to have only one. Now fixed!
        pos_index = ((np.array(pos[:self._num_dimensions-1]) - tmp_array[0, :]) / tmp_array[1, :] + 0.5).astype(np.int, copy=False)
        pos_index = pos_index.tolist()
        tmp_len = [self.fileid[self._data_name_in_file].shape[self._num_dimensions-1-i] for i in range(self._num_dimensions) if i!=dir]
        for i in range(self._num_dimensions-1):
            if tmp_len[i]<=pos_index[i]:
                print('Warning: pos[{0}] is larger than the upper bundary! Force slicing at the upper bundary.'.format(i))
                pos_index[i] = tmp_len[i]-1
            elif 0>pos_index[i]:
                print('Warning: pos[{0}] is smaller than the lower bundary! Force slicing at the lower bundary.'.format(i))
                pos_index[i] = 0
        new_shape = self.fileid[self._data_name_in_file].shape[self._num_dimensions-1-dir]
        self._data = np.zeros(new_shape, dtype=float_type)
        pos_index.insert(dir, slice(None, None, None))
        slice_tuple = tuple(pos_index)[::-1]
        self.fileid[self._data_name_in_file].read_direct(self._data, source_sel=slice_tuple)
        self._axis_slices = [slice(self._axis_range[0, dir], self._axis_range[1, dir], self._cell_size[dir]),]#this can be simplified to self._axis_slices = [self._axis_slices[dir],] but bug test is required
        self._axis_labels = [self._axis_labels_original[dir],]
        self._axis_units = [self._axis_units_original[dir],]
        tmp_list = [self._axis_labels_original[i] for i in range(self._num_dimensions) if i!=dir]
        tmp_list = ['{0} = {1}'.format(tmp_list[i], pos[i]) for i in range(self._num_dimensions-1)]
        self._fig_title = 't = {0:.2f}, lineout at '.format(self.fileid.attrs['TIME'][0]) + ', '.join(tmp_list)

################################method data_FFT1d################################
    def data_FFT1d(self):
        '''Do fast Fourier transform for a 1D data'''
        if self._data.ndim!=1:
            raise RuntimeError('Data is not one dimensional! The OutFile.data_FT1d() method cannot proceed.')
        #self._data = np.absolute(np.fft.rfft(self._data))*(4.0 * pi / self._data.size / 84.7)
        self._data = np.absolute(np.fft.rfft(self._data))*((self._axis_slices[0].stop-self._axis_slices[0].start)/self._axis_slices[0].step/np.sqrt(2.0 * pi))
        delta_frequency = 2.0 * pi/(self._axis_slices[0].stop-self._axis_slices[0].start)
        self._axis_slices = [slice(0, delta_frequency+pi/self._axis_slices[0].step, delta_frequency)]
        self._axis_labels = ['$\omega$', '$E_{y}(\omega)$']
        self._axis_units = ['$\omega_p$', 'arb. unit']

################################method data_STFT################################
    def data_STFT(self):
        '''Do short time Fourier transform for the data'''
        if self._data.ndim!=1:
            raise RuntimeError('Data is not one dimensional! The OutFile.data_STFT() method cannot proceed.')
        one_over_dz = 1./self._axis_slices[0].step
        sample_freq = 2*pi * one_over_dz
        f, t, self._data = signal.spectrogram(self._data, one_over_dz, nperseg = 2*int(one_over_dz), noverlap = int(one_over_dz / 1.))
        self._axis_slices=[slice(t[0]+self._axis_slices[0].start, t[-1]+t[1]-t[0]+self._axis_slices[0].start, t[1]-t[0]), slice(2*pi * f[0], 2*pi * (f[-1]+f[1]-f[0]), 2*pi * (f[1]-f[0]))]
        self._axis_labels = [self._axis_labels_original[0], 'k']
        self._axis_units = [self._axis_units_original[0], '$k_p$']
        self._fig_title = 't = {0:.2f}, spectrogram'.format(self.fileid.attrs['TIME'][0])
        #print(self._axis_slices)
        #print('t={0}:{1}:{2}'.format(t[0],t[-1]+t[1]-t[0],t[1]-t[0]))
        #print(type(Sxx))
        #plt.pcolormesh(t, f, Sxx, cmap=plt.cm.jet)#, shading='gouraud')
        #plt.ylabel('Frequency [Hz]')
        #plt.xlabel('{0} [{1}]'.format(self._axis_labels[0], self._axis_units[0]))
        #color_bar = plt.colorbar()
        #plt.show()

################################method plot_data_line################################
    def plot_data_line(self, h_fig=None, h_ax=None, semilogy=False, linestyle='', if_xlabel=True, if_ylabel=True, if_title=True, multiple=1., offset=0.):
        '''Plot 1D data in a line. h_ax is the handle of the axis to be plotted on. If h_ax=None, a new figure and axis is created.'''
        if self._data.ndim!=1:
            raise RuntimeError('Data is not one dimensional! The OutFile.plot_data_line() method cannot proceed.')
        if h_fig is None:
            h_fig = plt.figure()
        if h_ax is None:
            h_ax = h_fig.add_subplot(111)
        if semilogy:
            plotfunc=h_ax.semilogy
        else:
            plotfunc=h_ax.plot
        plotfunc(np.mgrid[self._axis_slices[0]], (self._data*multiple)+offset, linestyle)
        if if_xlabel:
            h_ax.set_xlabel('{0} [{1}]'.format(self._axis_labels[0], self._axis_units[0]))
        #h_ax.set_ylabel(self._field_names[self._data_name_in_file])
        if if_ylabel:
            if len(self._axis_labels)>=2:
                h_ax.set_ylabel('{0} [{1}]'.format(self._axis_labels[1], self._axis_units[1]))
            #if len(self._axis_labels)<2, calling this will cause an error
        if if_title:
            h_ax.set_title(self._fig_title)
        return h_fig, h_ax

################################method pcolor_data_2d################################
    def pcolor_data_2d(self, h_fig=None, h_ax=None, if_colorbar=True, colorbar_orientation='vertical', if_log_colorbar=False, vmin=None, vmax=None, cmap=plt.cm.jet, alpha=None):
        '''Plot 2D data in as pcolor'''
        if self._data.ndim!=2:
            raise RuntimeError('Data is not two dimensional! The OutFile.pcolor_data_2d() method cannot proceed.')
        y_spread, x_spread = np.mgrid[self._axis_slices[1], self._axis_slices[0]]
        if h_fig is None:
            h_fig = plt.figure()
        if h_ax is None:
            h_ax = h_fig.add_subplot(111)
        if if_log_colorbar:
            h_plot = h_ax.pcolormesh(x_spread, y_spread, np.absolute(self._data), norm=LogNorm(vmin=vmin, vmax=vmax), cmap=cmap, alpha=alpha, antialiased=True)#, shading='gouraud')
        else:
            h_plot = h_ax.pcolormesh(x_spread, y_spread, self._data, vmin=vmin, vmax=vmax, cmap=cmap, alpha=alpha)#, shading='gouraud')
        h_ax.set_xlabel('{0} [{1}]'.format(self._axis_labels[0], self._axis_units[0]))
        h_ax.set_ylabel('{0} [{1}]'.format(self._axis_labels[1], self._axis_units[1]))
        if if_colorbar:
            self._color_bar = plt.colorbar(h_plot, ax=h_ax, orientation=colorbar_orientation)
            self._color_bar.set_label(self._field_names[self._data_name_in_file])
        h_ax.set_title(self._fig_title)
        return h_fig, h_ax

################################method plot_raw_hist_p1################################
    def plot_raw_hist_p1(self, h_fig=None, h_ax=None, num_bins=256, range_max=None, range_min=None):
        '''Plot histogram from p1 raw data.'''
        if h_fig is None:
            h_fig = plt.figure()
        if h_ax is None:
            h_ax = h_fig.add_subplot(111)
        if range_max is None:
            range_max = self._raw_p1.max()
        if range_min is None:
            range_min = self._raw_p1.min()
        '''dn = (range_max - range_min)/num_bins
        hist = np.zeros(num_bins+1, dtype=float_type)
        for i in range(self._raw_p1.size):
            bar_num = int((self._raw_p1[i]-range_min)//dn)
            if bar_num in range(num_bins):
                hist[bar_num]+=abs(self._raw_q[i])'''
        weights=np.absolute(self._raw_q)
        hist, bin_edges = np.histogram(self._raw_p1, num_bins, (range_min, range_max), weights=weights)
        bin_edges = bin_edges[0:-1]
        #h_ax.plot(np.arange(range_min, range_max+dn, dn, dtype=float_type), hist)
        #average = np.average(self._raw_p1, weights=weights)
        p = (np.average((np.absolute(self._raw_p1)), weights=weights))
        h_ax.plot(bin_edges, hist)
        h_ax.set_xlabel('$p_z$ [$m_ec$]')
        h_ax.set_ylabel('Counts [arbt. units]')
        h_ax.set_title('$p_z$ histogram. $\\|p_z\\|$ average = {0:.2f}'.format(p))
        return h_fig, h_ax, bin_edges, hist

################################method plot_raw_hist_gamma################################
    def raw_hist_gamma(self, num_bins=256, range_max=None, range_min=None):
        '''Get histogram from ene raw data +1 = gamma.'''
        if range_max is None:
            range_max = self._raw_ene.max()+1.
        if range_min is None:
            range_min = self._raw_ene.min()+1.
        weights=np.absolute(self._raw_q)
        hist, bin_edges = np.histogram(self._raw_ene+1., num_bins, (range_min, range_max), weights=weights)
        bin_edges = bin_edges[0:-1]
        gamma_ave = (np.average(self._raw_ene+1., weights=weights))
        return bin_edges, hist, gamma_ave

################################method plot_raw_hist_gamma################################
    def plot_raw_hist_gamma(self, h_fig=None, h_ax=None, num_bins=256, range_max=None, range_min=None):
        '''Plot histogram from ene raw data +1 = gamma.'''
        if h_fig is None:
            h_fig = plt.figure()
        if h_ax is None:
            h_ax = h_fig.add_subplot(111)
        bin_edges, hist, gamma_ave = self.raw_hist_gamma(num_bins, range_max, range_min)
        h_ax.plot(bin_edges, hist)
        h_ax.set_xlabel('$\\gamma$')
        h_ax.set_ylabel('Counts [arbt. units]')
        h_ax.set_title('$\\gamma$ histogram. $\\gamma$ average = {0:.2f}'.format(gamma_ave))
        return h_fig, h_ax, bin_edges, hist

################################method plot_data################################
    def plot_data(self, *args, **kwargs):
        '''Automatically detect the dimension of data and choose the proper plot method'''
        plot_methods = (self.plot_data_line, self.pcolor_data_2d)
        return plot_methods[self._data.ndim-1](*args, **kwargs)

################################method fit_for_W################################
    def fit_for_W(self, *args, **kwargs):
        '''Automatically detect the dimension of data and choose the proper fit dimension for W'''
        fit_methods = (self.fit_for_W_1d, self.fit_for_W_2d)
        return fit_methods[self._data.ndim-1](*args, **kwargs)

################################method fit_for_W_2d################################
    def fit_for_W_2d(self, h_fig=None, h_ax=None, guess_values=None):
        '''Do fitting for W. Usually use read_data_project(if_bas=True) to read data absolut value projection along 0th direction before this method.'''
        #self.read_data_project(dir = 0, if_abs = True)
        max_value = np.amax(self._data)
        min_value = np.amin(self._data)
        if guess_values is None:
            #guess_values are [amplitude, x0, y0, sigma, offset]
            guess_values = [max_value-min_value, (self._axis_range[1, 1]+self._axis_range[0, 1])*0.5, (self._axis_range[1, 2]+self._axis_range[0, 2])*0.5, (self._axis_range[1, 1]-self._axis_range[0, 1])*0.2, min_value]
        y_spread, x_spread = np.mgrid[self._axis_slices[1], self._axis_slices[0]]
        try:
            popt, pcov = tdgf.Fit2DGauss_simple(x_spread, y_spread, self._data, guess_values)
        except RuntimeError:
            guess_values[0] = -guess_values[0]
            guess_values[-1] = max_value
            popt, pcov = tdgf.Fit2DGauss_simple(x_spread, y_spread, self._data, guess_values)
        self._W = popt[3]*np.sqrt(2)
        self._a = popt[0]-popt[4]
        #plot contour
        if h_fig is None:
            h_fig = plt.figure()
        if h_ax is None:
            h_ax = h_fig.add_subplot(111)
        self.pcolor_data_2d(h_fig, h_ax)
        h_plot = h_ax.contour(x_spread, y_spread, self._data, [popt[0]/const_e+popt[-1]], colors='w')
        theta = np.arange(0., 2*pi, 0.01)
        h_ax.plot(self._W*np.cos(theta)+popt[1], self._W*np.sin(theta)+popt[2], 'k')
        h_ax.contour(x_spread, y_spread, self._data, [popt[0]/(const_e**0.25)+popt[-1]], colors='w')
        theta = np.arange(0., 2*pi, 0.01)
        h_ax.plot(self._W*np.cos(theta)*0.5+popt[1], self._W*np.sin(theta)*0.5+popt[2], 'k')
        return popt, h_fig, h_ax

################################method fit_for_W_1d################################
    def fit_for_W_1d(self, h_fig=None, h_ax=None, guess_values=None):
        '''Do fitting for W in quasi-3D modes osiris simulation.'''
        #guess_values are [amplitude, sigma]
        if guess_values is None:
            max_value = np.amax(self._data)
            min_value = 0.
            guess_values = [max_value-min_value, (self._axis_range[1, 1]-self._axis_range[0, 1])*0.2]
        r_spread = np.mgrid[self._axis_slices[0]]
        popt, pcov = tdgf.Fit_Gauss_simple(r_spread, self._data, guess_values)
        self._W = popt[1]*np.sqrt(2)
        self._a = popt[0]
        #plot contour
        if h_fig is None:
            h_fig = plt.figure()
        if h_ax is None:
            h_ax = h_fig.add_subplot(111)        
        self.plot_data(h_fig, h_ax)
        h_plot = h_ax.plot(r_spread, tdgf.Gaussian_simple(r_spread,popt[0],popt[1]), 'r--')
        return popt, h_fig, h_ax

################################method start_index################################
    def start_index(self, limit = 1.e-2):
        '''Find the start index when the data abslute value > limit.'''
        if self._data.ndim!=1:
            raise RuntimeError('Data is not one dimensional! The OutFile.start_index() method cannot proceed.')
        for i in range(self._data.size-1, 0, -1):
            if np.absolute(self._data[i])>limit:
                return i
        raise RuntimeError('The OutFile.start_index() method cannot find a start index with abs(self._data[index])>limit.')

################################method next_local_max################################
    def next_local_max(self, start_index=None, local_range = 6.28):
        '''Find the next local maximum point for a 1D data from the start_index to 0 index. If the data is not 1D, raise an exception. local_range is the length that local maximum is taking.'''
        if self._data.ndim!=1:
            raise RuntimeError('Data is not one dimensional! The OutFile.next_local_max() method cannot proceed.')
        if None==start_index:
            start_index = self._data.size
        till_index = start_index - int(local_range / self._axis_slices[0].step)
        return_ind=np.argmax(self._data[till_index:start_index])+till_index
        return return_ind

################################method next_local_min################################
    def next_local_min(self, start_index=None, local_range = 6.28):
        '''Find the next local minimum point for a 1D data from the start_index to 0 index. If the data is not 1D, raise an exception. local_range is the length that local minimum is taking.'''
        if self._data.ndim!=1:
            raise RuntimeError('Data is not one dimensional! The OutFile.next_local_min() method cannot proceed.')
        if None==start_index:
            start_index = self._data.size
        till_index = start_index - int(local_range / self._axis_slices[0].step)
        return_ind=np.argmin(self._data[till_index:start_index])+till_index
        return return_ind
#not completed

################################method zero_point################################
    def zero_point(self, start, stop):
        '''Find the zero point for a 1D data from index of start to stop. If the data is not 1D, raise an exception. If self._data[start] andself._data[stop] have the same sign, raise an exception.'''
        if self._data.ndim!=1:
            raise RuntimeError('Data is not one dimensional! The OutFile.zero_poin() method cannot proceed.')
        if np.sign(self._data[start])*np.sign(self._data[stop])>0:
            raise RuntimeError('self._data[start] andself._data[stop] have the same sign! The OutFile.zero_poin() method cannot take a zero point between them.')
        ind1 = start
        ind2 = stop
        while ind2-ind1 > 1:
            ind_mid = (ind1 + ind2) // 2
            if np.sign(self._data[ind1])*np.sign(self._data[ind_mid]) > 0:
                ind1 = ind_mid
            else:
                ind2 = ind_mid
        #do linear interpolation between ind1 and ind2 for the zero point location
        x1 = self._axis_slices[0].start + self._axis_slices[0].step*ind1
        x2 = x1 + self._axis_slices[0].step#or x2 = self._axis_slices[0].start + self._axis_slices[0].setp*ind2
        y1 = self._data[ind1]
        y2 = self._data[ind2]
        return (x1*y2-x2*y1)/(y2-y1)
#not completed

if __name__ == '__main__':
    out_num=990
    def tmp_3D():
        h_fig = plt.figure(figsize=(10,8))
        file1 = OutFile(path='/home/zming/simulations/os2D/Hi_beam3D52',field_name='charge',average='',spec_name='beam',out_num=out_num)
        h_ax = h_fig.add_subplot(221)
        file1.open()
        file1.read_data_slice(dir=1)
        file1.plot_data(h_fig, h_ax, vmax=0., vmin=-1.)
        file1.close()
        file1.spec_name='plasma'
        h_ax = h_fig.add_subplot(222)
        file1.open()
        file1.read_data_slice(dir=1)
        file1.plot_data(h_fig, h_ax)#, vmax=5., vmin=-5.)
        file1.close()

        file1.field_name='raw'
        file1.spec_name='driver'
        file1.open()
        file1.read_raw_q()
        print('driver charge = {} pC'.format(file1.calculate_q_pC(4.9e16, 10./256, 8./256, 8./256)))
        file1.read_raw_p1()
        file1.read_raw_p2()
        file1.read_raw_p3()
        file1.read_raw_ene()
        h_ax = h_fig.add_subplot(223)
        file1.plot_raw_hist_gamma(h_fig, h_ax)
        file1.close()

        file1.field_name='raw'
        file1.spec_name='trailer'
        file1.open()
        file1.read_raw_q()
        print('trailer charge = {}'.format(file1.calculate_q_pC(4.9e16, 10./256, 8./256, 8./256)))
        file1.read_raw_p1()
        file1.read_raw_p2()
        file1.read_raw_p3()
        file1.read_raw_ene()
        h_ax = h_fig.add_subplot(224)
        file1.plot_raw_hist_gamma(h_fig, h_ax)
        file1.close()
        plt.tight_layout()
        plt.show()
    def tmp_2D():
        h_fig = plt.figure(figsize=(16.5,11))
        file1 = OutFile(path='/home/zming/simulations/os2D/os_res2D81',field_name='charge',spec_name='e_He',out_num=out_num)
        h_ax = h_fig.add_subplot(221)
        file1.open()
        file1.read_data()
        file1.plot_data(h_fig, h_ax)#, if_log_colorbar=True)
        file1.close()
        h_ax = h_fig.add_subplot(222)
        file1.spec_name='e_He'
        #file1.field_name='p1x1'
        file1.open()
        file1.read_data()
        file1.plot_data(h_fig, h_ax)#, if_log_colorbar=True)#, vmax=1.)
        file1.close()
        h_ax = h_fig.add_subplot(223)
        file1.field_name = 'e1'
        file1.open()
        file1.read_data()
        file1.plot_data(h_fig, h_ax)#, vmax=3., vmin=-3.)
        file1.close()
        '''
        h_ax = h_fig.add_subplot(224)
        file1.field_name = 'psi'
        file1.open()
        file1.read_data()
        print('Delta psi = {0}'.format(file1._data.max()-file1._data.min()))
        file1.plot_data(h_fig, h_ax)#, vmax=5., vmin=-1.)
        file1.close()'''
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
    def tmp_1D():
        file1 = OutFile(path='/home/zming/simulations/os2D/os_laser1D5',field_name='e3',spec_name='e',out_num=out_num)
        file1.open()
        file1.read_data()
        h_fig = plt.figure()
        h_ax = h_fig.add_subplot(311)
        file1.plot_data(h_fig, h_ax)
        #plt.show(block=True)
        file1.data_STFT()
        h_ax = h_fig.add_subplot(312)
        file1.plot_data(h_fig, h_ax, if_log_colorbar=True)
        file1.close()
        file1.field_name = 'e1'
        file1.open()
        file1.read_data()
        h_ax = h_fig.add_subplot(313)
        file1.plot_data(h_fig, h_ax)
        file1.close()
        plt.show(block=True)
    def test():
        semilogy=True
        h_fig = plt.figure(figsize=(10,8))
        file1 = OutFile(path='/home/zming/simulations/os2D/os_PT3D3',field_name='e3',average='-savg',spec_name='plasma',out_num=out_num)
        h_ax = h_fig.add_subplot(111)
        file1.open()
        file1.read_data_lineout()
        file1.data_FFT1d()
        file1.plot_data(h_fig, h_ax, semilogy=semilogy, linestyle='k-')
        file1.close()
        file1.out_num=10
        file1.open()
        file1.read_data_lineout()
        file1.data_FFT1d()
        file1.plot_data(h_fig, h_ax, semilogy=semilogy, linestyle='r-')
        file1.close()
        file1.out_num=20
        file1.open()
        file1.read_data_lineout()
        file1.data_FFT1d()
        file1.plot_data(h_fig, h_ax, semilogy=semilogy, linestyle='b-')
        file1.close()
        plt.tight_layout()
        plt.show()
        
    #test()
    #for out_num in range(43,88):
    #    tmp_3D()
    #    plt.savefig('/home/zming/simulations/os2D/os_PT3D4/plots/{}.png'.format(out_num))
    #    plt.close()
    #cyl_m_2D()
    tmp_3D()

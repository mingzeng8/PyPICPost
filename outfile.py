import numpy as np
import h5py
import matplotlib.pyplot as plt
import TwoDGaussianFit as tdgf
from scipy import signal, pi
from scipy import e as const_e
from scipy.constants import e as e_charge
from scipy.signal import find_peaks
from matplotlib.colors import LogNorm
from glob import glob
import os, warnings
try: import parse
except ImportError: warnings.warn('Cannot import lib \'parse\'. num_list cannot be activated.')
import my_cmap

float_type=np.float64

class OutFile:
    def __init__(self, code_name = 'osiris', path = '.', out_type = None, field_name = 'e3', spec_name = '', use_num_list = False, out_num = 0, average='', fld_slice=None, cyl_m_num = 0, cyl_m_re_im='re'):
##value digit_num##
        self._accepted_code_name = {'osiris', 'quickpic', 'hipace'}
        self.code_name = code_name.lower()
        self.path = path
        self.spec_name = spec_name
        self._field_name_to_out_type = {'psi':'FLD', 'e1':'FLD', 'e2':'FLD', 'e3':'FLD', 'e3_cyl_m':'FLD_CYL_M', 'b1':'FLD', 'b2':'FLD', 'b3':'FLD', 'j1':'DENSITY', 'ene':'DENSITY', 'charge':'DENSITY', 'ion_charge':'ION', 'p1x1':'PHA', 'p2x2':'PHA', 'raw':'RAW', 'ExmBy':'FLD', 'Ez':'FLD', 'EypBx':'FLD', 'tracks':'TRACKS'}
        #self.out_type = out_type
        self.field_name = field_name
        self.average = average
        # fld_slice can be 1, 2 or 3 (for QuickPIC)
        self.fld_slice = fld_slice
        #2D cylindrical modes for fields
        self.cyl_m_num = cyl_m_num
        self.cyl_m_re_im = cyl_m_re_im

        #self._num_dimensions = 0
        ''' This is moved to self.set_code_name() so that self._axis_labels_original is modified every time code name is changed.
        if self.code_name in {'quickpic', 'hipace'}: self._axis_labels_original = ('$\\xi$', 'x', 'y', '$p_z$', '$p_x$', '$p_y$')
        else: self._axis_labels_original = ('z', 'x', 'y', '$p_z$', '$p_x$', '$p_y$')'''
        self._axis_units_original = ('$c / \\omega_p$',)*3 + ('$m_ec$',)*3
        #self._axis_units_original = ('$k_p^{-1}$',)*3 + ('$m_ec$',)*3
        self._field_names = {'psi':'$\psi$', 'e1':'$E_z$', 'e2':'$E_x$', 'e3':'$E_y$', 'e3_cyl_m':'$E_y$', 'b1':'$B_z$', 'b2':'$B_x$', 'b3':'$B_y$', 'j1':'$J_z$', 'ene':'$E_k / n_p m_e c^2$', 'charge':'$\\rho$', 'ion_charge':'$\\rho$', 'p1x1':'$p_1x_1$ [arb. units]', 'p2x2':'$p_2x_2$ [arb. units]', 'beam_charge':'$\\rho_b$', 'plasma_charge':'$\\rho_e$', 'ExmBy':'$E_x-B_y$', 'Ez':'$E_z$', 'EypBx':'$E_y+B_x$'\
        # Some field naming problem in new HiPACE
        , 'plasma_electrons': '$\\rho_e$', 'driver': '$\\rho_d$', 'drive_beam': '$\\rho_d$', 'trailer': '$\\rho_t$'\
        # Some field naming in QuickPIC
        , 'ez':'$E_z$', 'ex':'$E_x$', 'ey':'$E_y$', 'bz':'$B_z$', 'bx':'$B_x$', 'by':'$B_y$', 'charge_slice_xz':'$\\rho$'}
        self._accepted_out_type = {'DENSITY', 'FLD', 'FLD_CYL_M', 'ION', 'PHA', 'RAW', 'TRACKS'}
        # if use_num_list = True, the actual out_num used is self._avail_num_list[out_num]
        self.use_num_list = use_num_list
        self.out_num = out_num


################################property code_name################################
    def get_code_name(self):
        return self._code_name

    def set_code_name(self, value):
        if value not in self._accepted_code_name:
            raise ValueError('code_name \'{0}\' not implemented!'.format(value))
        self._code_name = value
        if value=='quickpic': self._axis_labels_original = ('$\\xi$', 'x', 'y', '$p_z$', '$p_x$', '$p_y$')
        elif value=='hipace': self._axis_labels_original = ('$\\zeta$', 'x', 'y', '$p_z$', '$p_x$', '$p_y$')
        else: self._axis_labels_original = ('z', 'x', 'y', '$p_z$', '$p_x$', '$p_y$')
        if 'quickpic' == self.code_name: self.digit_num = 8
        else: self.digit_num = 6
        ''' Deprecated
        # Set the AXIS group and data matrix axis sequence of the h5 file for different codes
        if value in {'osiris'}:
            # OSIRIS 3D h5 file has data axes [y,x,z] and AXIS group [z,x,y].
            # OSIRIS 2D h5 file has data axes [x,z] and AXIS group [z,x].
            self._h5_data_seq = [2,1,0]
            self._h5_AXIS_seq = [0,1,2]
        elif value in {'quickpic'}:
            # QuickPIC h5 file has data axes [z,y,x] and AXIS group [x,y,z].
            self._h5_data_seq = [0,2,1]
            self._h5_AXIS_seq = [2,0,1]
        elif value in {'hipace'}:
            # HiPACE h5 file has data axes [z,x,y] and AXIS group [z,x,y].
            self._h5_data_seq = [0,1,2]
            self._h5_AXIS_seq = [0,1,2]'''

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

################################property num_of_zeros################################
    def get_num_of_zeros(self):
        return self._num_of_zeros

    def set_num_of_zeros(self, **kwargs):
        raise RuntimeError('Please do not set num_of_zeros directly!. It will be automatically set when setting out_num.')
    num_of_zeros = property(get_num_of_zeros, set_num_of_zeros)

################################property use_num_list################################
    def get_use_num_list(self):
        return self._use_num_list

    def set_use_num_list(self, value):
        if isinstance(value, bool):
            self._use_num_list = value
            if value: self.reset_avail_num_list()
            else: pass #maybe del self._avail_num_list is better
        else: raise TypeError("set_use_num_list should recieve a boolean type!")

    use_num_list = property(get_use_num_list, set_use_num_list)

################################property out_num################################
    def get_out_num(self):
        return self._out_num

    def set_out_num(self, value):
        if self.use_num_list:
            #try: self._avail_num_list
            #except: self.reset_avail_num_list()
            if value not in range(len(self._avail_num_list)):
                #print('self._avail_num_list is ({})'.format(self._avail_num_list))
                raise KeyError('Using num_list. out_num = {}, not in range({})!'.format(value, len(self._avail_num_list)))
            self.actual_num = self._avail_num_list[value]
        else:
            if value not in range(int(10**self.digit_num)):
                raise ValueError('out_num \'{0}\' not in range!'.format(value))
            self.actual_num = value
        self._out_num = value
##value num_of_zeros##
        if 0 == self.actual_num:
            self._num_of_zeros = self.digit_num - 1
        else:
            self._num_of_zeros = self.digit_num - 1 - int(np.log10(self.actual_num))

    out_num = property(get_out_num, set_out_num)

################################property average################################
    def get_average(self):
        return self._average

    def set_average(self, value):
        if value not in ('','-savg'):
            raise NotImplementedError('average = \'{0}\' not implemented!'.format(value))
        self._average = value

    average = property(get_average, set_average)

################################ prefix_filename ################################
    def get_prefix_filename(self):
        self.reset_prefix_filename()
        return self._prefix_filename

    def reset_prefix_filename(self):
        if 'osiris' == self.code_name:
            main_folder_path = '{}/MS/{}'.format(self.path, self.out_type)
            if 'DENSITY' == self._out_type:
                self._prefix_filename = '{0}/{1}/{2}{3}/{2}{3}-{1}-'.format(main_folder_path, self.spec_name, self.field_name, self.average)
            elif 'ION' == self._out_type:
                self._prefix_filename = '{0}/{1}/{2}/{2}-{1}-'.format(main_folder_path, self.spec_name, self.field_name)
            elif 'FLD' == self._out_type:
                self._prefix_filename = '{0}/{1}{2}/{1}{2}-'.format(main_folder_path, self.field_name, self.average)
            elif 'FLD_CYL_M' == self._out_type:
                self._prefix_filename = '{0}/MODE-{1}-{2}/{3}/{3}-{1}-{4}-'.format(main_folder_path, self.cyl_m_num, self.cyl_m_re_im.upper(), self.field_name, self.cyl_m_re_im.lower())
            elif 'PHA' == self._out_type:
                self._prefix_filename = '{0}/{1}/{2}/{1}-{2}-'.format(main_folder_path, self.field_name, self.spec_name)
            elif 'RAW' == self._out_type:
                self._prefix_filename = '{0}/{1}/RAW-{1}-'.format(main_folder_path, self.spec_name)
            elif 'TRACKS' == self._out_type:
                self._prefix_filename = '{0}/{1}-tracks.h5'.format(main_folder_path, self.spec_name)
        elif 'quickpic' == self.code_name:
            field_name_dict = {'charge':'Charge', 'e1':'Ez', 'e2':'Ex', 'e3':'Ey', 'b1':'Bz', 'b2':'Bx', 'b3':'By', 'psi':'Psi'}
            # For slice outputs of quickpic
            slc_dic = {1:'xz', 2:'yz', 3:'xy'}
            if 'DENSITY' == self._out_type:
                if self.fld_slice is None:
                    fld_name1 = field_name_dict[self.field_name]
                    fld_name2 = self.field_name
                else:
                    fld_name1 = '{}_slice_000{}'.format(field_name_dict[self.field_name], self.fld_slice)
                    fld_name2 = '{}_slice_{}'.format(self.field_name, slc_dic[self.fld_slice])
                self._prefix_filename = '{0}/{1}/{2}/{3}_'.format(self.path, self.spec_name, fld_name1, fld_name2)
            elif 'FLD' == self._out_type:
                if self.fld_slice is None:
                    fld_name1 = field_name_dict[self.field_name]
                    fld_name2 = field_name_dict[self.field_name].lower()
                else:
                    fld_name1 = '{}_slice000{}'.format(field_name_dict[self.field_name], self.fld_slice)
                    fld_name2 = '{}slice{}'.format(field_name_dict[self.field_name].lower(), slc_dic[self.fld_slice])
                self._prefix_filename = '{0}/Fields/{1}/{2}_'.format(self.path, fld_name1, fld_name2)
            elif 'RAW' == self._out_type:
                self._prefix_filename = '{0}/{1}/Raw/raw_'.format(self.path, self.spec_name)
        elif 'hipace' == self.code_name:
            main_folder_path = '{}/DATA'.format(self.path)
            if 'DENSITY' == self._out_type:
                #self._prefix_filename = '{0}/density_{1}_{2}_'.format(main_folder_path, self.spec_name, self.field_name)
                # Changed according to new HiPACE
                self._prefix_filename = '{0}/density_{1}_'.format(main_folder_path, self.spec_name)
            elif 'FLD' == self._out_type:
                # HiPACE has different field names
                hi_fld_dict = {'e1':'Ez', 'e2':'Ex', 'e3':'Ey', 'psi':'Psi'}
                self._prefix_filename = '{0}/field_{1}_'.format(main_folder_path, hi_fld_dict[self.field_name])
            elif 'RAW' == self._out_type:
                self._prefix_filename = '{0}/raw_{1}_'.format(main_folder_path, self.spec_name)
        else:
            raise NotImplementedError('Code name {} not implemented!'.format(self.code_name))

################################ avail_num_list ################################
# obtain the number list from the simulation output folder
    def get_avail_num_list(self):
        self.reset_avail_num_list()
        return self._avail_num_list

    def reset_avail_num_list(self):
        self.reset_prefix_filename()
        format_string=self._prefix_filename+'{}.h5'
        file_list = glob(format_string.format('*'))
        try: num_list=[int(parse.parse(format_string, file_list[i])[0]) for i in range(len(file_list))]
        except ValueError: num_list=[float(parse.parse(format_string, file_list[i])[0]) for i in range(len(file_list))]
        if len(num_list)<1: raise RuntimeError('num_list is empty! Check the file path. format_string = {}'.format(format_string))
        num_list.sort()
        self._avail_num_list = num_list

################################property path_filename################################
    def get_path_filename(self):
        tmp_filename = self.get_prefix_filename()
        if 'TRACKS' != self._out_type:
            tmp_filename = '{}{}{}.h5'.format(tmp_filename, '0'*self.num_of_zeros, self.actual_num)
        return tmp_filename

    def set_path_filename(self):
        raise RuntimeError('You cannot set path_filename directly!')

    path_filename = property(get_path_filename, set_path_filename)

################################property num_dimensions################################
    def get_num_dimensions(self):
        return self._num_dimensions

    def set_num_dimensions(self, value):
        raise RuntimeError('You cannot set num_dimensions directly! This number is set by self.open() automatically.')

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
        return self._time

    def set_time(self, value):
        raise RuntimeError('You cannot set time directly!')

    time = property(get_time, set_time)

################################method open################################
    def open(self, filename=None, t_offset = None, cell_size_qp_raw = None):
        '''
            The normal use is to let the code determin file name according to the PIC code's file name regulation. In some cases one may need to set file name manually, so one may provide filename in the argument.
        filename: None or string.
                  The path to an output h5 file. If None, automatically set the path.
        t_offset: None or float.
                  If not None, self._time will have an offset respect to the time read from attribute 'TIME' of the h5 file.
        cell_size_qp_raw: None or a array of 3 elements.
                  The cell size, if the code is QuickPIC and the file type is RAW. Currently QuickPIC RAW output does not have attributes for the cell sizes.
        '''
        if self._out_type not in self._accepted_out_type: raise NotImplementedError('Out type {} not implemented!'.format(self._out_type))
        if filename is None: filename=self.path_filename
        if not os.path.isfile(filename): raise FileNotFoundError('File {} not found!'.format(filename))
        self.fileid = h5py.File(filename,'r')
        # Do not obtain the grid information if the type is TRACKS
        if self._out_type == 'TRACKS': pass
        else:
        # Obtain the grid information
##value _num_dimensions##
            if self._out_type != 'RAW':
                #try to read the dimension from 'AXIS'
                #for the code osiris, when plotting phasespace, the data matrix size may be different from the attribute NX
                self._num_dimensions = len(self.fileid['AXIS'].keys())
                ''' Deprecated
                if self.code_name in {'osiris'}:
                    # For OSIRIS, modifiy self._h5_data_seq for 2d and 1d cases
                    if 2==self._num_dimensions: self._h5_data_seq = [1,0]
                    elif 1==self._num_dimensions: self._h5_data_seq = [0]'''
##value _axis_range##
                self._axis_range = np.zeros((2,self._num_dimensions), dtype=float_type)
                self._axis_labels = []
                self._axis_units = []
                for i in range(self._num_dimensions):
                    self.fileid['AXIS/AXIS{0}'.format(i+1)].read_direct(self._axis_range, np.s_[:], np.s_[:,i])
                    self._axis_labels.append('${}$'.format(self.fileid['AXIS/AXIS{0}'.format(i+1)].attrs.get('NAME')[0].decode("utf-8")))
                    self._axis_units.append('${}$'.format(self.fileid['AXIS/AXIS{0}'.format(i+1)].attrs.get('UNITS')[0].decode("utf-8")))
            else:
                if self.code_name != 'quickpic':
                    # For RAW data, 'AXIS' does not exist. Read simulation box information from the attribute
                    xmax = self.fileid.attrs.get('XMAX')
                    xmin = self.fileid.attrs.get('XMIN')
                    nx = self.fileid.attrs.get('NX')
                    self._num_dimensions = len(xmax)
##value _axis_range##
                    self._axis_range = np.zeros((2,self._num_dimensions), dtype=float_type)
                    for i in range(self._num_dimensions):
                        self._axis_range[0,i] = xmin[i]
                        self._axis_range[1,i] = xmax[i]
                else: self._num_dimensions = 3
            for i in self.fileid.keys():
                if i!='AXIS':
                    break
##value _data_name_in_file##
            self._data_name_in_file = i
            try:
                nx
            except NameError:
                #if nx not defined, get nx from the size of the matrix. In case of code osiris, the data matrix size may be different from the attribute when plotting phasespace.
                nx = np.flip(self.fileid[self._data_name_in_file].shape)
##value _cell_size##
            if self.code_name == 'quickpic' and self._out_type == 'RAW': self._cell_size = cell_size_qp_raw
            else: self._cell_size = (self._axis_range[1, :] - self._axis_range[0, :]) / nx
            self._time = self.fileid.attrs['TIME'][0]
            if t_offset is not None: self._time += t_offset

################################method close################################
    def close(self):
        self.fileid.close()

################################method read_raw_tag################################
    def read_raw_tag(self):
        '''In osiris the tag for one macro particle has two numbers:
           the first is the process index, and the second is the particle index.'''
        self._raw_tag = np.zeros(self.fileid['tag'].shape, dtype=int)
        self.fileid['tag'].read_direct(self._raw_tag)

################################method read_raw_q################################
    def read_raw_q(self):
        self._raw_q = np.zeros(self.fileid['q'].shape, dtype=float_type)
        self.fileid['q'].read_direct(self._raw_q)
        return self._raw_q

################################method read_raw_x1################################
    def read_raw_x1(self):
        if self.code_name == 'quickpic': self._raw_x1 = self.fileid['x3'][()]
        else: self._raw_x1 = self.fileid['x1'][()]
        return self._raw_x1

################################method read_raw_x2################################
    def read_raw_x2(self):
        if self.code_name == 'quickpic': self._raw_x2 = self.fileid['x1'][()]
        else: self._raw_x2 = self.fileid['x2'][()]
        return self._raw_x2

################################method read_raw_x3################################
    def read_raw_x3(self):
        if self.code_name == 'quickpic': self._raw_x3 = self.fileid['x2'][()]
        else: self._raw_x3 = self.fileid['x3'][()]
        return self._raw_x3

################################method read_raw_p1################################
    def read_raw_p1(self):
        if self.code_name == 'quickpic': self._raw_p1 = self.fileid['p3'][()]
        else: self._raw_p1 = self.fileid['p1'][()]
        return self._raw_p1

################################method read_raw_p2################################
    def read_raw_p2(self):
        if self.code_name == 'quickpic': self._raw_p2 = self.fileid['p1'][()]
        else: self._raw_p2 = self.fileid['p2'][()]
        return self._raw_p2

################################method read_raw_p3################################
    def read_raw_p3(self):
        if self.code_name == 'quickpic': self._raw_p3 = self.fileid['p2'][()]
        else: self._raw_p3 = self.fileid['p3'][()]
        return self._raw_p3

################################method read_raw_ene################################
    def read_raw_ene(self, ene_key_warning=True):
        try:
            self._raw_ene = np.zeros(self.fileid['ene'].shape, dtype=float_type)
            self.fileid['ene'].read_direct(self._raw_ene)
        except KeyError:
            if ene_key_warning:
            #by default, if 'ene' key does not exist, print the following warning message.
            #but one can explicitly silence this message by setting ene_key_warning=False
                warnings.warn('Warning! Key \'ene\' does not exist in particle raw data! Reading p1, p2, p3 and doing ene=sqrt(p1^2+p2^2+p3^2)-1 instead. Please make sure p1, p2, p3 are read before this!')
            self._raw_ene = np.sqrt(np.square(self._raw_p1)+np.square(self._raw_p2)+np.square(self._raw_p3))-1.
        return self._raw_ene

################################method select_raw_data################################
# Select particles according to the raw data
# x1_low, x1_up: lower and upper limits of x1.
# ...
# r_low, r_up: lower and upper limits of radius of position, i.e. sqrt(x2^2+x3^2).
# sample_size: select a random sample of size sample_size from selected macro particles according to the former conditions. sample_size overrides sample_rate.
# sample_rate: a number between 0 and 1. sample_size = int(sample_rate * number of remaining particles from former range selections). If sample_size is not None, this parameter is ignored.
# set self._raw_select_index with the selection index array, and also return this array.
    def select_raw_data(self, x1_low=None, x1_up=None, x2_low=None, x2_up=None, x3_low=None, x3_up=None, p1_low=None, p1_up=None, p2_low=None, p2_up=None, p3_low=None, p3_up=None, ene_low=None, ene_up=None, r_low=None, r_up=None, sample_size=None, sample_rate=None):
        #get number of particles
        try: n_part = self._raw_tag.size // 2
        except AttributeError:
            try: n_part = self._raw_q.size
            except AttributeError:
                try: n_part = self._raw_x1.size
                except AttributeError:
                    try: n_part = self._raw_p1.size
                    except AttributeError:
                        try: n_part = self._raw_x2.size
                        except AttributeError:
                            try: n_part = self._raw_p2.size
                            except AttributeError:
                                try: n_part = self._raw_x3.size
                                except AttributeError:
                                    try: n_part = self._raw_p3.size
                                    except AttributeError:
                                        n_part = self._raw_ene.size
        #select_list is a boolean list, if one of its element is True, the corresponding macro particle is selected
        select_list = np.full(n_part, True, dtype=bool)
        if x1_low is not None:
            select_list = self._raw_x1 > x1_low
        if x1_up is not None:
            select_list = (self._raw_x1 < x1_up) & (select_list)
        if x2_low is not None:
            select_list = (self._raw_x2 > x2_low) & (select_list)
        if x2_up is not None:
            select_list = (self._raw_x2 < x2_up) & (select_list)
        if x3_low is not None:
            select_list = (self._raw_x3 > x3_low) & (select_list)
        if x3_up is not None:
            select_list = (self._raw_x3 < x3_up) & (select_list)
        if r_low is not None:
            r_square = np.square(self._raw_x3)+np.square(self._raw_x2)
            select_list = (r_square > (r_low*r_low)) & (select_list)
        if r_up is not None:
            try: r_square
            except UnboundLocalError: r_square = np.square(self._raw_x3)+np.square(self._raw_x2)
            select_list = (r_square < (r_up*r_up)) & (select_list)

        if p1_low is not None:
            select_list = (self._raw_p1 > p1_low) & (select_list)
        if p1_up is not None:
            select_list = (self._raw_p1 < p1_up) & (select_list)
        if p2_low is not None:
            select_list = (self._raw_p2 > p2_low) & (select_list)
        if p2_up is not None:
            select_list = (self._raw_p2 < p2_up) & (select_list)
        if p3_low is not None:
            select_list = (self._raw_p3 > p3_low) & (select_list)
        if p3_up is not None:
            select_list = (self._raw_p3 < p3_up) & (select_list)

        if ene_low is not None:
            select_list = (self._raw_ene > ene_low) & (select_list)
        if ene_up is not None:
            select_list = (self._raw_ene < ene_up) & (select_list)

        #debug: for predictable samples
        np.random.seed(0)
        #select_index_array is an array containing the index of selected macro particles
        #np.nonzero() returns a tuple.
        select_index_array = np.nonzero(select_list)[0]
        if sample_rate is not None:
            if sample_size is None: sample_size = int(np.size(select_index_array) * sample_rate)
        if sample_size is not None:
            select_index_array = np.random.choice(select_index_array, sample_size, replace=False)
        self._raw_select_index = select_index_array
        return select_index_array

################################method calculate_q_pC################################
# calculate the charge in unit of pico-coulumb
# please call read_raw_q() before this function
# if_select = False: use all the macro particles
# if_select = True: only use the macro particles with index in self._raw_select_index
    def calculate_q_pC(self, n0_per_cc, if_select = False):
        '''
        n0_per_cc: reference density in simulation, in unit of per centimeter cube
        when if_select = True, only calculate the selected macro particles according to self._raw_select_index.
        The formula for charge in pC: Q[pC] = sum_q*n0*V_cell_norm*k0**-3*e_charge*10**12
        and k0=(4*pi*r_e*n0)**0.5, where r_e is classical electron radius.
        So Q[pC] = sum_q*V_cell_norm*(4pi*r_e)**-1.5*e_charge*10**12/n0**0.5,
        where n0 in unit of cm**-3.
        e_charge*(4pi*r_e)**-1.5*10**12 = 24043512116.12064
        or (4pi*r_e)**-1.5*10**12 = 1.5006780029105165e+29
        Finally Q[pC] = sum_q*V_cell_norm/n0**0.5*24043512116.12064
        '''
        q_array = self._raw_q
        if if_select:
            try:
                if 0==len(self._raw_select_index):
                    print("Warning: no particle is selected! Charge is set to 0.")
                    return 0.0
                q_array = q_array[self._raw_select_index]
            except: warnings.warn('Particle select condition is not valid! All particles are used.')
        #normalized cell volume
        cell_volume_norm = 1.0
        for i in range(self.num_dimensions):
            cell_volume_norm = cell_volume_norm*self._cell_size[i]
        if 3>self.num_dimensions:
            print('Warning! Similation is in {} dimensional. Charge calculation may not be correct.'.format(self.num_dimensions))
        return np.sum(q_array)*cell_volume_norm/np.sqrt(n0_per_cc)*24043512116.12064

################################method calculate_norm_rms_emittance_um################################
# calculate the normalized rms emittance in unit of micrometer radian
# please call read_raw_x2(), read_raw_p2() and/or read_raw_x3(), read_raw_p3()
# and/or read_raw_x1(), read_raw_p1()
# and read_raw_q() before this function
    def calculate_norm_rms_emittance_um(self, n0_per_cc, directions=(2,), if_select = False):
        '''
        n0_per_cc: reference density in simulation, in unit of per centimeter cube
        one_over_k0_um: one over k0, in unit of micrometer
        one_over_k0_um = (4pi*r_e_in_um*n0_per_cc*10^-12)^-0.5 = (3.541e-20*n0_per_cc)^-0.5
        when if_select = True, only calculate the selected macro particles according to self._raw_select_index.
        return: 2 elements. The first is normalized emittance list in each of the given directions. The second is a list of 3 elements: the 3 Courant-Snyder parameters in each of the given directions.
        return all 0 if no particle is selected.
        '''
        one_over_k0_um_square = 2.8239587227915743e19/n0_per_cc
        one_over_k0_mm = np.sqrt(one_over_k0_um_square)/1.e3
        norm_emittances=np.zeros(len(directions))

        # Courant-Snyder parameters
        alphas=np.zeros(len(directions))
        betas=np.zeros(len(directions))
        gammas=np.zeros(len(directions))

        weights = np.absolute(self._raw_q)
        p1_array=self._raw_p1
        if if_select:
            try:
                if 0==len(self._raw_select_index):
                    print("Warning: no particle is selected! Emittance is set to 0.")
                    return emittances, [alphas, betas, gammas]
                weights = weights[self._raw_select_index]
                p1_array = p1_array[self._raw_select_index]
            except: warnings.warn('Particle select condition is not valid! All particles are used.')
        sum_weight=np.sum(weights)
        mean_p1=np.sum(np.multiply(p1_array, weights))/sum_weight
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
            if if_select:
                try:
                    x = x[self._raw_select_index]
                    p = p[self._raw_select_index]
                except: warnings.warn('Particle select condition is not valid! All particles are used.')
            # Centering
            x = x - np.average(x)
            p = p - np.average(x)
            # The geometric emittance calculation refer to
            # http://nicadd.niu.edu/~syphers/uspas/2018w/some-notes-on-ellipses.html
            # unnormalize x to mm
            x=x*one_over_k0_mm
            # calculate x prime in mrad
            y=p/p1_array*1.e3
            s11=np.sum(np.multiply(np.square(x), weights))/sum_weight
            s22=np.sum(np.multiply(np.square(y), weights))/sum_weight
            s12=np.sum(np.multiply(np.multiply(x, y), weights))/sum_weight
            eps_pi=np.sqrt(s11*s22-np.square(s12))
            # alpha does not have unit
            alphas[i] = -s12/eps_pi
            # beta in unit of meter
            betas[i]  =  s11/eps_pi
            # gamma in unit of 1/meter
            gammas[i] =  s22/eps_pi
            norm_emittances[i] = eps_pi*mean_p1
        return norm_emittances, [alphas, betas, gammas]

################################method save_tag_file################################
# Save tage file for particle tracking
# in a tag file of OSIRIS, the ! symble is for comment;
# the first uncommented line is the number of tracked particles;
# the following lines are index of the tracked particles.
# Each particle has two index: the first is the process index (where it is born)
# the second is the particle index in this process.
# User has to make sure that the raw data has been read before calling this.
    def save_tag_file(self, tag_file_name="particles.tags", path=None, if_select = False):
        if path is None:
            path = self.path
        if if_select:
            try:
                if 0==len(self._raw_select_index):
                    print("Warning: no particle is selected! Return without saving a tag file.")
                    return 1
                tag_array = self._raw_tag[self._raw_select_index]
            except: warnings.warn('Particle select condition is not valid! All particles are used.')
        else: tag_array = self._raw_tag
        with open("{}/{}".format(path, tag_file_name),'w',encoding = 'utf-8') as h_file:
        #h_file.close() is guaranteed
            n_particles = tag_array.shape[0]
            h_file.write("{}\n".format(n_particles))
            for i in range(n_particles):
                h_file.write("{}, {}\n".format(tag_array[i,0], tag_array[i,1]))
        return 0

################################method read_data################################
    def read_data(self):
        self._data = np.zeros(self.fileid[self._data_name_in_file].shape, dtype=float_type)
        self.fileid[self._data_name_in_file].read_direct(self._data)
        self._axis_slices = [slice(self._axis_range[0, i], self._axis_range[1, i], self._cell_size[i]) for i in range(self._num_dimensions)]
        if self.field_name in {'p1x1', 'p2x2', 'p3,x3'}:
            self._axis_labels = [self._axis_labels_original[0], self._axis_labels_original[3]]
            self._axis_units = [self._axis_units_original[0], self._axis_units_original[3]]
        else:
            pass
            # self._axis_labels and self._axis_units already set in open()
            #self._axis_labels = [self._axis_labels_original[i] for i in range(self._num_dimensions)]
            #self._axis_units = [self._axis_units_original[i] for i in range(self._num_dimensions)]
        self._fig_title = 't = {0:.2f}'.format(self.time)

################################method dirs_in_AXIS_data################################
    def dirs_in_AXIS_data(self, dir):
        '''
            Getting dir_in_AXIS and dir_in_data according to dir.
            OSIRIS 3D h5 file has data axes [y,x,z] and AXIS group [z,x,y].
            OSIRIS 2D h5 file has data axes [x,z] and AXIS group [z,x].
            QuickPIC 3D h5 file has data axes [z,y,x] and AXIS group [x,y,z].
            QuickPIC 2D h5 file has: x-z slices data axes [z,x] and AXIS group [x,z];
                                     y-z slices data axes [z,y] and AXIS group [y,z].
            HiPACE h5 file has data axes [z,x,y] and AXIS group [z,x,y].
        '''
        if self.code_name in {'osiris'}:
            # OSIRIS h5 file has data axes [y,x,z] and AXIS group [z,x,y].
            dir_in_AXIS = dir
            dir_in_data = self.num_dimensions-1-dir
        elif self.code_name in {'quickpic'}:
            # QuickPIC h5 file has data axes [z,y,x] and AXIS group [x,y,z].
            if self.num_dimensions ==3:
                dir_in_AXIS = (dir+2)%3
                qp_data_axes = np.array((0,2,1))
                dir_in_data = qp_data_axes[dir]
            else: # 2 dimensional data
                dir_in_AXIS = 1-dir
                dir_in_data = dir
        elif self.code_name in {'hipace'}:
            # HiPACE h5 file has data axes [z,x,y] and AXIS group [z,x,y].
            dir_in_AXIS = dir
            dir_in_data = dir
        # Deprecated
        # return self._h5_AXIS_seq[dir], self._h5_data_seq[dir]
        return dir_in_AXIS, dir_in_data

################################method read_data_slice################################
    def read_data_slice(self, dir = 2, pos = None):
        '''dir is the direction perpendicular to the slice plane.
           dir = 0 corresponds to longidudinal direction (z).
           dir = 1 corresponds to the first transverse direction (x).
           dir = 2 corresponds to the second transverse direction (y).
           OSIRIS h5 file has data axes [y,x,z] and AXIS group [z,x,y].
           QuickPIC h5 file has data axes [z,y,x] and AXIS group [x,y,z].
           HiPACE h5 file has data axes [z,x,y] and AXIS group [z,x,y].
        '''
        if 3 > self._num_dimensions:
            raise RuntimeError('Method OutFile.read_data_slice() cannot work on data with dimension number < 3!')
        if dir not in range (3):
            raise ValueError('Slice direction should be 0, 1 or 2!')

        # Setting dir_in_AXIS and dir_in_data
        dir_in_AXIS, dir_in_data = self.dirs_in_AXIS_data(dir)

        # Determine the slice position
        slice_dir_len = self.fileid[self._data_name_in_file].shape[dir_in_data]
        if pos is None:
            #set pos at the middle of the box
            pos = (self._axis_range[0, dir_in_AXIS] + self._axis_range[1, dir_in_AXIS])/2.
            pos_index = int(slice_dir_len/2.)
        else:
            #get the index of the nearest grid. +0.5 here has a similar effect of rounding.
            pos_index = int((pos - self._axis_range[0, dir_in_AXIS]) / self._cell_size[dir_in_AXIS] + 0.5)
            if slice_dir_len<=pos_index:
                print('Warning: pos is larger than the upper bundary! Force slicing at the upper bundary.')
                pos_index = slice_dir_len-1
                pos = self._axis_range[1, dir_in_AXIS]
            elif 0>pos_index:
                print('Warning: pos is smaller than the lower bundary! Force slicing at the lower bundary.')
                pos_index = 0
                pos = self._axis_range[0, dir_in_AXIS]
        # Read date matrix from h5 file
        new_shape = [self.fileid[self._data_name_in_file].shape[i] for i in range(3) if i != dir_in_data]
        slice_tuple = tuple([slice(None, None, None) if i != dir_in_data else pos_index for i in range(3)])
        self._data = np.zeros(new_shape, dtype=float_type)
        self.fileid[self._data_name_in_file].read_direct(self._data, source_sel=slice_tuple)
        if self.code_name=='hipace' or (self.code_name=='quickpic' and dir!=0):
            # Transpose if because the data axes order is reversed compared to OSIRIS, for all dir cases of HiPACE and dir==1, 2 cases of QuickPIC
            self._data = np.transpose(self._data)

        # Determine axes perpendicular to the dir
        self._axis_slices = [slice(self._axis_range[0, i], self._axis_range[1, i], self._cell_size[i]) for i in range(3) if i!=dir_in_AXIS]
        if self.code_name=='quickpic' and dir!=0:
            # Flip for QuickPIC with dir==1 or 2 because it has different order in AXIS
            self._axis_slices = np.flip(self._axis_slices)
        # self._axis_labels_original and self._axis_units_original are defined in thie class, independent of the code.
        self._axis_labels = [self._axis_labels_original[i] for i in range(self._num_dimensions) if i!=dir]
        self._axis_units = [self._axis_units_original[i] for i in range(self._num_dimensions) if i!=dir]
        self._fig_title = 't = {0:.2f}, slice at {1} = {2}'.format(self.time, self._axis_labels_original[dir], pos)

################################method read_data_project################################
    def read_data_project(self, *args, **kwargs):
        '''Read data, automatically detect the dimension of data and choose the proper project method'''
        self.read_data()
        project_methods = (self.data_project2d, self.data_project3d)
        return project_methods[self._data.ndim-2](*args, **kwargs)

################################method data_project################################
    def data_project(self, *args, **kwargs):
        '''Automatically detect the dimension of data and choose the proper project method'''
        project_methods = (self.data_project2d, self.data_project3d)
        return project_methods[self._data.ndim-2](*args, **kwargs)

################################method data_project3d################################
    def data_project3d(self, dir = 2, if_abs = False, if_square = False):
        '''
            Taking projection and devide by number of grid.
            dir is the direction the summation is taking through.
            dir = 0 corresponds to longidudinal direction (z).
            dir = 1 corresponds to the first transverse direction (x).
            dir = 2 corresponds to the second transverse direction (y).
            OSIRIS h5 file has data axes [y,x,z] and AXIS group [z,x,y].
            QuickPIC h5 file has data axes [z,y,x] and AXIS group [x,y,z].
            HiPACE h5 file has data axes [z,x,y] and AXIS group [z,x,y].
        '''
        if dir not in range (3):
            raise ValueError('Project direction should be 0, 1 or 2!')

        self._axis_labels = [self._axis_labels[i] for i in range(len(self._axis_labels)) if i!=dir]
        self._axis_units = [self._axis_units[i] for i in range(len(self._axis_units)) if i!=dir]
        self._fig_title = 't = {0:.2f}, {1}projection'.format(self.time, 'absolute value ' if if_abs else '')

        # Setting dir_in_AXIS and dir_in_data
        dir_in_AXIS, dir_in_data = self.dirs_in_AXIS_data(dir)

        if if_abs:
            self._data = np.absolute(self._data)
        if if_square:
            self._data = np.square(self._data)
        self._data = np.sum(self._data, axis = dir_in_data)/self.fileid[self._data_name_in_file].shape[dir_in_data]
        if self.code_name in {'quickpic', 'hipace'}:
            # Transpose if because the data axes order is reversed compared to OSIRIS
            self._data = np.transpose(self._data)

        self._axis_slices = [self._axis_slices[i] for i in range(3) if i!=dir_in_AXIS]
        if self.code_name in {'quickpic'}:
            # Flip for QuickPIC because it has different order in AXIS
            self._axis_slices = np.flip(self._axis_slices)

################################method data_project2d################################
    def data_project2d(self, dir = 1, if_abs = False, if_square = False):
        '''dir is the direction the summation is taking through.
           Currently the project dir is correct for OSIRIS.
           For other codes this shall be also correct if the data is taken by self.read_data_slice(), or self.read_data() followed by self.data_project3d()
        '''
        if dir not in range (2):
            raise ValueError('Project direction should be 0 or 1!')
        self._axis_slices = [self._axis_slices[1-dir]]
        self._axis_labels = [self._axis_labels_original[1-dir], None]
        self._axis_units = [self._axis_units_original[1-dir], None]
        self._fig_title = 't = {0:.2f}, {1}{2}projection'.format(self.time, 'absolute value ' if if_abs else '', 'squared ' if if_square else '')
        if if_abs:
            self._data = np.absolute(self._data)
        if if_square:
            self._data = np.square(self._data)
        #self._data = np.sum(self._data, axis = 1-dir)/self.fileid[self._data_name_in_file].shape[1-dir]
        self._data = np.sum(self._data, axis = 1-dir)/self._data.shape[1-dir] # Need debugging

################################method data_center_of_mass2d################################
    def data_center_of_mass2d(self, dir = 0, if_abs = True, if_square = False, weigh_threshold=0.):
        '''
            get a line of the center of mass of a 2D data along "dir" direction.
            self._data will be transformed to a 1D numpy array.
        '''
        if 2 != self._data.ndim: raise RuntimeError('Data not 2 dimensional! Exit.')
        if dir not in range (2):
            raise ValueError('Direction should be 0 or 1!')
        if 1==dir: weights = self._data
        else: weights = np.transpose(self._data)
        if if_square:
            weights = np.square(weights)
        if if_abs:
            weights = np.absolute(weights)
        transverse_mgrid = np.mgrid[self._axis_slices[1-dir]]
        # Find the bundary of the data region
        weights_sum = np.sum(weights, axis=1)
        weights_sum_max = weights_sum.max()
        if weights_sum_max>weigh_threshold:
            for i_bund_left in range(weights_sum.shape[0]):
                if weights_sum[i_bund_left]>weigh_threshold: break
            for i_bund_right in range(weights_sum.shape[0], 1, -1):
                if weights_sum[i_bund_right-1]>weigh_threshold: break
        else: raise RuntimeError("Weights maximum = {} < threshold. Cannot proceed.".format(weights_sum_max))
        center_of_mass_len = i_bund_right-i_bund_left
        center_of_mass = np.zeros(center_of_mass_len)
        for i in range(center_of_mass_len):
            # The weight may still have small values in the range between i_bund_left and i_bund_right
            if weights_sum[i+i_bund_left]>0: center_of_mass[i] = np.average(transverse_mgrid, weights=weights[i+i_bund_left])
        self._data = center_of_mass
        #print(self._data)
        whole_slice = self._axis_slices[dir]
        self._axis_slices = [np.s_[(whole_slice.start+whole_slice.step*i_bund_left):(whole_slice.start+whole_slice.step*i_bund_right):whole_slice.step]]
        #self._axis_labels = [self._axis_labels_original[dir], self._axis_labels_original[1-dir]]
        #self._axis_units = [self._axis_units_original[dir], self._axis_units_original[1-dir]]
        if 1==dir:
            self._axis_labels = self._axis_labels[::-1]
            self._axis_units = self._axis_units[::-1]
        self._fig_title = 't = {0:.2f}, centre of mass along {1} direction'.format(self.time, self._axis_labels_original[dir])
        return i_bund_left, center_of_mass_len

################################method data2D_slice_spread################################
    def data2D_slice_spread(self, dir = 0, method = 'rms', weigh_threshold=0., relative_spread= False):
        '''
            Get a slice spread of a 2D data along "dir" direction.
            self._data will be transformed to a 1D numpy array.
            method: can be either 'rms' (root-mean-square) or 'lfit' (fit by Lorentz function)
        '''
        if 1==dir: weights = self._data
        else: weights = np.transpose(self._data)
        transverse_mgrid = np.mgrid[self._axis_slices[1-dir]]
        if 'rms' == method:
            # Get the weighted average
            i_bund_left, center_of_mass_len = self.data_center_of_mass2d(dir = dir, if_abs = True, if_square = False, weigh_threshold=weigh_threshold)
            # self._data is now slice average
            #print(weights[i_bund_left:i_bund_right].shape)
            #spread_array = np.average(np.square(transverse_mgrid), axis=0, weights=weights[i_bund_left:i_bund_right]) - np.square(self._data)
            spread_array = np.zeros(center_of_mass_len)
            square_transverse_mgrid = np.square(transverse_mgrid)
            for i in range(center_of_mass_len):
                spread_array[i] = np.average(square_transverse_mgrid, weights=weights[i+i_bund_left]) - self._data[i]**2
        else:
            # Not implemented
            raise NotImplementedError('method {} not implemented yet.'.format(method))
        # If relative_spread, transform to relative spread
        if relative_spread: spread_array/=self._data
        self._data = spread_array
        self._fig_title = 't = {0:.2f}, slice spread'.format(self.time)
        self._axis_labels[1] = '$\\Delta$'+self._axis_labels[1]

################################method read_data_lineout################################
    def read_data_lineout(self, dir = 0, pos = None):
        '''
            Read lineout of the data and save to self._data. Also modify self._axis_slices, self._axis_labels, self._axis_units and self._fig_title.
            dir : Integer.
                  Can be [0, 1] in 2D case, and [0, 1, 2] in 3D case.
            pos : None or a turple of two floats.
                  The position to take lineout.
                  If pos is None, set the position at the center of the simulation box.
            OSIRIS h5 file has data axes [y,x,z] and AXIS group [z,x,y].
            QuickPIC h5 file has data axes [z,y,x] and AXIS group [x,y,z].
            HiPACE h5 file has data axes [z,x,y] and AXIS group [z,x,y].
        '''
        if 2 > self._num_dimensions:
            raise RuntimeError('Method OutFile.read_data_lineout() cannot work on data with dimension number < 2!')
        if dir not in range(self._num_dimensions):
            raise ValueError('Lineout direction should be in range({0})!'.format(self._num_dimensions))

        # Setting dir_in_AXIS and dir_in_data
        dir_in_AXIS, dir_in_data = self.dirs_in_AXIS_data(dir)

        # box_info_array contains the information of the lower boundary position, cell size and upper boundary position in the directions except the dir direction.
        box_info_array = np.array([[self._axis_range[0, i], self._cell_size[i], self._axis_range[1, i],] for i in range(self._num_dimensions) if i!=dir_in_AXIS])
        # For QuickPIC, if dir==1 or 2, the AXIS should be flipped so that the AXIS is in ascending order
        if self.code_name=='quickpic' and dir!=0: box_info_array = np.flipud(box_info_array)
        if pos is None:
            # Set pos to the center of box
            pos = (box_info_array[:,0]+box_info_array[:,2])/2
        # Calculate the index of the nearest grid of the lineout position. +0.5 here has a similar effect of rounding.
        pos_index = ((np.array(pos[:self._num_dimensions-1]) - box_info_array[:, 0]) / box_info_array[:, 1] + 0.5).astype(np.int, copy=False)
        pos_index = pos_index.tolist()
        # Correct pos to the cell position.
        pos = box_info_array[:,0] + box_info_array[:,1]*np.array(pos_index)

        # data size in the directions perpendicular to dir
        size_perp_dir = [self.fileid[self._data_name_in_file].shape[i] for i in range(self._num_dimensions) if i!=dir_in_data]
        # For OSIRIS, size_perp_dir should be flipped to be ascending
        if self.code_name=='osiris': size_perp_dir = size_perp_dir[::-1]
        # For QuickPIC, if dir==0, size_perp_dir should be flipped to be ascending
        elif self.code_name=='quickpic' and 0==dir: size_perp_dir = size_perp_dir[::-1]

        # If pos_index out of the array range, limit it to the boundary.
        for i in range(self._num_dimensions-1):
            if size_perp_dir[i]<=pos_index[i]:
                print('Warning: pos[{0}] is larger than the upper bundary! Force slicing at the upper bundary.'.format(i))
                pos_index[i] = size_perp_dir[i]-1
            elif 0>pos_index[i]:
                print('Warning: pos[{0}] is smaller than the lower bundary! Force slicing at the lower bundary.'.format(i))
                pos_index[i] = 0

        new_shape = self.fileid[self._data_name_in_file].shape[dir_in_data]
        # To insert the dir_in_data to the pos_index for h5 slicing.
        # Make pos_index to have the same order as data matrix in h5 file. For OSIRIS and QuickPICin the dir==0 case, flip the position index first.
        if self.code_name=='osiris' or (self.code_name=='quickpic' and 0==dir): pos_index = pos_index[::-1]
        pos_index.insert(dir_in_data, slice(None, None, None))
        slice_tuple = tuple(pos_index)

        self._data = np.zeros(new_shape, dtype=float_type)
        self.fileid[self._data_name_in_file].read_direct(self._data, source_sel=slice_tuple)
        self._axis_slices = [slice(self._axis_range[0, dir_in_AXIS], self._axis_range[1, dir_in_AXIS], self._cell_size[dir_in_AXIS]),]

        self._axis_labels = [self._axis_labels_original[dir], self._field_names[self.field_name]]
        self._axis_units = [self._axis_units_original[dir], '']
        tmp_list = [self._axis_labels_original[i] for i in range(self._num_dimensions) if i!=dir]
        tmp_list = ['{0} = {1}'.format(tmp_list[i], pos[i]) for i in range(self._num_dimensions-1)]
        self._fig_title = 't = {0:.2f}, lineout at '.format(self.time) + ', '.join(tmp_list)

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
        self._fig_title = 't = {0:.2f}, spectrogram'.format(self.time)

################################method data_profile1d################################
    def data_profile1d(self, peak_height_min_ratio=0.00001):
        '''Turn 1-D high-frequency signal (eg. laser) to its profile.'''
        if self._data.ndim!=1:
            raise RuntimeError('Data is not one dimensional! The OutFile.profile1d() method cannot proceed.')
        self._data = np.abs(self._data)
        # Longitudinal grid points
        lon_grid = np.mgrid[self._axis_slices[0]]
        lon_len_minus1 = len(lon_grid)-1
        peaks_ind=find_peaks(self._data, height=0.00001*np.average(self._data))[0]
        # Adding points at the beginning and ending of the peak_ind, preparing for interpolation
        if 0<peaks_ind[0]: peaks_ind=np.append(0,peaks_ind)
        if lon_len_minus1>peaks_ind[-1]: peaks_ind=np.append(peaks_ind, lon_len_minus1)
        self._data = np.interp(lon_grid, lon_grid[peaks_ind], self._data[peaks_ind])

################################method data_profile2d################################
    def data_profile2d(self, peak_height_min_ratio=0.00001, dir=0):
        '''Turn 2-D high-frequency signal (eg. laser) to its profile.
           dir=0 means the high frequency signal is along 0th direction.
        '''
        if self._data.ndim!=2:
            raise RuntimeError('Data is not two dimensional! The OutFile.profile2d() method cannot proceed.')
        self._data = np.abs(self._data)
        # Longitudinal grid points
        lon_grid = np.mgrid[self._axis_slices[dir]]
        lon_len_minus1 = len(lon_grid)-1
        abs_avg = np.average(self._data)
        for transverse_ind in range(self._data.shape[dir]):
            peaks_ind=find_peaks(self._data[transverse_ind], height=0.00001*abs_avg)[0]
            if len(peaks_ind)<1: continue
            # Adding points at the beginning and end of the peak_ind, preparing for interpolation
            if 0<peaks_ind[0]: peaks_ind=np.append(0,peaks_ind)
            if lon_len_minus1>peaks_ind[-1]: peaks_ind=np.append(peaks_ind, lon_len_minus1)
            self._data[transverse_ind] = np.interp(lon_grid, lon_grid[peaks_ind], self._data[transverse_ind, peaks_ind])
        # not yet finished. Currently this can only work with dir=0.

################################method plot_data_line################################
    def plot_data_line(self, h_fig=None, h_ax=None, if_flip_xy=False, semilogy=False, if_xlabel=None, if_ylabel=None, if_title=None, multiple=1., offset=0., if_z2zeta=False, **kwargs):
        '''Plot 1D data in a line. h_ax is the handle of the axis to be plotted on. If h_ax=None, a new figure and axis is created.
        if_z2zeta: boolean
                   If True, offset x axis by -t (thus it is effectively zeta=z-t)
        '''
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
        if if_z2zeta:
            x_slice = slice(self._axis_slices[0].start-self.time, self._axis_slices[0].stop-self.time, self._axis_slices[0].step)
        else:
            x_slice = self._axis_slices[0]
        if if_flip_xy:
            plotfunc((self._data*multiple)+offset, np.mgrid[x_slice], **kwargs)
            set_xlabel = h_ax.set_ylabel
            set_ylabel = h_ax.set_xlabel
            get_xlabel = h_ax.get_ylabel
            get_ylabel = h_ax.get_xlabel
        else:
            plotfunc(np.mgrid[x_slice], (self._data*multiple)+offset, **kwargs)
            set_xlabel = h_ax.set_xlabel
            set_ylabel = h_ax.set_ylabel
            get_xlabel = h_ax.get_xlabel
            get_ylabel = h_ax.get_ylabel

        if if_xlabel is None:
            # Default if_xlabel to True if get_xlabel is empty
            if '' == get_xlabel(): if_xlabel = True
        if if_ylabel is None:
            # Default if_ylabel to True if get_ylabel is empty
            if '' == get_ylabel(): if_ylabel = True
        if if_title is None:
            # Default if_title to True if get_title is empty
            if '' == h_ax.get_title(): if_title = True

        if if_xlabel:
            set_xlabel('{0} [{1}]'.format(self._axis_labels[0], self._axis_units[0]))
        if if_ylabel:
            set_ylabel('{0} [{1}]'.format(self._axis_labels[1], self._axis_units[1]))
        if if_title:
            h_ax.set_title(self._fig_title)
        return h_fig, h_ax

################################method pcolor_data_2d################################
    def pcolor_data_2d(self, h_fig=None, h_ax=None, if_colorbar=True, colorbar_orientation='vertical', if_log_colorbar=False, vmin=None, vmax=None, cmap=plt.cm.jet, alpha=None, if_z2zeta=False, if_transpose=None, **kwargs):
        '''Plot 2D data in as pcolor
        if_z2zeta: boolean
                   If True, offset x axis by -t (thus it is effectively zeta=z-t)
        if_transpose: boolean or None
                   If None, transpose the plot for QuickPIC by default, and do not transpose for OSIRIS and HiPACE'''
        if self._data.ndim!=2:
            raise RuntimeError('Data is not two dimensional! The OutFile.pcolor_data_2d() method cannot proceed.')
        if if_transpose is None:
            if self.code_name == 'quickpic' and self.out_type != 'RAW': if_transpose = True
            else: if_transpose = False
        if if_transpose:
            self._data = np.transpose(self._data)
            self._axis_slices = np.flip(self._axis_slices)
            self._axis_labels = np.flip(self._axis_labels)
            self._axis_units = np.flip(self._axis_units)
        if if_z2zeta:
            x_slice = slice(self._axis_slices[0].start-self.time, self._axis_slices[0].stop-self.time, self._axis_slices[0].step)
        else:
            x_slice = self._axis_slices[0]
        y_spread, x_spread = np.mgrid[self._axis_slices[1], x_slice]
        if h_fig is None:
            h_fig = plt.figure()
        if h_ax is None:
            h_ax = h_fig.add_subplot(111)
        if if_log_colorbar:
            h_plot = h_ax.pcolormesh(x_spread, y_spread, np.absolute(self._data), norm=LogNorm(vmin=vmin, vmax=vmax), cmap=cmap, alpha=alpha, **kwargs)
        else:
            h_plot = h_ax.pcolormesh(x_spread, y_spread, self._data, vmin=vmin, vmax=vmax, cmap=cmap, alpha=alpha, **kwargs)
        if self._axis_labels[0] is not None:
            xlabel = self._axis_labels[0]
            if self._axis_units[0] is not None:
                xlabel = xlabel + ' [{}]'.format(self._axis_units[0])
            h_ax.set_xlabel(xlabel)
        if self._axis_labels[1] is not None:
            ylabel = self._axis_labels[1]
            if self._axis_units[1] is not None:
                ylabel = ylabel + ' [{}]'.format(self._axis_units[1])
            h_ax.set_ylabel(ylabel)
        if if_colorbar:
            self._color_bar = plt.colorbar(h_plot, ax=h_ax, orientation=colorbar_orientation)
            try: self._color_bar.set_label(self._field_names[self._data_name_in_file])
            except: self._color_bar.set_label('Counts [arb. units]')
        h_ax.set_title(self._fig_title)
        return h_fig, h_ax

################################method isosurface_data_3d################################
    def isosurface_data_3d(self, h_fig=None, h_ax=None, **kwargs):
        '''Plot 3D data as isosurface.'''
        if self._data.ndim!=3:
            raise RuntimeError('Data is not 3 dimensional! The OutFile.isosurface_data_3d() method cannot proceed.')
        raise NotImplementedError('3D data plot not implemented')

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

################################method raw_hist_gamma################################
    def raw_hist_gamma(self, num_bins=256, range_max=None, range_min=None, if_select = False):
        '''Get histogram from ene raw data +1 = gamma.'''
        weights=np.absolute(self._raw_q)
        gamma = self._raw_ene+1.
        if if_select:
            try:
                weights = weights[self._raw_select_index]
                gamma = gamma[self._raw_select_index]
            except: warnings.warn('Particle select condition is not valid! All particles are used.')
        if range_max is None:
            range_max = gamma.max()
        if range_min is None:
            range_min = gamma.min()
        hist, bin_edges = np.histogram(gamma, num_bins, (range_min, range_max), weights=weights)
        bin_edges = bin_edges[0:-1]
        return bin_edges, hist

################################method raw_hist_current################################
    def raw_hist_current(self, n0_per_cc, dir=1, num_bins=256, range_max=None, range_min=None, if_select = False, if_abs=True):
        '''Get histogram for current distribution.
        dir : integer
              The longitudinal direction. 1 for z, 2 for x, 3 for y.
        n0_per_cc : float
              the normalization density in cm^-3'''
        if if_abs: q=np.absolute(self._raw_q)
        else: q=self._raw_q
        if 1==dir: z=self._raw_x1
        elif 2==dir: z=self._raw_x2
        elif 3==dir: z=self._raw_x3
        else: raise KeyError('dir should be 1, 2 or 3.')
        if if_select:
            try:
                q = q[self._raw_select_index]
                z = z[self._raw_select_index]
            except: warnings.warn('Particle select condition is not valid! All particles are used.')
        if range_max is None:
            range_max = z.max()
        if range_min is None:
            range_min = z.min()
        hist, bin_edges = np.histogram(z, num_bins, (range_min, range_max), weights=q, density=False)
        bin_edges = bin_edges[0:-1]

        # Transform charge to Coulomb
        cell_volume_norm = 1.0
        for i in range(self.num_dimensions):
            cell_volume_norm = cell_volume_norm*self._cell_size[i]
        if 3>self.num_dimensions:
            print('Warning! Similation is in {} dimensional. Charge calculation may not be correct.'.format(self.num_dimensions))
        hist *= (cell_volume_norm/np.sqrt(n0_per_cc)*2.404351211612064e-2)

        # Transform to Amper
        one_over_k0_m_square = 2.8239587227915743e7/n0_per_cc
        hist /= (bin_edges[1]-bin_edges[0])*one_over_k0_m_square**0.5/2.99792458e8

        return bin_edges, hist

################################method plot_raw_hist_gamma################################
    def plot_raw_hist_gamma(self, h_fig=None, h_ax=None, num_bins=256, range_max=None, range_min=None, if_select = False, **kwargs):
        '''Plot histogram from ene raw data +1 = gamma.'''
        if h_fig is None:
            h_fig = plt.figure()
        if h_ax is None:
            h_ax = h_fig.add_subplot(111)
        bin_edges, hist = self.raw_hist_gamma(num_bins, range_max, range_min, if_select = if_select)
        h_ax.plot(bin_edges, hist, **kwargs)
        h_ax.set_xlabel('$\\gamma$')
        h_ax.set_ylabel('Counts [arbt. units]')
        h_ax.set_title('$t = {0:.2f}$'.format(self.time))
        return h_fig, h_ax, bin_edges, hist

################################method plot_tracks################################
    def plot_tracks(self, h_fig=None, h_ax=None, x_quant = b'x1', y_quant = b'x2', mwin_v = 0., mwin_t_offset = 0.):
        '''Plot tracks from tracking data'''
        if h_fig is None:
            h_fig = plt.figure()
        if h_ax is None:
            h_ax = h_fig.add_subplot(111)
        n_tracks = self.fileid.attrs.get("NTRACKS")[0]
        # in OSIRIS track files, QUANTS attribute is a list of string telling which quantities are recorded in data. The 0th element in QUANTS should be ignored.
        # the quants in OSIRIS output are bytes, thus have prefix b.
        quants = self.fileid.attrs.get("QUANTS")
        t_quant_ind = np.where(quants == b't')[0][0]-1
        x2_quant_ind = np.where(quants == b'x2')[0][0]-1
        x3_quant_ind = np.where(quants == b'x3')[0][0]-1
        x_quant_ind = np.where(quants == x_quant)[0][0]-1
        y_quant_ind = np.where(quants == y_quant)[0][0]-1
        itermap = np.zeros(self.fileid['itermap'].shape, dtype=int)
        self.fileid['itermap'].read_direct(itermap)
        data = np.zeros(self.fileid['data'].shape, dtype=float_type)
        self.fileid['data'].read_direct(data)
        seperation_pointer = 0
        print(np.where(1 == itermap[:,0])[0])
        print(len(itermap))
        #for i in range(n_tracks):
        for i in range(100):
            #h_ax.plot(data[seperation_pointer:itermap[i,1]+seperation_pointer, x_quant_ind], data[seperation_pointer:itermap[i,1]+seperation_pointer, y_quant_ind])
            # select condition
            if (data[itermap[i,1]+seperation_pointer-1, x2_quant_ind]**2+data[itermap[i,1]+seperation_pointer-1, x3_quant_ind]**2)<0.0289:
                # replace x1 by x1-mwin_v*t+mwin_t_offset
                h_ax.plot(data[seperation_pointer:itermap[i,1]+seperation_pointer, x_quant_ind]-data[seperation_pointer:itermap[i,1]+seperation_pointer, t_quant_ind]*mwin_v+mwin_t_offset, data[seperation_pointer:itermap[i,1]+seperation_pointer, y_quant_ind])
            seperation_pointer += itermap[i,1]
        #h_ax.set_xlabel('$\\gamma$')
        #h_ax.set_ylabel('Counts [arbt. units]')
        #h_ax.set_title('$t = {0:.2f}$'.format(self.time))
        return h_fig, h_ax

################################method raw_hist2D################################
    def raw_hist2D(self, dims='p1x1', num_bins=128, range=None, if_select = False):
        '''Generate 2D histogram for phasespace from raw data.
           One does not have to read raw data before calling this.
           dims can be combinations of 'x1', 'x2', 'x3', 'p1', 'p2', 'p3'.'''
        weights=np.absolute(self.read_raw_q())
        if weights.size<2:
            warnings.warn("No particle contained in the RAW file! Skipping...")
            return
        y_type = dims[0]
        y_dir = int(dims[1])-1
        x_type = dims[2]
        x_dir = int(dims[3])-1
        type_tuple = ('x','p')
        if (y_type not in type_tuple) or (x_type not in type_tuple):
            raise NotImplementedError('dims {} not implemented!'.format(dims))
        if (y_dir > 2) or (x_dir > 2):
            raise NotImplementedError('direction in dims should be in (1, 2, 3)!')
        read_raw_tuple = ((self.read_raw_x1(), self.read_raw_x2(), self.read_raw_x3()), (self.read_raw_p1(), self.read_raw_p2(), self.read_raw_p3()))
        #raw_tuple = ((self._raw_x1, self._raw_x2, self._raw_x3), (self._raw_p1, self._raw_p2, self._raw_p3))
        label_tuple = (('$k_p z$', '$k_p x$', '$k_p y$'), ('$p_z / m_ec$', '$p_x / m_ec$', '$p_y / m_ec$'))
        y_type_ind = type_tuple.index(y_type)
        y_array = read_raw_tuple[y_type_ind][y_dir]
        x_type_ind = type_tuple.index(x_type)
        x_array = read_raw_tuple[x_type_ind][x_dir]
        self._axis_labels = [label_tuple[x_type_ind][x_dir], label_tuple[y_type_ind][y_dir]]
        if if_select:
            try:
                weights = weights[self._raw_select_index]
                x_array = x_array[self._raw_select_index]
                y_array = y_array[self._raw_select_index]
            except: warnings.warn('Particle select condition is not valid! All particles are used.')
        self._data, yedges, xedges = np.histogram2d(y_array, x_array, bins=num_bins, range=range, weights=weights)
        #self._data=np.transpose(self._data)
        self._axis_slices = [slice(xedges[0], xedges[-1], xedges[1]-xedges[0]), slice(yedges[0], yedges[-1], yedges[1]-yedges[0])]
        self._axis_units = [None,None]
        self._fig_title = 't = {0:.2f}, phasespace'.format(self.time)
        return

################################method plot_raw_hist2D################################
    def plot_raw_hist2D(self, h_fig=None, h_ax=None, dims='p1x1', num_bins=128, range=None, if_select = False, **kwargs):
        '''Plot 2D histogram for phasespace from raw data.
           Please make sure the corresponding raw data is read before calling this.
           dims can be combinations of 'x1', 'x2', 'x3', 'p1', 'p2', 'p3'.'''
        self.raw_hist2D(dims, num_bins, range, if_select)
        h_fig, h_ax = self.plot_data(h_fig, h_ax, **kwargs)
        return h_fig, h_ax

################################method raw_mean_rms_ene################################
    def raw_mean_rms_ene(self, if_select = False):
        '''Return weighted mean value and RMS spread of ene.
           If if_select and self._raw_select_index is valid, use selection of macroparticles.
           Otherwise use all the macroparticles.'''
        weights=np.absolute(self._raw_q)
        ene = self._raw_ene
        if if_select:
            try:
                weights = weights[self._raw_select_index]
                ene = ene[self._raw_select_index]
            except: warnings.warn('Particle select condition is not valid! All particles are used.')
        ene_avg, sum_weights = np.average(ene, weights=weights, returned=True)
        ene_rms_spread = np.sqrt(np.sum(np.square(ene-ene_avg)*weights)/sum_weights)
        return ene_avg, ene_rms_spread

################################method plot_data################################
    def plot_data(self, *args, **kwargs):
        '''Automatically detect the dimension of data and choose the proper plot method'''
        plot_methods = (self.plot_data_line, self.pcolor_data_2d, self.isosurface_data_3d)
        h_fig, h_ax = plot_methods[self._data.ndim-1](*args, **kwargs)
        h_ax.minorticks_on()
        return h_fig, h_ax

################################method fit_for_W################################
    def fit_for_W(self, *args, **kwargs):
        '''Automatically detect the dimension of data and choose the proper fit dimension for W'''
        fit_methods = (self.fit_for_W_1d, self.fit_for_W_2d)
        return fit_methods[self._data.ndim-1](*args, **kwargs)

################################method fit_for_W_2d################################
    def fit_for_W_2d(self, h_fig=None, h_ax=None, guess_values=None):
        '''Do fitting for W. If 'abs' == data_preprocess, use read_data_project(if_abs=True) to read data absolut value projection along 0th direction before this method. If 'square' == data_preprocess, use read_data_project(if_square=True) instead, and fited W should be .'''
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
            raise RuntimeError('self._data[start] and self._data[stop] have the same sign! The OutFile.zero_poin() method cannot take a zero point between them.')
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
        file1 = OutFile(code_name='osiris',path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/os_beam3D/os_beam3D262',field_name='raw',spec_name='He_e',use_num_list=True,out_num=50)
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
        h_fig, h_ax = file1.plot_raw_hist2D(dims='p2x2', if_select=True, if_log_colorbar=True)
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
        file1 = OutFile(code_name='osiris',path='',field_name='raw',out_num=0)
        file1.open(filename='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/hi_qp_compare/hi3/trailer.h5')
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
        h_fig, h_ax = file1.plot_raw_hist2D(dims='p1x1', if_select=True, if_log_colorbar=True, if_colorbar=True,range=[[19569.4, 19569.6], [5.,11.]])
        file1.close()
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
        file1 = OutFile(code_name='quickpic',path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/qp_hi_compare/qp1',field_name='raw',spec_name='Beam0003', out_num=6900)
        file1.open(cell_size_qp_raw = np.array([0.046875, 0.046875, 0.025390625]))
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
        h_fig, h_ax = file1.plot_raw_hist2D(dims='p1x1', if_select=True, if_log_colorbar=True, if_colorbar=True)
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

    test_os()
    test_hi()
    test_qp()
    plt.show()
    #beam_measure()
    #beam_measure_hi()
    #hi_beam3D()
    #beam_measure_beam3D()
    #tracking()
    #select_tag()
    #tmp_hipace3D()
    #scatter()
    #tmp_3D()
    #X1plot()
    #oblique_laser_3D()
    #test_profile()

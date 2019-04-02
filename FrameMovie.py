import outfile
import my_cmap
import os
import matplotlib.pyplot as plt
import numpy as np

class Frames:
    def __init__(self, simulation_path, frame_folder, dirver_type = 0, start_num = 0, stride_num=1, count_num=1, dirver_spec_name='driver', background_spec_name='e', trail_spec_name=None, if_e1=False, if_psi=False, dir=2, save_type='png'):
        self.simulation_path = simulation_path
        self.frame_path = simulation_path+'/'+frame_folder
        self.start_num = start_num
        self.stride_num = stride_num
        self.count_num = count_num
        self.dirver_type = dirver_type
        self.dirver_spec_name = dirver_spec_name
        self.background_spec_name = background_spec_name
        self.trail_spec_name = trail_spec_name
        self.if_e1 = if_e1
        self.if_psi = if_psi
        self.dir = dir
        self.save_type = save_type

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

################################property dirver_type################################
    def get_dirver_type(self):
        return self._dirver_type

    def set_dirver_type(self, value):
        self._dirver_type = value

    dirver_type = property(get_dirver_type, set_dirver_type)

################################method plot_beam_driven################################
#plot one frame for beam driven cases
    def plot_beam_driven(self, out_num):
        h_fig = plt.figure(figsize=(8,4.5))
        file1 = outfile.OutFile(path=self.simulation_path, field_name='charge', average='', spec_name=self.background_spec_name, out_num=out_num)
        h_ax = h_fig.add_subplot(111)
        h_ax.set_aspect('equal','box')
        file1.open()
        file1.read_data_slice(dir=self.dir)
        file1.plot_data(h_fig, h_ax, vmin=-4., vmax=0, cmap='gray')
        file1._color_bar.set_label('$\\rho_e$')
        file1.close()

        file1.spec_name=self.dirver_spec_name
        file1.open()
        file1.read_data_slice(dir=self.dir)
        file1.plot_data(h_fig, h_ax, vmin=-5., vmax=0, cmap=my_cmap.cmap_higher_range_transparent(plt.cm.hot))
        file1._color_bar.set_label('$\\rho_d$')
        file1.close()

        if self.trail_spec_name is not None:
            file1.spec_name=self.trail_spec_name
            file1.open()
            file1.read_data_slice(dir=self.dir)
            file1.plot_data(h_fig, h_ax, vmin=-0.2, vmax=0, cmap=my_cmap.cmap_higher_range_transparent())
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
        file1 = outfile.OutFile(path=self.simulation_path, field_name='charge', average='-savg', spec_name='plasma', out_num=out_num)
        h_ax = h_fig.add_subplot(111)
        h_ax.set_aspect('equal', 'box')
        plt.ylim(-10,10)
        file1.open()
        file1.read_data_slice(dir=self.dir)
        file1.plot_data(h_fig, h_ax, cmap='gray', vmin=-5.)
        file1._color_bar.set_label('$\\rho_e$')
        file1.close()

        file1.field_name='e3'
        file1.open()
        file1.read_data_slice(dir=self.dir)
        file1.plot_data(h_fig, h_ax, cmap=my_cmap.cmap_middle_range_transparent())
        file1._color_bar.set_label('$E_y$')
        file1.close()
        plt.tight_layout()
        return h_fig

################################method plot_save################################
#plot one frame and save using the method either plot_save_beam_driven or plot_save_laser_driven, depends on the "dirver_type" property (0 for beam driver and 1 for laser driver)
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
            methods = (self.plot_beam_driven, self.plot_laser_driven)
            h_fig = methods[self.dirver_type](*args, **kwargs)
            plt.savefig(save_file_name, format=self.save_type)
            plt.close(h_fig)

################################method save_frames################################
#save all frames
    def save_frames(self):
        print('Working on simulation \'{0}\' and saving frames at \'{1}\'.'.format(self.simulation_path, self.frame_path))
        try:
            for i in range(self.start_num, self.start_num+self.count_num*self.stride_num, self.stride_num):
                self.plot_save(i)
        except IOError:
            print('Iteration stops at frame number {0}.'.format(i))

################################method make_movie################################
#make movie based on the saved frames
#    def make_movie(self):
#        print('Making movie.mpg under \'{0}\' - this may take a while.'.format(self.frame_path))
#        subprocess.call('mencoder \'{0}/*.png\' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o {0}/movie.mpg'.format(self.frame_path), shell=True)

if __name__ == '__main__':
    #frame1 = Frames('/home/zming/simulations/os2D/os_PT3D4','/home/zming/simulations/os2D/os_PT3D4/movie1', 1, start_num = 63, count_num=1000)
    frame1 = Frames('/home/zming/simulations/os2D/os_beam3D52','movie2', 0, start_num = 0, count_num=999, trail_spec_name='He_e', if_e1=True, if_psi=True, dir=2)
    frame1.save_frames()
    #frame1.make_movie()

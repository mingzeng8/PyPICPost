import outfile
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

class Frames:
    def __init__(self, simulation_path, frame_path, dirver_type = 0, start_num = 0, stride_num=1, count_num=1):
        self.simulation_path = simulation_path
        self.frame_path = frame_path
        self.start_num = start_num
        self.stride_num = stride_num
        self.count_num = count_num
        self.dirver_type = dirver_type

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
                    raise IOError('\'{0}\' is not empty! Use a empty folder instead.'.format(tmp_path))
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

################################method plot_save_beam_driven################################
#plot one frame and save for beam driven cases
    def plot_save_beam_driven(self, out_num):
        h_fig = plt.figure(figsize=(8,4.5))
        file1 = outfile.OutFile(path=self.simulation_path, field_name='charge', average='', spec_name='e', out_num=out_num)
        h_ax = h_fig.add_subplot(111)
        file1.open()
        file1.read_data_project(dir=2)
        file1.plot_data(h_fig, h_ax, cmap='gray')
        file1._color_bar.set_label('$\\rho_e$')
        file1.close()

        file1.spec_name='driver'
        file1.open()
        file1.read_data_project(dir=2)
        cmap=plt.cm.jet
        my_cmap = cmap(np.arange(cmap.N))
        transparency_start=int(cmap.N*0.8)
        my_cmap[transparency_start:cmap.N,-1] = np.linspace(1, 0, cmap.N-transparency_start)
        my_cmap = ListedColormap(my_cmap)
        file1.plot_data(h_fig, h_ax, cmap=my_cmap)
        file1._color_bar.set_label('$\\rho_d$')
        file1.close()

        plt.savefig('{0}/{1}.png'.format(self.frame_path, out_num), format='png')
        plt.close(h_fig)

################################method plot_save_laser_driven################################
#plot one frame and save for laser driven cases
    def plot_save_laser_driven(self, out_num):
        h_fig = plt.figure(figsize=(6.5,5))
        file1 = outfile.OutFile(path=self.simulation_path, field_name='charge', average='-savg', spec_name='plasma', out_num=out_num)
        h_ax = h_fig.add_subplot(111)
        h_ax.set_aspect('equal', 'box')
        plt.ylim(-10,10)
        file1.open()
        file1.read_data_slice(dir=2)
        file1.plot_data(h_fig, h_ax, cmap='gray', vmin=-5.)
        file1._color_bar.set_label('$\\rho_e$')
        file1.close()

        file1.field_name='e3'
        file1.open()
        file1.read_data_slice(dir=2)
        cmap=plt.cm.bwr
        my_cmap = cmap(np.arange(cmap.N))
        transparency_x=[int(cmap.N*0.3), int(cmap.N*0.45), int(cmap.N*0.55), int(cmap.N*0.7)]
        my_cmap[transparency_x[0]:transparency_x[1],-1] = np.linspace(1, 0, transparency_x[1]-transparency_x[0])
        my_cmap[transparency_x[1]:transparency_x[2],-1] = np.zeros(transparency_x[2]-transparency_x[1])
        my_cmap[transparency_x[2]:transparency_x[3],-1] = np.linspace(0, 1, transparency_x[3]-transparency_x[2])
        my_cmap = ListedColormap(my_cmap)
        file1.plot_data(h_fig, h_ax, cmap=my_cmap)
        file1._color_bar.set_label('$E_y$')
        file1.close()
        plt.tight_layout()
        plt.savefig('{0}/{1}.png'.format(self.frame_path, out_num), format='png')
        plt.close(h_fig)

################################method plot_save################################
#plot one frame and save using the method either plot_save_beam_driven or plot_save_laser_driven, depends on the "dirver_type" property (0 for beam driver and 1 for laser driver)
    def plot_save(self, *args, **kwargs):
        methods = (self.plot_save_beam_driven, self.plot_save_laser_driven)
        return methods[self.dirver_type](*args, **kwargs)

################################method save_frames################################
#save all frames
    def save_frames(self):
        print('Working on simulation \'{0}\' and saving frames at \'{1}\'.'.format(self.simulation_path, self.frame_path))
        try:
            for i in range(self.start_num, self.start_num+self.count_num*self.stride_num, self.stride_num):
                print(i)
                self.plot_save(i)
        except IOError:
            print('Iteration stops at frame number {0}.'.format(i))

################################method make_movie################################
#make movie based on the saved frames
#    def make_movie(self):
#        print('Making movie.mpg under \'{0}\' - this may take a while.'.format(self.frame_path))
#        subprocess.call('mencoder \'{0}/*.png\' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o {0}/movie.mpg'.format(self.frame_path), shell=True)

if __name__ == '__main__':
    frame1 = Frames('/Path/to/OSIRIS/running','/Path/to/frames/saving', 1, start_num = 0, count_num=1000)
    frame1.save_frames()
    #frame1.make_movie()

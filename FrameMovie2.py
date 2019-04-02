'''Similar to FrameMovie.py but output two simulations in one frame'''
from FrameMovie import *

class Frames2(Frames):
    def __init__(self, simulation_path, simulation_path2, frame_path, dirver_type = 0, start_num = 0, stride_num=1, count_num=1):
        Frames.__init__(self, simulation_path, frame_path, dirver_type, start_num, stride_num, count_num)
        self.simulation_path2 = simulation_path2

################################property simulation_path2################################
    def get_simulation_path2(self):
        return self._simulation_path2

    def set_simulation_path2(self, value):
        if os.path.isdir(value):
            self._simulation_path2 = value
        else:
            raise IOError('\'{0}\' is not a folder or does not exist!'.format(value))

    simulation_path2 = property(get_simulation_path2, set_simulation_path2)

################################method plot_save_beam_driven2################################
#plot one frame with 2 simulations and save for beam driven cases
    def plot_save_beam_driven2(self, out_num):
        h_fig = plt.figure(figsize=(5,2))
        file1 = outfile.OutFile(path=self.simulation_path, field_name='charge', average='', spec_name='e', out_num=out_num)
        h_ax = h_fig.add_subplot(121)
        h_ax.set_aspect('equal', 'box')
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

        file1 = outfile.OutFile(path=self.simulation_path2, field_name='charge', average='', spec_name='e', out_num=out_num)
        h_ax = h_fig.add_subplot(122)
        h_ax.set_aspect('equal', 'box')
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

        plt.tight_layout()
        plt.savefig('{0}/{1}.png'.format(self.frame_path, out_num), format='png')
        plt.close(h_fig)

################################method plot_save_laser_driven2################################
#plot one frame with 2 simulations and save for laser driven cases
    def plot_save_laser_driven2(self, out_num):
        h_fig = plt.figure(figsize=(5,2))
        file1 = outfile.OutFile(path=self.simulation_path, field_name='charge', average='-savg', spec_name='plasma', out_num=out_num)
        h_ax = h_fig.add_subplot(121)
        h_ax.set_aspect('equal', 'box')
        plt.ylim(-10,10)
        file1.open()
        file1.read_data_slice(dir=2)
        file1.plot_data(h_fig, h_ax, cmap='gray', vmin=-5., if_colorbar=False)
        #file1._color_bar.set_label('$\\rho_e$')
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
        file1._color_bar.set_label('$E_y [k_p m_e c^2/e]$')
        h_ax.set_title('With plasma lens')
        file1.close()

        file1 = outfile.OutFile(path=self.simulation_path2, field_name='charge', average='-savg', spec_name='plasma', out_num=out_num)
        h_ax = h_fig.add_subplot(122)
        h_ax.set_aspect('equal', 'box')
        plt.ylim(-10,10)
        file1.open()
        file1.read_data_slice(dir=2)
        file1.plot_data(h_fig, h_ax, cmap='gray', vmin=-5., if_colorbar=False)
        #file1._color_bar.set_label('$\\rho_e$')
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
        file1._color_bar.set_label('$E_y [k_p m_e c^2/e]$')
        h_ax.set_title('W/o plasma lens')
        file1.close()

        plt.tight_layout()
        plt.savefig('{0}/{1}.png'.format(self.frame_path, out_num), format='png')
        plt.close(h_fig)

################################method plot_save2################################
#plot one frame with 2 simulations and save using the method either plot_save_beam_driven2 or plot_save_laser_driven2, depends on the "dirver_type" property (0 for beam driver and 1 for laser driver)
    def plot_save2(self, *args, **kwargs):
        methods = (self.plot_save_beam_driven2, self.plot_save_laser_driven2)
        return methods[self.dirver_type](*args, **kwargs)

################################method save_frames2################################
#save all frames with 2 simulations
    def save_frames2(self):
        print('Working on simulation \'{0}\' and \'{1}\', saving frames at \'{2}\'.'.format(self.simulation_path, self.simulation_path2, self.frame_path))
        try:
            for i in range(self.start_num, self.start_num+self.count_num*self.stride_num, self.stride_num):
                print(i)
                self.plot_save2(i)
        except IOError:
            print('Iteration stops at frame number {0}.'.format(i))

################################method make_movie################################
#make movie based on the saved frames
#    def make_movie(self):
#        print('Making movie.mpg under \'{0}\' - this may take a while.'.format(self.frame_path))
#        subprocess.call('mencoder \'{0}/*.png\' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o {0}/movie.mpg'.format(self.frame_path), shell=True)

if __name__ == '__main__':
    frame1 = Frames2('/home/zming/simulations/os2D/os_PT3D3','/home/zming/simulations/os2D/os_PT3D4','/home/zming/simulations/os2D/os_PT3D3/movie2', 1, start_num = 0, count_num=146)
    frame1.save_frames2()

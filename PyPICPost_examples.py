from outfile import *

def os_pond_scatter3D():
    file1 = OutFile(code_name='osiris',path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/os_pond_scatter3D1', field_name='raw', spec_name='e', out_num=5)
    file1.open()
    file1.read_raw_q()
    file1.read_raw_x1()
    file1.read_raw_x2()
    file1.read_raw_x3()
    file1.read_raw_p1()
    file1.read_raw_p2()
    file1.read_raw_p3()
    file1.read_raw_ene()
    file1.close()
    file1.select_raw_data(ene_low=0.1)
    plt.scatter(file1._raw_x1[file1._raw_select_index], file1._raw_x2[file1._raw_select_index])
    plt.show()

if __name__ == '__main__':
    os_pond_scatter3D()

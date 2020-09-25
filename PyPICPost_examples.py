import outfile
import matplotlib.pyplot as plt
import numpy as np

def qp_hi_compare_Ez(h_fig, h_axs, out_num=1):
    file1 = outfile.OutFile(code_name='quickpic',path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/qp_hi_compare/qp5', field_name='e1', use_num_list=True, out_num=out_num, fld_slice=1)
    file1.open()
    file1.read_data()
    file1.plot_data(h_fig, h_axs[0])
    data_qp = file1._data
    file1.read_data_lineout()
    lineout_qp = file1._data
    file1.close()

    file1 = outfile.OutFile(code_name='hipace', path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/hi_qp_compare/hi5', field_name='e1', use_num_list=True, out_num=out_num)
    file1.open()
    file1.read_data_slice()
    file1.plot_data(h_fig, h_axs[1])

    file1._data = file1._data - np.fliplr(data_qp)
    file1.plot_data(h_fig, h_axs[2])
    
    file1.read_data_lineout()
    file1.plot_data(h_fig, h_axs[2], c='k', ls='--', label='HiPACE')
    file1._data = np.flip(lineout_qp)
    file1.plot_data(h_fig, h_axs[2], c='r', ls=':', label='QuickPIC')
    file1.close()
    h_axs[2].legend()
    h_axs[2].set_title('Difference')

def qp_hi_compare_pzz(h_fig, h_axs, out_num=1):
    file1 = outfile.OutFile(code_name='hipace', path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/hi_qp_compare/hi4', field_name='raw', spec_name='driver', use_num_list=True, out_num=out_num)
    file1.open()
    file1.read_raw_p1()
    file1.read_raw_x1()
    ax_range = [[np.min(file1._raw_p1), np.max(file1._raw_p1)], [np.min(file1._raw_x1), np.max(file1._raw_x1)]]
    h_fig, h_ax = file1.plot_raw_hist2D(h_fig, h_axs[1], dims='p1x1', range=ax_range, if_log_colorbar=True, if_colorbar=True)
    hi_data = file1._data
    file1.close()

    file1 = outfile.OutFile(code_name='quickpic',path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/qp_hi_compare/qp4', field_name='raw', spec_name='Beam0001', use_num_list=True, out_num=out_num)
    file1.open(cell_size_qp_raw = np.array([0.03125, 0.03125, 0.01953125]))
    file1.read_raw_q()
    file1.read_raw_x1()
    file1.read_raw_p1()
    file1._raw_x1 = 10.-file1._raw_x1
    h_fig, h_ax = file1.plot_raw_hist2D(h_fig, h_axs[0], dims='p1x1', range=ax_range, if_reread=False, if_log_colorbar=True, if_colorbar=True)
    file1.close()

    file1._data = hi_data - file1._data
    file1.plot_data(h_fig, h_axs[2])

    file1.close()
    h_axs[2].set_title('Difference')

if __name__ == '__main__':
    h_fig, h_axs = plt.subplots(2, 3, figsize=[12.5, 7.])
    qp_hi_compare_pzz(h_fig, h_axs[0], out_num=0)
    qp_hi_compare_pzz(h_fig, h_axs[1], out_num=50)
    plt.tight_layout()
    plt.show()

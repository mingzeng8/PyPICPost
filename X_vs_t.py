import outfile
import numpy as np
import matplotlib.pyplot as plt

# plot what ever parameter vs time

def get_field_peak(outfile_object, if_abs=True):
    ''' Get the maximum of a field from outfile_object.
        outfile_object: a outfile.OutFile() object.
        if_abs: Boolean. If Ture, use absolute value of data.
    '''
    outfile_object.read_data_lineout()
    if if_abs:
        return np.max(np.abs(outfile_object._data))
    else:
        return np.max(outfile_object._data)

def plot_field_peak_vs_t(path, code_name='osiris', field_name='e3', if_abs=False, use_num_list=False, start=0, count=0, stride=1, t_offset = 0., h_ax=None, field_abs=True, **kwargs):
    '''
    
    '''
    t_array = list()
    peak_array = list()
    outfile_object = outfile.OutFile(code_name=code_name, path=path, field_name=field_name, use_num_list=use_num_list)
    for i in range(start, start+count*stride, stride):
        try:
            outfile_object.out_num=i
            outfile_object.open()
        except KeyError:
            print('Reaching the end of number list. Finishing...')
            break
        print('File number {}'.format(outfile_object.actual_num))
        peak = get_field_peak(outfile_object, if_abs=if_abs)
        outfile_object.close()
        t_array.append(outfile_object.time+t_offset)
        peak_array.append(peak)
    if h_ax is None: _, h_ax =  plt.subplots()
    h_ax.plot(t_array, peak_array, **kwargs)
    plt.xlabel('$t$')
    plt.minorticks_on()
    return h_ax

if __name__ == '__main__':
    h_ax = plot_field_peak_vs_t(code_name='fbpic', path='/home/zengming/project/FB/PTE_k28/diag_k=28_a0=15.7_w0=2.4_d=0', field_name='e2', if_abs=True, use_num_list=True, count=999, color='k', label='$d=0$')
    plot_field_peak_vs_t(code_name='fbpic', path='/home/zengming/project/FB/PTE_k28/diag_k=28_a0=15.7_w0=2.4_d=56', field_name='e2', if_abs=True, use_num_list=True, count=999, h_ax = h_ax, color='r', label='$d=56$')
    plot_field_peak_vs_t(code_name='fbpic', path='/home/zengming/project/FB/PTE_k28/diag_k=28_a0=15.7_w0=2.4_d=60', field_name='e2', if_abs=True, use_num_list=True, count=999, h_ax = h_ax, color='b', label='$d=60$')
    plot_field_peak_vs_t(code_name='fbpic', path='/home/zengming/project/FB/PTE_k28/diag_k=28_a0=15.7_w0=2.4_d=70', field_name='e2', if_abs=True, use_num_list=True, count=999, h_ax = h_ax, color='g', label='$d=70$')
    plot_field_peak_vs_t(code_name='fbpic', path='/home/zengming/project/FB/PTE_k28/diag_k=28_a0=15.7_w0=2.4_d=100', field_name='e2', if_abs=True, use_num_list=True, count=999, h_ax = h_ax, color='y', label='$d=100$')
    '''h_ax = plot_field_peak_vs_t(code_name='fbpic', path='/home/zengming/project/FB/PTE/diag_k=25_a0=20.3_w0=2.7_d=0', field_name='e2', if_abs=True, use_num_list=True, count=999, color='k', label='$d=0$')
    plot_field_peak_vs_t(code_name='fbpic', path='/home/zengming/project/FB/PTE/diag_k=25_a0=20.3_w0=2.7_d=80', field_name='e2', if_abs=True, use_num_list=True, count=999, h_ax = h_ax, color='r', label='$d=80$')
    h_ax = plot_field_peak_vs_t(code_name='fbpic', path='/home/zengming/project/FB/PTE/diag_k=25_a0=20.3_w0=2.7_d=50', field_name='e2', if_abs=True, use_num_list=True, count=999, color='y', label='$d=50$')
    plot_field_peak_vs_t(code_name='fbpic', path='/home/zengming/project/FB/PTE/diag_k=25_a0=20.3_w0=2.7_d=60', field_name='e2', if_abs=True, use_num_list=True, count=999, h_ax = h_ax, color='k', label='$d=60$')
    plot_field_peak_vs_t(code_name='fbpic', path='/home/zengming/project/FB/PTE/diag_k=25_a0=20.3_w0=2.7_d=70', field_name='e2', if_abs=True, use_num_list=True, count=99, h_ax = h_ax, color='g', label='$d=70$')
    plot_field_peak_vs_t(code_name='fbpic', path='/home/zengming/project/FB/PTE/diag_k=25_a0=20.3_w0=2.7_d=80', field_name='e2', if_abs=True, use_num_list=True, count=999, h_ax = h_ax, color='r', label='$d=80$')
    plot_field_peak_vs_t(code_name='fbpic', path='/home/zengming/project/FB/PTE/diag_k=25_a0=20.3_w0=2.7_d=90', field_name='e2', if_abs=True, use_num_list=True, count=999, h_ax = h_ax, color='b', label='$d=90$')'''
    plt.xlabel('$t$ [s]')
    plt.ylabel('$E_L$ [V/m]')
    plt.legend()
    plt.tight_layout()
    plt.show()

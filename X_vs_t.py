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

def plot_field_peak_vs_t(path, code_name='osiris', field_name='e3', if_abs=False, use_num_list=False, start=0, count=0, stride=1, t_offset = 0., h_ax=None, field_abs=True, linestyle='-'):
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
    h_ax.plot(t_array, peak_array, color='k', linestyle = linestyle)
    plt.xlabel('$t$')
    plt.minorticks_on()
    return h_ax

if __name__ == '__main__':
    plot_field_peak_vs_t(code_name='fbpic', path='/home/ming/mnt/CNG12/FB/DCII/diag_k=35_a0=1.5_w0=8.9_dphase=0.0_1mmdope', field_name='e2', if_abs=True, use_num_list=True, count=9, linestyle='-')
    plt.tight_layout()
    plt.show()

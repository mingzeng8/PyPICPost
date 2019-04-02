import outfile
import numpy as np
import matplotlib.pyplot as plt
import time

def plot_W_vs_t(path, field_name='e3', species_name='', average='', cyl_m='False', start=0, count=0, stride=1):
    '''Plot the W evolution'''
    if not isinstance(path, str):
        raise TypeError('The path should be a string!')
    if not isinstance(field_name, str):
        raise TypeError('The field name should be a string!')
    if not isinstance(species_name, str):
        raise TypeError('The species name should be a string!')
    if not isinstance(start, int):
        raise TypeError('The start should be a integer!')
    if not isinstance(count, int):
        raise TypeError('The count should be a integer!')
    if not isinstance(stride, int):
        raise TypeError('The stride should be a integer!')
    if cyl_m:
        field_name=field_name+'_cyl_m'
    os_file = outfile.OutFile(path=path, field_name=field_name, average=average, cyl_m_num=1, cyl_m_re_im='re')
    if start not in range(10**os_file.digit_num):
        raise ValueError('start = {0} not in range!'.format(start))
    if count not in range(10**os_file.digit_num):
        raise ValueError('count = {0} not in range!'.format(count))
    if stride not in range(10**os_file.digit_num):
        raise ValueError('stride = {0} not in range!'.format(stride))

    W_array = list()
    t_array = list()
    a_array = list()
    try:
        popt = None
        for i in range(start, start+count*stride, stride):
            os_file.out_num = i
            try:
                os_file.open()
            except IOError:
                print('Warning: unable to open file \'{0}\'. Iteration breaks.'.format(os_file.path_filename))
                break
            os_file.read_data_project(if_abs = True)
            popt, h_fig, h_ax = os_file.fit_for_W(guess_values=popt)
            W_array.append(os_file._W)
            t_array.append(os_file.time)
            a_array.append(os_file._a)
            plt.show(block=False)
            time.sleep(1)
            plt.close(h_fig)
            #h_ax=None
            os_file.close()
    except KeyboardInterrupt:
        print('Keyboard interruption occurs at file \'{0}\'. Iteration breaks.'.format(os_file.path_filename))

    W_max = max(W_array)
    W_max_i = W_array.index(W_max)
    print('W_max = {0} at index {1}, t = {2}'.format(W_max, W_max_i, t_array[W_max_i]))
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(t_array, W_array)
    plt.xlabel('t [$\\omega_p$]')
    plt.ylabel('w [$c / \\omega_p$]')
    return ([[i*stride+start, t_array[i], W_array[i], a_array[i]] for i in range(len(t_array))], fig1)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Measure the laser size W vs t.')
    parser.add_argument('-s', '--savg', help="use space average [default false]", action='store_true', default=False)
    parser.add_argument('-c', '--cyl_m', help="use 2D cylindrical modes [default false]", action='store_true', default=False)
    parser.add_argument("fnum", help="plot w vs z for a given simulation",type=str)
    args=parser.parse_args()
    fnum=args.fnum
    if parser.parse_args().savg:
        average='-savg'
    else:
        average=''
    #prefix='PT3D'
    prefix='laser3D'
    print('Working on os_'+prefix+fnum+'...')
    #if fnum[-1]=='h':
    #    average='-savg'
    t_W_array, h_fig = plot_W_vs_t('/home/zming/simulations/os2D/os_'+prefix+fnum, 'e3', average=average, cyl_m=parser.parse_args().cyl_m, start=0, count=600, stride=1)
    print(t_W_array)
    plot_save_name='/home/zming/simulations/os2D/os_{0}{1}/{0}{1}.png'.format(prefix,fnum)
    plt.savefig(plot_save_name)
    print('Plot saved at '+plot_save_name)
    data_save_name='/home/zming/simulations/os2D/os_{0}{1}/{0}{1}.data'.format(prefix,fnum)
    np.savetxt(data_save_name, np.array(t_W_array), delimiter=',', fmt='%.5f')
    print('Data saved at '+data_save_name)
    plt.show(block=True)
    

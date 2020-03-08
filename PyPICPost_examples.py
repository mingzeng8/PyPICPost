from outfile import *

def os_pond_scatter3D():
    file1 = OutFile(code_name='osiris',path='/beegfs/desy/group/fla/plasma/OSIRIS-runs/2D-runs/MZ/os_pond_scatter3D/os_pond_scatter3D2', field_name='raw', spec_name='e', out_num=9)
    file1.open()
    file1.read_raw_q()
    file1.read_raw_x1()
    file1.read_raw_x2()
    file1.read_raw_x3()
    file1.read_raw_p1()
    file1.read_raw_p2()
    #file1.read_raw_p3()
    file1.read_raw_ene()
    file1.close()
    file1.select_raw_data(ene_low=0.1, r_low=1.6, x1_low=-2.)
    plt.scatter(file1._raw_x1[file1._raw_select_index], file1._raw_x2[file1._raw_select_index])
    # Get the filtered particles
    weight = np.abs(file1._raw_q[file1._raw_select_index])
    p_z_prime = file1._raw_p1[file1._raw_select_index]
    p_x_prime = file1._raw_p2[file1._raw_select_index]
    #p_y_prime = file1._raw_p3[file1._raw_select_index]
    ene = file1._raw_ene[file1._raw_select_index]
    p_square = np.square(ene+1.)-1.
    # Generate histogram for N (number of particles) vs. phi_M and theta
    phi_M = np.linspace(.55, .9, 64)
    theta_degree = np.linspace(20., 140., 64)
    # Transform theta from degree to rad
    theta = theta_degree * (np.pi/180)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    N = np.zeros((phi_M.size, theta.size))
    pzmax=0.
    for theta_ind in range(theta.size):
        # Rotate axis for particle momentum
        p_z = cos_theta[theta_ind]*p_z_prime + sin_theta[theta_ind]*p_x_prime
        if p_z.max() > pzmax: pzmax = p_z.max()
        p_perp_square = p_square - np.square(p_z)
        for phi_M_ind in range(phi_M.size):
            # Same as p_z >= (1.-np.square(phi_M[phi_M_ind])+p_perp_square)/phi_M[phi_M_ind]/2 but faster
            #condition = p_z >= (1.-np.square(phi_M[phi_M_ind])+p_perp_square)/phi_M[phi_M_ind]/2
            condition = ((p_z - p_perp_square/(phi_M[phi_M_ind]*2)) >= (1.-np.square(phi_M[phi_M_ind]))/phi_M[phi_M_ind]/2)
            #condition = p_z >= 1.-phi_M[phi_M_ind]
            #if any(condition): print("True")
            N[phi_M_ind, theta_ind] = np.sum(weight[condition])
    #print(pzmax)
    #h_fig = plt.figure()
    #plt.imshow(N)
    h_fig = plt.figure()
    h_ax = h_fig.add_subplot(111)
    h_plot = h_ax.pcolormesh(theta_degree, phi_M, N, cmap=my_cmap.cmap_lower_range_transparent(transparency_transition_region=[0.,0.02]))
    h_ax.set_ylabel('$\\phi_M$')
    h_ax.set_xlabel('$\\theta$ [$^{{\\circ}}$]')
    h_cb = plt.colorbar(h_plot, ax=h_ax)
    h_cb.set_label('$N$ [arb. unit]')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    os_pond_scatter3D()

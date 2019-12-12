from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import numpy as np

def cmap_higher_range_transparent(original_cmap=plt.cm.jet, transparency_transition_region=[0.8,1.0]):
    ''' creat a colormap based on original_cmap, with a linear transparency transition from (transparency_start_ratio*total range) to the higher end.'''
    my_cmap = original_cmap(np.arange(original_cmap.N))
    transition_ind=[int(transparency_transition_region[0]*original_cmap.N), int(transparency_transition_region[1]*original_cmap.N)]
    my_cmap[transition_ind[1]:,-1] = 0.
    my_cmap[transition_ind[0]:transition_ind[1],-1] = np.linspace(1, 0, transition_ind[1]-transition_ind[0])
    return ListedColormap(my_cmap)

def cmap_middle_range_transparent(original_cmap=plt.cm.bwr, transparency_start_ratio=(0.05, 0.2)):
    ''' creat a colormap based on original_cmap, with linear transparency transitions from ((0.5-transparency_start_ratio[1])*total range) to ((0.5-transparency_start_ratio[0])*total range) and from ((0.5+transparency_start_ratio[0])*total range) to ((0.5+transparency_start_ratio[1])*total range).'''
    my_cmap = original_cmap(np.arange(original_cmap.N))
    transparency_index=[int((0.5-transparency_start_ratio[1])*original_cmap.N), int((0.5-transparency_start_ratio[0])*original_cmap.N), int((0.5+transparency_start_ratio[0])*original_cmap.N), int((0.5+transparency_start_ratio[1])*original_cmap.N)]
    my_cmap[transparency_index[0]:transparency_index[1],-1] = np.linspace(1, 0, transparency_index[1]-transparency_index[0])
    my_cmap[transparency_index[1]:transparency_index[2],-1] = np.zeros(transparency_index[2]-transparency_index[1])
    my_cmap[transparency_index[2]:transparency_index[3],-1] = np.linspace(0, 1, transparency_index[3]-transparency_index[2])
    return ListedColormap(my_cmap)

def cmap_lower_range_transparent(original_cmap=plt.cm.jet, transparency_transition_region=[0.,0.2]):
    ''' creat a colormap based on original_cmap, with a linear transparency transition from (transparency_transition_region[0]*total range) to (transparency_transition_region[1]*total range).
        One has to make sure transparency_transition_region[0]<transparency_transition_region[1], and both of these two elemenst are between 0 and 1.
    '''
    my_cmap = original_cmap(np.arange(original_cmap.N))
    transition_ind=[int(transparency_transition_region[0]*original_cmap.N), int(transparency_transition_region[1]*original_cmap.N)]
    my_cmap[0:transition_ind[0],-1] = 0. # same effect as ...=np.zerps(transition_ind[0], dtype=float)
    my_cmap[transition_ind[0]:transition_ind[1],-1] = np.linspace(0, 1, transition_ind[1]-transition_ind[0])
    return ListedColormap(my_cmap)

def combine_cmaps(cmap1, cmap2, ratio = 0.5, N = 256):
    ''' linearly combine two colormaps. ratio determines the portion of the two colormaps. N is the total color bins of the output colormap.'''
    if ratio <= 0. or ratio >= 1.:
        raise ValueError("Ratio should be in the range (0., 1.).")
    N1 = int(N*ratio)
    N2 = N - N1
    colors1 = cmap1(np.linspace(0., 1., N1))
    colors2 = cmap2(np.linspace(0., 1., N2))
    colors = np.vstack((colors1, colors2))
    return ListedColormap(colors)

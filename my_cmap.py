from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import numpy as np

def cmap_higher_range_transparent(original_cmap=plt.cm.jet, transparency_start_ratio=0.8):
    ''' creat a colormap based on original_cmap, with a linear transparency transition from (transparency_start_ratio*total range) to the higher end.'''
    my_cmap = original_cmap(np.arange(original_cmap.N))
    transparency_start=int(transparency_start_ratio*original_cmap.N)
    my_cmap[transparency_start:original_cmap.N,-1] = np.linspace(1, 0, original_cmap.N-transparency_start)
    return ListedColormap(my_cmap)

def cmap_middle_range_transparent(original_cmap=plt.cm.bwr, transparency_start_ratio=(0.05, 0.2)):
    ''' creat a colormap based on original_cmap, with linear transparency transitions from ((0.5-transparency_start_ratio[1])*total range) to ((0.5-transparency_start_ratio[0])*total range) and from ((0.5+transparency_start_ratio[0])*total range) to ((0.5+transparency_start_ratio[1])*total range).'''
    my_cmap = original_cmap(np.arange(original_cmap.N))
    transparency_index=[int((0.5-transparency_start_ratio[1])*original_cmap.N), int((0.5-transparency_start_ratio[0])*original_cmap.N), int((0.5+transparency_start_ratio[0])*original_cmap.N), int((0.5+transparency_start_ratio[1])*original_cmap.N)]
    my_cmap[transparency_index[0]:transparency_index[1],-1] = np.linspace(1, 0, transparency_index[1]-transparency_index[0])
    my_cmap[transparency_index[1]:transparency_index[2],-1] = np.zeros(transparency_index[2]-transparency_index[1])
    my_cmap[transparency_index[2]:transparency_index[3],-1] = np.linspace(0, 1, transparency_index[3]-transparency_index[2])
    return ListedColormap(my_cmap)

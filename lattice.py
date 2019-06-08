import numpy as np
import matplotlib.pyplot as plt

# useful dicts
planar = {'r': '1000',
          'u': '0100',
          'l': '0010',
          'd': '0001'}

k_dict = {'r': 0,
          'u': 1,
          'l': 2,
          'd': 3}

move_dict = {'r': [+1, 0],
             'u': [ 0,+1],
             'l': [-1, 0],
             'd': [ 0,-1]}

mark_dict = {'r': '>',
             'u': '^',
             'l': '<',
             'd': 'v'}


# `turns`: string e.g. "lurrddr"
# `pep`: string e.g. "PSVK"
def plot_lattice(turns, pep=None):
    
    coords = np.zeros((len(turns)+1,2), dtype=int)
    coords[0] = [0,0]
    for t,turn in enumerate(turns):
        coords[t+1] = coords[t] + move_dict[turn]

    n_coords = len(coords)
    max_v = abs(coords).max()+1
    plt.figure(figsize=(n_coords*1.1,n_coords*1.1))
    plt.xlim(-n_coords,n_coords)
    plt.ylim(-n_coords,n_coords)
    plt.hlines(np.arange(-n_coords,n_coords,1),-n_coords,n_coords, color='lightgrey',linestyle='-.', lw=.7)
    plt.vlines(np.arange(-n_coords+1,n_coords,1),-n_coords,n_coords, color='lightgrey',linestyle='-.', lw=.7)
    plt.axis('off')
    for c,coord in enumerate(coords):
        if c < len(coords)-1:
            t = turns[c]
            plt.arrow(*coord,*move_dict[t],
                      width=.13,head_width=0,length_includes_head=True,
                      color='lightgrey')
            plt.scatter(coords[c,0]+move_dict[t][0]*.07,coords[c,1]+move_dict[t][1]*.07,
                        s=350, facecolors='lightgrey', marker=mark_dict[t],facecolor='lightgrey')
        plt.scatter(coords[-1,0],coords[-1,1],
                        s=250, marker='o', facecolor='lightgrey')
            
        plt.annotate(s=c, xy=coord+[-.12,.35],
                     horizontalalignment='right',verticalalignment='center', 
                     color='grey',fontfamily='monospace')
        if pep:
            plt.text(coord[0], coord[1], pep[c],
                     horizontalalignment='center',verticalalignment='center',
                     color='darkblue', fontfamily='monospace')
        
        def space(ins):
            return ' '.join([str(x).upper() for x in ins])
        for i,x in enumerate([pep, range(n_coords), turns]):
            plt.text(n_coords,n_coords-i/1.5, space(x), fontfamily='monospace', color='darkblue')
        
    plt.show()
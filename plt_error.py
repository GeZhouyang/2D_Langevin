import sys
import numpy as np
from matplotlib import pyplot as plt

params = {'text.usetex' : True,
          'font.size' : 12.5,
          'font.family' : 'lmodern',
          'text.latex.unicode' : True}

plt.rcParams['text.latex.preamble']=[r'\usepackage{lmodern}']
plt.rcParams.update(params)
    

#---- define a function to plot data

def plot_data(ax, t,x):

    # data
    c0 = ax.loglog(t, x, 'o', mec='C0', mfc='C0')

    # slopes
    xx = np.array([1e-2,1e1])
    c1 = ax.plot(xx, 0.005*xx, '--')
    c2 = ax.plot(xx, 0.015*np.sqrt(xx), '--')

    ax.set_xlim(3e-2, 3e0)
    ax.set_ylim(2e-4, 3e-2)

    ax.set_xlabel(r'$\Delta t$')
    ax.set_ylabel(r'$\langle r \rangle$')

    ax.get_xaxis().set_tick_params(which='both', direction='in')
    ax.get_yaxis().set_tick_params(which='both', direction='in')

    ax.legend([r'Data',r'$\propto \Delta t$',r'$\propto \sqrt{\Delta t}$'],
              frameon=False, loc=2, numpoints=1, fontsize=11)
    
    return ax


#---- Script

if __name__ == '__main__':

    #---> Load data

    dt   = np.loadtxt('data.txt', usecols=(0,), skiprows=0)
    disp = np.loadtxt('data.txt', usecols=(1,), skiprows=0)

    #---> Create the figure
    
    fig, ax = plt.subplots(1, 1, sharey=True)
    fig.set_size_inches(8, 8)
    
    ax = plot_data(ax, dt,disp)
        
    plt.savefig('error.png', bbox_inches='tight', transparent=False)
    plt.show()

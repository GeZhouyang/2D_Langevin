import sys
import numpy as np
from matplotlib import pyplot as plt

params = {'text.usetex' : True,
          'font.size' : 12.5,
          'font.family' : 'lmodern',
          'text.latex.unicode' : True}

plt.rcParams['text.latex.preamble']=[r'\usepackage{lmodern}']
plt.rcParams.update(params)

pi  = np.pi

def plot_traj(ax,x,y):

    il = len(x)
    for i in range(il):
        ax.plot(x[i],y[i])

    xy = 1.5
    ax.set_xlim(-xy, xy)
    ax.set_ylim(-xy, xy)
    
    ax.set_xlabel(r'$x/a$')
    ax.set_ylabel(r'$y/a$')

    ax.get_xaxis().set_tick_params(which='both', direction='in')
    ax.get_yaxis().set_tick_params(which='both', direction='in')
    
    return ax

def plot_scat(ax, x,y):

    ax.scatter(x,y, alpha=0.5)

    mm = 20
    ax.set_xlim(-mm, mm)
    ax.set_ylim(-mm, mm)
    
    ax.set_xlabel(r'$x/a$')
    ax.set_ylabel(r'$y/a$')
    
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    ax.get_yaxis().set_tick_params(which='both', direction='in')

    #plt.axis('off')

    return ax

def plot_hist(ax, r,Nx,t):

    mm = 8
    bin_size = 0.1
    
    ax.hist(r, bins=np.arange(0,mm,bin_size))

    #-- analytical solution of the diffusion eqn --#
    
    if t != 0.:
        
        rr = np.arange(bin_size/2.,mm+bin_size/2.,bin_size)  # same bin but centered value (more accurate)
        nn = Nx/(4.*pi*D_coef*t)*np.exp( - (a*rr)**2/(4.*D_coef*t) )  # non-dimensionalize x by a

        for i in range(len(rr)):
            nn[i] *= pi*2.*rr[i]*bin_size*a**2  # multiply by the area to get the count
            
        ax.plot(np.insert(rr,0,0), np.insert(nn,0,0)) # insert 0 in the beginning

    ax.set_xlim(0, mm)
    ax.set_ylim(0, Nx/5)
    
    ax.set_xlabel(r'$r/a$')
    ax.set_ylabel(r'$N$')
    
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    ax.get_yaxis().set_tick_params(which='both', direction='in')

    return ax


def ensemble_average(x_all,Nx):
    
    return np.sum(x_all, axis=0)/Nx


def solve_langevin(x0,y0,G_scale,Nt):

    method_num = 1
    
    mu    = 0.  # mean
    sigma = 1.  # standard deviation

    if method_num == 1: # 2D Gaussian step
        Gstep2d = np.random.multivariate_normal([mu,mu],[[sigma,0],[0,sigma]],Nt)*G_scale
    
    xt,yt = x0,y0
    xt_all,yt_all = [x0],[y0]  # initialize the position lists
    
    for i in range(Nt):  # integrate in *time*

        if method_num == 0:  # faster but less accurate

            Gstep = np.random.normal(mu,sigma)*G_scale    # Gaussian step
            angle = np.random.uniform(0.,2.*pi)           # random angle (uniform)
        
            xt += Gstep*np.cos(angle)       # random walk (x-dir)
            yt += Gstep*np.sin(angle)       # random walk (y-dir)

        elif method_num == 1:  # accurate

            xt += Gstep2d[i][0]       # random walk (x-dir)
            yt += Gstep2d[i][1]       # random walk (y-dir)
        
        xt_all.append(xt)  
        yt_all.append(yt)

    return [xt_all,yt_all]



def main_compute(G_scale,Nt,t, x0,y0,Nx):

    write_disp_to_file  = False
    plot_trajectories   = False
    plot_scat_snapshots = True
    plot_scat_meansnaps = False
    plot_hist_snapshots = False

    x_all,y_all = [],[] # list of samples

    for i in range(Nx):  # sampling the *ensemble*
        
        [xt_all,yt_all] = solve_langevin(x0,y0,G_scale,Nt)
        
        x_all.append(xt_all)
        y_all.append(yt_all)

    # list of averaged positions in time
    x_avr = ensemble_average(x_all,Nx)
    y_avr = ensemble_average(y_all,Nx)

    disp_last_avr = np.sqrt(x_avr[-1]**2+y_avr[-1]**2)
    print 'Averaged final displacement r = ', disp_last_avr

    if write_disp_to_file:

        f=open("data.txt", "a+")  # append to file
        f.write( '\n%0.2e %0.6e' % (dt, disp_last_avr) )
        f.close()
        
    if plot_trajectories:

        fig, ax = plt.subplots(1, 1, sharey=True)
        fig.set_size_inches(8,8)
    
        ax = plot_traj(ax, x_all,y_all)

        plt.savefig('trajectory_one.png', bbox_inches='tight', transparent=False)
        plt.show()

    if plot_scat_snapshots or plot_scat_meansnaps or plot_hist_snapshots:
        
        t_inserted = np.insert(t,0,0)  # insert 0 to the beginning of t
        x_transpose = map(list, zip(*x_all))
        y_transpose = map(list, zip(*y_all))

        for i in range(len(t_inserted)):

            fig, ax = plt.subplots(1, 1, sharey=True)
            fig.set_size_inches(4,4)#(8,8)

            if plot_scat_snapshots:

                ax = plot_scat(ax, x_transpose[i],y_transpose[i])
                plt.savefig('scat-'+('%0.2e' % t[0])+'_only1/scat'+str(i).zfill(6)+'.png',
                            bbox_inches='tight', transparent=False)

            if plot_scat_meansnaps:

                ax = plot_scat(ax, x_avr[i],y_avr[i])
                plt.savefig('scat-'+('%0.2e' % t[0])+'_mean/scat'+str(i).zfill(6)+'.png',
                            bbox_inches='tight', transparent=False)        

            if plot_hist_snapshots:
        
                disp = []
                for j in range(Nx):
                    bufx = x_transpose[i]
                    bufy = y_transpose[i]
            
                    r0 = np.sqrt(bufx[j]**2 + bufy[j]**2)
                    disp.append(r0)
            
                ax = plot_hist(ax, disp,Nx,t_inserted[i])
                plt.savefig('hist-'+('%0.2e' % t[0])+'/hist'+str(i).zfill(6)+'.png',
                            bbox_inches='tight', transparent=False)       
    

    return [x_all,y_all]


#---- Main script ----#


if __name__ == '__main__':
    
    #--- init parameters ---#

    x0,y0 = 0.,0.   # starting position [1]
    t0 = 0.         # starting time [s]
    t1 = 20.        # end time [s]
    Nx = 5000      # number of samples (the ensemble)

    T   = 300.      # temperature (Kelvin)
    eta = 1e-3      # dynamic viscosity of water [Pa*s]
    a   = 1e-6      # particle radius [m]
    k_B = 1.38e-23  # Boltzmann const [m2*kg/s^2/K] 

    alpha = 6.*pi*eta*a    # Stokes drag coef
    D_coef = k_B*T/alpha   # diffusion coef
    
    #dtl = [2,1,5e-1,2e-1,1e-1,5e-2,2e-2]   # time step
    dtl = [2e-1]
    
    #--- main solve ---#
    
    print '\nSolving the Langevin equation at T =',T, '[K],'
    print 'corresponding to diffusion coef = ',D_coef, '[m^2/s].'
        
    for dt in dtl:

        # G_scale is the (non-dimensional) standard deviation
        G_scale = np.sqrt(2.*D_coef*dt)/a  

        print '\nTime step dt = ',dt
        print 'The scale of the Gaussian = ',G_scale

        ##t1 = dt  # only if testing one step
        Nt = int((t1-t0)/dt)            # total number of time steps
        t = np.linspace(dt,t1, num=Nt)  # list of discrete time

        [x_all,y_all] = main_compute(G_scale,Nt,t, x0,y0,Nx)
        

    print '\nDone.'

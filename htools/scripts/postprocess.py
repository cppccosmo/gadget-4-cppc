import h5py
import glob
import numpy as np
import matplotlib.pyplot as plt
import cmasher as cmr
import MAS_library as MASL

class Snap():
    def __init__(self, fname):
        self.fname = fname
        self.f = h5py.File(fname+'.hdf5','r')
        self.z = self.f['Header'].attrs['Redshift']
        self.a = self.f['Header'].attrs['Time']
        self.L = self.f['Header'].attrs['BoxSize']
        self.grid  = self.f['Parameters'].attrs['GridSize']

    @property
    def redshift(self):
        return self.z

    @property
    def boxsize(self):
        return self.L

    def get_vels(self, ptype, n=1000):
        vels = self.f['PartType%d'%ptype]
        vels = vels['Velocities'][:n]
        absv = np.sqrt(vels[:,0]**2+vels[:,1]**2+vels[:,2]**2)
        return absv*self.a**(3/2)

    def get_delta(self, ptype, axis=0, threads=4, verbose=True):
        ptypes    = [ptype]
        MAS       = 'CIC'
        verbose   = verbose
        threads   = threads
        delta     = MASL.density_field_gadget(self.fname, ptypes, self.grid, MAS, do_RSD=False, axis=axis, verbose=verbose)
        return delta

    def projection(self, ptype, axis=0, threads=4, verbose=False, vmin=-1,vmax=2):
        ptypes    = [ptype]
        MAS       = 'CIC'
        verbose   = verbose
        threads   = threads
        delta     = MASL.density_field_gadget(self.fname, ptypes, self.grid, MAS, do_RSD=False, axis=axis, verbose=verbose)
        fig,ax = plt.subplots(1,1,figsize=(13,12))
        im = ax.imshow(np.log10(np.mean(delta,axis=axis)),
               cmap='cmr.pride',origin='lower',
               extent=[-self.L/2,self.L/2,-self.L/2,self.L/2],vmin=vmin,vmax=vmax)
        plt.colorbar(im)
        ax.set_title(r'$z=%d$'%(round(self.z)))
        ax.set_xlabel(r'$x~({\rm Mpc}/h)$')
        ax.set_ylabel(r'$y~({\rm Mpc}/h)$')
        plt.show()

class PowerSpectra:
    def __init__(self, sdir, verbose=False):
        self.mf_dir = sdir + '/output'
        self.c1_dir = sdir + '/output_c1'
        self.c2_dir = sdir + '/output_c2'
        self.c3_dir = sdir + '/output_c3'
        self.c4_dir = sdir + '/output_c4'
        self.Ntau = 20
        self.Nf = len([f for f in glob.iglob(sdir +'/output/snapshot_*', recursive=True)])
        self.vlist  = np.loadtxt(sdir+'/vel_list_full.txt')
        if verbose:
            print('Found %d spectra to analyse'%self.Nf)

        # Multifluid 
        de_dir = self.mf_dir + '/neutrino_stream_data/neutrino_delta_stream_'
        th_dir = self.mf_dir + '/neutrino_stream_data/neutrino_theta_stream_'

        de = [np.loadtxt(de_dir + '%.3d.csv'%i ,delimiter=',',usecols=range(0,self.Ntau+1)) for i in range(self.Nf)]
        th = [np.loadtxt(th_dir + '%.3d.csv'%i ,delimiter=',',usecols=range(0,self.Ntau+1)) for i in range(self.Nf)]

        self.k_mf = de[0][:,0]
        self.de_mf = de
        self.th_mf = th

        # CB particles 

    @property
    def mf_delta_flow(self, snap, alpha):
        return self.de_mf[snap][:, alpha+1]

    @property
    def mf_delta_flow0(self, alpha):
        return self.de_mf[-1][:, alpha+1]

    @property
    def mf_theta_flow(self, snap, alpha):
        return self.th_mf[snap][:, alpha+1]

    @property
    def mf_theta_flow0(self, alpha):
        return self.th_mf[-1][:, alpha+1]









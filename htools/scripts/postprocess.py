import os, glob
import h5py
import numpy as np
import matplotlib.pyplot as plt
import cmasher as cmr
import Pk_library as PKL
import MAS_library as MASL

THREADS = len(os.sched_getaffinity(0))

class Snap():
    def __init__(self, fname):
        self.fname = fname
        self.f = h5py.File(fname+'.hdf5','r')
        self.z = self.f['Header'].attrs['Redshift']
        self.a = self.f['Header'].attrs['Time']
        self.L = self.f['Header'].attrs['BoxSize']
        self.grid = int(self.f['Config'].attrs['PMGRID'])

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

    def get_delta(self, ptype, grid=512, axis=0, threads=THREADS, verbose=True):
        ptypes    = ptype
        MAS       = 'CIC'
        verbose   = verbose
        threads   = threads
        delta     = MASL.density_field_gadget(self.fname, ptypes, grid, MAS, do_RSD=False, axis=axis, verbose=verbose)
        return delta/np.mean(delta)

    def projection(self, ptype, axis=0, threads=4, verbose=False, vmin=-1,vmax=3):
        delta = self.get_delta(ptype, axis=axis, threads=threads, verbose=verbose)
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
    def __init__(self, sdir):
        self.sdir   = sdir
        self.mf_dir = sdir + '/output'
        self.gs_dir = sdir + '/output_gse'
        self.c1_dir = sdir + '/output_c1'
        self.c2_dir = sdir + '/output_c2'
        self.c3_dir = sdir + '/output_c3'
        self.c4_dir = sdir + '/output_c4'

        self.L = int(np.loadtxt(sdir + '/output/powerspecs/powerspec_000.txt', max_rows=3)[-1])

        c0   = int(np.loadtxt(sdir + '/output/powerspecs/powerspec_000.txt', max_rows=2)[-1])
        self.k_c0   = np.loadtxt(sdir + '/output/powerspecs/powerspec_000.txt', skiprows=5, max_rows=c0)[:,0]

        try:
            c1   = int(np.loadtxt(sdir + '/output_c1/powerspecs/powerspec_000.txt', max_rows=2)[-1])
            self.k_c1   = np.loadtxt(sdir + '/output_c1/powerspecs/powerspec_000.txt', skiprows=5, max_rows=c1)[:,0]
        except:
            print("Missing conversion c1")
        try:
            c2   = int(np.loadtxt(sdir + '/output_c2/powerspecs/powerspec_000.txt', max_rows=2)[-1])
            self.k_c2   = np.loadtxt(sdir + '/output_c2/powerspecs/powerspec_000.txt', skiprows=5, max_rows=c2)[:,0]
        except:
            print("Missing conversion c2")

        try:
            c3   = int(np.loadtxt(sdir + '/output_c3/powerspecs/powerspec_000.txt', max_rows=2)[-1])
            self.k_c3   = np.loadtxt(sdir + '/output_c3/powerspecs/powerspec_000.txt', max_rows=2)[-1]
        except:
            print("Missing conversion c3")

        try:
            c4   = int(np.loadtxt(sdir + '/output_c4/powerspecs/powerspec_000.txt', max_rows=2)[-1])
            self.k_c4   = np.loadtxt(sdir + '/output_c4/powerspecs/powerspec_000.txt', max_rows=2)[-1]
        except:
            print("Missing conversion c4")

        # Multifluid 
        self.Ntau = 20
        self.Nf = len([f for f in glob.iglob(sdir +'/output/snapshot_*', recursive=True)])
        self.vlist  = np.loadtxt(sdir+'/vel_list_full.txt')

        de_dir = self.mf_dir + '/neutrino_stream_data/neutrino_delta_stream_'
        th_dir = self.mf_dir + '/neutrino_stream_data/neutrino_theta_stream_'

        de = [np.loadtxt(de_dir + '%.3d.csv'%i ,delimiter=',',usecols=range(0,self.Ntau+1)) for i in range(self.Nf)]
        th = [np.loadtxt(th_dir + '%.3d.csv'%i ,delimiter=',',usecols=range(0,self.Ntau+1)) for i in range(self.Nf)]

        self.k_mf = de[0][:,0]
        self.de_mf = de
        self.th_mf = th
        self.De_mf  = np.mean([self.mf_De_flow0(i) for i in range(self.Ntau)],axis=0)
        self.De_mf_i  = np.mean([self.mf_De_flow0_i(i) for i in range(self.Ntau)],axis=0)

    # MF delta, theta, Delta^2
    def mf_delta_flow(self, snap, alpha):
        return self.de_mf[snap][:, alpha+1]

    def mf_delta_flow0(self, alpha):
        return self.de_mf[-1][:, alpha+1]

    def mf_theta_flow(self, snap, alpha):
        return self.th_mf[snap][:, alpha+1]

    def mf_theta_flow0(self, alpha):
        return self.th_mf[-1][:, alpha+1]

    def mf_De_flow(self, snap, alpha):
        return self.k_mf**3/(2*np.pi**2)*self.mf_delta_flow(snap, alpha)**2

    def mf_De_flow0(self, alpha):
        return self.k_mf**3/(2*np.pi**2)*self.mf_delta_flow0(alpha)**2

    def mf_De_flow_i(self, snap, alpha):
        return np.interp(self.k_c0, self.k_mf, self.mf_De_flow(snap, alpha))

    def mf_De_flow0_i(self, alpha):
        return np.interp(self.k_c0, self.k_mf, self.mf_De_flow0(alpha))

    # Particle type generic
    def ptype_De(self, conv, snap, ptype, smooth=False, grid=1024, gse=False):
        if conv == 0:
            if gse:
                print("GSE spectrum")
                ddir = self.gs_dir
            ddir = self.mf_dir
            f = np.loadtxt(ddir + '/powerspecs/powerspec_%.3d.txt'%snap, max_rows=5)
            z = int(np.round(1/f[0]-1))
            rows = int(f[1])
            print('Type 1 (CDM) spectrum at z=%d'%z)
            if smooth:
                snapshot = Snap(self.sdir + '/output/snapshot_%.3d'%snap)
                delta  = snapshot.get_delta([ptype], grid=grid, threads=THREADS, verbose=False)
                Pk = PKL.Pk(delta, self.L, axis=0, MAS="CIC", threads=THREADS, verbose=False)
                return Pk.k3D, Pk.Pk[:,0]*Pk.k3D**3/(2*np.pi**2), np.sqrt(Pk.Pk[:,0])
            else:
                f = np.loadtxt(ddir + '/powerspecs/powerspec_%.3d.txt'%snap, skiprows=5, max_rows=rows)
                return f[:,0], f[:,1], f[:,2]
        else:
            ddir = self.sdir + '/output_c%d'%conv
            f = np.loadtxt(ddir + '/powerspecs/powerspec_type%d_%.3d.txt'%(ptype, snap), max_rows=5)
            z = int(np.round(1/f[0]-1))
            rows = int(f[1])
            print('Type %d spectrum at z=%d'%(ptype,z))
            if smooth:
                snapshot = Snap(self.sdir + '/output_c%d/snapshot_%.3d'%(conv,snap))
                delta  = snapshot.get_delta([ptype], grid=grid, threads=THREADS, verbose=False)
                Pk = PKL.Pk(delta, self.L, axis=0, MAS="CIC", threads=THREADS, verbose=False)
                return Pk.k3D, Pk.Pk[:,0]*Pk.k3D**3/(2*np.pi**2), np.sqrt(Pk.Pk[:,0])
            else:
                f = np.loadtxt(ddir + '/powerspecs/powerspec_type%d_%.3d.txt'%(ptype, snap), skiprows=5, max_rows=rows)
                return f[:,0], f[:,1], f[:,2]

    def De_hdm_hy(self, deg_list):
        pass


    def P_hdm(self, deg_list):
        Nco = len(deg_list)
        fu = 0
        for i in range(Nco):
            fu += deg_list[i]
        p_hy = np.mean([np.loadtxt(self.sdir+'/postprocess/ps_hy_type%d_raw.txt'%(i+2)) for i in range(4)],axis=0)
        k    = p_hy[:,0] # Normalise to raw spectrum k
        conv = p_hy[:,1] # Dimensionless Delta^2
        unc_raw = np.mean([self.mf_De_flow0(j) for j in range(fu, self.Ntau)],axis=0)
        unc  = np.interp(k, self.k_mf, unc_raw)
        return np.average(np.array([conv,unc]), axis=0, weights=[fu,self.Ntau-fu])




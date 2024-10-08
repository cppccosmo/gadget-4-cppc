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
        if fu == self.Ntau:
            res = conv
        else:
            unc_raw = np.mean([self.mf_De_flow0(j) for j in range(fu, self.Ntau)],axis=0)
            unc  = np.interp(k, self.k_mf, unc_raw)
            res = np.average(np.array([conv,unc]), axis=0, weights=[fu,self.Ntau-fu])
        return res


def MF_cb(ddir):
    PS = PowerSpectra(ddir)
    f = np.sort([f for f in glob.iglob(ddir+'/output/powerspecs/powerspec_*.txt', recursive=True)])[-1]
    max_rows = int(np.loadtxt(f, max_rows=2)[1])
    f_mf_cb = np.loadtxt(f, skiprows=5, max_rows=max_rows)
    return f_mf_cb[:,0], f_mf_cb[:,1]

def HY_cb(ddir):
    PS = PowerSpectra(ddir)
    fo = np.sort([f for f in glob.iglob(ddir+'/output_c*', recursive=True)])[-1]
    snap = int(np.sort([f for f in glob.iglob(fo+'/powerspecs/powerspec_type1_*.txt', recursive=True)])[-1].split(".")[0][-1])
    max_rows = int(np.loadtxt(fo+'/powerspecs/powerspec_type1_%.3d.txt'%snap, max_rows=2)[1])
    f_hy_cb = np.loadtxt(fo+'/powerspecs/powerspec_type1_%.3d.txt'%snap,skiprows=5,max_rows=max_rows)
    return f_hy_cb[:,0], f_hy_cb[:,1]

def MF_m(ddir, omega_cb, omega_hdm, Ntau=20):
    PS = PowerSpectra(ddir)
    f_mf_cb = np.sort([f for f in glob.iglob(ddir+'/output/powerspecs/powerspec_*.txt', recursive=True)])[-1]
    max_rows = int(np.loadtxt(f_mf_cb, max_rows=2)[1])
    mf_cb = np.loadtxt(f_mf_cb,skiprows=5,max_rows=max_rows)
    snap_mf = int(np.sort([f for f in glob.iglob(ddir+'/output/powerspecs/powerspec_*.txt', recursive=True)])[-1].split(".")[0][-1])
    k = mf_cb[:,0]
    mf_De_cb = mf_cb[:,1]
    mf_De_hdm = np.interp(k, PS.k_mf, np.mean([PS.mf_De_flow(snap_mf, l) for l in range(Ntau)], axis=0))
    de = np.average(np.array([mf_De_cb, mf_De_hdm]), axis=0, weights=[omega_cb, omega_hdm])
    return k, de

def HY_m(ddir, omega_cb, omega_hdm, fu, Ntau=20):
    PS = PowerSpectra(ddir)
    k = PS.k_c0
    fo = np.sort([f for f in glob.iglob(ddir+'/output_c*', recursive=True)])[-1]
    fon = int(fo[-1])
    snap = int(np.sort([f for f in glob.iglob(fo+'/powerspecs/powerspec_type1_*.txt', recursive=True)])[-1].split(".")[0][-1])
    snap_mf = int(np.sort([f for f in glob.iglob(ddir+'/output/powerspecs/powerspec_*.txt', recursive=True)])[-1].split(".")[0][-1])
    max_rows = int(np.loadtxt(fo+'/powerspecs/powerspec_type1_%.3d.txt'%snap, max_rows=2)[1])
    hy_De_cb = np.loadtxt(fo+'/powerspecs/powerspec_type1_%.3d.txt'%snap,skiprows=5,max_rows=max_rows)[:,1]
    hy_De_hdm_p = np.mean([PS.ptype_De(fon, snap, i+2) for i in range(fon)], axis=0)[1]
    if fu == Ntau:
        hy_De_m     = np.average(np.array([hy_De_cb, hy_De_hdm_p]), axis=0, weights=[omega_cb, omega_hdm])
    else:
        hy_De_hdm_l = np.interp(k, PS.k_mf, np.mean([PS.mf_De_flow(snap_mf, l) for l in range(fu, Ntau)], axis=0))
        hy_De_m     = np.average(np.array([hy_De_cb, hy_De_hdm_p, hy_De_hdm_l]), axis=0, weights=[omega_cb, omega_hdm*fu/Ntau, omega_hdm*(1-fu/Ntau)])
    return k, hy_De_m

### Gauss-Laguerre functions

def Pcb(ddir, file):
    boxsize  = int(np.loadtxt(ddir+'/powerspecs/powerspec_%.3d.txt'%file,max_rows=3)[2])
    max_rows = int(np.loadtxt(ddir+'/powerspecs/powerspec_%.3d.txt'%file,max_rows=2)[1])
    f = np.loadtxt(ddir+'/powerspecs/powerspec_%.3d.txt'%file,skiprows=5, max_rows=max_rows)
    return f[:,0], f[:,1], f[:,2]*boxsize**3

def Pcb_hy(ddir, file):
    boxsize  = int(np.loadtxt(ddir+'/powerspecs/powerspec_type1_%.3d.txt'%file,max_rows=3)[2])
    max_rows = int(np.loadtxt(ddir+'/powerspecs/powerspec_type1_%.3d.txt'%file,max_rows=2)[1])
    f = np.loadtxt(ddir+'/powerspecs/powerspec_type1_%.3d.txt'%file,skiprows=5, max_rows=max_rows)
    return f[:,0], f[:,1], f[:,2]*boxsize**3

def Pnu_mflr(ddir, file, omega):

    # nu part
    ncols = len(open(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file).readline().split(',')) - 1
    k_mf = np.loadtxt(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file, delimiter=',',usecols=range(1))
    delta_nu_mf = np.average(np.loadtxt(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file, delimiter=',',usecols=range(1,ncols)),axis=1,weights=omega)
    max_rows = int(np.loadtxt(ddir+'/powerspecs/powerspec_%.3d.txt'%file,max_rows=2)[1])
    k = np.loadtxt(ddir+'/powerspecs/powerspec_%.3d.txt'%file,skiprows=5, max_rows=max_rows)[:,0]
    delta_nu_mfi = np.interp(k, k_mf, delta_nu_mf)
    Pnu_mf = delta_nu_mfi**2
    return k, Pnu_mf, Pnu_mf*k**3/(2*np.pi**2)

def Pm_mflr(ddir, file, omega, omega_cb):

    # nu part
    ncols = len(open(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file).readline().split(',')) - 1
    k_mf = np.loadtxt(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file, delimiter=',',usecols=range(1))
    delta_nu_mf = np.average(np.loadtxt(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file, delimiter=',',usecols=range(1,ncols)),axis=1,weights=omega)
    max_rows = int(np.loadtxt(ddir+'/powerspecs/powerspec_%.3d.txt'%file,max_rows=2)[1])
    boxsize  = int(np.loadtxt(ddir+'/powerspecs/powerspec_%.3d.txt'%file,max_rows=3)[2])
    k = np.loadtxt(ddir+'/powerspecs/powerspec_%.3d.txt'%file,skiprows=5, max_rows=max_rows)[:,0]
    delta_nu_mfi = np.interp(k, k_mf, delta_nu_mf)

    # cdm part
    f_cb = np.loadtxt(ddir+'/powerspecs/powerspec_%.3d.txt'%file,skiprows=5, max_rows=max_rows)
    delta_cb = np.sqrt(f_cb[:,2])*boxsize**(3/2)

    # tot
    omega_tot = omega_cb+np.sum(omega)
    f_cb = omega_cb/omega_tot
    f_nu = np.sum(omega)/omega_tot
    delta_tot = f_cb*delta_cb + f_nu*delta_nu_mfi
    Pm_mflr = delta_tot**2
    return k, Pm_mflr, Pm_mflr*k**3/(2*np.pi**2)#/boxsize**3


def Pnu_hyb(ddir, file, nconv, omega):

    # nu part
    ncols = len(open(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file).readline().split(',')) - 1
    k_mf = np.loadtxt(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file, delimiter=',',usecols=range(1))
    delta_nu_mf = np.average(np.loadtxt(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file, delimiter=',',usecols=range(nconv+1,ncols)),axis=1,weights=omega[nconv:])
    max_rows = int(np.loadtxt(ddir+'/powerspecs/powerspec_type1_%.3d.txt'%file,max_rows=2)[1])
    boxsize = int(np.loadtxt(ddir+'/powerspecs/powerspec_type1_%.3d.txt'%file,max_rows=3)[2])
    nu_parts = [np.loadtxt(ddir+'/powerspecs/powerspec_type%d_%.3d.txt'%(j, file),skiprows=5, max_rows=max_rows) for j in range(2, nconv+2)]
    delta_nu_parts = [np.sqrt(nu_parts[k][:,2])*boxsize**(3/2) for k in range(nconv)]
    k = nu_parts[0][:,0]
    delta_nu_mfi = np.interp(k, k_mf, delta_nu_mf)
    omega_hy_tot = np.sum(omega[:nconv])
    omega_lr_tot = np.sum(omega)-omega_hy_tot
    om_list = [omega[i] for i in range(nconv)]
    om_list.append(omega_lr_tot)
    de_list = delta_nu_parts
    de_list.append(delta_nu_mfi)
    delta_nu = np.average(np.array(de_list), axis=0, weights=np.array(om_list))
    Pnu_hyb = delta_nu**2
    return k, Pnu_hyb, Pnu_hyb*k**3/(2*np.pi**2)

def Pm_gse(ddir, file, tau, omega, omega_cb, mass):
    max_rows = int(np.loadtxt(ddir+'/powerspecs/powerspec_%.3d.txt'%file,max_rows=2)[1])
    boxsize = int(np.loadtxt(ddir+'/powerspecs/powerspec_%.3d.txt'%file,max_rows=3)[2])

    f_cb = np.loadtxt(ddir+'/powerspecs/powerspec_%.3d.txt'%file,skiprows=5, max_rows=max_rows)
    k = f_cb[:,0]
    delta_cb = np.sqrt(f_cb[:,2])*boxsize**(3/2)
    delta_nu = np.zeros_like(delta_cb)

    omega_tot = omega_cb+np.sum(omega)
    sumN = 0
    for i in range(len(omega)):
        kfs = np.sqrt(1.5)*1e2/299792.0*mass/tau[i]*np.sqrt(omega_tot)
        w = omega[i]/omega_tot
        sumN += (kfs**2/(k**2+k*kfs+kfs**2)-1)*w

    f_cb = omega_cb/omega_tot
    delta_nu = sumN*delta_cb
    f_nu = np.sum(omega)/omega_tot
    delta_tot = (1+sumN)*delta_cb
    #delta_tot = f_cb*delta_cb+f_nu*delta_nu
    Pm = delta_tot**2
    return k, Pm, Pm*k**3/(2*np.pi**2)#/boxsize**3




def Pm_hyb(ddir, file, nconv, omega, omega_cb):

    # nu part
    ncols = len(open(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file).readline().split(',')) - 1
    k_mf = np.loadtxt(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file, delimiter=',',usecols=range(1))
    delta_nu_mf = np.average(np.loadtxt(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file, delimiter=',',usecols=range(nconv+1,ncols)),axis=1,weights=omega[nconv:])
    max_rows = int(np.loadtxt(ddir+'/powerspecs/powerspec_type1_%.3d.txt'%file,max_rows=2)[1])
    boxsize = int(np.loadtxt(ddir+'/powerspecs/powerspec_type1_%.3d.txt'%file,max_rows=3)[2])
    nu_parts = [np.loadtxt(ddir+'/powerspecs/powerspec_type%d_%.3d.txt'%(j, file),skiprows=5, max_rows=max_rows) for j in range(2, nconv+2)]
    delta_nu_parts = [np.sqrt(nu_parts[k][:,2])*boxsize**(3/2) for k in range(nconv)]
    k = nu_parts[0][:,0]
    delta_nu_mfi = np.interp(k, k_mf, delta_nu_mf)
    omega_hy_tot = np.sum(omega[:nconv])
    omega_lr_tot = np.sum(omega)-omega_hy_tot
    om_list = [omega[i] for i in range(nconv)]
    om_list.append(omega_lr_tot)
    de_list = delta_nu_parts
    de_list.append(delta_nu_mfi)
    delta_nu = np.average(np.array(de_list), axis=0, weights=np.array(om_list))

    # cdm part
    f_cb = np.loadtxt(ddir+'/powerspecs/powerspec_type1_%.3d.txt'%file,skiprows=5, max_rows=max_rows)
    delta_cb = np.sqrt(f_cb[:,2])*boxsize**(3/2)

    # tot
    omega_tot = omega_cb+np.sum(omega)
    f_cb = omega_cb/omega_tot
    f_nu = np.sum(omega)/omega_tot
    delta_tot = f_cb*delta_cb + f_nu*delta_nu
    Pm_mflr = delta_tot**2
    return k, Pm_mflr, Pm_mflr*k**3/(2*np.pi**2)#/boxsize**3

def Pm_hyb_tot(ddir, file, nconv, omega, omega_cb):

    # nu part
    ncols = len(open(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file).readline().split(',')) - 1
    k_mf = np.loadtxt(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file, delimiter=',',usecols=range(1))
    delta_nu_mf = np.average(np.loadtxt(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file, delimiter=',',usecols=range(nconv+1,ncols)),axis=1,weights=omega[nconv:])
    max_rows = int(np.loadtxt(ddir+'/powerspecs/powerspec_type1_%.3d.txt'%file,max_rows=2)[1])
    boxsize = int(np.loadtxt(ddir+'/powerspecs/powerspec_type1_%.3d.txt'%file,max_rows=3)[2])

    parts = np.loadtxt(ddir+'/powerspecs/powerspec_%.3d.txt'%(file),skiprows=5, max_rows=max_rows)
    delta_parts = np.sqrt(parts[:,2])*boxsize**(3/2)
    k = parts[:,0]
    delta_nu_mfi = np.interp(k, k_mf, delta_nu_mf)
    omega_tot = omega_cb+np.sum(omega)
    f_cb = omega_cb/omega_tot
    f_nu_hyb = np.sum(omega[:nconv])/omega_tot
    f_nu_lr = 1-f_cb-f_nu_hyb

    delta_tot = (f_cb+f_nu_hyb)*delta_parts+f_nu_lr*delta_nu_mfi

    Pm_mflr = delta_tot**2
    return k, Pm_mflr, Pm_mflr*k**3/(2*np.pi**2)

def Pnu_hyb_tot(ddir, file, nconv, omega):

    # nu part
    ncols = len(open(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file).readline().split(',')) - 1
    k_mf = np.loadtxt(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file, delimiter=',',usecols=range(1))
    delta_nu_mf = np.average(np.loadtxt(ddir+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%file, delimiter=',',usecols=range(nconv+1,ncols)),axis=1,weights=omega[nconv:])
    max_rows = int(np.loadtxt(ddir+'/powerspecs/powerspec_type1_%.3d.txt'%file,max_rows=2)[1])
    boxsize = int(np.loadtxt(ddir+'/powerspecs/powerspec_type1_%.3d.txt'%file,max_rows=3)[2])

    nu_parts = np.loadtxt(ddir+'/powerspecs/powerspec_hdm_%.3d.txt'%(file),skiprows=5, max_rows=max_rows)
    delta_nu_parts = np.sqrt(nu_parts[:,2])*boxsize**(3/2) 
    k = nu_parts[:,0]
    delta_nu_mfi = np.interp(k, k_mf, delta_nu_mf)
    omega_hy_tot = np.sum(omega[:nconv])
    omega_lr_tot = np.sum(omega)-omega_hy_tot
    delta_nu = np.average(np.array([delta_nu_parts, delta_nu_mfi]), axis=0, weights=np.array([omega_hy_tot, omega_lr_tot]))
    Pnu_hyb = delta_nu**2
    return k, Pnu_hyb, Pnu_hyb*k**3/(2*np.pi**2)



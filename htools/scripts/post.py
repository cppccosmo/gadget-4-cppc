import sys
import argparse
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import Pk_library as PKL
import MAS_library as MASL


simfolder = sys.argv[1]
fname     = sys.argv[2]
ptype     = int(sys.argv[3])



def PowerSpectrum(fname,ptype,axis=0,verb=True):
    '''
    Computes the power spectrum of the divergence of a 3D velocity field
    As input pass the snapshot name, without the .hdf5 extension, and the particle type
    Returns k,Pk,Nmodes
    '''
    snap    = str(simfolder)+'/'+str(fname)
    f       = h5.File(str(simfolder)+'/'+str(fname)+'.hdf5','r')
    N       = f['Header'].attrs['NumPart_Total'][ptype]
    n       = int(np.ceil(N**(1/3)))
    BoxSize = f['Header'].attrs['BoxSize']
    grid    = f['Parameters'].attrs['GridSize']
    MAS     = 'CIC'  
    axis    = 0
    threads = 16
    ptypes  = [ptype] 


    delta = MASL.density_field_gadget(snap, ptypes, grid, MAS, do_RSD=False, axis=axis, verbose=verb)
    Pk    = PKL.Pk(delta, BoxSize, axis, MAS, threads, verbose=verb)
    k       = Pk.k3D
    Pk0     = Pk.Pk[:,0] #monopole
    Pk2     = Pk.Pk[:,1] #quadrupole
    Pk4     = Pk.Pk[:,2] #hexadecapole
    Pkphase = Pk.Pkphase #power spectrum of the phases
    Nmodes  = Pk.Nmodes3D
    np.savetxt("postprocess/spectra/"+str(simfolder)+"_ps.txt",np.column_stack((k,Pk0,Pk2,Pk4,Pkphase,Nmodes)))
    #return k,Pk0,Pk2,Pk4,Pkphase,Nmodes
    return 


def PS_theta(fname,ptype):
    '''
    Computes the power spectrum of the divergence of a 3D velocity field
    As input pass the snapshot name, without the .hdf5 extension, and the particle type
    Returns k,Pk,Nmodes
    '''
    f       = h5.File(fname+'.hdf5','r')
    N       = f['Header'].attrs['NumPart_Total'][ptype]
    n       = int(np.ceil(N**(1/3)))
    BoxSize = f['Header'].attrs['BoxSize']
    grid    = f['Parameters'].attrs['GridSize']
    MAS     = 'CIC'  
    axis    = 0
    threads = 16

    Vx = np.array(f['PartType%d/Velocities'%ptype][:,0]).reshape((n,n,n)).astype('float32')
    Vy = np.array(f['PartType%d/Velocities'%ptype][:,1]).reshape((n,n,n)).astype('float32')
    Vz = np.array(f['PartType%d/Velocities'%ptype][:,2]).reshape((n,n,n)).astype('float32')

    k, Pk, Nmodes = PKL.Pk_theta(Vx,Vy,Vz,BoxSize,axis,MAS,threads)
    return k,Pk,Nmodes

def projection(simfol,snap_num,patype,axis,mtype):
    snapshot  = str(simfol)+'/snapshot_%3d'%snap_num
    f = h5.File(str(snapshot)+'.hdf5','r')
    z = f['Header'].attrs['Redshift']
    L = f['Header'].attrs['BoxSize']
    grid  = f['Parameters'].attrs['GridSize']
    BoxSize = L
    ptypes  = patype                   
    MAS       = 'CIC'               
    verbose   = True
    threads   = 4
    delta     = MASL.density_field_gadget(snapshot, ptypes, grid, MAS, do_RSD=False, axis=0, verbose=verbose)
    fig,ax = plt.subplots(1,1,figsize=(13,12))
    im = ax.imshow(np.log10(np.mean(delta,axis=axis)),cmap='cmr.%s'%mtype,origin='lower',extent=[-L/2,L/2,-L/2,L/2])
    plt.colorbar(im)
    ax.set_xlabel(r'$x~({\rm Mpc}/h)$')
    ax.set_ylabel(r'$y~({\rm Mpc}/h)$')
    print('Saving plot in %s'%(str(simfol)+'/plots'))
    if len(ptype) > 1:
        fig.savefig(str(simfol)+'/plots/'+'ptype_%d-%d_z%.1f.pdf'%(ptype[0],ptype[-1],z),bbox_inches='tight')
    else:
        fig.savefig(str(simfol)+'/plots/'+'ptype_%d_z%.1f.pdf'%(ptype[0],z),bbox_inches='tight')

    
#axis = 0
#cmap = 'amber'
#patype = [1]  
#projection(args.simfol,args.snap_num,patype,axis,cmap)

if __name__ == '__main__':
    PowerSpectrum(fname,ptype)
    #projection(simfol,snap_num)

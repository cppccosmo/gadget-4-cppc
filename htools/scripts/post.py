import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import Pk_library as PKL
import MAS_library as MASL


parser = argparse.ArgumentParser()

parser.add_argument('-pt',"--part_type",dest="ptype",type=int)

patype = [ptype]

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
    print("type:", patype)
    #projection(args.simfol,args.snap_num,patype,axis,cmap)

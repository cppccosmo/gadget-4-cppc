import numpy as np


'''
TO ORGANIZE
'''

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



def PS_momentum(fname,ptype):
    '''
    Computes the density constrast and momentum auto- and cross-power spectra
    As input pass the snapshot name, without the .hdf5 extension, and the particle type
    Returns k,Pk_dd,Pk_tt,Pk_dt,Nmodes, where:
            Pk_dd is the standard power spectrum in (Mpc/h)^3
            Pk_tt is the momentum auto-power spectrum in (km/s)^2*(Mpc/h)^3
            Pk_dt is the density-momentum cross-power spectrum in (km/s)*(Mpc/h)^3
    '''
    f       = h5.File(fname+'.hdf5','r')
    N       = f['Header'].attrs['NumPart_Total'][ptype]
    n       = int(np.ceil(N**(1/3)))
    BoxSize = f['Header'].attrs['BoxSize']
    grid    = f['Parameters'].attrs['GridSize']
    MAS     = 'CIC'  
    axis    = 0
    threads = 16

    delta = MASL.density_field_gadget(fname, [ptype], grid, MAS, do_RSD=False, axis=0, verbose=True)
    delta /= np.mean(delta, dtype=np.float32);  delta -= 1.0

    Vx = np.array(f['PartType%d/Velocities'%ptype][:,0]).reshape((n,n,n)).astype('float32')
    Vy = np.array(f['PartType%d/Velocities'%ptype][:,1]).reshape((n,n,n)).astype('float32')
    Vz = np.array(f['PartType%d/Velocities'%ptype][:,2]).reshape((n,n,n)).astype('float32')

    k, Pk_dd, Pk_tt, Pk_dt, Nmodes = PKL.XPk_dv(delta, Vx, Vy, Vz, BoxSize, axis, MAS, threads)
    return k,Pk_dd,Pk_tt,Pk_dt,Nmodes

    

def total_Pnu(streams,converted_list,converted_pk_list, Nstreams):
    '''
    Computes the neutrino delta power spectrum ...
    TO FINISH
    '''
    assert len(converted_list) == len(converted_pk_list)
    num_particle_conv = len(converted_list)
    num_streams_in_particle = np.ones(num_particle_conv)
    
    all_converted_streams = []
    all_remaining_streams = list(range(1,Nstreams+1))
    k_streams = streams[:,0]
    for i in range(num_particle_conv):
        #print(i)
        num_streams_in_particle *= len(converted_list[i])
        
        for j in range(len(converted_list[i])):
            all_converted_streams.append(converted_list[i][j])
            
    for stream in all_converted_streams:
        all_remaining_streams.remove(stream)
        
    stream_delta_sum = np.zeros(streams.shape[0])
    for stream in all_remaining_streams:
        stream_delta_sum += streams[:,stream]
    delta_converted_sum = np.zeros(converted_pk_list[0].k.shape)
    k_converted=converted_pk_list[0].k
    for i in range(num_particle_conv):
        Pk = converted_pk_list[i]
        assert np.all(Pk.k==k_converted)
        delta_converted_sum += np.sqrt(Pk.Pk)*num_streams_in_particle[i]
        
    #total_delta = (delta_converted_sum + np.interp(k_converted,k_streams,stream_delta_sum))/Nstreams
    total_delta = (stream_delta_sum + np.interp(k_streams,k_converted,delta_converted_sum))/Nstreams
    
    returnarray = np.zeros([total_delta.size,2])
    #returnarray[:,0] = k_converted
    returnarray[:,0] = k_streams
    returnarray[:,1] = total_delta**2
    
    return returnarray


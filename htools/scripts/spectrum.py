import numpy as np

def total_Pnu(streams,converted_list,converted_pk_list, Nstreams):
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



def velocity_PS():


import os, sys
import numpy as np

root      = sys.argv[1]
snap      = int(sys.argv[2])
Nstreams  = int(sys.argv[3])
in_stream = int(sys.argv[4])
conv      = int(sys.argv[5])

def load_powerspec(powerspecfile):
    header_length = 5
    with open(powerspecfile) as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]

    numlines = int(lines[1])
    boxlength = float(lines[2])

    datalines = lines[header_length:header_length+numlines]
    data = np.zeros([numlines,5])

    for i, line in enumerate(datalines):
        data[i,:] = np.fromstring(line, dtype=float, sep=' ')

    k = data[:,0]
    pk = data[:,2]*boxlength**3
    datarray = np.array([k,pk])

    return datarray.transpose()

def hyb_restart(root,snap,Nstreams,stream):
    '''
    '''
    proot = os.path.abspath(os.path.join(root, os.pardir))
    ddelta = np.loadtxt(root+'/neutrino_stream_data/neutrino_delta_stream_%.3d.csv'%snap,delimiter=',',usecols=range(0,stream*Nstreams+1))
    dtheta = np.loadtxt(root+'/neutrino_stream_data/neutrino_theta_stream_%.3d.csv'%snap,delimiter=',',usecols=range(0,stream*Nstreams+1))
    ps_cb = load_powerspec(root+'/powerspecs/powerspec_%.3d.txt'%snap)

    if Nstreams == 1:
        print('Converting stream %d'%stream)
        delta = np.column_stack([ddelta[:,0],ddelta[:,stream]**2])
        growth = np.column_stack([ddelta[:,0],dtheta[:,stream]/ddelta[:,stream]])
    else:
        print('Converting %d streams starting from stream %d'%(Nstreams,stream))
        delta_s = np.sum(ddelta[:,stream:stream+Nstreams],axis=1)/Nstreams
        theta_s = np.sum(dtheta[:,stream:stream+Nstreams],axis=1)/Nstreams

        delta = np.column_stack([ddelta[:,0],delta_s**2])
        growth = np.column_stack([ddelta[:,0],theta_s/delta_s])

    # Save param files 
    np.savetxt(proot+'/PSFile_c%.3d.txt'%conv,delta,fmt='%.6e')
    np.savetxt(proot+'/GRFile_c%.3d.txt'%conv,growth,fmt='%.6e')
    np.savetxt(proot+'/CBFile_c%.3d.txt'%conv,ps_cb,fmt='%.6e')

if __name__ == "__main__":
    hyb_restart(root,snap,Nstreams,in_stream)


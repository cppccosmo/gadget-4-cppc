import argparse
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument('-pk',"--out_pk",dest="out_pk",type=str,required=False)
parser.add_argument('-gr',"--out_gr",dest="out_gr",type=str,required=False)
parser.add_argument('-cb',"--out_cb",dest="out_cb",type=str,required=False)
parser.add_argument('-pp',"--pp_index",dest="pp_index",type=int,required=False)
parser.add_argument('-sa',"--snap_num",dest="snap_num",type=int,required=False)
parser.add_argument('-id',"--input_delta",dest="in_de",type=str,required=False)
parser.add_argument('-it',"--input_theta",dest="in_th",type=str,required=False)
parser.add_argument('-icb',"--input_cb",dest="in_cb",type=str,required=False)
parser.add_argument('-str',"--stream",dest="stream",type=int,required=False)
parser.add_argument('-Ns',"--Nstreams",dest="Nstreams",type=int,required=False)

args = parser.parse_args()

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

def restart_hybrid(simfol,out_pk,out_gr,out_cb,input_delta,input_theta,input_pcb,streamno,Nstreams):
    # Prepare PowerSpectrumFile
    delta_nu_dat = np.loadtxt(input_delta,delimiter=',',usecols=range(0,Nstreams+1))
    k_delta = delta_nu_dat[:,0]
    delta = delta_nu_dat[:,streamno]**2
    data_delta = np.array([k_delta,delta])
    np.savetxt(simfol+'/'+f'{out_pk}.txt',data_delta.transpose())
    # Prepare GrowthRateFile
    theta_nu_dat = np.loadtxt(input_theta,delimiter=',',usecols=range(0,Nstreams+1))
    growthrate = theta_nu_dat[:,streamno]/delta_nu_dat[:,streamno]
    data_growthrate = np.array([k_delta,growthrate])
    np.savetxt(simfol+'/'+f'{out_gr}.txt',data_growthrate.transpose())
    # Prepare CBPowerSpectrumFile
    pcbdat = load_powerspec(input_pcb)
    np.savetxt(simfol+'/'+f'{out_cb}.txt',pcbdat)

def restart_hybrid_multi(simfol,out_pk,out_gr,out_cb,input_delta,input_theta,input_pcb,streamno,Nstreams):
    # Prepare PowerSpectrumFile
    delta_nu_dat = np.loadtxt(input_delta,delimiter=',',usecols=range(0,51)) # To fix number
    k_delta = delta_nu_dat[:,0]
    delta = np.sum(delta_nu_dat[:,streamno:streamno+Nstreams],axis=1)/Nstreams
    delta_sq = delta**2
    data_delta = np.array([k_delta,delta_sq])
    np.savetxt(simfol+'/'+f'{out_pk}.txt',data_delta.transpose())
    # Prepare GrowthRateFile
    theta_nu_dat = np.loadtxt(input_theta,delimiter=',',usecols=range(0,51)) # To fix number
    theta = np.sum(theta_nu_dat[:,streamno:streamno+Nstreams],axis=1)/Nstreams
    growthrate = theta/delta
    data_growthrate = np.array([k_delta,growthrate])


if args.Nstreams > 1:
    restart_hybrid_multi(args.simfol,args.out_pk,args.out_gr,args.out_cb,args.in_de,args.in_th,args.in_cb,args.stream,args.Nstreams)
else:
    restart_hybrid(args.simfol,args.out_pk,args.out_gr,args.out_cb,args.in_de,args.in_th,args.in_cb,args.stream,args.Nstreams)


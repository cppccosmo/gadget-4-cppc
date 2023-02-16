import sys
import argparse
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument('-ic',"--in_class",dest="in_class",type=str)
parser.add_argument('-im',"--in_mf",dest="in_mf",type=str)
parser.add_argument('-op',"--out_pk",dest="out_pk",type=str)
parser.add_argument('-og',"--out_gr",dest="out_gr",type=str)

args = parser.parse_args()

def scale_back(fclass,fmflr,f_pk,f_f):
    k,Pk_class = np.loadtxt(f"{fclass}.dat",unpack=True)
    mf_data    = np.loadtxt(f"{fmflr}.txt",skiprows=4)
    Pk_i       = interp1d(k,Pk_class,kind='cubic')
    np.savetxt(str(f_pk)+'.txt',np.column_stack((mf_data[:,1],Pk_i(mf_data[:,1])*mf_data[:,2]**2)))
    np.savetxt(str(f_f)+'.txt', np.column_stack((mf_data[:,1],mf_data[:,3])))

def main():
    scale_back(args.in_class,args.in_mf,args.out_pk,args.out_gr)

if __name__ == "__main__":
    main()
















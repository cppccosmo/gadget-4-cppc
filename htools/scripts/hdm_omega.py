import sys
import numpy as np

switch        = int(sys.argv[1])
total_streams = int(sys.argv[2])
omega_tot     = float(sys.argv[3])
deg           = int(sys.argv[4])
curr_type     = int(sys.argv[5])
stream        = int(sys.argv[6])

def create_budget(o_tot, deg):
    omega_part = 0
    omega_lin = o_tot
    stream_omega = o_tot/total_streams
    opart_list = []
    olin_list = []
    for i in range(4):
        omega_part += stream_omega*deg
        omega_lin -= stream_omega*deg
        opart_list.append(omega_part)
        olin_list.append(omega_lin)
    np.savetxt('hdm_omega.txt',np.column_stack([olin_list, opart_list]), fmt='%.5f %.5f')

def read_budget(curr):
    olin, opart = np.loadtxt('hdm_omega.txt', unpack=True)
    c_olin = olin[curr-1]
    c_opart = opart[curr-1]
    print(c_olin, c_opart)

def read_omega(i):
    opart = np.loadtxt('omega_table.txt')
    c_opart = np.sum(opart[i])
    c_olin = omega_tot - cumpart
    print(c_olin, c_opart)

if __name__ == "__main__":
    if switch == 0:
        create_budget(omega_tot, deg)
    elif switch ==1:
        #read_budget(curr_type)
        read_omega(stream)




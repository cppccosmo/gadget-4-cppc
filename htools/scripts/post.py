import sys, glob
import postprocess as p
#import visualisation as v
import numpy as np

opt       = int(sys.argv[1])
sim       = sys.argv[2]

sim_dir = '/scratch/xu32/gp5547/hybrid/' + sim
pp_dir = sim_dir + '/postprocess'
plot_dir = sim_dir + '/postprocess/plots'

# Create power spectra
PS = p.PowerSpectra(sim_dir)
#
#k = PS.k_c0
#vlist = PS.vlist

def dump_spectra(grid=1024):
    output_mf = sim_dir+'/output'
    output_gse = sim_dir+'/output_gse'
    output_hy = np.sort([f for f in glob.iglob(sim_dir+'/output_c*', recursive=True)])[-1]
    Nconv = int(output_hy[-1])
    deg = len(np.loadtxt(sim_dir+'/vel_list_%d.txt'%Nconv))


    # MF 
    snap = int(np.sort([f for f in glob.iglob(output_mf+'/snapshot_*', recursive=True)])[-1].split(".")[0][-1])
    mf_x1, mf_y1, mf_z1 = PS.ptype_De(0, snap, 1)
    mf_x2, mf_y2, mf_z2 = PS.ptype_De(0, snap, 1, smooth=True, grid=grid)
    np.savetxt(pp_dir+"/ps_mf_type1_raw.txt", np.column_stack([mf_x1,mf_y1,mf_z1]))
    np.savetxt(pp_dir+"/ps_mf_type1.txt", np.column_stack([mf_x2,mf_y2,mf_z2]))

    # GSE
    snap = int(np.sort([f for f in glob.iglob(output_gse+'/snapshot_*', recursive=True)])[-1].split(".")[0][-1])
    g_x1, g_y1, g_z1 = PS.ptype_De(0, snap, 1, gse=True)
    g_x2, g_y2, g_z2 = PS.ptype_De(0, snap, 1, gse=True, smooth=True, grid=grid)
    np.savetxt(pp_dir+"/ps_gse_type1.txt",     np.column_stack([g_x1,g_y1,g_z1]))
    np.savetxt(pp_dir+"/ps_gse_type1_raw.txt", np.column_stack([g_x2,g_y2,g_z2]))

    # Hybrid 
    snap = int(np.sort([f for f in glob.iglob(output_hy+'/snapshot_*', recursive=True)])[-1].split(".")[0][-1])
    tot = []

    for i in range(Nconv+1):
        x1, y1, z1 = PS.ptype_De(Nconv, snap, i+1, smooth=True, grid=grid)
        x2, y2, z2 = PS.ptype_De(Nconv, snap, i+1)
        np.savetxt(pp_dir+"/ps_hy_type%d.txt"%(i+1), np.column_stack([x1,y1,z1]))
        np.savetxt(pp_dir+"/ps_hy_type%d_raw.txt"%(i+1), np.column_stack([x2,y2,z2]))

def dump_deltas():
    pass

def plot_something():
    pass

if __name__ == "__main__":
    if opt == 0:
        dump_spectra()
    #compute_deltas()
    #save_smoothed_spectra()




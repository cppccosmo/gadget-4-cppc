import postprocess as p
import numpy as np

sim       = sys.argv[1]
Nconv     = int(sys.argv[2])
deg       = int(sys.argv[3])
snap_list = list(sys.argv[4].split(" "))

sim_dir = '/scratch/xu32/gp5547/hybrid' + sim
spectra_dir = sim + '/postprocess/spectra'
plot_dir = sim + '/postprocess/plots'

# Create power spectra
PS = p.PowerSpectra(sim_dir)

k = PS.k_c0
vlist = PS.vlist

def compute_deltas():
    pass

def save_smoothed_spectra():
    for i in range(Nconv):
        conv = i+1; ptype = i+2
        snapshot = snpa_list[i]
        print("Computing spectrum_c%d_snap%.3d_type%d ... "%(conv, snapshot, ptype))
        sp = PS.ptype_De(conv, snapshot, ptype, smooth=True, grid=1024))
        np.save(spectra_dir+"/spectrum_c%d_snap%.3d_type%d"%(conv, snapshot, ptype),sp)

def plot_something():
    pass

if __name__ == "__main__":
    compute_deltas()
    save_smoothed_spectra()




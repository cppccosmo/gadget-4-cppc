import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-in',"--in_file",dest="ini_file",type=str)
parser.add_argument('-ou',"--out_file",dest="muflr_file",type=str)
parser.add_argument('-pf',"--print_flag",dest="pr_flag",type=int)
parser.add_argument('-zi',"--zi",dest="zi",type=float)
parser.add_argument('-zf',"--zf",dest="zf",type=float)

args = parser.parse_args()

def create_mflr_file(ini_file,muflr_file,pr_flag,zi,zf):
    d = {}
    with open(f"{ini_file}.ini") as f:
        for line in f:
            if "#" in line:
                continue
            if line.startswith('root'):
                break
            d[line.split()[0]] = line.split()[-1]

    h        = d['h']
    Omega_m  = d['Omega_m']
    Omega_b  = d['Omega_b']
    Omega_nu = d['Omega_ncdm']
    sigma8   = d['sigma8']
    n_s      = d['n_s']

    T_cmb_k  = 2.7255
    w0       = -1
    wa       = 0
    cs_1     = 0
    cs_2     = 0
    print_r  = pr_flag
    notu     = 0
    z_i      = zi
    nrs      = 1
    z_f      = zf
    cambf    = '/srv/scratch/cppcnbody/MuFLR/nu05_transfer_z0.dat'
    nu_appr  = 1

    fout = open(f'{muflr_file}.txt','w+')
    fout.write(str(n_s)+" # n_s: scalar spectral index \n")
    fout.write(str(sigma8)+"\n")
    fout.write(str(h)+"\n")
    fout.write(str(Omega_m)+"\n")
    fout.write(str(Omega_b)+"\n")
    fout.write(str(Omega_nu)+"\n")
    fout.write(str(T_cmb_k)+"\n")
    fout.write(str(w0)+"\n")
    fout.write(str(wa)+"\n")
    fout.write(str(cs_1)+"\n")
    fout.write(str(cs_2)+"\n")
    fout.write(str(print_r)+"\n")
    fout.write(str(notu)+"\n")
    fout.write(str(z_i)+"\n")
    fout.write(str(nrs)+"\n")
    fout.write(str(z_f)+"\n")
    fout.write(str(cambf)+"\n")
    fout.write(str(nu_appr)+"\n")
    fout.close()

def main():
    create_mflr_file(args.ini_file,args.muflr_file,args.pr_flag,args.zi,args.zf)

if __name__ == "__main__":
    main()


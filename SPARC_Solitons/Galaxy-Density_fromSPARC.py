#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.integrate import quad
import sys


# parameters 
GMsolar = 1.32e11 # G*Msolar in km^3/s^2
kpcinkm = 3.086e16 # 1 kpc in km
epsilon = 0.01 # parameter for derivative

mDM = 1.0 # times 10e-22 eV

rescalerho = 2.478e-24 / mDM / mDM # equal to (Msun/kpc^3)/(Mpl m/sqrt(4pi))^2, with m=10e-22 eV, used to rescale mass density
rescaler = 15615.2 / mDM# equal to kpc*m with m=10e-22 eV, used to rescale radius


# files
densfile_folder = 'BulgeDiskDec' # xxx.dens files should be in this folder
rotmodfile_folder = 'Rotmod' # xxx_rotmod.dat files should be in this folder
sparcfile_name = 'SPARC_Lelli2016c.mrt.txt'


########################################################################
# calculate density

# this function returns density
#
# r        : list of radius in kpc
# rhobulge : list of rho_bulge in Msolar/kpc^3 at each radius
# rhodisk  : list of rho_disk in Msolar/kpc^3 at each radius (z=0)
# zd       : rhodisk(z) = rhodisk(z=0) x exp(-|z|/zd)

def getrho(galaxyname):
    print("analyzing " + galaxyname + "...")

    # read Rd
    Rd = 0.
    ld = open(sparcfile_name)
    lines = ld.readlines()
    ld.close()
    flag = 0
    angle = 0.
    for line in lines:
        if line.find(galaxyname)>=0:
            items = line[:-1].split()
            Rd = float(items[9])
            angle = float(items[5])
            flag=1
    if(flag==0):
        print('no info : ' + galaxyname)
        return [[0.],[0.],[0.],1.]

    
    
    densfile_name = densfile_folder + '/' + galaxyname + '.dens'
    rotmodfile_name = rotmodfile_folder + '/' + galaxyname + '_rotmod.dat'
    
    # read .dens file
    rad_dens = []
    sbdisk_dens = []
    sbbul_dens = []
    for line in open(densfile_name,'r'):
        if(line[0] != '#'):
            items = line[:-1].split('\t')
            if(len(items)>=3):
                rad_dens.append(float(items[0]))
                sbdisk_dens.append(float(items[1]))
                sbbul_dens.append(float(items[2]))
    
    # read _rotmod.dat file
    rad_rotmod = []
    vgas_rotmod = []
    vdisk_rotmod = []
    vbul_rotmod = []
    sbdisk_rotmod = []
    sbbul_rotmod = []
    for line in open(rotmodfile_name,'r'):
        if(line[0] != '#'):
            items = line[:-1].split('\t')
            if(len(items)>=3):
                rad_rotmod.append(float(items[0]))
                vgas_rotmod.append(float(items[3]))
                vdisk_rotmod.append(float(items[4]))
                vbul_rotmod.append(float(items[5]))
                sbdisk_rotmod.append(float(items[6]))
                sbbul_rotmod.append(float(items[7]))
    """    
    # calculate r*v^2 in [kpc * km^2/s^2]
    rvv_rad = [0.]
    rvv_gas = [0.]
    rvv_disk = [0.]
    rvv_bul = [0.]
    for r,vgas,vdisk,vbul in zip(rad_rotmod, vgas_rotmod, vdisk_rotmod, vbul_rotmod):
        rvv_rad.append(r)
        rvv_gas.append(r*vgas*vgas)
        rvv_disk.append(r*vdisk*vdisk)
        rvv_bul.append(r*vbul*vbul)
    
    rvvfunc_gas = interp1d(rvv_rad, rvv_gas, kind='cubic')
    rvvfunc_disk = interp1d(rvv_rad, rvv_disk, kind='cubic')
    rvvfunc_bul = interp1d(rvv_rad, rvv_bul, kind='cubic')
    """
    
    
    ########################################################################
    # density of bulge
    """
    # bulge density from rotation curve
    rvvderiv_bul = lambda x : ( rvvfunc_bul(x+epsilon) - rvvfunc_bul(x) ) / epsilon
    rhofromv_bul = lambda rad : 1./GMsolar*kpcinkm/4./np.pi/rad/rad * rvvderiv_bul(rad)
    """    
    
    # density of bulge from photometry
    # solving sbbul(r) = \int dz \rho( \sqrt(r^2+z^2) )
    radnew = []
    rhonew = []
    for i in range(len(rad_dens)):
        r0 = rad_dens[-1-i]
        if(len(radnew)==0):
            rhonew.append(0.)
        else:
            r1 = rad_dens[-i]
            b = rhonew[-1] # b = rho(r1), a = rho(r0)
            c = (quad(lambda z: b * (r0 - np.sqrt(r0**2+z**2))/(r0-r1), 0., np.sqrt(r1**2 - r0**2) ))[0]
            k = (quad(lambda z: (np.sqrt(r0**2+z**2)-r1)/(r0-r1), 0., np.sqrt(r1**2 - r0**2) ))[0]
            if(len(radnew)>1):
                rhotmp = interp1d(radnew,rhonew,kind='linear')
                #c += (quad(lambda z: rhotmp(np.sqrt(r0**2+z**2)), np.sqrt(r1**2 - r0**2), np.sqrt(max(radnew)**2 - r0**2), points=radnew ))[0]
                c += (quad(lambda z: rhotmp(np.sqrt(r0**2+z**2)), np.sqrt(r1**2 - r0**2), np.sqrt(max(radnew)**2 - r0**2)))[0]
            rhonew.append( (sbbul_dens[-1-i]-c)/k )
        radnew.append(r0)
    
    rhonew_bul = interp1d(radnew,[1e6*x for x in rhonew], kind='linear')

    # check the consistency
    #for r,s in zip(rad_dens, sbbul_dens):
    #    print(r, s, quad(lambda z: 1e-6*rhonew_bul(np.sqrt(z**2+r**2)), 0., np.sqrt(max(rad_dens)**2-r**2)) )
    
    ########################################################################
    # density of disk

    # see p.6 in 1606.09251
    zd = 0.196 * Rd * 0.633
    
    rho_disk_dens = interp1d(rad_dens, [sb/2./zd for sb in sbdisk_dens], kind='cubic')
    rho_disk_rotmod = interp1d(rad_rotmod, [sb/2./zd for sb in sbdisk_rotmod], kind='cubic')
   
    # normalization between surface brightness in .dens file and _rotmod.dat file 
    rref = ( max(min(rad_dens),min(rad_rotmod)) + min(max(rad_dens),max(rad_rotmod)) )/2.
    factor = rho_disk_rotmod(rref) / rho_disk_dens(rref)
    #print('\t:', factor, angle, np.cos(angle*np.pi/180.), np.sin(angle*np.pi/180.))
    
    ########################################################################
    # density of gas
    # ?????
    """
    ########################################################################
    # plot for check
    rads = np.linspace(max(rvv_rad)*0.01, max(rvv_rad)*0.99, 100 )
    plt.plot(rads, [rhofromv_bul(r) for r in rads], label='Bulge density')
    plt.plot(rad_dens, [rhonew_bul(r) for r in rad_dens], marker='o', label='Bulge density')
    rads2 = np.linspace(min(rad_rotmod), max(rad_rotmod), 100)
    plt.plot(rads2, [rho_disk_rotmod(r) for r in rads2], label='Disk density at z=0')
    #rads3 = np.linspace(min(rad_dens), max(rad_dens), 100)
    plt.plot(rad_dens, [factor*rho_disk_dens(r) for r in rad_dens], marker='o', label='Disk density at z=0')
    plt.legend()
    plt.xlabel('radius [kpc]')
    plt.ylabel('density [$M_\odot$/kpc$^3$]')
    plt.yscale('log')
    plt.xscale('log')
    plt.savefig(galaxyname + '_check.pdf')
    plt.clf()
    """

    return [rad_dens, [rhonew_bul(r) for r in rad_dens], [1e6*factor*s for s in sbdisk_dens], zd]

########################################################################


# list of galaxies in SPARC data
names=["CamB", "D512-2", "D564-8", "D631-7", "DDO064", "DDO154", "DDO161", "DDO168", "DDO170", "ESO079-G014", "ESO116-G012", "ESO444-G084", "ESO563-G021", "F561-1", "F563-1", "F563-V1", "F563-V2", "F565-V2", "F567-2", "F568-1", "F568-3", "F568-V1", "F571-8", "F571-V1", "F574-1", "F574-2", "F579-V1", "F583-1", "F583-4", "IC2574", "IC4202", "KK98-251", "NGC0024", "NGC0055", "NGC0100", "NGC0247", "NGC0289", "NGC0300", "NGC0801", "NGC0891", "NGC1003", "NGC1090", "NGC1705", "NGC2366", "NGC2403", "NGC2683", "NGC2841", "NGC2903", "NGC2915", "NGC2955", "NGC2976", "NGC2998", "NGC3109", "NGC3198", "NGC3521", "NGC3726", "NGC3741", "NGC3769", "NGC3877", "NGC3893", "NGC3917", "NGC3949", "NGC3953", "NGC3972", "NGC3992", "NGC4010", "NGC4013", "NGC4051", "NGC4068", "NGC4085", "NGC4088", "NGC4100", "NGC4138", "NGC4157", "NGC4183", "NGC4214", "NGC4217", "NGC4389", "NGC4559", "NGC5005", "NGC5033", "NGC5055", "NGC5371", "NGC5585", "NGC5907", "NGC5985", "NGC6015", "NGC6195", "NGC6503", "NGC6674", "NGC6789", "NGC6946", "NGC7331", "NGC7793", "NGC7814", "PGC51017", "UGC00128", "UGC00191", "UGC00634", "UGC00731", "UGC00891", "UGC01230", "UGC01281", "UGC02023", "UGC02259", "UGC02455", "UGC02487", "UGC02885", "UGC02916", "UGC02953", "UGC03205", "UGC03546", "UGC03580", "UGC04278", "UGC04305", "UGC04325", "UGC04483", "UGC04499", "UGC05005", "UGC05253", "UGC05414", "UGC05716", "UGC05721", "UGC05750", "UGC05764", "UGC05829", "UGC05918", "UGC05986", "UGC05999", "UGC06399", "UGC06446", "UGC06614", "UGC06628", "UGC06667", "UGC06786", "UGC06787", "UGC06818", "UGC06917", "UGC06923", "UGC06930", "UGC06973", "UGC06983", "UGC07089", "UGC07125", "UGC07151", "UGC07232", "UGC07261", "UGC07323", "UGC07399", "UGC07524", "UGC07559", "UGC07577", "UGC07603", "UGC07608", "UGC07690", "UGC07866", "UGC08286", "UGC08490", "UGC08550", "UGC08699", "UGC08837", "UGC09037", "UGC09133", "UGC09992", "UGC10310", "UGC11455", "UGC11557", "UGC11820", "UGC11914", "UGC12506", "UGC12632", "UGC12732", "UGCA281", "UGCA442", "UGCA444"]

for name in ['UGC02953']:#'UGC01281' 'UGC01281', UGC02953                      'NGC1705','UGC02953','UGC05253','UGC06786']: #'CamB','NGC0891'...... 'UGC01281','UGC07524','NGC1705'
#for name in names:
    # r        : list of radius in kpc
    # rhobulge : list of rho_bulge in Msolar/kpc^3 at each radius
    # rhodisk  : list of rho_disk in Msolar/kpc^3 at each radius (z=0)
    # zd       : rhodisk(z) = rhodisk(z=0) x exp(-|z|/zd)
    r,rhobulge,rhodisk,zd = getrho(name)
    
    rdim = [item*rescaler for item in r]
    rhobulge_dim = [item*rescalerho for item in rhobulge]
    rhodisk_dim = [item*rescalerho for item in rhodisk]
    zd_dim = zd*rescaler
    
    f = open('Fits/' + name + '_rho.txt','w') # + name + '/' 
    for i in range(len(rhobulge_dim)):
        f.write("%s %s %s %s\n" % (rdim[i] , rhobulge_dim[i] , rhodisk_dim[i] , zd_dim))
    f.close()
    
    #r_1D = np.arange(r[0],r[-1],1e-2)
    #z_1D = np.arange(r[0],r[-1],1e-2)

#rhobulge_1D = interp1d(rdim,[x for x in rhobulge_dim], kind='linear', fill_value='extrapolate') # function of r
#rhodisk_1D = interp1d(rdim,[x for x in rhodisk_dim], kind='linear', fill_value='extrapolate') # function of \rho
    
    #for i in range(len(rdim)):
    #    for j in range(len(rdim)):
    #        tottable[i][j] = rhobulge_1D(np.sqrt(rdim[j]**2 + rdim[i]**2)) + rhodisk_1D(rdim[i])*np.exp(-abs(rdim[j])/zd)
    #[rhobulge_1D(np.sqrt(z**2 + r**2)) + rhodisk_1D(r)*np.exp(-abs(z)/zd) for r,z in itertools.product(rdim,rdim)]
    #print(tottable[0][0])
    #rhototal = interp2d(rdim,rdim,tottable,kind='linear')
    #func = rhototal(r_1D,z_1D)
    #rhobulge_1D(np.sqrt(z**2 + r**2)) + rhodisk_dim(r)*np.exp(-abs(z)/zd)

#plt.plot(rdim, [rhototal(x,0) for x in rdim], marker='o', label='z=0')
#    plt.plot(rdim, [rhototal(x,50.0) for x in rdim], marker='o', label='z=1')
#    #plt.plot(r, rhodisk, marker='o', label='Disk density (z=0)')
#    #plt.plot(r, [r*np.exp(-1./zd) for r in rhodisk], marker='o', label='Disk density (z=1kpc)')
#    plt.legend()
#    plt.xlabel('radius [1/m]')
#    plt.ylabel('density [$M_P^2 m^2 / 4\pi$]')
    #    plt.yscale('log')
    #    plt.xscale('log')
    #plt.title(name)
    #plt.savefig(name + '_density_JE.pdf')
    #plt.clf()


    # make a sample plot
    plt.plot(r, rhobulge, marker='o', label='Bulge density')
    plt.plot(r, rhodisk, marker='o', label='Disk density (z=0)')
    plt.plot(r, [r*np.exp(-1./zd) for r in rhodisk], marker='o', label='Disk density (z=1kpc)')
    plt.legend()
    plt.xlabel('radius [kpc]')
    plt.ylabel('density [$L_\odot$/kpc$^3$]')
    plt.yscale('log')
    plt.xscale('log')
    plt.title(name)
    plt.savefig(name + '_density.pdf')
    plt.clf()


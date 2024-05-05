The Spitzer Photometry and Accurate Rotation Curves (SPARC) database contains luminosity data and mass models for 175 galaxies of spiral and irregular type; the data can be found [here](http://astroweb.cwru.edu/SPARC/)

This project contains two tools:
1. **Galaxy-Density_fromSPARC.py** : Reads data for a given galaxy from SPARC and reconstructs its baryonic density function, both for the (spherically-symmetric) bulge and (cylindrically-symmetric) disk components. 

*Inputs:* To run this tool, one must download three datasets from SPARC (link above):
   * 'Galaxy Sample': Save to text file titled `SPARC_Lelli2016c.mrt.txt`
   * 'Bulge-Disk Decompositions': Unzip and save to directory titled `BulgeDiskDec`
   * 'Newtonian Mass Models': Unzip and save to directory titled `Rotmod`

One must also pass a galaxy name, e.g. `UGC02953`, as an argument when calling the Python script.

*Outputs:* Returns a figure illustrating the bulge and disk densities as a function of radial distance from the center of the galaxy, e.g. `UGC02953_density.pdf`

2. **ScrodingerPoissonSolver_Disk.cc** : Employs a successive-overrelaxation (SOR) method to derive the ground-state solution of the two-dimensional Schrodinger-Poisson system, in the presence of a (cylindrically-symmetric) external potential, e.g. arising from the presence of some disc galaxy in the SPARC sample.
 
*Inputs:* To run this tool, one must read in a density profile in cylindrical `r` and `z` coordinates, $\rho(r,z)$, in a .txt file organized as `[density    r    z]`. The units are assumed to be solar masses and kpc.
   
*Outputs:* A file containing the density profile of the cylidrical soliton / ground-state solution, titled `solution.dat`

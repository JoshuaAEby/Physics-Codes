Inelastic dark matter (DM) is a class of models in which two dark-sector particles which are narrowly-split in mass, which we call the ground and excited states. If they can interact with nuclei, such particles can upscatter (i.e. excite from the lighter ground state to the heavier excited state) when they interact with heavy elements in the bulk of the Earth. These excited states are metastable, and after traveling some distance they decay, often emitting a photon which can be searched for in experiment. This work determines the rate of DM upscatter and subsequent decay, and outputs the rate of signal photons as a function of time from this series of processes.

The DM-nuclear interaction chosen here is that of a magnetic transition dipole operator, defined by its coupling `gM`. For details, see [our paper](https://arxiv.org/abs/2312.08478)

This project contains four files:
* `Constants.c`: contains global constants used in the execution.
* `isodata.txt`: contains information on stable isotopes in the bulk of Earth, including abundance (in core, mantle, and crust), spin, and isotopic fraction.
* `Functions.h`: header file which contains functions, most notably
  * `sixDim_Integrand`, which outputs the signal rate of photons in the experiment for a given incoming DM velocity and scatter site in Earth;
  * `MC_Ints`, which runs a Monte Carlo sampling of velocites and scatter sites to approximate the 6-dimensional integral to give the total signal rate.
* `iDM.c`, which is the executable file written in C.

**Inputs:** The execution of this script in the terminal should be done through a command like so:
* `gcc -Wall -I/sw/include/ -c iDM.c; gcc -L/sw/lib iDM.o -lgsl -lm -lgslcblas; mv a.out iDM.out`, which produces an executable file `iDM.out`
* `./iDM.out <det> <elt> <A> <gM> <mDM> <delta> <lifemult> <Month> <logprec>`, where
  * `<det>` is one of the detector options defined in `Constants.c`, labeled by their first letter (e.g. `B` for Borexino at Gran Sasso, `S` for SUPL at Stawell Laboratory, etc.). One can add additional options to the file if needed.
  * `<elt>` is the element being targeted (e.g. `Al` for Aluminum, `Fe` for Iron, etc.).
  * `<A>` is the mass number of the isotope being targeted (e.g. `13` for 13Al, `56` for 56Fe, etc.).
    * The combination of `<elt>` and `<A>` determines which row of `isodata.txt` is read to the program, which therefore will terminate if the value is not in the table.
  * `<gM>` is the dimensionless coupling of the magnetic transition dipole operator (generally $\ll 1$).
  * `<mDM>` is the mass of the ground-state DM particle in GeV. (Benchmark value is 1000.)
  * `<delta>` is the mass splitting between the two DM states in keV (generally between $1-1000\,{\rm keV}$).
  * `<lifemult>` is a multiplier on the lifetime of the excited DM state. Normally the lifetime is derived from the operator under investigation, but this options lets one multiply this result by a constant equal to `lifemult` value. (Generally set to unity.)
  * `<Month>` is the month of the year (between 1 and 12, corresponding to January through December) being investigated; the signal is highly oscillatory and can vary significantly over the course of a day / week / month.
  * `<logprec>` is the log of the number of Monte Carlo calls, which largely determines the precision of the calculation and also the runtime. As an example, on my Macbook (Apple M2, 8GB memory, on a single core), a run with `logprec=7.4` will take several hours. Generally a higher `logprec` will reduce noise in the output, and too small a value will result in systematic error in the output. I have found a rule of thumb that seems to work, which is that for `logprec > 6.8` (roughly) the overall scale of the signal output converges, with noise that diminishes at larger values.
 
**Outputs:** The calculation will output a file with title `MiDM_det=<det>_<A><elt>_mDM=<mDM>_delta=<delta>_tMult=<lifemult>_TIME=2000-<month>-1_prec<logprec>_gM=<gM>`, e.g. `MiDM_det=B_29Si_mDM=1000_delta=15_tMult=1.0_TIME=2000-1-1_prec6.0_gM=0.141`

# Physics-Codes

This repository contains three projects:
* **AxionStars_and_Gravitational Atoms** analyzes the structure of and constraints on gravitationally-bound states of ultralight dark matter (also called *axions*). The two possible states, *axion stars* and *gravitational atoms*, are defined by the source of the gravitational force which holds the object together (self-gravity and the gravity of an external astrophysical object, respectively). It contains:
  * a Mathematica file containing the basic structure functions (e.g. mass, radius, and density) and current observational constraints on the density of gravitational atoms inside the Solar System.
* **Inelastic-Dark-Matter** allows one to determine the rate of signal photons in a two-state dark sector which upscatters against elements in the Earth and decays to a photon in a detector volume. It contains:
  * a C-based Monte-Carlo integration routine using gcc, and
  * auxiliary files required to run the integrator.
* **SPARC_Solitons** is used to analyze the density profile of disk galaxies and use their density profile to deduce the ground-state configuration from the Schrodinger-Poisson equations. It contains:
  * a Python script for extracting a two-dimensional (cylindrically-symmetric) density profile for a galaxy in the SPARC sample of galaxies using photometric data, and
  * a C++-based integrator using a successive over-relaxation method to simultaneously solve the Schrodinger-Poisson equations in cylindrical symmetry.

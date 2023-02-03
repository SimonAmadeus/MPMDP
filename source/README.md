# Multi purpose molecular dynamics package (MPMDP).

## Features:

- Harmonic bonds
- Harmonic angles
- Lennar-Jones pair interactions
- Hybrid particle-field molecular dynamics 
-> with classical interactions and Gaussian convolution of atomic denisities with FFT
- Shear flow (Non-eq.)


### In progress:

- Different particle types
- DPD
- Lowe-Andersen thermostat
- Multi particle collision dynamics

## Tools

- RDF
- MSD
- E2E autocorrelation & relaxation time

## Running the code:

### Julia packages to install before running:
import Pkg; Pkg.add("Formatting"); Pkg.add("Distributions"); Pkg.add("ProgressMeter"); Pkg.add("Debugger"); Pkg.add("Statistics"); Pkg.add("LinearAlgebra"); Pkg.add("Random"); Pkg.add("FFTW")

### Runing the code.
You will need an input file, typically called input.data and a control file, typically named run.jl. Run with the folloqing command:

julia run.jl

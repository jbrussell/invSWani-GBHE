# invSWani-GBHE
Invert Rayleigh and Love wave phase velocities for azimuthal anisotropy parameters G, B, H, and E using MINEOS to calculate kernels.

## Contents
- ./FORTRAN : contains all Fortran binaries required to build MINEOS and idagrn6
- ./run_MINEOS : contains the MATLAB wrappers for running MINEOS to build the mode tables and invert for depth dependent azimuthal anisotropy.

## Getting Started

Must have installed
- gfortran (other Fortran compilers might work but have not been tested)
- MATLAB

### For compiling fortran binaries, see [./FORTRAN/README](https://github.com/jbrussell/MINEOS_synthetics/blob/master/FORTRAN/README)


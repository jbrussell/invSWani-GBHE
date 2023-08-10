Steps for running inversion:

Modify parameter_FRECHET.m
- set param.DATAmat_path to point to your dataset
- set param.CARDID to your desired MINEOS card
- set SONLY=1 or TONLY=1 to calculate spheroidal or toroidal modes, respectively (one of them must be set to 0)

a1_run_mineos_check.m -> calculate mode sets
- Run once with SONLY = 1 in parameter_FRECHET to generate spheroidal modes
- Run once with TONLY = 1 in parameter_FRECHET to generate toroidal modes

a2_mk_ACFLNkernels_love_all.m -> loop through all Rayleigh and Love data and generate anisotropy kernels

b1_plot_kernels_scaled_S0S1T0_GBHE.m -> plot the kernels with their proper scalings

b2_InvGBHE_fxnize_bootstrap -> run bootstrap inversion of Rayleigh and Love anisotropy strength and direction for depth-dependent parameters G, B, H, and E strength and direction. The bootstrap procedure randomly perturbs each datapoint within its 1-sigma uncertainty bounds and is repeated par.nbs times.

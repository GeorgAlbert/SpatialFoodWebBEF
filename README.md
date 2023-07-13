# SpatialFoodWebBEF
Code used to generate and simulate spatial food webs used in Albert et al. (2023) Animal and plant space-use drive plant diversity-productivity relationships. Ecology Letters

Original simulations were done in Julia version 1.6.1. Code was last tested in Julia version 1.8. 

Especially in more complex setups (i.e. nested food webs), simulations can take very long. Therfore, output is saved regularly, allowing that simulations can be interrupted and started again later. Code was designed to run on a cluster with multiple simulations running in parallel.

The general workflow:
1) If necessary, install required packages using 0_package_install.jl.
2) Create parameter sets for simulations using 2_create_parameters.jl.
3) Simulate food web dynamics using 3_ODEsimulate.jl. This can be done directly, specifying for which parameter set simulations should be done when defining args (line 48). An example SLURM-script is provided (submit.sh) when using a cluster.
4) Extract data from saved simulation output using 4_extract.jl.
5) Data is processed in R using datahandling.R. Figures need to be finalized using other programs (e.g. inkscape, gimp, powerpoint, ...)

For further details, see comments in code, check out the paper, or contact me. 

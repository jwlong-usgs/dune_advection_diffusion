# dune_advection_diffusion
Development codes for modeling storm-induced dune evolution

DUNE DIFFUSION README

- This set of codes is being established in order to develop codes capable of simulating dune evolution during storm events.

- The goal of the model is to efficiently compute the morphological change during the regimes identified by Sallenger 2000 (collision, overwash, inundation).

- Initially this is a framework that relies on different models to simulated the different regimes with a wrapper code that checks conditions at each time step and calls the appropriate subroutines.

Current codes include:

 - wrapper.m:  This is the framework code that loads the dune information and loops through each time step to call the required model components.
 - LEH04MainProgram.m:  Original code to simulate dune erosion provided by Meg Palmsten.  Code is based on the LEH04 model with updates published in Palmsten and Holman 2011
 - LEH04MainProgram_v2.m:  Updated version of the code that takes out figure plotting and other loops so it can be used in the overall model framework.
 - LEH04_notime.m:  Updated version of the original code that takes out time loop so it can be used in the overall model framework.
 - find_dlow_dhigh.m:  Identifies the dune toe and dune crest.  Within the framework this is called each time step and compared with the wave levels to identify the current regime.
 - extreme.m:  Function called by find_dlow_dhigh.m to help find the dune crest.
 - callLEH.m:  Testing code that calls LEH04_notime.m
 - dune_diffusion.m:  Initial code of the diffusion equation to diffuse the dune face during dune erosion.
 - dune_advection_diffusion.m:  Initial code of the advection-diffusion equation to move the dune landward and diffuse it during overwash and inundation.

** many parameters and coding efficiencies still need tested and implemented.

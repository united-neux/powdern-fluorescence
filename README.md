# powdern-flourescence
The objectives of this project are:

1. [ ] add a powder-scattering option to the flourescence-sample component based on xraylib
2. [ ] clean up code - with the aim of creating a generalized library that can be called for the sample model.
  This should be based on xraylib as well
  1. [ ] identify compounds (or just define a flag for that)
  2. [ ] create fucntion for retrieving cross sections (mapping to xraylib)
  3. [ ] functions for choosing between the set of enabled processes.

The general flow should be:
1. intersect the sample object, and move to it if necessary
2. get the overall cross-section for all processes
3. two branches: tunnel or "scatter" - pick one.
4. If tunnel:
  1. move to end and exit component
5. If scatter: 
  1. pick one and get new direction/photon - return to 1.

- Scattering processes should be created as a static list (this will allow for GPU-mode), with a fixed set of outcomes.
- Each process will then have an enabled flag that determines if it is possible to reach it or no.

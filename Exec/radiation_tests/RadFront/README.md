# RadFront

This setup is based on the Optically-Thin Streaming of Radiation front
described on the Castro paper II but with some differences:

1. Instead of defining the B.C. by Er, I'm using the corresponding value
   for Fr, assuming that for optically thin Fr=-cEr

2. There is no filled region with radiation

3. I'm using do_hydro=1 instead of having the hyperbolic update off. The medium
   has a very low density and temperature

4. I'm using the same value for Rosseland and Planck coeff.

The 1d version, has a domain size of 100 cm, whereas the 2d has both x and y of ~10^10 cm.
The opacities are changed to have similar optical depths along the domain for both setups.
In both cases it is seen that the radiation front propagates at speed of light. 


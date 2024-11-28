Python and fortran code to reproduce Philips 1956 model (https://doi.org/10.1002/qj.49708235202). Written by Martin Dix (@MartinDix)

- `msq_rand.f90` is an implementation the same random number generator that Phillips used, but has a short period. Allows the use of the same RNG in both fortran and python, and the original implementation in Philips (1956)
- `phillips_model.py` is the python script for producing an animation of a run
- `phillips.f90` is the fortran version of the same


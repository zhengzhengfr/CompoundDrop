# PNAS-Compound-drop-code
Numerical codes for compound drop impact

This is the numerical codes (based on the Basilisk open-source software, basilisk.fr) used to produce the numerical results presented in the PNAS paper.
Note that the Basilisk should be installed before running these codes.
Basically, to run in parallel, type these line one by one (27 means number of cores used):

qcc -source -D_MPI=1 0impactCompoundSolid-water.c 
mpicc -O2 -Wall -std=c99 -D_MPI=1 _0impactCompoundSolid-water.c -o cs -lm
mpirun -np 27 ./cs

And to run in normal conditions:
qcc -O2 -Wall 0impactCompoundSolid-water.c -o cs -lm -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
./cs
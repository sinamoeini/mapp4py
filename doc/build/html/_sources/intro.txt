
**************
 Introduction
**************

MAPP is a platform written in C++ for atomistic simulation of metals and alloys. It is equipped with most of popular forcefields designed for metals. In addition to being able to perform energy minimization and molecular dynamics (MD), it can also do diffusive molecular dynamics (DMD) simulations[?]. Since in addition to displacive mechanisms DMD is capable of addressing diffusion of atoms, in general it can achieve much longer timescales. For now MAPP is parallelized using message passing interface (MPI). In the “near future” we will release a hybrid version of MAPP using CUDA that will utilize the GPU cores as well as CPU.

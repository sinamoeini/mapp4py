************************
 Performing Simulations
************************

Whether we want to perform a molecular dynamics, energy minimization,
or any other kind of simulation one needs define the system of interest.
In other words inital conditions of the problem must be given. In
addition, prior to running a simulation the force field (inter-atomistic
potential) needs to be given as well namely, the governing equations
of the problem.

.. autosummary::
   :toctree: generated/

   mapp.md.atoms
   mapp.dmd.atoms

class :class:`mapp.md.atoms` and :class:`mapp.dmd.atoms` are containers
for such data, they store the configuration of the system (position of the
atoms and etc.) and force field kind and related parameters


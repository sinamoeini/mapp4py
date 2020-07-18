**********************
Performing Simulations
**********************

Whether we want to perform a molecular dynamics, energy minimization, or any other kind of simulation one needs define the system of interest, prior to execution. In other words inital conditions of the problem must be given. In addition, the force field (inter-atomistic potential) namely, the governing equations of the problem, needs to be given as well.

.. autosummary::
   :toctree: generated/

   mapp4py.md.atoms
   mapp4py.dmd.atoms


class :class:`mapp4py.md.atoms` and :class:`mapp4py.dmd.atoms` are containers for such information, they store the configuration of the system (per atom quantities such as position inside the unitcell, and etc.) and force field kind and relevant parameters


Defining Initial Value
**********************

One of the common ways to give present a system configuration is using a file. There are some static methods provided in both containers that can be utilised to present construct containers from a file. These static method names start with :code:`import_`

MD
--
.. autosummary::
   :toctree: generated/

   mapp4py.md.atoms.import_cfg

DMD
---
.. autosummary::
   :toctree: generated/

   mapp4py.dmd.atoms.import_cfg

Defining Governing Equations
****************************

explain defining forcefield


MD
--
.. autosummary::
   :toctree: generated/
   
   mapp4py.md.atoms.ff_lj
   mapp4py.md.atoms.ff_fs
   mapp4py.md.atoms.ff_eam_funcfl
   mapp4py.md.atoms.ff_eam_setfl
   mapp4py.md.atoms.ff_eam_fs


DMD
---
.. autosummary::
   :toctree: generated/
   
   mapp4py.dmd.atoms.ff_eam_funcfl
   mapp4py.dmd.atoms.ff_eam_setfl
   mapp4py.dmd.atoms.ff_eam_fs




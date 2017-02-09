
import numpy
import mapp
from math import *

mapp.pause_slave_out();


"""
#Kob-Anderson
sim=mapp.cfg_md("configs/KA.cfg")
sim.ff_lj(sigma=[[1.0],[0.8,0.88]],eps=[[1.0],[1.5,0.5]],r_c=[[2.5],[2.0,2.2]],shift=False)
sim.kB=1.0;
sim.create_T(1.0,8569643);
dt=0.05**2/sqrt(3.0);


nst=mapp.md_nst(sim,[[0.0],[1.0,0.0],[0,0,0.0]],0.01,dt);
nst.t_relax_baro=dt*10000
nst.S_dof=[[False],[True,False],[False,False,False]];
nst.nreset=1000;
nst.ntally=100;
nst.run(1000)

nvt=mapp.md_nvt(sim,1.0,dt);
nvt.ntally=100;
nvt.run(1000)

min=mapp.min_cg(sim);
min.run(1000)
nvt.run(1000)
"""


"""
nst=mapp.md_nst(sim,[[0.0],[1.0,0.0],[0,0,0.0]],0.01,dt);
nst.t_relax_baro=dt*10000
nst.S_dof=[[False],[True,False],[False,False,False]];
nst.nreset=1000;
nst.ntally=100;
nst.run(10000)
"""


"""
#Kob-Anderson
sim=mapp.cfg_md("configs/KA.cfg",mpi0)
sim.ff_lj(sigma=[[1.0],[0.8,0.88]],eps=[[1.0],[1.5,0.5]],r_c=[[2.5],[2.0,2.2]],shift=False)


sim=mapp.cfg_md("configs/Cementite.cfg")
# {0:Fe, 1:C}
sim.ff_fs(A=[1.8289905,2.9588787],t1=[[1.0,10.024001],[10.482408,0.0]],
          t2=[[0.504238,1.638980],[3.782595,-7.329211]],
          k1=[[1.237115],[8.972488,22.061824]],
          k2=[[-0.35921],[-4.086410,-17.468518]],
          k3=[[-0.038560],[1.483233,4.812639]],
          r_c_phi=[[3.40],[2.468801,2.875598]],
          r_c_rho=[[3.569745],[2.545937,2.892070]])



sim=mapp.cfg_md("configs/FeBCC.cfg")
# {0:Fe}
sim.ff_fs(A=[1.8289905],t1=[[1.0]],
    t2=[[0.504238]],
    k1=[[1.237115]],
    k2=[[-0.35921]],
    k3=[[-0.038560]],
    r_c_phi=[[3.40]],
    r_c_rho=[[3.569745]])











sim.kB=1.0;
sim.create_T(1.0,57569876);
dt=0.05**2/sqrt(3.0);
# NVT
nvt=mapp.md_nvt(sim,0.01,dt);
nvt.ntally=10;
nvt.run(100)


# NST
nst=mapp.md_nst(sim,[[0.0],[1.0,0.0],[0,0,0.0]],0.01,dt);
nst.t_relax_baro=dt*10000
nst.S_dof=[[False],[True,False],[False,False,False]];
nst.nreset=1000;
nst.ntally=100;
nst.run(10000)


muvt=mapp.md_muvt(sim,1.0,1.0,dt,"P",2411895);
muvt.gas_element="Ni"
muvt.nevery=10;
muvt.nattempts=10;
muvt.ntally=1;
muvt.run(100);

min=mapp.min_cg(sim);
#min.H_dof=[[True],[False,True],[False,False,True]]
#min.affine=True;
min.run(100)
print sim.H
"""









"""
sim=mapp.cfg_dmd("configs/Cu-DMD.cfg")
sim.ff_eam_setfl("potentials/Cu_mishin.eam.alloy")
min=mapp.min_cg(sim);
min.H_dof=[[True],[False,True],[False,False,True]]
min.affine=True;
min.ntally=1;
min.run(100)

min.H_dof=[[False],[False,False],[False,False,False]]
min.affine=False;
min.run(100)
"""





#delta=0.05;
#strain=[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]

#sim=mapp.cfg_dmd(5,"configs/Cu-DMD.cfg")
#sim.ff_eam_setfl("potentials/Cu_mishin.eam.alloy")

sim=mapp.cfg_dmd(3,"configs/FeH-DMD.cfg")
sim.ff_eam_fs("potentials/FeH.eam.fs")
min=mapp.min_cg_dmd(sim);
min.ntally=1;
min.max_dalpha=0.5;
min.run(1000)





"""
min=mapp.min_cg(sim);
min.H_dof=[[True],[False,True],[False,False,True]]
min.affine=True;
print sim.H
min.run(1000)
strain[0][0]-=delta;
strain[1][1]-=delta;
strain[2][2]-=delta;
sim.strain(strain)
min.run(1000)
"""

"""
for i in xrange(0,100):
    strain[0][0]-=delta;
    strain[1][1]-=delta;
    strain[2][2]-=delta;
    sim.strain(strain)
    min.run(0)
"""



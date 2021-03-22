import os, sys, time
import numpy as np
import mapp4py
from mapp4py import md

#####################################################

boltz = 8.617330350e-5
planck = 4.13566766225 * 0.1 * np.sqrt(1.60217656535/1.66053904020)

#####################################################
from mapp4py import mpi
if mpi().rank!=0:
    with open(os.devnull, 'w') as f:
        sys.stdout = f;
################################################################
def mu(p,T):
   return -2.37+0.0011237850013293155*T+0.00004308665175*T*np.log(p)-0.000193889932875*T*np.log(T);

################################################################

muvt = md.muvt(mu(1.0e-2,300.0),300.0,0.1,'H',756055104);
muvt.nevery = 1;
muvt.ntally = 0;

################################################################

nvt=md.nvt(300.0,0.1);
nvt.ntally = 0;

################################################################
sim = md.atoms.import_cfg("configs/Fe_300K.cfg");
sim *= [64, 64, 64];
sim.add_elem("H",1.007940)
sim.ff_eam_setfl("potentials/FeH-sina.eam.alloy");
sim.hP = planck
sim.kB = boltz

np.random.seed(sim.comm_rank+7236498);
frac = 0.001;
def remove_random(id):
    if np.random.random() < frac:
        return False;
    return True;

sim.do(remove_random)

sim.create_temp(300.0,8569643);

muvt.nattempts = 1024*1000;
nsteps = 1;


start = time.time()
muvt.run(sim,nsteps);
if sim.comm_rank == 0:
    print ("time elapsed: %lf seconds" % (time.time()-start))




import time 
import os
import subprocess
proc = subprocess.Popen('rm -rf dumps/*', shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)

from math import sqrt
boltz = 8.617330350e-5
planck = 4.13566766225 * 0.1 * sqrt(1.60217656535/1.66053904020)



import mapp
from mapp import md
mapp.pause_slave_out();

################################################################

min=md.min_cg(e_tol=1.0e-8);
min.ntally=10000;
min.ls=mapp.ls_brent();

################################################################

muvt=md.muvt(-2.33,300,0.1,'H',2411895);
muvt.nevery=100;
muvt.nattempts=10000;
muvt.ntally=1000;
muvt.export=md.export_cfg('dumps/dump',10000)

################################################################



sim=md.atoms.import_cfg("configs/Fe.cfg")
sim.add_elem('H',1.007940)
sim.ff_eam_fs("potentials/FeH.eam.fs")
sim.hP=planck
sim.kB=boltz


min.run(sim,500000)

min.H_dof=[[True],[False,False],[False,False,True]]
min.affine=True

min.run(sim,500000)

sim.create_temp(300.0,8569643);



sim.step=0

start = time.time()
muvt.run(sim,50000000);
print "time elapsed: %lf seconds" % (time.time()-start)


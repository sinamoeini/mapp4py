
#######################################################################################
#    [M] = a.m.u
#    [L] = Angs
#    [U] = [ML^2/T^2] = eV
#
#    therefore
#    [T] = sqrt(a.m.u/eV)*Angs
#
#
#    a.m.u = 1.66053904020e-27 Kg
#    Angs = 1.0e-10 m
#    eV = 1.60217656535e-19 Kg m^2/s^2
#    therefore
#    [T] = sqrt(1.66053904020/1.60217656535)*1.0e-14 s
#
#    kB = 8.617330350e-5 eV/K
#    h = 4.13566766225e-15 eVs
#    h = 4.13566766225 * 0.1 * sqrt(1.60217656535/1.66053904020)  sqrt(eV*a.m.u)*Angs
#######################################################################################

import time 
import os
import subprocess
proc = subprocess.Popen('rm -rf dumps/*', shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)

from math import sqrt
from math import exp
from math import log
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

nan=float('nan')
nst=md.nst([[0.0],[nan,0.0],[nan,nan,nan]],300,0.1)
nst.ntally=10000
nst.export=md.export_cfg('dumps/relax', 100000)

################################################################

sim=md.atoms.import_cfg("configs/poly.cfg")
sim.comm_dims=[8,4,4]
sim.add_elem('H',1.007940)
sim.ff_eam_setfl("potentials/NiAlH_jea.eam.alloy")
sim.hP=planck
sim.kB=boltz

################################################################

def get_muvt(mu, nt):
    muvt=md.muvt(mu,300.0,0.1,'H',2411895);
    muvt.nevery=1000;
    muvt.nattempts=10000;
    muvt.export=md.export_cfg('dumps/dump', nt)
    return muvt;

################################################################

min.run(sim,500000)
sim.create_temp(300.0,8569643);
sim.step=0
nst.run(sim,500000)
sim.step=0




mu0=-2.547;
mu1=-2.447;
no=50;


beta=1.0/(boltz*300.0);
X=exp(2.0*beta*(mu1-mu0))-1.0;
d=1.0/(no*1.0);
k=0.0;
for i in range(0,no):
    mu=mu0+0.5*log(1.0+k*X)/beta;
    muvt=get_muvt(mu, 20000)
    muvt.run(sim,200000);
    k+=d;

for i in range(no,-1,-1):
    mu=mu0+0.5*log(1.0+k*X)/beta;
    muvt=get_muvt(mu, 20000)
    muvt.run(sim,200000);
    k-=d;



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
boltz = 8.617330350e-5
planck = 4.13566766225 * 0.1 * sqrt(1.60217656535/1.66053904020)



import mapp4py
from mapp4py import dmd



sim=dmd.atoms.import_cfg(5,'configs/Al_DMD.cfg')


def fix_dof(id,x_dof,c_dof,alpha_dof,alpha):
    if id>3197:
        x_dof[0]=False;
        x_dof[1]=False;
        x_dof[2]=False;
        c_dof[0]=False;

sim.do(fix_dof)

sim.ff_eam_setfl("potentials/NiAl.eam.alloy",[3.0],[2.0])
sim.hP=planck
sim.kB=boltz
sim.temp=300.0


min=dmd.min_cg();
min.ntally=1000;
min.e_tol=1.0e-10;
min.ls=mapp4py.ls_brent();
min.ntally=0;
#min.export=dmd.export_cfg('dumps/dump',1000,sort=True,extra_vecs=['f','f_alpha'])
min.run(sim,10000)

def fix_dof(id,alpha_dof,c):
    if id>3197:
        alpha_dof[0]=False;
    if id==0:
        c[0]=0.0;

sim.do(fix_dof)

chem=dmd.bdf();
chem.a_tol=1.0e-6;
chem.ntally=1;
chem.max_ngmres_iters=20;
chem.max_nnewton_iters=100;
chem.max_nsteps=10000;
chem.nreset=1;
chem.export=dmd.export_cfg('dumps/dump',10,sort=True,extra_vecs=['c_d'])

sim.step=0
chem.run(sim,10.0**(30));







.. MAPP documentation master file, created by
   sphinx-quickstart on Tue Feb 14 00:27:45 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MAPP's Documentation!
================================

:Page Status: Incomplete
:Last Reviewed: |today|

MAPP is a parallel atomistic simulation package written entirely in C++, it is presented in form of a `python <http://www.python.org>`_ package in order to facilitate usage. Some of the main features of MAPP are:

  * Molecular Dynamics (MD)
  * Diffusive Molecular Dynamics (DMD) :cite:`ju_li_diffusive_2011` canonical ensemble
  * DMD isothermoal-isostress ensemble
  * Sophisticated energy minimization (CG and l-BFGS) and line search (brent, golden section and back track) algorithms
  * Parallel Grand Canonical Monte Carlo (pGCMC) for many body potentials such as EAM


At the moment the paralleization algorithms are done using standard `message passing interface (MPI) library <https://en.wikipedia.org/wiki/Message_Passing_Interface>`_. mapp4py package is consist of two modules

.. autosummary::
   :toctree: generated/

   mapp4py
   mapp4py.md
   mapp4py.dmd




.. toctree::
   :maxdepth: 2

   install
   simul
   md
   dmd
   min
   ff

References
----------
.. bibliography:: refs.bib
   :filter: docname in docnames
   :style: unsrt

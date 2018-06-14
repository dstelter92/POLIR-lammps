/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(POLIR/CHARGE/ATOM,ComputePolirChargeAtom)

#else

#ifndef LMP_COMPUTE_POLIR_CHARGE_ATOM_H
#define LMP_COMPUTE_POLIR_CHARGE_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePolirChargeAtom : public Compute {
 public:
  ComputePolirChargeAtom(class LAMMPS *, int, char **);
  ~ComputePolirChargeAtom();
  void init();
  void compute_peratom();
  double memory_usage();

 protected:
  int me,np;
  int ifix,icompute;
  int nmax,nbonds;
  int *typeH,*typeO;
  double *c1,*c3,*qeH,*reOH;
  double *qpolir;
  double **indx2bond;
  double *lbond;
  double *qH;
  double *qO;

  char *fix_polir,*compute_lbond;
  
  void allocate();
  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Per-atom charge was not tallied on needed timestep

You are using a thermo keyword that requires potentials to
have tallied charge, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

*/

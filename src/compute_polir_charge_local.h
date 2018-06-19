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

ComputeStyle(polir/charge/local,ComputePolirChargeLocal)

#else

#ifndef LMP_COMPUTE_POLIR_CHARGE_LOCAL_H
#define LMP_COMPUTE_POLIR_CHARGE_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePolirChargeLocal : public Compute {
 public:
  ComputePolirChargeLocal(class LAMMPS *, int, char **);
  ~ComputePolirChargeLocal();
  void init();
  void compute_local();
  void allocate();
  double memory_usage();

 private:
  int nlocal,typeH,typeO;
  double c1,c3,qeH,reOH;
  double *qpolir;
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

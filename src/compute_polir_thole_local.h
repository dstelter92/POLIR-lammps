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

ComputeStyle(POLIR/THOLE/LOCAL,ComputePolirTholeLocal)

#else

#ifndef LMP_COMPUTE_POLIR_THOLE_LOCAL_H
#define LMP_COMPUTE_POLIR_THOLE_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePolirTholeLocal : public Compute {
 public:
  ComputePolirTholeLocal(class LAMMPS *, int, char **);
  ~ComputePolirTholeLocal();
  void init();
  void init_list(int, class NeighList *);
  void compute_local();
  double memory_usage();

 private:
  class NeighList *list;
  int nvalues,nmax;
  double A_const,alphaH,alphaO;

  double **thole;
  double *damping;
  int *typeH,*typeO;
  double *CD_intra_OH,*CD_intra_HH,*DD_intra_OH,*DD_intra_HH;
  double *CC_inter,*CD_inter,*DD_inter;

  char *fix_polir;

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

E: Compute bond/local used when bonds are not allowed

The atom style does not support bonds.

E: Invalid keyword in compute bond/local command

Self-explanatory.

E: No bond style is defined for compute bond/local

Self-explanatory.

E: Sanity check on 3 energy components failed

UNDOCUMENTED

*/

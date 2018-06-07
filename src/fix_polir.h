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

#ifdef FIX_CLASS

FixStyle(polir,FixPolir)

#else

#ifndef LMP_FIX_POLIR_H
#define LMP_FIX_POLIR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPolir : public Fix {
 public:
  FixPolir(class LAMMPS *, int, char **);
  ~FixPolir();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void pre_force(int);
  void post_force(int);
  void end_of_step();

 private:
  Compute *compute_pca;

 protected:
  int typeH;                        // type index of H
  int typeO;                        // type index of O
  double qeH;                       // equilibrium charge
  double reOH;                      // equilibrium O-H bond length
  double c1;                        // 1st coeff for charge-dep bonds
  double c3;                        // 2nd coeff for charge-dep bonds
  
  int q_compute_id;
  char * id_q;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/

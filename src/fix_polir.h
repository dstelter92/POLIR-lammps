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
  void allocate();
  double memory_usage();

 private:
  Compute *compute_pca;
  Compute *compute_thole;

 protected:
  // INPUT VARIABLES
  int typeH;                        // type index of H
  int typeO;                        // type index of O
  double qeH;                       // equilibrium charge
  double reOH;                      // equilibrium O-H bond length
  double uH;                        // H dipole along HOH bisector pointing to H
  double uO;                        // O dipole in HOH plane along O->H bond
  double alphaH;                    // polarizibility for type H
  double alphaO;                    // polarizibility for type O

  // POLIR PARAMETERS
  double c1;              // 1st coeff for charge-dep bonds
  double c3;              // 2nd coeff for charge-dep bonds
  double CD_intra_OH;     // Thole damping coeff, charge-dipole intramolecular O-H
  double CD_intra_HH;     // Thole damping coeff, charge-dipole intramolecular H-H
  double DD_intra_OH;     // Thole damping coeff, dipole-dipole intramolecular O-H
  double DD_intra_HH;     // Thole damping coeff, dipole-dipole intramolecular H-H
  double CC_inter;        // Thole damping coeff, charge-charge intermolecular
  double CD_inter;        // Thole damping coeff, charge-dipole intermolecular
  double DD_inter;        // Thole damping coeff, dipole-dipole intermolecular
  
  int polir_output,count;
  int q_compute_id, thole_compute_id;
  char *id_q,*id_thole;
  char *stypeH,*stypeO,*sqeH,*sreOH,*suH,*suO,*salphaH,*salphaO;

  double *global_vector;
  double *charges;

  int me_universe,nworlds,iworld,me;
  int nlocal,nmax;
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
